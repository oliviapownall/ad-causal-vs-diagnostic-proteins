# Loading libraries
library(data.table)
library(ieugwasr)
library(dplyr)

# Read in the dataset with rsIDs
snps_cis <- fread("snps_cis_with_rsids.csv")

# Assessing duplicates: find all protein-SNP pairs that appear more than once - it will be the allele differences why there are many
multi_allelic_counts <- snps_cis[, .N, by = .(source_protein, rsid)][N > 1]

cat("Number of SNP-Protein pairs with multiple versions:", nrow(multi_allelic_counts), "\n")
cat("Total rows involved in these duplicates:", sum(multi_allelic_counts$N), "\n")

# We need to only keep the most significant association
# The reason we are not specifically just picking the alleles that are in the outcome data is to avoid cherry-picking which could bias results

# Sorting by protein, then rsid, then pval (largest log p-value first)
setorder(snps_cis, source_protein, rsid, -LOG10P)

# Keeping only the first occurrence for each SNP-Protein pair (aka most significant due to the sorting)
# Since we sorted by pval, the 'first' one is the most significant
snps_unique <- unique(snps_cis, by = c("source_protein", "rsid"))

# Double checking all duplicates have been removed
multi_allelic_counts <- snps_unique[, .N, by = .(source_protein, rsid)][N > 1]

cat("Number of SNP-Protein pairs with multiple versions:", nrow(multi_allelic_counts), "\n")
cat("Total rows involved in these duplicates:", sum(multi_allelic_counts$N), "\n")

# Assessing SNPs and proteins remaining after removing duplicate protein-SNP pairs
cat("SNPs remaining after duplicate removal:", nrow(snps_unique), "\n")
cat("Proteins remaining after duplicate removal:", uniqueN(snps_unique$source_protein), "\n")

## Excluding MHC region due to extensive linkage disequilibrium
snps_unique <- snps_unique %>%
  filter(!(CHROM == 6 & GENPOS >= 26000000 & GENPOS <= 34000000))

n_after <- nrow(snps_unique)

cat('SNPs remaining after MHC region exclusion:', n_after, '\n')
proteins_after <- length(unique(snps_unique$source_protein))
cat('Proteins remaining:', proteins_after, '\n')

# Creates a smaller dataset with just the columns that ld_clump needs
snps_formatted <- snps_unique[, .(
  rsid,
  pval = 10^(-LOG10P),
  source_protein
)]

# Clumping per protein

# Run the clumping and store results in a list column
clump_list <- snps_formatted[, .(
  result = list(
    tryCatch(
      {
        # Skip clumping if only 1 SNP exists
        if (.N == 1) {
          .SD 
        } else {
          # If the protein has 2 or more SNPs, attempt clumping
          ld_clump(
            dat = .SD, # specifies for this set of data (found by 'by = source_protein')
            clump_r2 = 0.001, # for any two SNPs with r^2 >0.001 ,the less significant one is dropped
            clump_kb = 10000, # the window size for location-based LD (1Mb)
            clump_p1 = 1, # keep all SNPs as potential SNP chosen to represent a region
            clump_p2 = 1, # keep all SNPs for comparison - this can be done. as we have already filtered for significant SNPs
            plink_bin = "plink", # link to how to access the reference set (downloaded previously)
            bfile = "EUR" # European reference set of correlations
          )
        }
      },
      error = function(e) {
        message("Skipping protein ", .BY, ": ", conditionMessage(e)) # Sometimes it can fail becasue the SNPs for one protein are not in the reference panel
        return(NULL) 
      }
    )
  )
), by = source_protein]

# Name the list items so rbindlist knows which protein is which
# We use setNames to pair the results with the protein names
named_results <- setNames(clump_list$result, clump_list$source_protein)

# Use idcol = "source_protein" to bring that column back to life
snps_clumped <- rbindlist(named_results, fill = TRUE, idcol = "source_protein")

# Assessing results
cat("SNPs remaining after clumping:", nrow(snps_clumped), "\n")
cat("Proteins remaining after clumping:", uniqueN(snps_clumped$source_protein), "\n")

# Merge back with original dataset
final_instruments <- merge(
  snps_clumped[, .(rsid, source_protein)],
  snps_unique,
  by = c("rsid", "source_protein")
)

# Removing the not needed columns
prot_final <- final_instruments[, c('rsid', 'BETA', 'LOG10P', 'SE', 'ALLELE1', 'ALLELE0', 'A1FREQ', 'source_protein')]

# renaming to what TwoSampleMR will expect
prot_final <- prot_final %>%
  rename(
    SNP = rsid,
    `pval(log)` = LOG10P,
    beta = BETA,
    se = SE,
    effect_allele = ALLELE1,
    other_allele = ALLELE0,
    eaf = A1FREQ,
    exposure = source_protein
  )

# Assessing for SNP-protein duplicates
unique_pairs_with_dupes <- prot_final %>%
  group_by(SNP, exposure) %>%
  filter(n() > 1) %>%
  summarise(n_snps = n_distinct(SNP)) %>%
  pull(n_snps)

cat("Number of SNP-Protein pairs that have duplicates:", unique_pairs_with_dupes)
# No duplicates

# Save final instruments - exposure dataset
fwrite(final_instruments, "FinalSNPsProteins.csv")
