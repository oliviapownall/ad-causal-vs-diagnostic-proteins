## Had to filter AD_SNPs in the terminal- ran the following: 

## Extracting rsIDs from the protein cis-pQTL file (column 18):
# awk -F',' 'NR>1 {print $18}' snps_genetic_coloc_rsids.csv | sort -u > protein_rsids.txt

# Stream-filter the Bellenguez AD GWAS file (Variant_id (column 1 in AD file) is matched against the protein rsID list)
# awk -F'\t' '
#   NR==FNR {ids[$1]; next}        # First file: load protein rsIDs into hash
#   FNR==1 {print; next}            # Print AD header row
#   $1 in ids                       # Print AD rows where variant_id matches
# ' protein_rsids.txt AD_SNPs.tsv > ad_filtered.tsv
# Output: ad_filtered.tsv

## Importing libraries
library(data.table)
library(dplyr)

## Load the AD-filtered file
ad <- fread("ad_filtered.tsv")

## Assessing number of SNPs remaining
cat("AD rows after filtering:", nrow(ad))

## Adding protein data back onto ad dataset
# Reading files
prot <- fread("snps_genetic_coloc_rsids.csv")
ad_merged <- merge(
  ad,
  prot[, .(rsid, source_protein)],
  by.x = "variant_id",
  by.y = "rsid",
)

# Checking all proteins remain
cat("Proteins remaining after filtering:", length(unique(ad_merged$source_protein)))

# Creating a df of separate SNP counts for each protein
count(ad_merged, source_protein)