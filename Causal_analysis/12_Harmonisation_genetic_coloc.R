# Loading libraries
library(TwoSampleMR)
library(dplyr)
library(data.table)

# Reading in datasets
exposure_raw <- read.csv('snps_genetic_coloc_rsids.csv')
outcome_raw <- read.delim('ad_filtered_genetic_coloc.tsv')



# Formatting data to run under the harmonise_data function
# Exposure
exposure_dat <- format_data(
  exposure_raw,
  type = "exposure",
  snp_col = "rsid",             
  beta_col = "BETA",               
  se_col = "SE",                   
  effect_allele_col = "ALLELE1",   
  other_allele_col = "ALLELE0",    
  eaf_col = "A1FREQ",         
  pval_col = "LOG10P",   
  log_pval = TRUE, 
  phenotype_col = "source_protein"
)

# Outcome
outcome_dat <- format_data(
  outcome_raw,
  type = "outcome",
  snp_col = "variant_id",              
  beta_col = "beta",
  se_col = "standard_error",          
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "effect_allele_frequency",
  pval_col = "p_value"         
)

# Harmonising datasets -> handles palindromic SNPs

harmonised <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
  action = 2  # Removes palindromic SNPs with mean allele frequency of > 0.42
)

# Filtering for only the SNPs the function has coded to keep

harmonised_final <- harmonised[harmonised$mr_keep == TRUE, ]

cat("Total SNPs:", nrow(harmonised_final))
cat("Proteins remaining after filtering:", length(unique(harmonised_final$exposure)))
count(harmonised_final, exposure)

# Writing to a file

fwrite(harmonised_final, 'harmonised_genetic_coloc.csv')
