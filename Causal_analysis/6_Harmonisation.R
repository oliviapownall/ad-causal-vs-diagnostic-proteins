# Loading libraries
library(TwoSampleMR)
library(dplyr)

# Reading in datasets
exposure_raw <- read.csv('FinalSNPsProteins.csv')
outcome_raw <- read.csv('finalSNPsAD.csv')

# Formatting data to run under the harmonise_data function

# Exposure
exposure_dat <- format_data(
  exposure_raw,
  type = "exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf",
  phenotype_col = "exposure"
)

# Outcome
outcome_dat <- format_data(
  outcome_raw,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  eaf_col = "eaf"
)

# Harmonising datasets -> handles palindromic SNPs

harmonised <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outcome_dat,
  action = 2  # Removes palindromic SNPs with mean allele frequency of > 0.42
)

# Filtering for only the SNPs the function has coded to keep

harmonised_final <- harmonised[harmonised$mr_keep == TRUE, ]

cat("SNPs after harmonisation:", nrow(harmonised_final), "\n")
cat("Proteins remaining after harmonisation:", uniqueN(harmonised_final$exposure), "\n")

# Writing to a file

fwrite(harmonised_final, "harmonised_dataset.csv")
