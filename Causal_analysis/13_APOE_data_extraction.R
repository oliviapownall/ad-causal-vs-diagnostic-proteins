# Loading libraries
library(coloc)
library(dplyr)
library(data.table)
library(ieugwasr)
library(TwoSampleMR)

# Opening GWAS using token
Sys.setenv(OPENGWAS_JWT = "")

# Accessing Kunkle 2019 for the APOE region

kunkle_apoe <- associations(
  variants = "19:44400000-46400000",  # GRCh37 coordinates (Because it's 2019 it uses an older dataset)
  id = "ieu-b-2"  # Kunkle 2019 Stage 1 on IEU OpenGWAS
)

cat("Kunkle SNPs in APOE region:", nrow(kunkle_apoe), "\n")
cat("rs429358 present?", "rs429358" %in% kunkle_apoe$rsid, "\n")
cat("rs7412 present?", "rs7412" %in% kunkle_apoe$rsid, "\n")

# Reformat kunkle to match my previously made pipeline
kunkle_clean <- kunkle_apoe %>%
  transmute(
    SNP = rsid,
    beta = beta,
    se = se,
    effect_allele = ea,
    other_allele = nea,
    eaf = eaf,
    pval = p
  )

# Reading in exposure dataset
exposure_raw <- read.csv('snps_genetic_coloc_rsids.csv')

# Subsetting to just APOE and APOC1 (APOE gene associated proteins)
exposure_apoe <- exposure_raw[exposure_raw$source_protein %in% c("APOE", "APOC1"), ]


# Reformatting (same as before)
exposure_dat_apoe <- format_data(
  exposure_apoe,
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

outcome_dat_kunkle <- format_data(
  kunkle_clean,
  type = "outcome",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval"
)

# Harmonising dataset
harmonised_kunkle <- harmonise_data(
  exposure_dat = exposure_dat_apoe,
  outcome_dat = outcome_dat_kunkle,
  action = 2
)

# Filtering for those to keep
harmonised_kunkle_final <- harmonised_kunkle[harmonised_kunkle$mr_keep == TRUE, ]

# Sanity check
cat("APOE SNPs after Kunkle harmonisation:", 
    sum(harmonised_kunkle_final$exposure == "APOE"))
cat("APOC1 SNPs after Kunkle harmonisation:", 
    sum(harmonised_kunkle_final$exposure == "APOC1"))
cat("rs429358 present in harmonised data?", 
    "rs429358" %in% harmonised_kunkle_final$SNP)
  
  
fwrite(harmonised_kunkle_final, "harmonised_kunkle.csv")
  

  
  