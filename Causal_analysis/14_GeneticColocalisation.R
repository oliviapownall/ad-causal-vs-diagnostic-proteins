# Loading libraries
library(coloc)
library(dplyr)
library(data.table)
library(ggplot2)
library(patchwork)

# reading in harmonised dataset
harmonised <- read.csv('harmonised_genetic_coloc.csv')

# Function that will run coloc for one protein
run_coloc_protein <- function(protein_data, protein_name,
                              n_protein = 54219,        # UKB-PPP sample size
                              n_cases = 111326,          # Bellenguez cases
                              n_controls = 677663) {     # Bellenguez controls
  # Compute MAF (coloc needs MAF, not EAF)
  protein_data$maf <- pmin(protein_data$eaf.exposure, 
                           1 - protein_data$eaf.exposure)
  
  # Sample size and case proportion for AD outcome
  n_total_outcome <- n_cases + n_controls
  case_proportion <- n_cases / n_total_outcome
  
  # Build coloc input list for protein
  D1 <- list(
    beta    = protein_data$beta.exposure,
    varbeta = protein_data$se.exposure^2,
    snp     = protein_data$SNP,
    type    = "quant",
    N       = n_protein,
    MAF     = protein_data$maf
  )
  
  # Build coloc input list for AD
  D2 <- list(
    beta    = protein_data$beta.outcome,
    varbeta = protein_data$se.outcome^2,
    snp     = protein_data$SNP,
    type    = "cc",
    s       = case_proportion,
    N       = n_total_outcome,
    MAF     = protein_data$maf
  )
  
  # Run coloc.abf (default priors)
  res <- coloc.abf(dataset1 = D1, 
                   dataset2 = D2,
                   p1  = 1e-4, 
                   p2  = 1e-4, 
                   p12 = 1e-5)
  
  # Find lead pQTL (smallest exposure p-value) for documentation
  lead_idx <- which.min(protein_data$pval.exposure)
  
  # Return summary
  data.frame(
    protein           = protein_name,
    lead_cis_pqtl     = protein_data$SNP[lead_idx],
    n_snps_in_coloc   = nrow(protein_data),
    PP.H0             = res$summary["PP.H0.abf"],
    PP.H1             = res$summary["PP.H1.abf"],
    PP.H2             = res$summary["PP.H2.abf"],
    PP.H3             = res$summary["PP.H3.abf"],
    PP.H4             = res$summary["PP.H4.abf"],
    interpretation    = case_when(
      res$summary["PP.H4.abf"] >= 0.80 ~ "Strong colocalisation (H4 ≥ 0.80)",
      res$summary["PP.H4.abf"] >= 0.60 ~ "Moderate colocalisation (H4 ≥ 0.60)",
      res$summary["PP.H3.abf"] >= 0.80 ~ "Distinct causal variants (H3 ≥ 0.80)",
      TRUE                              ~ "Inconclusive"
    ),
    coloc_status      = ifelse(nrow(protein_data) < 50, 
                               "warn_low_n_snps", "ok"),
    stringsAsFactors  = FALSE
  )
}

# Grabbing proteins
proteins <- unique(harmonised$exposure)

# creating empty list to store results in?
results_list <- list()

# Looping through proteins and runnign function
for (protein in proteins) {
  cat("Processing", protein)
  
  # Subset to this protein
  protein_data <- harmonised[harmonised$exposure == protein, ]
  
  # Run coloc 
  result <- run_coloc_protein(protein_data, protein)
  results_list[[protein]] <- result
  
}

coloc_results <- bind_rows(results_list) %>%
  arrange(desc(PP.H4))

# Write up results
fwrite(coloc_results, "genetic_coloc_results.csv")

### Further assessment of APOE

ad <- fread("ad_filtered_genetic_coloc.tsv")

# Filter to APOE locus (chr19, ±1 Mb of rs429358 at 44908684)
ad_apoe <- ad[ad$chromosome == 19 & 
                ad$base_pair_location >= 43900000 & 
                ad$base_pair_location <= 45900000, ]

cat("Bellenguez SNPs in APOE region:", nrow(ad_apoe), "\n")

# Assessing for known APOE isoform coders in data
"rs429358" %in% ad$variant_id
"rs7412" %in% ad$variant_id

# Assessing what was counted as the strongest pQTL score
ad_apoe[order(p_value)][1:10]



##################### RERUNNING COLOC WITH KUNKLE DATA ####################################

# Reading in dataset
harmonised_kunkle_final <- read.csv('harmonised_kunkle.csv')

# Specifying cases and controls
n_cases_kunkle <- 21982
n_controls_kunkle <- 41944

results_kunkle <- list()
for (protein in c("APOE", "APOC1")) {
  protein_data <- harmonised_kunkle_final[harmonised_kunkle_final$exposure == protein, ]
  
  result <- run_coloc_protein(
    protein_data, protein,
    n_protein = 54219,
    n_cases = n_cases_kunkle,
    n_controls = n_controls_kunkle
  )
  
  results_kunkle[[protein]] <- result
  cat(protein, ": PP.H4 =", round(result$PP.H4, 3), 
      "| PP.H3 =", round(result$PP.H3, 3), 
      "| n_snps =", result$n_snps_in_coloc, "\n")
}

coloc_kunkle_results <- bind_rows(results_kunkle)
print(coloc_kunkle_results)

# Writing up results 

fwrite(coloc_kunkle_results, "genetic_coloc_results_apoe.csv")