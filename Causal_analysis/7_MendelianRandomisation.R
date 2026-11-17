# Loading libraries
library(TwoSampleMR)
library(dplyr)

# Reading in dataset
harmonised_dataset <- read.csv('harmonised_dataset.csv')

## Steiger filtering: confirms direction is protein -> AD, not reverse

# Adding sample sizes onto the exposure and outcome datasets for test
harmonised_dataset$samplesize.exposure <- 54219 
harmonised_dataset$samplesize.outcome <- 788989 # got these from the papers

# Steiger filtering
steiger_results <- directionality_test(harmonised_dataset) # assesses whether a SNP explains more variance in the outcome data than in the exposure data. If that is the case, the causal effect might be the other way around
print(steiger_results)

# Filter to keep only correctly-directed
harmonised_filtered <- harmonised_dataset %>%
  left_join(steiger_results %>% select(id.exposure, correct_causal_direction),
            by = 'id.exposure') %>%
  filter(correct_causal_direction == TRUE) 

# Assessing how many SNPs/proteins have been lost
cat('SNPs after Steiger filtering:', nrow(harmonised_filtered))
cat('Proteins remaining after Steiger filtering:', length(unique(harmonised_filtered$exposure)))

# Documenting how many proteins only have 1 SNP -> means that the sensitivity analyses don't really work
instrument_count <- harmonised_filtered %>%
  filter(mr_keep == TRUE) %>%
  group_by(exposure) %>%
  summarise(n_snps = n()) %>%
  mutate(single_instrument = n_snps == 1)

cat('Proteins with only 1 SNP:', sum(instrument_count$single_instrument))

# Calculating mean, iqr and range of snps per protein
snp_summary <- harmonised_filtered %>%
  filter(mr_keep == TRUE) %>%
  group_by(exposure) %>%
  summarise(n_snps = n(), .groups = 'drop')

cat('Median SNPs per protein:', median(snp_summary$n_snps))
cat('IQR:', quantile(snp_summary$n_snps, c(0.25, 0.75)))
cat('Range:', min(snp_summary$n_snps), max(snp_summary$n_snps))

# Calculating F stat summary stats
harmonised_filtered <- harmonised_filtered %>%
  mutate(f_stat = (beta.exposure / se.exposure)^2)

mean_f <- mean(harmonised_filtered$f_stat[harmonised_filtered$mr_keep == TRUE])
min_f <- min(harmonised_filtered$f_stat[harmonised_filtered$mr_keep == TRUE])
max_f <- max(harmonised_filtered$f_stat[harmonised_filtered$mr_keep == TRUE])

cat('Mean F-statistic:', round(mean_f, 1), '\n')
cat('Range:', round(min_f, 1), '–', round(max_f, 1), '\n')

## Calculating r squared
# getting r squared from stats
harmonised_filtered <- harmonised_filtered %>%
  mutate(r2_snp = 2 * eaf.exposure * (1 - eaf.exposure) * beta.exposure^2)

# total r squared per protein (sum across that protein's instruments)
r2_per_protein <- harmonised_filtered %>%
  filter(mr_keep == TRUE) %>%
  group_by(exposure) %>%
  summarise(total_r2 = sum(r2_snp), .groups = 'drop')

# Summary stats
mean_r2 <- mean(r2_per_protein$total_r2) * 100
min_r2  <- min(r2_per_protein$total_r2) * 100
max_r2  <- max(r2_per_protein$total_r2) * 100

cat('Mean R^2 per protein:', round(mean_r2, 2), '%\n')
cat('Range:', round(min_r2, 2), '% –', round(max_r2, 2), '%\n')


## Function that runs IVW, MR-Egger, and Weighted Median by default - MAIN MR TEST
# handles the mr_keep column internally
mr_results <- mr(harmonised_filtered)

# Printing results 
print(mr_results)

## Multiple testing correction for the IVW results
# Combine IVW and Wald ratio into primary results
primary_results <- mr_results %>% 
  filter(method %in% c('Inverse variance weighted', 'Wald ratio')) # keeping only these  two methods

# Apply multiple testing correction across all proteins tested
primary_results$pval_fdr <- p.adjust(primary_results$pval, method = 'BH') #BH = Benjamin-Hochberg

# Count results
n_tested <- length(unique(primary_results$exposure))
n_uncorrected <- length(unique(primary_results$exposure[primary_results$pval < 0.05]))
n_fdr <- length(unique(primary_results$exposure[primary_results$pval_fdr < 0.05]))

# Reporting results
cat('Proteins tested:', n_tested, '\n')
cat('  - via IVW (2+ SNPs):', sum(primary_results$method == 'Inverse variance weighted'), '\n')
cat('  - via Wald ratio (1 SNP):', sum(primary_results$method == 'Wald ratio'), '\n')
cat('Significant at p < 0.05 (uncorrected):', n_uncorrected, '\n')
cat('Significant after FDR correction:', n_fdr, '\n')

## Cochran's Q heterogeneity test → tests for consistency within SNPs
heterogeneity <- mr_heterogeneity(harmonised_filtered)

# Reporting results
cat('Proteins with significant heterogeneity (Q p < 0.05):', 
    sum(heterogeneity$Q_pval < 0.05, na.rm = TRUE), '\n')

## MR-Egger intercept test — identifies pleiotropy
pleiotropy <- mr_pleiotropy_test(harmonised_filtered)

# Reporting results
cat('Proteins with significant directional pleiotropy (Egger intercept p < 0.05):', 
    sum(pleiotropy$pval < 0.05, na.rm = TRUE), '\n')

## Testing for consistency across causal methods from mr analysis
consistency <- mr_results %>%
  filter(method %in% c('Inverse variance weighted', 'MR Egger', 'Weighted median')) %>%
  group_by(exposure) %>%
  summarise(
    n_methods_tested = sum(!is.na(pval)),
    n_methods_significant = sum(pval < 0.05, na.rm = TRUE),
    all_consistent = (n_methods_tested >= 3) & (n_methods_significant == n_methods_tested)
  )

# Reporting results
cat('Proteins with consistently significant results:', 
    sum(consistency$all_consistent == TRUE), '\n')

## Combining all criteria
high_conf <- primary_results %>%
  filter(pval_fdr < 0.05) %>% # FDR p-value significant
  left_join(instrument_count, by = 'exposure') %>% # merging with dataset documenting one SNP proteins
  left_join(consistency %>% # merging with dataset documenting consistent results
              select(exposure, all_consistent), by = 'exposure') %>% # keeping only columns needed
  left_join(heterogeneity %>% filter(method == 'Inverse variance weighted') %>% # merging with dataset documenting heterogeneity of IVW
              select(exposure, Q_pval), by = 'exposure') %>% # keeping only columns needed
  left_join(pleiotropy %>% # merging with dataset documenting pleiotropy
              select(exposure, egger_pval = pval), by = 'exposure') # keeping only columns needed

# Identifying which criteria each protein passes in new columns in the dataset
high_conf$pass_multi_snp <- high_conf$n_snps >= 2 # Multiple SNPs
high_conf$pass_consistency <- high_conf$all_consistent == TRUE | is.na(high_conf$all_consistent) # Consistently passes MR
high_conf$pass_heterogeneity <- high_conf$Q_pval > 0.05 | is.na(high_conf$Q_pval) # Not hetergenous
high_conf$pass_pleiotropy <- high_conf$egger_pval > 0.05 | is.na(high_conf$egger_pval) # No pleiotropy

# Reporting results
cat('FDR-significant:', nrow(high_conf), '\n')
cat('Multi-SNP:', sum(high_conf$pass_multi_snp), '\n')
cat('Consistent across methods:', sum(high_conf$pass_consistency, na.rm = TRUE), '\n')
cat('No heterogeneity:', sum(high_conf$pass_heterogeneity, na.rm = TRUE), '\n')
cat('No directional pleiotropy:', sum(high_conf$pass_pleiotropy, na.rm = TRUE), '\n')

# Calculating OR instead as more interpretable
high_conf$OR <- exp(high_conf$b)
high_conf$OR_lower <- exp(high_conf$b - 1.96 * high_conf$se)
high_conf$OR_upper <- exp(high_conf$b + 1.96 * high_conf$se)

## Assessing number of OR effects
cat('Risk increasing direction effect:', sum(high_conf$OR > 1))
cat('Risk decreasing direction effect:', sum(high_conf$OR < 1))

## Assessing biggest OR effects
biggest_OR <- which.max(ad$OR)
cat('Highest ')

# Saving results
fwrite(high_conf, "mr_results.csv")
fwrite(mr_results, "mr_results_full.csv")

