# ad-causal-vs-diagnostic-proteins
Code for 'Mendelian Randomization and Explainable Machine Learning Show Divergent Causal and Diagnostic Value Among Plasma Proteins in Alzheimer’s Disease': Mendelian randomisation of 1,875 plasma proteins in Alzheimer’s disease (UKB-PPP), with explainable machine learning validation in the GNPC cohort.

## Overview

Plasma proteins with genetic evidence of a causal role in Alzheimer's disease (AD)
are often prioritised as both drug targets and blood-based biomarkers. This project
tests whether causal evidence predicts diagnostic value, in two stages:

1. **Causal analysis:** two-sample Mendelian randomisation (MR) of 1,875 plasma
   proteins, using cis-pQTLs from the UK Biobank Pharma Proteomics Project
   (UKB-PPP) and AD GWAS summary statistics (Bellenguez et al., 2022), followed by
   Bayesian colocalisation. MR-significant proteins are grouped into three tiers of
   causal confidence.
2. **Discriminative analysis:** the diagnostic value of each tier is evaluated in an
   independent cohort from the Global Neurodegeneration Proteomics Consortium (GNPC),
   using machine learning and SHapley Additive exPlanations (SHAP).

## Repository structure

```
├── Causal_analysis/          # MR and colocalisation (summary-level data)
└── Discriminative_analysis/  # Machine learning in GNPC (individual-level data)
```

No data are included in this repository. See [Data availability](#data-availability).

## Causal analysis

Run in numbered order:

| Step | File | Description |
|---|---|---|
| 1 | `1_extracting_SNPs.ipynb` | Downloads UKB-PPP (European ancestry) protein GWAS from Synapse and extracts genome-wide significant SNPs |
| 2 | `2_ExtractingCispQTLs.R` | Restricts to cis-pQTLs (within 1 Mb of the protein-coding gene) |
| 3 | `3_matching_rsIDs.ipynb` | Maps variant positions to rsIDs via the Ensembl REST API |
| 4 | `4_Clumping.R` | LD clumping using PLINK and the 1000 Genomes European reference panel |
| 5 | `5_extracting_AD_SNPs.ipynb` | Extracts the instrument SNPs from the AD GWAS |
| 6 | `6_Harmonisation.R` | Harmonises exposure and outcome alleles |
| 7 | `7_MendelianRandomisation.R` | Runs MR (IVW, MR-Egger, weighted median, simple and weighted mode) and sensitivity analyses |
| 8 | `8_GO_analysis.md` | Gene Ontology enrichment (run via the PANTHER web interface; settings documented) |
| 9 | `9_Genetic_coloc_extraction.ipynb` | Extracts regional SNPs for the MR-significant proteins |
| 10 | `10_Matching_rsIDs_genetic_coloc.ipynb` | Maps regional SNPs to rsIDs |
| 11 | `11_ExtractingADSNPsGeneColoc.r` | Extracts matching regional SNPs from the AD GWAS |
| 12 | `12_Harmonisation_genetic_coloc.R` | Harmonises regional data for colocalisation |
| 13 | `13_APOE_data_extraction.R` | Extracts the APOE region from Kunkle et al. (2019), as most APOE variants are absent from Bellenguez et al. (2022) |
| 14 | `14_GeneticColocalisation.R` | Bayesian colocalisation (coloc) |
| 15 | `15_Result_figures.ipynb` | Figures for the causal analysis |

## Discriminative analysis

| Step | File | Description |
|---|---|---|
| 1 | `1_creating_demographic_cohort.ipynb` | Defines AD cases and cognitively normal controls (MMSE ≥ 27 and CDR = 0, following An et al., 2026); excludes participants with other neurological conditions |
| 2 | `2_creating_proteomics_cohort.ipynb` | Extracts MR-identified proteins from the GNPC SomaScan data |
| 3 | `3_creating_combined_cohort.ipynb` | Merges demographic and protein data; keeps baseline visits only |
| 3.1 | `3.1_descriptives.ipynb` | Cohort descriptive statistics |
| 4 | `4_feature_selection.ipynb` | Selects one aptamer per protein and adds age², age × sex terms |
| 5 | `5_model_creation(all_tiers).ipynb` | Models using Tiers 1 + 2 + 3 |
| 5.1 | `5.1_model_creation(tier_1_and_2).ipynb` | Models using Tiers 1 + 2 |
| 5.2 | `5.2_model_creation(tier_1).ipynb` | Models using Tier 1 only |
| 6 | `6_paper_figures.ipynb` | Figures for the discriminative analysis |

## Data availability

| Data | Source | Access |
|---|---|---|
| UKB-PPP protein GWAS (Sun et al., 2023) | [Synapse project syn51364943](https://www.synapse.org/Synapse:syn51364943); European-ancestry folder `syn51365303` | Free, Synapse account required |
| AD GWAS (Bellenguez et al., 2022) | [GWAS Catalog GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158) | Public |
| AD GWAS, APOE region (Kunkle et al., 2019) | [IEU OpenGWAS `ieu-b-2`](https://gwas.mrcieu.ac.uk/datasets/ieu-b-2/) | Free, OpenGWAS token required |
| Olink protein map (`olink_protein_map_3k_v1.tsv`) | UKB-PPP metadata, [Synapse project syn51364943](https://www.synapse.org/Synapse:syn51364943) | Free, Synapse account required |
| 1000 Genomes European LD reference panel | [MRC IEU (`1kg.v3.tgz`)](http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz); uses the `EUR` files | Public |
| GNPC V1 harmonised dataset (Imam et al., 2025) | [AD Data Initiative, AD Workbench](https://www.alzheimersdata.org/gnpc-v1-harmonized-dataset) | Free to qualified researchers, by application |

## Requirements

- Python 3.11.6: pandas, numpy, scipy, statsmodels, scikit-learn (1.8.0),
  scikit-optimize, shap, matplotlib, seaborn, synapseclient, requests
- R 4.5.2: TwoSampleMR (0.7.1), coloc (5.2.3), ieugwasr, data.table, dplyr,
  ggplot2, patchwork, httr, jsonlite, rtracklayer (Bioconductor)
- PLINK 1.9

## Citation

If you use this code, please cite:



## Licence

This code is released under the MIT Licence. The licence covers the code only; each
dataset remains subject to its own terms of use.

## Contact

Olivia Pownall, Imperial College London — olivia.pownall25@imperial.ac.uk
