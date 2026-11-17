## Calling libraries
library(data.table) # needed for fread -> helps handle large data
library(ieugwasr) # provides the clumping function
library(rtracklayer) # provides a chain file that can translate GRCh38 (original SNP dataset) to GRCh37 (1000 Genomes EUR reference panel -> correlation dataset)
library(httr)
library(jsonlite)

## Loading in dataset
snps <- fread("significant_snps.csv") # This is a more efficient way to read in such a big file
head(snps)

# Because the kernel kept restarting in the extraction, it is not implausible that we have ended up with some duplicates
# Checking for duplicates: Completely identical rows
identical_dupes <- nrow(snps) - nrow(distinct(snps))
print(paste("Completely identical duplicate rows:", identical_dupes))


## Assessing how many proteins we've gotten from the initial filter
unique_proteins <- unique(snps$source_protein)
length(unique_proteins)
#2620 proteins that have at least one SNP that meets the significance and strength threshold

# Looking at the number of SNPs for each protein
snps[, .N, by = source_protein][order(-N)]
# Big variation in how many SNPs each protein has

## FILTERING TO CIS-pQTLS WITHIN 1MB OF ENCODING GENE

## Read in annotation file
protein_info <- fread("olink_protein_map_3k_v1.tsv") # File contains location information about where each protein-encoding gene is

## Check for duplicates in the annotation file
protein_info[duplicated(protein_info$HGNC.symbol), HGNC.symbol]
# Some proteins are measured in multiple Olink panels (e.g. neurology and immunology) so need to ensure only one remains

## Keep one entry per gene symbol
protein_info_unique <- unique(protein_info[, .(HGNC.symbol, chr, gene_start, gene_end)]) # only keeping relevant location-related columns

## Merges gene coordinates onto the SNP data by matching the protein name
snps <- merge(snps, protein_info_unique,
              by.x = "source_protein",
              by.y = "HGNC.symbol")

## Cis filter - SNPs within 1Mb of encoding gene
snps_cis <- snps[
  CHROM == chr & # Must be on the same chromosone as the encoding gene
    GENPOS >= (gene_start - 1e6) & # position must be not more than 1Mb before the gene start
    GENPOS <= (gene_end + 1e6) # position must be not more than 1Mb after the gene end
]

# Checks remaining SNPs and proteins
cat("SNPs remaining after cis filter:", nrow(snps_cis), "\n")
cat("Proteins remaining after cis filter:", uniqueN(snps_cis$source_protein), "\n")

# Write to file
fwrite(snps_cis, "snps_cis")
