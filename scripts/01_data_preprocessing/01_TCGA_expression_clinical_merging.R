## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 01 | Data preprocessing (TCGA expression + clinical data, validation cohorts)
## Script  : 01_TCGA_expression_clinical_merging.R
## Purpose : Merge TCGA-HNSC gene expression with curated clinical data for the candidate gene
##           sets (log2FC threshold 0.5 and 1): Ensembl-to-symbol mapping (biomaRt), log2
##           transformation, sample matching and renaming of clinical columns.
## Input   : data/tcga/expData.txt, data/tcga/clinical.txt
## Output  : results/01_data_preprocessing/merged_data_final_threshold_{0_5,1}.txt
## =============================================================================

## ---- Packages ---------------------------------------------------------------

library(scales)
library(readr)
library(survival)
library(caret)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(limma)
library(edgeR)
library(tibble)
library(stringr)
library(Matrix)
library(glmnet)
library(ggpubr)
library(survminer)
library(devtools)
library(AnnotationDbi)
library(MASS)
library(tidyr)
library(cowplot)
library(pheatmap)
library(grid)
library(leaps)
library(GGally)
library(corrplot)
library(broom)
library(biomaRt)

## ---- Paths ------------------------------------------------------------------
## Run from the repository root. Inputs are read from data/tcga/
## (see data/README.md); outputs are written to results/01_data_preprocessing/.
data_dir <- file.path(getwd(), "data", "tcga")
out_dir  <- file.path(getwd(), "results", "01_data_preprocessing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## Sub data
file_name <- file.path(data_dir, "expData.txt")
expData <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

file_name <- file.path(data_dir, "clinical.txt")
clinical <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

# Remove the last character "-01A" from all patient IDs
colnames(expData) <- gsub("-\\d+A$", "", colnames(expData))

# Print the updated column names
print(colnames(expData))

# Extract the Ensembl ID column from expData
ensemblIDs <- expData$Ensembl_ID
# Remove the Ensembl ID column from expData
expData <- expData[, -1]
# Set rownames of expData to be ensemblIDs
rownames(expData) <- ensemblIDs

# Remove transcript ID from row names
rownames(expData) <- sub("\\..*", "", rownames(expData))
rownames(expData)

selected_genes <- c("MAD1L1", "CYP26B1", "ITGA3", "ITGAX", "MMP25", "SERPINA1",
                    "PLAUR", "VCAN", "FAS", "CD44", "TMSB10", "CCL5", "TRIB2",
                    "PLK3", "TBC1D1", "PDE4A", "PKM", "ATP2B1", "MAP3K9", "MAP4K4",
                    "NTRK2", "PTGS2", "CA12", "FSCN1", "RAC2", "MCAM", "NFKB2",
                    "IL4R", "ITGB5", "PLOD1", "CD59", "MMP2", "PTK6", "ICAM1",
                    "LAMB1", "DPYSL2", "PTGS1", "PGF", "HCN2", "LGALS1", "ITGA5",
                    "PROCR", "NFATC2", "SLCO4A1", "ICAM3", "SLC17A9", "JAG1", "SMS",
                    "NPR2", "KCNN4", "SQLE", "PLAT", "GPI", "CAV1", "FGFR1",
                    "ISG15", "DKK1", "RPS6KA2", "LRRC59", "PMP22", "HMGCS1",
                    "PANX1", "ASIC1", "SH2B3", "CD83", "AXL", "VEGFA", "PAX8",
                    "HBEGF", "LIFR", "PIM1", "ALOX15B", "FN1", "ODC1", "SDC1",
                    "ERRFI1", "GRK5", "PDK1", "PRKCA", "NRP2", "TNFAIP3", "CCN2",
                    "PRSS3", "GPR68", "ACAT2", "TNFRSF10B", "S1PR3", "INHBA",
                    "P4HA1", "TUBA1B", "GRM4", "EREG", "IL1B", "CSF2", "MET",
                    "GLIPR1", "PCSK9", "TEP1", "TPBG", "DPP6", "DNMT1", "SCD", "BHLHE40", "LDHA", "FST",
                    "COL4A2", "AREG", "IL1A", "CD63", "DYRK3", "NUDT15", "CTSV",
                    "TCF19", "ARRB1", "MLKL", "ETV1", "SLC22A3", "EDNRA", "OASL",
                    "CDH11", "SLC16A3", "SAE1", "PADI1", "MCL1", "PPP1R18",
                    "COL4A1", "FOSL1", "HTR7", "HYOU1", "ITGB1", "APLN", "ETV5",
                    "ADAM8", "ZFP36L2", "ADAMTS1", "WNT7A", "SLC16A1", "KCNMA1",
                    "NRG1", "HK2", "GNE", "PDE9A", "TGM2", "AK4", "LDLR", "SGMS2",
                    "ANXA5", "F2RL1", "PLAU", "FYN", "E2F7", "IFI27", "PPIB", "CNPY4", "TUBA1A",
                    "PIK3CD", "TK1", "DDIT4", "MGAT2", "FEN1", "VEGFC", "MMP9",
                    "PTAFR", "FOS", "S1PR1", "CHST11", "IL7R", "SLFN11", "DHCR7",
                    "PC", "EGLN3", "ETV4", "TUBB6", "CDK5R1", "JUN", "CALR",
                    "LOXL2", "F2R", "ANXA2", "GJC1", "GPR39", "SMTN", "SLITRK6",
                    "STING1", "KCNQ5", "TEAD1", "NTSR1", "SERPINE1", "GJB4", "TUBB",
                    "PPIA", "FUT4", "S100A2", "LAMB3", "DUSP6", "SERPINB2", "LPAR1",
                    "PNP", "NT5E", "TGFA", "THBS1", "VIM")


ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'),
                   filters='hgnc_symbol',
                   values=selected_genes,
                   mart=ensembl)

# Print the result
gene_info

# Transform the gene expression data into a wide format
ensembl_ids <- gene_info$ensembl_gene_id

# Find the common Ensembl IDs between ensembl_ids and expData
common_ensembl_ids <- intersect(ensembl_ids, rownames(expData))

# Subset expData_se using common_ensembl_ids
expData_common <- expData[common_ensembl_ids,]
# Print the dimensions of the subsetted data
print(dim(expData_common))

## Convert ensembl id to gene symbol
# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Extract the Ensembl IDs from the row names of the transposed matrix
ensembl_ids <- rownames(expData_common)

# Get gene symbols for the Ensembl IDs
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = ensembl_ids,
                   mart = ensembl)

# Create a named vector to map Ensembl IDs to gene symbols
ensembl_to_symbol <- setNames(gene_info$hgnc_symbol, gene_info$ensembl_gene_id)

# Replace the row names (Ensembl IDs) with the gene symbols
rownames(expData_common) <- ensembl_to_symbol[ensembl_ids]

##clean expData
head(expData_common)
any(is.na(expData_common))

expData_common <- apply(expData_common, 2, function(x) gsub("\\,", "", x))

clean_expData_common <- apply(expData_common, 2, function(x) as.numeric(gsub("\\.", "", x)))
clean_expData_common <- clean_expData_common[complete.cases(clean_expData_common), ]
clean_expData_common[is.na(clean_expData_common)] <- 0

clean_expData_common_log2 <- log2(clean_expData_common)

## transpose clean_expData_common_log2
log2_expData_t <- t(clean_expData_common_log2)


## clinical
# Extract the Sample_ID column from clinicalData
sampleIDs <- clinical$Sample_ID
duplicated(sampleIDs)
rownames(clinical) <- sampleIDs

# intersect clinical and filtered_expData_t to have a common sample_ids
common_samples <- intersect(rownames(clinical), rownames(log2_expData_t))

## subset clinical and log2_expData_t based on common samples
clinical_common <- clinical[common_samples, ]
expData_common <- log2_expData_t[rownames(clinical_common), ]

## merge the filtered_expData_t_common and clinical_common
merged_data <- cbind(expData_common, clinical_common)
names(merged_data)

## assign the colnmaes of the merged data
colnames(merged_data) <- ensembl_to_symbol[ensembl_ids]
## rename last 9 column of the merged_data
names(merged_data)[(ncol(merged_data)-8):ncol(merged_data)] <- c("Age",
                                                                 "Aneuploidy_score", "Neoplasm_histologic_grade",
                                                                 "Disease_specific_survival_Months",
                                                                 "OS_Months", "Progress_free_survival_Months", "Vital_status",
                                                                 "Radiation_therapy", "HPV_status")

##there is column named called 'NA' which contains Patients Id, so lets remove it
merged_data <- subset(merged_data, select = -`NA`)

write.table(merged_data, file = "merged_data_final_threshold_0_5.txt", sep = "\t", quote = FALSE, row.names = TRUE)


## For throshold 1
## Sub data
file_name <- file.path(data_dir, "expData.txt")
expData <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

file_name <- file.path(data_dir, "clinical.txt")
clinical <- read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

# Remove the last character "-01A" from all patient IDs
colnames(expData) <- gsub("-\\d+A$", "", colnames(expData))

# Print the updated column names
print(colnames(expData))

# Extract the Ensembl ID column from expData
ensemblIDs <- expData$Ensembl_ID
# Remove the Ensembl ID column from expData
expData <- expData[, -1]
# Set rownames of expData to be ensemblIDs
rownames(expData) <- ensemblIDs

# Remove transcript ID from row names
rownames(expData) <- sub("\\..*", "", rownames(expData))
rownames(expData)

selected_genes <- c("MAD1L1", "CYP26B1", "ITGA3", "ITGAX", "SERPINA1", "PLAUR",
                    "VCAN", "FAS", "TMSB10", "CCL5", "TRIB2", "PLK3", "PKM",
                    "ATP2B1", "MAP3K9", "MAP4K4", "NTRK2", "PTGS2", "CA12",
                    "FSCN1", "RAC2", "MCAM", "IL4R", "PTK6", "ICAM1", "DPYSL2",
                    "PGF", "LGALS1", "ITGA5", "PROCR", "NFATC2", "SLCO4A1", "ICAM3",
                    "SLC17A9", "KCNN4", "SQLE", "PLAT", "ISG15", "HMGCS1", "ASIC1",
                    "VEGFA", "PAX8", "HBEGF", "ALOX15B", "ODC1", "CCN2", "PRSS3",
                    "GPR68", "ACAT2", "S1PR3", "INHBA", "P4HA1", "TUBA1B", "GRM4", "EREG",
                    "IL1B", "CSF2", "GLIPR1", "PCSK9", "TPBG", "SCD", "BHLHE40",
                    "FST", "COL4A2", "AREG", "IL1A", "TCF19", "ETV1", "SLC22A3",
                    "EDNRA", "OASL", "CDH11", "SLC16A3", "PADI1", "COL4A1", "FOSL1",
                    "APLN", "ETV5", "ADAM8", "WNT7A", "KCNMA1", "NRG1", "HK2", "TGM2",
                    "AK4", "LDLR", "PLAU", "E2F7", "IFI27", "TK1", "VEGFC", "MMP9",
                    "FOS", "S1PR1", "CHST11", "IL7R", "DHCR7", "PC", "EGLN3", "ETV4",
                    "LOXL2", "SMTN", "NTSR1", "SERPINE1", "S100A2",
                    "LAMB3", "DUSP6", "NT5E", "TGFA", "THBS1", "VIM", "GALNT14",
                    "CTSS", "HLA-DQB1", "COL6A3", "CYP51A1", "TNFRSF12A", "TBXA2R",
                    "TBXAS1", "CACNB1", "SREBF1", "IL12RB2", "IL11", "CHEK1",
                    "CYP27B1", "MT2A", "OPRL1", "TRPC6", "SLC5A2", "SYT2", "EPOR",
                    "H2AX", "RCSD1", "GJC2", "TNFRSF1B", "RAD51", "PKMYT1", "AURKB",
                    "BIRC5", "E2F1", "NUDT1", "EGR2", "CYP2J2", "CDC25A", "CXCL8",
                    "FOXC2", "TYMS", "FUT2", "DNER", "SH2D5", "S100A4", "BLM",
                    "S1PR2")

ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id'),
                   filters='hgnc_symbol',
                   values=selected_genes,
                   mart=ensembl)

# Print the result
gene_info

# Transform the gene expression data into a wide format
ensembl_ids <- gene_info$ensembl_gene_id

# Find the common Ensembl IDs between ensembl_ids and expData
common_ensembl_ids <- intersect(ensembl_ids, rownames(expData))

# Subset expData_se using common_ensembl_ids
expData_common <- expData[common_ensembl_ids,]
# Print the dimensions of the subsetted data
print(dim(expData_common))

## Convert ensembl id to gene symbol
# Connect to the Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Extract the Ensembl IDs from the row names of the transposed matrix
ensembl_ids <- rownames(expData_common)

# Get gene symbols for the Ensembl IDs
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                   filters = 'ensembl_gene_id',
                   values = ensembl_ids,
                   mart = ensembl)

# Create a named vector to map Ensembl IDs to gene symbols
ensembl_to_symbol <- setNames(gene_info$hgnc_symbol, gene_info$ensembl_gene_id)

# Replace the row names (Ensembl IDs) with the gene symbols
rownames(expData_common) <- ensembl_to_symbol[ensembl_ids]

##clean expData
head(expData_common)
any(is.na(expData_common))

expData_common <- apply(expData_common, 2, function(x) gsub("\\,", "", x))

clean_expData_common <- apply(expData_common, 2, function(x) as.numeric(gsub("\\.", "", x)))
clean_expData_common <- clean_expData_common[complete.cases(clean_expData_common), ]
clean_expData_common[is.na(clean_expData_common)] <- 0

clean_expData_common_log2 <- log2(clean_expData_common)

## transpose clean_expData_common_log2
log2_expData_t <- t(clean_expData_common_log2)


## clinical
# Extract the Sample_ID column from clinicalData
sampleIDs <- clinical$Sample_ID
duplicated(sampleIDs)
rownames(clinical) <- sampleIDs

# intersect clinical and filtered_expData_t to have a common sample_ids
common_samples <- intersect(rownames(clinical), rownames(log2_expData_t))

## subset clinical and log2_expData_t based on common samples
clinical_common <- clinical[common_samples, ]
expData_common <- log2_expData_t[rownames(clinical_common), ]

## merge the filtered_expData_t_common and clinical_common
merged_data <- cbind(expData_common, clinical_common)
names(merged_data)

## assign the colnmaes of the merged data
colnames(merged_data) <- ensembl_to_symbol[ensembl_ids]
## rename last 5 column of the merged_data
names(merged_data)[(ncol(merged_data)-7):ncol(merged_data)] <- c("Age",
                                                                 "Aneuploidy_score", "Neoplasm_histologic_grade",
                                                                 "OS_Months", "Vital_status", "Progress_free_survival_Months",
                                                                 "Radiation_therapy", "HPV_status")

##there is column named called 'NA' which contains Patients Id, so lets remove it
merged_data <- subset(merged_data, select = -`NA`)

write.table(merged_data, file = "merged_data_final_threshold_1.txt", sep = "\t", quote = FALSE, row.names = TRUE)
