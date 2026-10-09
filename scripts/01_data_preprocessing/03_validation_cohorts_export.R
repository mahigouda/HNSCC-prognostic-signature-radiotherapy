## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 01 | Data preprocessing (TCGA expression + clinical data, validation cohorts)
## Script  : 03_validation_cohorts_export.R
## Purpose : Export expression and phenotype tables of the two independent oral cavity SCC
##           cohorts (MDACC: GSE42743, FHCRC: GSE41613) from R data objects to text files.
## Input   : data/validation_cohorts/*.Rdata
## Output  : results/01_data_preprocessing/{exp,pd}_{mdacc,fhcrc}.txt
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
## Run from the repository root. Inputs are read from data/validation_cohorts/
## (see data/README.md); outputs are written to results/01_data_preprocessing/.
data_dir <- file.path(getwd(), "data", "validation_cohorts")
out_dir  <- file.path(getwd(), "results", "01_data_preprocessing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

load(file.path(data_dir, "GSE42743&41613.Rdata"))
ls()

## Export files as txt
write.table(exp_mdacc, file = "exp_mdacc.txt", sep = "\t", row.names = TRUE)
write.table(pd_mdacc, file = "pd_mdacc.txt", sep = "\t", row.names = TRUE)

write.table(exp_fhcrc, file = "exp_fhcrc.txt", sep = "\t", row.names = TRUE)
write.table(pd_fhcrc, file = "pd_fhcrc.txt", sep = "\t", row.names = TRUE)

load(file.path(data_dir, "Mdacc_cohort.Rdata"))
load(file.path(data_dir, "Fhcrc_cohort.Rdata"))
