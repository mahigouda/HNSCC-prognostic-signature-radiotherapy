## Install the R / Bioconductor packages used in this repository.
cran <- c("survival", "survminer", "forestplot", "forestmodel", "car", "lmtest", "caret",
          "glmnet", "leaps", "MASS", "broom", "plyr", "dplyr", "tidyr", "tibble", "purrr",
          "readr", "stringr", "scales", "ggplot2", "ggpubr", "ggfortify", "cowplot", "GGally",
          "corrplot", "pheatmap", "gplots", "gridExtra", "RColorBrewer", "circlize", "factoextra",
          "metafor", "gt", "knitr", "Seurat", "patchwork", "reshape2", "zoo")
bioc <- c("biomaRt", "AnnotationDbi", "SummarizedExperiment", "limma", "edgeR", "Biobase",
          "progeny", "GSVA", "ComplexHeatmap", "clusterProfiler", "enrichplot", "scater")
install.packages(setdiff(cran, rownames(installed.packages())))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(setdiff(bioc, rownames(installed.packages())), update = FALSE)
