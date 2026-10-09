## =============================================================================
## Project : Prognostic gene-expression signature for radiotherapy outcome in HPV-negative HNSCC
## Module  : 05 | Single-cell expression and pathway activity of signature genes (GSE103322)
## Script  : 01_signature_genes_single_cell_PROGENy.R
## Purpose : Seurat processing of malignant HNSCC single cells, expression of OS/PFS signature
##           genes (feature/dot plots) and correlation of signature genes with PROGENy
##           pathway activity scores (heatmaps with significance annotation).
## Input   : data/single_cell/data_malignant.txt (GSE103322, malignant cells)
## Output  : plots (interactive)
## Notes   : CTSV is not detected in this dataset (alias CTSL2).
## =============================================================================

## ---- Packages ---------------------------------------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(scater)
library(monocle)
library(reshape2)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(GSVA)
library(progeny)
library(IRanges)
library(S4Vectors)
library(gplots)
library(ComplexHeatmap)
library(grid)
library(circlize)
library(Biobase)
library(BiocGenerics)
library(gridExtra)
library(factoextra)

## ---- Paths ------------------------------------------------------------------
## Run from the repository root. Inputs are read from data/single_cell/
## (see data/README.md); outputs are written to results/05_single_cell/.
data_dir <- file.path(getwd(), "data", "single_cell")
out_dir  <- file.path(getwd(), "results", "05_single_cell")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)


# Load the data
data <- read.table(file.path(data_dir, "data_malignant.txt"), header = TRUE, row.names = 1, sep = "\t")
sce <- CreateSeuratObject(counts = data)

# Quality Control
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 9000)

# Preprocess the data
sce <- SCTransform(sce)
sce <- FindVariableFeatures(sce)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
# Run UMAP
sce <- RunUMAP(sce, dims = 1:10)
# Run t-SNE
sce <- RunTSNE(sce, dims = 1:10)

DimPlot(sce, reduction = "tsne", label = TRUE)

DimPlot(sce, reduction = "umap", label = TRUE)

# Define vector of samples to keep
samples_to_keep <- c("HNSCC16", "HNSCC17", "HNSCC18", "HNSCC20", "HNSCC22", "HNSCC25", "HNSCC26", "HNSCC28", "HNSCC5", "HNSCC6")

# Subset Seurat object to only include cells from these samples
sce_subset <- subset(sce, subset = orig.ident %in% samples_to_keep)

# Plot UMAP with labels
DimPlot(sce_subset, reduction = "umap", label = TRUE)

# Plot tsne with labels
DimPlot(sce_subset, reduction = "tsne", label = TRUE, pt.size = 2.5, repel = TRUE)

# Define the gene list
genes_of_interest <- c('ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB', 'PLAU',
                       'LRRC59', 'CTSV', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL', 'ADAM8', 'RAC2',
                       'CSF2', 'PPP1R18', 'HBEGF', 'ITGA5', 'LRRC59', 'PLAU', 'FOSL1', 'PLAUR')

# Remove duplicates
genes_of_interest <- unique(genes_of_interest)


genes_of_interest_1 <- c('ODC1', 'TUBB', 'ETV4', 'NUDT15')

FeaturePlot(sce, features = genes_of_interest_1, reduction = "tsne", pt.size = 1)

genes_of_interest_2 <- c('AK4', 'PPIB', 'PLAU', 'LRRC59')

FeaturePlot(sce, features = genes_of_interest_2, reduction = "tsne", pt.size = 1)

genes_of_interest_3 <- c('CTSV', 'HBEGF', 'PPIA', 'MAP4K4')

FeaturePlot(sce, features = genes_of_interest_3, reduction = "tsne", pt.size = 1)

genes_of_interest_4 <- c('AXL', 'ADAM8', 'RAC2', 'CSF2')

FeaturePlot(sce, features = genes_of_interest_4, reduction = "tsne", pt.size = 1)

genes_of_interest_5 <- c('PPP1R18', 'ITGA5', 'FOSL1', 'PLAUR')

FeaturePlot(sce, features = genes_of_interest_5, reduction = "tsne", pt.size = 1)

## Signature genes Dot plot
genes_of_interest_dot_plot <- c('ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB', 'PLAU',
                       'LRRC59', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL', 'ADAM8', 'RAC2',
                       'CSF2', 'PPP1R18', 'HBEGF', 'ITGA5', 'LRRC59', 'PLAU', 'FOSL1', 'PLAUR')
# Remove duplicates
genes_of_interest_dot_plot <- unique(genes_of_interest_dot_plot)
# Create the dot plot
dot_plot <- DotPlot(sce_subset, features = genes_of_interest_dot_plot) + RotatedAxis()

# Change the color of the plot
dot_plot <- dot_plot + scale_color_gradient(low = "blue", high = "red")

# Display the plot
print(dot_plot)


## PROGENY
# Extract the expression matrix
# Convert the expression data to a regular matrix
expr_matrix <- as.matrix(sce@assays$RNA@counts)

# Apply PROGENy
pathway_scores <- progeny(expr_matrix, scale=TRUE, organism="Human")

# Define the gene list
genes_of_interest <- c('ODC1', 'TUBB', 'ETV4', 'NUDT15', 'AK4', 'PPIB', 'PLAU',
                       'LRRC59', 'HBEGF', 'PPIA', 'MAP4K4', 'AXL', 'ADAM8', 'RAC2',
                       'CSF2', 'PPP1R18', 'HBEGF', 'ITGA5', 'LRRC59', 'PLAU', 'FOSL1', 'PLAUR')

missing_genes <- genes_of_interest[!(genes_of_interest %in% rownames(sce@assays$RNA@counts))]
print(missing_genes)

# Correlate the expression with the PROGENy scores
correlation_results <- sapply(genes_of_interest, function(gene) {
  gene_expression <- sce@assays$RNA@counts[gene, ]
  sapply(colnames(pathway_scores), function(pathway) {
    pathway_activity <- pathway_scores[, pathway]
    cor.test(gene_expression, pathway_activity)$estimate
  })
})

# View the correlation results
correlation_results

# Transpose the correlation matrix, as pheatmap expects genes (or other samples) to be rows
correlation_results_t <- t(correlation_results)

# Make a heatmap
pheatmap(correlation_results_t,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE)

## Progeny pathways dotplot
# Calculate average expression per pathway
avg_expression <- sapply(genes_of_interest, function(gene) {
  if (gene %in% rownames(sce@assays$RNA@counts)) {
    gene_expression <- sce@assays$RNA@counts[gene, ]
    sapply(colnames(pathway_scores), function(pathway) {
      pathway_cells <- which(pathway_scores[, pathway] > 0)
      mean(gene_expression[pathway_cells])
    })
  } else {
    return(rep(NA, length(colnames(pathway_scores))))
  }
})

# Create a data frame in the format required by DotPlot
dotplot_data <- data.frame(
  gene = rep(rownames(avg_expression), each = ncol(avg_expression)),
  pathway = rep(colnames(avg_expression), nrow(avg_expression)),
  expression = as.vector(avg_expression)
)

# Remove any rows with missing data
dotplot_data <- dotplot_data[!is.na(dotplot_data$expression), ]

## Dot plot using ggplot2
# Convert the average expression matrix to a data frame for ggplot
avg_expression_df <- as.data.frame(avg_expression)
avg_expression_df$gene <- rownames(avg_expression_df)

# Convert to long format
avg_expression_long <- tidyr::pivot_longer(avg_expression_df, -gene, names_to = "pathway", values_to = "expression")

# Remove NA values
avg_expression_long <- avg_expression_long[!is.na(avg_expression_long$expression), ]

# Create the plot
ggplot(avg_expression_long, aes(x = pathway, y = gene, size = expression, color = expression)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_color_gradient(low = "blue", high = "red")


## Correlation with pvalue
# Correlate the expression with the PROGENy scores
correlation_results <- sapply(genes_of_interest, function(gene) {
  if (gene %in% rownames(sce@assays$RNA@counts)) {
    gene_expression <- sce@assays$RNA@counts[gene, ]
    sapply(colnames(pathway_scores), function(pathway) {
      pathway_activity <- pathway_scores[, pathway]
      cor.test(gene_expression, pathway_activity)$estimate
    })
  } else {
    rep(NA, length(colnames(pathway_scores)))
  }
})

print(range(correlation_results, na.rm = TRUE))

# Calculate the p-values for the correlations
pvalues <- sapply(genes_of_interest, function(gene) {
  if (gene %in% rownames(sce@assays$RNA@counts)) {
    gene_expression <- sce@assays$RNA@counts[gene, ]
    sapply(colnames(pathway_scores), function(pathway) {
      pathway_activity <- pathway_scores[, pathway]
      cor.test(gene_expression, pathway_activity)$p.value
    })
  } else {
    rep(NA, length(colnames(pathway_scores)))
  }
})


# Generate the correlation heatmap
correlation_heatmap <- pheatmap(correlation_results,
                                cluster_rows = FALSE,
                                cluster_cols = FALSE,
                                display_numbers = TRUE,
                                number_format = "%.2f")

# Prepare the heatmap annotation based on p-value
significance_annotations <- ifelse(pvalues < 0.001, "***",
                                   ifelse(pvalues < 0.01, "**",
                                          ifelse(pvalues < 0.05, "*", "")))
# Remove .cor from pathway labels
rownames(correlation_results) <- sub("\\.cor", "", rownames(correlation_results))
# Create the heatmap again
heatmap_plot <- Heatmap(correlation_results, name = "correlation",
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_annotations[i, j], x, y, gp = gpar(fontsize = 10))
                        },
                        heatmap_legend_param = list(direction = "horizontal", at = c(-1, 0, 1),
                                                    labels = c("-1", "0", "1"), title = "Correlation"))
# Draw the heatmap
draw(heatmap_plot, heatmap_legend_side = "right")
# Create the legend for significance levels
grid.text("Significance levels: *** p<0.001, ** p<0.01, * p<0.05", x = unit(0.98, "npc"), y = unit(0.02, "npc"),
          just = "right", gp = gpar(fontsize = 10))
