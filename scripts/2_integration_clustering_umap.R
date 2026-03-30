# File: scripts/2_integration_clustering_umap.R
# Purpose: Run integration-ready preprocessing, PCA, clustering, and UMAP

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

set.seed(123)

# User inputs
input_rds  <- "results/seurat_qc_processed.rds"
output_rds <- "results/seurat_clustered_umap.rds"

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load processed object
if (!file.exists(input_rds)) {
  stop("Processed object not found: ", input_rds,
       "\nRun scripts/01_qc_preprocessing.R first.")
}

seurat_obj <- readRDS(input_rds)

# Optional batch correction / integration note

if ("batch" %in% colnames(seurat_obj@meta.data)) {
  message("Batch column detected. You can extend this script with Seurat integration if needed.")
}

# Dimensionality reduction
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

png("results/figures/elbow_plot.png", width = 800, height = 500)
print(ElbowPlot(seurat_obj, ndims = 30))
dev.off()

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

png("results/figures/umap_clusters.png", width = 900, height = 700)
print(
  DimPlot(seurat_obj, reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle("UMAP of Breast Cancer scRNA-seq Clusters")
)
dev.off()

saveRDS(seurat_obj, output_rds)

message("Saved clustered Seurat object to: ", output_rds)
