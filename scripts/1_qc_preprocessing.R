# File: scripts/01_qc_preprocessing.R
# Project: Single-Cell RNA-seq Analysis of Breast Cancer Tumor Microenvironment
# Purpose: Load data, compute QC metrics, filter cells, normalize, and save processed object

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

set.seed(123)

# inputs
input_rds  <- "data/breast_cancer_raw.rds"   # replace with your actual file
output_rds <- "results/seurat_qc_processed.rds"

dir.create("results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load object
if (!file.exists(input_rds)) {
  stop("Input file not found: ", input_rds,
       "\nPut your Seurat object at this path or update input_rds.")
}

seurat_obj <- readRDS(input_rds)

# QC metrics
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

png("results/figures/qc_violin_plots.png", width = 1200, height = 500)
print(VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
))
dev.off()

png("results/figures/qc_feature_scatter.png", width = 1000, height = 500)
print(
  CombinePlots(plots = list(
    FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt"),
    FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ))
)
dev.off()

# Filter cells
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 15
)

# Normalize and variable features
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))

# Save 
saveRDS(seurat_obj, output_rds)

message("Saved processed Seurat object to: ", output_rds)
