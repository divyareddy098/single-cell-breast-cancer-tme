# Purpose: Identify cluster markers, assign cell types, and save annotated object

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(123)

# User inputs
input_rds  <- "results/seurat_clustered_umap.rds"
output_rds <- "results/seurat_annotated.rds"

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load object
if (!file.exists(input_rds)) {
  stop("Clustered object not found: ", input_rds,
       "\nRun scripts/02_integration_clustering_umap.R first.")
}

seurat_obj <- readRDS(input_rds)

# Marker detection
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(markers, "results/tables/cluster_markers.csv", row.names = FALSE)

# cluster-based cell-type labels
seurat_obj$cell_type <- paste0("Cluster_", Idents(seurat_obj))

# Plot UMAP
png("results/figures/umap_celltypes.png", width = 1000, height = 700)
print(
  DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type", label = TRUE, repel = TRUE) +
    ggtitle("UMAP (Cluster-based Annotation - Placeholder)")
)
dev.off()

# Save 
saveRDS(seurat_obj, output_rds)

message("Saved annotated (placeholder) Seurat object to: ", output_rds)
message("Markers saved to: results/tables/cluster_markers.csv")
