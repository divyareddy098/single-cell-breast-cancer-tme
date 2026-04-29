# Seurat Analysis Pipeline

library(zellkonverter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)

# Load data
sce <- readH5AD("data/raw/breast_cancer_tme.h5ad")

# Convert to Seurat
seurat_obj <- as.Seurat(sce, counts = "X")


# QC summary
qc_summary <- data.frame(
  Metric = c("Total cells", "Total genes"),
  Value = c(ncol(seurat_obj), nrow(seurat_obj))
)

write.csv(qc_summary, "results/seurat_qc_summary.csv", row.names = FALSE)

# Normalize & PCA
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)


# Clustering & UMAP
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20, seed.use = 42)


# Plots

p1 <- DimPlot(seurat_obj, group.by = "celltype_major", label = TRUE) + NoLegend()
ggsave("figures/seurat_umap_celltype_major.png", p1, width = 8, height = 6, dpi = 300)

p2 <- DimPlot(seurat_obj, group.by = "subtype")
ggsave("figures/seurat_umap_subtype.png", p2, width = 8, height = 6, dpi = 300)

# Cell type counts
celltype_counts <- as.data.frame(table(seurat_obj$celltype_major))
colnames(celltype_counts) <- c("Cell_Type", "Cell_Count")
write.csv(celltype_counts, "results/seurat_celltype_counts.csv", row.names = FALSE)


# Subtype counts
subtype_counts <- as.data.frame(table(seurat_obj$subtype))
colnames(subtype_counts) <- c("Subtype", "Cell_Count")
write.csv(subtype_counts, "results/seurat_subtype_counts.csv", row.names = FALSE)

# Marker genes
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
write.csv(markers, "results/seurat_markers.csv")


# Save object
saveRDS(seurat_obj, "results/seurat_object.rds")

cat("Seurat pipeline complete\n")
