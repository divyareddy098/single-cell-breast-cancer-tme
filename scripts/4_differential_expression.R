# Purpose: Run differential expression between two chosen groups

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})


# User inputs
input_rds <- "results/seurat_annotated.rds"

group_column <- "cell_type"   # can also be "seurat_clusters"
ident_1 <- "Cluster_0"
ident_2 <- "Cluster_1"

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)


# Load object
if (!file.exists(input_rds)) {
  stop("Annotated object not found: ", input_rds,
       "\nRun scripts/03_celltype_annotation_markers.R first.")
}

seurat_obj <- readRDS(input_rds)

if (!(group_column %in% colnames(seurat_obj@meta.data))) {
  stop("Column not found in metadata: ", group_column)
}

Idents(seurat_obj) <- seurat_obj[[group_column, drop = TRUE]]

de_res <- FindMarkers(
  seurat_obj,
  ident.1 = ident_1,
  ident.2 = ident_2,
  logfc.threshold = 0.25,
  min.pct = 0.1
)

out_file <- paste0(
  "results/tables/differential_expression_",
  gsub(" ", "_", ident_1), "_vs_", gsub(" ", "_", ident_2), ".csv"
)

write.csv(de_res, out_file, row.names = TRUE)

message("DE results saved to: ", out_file)
