import scanpy as sc
import pandas as pd
from pathlib import Path

RAW = "data/raw/breast_cancer_tme.h5ad"

out_results = Path("results/markers")
out_results.mkdir(parents=True, exist_ok=True)

print("Loading dataset...")
adata = sc.read_h5ad(RAW)

# Use feature names if available
if "feature_name" in adata.var.columns:
    adata.var_names = adata.var["feature_name"].astype(str)
    adata.var_names_make_unique()

# Keep cells with major cell type labels
adata = adata[adata.obs["celltype_major"].notna()].copy()

print("Normalizing data for marker analysis...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print("Running marker gene analysis by major cell type...")
sc.tl.rank_genes_groups(
    adata,
    groupby="celltype_major",
    method="wilcoxon",
    n_genes=50
)

markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(out_results / "marker_genes_by_major_celltype.csv", index=False)

top10 = (
    markers
    .groupby("group", group_keys=False)
    .head(10)
)

top10.to_csv(out_results / "top10_marker_genes_by_major_celltype.csv", index=False)

print("Saved:")
print(" - results/markers/marker_genes_by_major_celltype.csv")
print(" - results/markers/top10_marker_genes_by_major_celltype.csv")

print("\nTop markers:")
print(top10[["group", "names", "scores", "pvals_adj", "logfoldchanges"]].head(30))
