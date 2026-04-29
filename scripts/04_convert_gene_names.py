import pandas as pd
import scanpy as sc
from pathlib import Path

RAW = "data/raw/breast_cancer_tme.h5ad"
MARKERS = "results/markers/top10_marker_genes_by_major_celltype.csv"

out_dir = Path("results/markers")
out_dir.mkdir(parents=True, exist_ok=True)

print("Loading dataset for gene mapping...")
adata = sc.read_h5ad(RAW)

# Create mapping ENSG -> gene symbol
if "feature_name" in adata.var.columns:
    mapping = dict(zip(adata.var_names, adata.var["feature_name"]))
else:
    raise ValueError("feature_name column not found in dataset")

print("Loading marker genes...")
df = pd.read_csv(MARKERS)

# Map gene names
df["gene_symbol"] = df["names"].map(mapping)

# Replace missing with original ENSG
df["gene_symbol"] = df["gene_symbol"].fillna(df["names"])

df.to_csv(out_dir / "top10_marker_genes_with_symbols.csv", index=False)

print("Saved:")
print(" - results/markers/top10_marker_genes_with_symbols.csv")

print("\nPreview:")
print(df.head(20))
