import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

RAW = "data/raw/breast_cancer_tme.h5ad"

out_data = Path("data/processed")
out_results = Path("results")
out_figures = Path("figures")

out_data.mkdir(parents=True, exist_ok=True)
out_results.mkdir(parents=True, exist_ok=True)
out_figures.mkdir(parents=True, exist_ok=True)

print("Loading dataset...")
adata = sc.read_h5ad(RAW)

print(adata)

# Basic QC summary
qc_summary = pd.DataFrame({
    "Metric": [
        "Total cells",
        "Total genes",
        "Median genes per cell",
        "Median UMIs per cell",
        "Median mitochondrial percent",
        "Number of donors",
        "Number of major cell types"
    ],
    "Value": [
        adata.n_obs,
        adata.n_vars,
        round(adata.obs["nFeature_RNA"].median(), 2),
        round(adata.obs["nCount_RNA"].median(), 2),
        round(adata.obs["percent_mito"].median(), 2),
        adata.obs["donor_id"].nunique(),
        adata.obs["celltype_major"].nunique()
    ]
})
qc_summary.to_csv(out_results / "qc_summary.csv", index=False)

print("\nQC summary:")
print(qc_summary)

# Cell type count summary
celltype_counts = (
    adata.obs["celltype_major"]
    .value_counts()
    .reset_index()
)
celltype_counts.columns = ["Cell_Type", "Cell_Count"]
celltype_counts.to_csv(out_results / "celltype_major_counts.csv", index=False)

# Minor cell type count summary
minor_counts = (
    adata.obs["celltype_minor"]
    .value_counts()
    .reset_index()
)
minor_counts.columns = ["Cell_Type_Minor", "Cell_Count"]
minor_counts.to_csv(out_results / "celltype_minor_counts.csv", index=False)

# Donor/subtype summary
donor_summary = (
    adata.obs[["donor_id", "subtype"]]
    .drop_duplicates()
    .sort_values("donor_id")
)
donor_summary.to_csv(out_results / "donor_subtype_summary.csv", index=False)

# Cell type composition by cancer subtype
composition = (
    adata.obs
    .groupby(["subtype", "celltype_major"], observed=True)
    .size()
    .reset_index(name="Cell_Count")
)

composition["Fraction"] = composition.groupby("subtype", observed=True)["Cell_Count"].transform(
    lambda x: x / x.sum()
)

composition.to_csv(out_results / "celltype_composition_by_subtype.csv", index=False)

# Cell type composition by donor
donor_composition = (
    adata.obs
    .groupby(["donor_id", "celltype_major"], observed=True)
    .size()
    .reset_index(name="Cell_Count")
)

donor_composition["Fraction"] = donor_composition.groupby("donor_id", observed=True)["Cell_Count"].transform(
    lambda x: x / x.sum()
)

donor_composition.to_csv(out_results / "celltype_composition_by_donor.csv", index=False)

# Use existing UMAP if available
if "X_umap" in adata.obsm:
    print("Using existing UMAP coordinates...")
    sc.pl.umap(
        adata,
        color=["celltype_major"],
        frameon=False,
        show=False,
        title="Breast Cancer TME Cell Types"
    )
    plt.savefig(out_figures / "umap_celltype_major.png", dpi=300, bbox_inches="tight")
    plt.close()

    sc.pl.umap(
        adata,
        color=["subtype"],
        frameon=False,
        show=False,
        title="Breast Cancer Subtypes"
    )
    plt.savefig(out_figures / "umap_subtype.png", dpi=300, bbox_inches="tight")
    plt.close()

# Cell type count barplot
plt.figure(figsize=(9, 5))
celltype_counts.plot(
    x="Cell_Type",
    y="Cell_Count",
    kind="bar",
    legend=False
)
plt.title("Major Cell Type Composition")
plt.xlabel("Major Cell Type")
plt.ylabel("Number of Cells")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(out_figures / "major_celltype_counts.png", dpi=300)
plt.close()

# TME composition stacked barplot
pivot = composition.pivot(
    index="subtype",
    columns="celltype_major",
    values="Fraction"
).fillna(0)

pivot.plot(
    kind="bar",
    stacked=True,
    figsize=(10, 6)
)
plt.title("Tumor Microenvironment Composition by Breast Cancer Subtype")
plt.xlabel("Breast Cancer Subtype")
plt.ylabel("Cell Type Fraction")
plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.savefig(out_figures / "tme_composition_by_subtype.png", dpi=300)
plt.close()

# Save lightweight processed object with metadata + embeddings
adata.write_h5ad(out_data / "breast_cancer_tme_processed_metadata.h5ad", compression="gzip")

print("\nSaved results:")
print(" - results/qc_summary.csv")
print(" - results/celltype_major_counts.csv")
print(" - results/celltype_minor_counts.csv")
print(" - results/donor_subtype_summary.csv")
print(" - results/celltype_composition_by_subtype.csv")
print(" - results/celltype_composition_by_donor.csv")

print("\nSaved figures:")
print(" - figures/umap_celltype_major.png")
print(" - figures/umap_subtype.png")
print(" - figures/major_celltype_counts.png")
print(" - figures/tme_composition_by_subtype.png")

print("\nDone.")
