import pandas as pd
from pathlib import Path

out_dir = Path("results/final")
out_dir.mkdir(parents=True, exist_ok=True)

qc = pd.read_csv("results/qc_summary.csv")
major = pd.read_csv("results/celltype_major_counts.csv")
composition = pd.read_csv("results/celltype_composition_by_subtype.csv")
markers = pd.read_csv("results/markers/top10_marker_genes_with_symbols.csv")

# Overall major cell type fractions
total_cells = major["Cell_Count"].sum()
major["Fraction"] = major["Cell_Count"] / total_cells
major.to_csv(out_dir / "major_celltype_fraction_summary.csv", index=False)

# Dominant cell type per subtype
dominant = (
    composition
    .sort_values(["subtype", "Fraction"], ascending=[True, False])
    .groupby("subtype", as_index=False)
    .head(3)
)
dominant.to_csv(out_dir / "top_celltypes_by_subtype.csv", index=False)

# Top marker symbols per major cell type
marker_summary = (
    markers
    .groupby("group")["gene_symbol"]
    .apply(lambda x: "; ".join(x.head(10)))
    .reset_index()
)
marker_summary.columns = ["Cell_Type", "Top_Marker_Genes"]
marker_summary.to_csv(out_dir / "top_marker_summary_by_celltype.csv", index=False)

# Simple immune/stromal/epithelial grouping
immune_types = ["T-cells", "B-cells", "Myeloid", "Plasmablasts"]
stromal_types = ["CAFs", "Endothelial", "PVL"]
epithelial_types = ["Cancer Epithelial", "Normal Epithelial"]

def category(celltype):
    if celltype in immune_types:
        return "Immune"
    if celltype in stromal_types:
        return "Stromal"
    if celltype in epithelial_types:
        return "Epithelial"
    return "Other"

composition["TME_Category"] = composition["celltype_major"].apply(category)

category_summary = (
    composition
    .groupby(["subtype", "TME_Category"], as_index=False)["Cell_Count"]
    .sum()
)
category_summary["Fraction"] = category_summary.groupby("subtype")["Cell_Count"].transform(
    lambda x: x / x.sum()
)
category_summary.to_csv(out_dir / "tme_category_composition_by_subtype.csv", index=False)

print("Final TME summary complete")
print("Saved:")
print(" - results/final/major_celltype_fraction_summary.csv")
print(" - results/final/top_celltypes_by_subtype.csv")
print(" - results/final/top_marker_summary_by_celltype.csv")
print(" - results/final/tme_category_composition_by_subtype.csv")
