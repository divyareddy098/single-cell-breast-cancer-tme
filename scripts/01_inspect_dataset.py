import scanpy as sc

adata = sc.read_h5ad("data/raw/breast_cancer_tme.h5ad")

print(adata)

print("\nObservation columns:")
print(list(adata.obs.columns))

print("\nVariable columns:")
print(list(adata.var.columns))

print("\nFirst obs rows:")
print(adata.obs.head())
