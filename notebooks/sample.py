import scanpy as sc

adata = sc.read_h5ad("data/external/single_cell_atlases/tms_kidney_female_ALLDATASETS_innerGenes.h5ad")

print("X dtype:", adata.X.dtype, "shape:", adata.shape)
print("layers:", list(adata.layers.keys()))
print("raw exists:", adata.raw is not None)
if adata.raw is not None:
    print("raw shape:", adata.raw.shape)
