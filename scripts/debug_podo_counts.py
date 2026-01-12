#!/usr/bin/env python3
"""
Simple diagnostic for podocyte expression.
"""
import scanpy as sc
import numpy as np
import pandas as pd

h5ad_path = "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"

print("Loading H5AD...")
adata = sc.read_h5ad(h5ad_path)

print(f"\n=== DATA OVERVIEW ===")
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"X dtype: {adata.X.dtype}")
print(f"Has .raw: {adata.raw is not None}")

# Check X values
X = adata.X
if hasattr(X, "data"):
    x_min, x_max = X.data.min(), X.data.max()
else:
    x_min, x_max = X.min(), X.max()
print(f"X min: {x_min:.4f}, X max: {x_max:.4f}")

# Check raw if exists
if adata.raw is not None:
    rawX = adata.raw.X
    if hasattr(rawX, "data"):
        r_min, r_max = rawX.data.min(), rawX.data.max()
    else:
        r_min, r_max = rawX.min(), rawX.max()
    print(f"raw.X min: {r_min:.4f}, raw.X max: {r_max:.4f}")
    print(f"raw.X dtype: {rawX.dtype}")

print(f"\n=== CELL TYPE LABELS ===")
print(adata.obs["free_annotation"].value_counts()[:20])

# Find podocyte cells
podo_idx = adata.obs["free_annotation"].str.contains("podocyte", case=False, na=False)
n_podo = podo_idx.sum()
print(f"\n=== PODOCYTES ===")
print(f"Podocyte cells: {n_podo}")

if n_podo > 0:
    # Check donor distribution of podocytes
    print(f"\nPodocyte donor distribution:")
    print(adata.obs.loc[podo_idx, "donor_id"].value_counts())

# Podocyte markers
markers = ["Nphs1", "Nphs2", "Wt1", "Synpo", "Podxl"]
marker_ens = {
    "Nphs1": "ENSMUSG00000006649",
    "Nphs2": "ENSMUSG00000026602", 
    "Wt1": "ENSMUSG00000016458",
    "Synpo": "ENSMUSG00000043079",
    "Podxl": "ENSMUSG00000029718",
}

print(f"\n=== MARKER GENES ===")
# Check in raw if available
if adata.raw is not None:
    var_names = list(adata.raw.var_names)
else:
    var_names = list(adata.var_names)

for name, ens_id in marker_ens.items():
    found = ens_id in var_names
    print(f"{name} ({ens_id}): {'FOUND' if found else 'NOT FOUND'}")

# Compare expression
if adata.raw is not None and n_podo > 0:
    print(f"\n=== MARKER EXPRESSION (from .raw) ===")
    raw_adata = adata.raw.to_adata()
    
    podo_barcodes = adata.obs_names[podo_idx].tolist()
    non_podo_barcodes = adata.obs_names[~podo_idx].tolist()
    
    for name, ens_id in marker_ens.items():
        if ens_id not in raw_adata.var_names:
            continue
            
        gene_data = raw_adata[:, ens_id].X
        if hasattr(gene_data, "toarray"):
            gene_data = gene_data.toarray().ravel()
        else:
            gene_data = np.array(gene_data).ravel()
        
        # Match to cell barcodes
        podo_expr = gene_data[podo_idx.values] if hasattr(podo_idx, 'values') else gene_data[np.array(podo_idx)]
        non_podo_expr = gene_data[~podo_idx.values] if hasattr(podo_idx, 'values') else gene_data[~np.array(podo_idx)]
        
        print(f"\n{name}:")
        print(f"  Podocytes ({len(podo_expr)} cells): mean={podo_expr.mean():.2f}, max={podo_expr.max():.0f}, %expressing={(podo_expr>0).mean()*100:.1f}%")
        print(f"  Non-podo ({len(non_podo_expr)} cells): mean={non_podo_expr.mean():.4f}, max={non_podo_expr.max():.0f}")

print("\n=== SUMMARY ===")
# Key question: is expression data actually raw counts?
if adata.raw is not None:
    rawX = adata.raw.X
    sample = rawX[0:1000, 0:100]
    if hasattr(sample, "toarray"):
        sample = sample.toarray()
    sample = sample[sample > 0]
    is_integer = np.allclose(sample, np.round(sample))
    print(f"Raw counts look like integers: {is_integer}")
    print(f"Sample raw values: {sample[:10]}")
