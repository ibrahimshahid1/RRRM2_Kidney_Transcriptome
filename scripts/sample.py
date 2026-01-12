import scanpy as sc
import numpy as np
import pandas as pd

adata = sc.read_h5ad("data/external/single_cell_atlases/tms_kidney_female_ALLDATASETS_innerGenes.h5ad")
print("n_vars:", adata.n_vars)
print("var_names head:", adata.var_names[:10].tolist())

# show helpful columns that might contain symbols
print("var columns:", adata.var.columns.tolist())

for col in ["gene_symbols", "gene_symbol", "feature_name", "name", "symbol"]:
    if col in adata.var.columns:
        print(f"{col} head:", adata.var[col].astype(str).head(10).tolist())

# quick heuristic
head = adata.var_names[:200].astype(str)
is_ens = (head.str.startswith("ENSMUSG").mean())
is_sym = (head.str.match(r"^[A-Za-z0-9\-]+$").mean())
print("looks like Ensembl:", is_ens, " looks like symbols:", is_sym)
# if X isn't log, make it usable for scoring
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# define signatures
sig = {
    "DCT": ["Slc12a3","Pvalb","Calb1","Trpv5","Fxyd2"],
    "Podocyte": ["Nphs1","Nphs2","Wt1","Synpo","Podxl"],
    "Endothelial": ["Pecam1","Kdr","Emcn","Robo4"],
    "Mesangial_Pericyte": ["Pdgfrb","Rgs5","Des"]
}

for name, genes in sig.items():
    genes_present = [g for g in genes if g in adata.var_names]
    print(name, "genes present:", len(genes_present), genes_present)
    sc.tl.score_genes(adata, genes_present, score_name=f"{name}_score")

# subset to the DCT-labeled cells (based on your obs column)
dct = adata[adata.obs["cell_type"].astype(str).str.contains("distal convoluted tubule", case=False, na=False)].copy()

# summarize scores
summary = dct.obs[[f"{k}_score" for k in sig]].describe(percentiles=[0.5,0.9,0.99]).T
print(summary)

# "Do we have obvious glomerular contamination?"
# check top-scoring "Podocyte" cells inside DCT
top = dct.obs.sort_values("Podocyte_score", ascending=False).head(30)
print(top[["Podocyte_score","DCT_score","Endothelial_score","Mesangial_Pericyte_score"]])
