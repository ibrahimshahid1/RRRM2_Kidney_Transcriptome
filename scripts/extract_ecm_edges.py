import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

repo = Path('/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome')
run_dir = repo / 'data/results/run_20260219_115338_2500g'
phase2_dir = repo / 'data/processed/networks/run_20260219_115338_2500g'

# Load genes
genes = [g.strip() for g in (phase2_dir / 'phase2_genes.txt').read_text().splitlines() if g.strip()]
gene_to_idx = {g: i for i, g in enumerate(genes)}

ecm_genes = ['Col1a1', 'Col3a1', 'Col4a1', 'Fn1', 'Mmp2']
ecm_indices = {g: gene_to_idx[g] for g in ecm_genes if g in gene_to_idx}
print(f"ECM gene indices: {ecm_indices}")

# Load edges
edge_i = np.load(phase2_dir / 'edge_i.npy')
edge_j = np.load(phase2_dir / 'edge_j.npy')

# Find edges involving ECM genes
ecm_idx_set = set(ecm_indices.values())
edge_mask = np.array([(i in ecm_idx_set) or (j in ecm_idx_set) for i, j in zip(edge_i, edge_j)])
ecm_edge_indices = np.where(edge_mask)[0]
print(f"Found {len(ecm_edge_indices)} edges involving ECM genes")

# Load LIONESS edge weights for those edges only
lioness_path = phase2_dir / 'lioness_edges.npy'
if not lioness_path.exists():
    lioness_path = phase2_dir / 'lioness_z_edges.npy'
W_full = np.load(lioness_path)
W_sub = W_full[:, ecm_edge_indices]

# Load and align metadata
meta = pd.read_csv(repo / 'data/processed/phase1_residuals/meta_phase1.tsv.gz', sep='\t', compression='gzip')
meta['Sample Name'] = meta[meta.columns[0]] if 'Sample Name' not in meta.columns else meta['Sample Name']
meta = meta.set_index('Sample Name')
samples = [s.strip() for s in (phase2_dir / 'lioness_samples.txt').read_text().splitlines() if s.strip()]
meta = meta.loc[samples]
mask_flt_gc = meta['EnvGroup'].isin(['FLT', 'GC'])
W_sub = W_sub[mask_flt_gc.values, :]
meta_sub = meta[mask_flt_gc].copy()

# Design matrix
flight = (meta_sub["EnvGroup"] == "FLT").astype(float).values
age = (meta_sub["Age"] == "OLD").astype(float).values
arm = (meta_sub["Arm"] == "ISS-T").astype(float).values
X = np.column_stack([
    np.ones(len(meta_sub)),
    age, arm, flight,
    age*flight, arm*flight, age*arm, age*arm*flight
])
term_idx = 5 # arm * flight

# Regression
XtX_inv = np.linalg.pinv(X.T @ X)
beta = XtX_inv @ X.T @ W_sub
b = beta[term_idx, :]
predicted = X @ beta
residuals = W_sub - predicted
df_resid = len(meta_sub) - X.shape[1]
var_resid = np.sum(residuals**2, axis=0) / df_resid
se = np.sqrt(var_resid * np.diag(XtX_inv)[term_idx])
t_stat = b / se
p_vals = 2 * (1 - stats.t.cdf(np.abs(t_stat), df_resid))

# BH correction on these specific edges
order = np.argsort(p_vals)
ranked = p_vals[order]
q_vals = ranked * len(p_vals) / np.arange(1, len(p_vals) + 1)
q_vals = np.minimum.accumulate(q_vals[::-1])[::-1]
q = np.empty_like(q_vals)
q[order] = q_vals

# Report top edges
results = []
for idx in np.where(q < 0.05)[0]:
    orig_edge_idx = ecm_edge_indices[idx]
    g1 = genes[edge_i[orig_edge_idx]]
    g2 = genes[edge_j[orig_edge_idx]]
    if g1 not in ecm_genes: g1, g2 = g2, g1
    results.append({
        'gene1': g1, 'gene2': g2,
        'beta_arm_flight': b[idx], 'p': p_vals[idx], 'q': q[idx]
    })

res_df = pd.DataFrame(results).sort_values('q')
print(f"\nTop significant rewiring edges (q < 0.05) for ECM genes (Arm x Flight interaction):")
print(res_df.head(20).to_string(index=False))

# What does a positive/negative interaction beta mean?
# flight=1(FLT), 0(GC); arm=1(ISS-T), 0(LAR).
# arm_flight interaction = (ISS-T FLT - ISS-T GC) - (LAR FLT - LAR GC)
# Positive: The increase in connectivity due to flight is GREATER on ISS-T than LAR.
# Negative: The increase in connectivity due to flight is WEAKER (or decreases) on ISS-T compared to LAR.
