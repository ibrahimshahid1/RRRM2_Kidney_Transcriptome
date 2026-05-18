"""Drill-down: which module does each of the 11 LIONESS candidates land in?"""

from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
LIONESS_RUN = REPO / "data/results/run_20260517_213205_2500g"
WGCNA_RUN = REPO / "data/results/run_20260505_remediated_2500g/wgcna"
OUTDIR = LIONESS_RUN / "module_convergence"

perm = pd.read_csv(LIONESS_RUN / "phase6_uncertainty" / "ISS_T_YNG_FLT_minus_GC_perm_pvals.tsv", sep="\t")
sig = perm[perm["q_BH_edge_sum"] < 0.05].copy()

assignments = pd.read_csv(WGCNA_RUN / "module_assignments.csv")
merged = sig.merge(assignments, on="gene", how="left")
merged = merged[["gene", "symbol", "module_color", "edge_sum_node_rewiring_obs", "p_perm_edge_sum", "q_BH_edge_sum"]]
merged = merged.sort_values("q_BH_edge_sum")

# Also flag the grey60-membership of each candidate
print("=== 11 LIONESS ISS-T-Young candidates -> WGCNA module ===")
print(merged.to_string(index=False))

merged.to_csv(OUTDIR / "lioness_11_candidates_module_membership.tsv", sep="\t", index=False)

print("\n=== Distribution across modules ===")
print(merged["module_color"].fillna("OFF_PANEL").value_counts())

# Reverse: which genes are in grey60? Are any LIONESS candidates among them?
grey60 = set(assignments[assignments["module_color"] == "grey60"]["gene"])
lioness_set = set(sig["gene"])
overlap = grey60 & lioness_set
print(f"\nGrey60 size: {len(grey60)}")
print(f"LIONESS candidates: {len(lioness_set)}")
print(f"Overlap LIONESS ∩ grey60: {len(overlap)}")
if overlap:
    print(assignments[assignments["gene"].isin(overlap)][["gene", "symbol", "module_color"]].to_string(index=False))

# Also: which are in purple?
purple = set(assignments[assignments["module_color"] == "purple"]["gene"])
overlap_purple = purple & lioness_set
print(f"\nPurple size: {len(purple)}")
print(f"Overlap LIONESS ∩ purple: {len(overlap_purple)}")
print(assignments[assignments["gene"].isin(overlap_purple)][["gene", "symbol", "module_color"]].to_string(index=False))
