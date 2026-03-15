#!/usr/bin/env python3
"""
Comprehensive comparison of two pipeline runs:
  run_20260312_203319_2500g  (Mar 12 — "NEW")
  run_20260226_132416_2500g  (Feb 26 — "OLD")

Sections covered:
  0. Run metadata & configuration
  1. Phase 0 – Deconvolution (cell-type proportions, pseudobulk recovery)
  2. Phase 1 – Residualization / gene counts
  3. Phase 2 – Network skeleton (edge & gene counts)
  4. Phase 3 – Node2Vec embeddings & seed stability
  5. Phase 3 – Procrustes rewiring (aggregated scores, top-50 genes)
  6. Phase 5 – Derived interaction / persistence metrics
  7. Phase 5 – Silent shifters (strict)
  8. Phase 6 – Edge regression (gene-level flight & interaction)
  9. Phase 6 – Permutation / bootstrap uncertainty
 10. Phase 7 – Gene-set enrichment / biological grounding
 11. Figure-level summary tables
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd

# ── paths ────────────────────────────────────────────────────────────────────
BASE = Path("/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome")
RUN_NEW = BASE / "data/results/run_20260312_203319_2500g"
RUN_OLD = BASE / "data/results/run_20260226_132416_2500g"

SEP = "=" * 88

def header(title: str) -> None:
    print(f"\n{SEP}")
    print(f"  {title}")
    print(SEP)


def safe_read(path: Path, **kw) -> pd.DataFrame | None:
    if path.exists():
        return pd.read_csv(path, **kw)
    return None


# ═══════════════════════════════════════════════════════════════════════════════
#  0.  RUN METADATA
# ═══════════════════════════════════════════════════════════════════════════════
header("0. RUN METADATA & CONFIGURATION")

for label, run in [("NEW (Mar 12)", RUN_NEW), ("OLD (Feb 26)", RUN_OLD)]:
    meta_path = run / "run_metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        print(f"\n  [{label}]")
        print(f"    run_id     : {meta['run_id']}")
        print(f"    timestamp  : {meta['timestamp']}")
        print(f"    git_commit : {meta['git_commit']}")
        print(f"    parameters : {meta['parameters']}")
        print(f"    phases     : {meta['phases']}")
        print(f"    skip_r     : {meta['skip_r']}")

# ═══════════════════════════════════════════════════════════════════════════════
#  1.  PHASE 0 — DECONVOLUTION
# ═══════════════════════════════════════════════════════════════════════════════
header("1. PHASE 0 — DECONVOLUTION (cell-type proportions)")

deconv_files = [
    "music_cluster_proportions.csv",
    "music_group_proportions.csv",
    "music_segment_direct_proportions_CLR.csv",
    "pseudobulk_recovery_comparison.csv",
    "pseudobulk_recovery_stats_P_rna.csv",
    "pseudobulk_recovery_stats_P_cell.csv",
    "pseudobulk_est_props.csv",
    "pseudobulk_normalize_comparison.csv",
    "validation_segment_correlations.csv",
    "immune_diagnostic.csv",
]

for fname in deconv_files:
    p_new = RUN_NEW / "deconvolution" / fname
    p_old = RUN_OLD / "deconvolution" / fname
    new_exists = p_new.exists()
    old_exists = p_old.exists()
    tag = ""
    if new_exists and not old_exists:
        tag = "  ← NEW-ONLY"
    elif old_exists and not new_exists:
        tag = "  ← OLD-ONLY"
    elif not new_exists and not old_exists:
        tag = "  ← MISSING BOTH"

    print(f"\n  File: {fname}{tag}")
    if new_exists:
        df = pd.read_csv(p_new)
        print(f"    NEW  shape={df.shape}  cols={list(df.columns[:8])}{'...' if len(df.columns)>8 else ''}")
    if old_exists:
        df = pd.read_csv(p_old)
        print(f"    OLD  shape={df.shape}  cols={list(df.columns[:8])}{'...' if len(df.columns)>8 else ''}")

# Compare deconvolution cluster proportions more deeply
for fname_base in ["music_cluster_proportions", "music_group_proportions",
                    "music_segment_direct_proportions_CLR"]:
    fname = f"{fname_base}.csv"
    p_new = RUN_NEW / "deconvolution" / fname
    p_old = RUN_OLD / "deconvolution" / fname
    if p_new.exists() and p_old.exists():
        df_n = pd.read_csv(p_new)
        df_o = pd.read_csv(p_old)
        # try numeric comparison
        num_n = df_n.select_dtypes(include="number")
        num_o = df_o.select_dtypes(include="number")
        if num_n.shape == num_o.shape and num_n.shape[0] > 0:
            diff = (num_n.values - num_o.values)
            print(f"\n  Δ {fname_base}:")
            print(f"    max |Δ| = {np.nanmax(np.abs(diff)):.6f}")
            print(f"    mean|Δ| = {np.nanmean(np.abs(diff)):.6f}")
            if np.nanmax(np.abs(diff)) == 0:
                print("    → IDENTICAL between runs.")
            else:
                print("    → VALUES DIFFER between runs.")
    elif p_new.exists():
        print(f"\n  {fname_base}: only in NEW run")

# pseudobulk recovery stats
for which in ["P_rna", "P_cell"]:
    fname = f"pseudobulk_recovery_stats_{which}.csv"
    p_new = RUN_NEW / "deconvolution" / fname
    if p_new.exists():
        df = pd.read_csv(p_new)
        print(f"\n  Pseudobulk recovery ({which}) — NEW run:")
        print(df.to_string(index=False, max_rows=20))

# ═══════════════════════════════════════════════════════════════════════════════
#  2.  PHASE 1 — RESIDUALIZATION / GENE COUNTS
# ═══════════════════════════════════════════════════════════════════════════════
header("2. PHASE 1 — RESIDUALIZATION / GENE & SAMPLE COUNTS")

for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    p1 = run / "phase1_residuals"
    rtech = p1 / "Rtech.tsv.gz"
    meta = p1 / "meta_phase1.tsv.gz"
    if rtech.exists():
        df = pd.read_csv(rtech, sep="\t", nrows=5)
        total_genes = pd.read_csv(rtech, sep="\t", usecols=[0]).shape[0]
        total_samples = df.shape[1] - 1  # first col is gene
        print(f"\n  [{label}] Rtech.tsv.gz")
        print(f"    genes   = {total_genes}")
        print(f"    samples = {total_samples}")
    elif label == "NEW":
        # NEW run has phase1; OLD may not
        print(f"\n  [{label}] Rtech.tsv.gz — not found (shared across runs)")
    if meta.exists():
        mdf = pd.read_csv(meta, sep="\t")
        print(f"    meta rows={mdf.shape[0]}, cols={list(mdf.columns)}")

# Gene overlap check
genes_new_path = RUN_NEW / "networks/phase2_genes.txt"
if genes_new_path.exists():
    genes_new = set(genes_new_path.read_text().strip().split("\n"))
    print(f"\n  Phase-2 gene set (NEW): {len(genes_new)} genes")

# DCT markers
for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    dct = run / "dct_markers/dct_marker_panel.txt"
    if dct.exists():
        markers = dct.read_text().strip().split("\n")
        print(f"  DCT marker panel [{label}]: {len(markers)} genes")

# ═══════════════════════════════════════════════════════════════════════════════
#  3.  PHASE 2 — NETWORK SKELETON
# ═══════════════════════════════════════════════════════════════════════════════
header("3. PHASE 2 — NETWORK SKELETON (edges & genes)")

for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    ei = run / "networks/edge_index.tsv"
    sk = run / "networks/skeleton_edges.tsv"
    if ei.exists():
        n_edges_ei = sum(1 for _ in open(ei)) - 1
        print(f"  [{label}] edge_index.tsv       : {n_edges_ei:,} edges")
    if sk.exists():
        n_edges_sk = sum(1 for _ in open(sk)) - 1
        print(f"  [{label}] skeleton_edges.tsv    : {n_edges_sk:,} edges")

    # regression contrasts
    contr = run / "networks/regression/contrasts.json"
    if contr.exists():
        c = json.loads(contr.read_text())
        print(f"  [{label}] contrasts: {list(c.get('effects', {}).keys())}")

# Compare limma regression results
header("3b. PHASE 2 — EDGE-LEVEL LIMMA REGRESSION")
contrasts = [
    "ISS_T_YNG_FLT_minus_GC",
    "ISS_T_OLD_FLT_minus_GC",
    "LAR_YNG_FLT_minus_GC",
    "LAR_OLD_FLT_minus_GC",
    "ISS_T_AgeDep_Flight",
    "LAR_AgeDep_Flight",
    "ISS_minus_LAR_YNG_Flight",
    "ISS_minus_LAR_OLD_Flight",
]

for c in contrasts:
    fname = f"{c}_limma.tsv"
    p_new = RUN_NEW / "networks/regression" / fname
    p_old = RUN_OLD / "networks/regression" / fname
    if not p_new.exists():
        continue
    df_n = pd.read_csv(p_new, sep="\t")
    has_old = p_old.exists()
    print(f"\n  Contrast: {c}")
    print(f"    NEW: {df_n.shape[0]:,} edges")
    # Significance summary
    if "adj.P.Val" in df_n.columns:
        sig_01 = (df_n["adj.P.Val"] < 0.01).sum()
        sig_05 = (df_n["adj.P.Val"] < 0.05).sum()
        sig_10 = (df_n["adj.P.Val"] < 0.10).sum()
        print(f"    NEW signif: FDR<0.01={sig_01:,}  FDR<0.05={sig_05:,}  FDR<0.10={sig_10:,}")
        print(f"    NEW logFC range: [{df_n['logFC'].min():.4f}, {df_n['logFC'].max():.4f}]")
        print(f"    NEW mean |logFC|: {df_n['logFC'].abs().mean():.4f}")

    if has_old:
        df_o = pd.read_csv(p_old, sep="\t")
        print(f"    OLD: {df_o.shape[0]:,} edges")
        if "adj.P.Val" in df_o.columns:
            sig_01_o = (df_o["adj.P.Val"] < 0.01).sum()
            sig_05_o = (df_o["adj.P.Val"] < 0.05).sum()
            sig_10_o = (df_o["adj.P.Val"] < 0.10).sum()
            print(f"    OLD signif: FDR<0.01={sig_01_o:,}  FDR<0.05={sig_05_o:,}  FDR<0.10={sig_10_o:,}")
            print(f"    OLD logFC range: [{df_o['logFC'].min():.4f}, {df_o['logFC'].max():.4f}]")
            print(f"    OLD mean |logFC|: {df_o['logFC'].abs().mean():.4f}")
        # logFC correlation
        if df_n.shape[0] == df_o.shape[0]:
            r = np.corrcoef(df_n["logFC"].values, df_o["logFC"].values)[0, 1]
            print(f"    logFC Pearson r = {r:.6f}")
            logfc_diff = (df_n["logFC"] - df_o["logFC"]).abs()
            print(f"    max |Δ logFC| = {logfc_diff.max():.6f}")
            print(f"    mean|Δ logFC| = {logfc_diff.mean():.6f}")

# ═══════════════════════════════════════════════════════════════════════════════
#  4.  PHASE 3 — NODE2VEC EMBEDDING STABILITY
# ═══════════════════════════════════════════════════════════════════════════════
header("4. PHASE 3 — NODE2VEC EMBEDDING SEED STABILITY")

for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    emb_sum = run / "phase3_embeddings/embedding_seed_summary.tsv"
    if emb_sum.exists():
        df = pd.read_csv(emb_sum, sep="\t")
        print(f"\n  [{label}] embedding_seed_summary.tsv:")
        print(df.to_string(index=False))

# Per-prediction seed stability
preds = [
    "Pred_YNG_ISS_T_FLT", "Pred_YNG_ISS_T_GC",
    "Pred_OLD_ISS_T_FLT", "Pred_OLD_ISS_T_GC",
    "Pred_YNG_LAR_FLT", "Pred_YNG_LAR_GC",
    "Pred_OLD_LAR_FLT", "Pred_OLD_LAR_GC",
]

print("\n  Per-prediction seed stability (median cosine similarity):")
for pred in preds:
    row = []
    for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
        ss = run / f"phase3_embeddings/{pred}/seed_stability.tsv"
        if ss.exists():
            df = pd.read_csv(ss, sep="\t")
            if "cosine_sim" in df.columns:
                row.append(f"{label}={df['cosine_sim'].median():.4f}")
            elif len(df.columns) >= 2:
                vals = df.iloc[:, 1].astype(float)
                row.append(f"{label}={vals.median():.4f}")
    print(f"    {pred:30s} {'  '.join(row)}")

# ═══════════════════════════════════════════════════════════════════════════════
#  5.  PHASE 3 — PROCRUSTES REWIRING
# ═══════════════════════════════════════════════════════════════════════════════
header("5. PHASE 3 — PROCRUSTES REWIRING SCORES")

for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    proc = run / "phase3_rewiring/procrustes_summary.tsv"
    if proc.exists():
        df = pd.read_csv(proc, sep="\t")
        print(f"\n  [{label}] procrustes_summary.tsv:")
        print(df.to_string(index=False))

# Aggregated rewiring comparisons
rewiring_contrasts = [
    "ISS_T_YNG_FLT_minus_GC",
    "ISS_T_OLD_FLT_minus_GC",
    "LAR_YNG_FLT_minus_GC",
    "LAR_OLD_FLT_minus_GC",
]

print("\n  Aggregated rewiring score distributions:")
for rc in rewiring_contrasts:
    fname = f"{rc}_rewiring_agg.tsv"
    rows_data = []
    for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
        p = run / f"phase3_rewiring/{fname}"
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            # find rewiring score column
            score_col = None
            for c in ["rewiring_mean", "mean_displacement", "rewiring_score", "displacement"]:
                if c in df.columns:
                    score_col = c
                    break
            if score_col is None and len(df.columns) >= 2:
                score_col = df.columns[1]
            if score_col:
                vals = df[score_col].dropna()
                rows_data.append(
                    f"    {label}: n={len(df):,}  mean={vals.mean():.4f}  "
                    f"std={vals.std():.4f}  max={vals.max():.4f}  "
                    f"median={vals.median():.4f}"
                )
    if rows_data:
        print(f"\n  {rc}:")
        for r in rows_data:
            print(r)

# Top-50 rewired genes comparison
print("\n  Top-50 rewired genes overlap (figures/):")
for rc in rewiring_contrasts:
    fname = f"{rc}_top50_rewired.tsv"
    p_new = RUN_NEW / f"figures/phase3_rewiring/{fname}"
    p_old = RUN_OLD / f"figures/phase3_rewiring/{fname}"
    if p_new.exists() and p_old.exists():
        df_n = pd.read_csv(p_new, sep="\t")
        df_o = pd.read_csv(p_old, sep="\t")
        gcol = df_n.columns[0]
        genes_n = set(df_n[gcol])
        genes_o = set(df_o[gcol])
        overlap = genes_n & genes_o
        pct = 100 * len(overlap) / max(len(genes_n), 1)
        print(f"    {rc}: {len(overlap)}/{len(genes_n)} overlap ({pct:.0f}%)")
        only_new = genes_n - genes_o
        only_old = genes_o - genes_n
        if only_new:
            print(f"      NEW-only: {sorted(only_new)[:10]}{'...' if len(only_new)>10 else ''}")
        if only_old:
            print(f"      OLD-only: {sorted(only_old)[:10]}{'...' if len(only_old)>10 else ''}")
    elif p_new.exists():
        print(f"    {rc}: only in NEW")

# rewiring summary (figures/)
for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    rs = run / "figures/phase3_rewiring/rewiring_summary.tsv"
    if rs.exists():
        df = pd.read_csv(rs, sep="\t")
        print(f"\n  [{label}] rewiring_summary.tsv (figures/):")
        print(df.to_string(index=False))

# ═══════════════════════════════════════════════════════════════════════════════
#  6.  PHASE 5 — DERIVED INTERACTION / PERSISTENCE
# ═══════════════════════════════════════════════════════════════════════════════
header("6. PHASE 5 — DERIVED INTERACTION & PERSISTENCE METRICS")

derived_files = [
    "ISS_T_interaction.tsv",
    "LAR_interaction.tsv",
    "ISS_minus_LAR_YNG_persistence.tsv",
    "ISS_minus_LAR_OLD_persistence.tsv",
]

for fname in derived_files:
    p_new = RUN_NEW / f"phase5_derived/{fname}"
    p_old = RUN_OLD / f"phase5_derived/{fname}"
    print(f"\n  {fname}:")
    for label, p in [("NEW", p_new), ("OLD", p_old)]:
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            # find score column
            num_cols = df.select_dtypes(include="number").columns.tolist()
            if num_cols:
                sc = num_cols[0]
                vals = df[sc].dropna()
                print(
                    f"    {label}: n={len(df):,}  {sc} mean={vals.mean():.4f}  "
                    f"std={vals.std():.4f}  [min={vals.min():.4f}, max={vals.max():.4f}]"
                )
    # correlation if both exist
    if p_new.exists() and p_old.exists():
        df_n = pd.read_csv(p_new, sep="\t")
        df_o = pd.read_csv(p_old, sep="\t")
        if df_n.shape == df_o.shape:
            num_n = df_n.select_dtypes(include="number")
            num_o = df_o.select_dtypes(include="number")
            for c in num_n.columns:
                if c in num_o.columns:
                    r = np.corrcoef(num_n[c].values, num_o[c].values)[0, 1]
                    d = (num_n[c] - num_o[c]).abs()
                    print(f"    Δ {c}: Pearson r={r:.6f}  max|Δ|={d.max():.6f}  mean|Δ|={d.mean():.6f}")

# ═══════════════════════════════════════════════════════════════════════════════
#  7.  PHASE 5 — SILENT SHIFTERS (STRICT)
# ═══════════════════════════════════════════════════════════════════════════════
header("7. PHASE 5 — SILENT SHIFTERS (strict)")

ss_contrasts = [
    "ISS_T_YNG_FLT_minus_GC",
    "ISS_T_OLD_FLT_minus_GC",
    "LAR_YNG_FLT_minus_GC",
    "LAR_OLD_FLT_minus_GC",
]

for rc in ss_contrasts:
    fname = f"{rc}_silent_shifters.tsv"
    print(f"\n  {rc}:")
    for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
        p = run / f"phase5_silent_shifters_strict/{fname}"
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            print(f"    [{label}] {len(df):,} silent shifters  cols={list(df.columns)}")
            if len(df) > 0:
                # show first few gene names
                gene_col = df.columns[0]
                print(f"      top genes: {list(df[gene_col].head(5))}")

    # overlap
    p_new = RUN_NEW / f"phase5_silent_shifters_strict/{fname}"
    p_old = RUN_OLD / f"phase5_silent_shifters_strict/{fname}"
    if p_new.exists() and p_old.exists():
        df_n = pd.read_csv(p_new, sep="\t")
        df_o = pd.read_csv(p_old, sep="\t")
        gc = df_n.columns[0]
        overlap = set(df_n[gc]) & set(df_o[gc])
        print(f"    Overlap: {len(overlap)} genes  ({100*len(overlap)/max(len(df_n),1):.0f}% of NEW)")
        new_only = set(df_n[gc]) - set(df_o[gc])
        old_only = set(df_o[gc]) - set(df_n[gc])
        if new_only:
            print(f"    NEW-only ({len(new_only)}): {sorted(new_only)[:8]}{'...' if len(new_only)>8 else ''}")
        if old_only:
            print(f"    OLD-only ({len(old_only)}): {sorted(old_only)[:8]}{'...' if len(old_only)>8 else ''}")

# Figures summary
for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    ss_sum = run / "figures/phase5_silent_shifters/silent_shifter_summary.tsv"
    if ss_sum.exists():
        df = pd.read_csv(ss_sum, sep="\t")
        print(f"\n  [{label}] silent_shifter_summary (figures/):")
        print(df.to_string(index=False))

# ═══════════════════════════════════════════════════════════════════════════════
#  8.  PHASE 6 — GENE-LEVEL EDGE REGRESSION
# ═══════════════════════════════════════════════════════════════════════════════
header("8. PHASE 6 — GENE-LEVEL EDGE REGRESSION")

phase6_files = [
    "gene_flight_effect.tsv",
    "gene_age_flight_interaction.tsv",
    "gene_arm_flight_interaction.tsv",
    "gene_age_arm_flight_3way.tsv",
]

for fname in phase6_files:
    p_new = RUN_NEW / f"phase6_regression/{fname}"
    p_old = RUN_OLD / f"phase6_regression/{fname}"
    print(f"\n  {fname}:")

    for label, p in [("NEW", p_new), ("OLD", p_old)]:
        if p.exists():
            df = pd.read_csv(p, sep="\t")
            print(f"    [{label}] {len(df):,} genes  cols={list(df.columns)}")
            # find p-value and FDR columns
            for pcol in ["p_value", "pvalue", "P.Value", "pval"]:
                if pcol in df.columns:
                    break
            else:
                pcol = None
            for qcol in ["q_BH", "fdr", "adj.P.Val", "qvalue", "padj"]:
                if qcol in df.columns:
                    break
            else:
                qcol = None

            if qcol:
                sig_01 = (df[qcol] < 0.01).sum()
                sig_05 = (df[qcol] < 0.05).sum()
                sig_10 = (df[qcol] < 0.10).sum()
                print(f"      FDR<0.01={sig_01}  FDR<0.05={sig_05}  FDR<0.10={sig_10}")

            # effect size column
            for ecol in ["beta", "coef", "logFC", "effect_size", "estimate"]:
                if ecol in df.columns:
                    break
            else:
                ecol = None
            if ecol:
                v = df[ecol].dropna()
                print(f"      {ecol}: mean={v.mean():.4f}  std={v.std():.4f}  [{v.min():.4f}, {v.max():.4f}]")

    # Cross-run correlation
    if p_new.exists() and p_old.exists():
        df_n = pd.read_csv(p_new, sep="\t")
        df_o = pd.read_csv(p_old, sep="\t")
        gcol_n = df_n.columns[0]
        gcol_o = df_o.columns[0]
        merged = df_n.merge(df_o, on=gcol_n, suffixes=("_new", "_old"))
        if len(merged) > 0:
            print(f"    Merged: {len(merged):,} genes in common")
            # find shared numeric columns
            for base_col in ["beta", "coef", "logFC", "effect_size", "estimate",
                             "p_value", "pvalue", "q_BH", "fdr"]:
                cn = f"{base_col}_new"
                co = f"{base_col}_old"
                if cn in merged.columns and co in merged.columns:
                    r = np.corrcoef(merged[cn].values, merged[co].values)[0, 1]
                    d = (merged[cn] - merged[co]).abs()
                    print(f"    {base_col}: Pearson r={r:.6f}  max|Δ|={d.max():.6f}  mean|Δ|={d.mean():.6f}")

# Figures-level regression summaries
for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    for figf in ["regression_summary.tsv", "significance_category_counts.tsv",
                 "flight_effect_top_significant.tsv",
                 "age_flight_interaction_top_significant.tsv"]:
        fp = run / f"figures/phase6_regression/{figf}"
        if fp.exists():
            df = pd.read_csv(fp, sep="\t")
            print(f"\n  [{label}] figures/phase6_regression/{figf}:")
            print(df.to_string(index=False, max_rows=30))

# ═══════════════════════════════════════════════════════════════════════════════
#  9.  PHASE 6 — PERMUTATION / BOOTSTRAP UNCERTAINTY
# ═══════════════════════════════════════════════════════════════════════════════
header("9. PHASE 6 — PERMUTATION & BOOTSTRAP UNCERTAINTY")

for rc in ss_contrasts:
    for suffix in ["perm_pvals", "bootstrap_ci"]:
        fname = f"{rc}_{suffix}.tsv"
        print(f"\n  {fname}:")
        for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
            p = run / f"phase6_uncertainty/{fname}"
            if p.exists():
                df = pd.read_csv(p, sep="\t")
                print(f"    [{label}] {len(df):,} genes  cols={list(df.columns[:6])}")
                if suffix == "perm_pvals":
                    for pc in ["perm_pval", "pval", "p_perm", "p_value"]:
                        if pc in df.columns:
                            sig = (df[pc] < 0.05).sum()
                            print(f"      perm p<0.05: {sig} genes ({100*sig/len(df):.1f}%)")
                            break
                elif suffix == "bootstrap_ci":
                    for lc, uc in [("ci_lo", "ci_hi"), ("lower", "upper"),
                                   ("ci_2.5", "ci_97.5"), ("ci_lower", "ci_upper")]:
                        if lc in df.columns and uc in df.columns:
                            excludes_zero = ((df[lc] > 0) | (df[uc] < 0)).sum()
                            print(f"      CI excludes 0: {excludes_zero} genes ({100*excludes_zero/len(df):.1f}%)")
                            break

# Figures uncertainty summary
for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
    fp = run / "figures/phase6_uncertainty/uncertainty_summary.tsv"
    if fp.exists():
        df = pd.read_csv(fp, sep="\t")
        print(f"\n  [{label}] uncertainty_summary.tsv:")
        print(df.to_string(index=False))

# ═══════════════════════════════════════════════════════════════════════════════
# 10.  PHASE 7 — GENE-SET ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════
header("10. PHASE 7 — GENE-SET ENRICHMENT / BIOLOGICAL GROUNDING")

for rc in ss_contrasts:
    for suffix in ["gene_set_enrichment", "gene_set_enrichment_significant"]:
        fname = f"{suffix}.tsv"
        print(f"\n  {rc} / {suffix}:")
        for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
            p = run / f"phase7_grounding/{rc}/{fname}"
            if p.exists():
                df = pd.read_csv(p, sep="\t")
                print(f"    [{label}] {len(df)} sets  cols={list(df.columns[:6])}")
                if len(df) > 0:
                    # show top enriched
                    for pc in ["p_value", "pvalue", "p_adj", "FDR", "padj", "fdr"]:
                        if pc in df.columns:
                            top = df.nsmallest(5, pc)
                            gene_set_col = df.columns[0]
                            for _, row in top.iterrows():
                                print(f"      {row[gene_set_col]:50s}  {pc}={row[pc]:.4e}")
                            break
            else:
                if suffix == "gene_set_enrichment_significant":
                    # only printed if exists
                    pass
                else:
                    print(f"    [{label}] NOT FOUND")

# Cross-run enrichment comparison
print("\n  Cross-run gene-set enrichment overlap (significant sets):")
for rc in ss_contrasts:
    fname = "gene_set_enrichment_significant.tsv"
    p_new = RUN_NEW / f"phase7_grounding/{rc}/{fname}"
    p_old = RUN_OLD / f"phase7_grounding/{rc}/{fname}"
    if p_new.exists() and p_old.exists():
        df_n = pd.read_csv(p_new, sep="\t")
        df_o = pd.read_csv(p_old, sep="\t")
        gc = df_n.columns[0]
        if gc in df_o.columns:
            sets_n = set(df_n[gc])
            sets_o = set(df_o[gc])
            overlap = sets_n & sets_o
            print(f"    {rc}: NEW={len(sets_n)}  OLD={len(sets_o)}  overlap={len(overlap)}")
            new_only = sets_n - sets_o
            old_only = sets_o - sets_n
            if new_only:
                print(f"      NEW-only: {sorted(new_only)[:5]}")
            if old_only:
                print(f"      OLD-only: {sorted(old_only)[:5]}")
    elif p_new.exists():
        df_n = pd.read_csv(p_new, sep="\t")
        print(f"    {rc}: NEW={len(df_n)} significant sets  (OLD has none)")
    elif p_old.exists():
        df_o = pd.read_csv(p_old, sep="\t")
        print(f"    {rc}: OLD={len(df_o)} significant sets  (NEW has none)")
    else:
        print(f"    {rc}: no significant sets in either run")

# ═══════════════════════════════════════════════════════════════════════════════
# 11.  FIGURES-LEVEL DECONVOLUTION SUMMARIES
# ═══════════════════════════════════════════════════════════════════════════════
header("11. FIGURES — DECONVOLUTION SUMMARIES")

for sumf in ["music_cluster_proportions_summary.tsv",
             "music_group_proportions_summary.tsv",
             "music_segment_direct_proportions_CLR_summary.tsv"]:
    print(f"\n  {sumf}:")
    for label, run in [("NEW", RUN_NEW), ("OLD", RUN_OLD)]:
        fp = run / f"figures/phase0_deconvolution/{sumf}"
        if fp.exists():
            df = pd.read_csv(fp, sep="\t")
            print(f"    [{label}]:")
            print(df.to_string(index=False))

# ═══════════════════════════════════════════════════════════════════════════════
# 12.  OVERALL STRUCTURAL DIFF — files unique to each run
# ═══════════════════════════════════════════════════════════════════════════════
header("12. STRUCTURAL DIFFERENCES — files unique to each run")

def rel_files(run: Path) -> set[str]:
    all_files = set()
    for root, dirs, files in os.walk(run):
        for f in files:
            if f.endswith(".npy") or f.endswith(".png") or f.endswith(".rds"):
                continue
            rel = os.path.relpath(os.path.join(root, f), run)
            all_files.add(rel)
    return all_files

files_new = rel_files(RUN_NEW)
files_old = rel_files(RUN_OLD)

only_new = sorted(files_new - files_old)
only_old = sorted(files_old - files_new)

print(f"\n  Files in BOTH runs: {len(files_new & files_old)}")
print(f"  Files ONLY in NEW:  {len(only_new)}")
for f in only_new:
    print(f"    + {f}")
print(f"  Files ONLY in OLD:  {len(only_old)}")
for f in only_old:
    print(f"    - {f}")

# ═══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 88}")
print("  COMPARISON COMPLETE")
print(f"{'=' * 88}\n")
