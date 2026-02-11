#!/usr/bin/env python3
"""
Phase 1.5: Dataset-Derived DCT Marker Discovery
=================================================

Identifies genes whose expression tracks the deconvolved DCT fraction,
after controlling for other nephron segments and experimental design.

KEY DESIGN DECISION:
  Uses VST-normalised expression (Y), NOT Rtech.
  Rtech has CLR segment proportions already regressed out (DESeq2.R line 224),
  so β_DCT ≈ 0 on Rtech. We need the pre-segment-residual expression.

Model (per gene g, sample i):
  Y_ig = α + β_DCT · CLR(DCT)_i + γ' · CLR(other)_i + θ_cell(i) + ε_ig

  where cell(i) = Age × Arm × EnvGroup (16 one-hot dummies, 15 df)

β_DCT > 0, q(BH) < 0.05, bootstrap_freq ≥ 0.70 → DCT marker gene

Usage:
    python scripts/discover_dct_markers.py \\
        --vst data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv \\
        --meta data/processed/phase1_residuals/meta_phase1.tsv.gz \\
        --clr data/processed/deconvolution/music_segment_direct_proportions_CLR.csv \\
        --outdir data/processed/dct_markers
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


# ── helpers ──────────────────────────────────────────────────────────────────

def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (NaN-safe)."""
    p = np.asarray(p, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)  # guard against NaN/Inf
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def build_design_matrix(meta: pd.DataFrame, clr: pd.DataFrame,
                        cell_cols: list[str],
                        tech_cols: list[str] | None = None,
                        drop_segment: str = "PT") -> tuple[np.ndarray, int]:
    """
    Build the full design matrix X [n_samples × p]:
      intercept | CLR(DCT) | CLR(other) | tech covariates | cell dummies

    Returns:
        X : (n, p) float64 design matrix
        dct_col_idx : column index of CLR(DCT) in X
    """
    n = len(meta)
    parts = [np.ones((n, 1))]  # intercept
    col_names = ["intercept"]

    # CLR(DCT) — the coefficient of interest
    dct_col_idx = 1
    dct_vals = clr[["DCT"]].values.astype(float)
    # Z-score CLR(DCT)
    dct_mean = np.nanmean(dct_vals)
    dct_std = np.nanstd(dct_vals)
    if dct_std > 1e-8:
        dct_vals = (dct_vals - dct_mean) / dct_std
    parts.append(dct_vals)
    col_names.append("CLR_DCT")

    # CLR(other segments), dropping one for compositionality
    other_cols = [c for c in clr.columns if c != "DCT" and c != drop_segment]
    other_vals = clr[other_cols].values.astype(float)
    # Z-score other CLR columns
    other_mean = np.nanmean(other_vals, axis=0)
    other_std = np.nanstd(other_vals, axis=0)
    # Avoid division by zero
    other_std[other_std < 1e-8] = 1.0
    other_vals = (other_vals - other_mean) / other_std
    parts.append(other_vals)
    col_names.extend([f"CLR_{c}" for c in other_cols])
    
    # Technical covariates (numeric: z-score; factor: one-hot drop-first)
    tech_added = []
    if tech_cols:
        for tc in tech_cols:
            if tc not in meta.columns:
                print(f"  WARNING: tech covariate '{tc}' not found in metadata, skipping")
                continue
            vals = meta[tc]
            if pd.api.types.is_numeric_dtype(vals):
                z = vals.values.astype(float)
                z = np.nan_to_num(z, nan=np.nanmedian(z))
                # Log-scale depth-like variables to tame leverage
                if tc.lower() in {"readdepth", "depth", "reads", "read_depth"}:
                    z = np.log10(np.maximum(z, 1.0))
                sd = z.std()
                if sd > 1e-8:
                    z = (z - z.mean()) / sd
                    parts.append(z.reshape(-1, 1))
                    col_names.append(f"tech_{tc}")
                    tech_added.append(tc)
            else:
                dummies = pd.get_dummies(vals, prefix=f"tech_{tc}",
                                        drop_first=True, dtype=float)
                if dummies.shape[1] > 0:
                    parts.append(dummies.values)
                    col_names.extend(dummies.columns.tolist())
                    tech_added.append(tc)
        if tech_added:
            print(f"  Technical covariates: {tech_added}")

    # Cell fixed effects (Age × Arm × EnvGroup)
    cell_labels = meta[cell_cols].astype(str).agg("|".join, axis=1)
    cell_dummies = pd.get_dummies(cell_labels, drop_first=True, dtype=float)
    parts.append(cell_dummies.values)
    col_names.extend(cell_dummies.columns.tolist())

    X = np.hstack(parts)
    print(f"  Design matrix: {X.shape[0]} samples × {X.shape[1]} predictors")
    print(f"  CLR(DCT) at column index {dct_col_idx}")
    print(f"  Segment CLR columns: DCT + {other_cols} (dropped: {drop_segment})")
    print(f"  Cell levels: {cell_dummies.shape[1] + 1}")
    return X.astype(np.float64), dct_col_idx, col_names


def ols_all_genes(X: np.ndarray, Y: np.ndarray,
                  dct_col: int) -> pd.DataFrame:
    """
    Vectorised OLS: solve X @ B = Y in one shot.

    Args:
        X:   (n, p) design matrix
        Y:   (n, G) expression matrix (samples × genes)
        dct_col: index of CLR(DCT) column in X

    Returns:
        DataFrame with columns: beta_dct, t_stat, p_value, partial_r2
    """
    n, p = X.shape
    G = Y.shape[1]
    rank = np.linalg.matrix_rank(X)

    # Solve via pseudoinverse (precompute once)
    X_pinv = np.linalg.pinv(X)             # (p, n)
    B = X_pinv @ Y                          # (p, G) — all betas at once
    residuals = Y - X @ B                   # (n, G)
    rss = np.sum(residuals ** 2, axis=0)    # (G,)
    df_resid = n - rank                     # use rank, not p

    # Standard error of β_DCT
    sigma2 = rss / max(df_resid, 1)         # (G,)
    # Variance of β̂ = diag(σ² (X'X)^{-1})
    XtX_inv = np.linalg.pinv(X.T @ X)      # (p, p)
    var_beta_dct = XtX_inv[dct_col, dct_col] * sigma2  # (G,)

    beta_dct = B[dct_col, :]               # (G,)
    se_dct = np.sqrt(np.maximum(var_beta_dct, 1e-30))
    t_stat = beta_dct / se_dct             # (G,)
    p_value = 2.0 * sp_stats.t.sf(np.abs(t_stat), max(df_resid, 1))

    # Partial R²: proper reduced model with DCT column DROPPED (not zeroed)
    X_red = np.delete(X, dct_col, axis=1)  # (n, p-1)
    B_red = np.linalg.pinv(X_red) @ Y
    rss_red = np.sum((Y - X_red @ B_red) ** 2, axis=0)
    partial_r2 = 1.0 - rss / np.maximum(rss_red, 1e-30)
    partial_r2 = np.clip(partial_r2, 0.0, 1.0)

    return pd.DataFrame({
        "beta_dct": beta_dct,
        "t_stat": t_stat,
        "p_value": p_value,
        "partial_r2": partial_r2,
    })


def specificity_filter(X: np.ndarray, Y: np.ndarray,
                       clr_col_indices: dict[str, int],
                       dct_col: int,
                       top_n: int = 3) -> np.ndarray:
    """
    For each gene, check if |β_DCT| is among the top-N segment coefficients.

    DCT is a minority cell type in bulk kidney, so requiring it to be THE
    largest beta is too stringent (PT/TAL_LOH dominate). Top-3 ensures
    DCT markers are reasonably specific without being impossible.

    Returns:
        Boolean array (G,): True if β_DCT ranks in top-N segments.
    """
    X_pinv = np.linalg.pinv(X)
    B = X_pinv @ Y  # (p, G)

    seg_cols = list(clr_col_indices.values())
    seg_betas = np.abs(B[seg_cols, :])  # (n_segments, G)

    dct_local_idx = seg_cols.index(dct_col)
    dct_abs = seg_betas[dct_local_idx, :]  # (G,)

    # Rank DCT among segments for each gene (higher rank = larger beta)
    # argsort gives ascending order; we want descending rank
    n_segs = seg_betas.shape[0]
    ranks = np.argsort(np.argsort(seg_betas, axis=0), axis=0)  # (n_segs, G)
    dct_rank = n_segs - ranks[dct_local_idx, :]  # 1 = largest, n = smallest

    return dct_rank <= top_n


def bootstrap_stability(X: np.ndarray, Y: np.ndarray, meta: pd.DataFrame,
                        cell_cols: list[str], dct_col: int,
                        n_boot: int, alpha: float,
                        rng: np.random.Generator) -> np.ndarray:
    """
    Stratified bootstrap within experimental cells.
    Counts gene as passing when BOTH β_DCT > 0 AND q < alpha,
    matching the actual marker selection criterion.
    Returns: (G,) array of fraction of bootstraps where gene passes.
    """
    G = Y.shape[1]
    pass_counts = np.zeros(G, dtype=int)

    cell_labels = meta[cell_cols].astype(str).agg("|".join, axis=1).values
    unique_cells = np.unique(cell_labels)
    cell_indices = {c: np.where(cell_labels == c)[0] for c in unique_cells}

    n, p = X.shape

    for b in range(n_boot):
        # Stratified resample: within each cell, sample with replacement
        boot_idx = []
        for cell, idx in cell_indices.items():
            boot_idx.append(rng.choice(idx, size=len(idx), replace=True))
        boot_idx = np.concatenate(boot_idx)

        X_b = X[boot_idx]
        Y_b = Y[boot_idx]

        # Recompute OLS for this resample
        try:
            X_pinv_b = np.linalg.pinv(X_b)
        except np.linalg.LinAlgError:
            continue
        rank_b = np.linalg.matrix_rank(X_b)
        B_b = X_pinv_b @ Y_b
        resid_b = Y_b - X_b @ B_b
        rss_b = np.sum(resid_b ** 2, axis=0)
        df_b = max(X_b.shape[0] - rank_b, 1)
        sigma2_b = rss_b / df_b

        XtX_inv_b = np.linalg.pinv(X_b.T @ X_b)
        var_beta_b = XtX_inv_b[dct_col, dct_col] * sigma2_b
        se_b = np.sqrt(np.maximum(var_beta_b, 1e-30))
        beta_b = B_b[dct_col, :]
        t_b = beta_b / se_b
        p_b = 2.0 * sp_stats.t.sf(np.abs(t_b), df_b)

        q_b = bh_fdr(p_b)
        # Must match full marker criterion: β > 0 AND q < alpha
        pass_counts += ((beta_b > 0) & (q_b < alpha)).astype(int)

    return pass_counts / n_boot


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Discover DCT marker genes via regression on VST expression"
    )
    ap.add_argument("--vst",
                    default=str(REPO_ROOT / "data/processed/vst_normalized"
                                / "GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"),
                    help="Path to VST expression CSV (genes × samples)")
    ap.add_argument("--meta",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals"
                                / "meta_phase1.tsv.gz"),
                    help="Path to metadata TSV (sample-level)")
    ap.add_argument("--clr",
                    default=str(REPO_ROOT / "data/processed/deconvolution/test16"
                                / "music_segment_direct_proportions_CLR.csv"),
                    help="Path to CLR proportions CSV")
    ap.add_argument("--cell_cols", default="Age,Arm,EnvGroup",
                    help="Comma-separated columns defining experimental cells")
    ap.add_argument("--tech_cols", default="",
                    help="Comma-separated technical covariate columns from metadata "
                         "(e.g. LibraryBatch,ReadDepth,rRNA). Numeric cols are z-scored, "
                         "categorical cols are one-hot encoded.")
    ap.add_argument("--drop_segment", default="PT",
                    help="Segment to drop from CLR for compositionality")
    ap.add_argument("--q_threshold", type=float, default=0.05,
                    help="BH FDR threshold for significance")
    ap.add_argument("--boot_n", type=int, default=200,
                    help="Number of bootstrap resamples for stability")
    ap.add_argument("--boot_freq", type=float, default=0.70,
                    help="Minimum bootstrap frequency to count as stable")
    ap.add_argument("--require_specific", action="store_true",
                    help="Require β_DCT to rank top-3 among segment betas "
                         "(off by default; DCT is ~6%% of kidney bulk, so "
                         "magnitude comparison is misleading)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/processed/dct_markers"),
                    help="Output directory")
    ap.add_argument("--diff_threshold", type=float, default=0.05,
                    help="Anti-confounding threshold: corr(expr, DCT) - max(corr(expr, TAL), corr(expr, CD)) > delta")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print("=" * 60)
    print("Phase 1.5: Dataset-Derived DCT Marker Discovery")
    print("=" * 60)

    # ── 1. Load VST expression (genes × samples) ────────────────────────
    print(f"\n1) Loading VST expression: {args.vst}")
    vst = pd.read_csv(args.vst, index_col=0)
    # Strip Ensembl versions
    vst.index = vst.index.str.replace(r"\.\d+$", "", regex=True)
    if vst.index.duplicated().any():
        # mean of duplicate Ensembl IDs (sum is wrong in VST space)
        vst = vst.groupby(vst.index).mean()
    print(f"   {vst.shape[0]} genes × {vst.shape[1]} samples")

    # ── 2. Load metadata ────────────────────────────────────────────────
    print(f"\n2) Loading metadata: {args.meta}")
    meta = pd.read_csv(args.meta, sep="\t")
    # Find sample column
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            meta = meta.set_index(col, drop=False)
            break

    # ── 3. Load CLR proportions ─────────────────────────────────────────
    print(f"\n3) Loading CLR proportions: {args.clr}")
    clr = pd.read_csv(args.clr, index_col=0)
    print(f"   Segments: {list(clr.columns)}")
    if "DCT" not in clr.columns:
        raise ValueError(f"CLR file must contain 'DCT' column. Found: {list(clr.columns)}")

    # ── 4. Align samples ────────────────────────────────────────────────
    print("\n4) Aligning samples across VST, metadata, CLR...")
    common = sorted(set(vst.columns) & set(meta.index) & set(clr.index))
    if len(common) < 10:
        raise ValueError(f"Only {len(common)} common samples found. Check sample IDs.")
    vst = vst[common]
    meta = meta.loc[common]
    clr = clr.loc[common]
    print(f"   Aligned: {len(common)} samples")

    # ── 5. Filter genes (CPM ≥ 1 in ≥ 20% of samples) ──────────────────
    # VST values are already normalised, but we should filter low-expression
    # Use a variance floor: keep genes with var > 0
    var = vst.var(axis=1)
    keep = var > 1e-8
    vst = vst.loc[keep]
    print(f"\n5) After low-variance filter: {vst.shape[0]} genes")

    # ── 6. Build design matrix ──────────────────────────────────────────
    cell_cols = [c.strip() for c in args.cell_cols.split(",")]
    tech_cols = [c.strip() for c in args.tech_cols.split(",") if c.strip()] if args.tech_cols else None
    print(f"\n6) Building design matrix...")
    X, dct_col, design_col_names = build_design_matrix(meta, clr, cell_cols, tech_cols, args.drop_segment)

    # Check rank
    rank = np.linalg.matrix_rank(X)
    print(f"   Rank: {rank}/{X.shape[1]} (should be full rank)")
    if rank < X.shape[1]:
        print("   WARNING: Design matrix is rank-deficient. Results may be unreliable.")

    # ── 7. Vectorised OLS for all genes ─────────────────────────────────
    print(f"\n7) Running OLS for {vst.shape[0]} genes...")
    Y = vst.values.T  # (n_samples, n_genes)
    gene_names = vst.index.tolist()

    results = ols_all_genes(X, Y, dct_col)
    results.index = gene_names
    results["gene"] = gene_names
    results["q_BH"] = bh_fdr(results["p_value"].values)

    # Marginal correlation with CLR(DCT) — makes Simpson's paradox legible
    from scipy.stats import pearsonr as _pearsonr
    
    # We need the original CLR values for correlation, not the z-scored ones from design matrix
    # But since pearsonr is invariant to linear scaling, it doesn't matter much.
    # However, let's use the raw input 'clr' dataframe to be transparent.
    dct_vec = clr["DCT"].values
    
    # Compute marginal correlations for DCT, TAL, and CD
    print(f"\n   Computing marginal correlations for DCT, TAL, CD...")
    marginal_r = np.array([_pearsonr(Y[:, g], dct_vec)[0] for g in range(Y.shape[1])])
    marginal_p = np.array([_pearsonr(Y[:, g], dct_vec)[1] for g in range(Y.shape[1])])
    results["pearson_r_vst_dct"] = marginal_r
    results["pearson_p_vst_dct"] = marginal_p
    
    # Strict Sign Consistency: β_DCT > 0 AND Corr(expr, DCT) > 0
    results["sign_consistent"] = (results["beta_dct"] > 0) & (results["pearson_r_vst_dct"] > 0)
    
    # Anti-confounding: Corr(DCT) >> Corr(TAL) and Corr(CD)
    # Identify TAL and CD columns. Adjust based on likely column names
    possible_tal = ["TAL", "TAL_LOH", "LOH"]
    possible_cd = ["CD", "CD_PC", "CD_IC", "CNT"] # CD might be split or combined
    
    tal_col = next((c for c in clr.columns if c in possible_tal), None)
    cd_col = next((c for c in clr.columns if c in possible_cd), None)
    
    if tal_col:
        print(f"   Using '{tal_col}' for TAL anti-confounding check")
        tal_vec = clr[tal_col].values
        r_tal = np.array([_pearsonr(Y[:, g], tal_vec)[0] for g in range(Y.shape[1])])
        results["pearson_r_vst_tal"] = r_tal
    else:
        print("   WARNING: Could not find TAL column. Skipping TAL anti-confounding.")
        results["pearson_r_vst_tal"] = -1.0 # dummy
        
    if cd_col:
        print(f"   Using '{cd_col}' for CD anti-confounding check")
        cd_vec = clr[cd_col].values
        r_cd = np.array([_pearsonr(Y[:, g], cd_vec)[0] for g in range(Y.shape[1])])
        results["pearson_r_vst_cd"] = r_cd
    else:
        print("   WARNING: Could not find CD column. Skipping CD anti-confounding.")
        results["pearson_r_vst_cd"] = -1.0
        
    # Delta check
    max_confound = np.maximum(results["pearson_r_vst_tal"], results["pearson_r_vst_cd"])
    results["anti_confounded"] = (results["pearson_r_vst_dct"] - max_confound) > args.diff_threshold
    
    results["sign_flip"] = (results["beta_dct"] > 0) & (results["pearson_r_vst_dct"] < 0)

    sig_pos = (results["beta_dct"] > 0) & (results["q_BH"] < args.q_threshold)
    print(f"   Significant positive β_DCT (q<{args.q_threshold}): {sig_pos.sum()} genes")

    # ── 8. DCT specificity filter ───────────────────────────────────────
    print(f"\n8) Computing DCT specificity (|β_DCT| top-3 among segments)...")
    # Build map of segment column indices
    seg_col_map = {"DCT": dct_col}
    other_cols = [c for c in clr.columns if c != "DCT" and c != args.drop_segment]
    for i, c in enumerate(other_cols):
        seg_col_map[c] = dct_col + 1 + i  # they follow DCT in the design matrix

    is_specific = specificity_filter(X, Y, seg_col_map, dct_col)
    results["dct_top3"] = is_specific
    sig_specific = sig_pos & is_specific
    print(f"   Of {sig_pos.sum()} significant, {sig_specific.sum()} are DCT top-3")

    # ── 9. Bootstrap stability ──────────────────────────────────────────
    print(f"\n9) Bootstrap stability ({args.boot_n} stratified resamples)...")
    boot_freq = bootstrap_stability(
        X, Y, meta, cell_cols, dct_col,
        n_boot=args.boot_n, alpha=args.q_threshold, rng=rng
    )
    results["bootstrap_freq"] = boot_freq
    stable = boot_freq >= args.boot_freq
    print(f"   Stable at ≥{args.boot_freq:.0%}: {stable.sum()} genes")

    # ── 10. Final panel ─────────────────────────────────────────────────
    # Core criteria: β_DCT > 0, q < threshold, bootstrap stable
    # NEW: Sign consistency AND Anti-confounding
    
    # 1. Significant positive beta
    c1 = (results["beta_dct"] > 0) & (results["q_BH"] < args.q_threshold)
    # 2. Bootstrap stability
    c2 = results["bootstrap_freq"] >= args.boot_freq
    # 3. Sign consistency (Marginal R > 0)
    c3 = results["sign_consistent"]
    # 4. Anti-confounding (Marginal R_DCT > R_TAL + delta)
    c4 = results["anti_confounded"]
    
    panel_mask = c1 & c2 & c3 & c4
    
    print(f"\n10) Filtering stats:")
    print(f"    Significant (q<{args.q_threshold}): {c1.sum()}")
    print(f"    Stable (freq>={args.boot_freq}): {c2.sum()}")
    print(f"    Sign Consistent (r>0): {c3.sum()}")
    print(f"    Anti-Confounded (δ>{args.diff_threshold}): {c4.sum()}")
    
    if args.require_specific:
        panel_mask = panel_mask & is_specific
        spec_label = ", DCT-dominant"
    else:
        spec_label = ""
        
    panel_genes = results.loc[panel_mask, "gene"].tolist()
    print(f"\n    FINAL DCT marker panel: {len(panel_genes)} genes")
    print(f"    (β>0, q<{args.q_threshold}, boot≥{args.boot_freq}, sign-consistent, anti-conf{spec_label})")

    # Sort by beta_dct descending
    results = results.sort_values("beta_dct", ascending=False)

    # ── 11. Save outputs ────────────────────────────────────────────────
    cols_order = ["gene", "beta_dct", "t_stat", "p_value", "q_BH",
                  "partial_r2", "dct_top3", "bootstrap_freq",
                  "pearson_r_vst_dct", "sign_consistent", "anti_confounded",
                  "pearson_r_vst_tal", "pearson_r_vst_cd",
                  "pearson_p_vst_dct", "sign_flip"]
    scores_path = outdir / "dct_marker_scores.tsv"
    results[cols_order].to_csv(scores_path, sep="\t", index=False)
    print(f"\n   Wrote full scores: {scores_path}")

    panel_path = outdir / "dct_marker_panel.txt"
    with open(panel_path, "w") as f:
        f.write("\n".join(sorted(panel_genes)) + "\n")
    print(f"   Wrote marker panel ({len(panel_genes)} genes): {panel_path}")

    # Save design matrix column list for reproducibility
    design_path = outdir / "design_cols.txt"
    with open(design_path, "w") as f:
        f.write("\n".join(design_col_names) + "\n")
    print(f"   Wrote design columns ({len(design_col_names)}): {design_path}")

    # ── 12. Quick sanity report ─────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    print(f"  Total genes tested:       {len(gene_names)}")
    print(f"  Significant (β>0, q<{args.q_threshold}): {sig_pos.sum()}")
    print(f"  + DCT top-3:              {sig_specific.sum()}")
    print(f"  + Bootstrap stable:       {panel_mask.sum()}")
    n_flips = results.loc[panel_mask, "sign_flip"].sum()
    print(f"  Sign flips (β>0, r<0):    {n_flips}/{panel_mask.sum()}")
    print(f"\n  Top 10 DCT markers by β_DCT:")
    top = results.loc[panel_mask].head(10)
    for _, row in top.iterrows():
        print(f"    {row['gene']:25s}  β={row['beta_dct']:+.4f}  "
              f"q={row['q_BH']:.2e}  R²={row['partial_r2']:.4f}  "
              f"boot={row['bootstrap_freq']:.2f}")

    print(f"\n[OK] Phase 1.5 complete → {outdir}")


if __name__ == "__main__":
    main()
