#!/usr/bin/env python3
"""Phase 1.5: Dataset-Derived DCT Marker Discovery"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from src.common import REPO_ROOT, find_sample_col, bh_fdr as _bh_fdr


# helpers

def bh_fdr(p: np.ndarray) -> np.ndarray:
    """BH-FDR with NaN/Inf guard (wraps src.common.bh_fdr)."""
    p = np.asarray(p, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    return _bh_fdr(p)


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


def bootstrap_marginal_stability(Y: np.ndarray, dct_vec: np.ndarray,
                                 meta: pd.DataFrame, cell_cols: list[str],
                                 n_boot: int, alpha: float,
                                 rng: np.random.Generator) -> np.ndarray:
    """
    Fallback bootstrap using marginal Pearson correlation instead of OLS.
    Used when the full OLS model is over-parameterised (≈0 genes pass q<α).

    For each bootstrap resample, compute Pearson r(gene, DCT) and test
    whether r > 0 with BH-corrected p < α.

    Returns: (G,) array of fraction of bootstraps where gene passes.
    """
    from scipy.stats import pearsonr as _pearsonr

    G = Y.shape[1]
    pass_counts = np.zeros(G, dtype=int)

    cell_labels = meta[cell_cols].astype(str).agg("|".join, axis=1).values
    unique_cells = np.unique(cell_labels)
    cell_indices = {c: np.where(cell_labels == c)[0] for c in unique_cells}

    for b in range(n_boot):
        boot_idx = []
        for cell, idx in cell_indices.items():
            boot_idx.append(rng.choice(idx, size=len(idx), replace=True))
        boot_idx = np.concatenate(boot_idx)

        dct_b = dct_vec[boot_idx]
        Y_b = Y[boot_idx]

        # Vectorised Pearson correlation
        n = len(boot_idx)
        dct_c = dct_b - dct_b.mean()
        Y_c = Y_b - Y_b.mean(axis=0)
        dct_ss = np.dot(dct_c, dct_c)
        Y_ss = np.sum(Y_c ** 2, axis=0)
        r_b = (dct_c @ Y_c) / np.sqrt(np.maximum(dct_ss * Y_ss, 1e-30))

        # t-test for correlation
        t_b = r_b * np.sqrt((n - 2) / np.maximum(1 - r_b ** 2, 1e-30))
        p_b = 2.0 * sp_stats.t.sf(np.abs(t_b), n - 2)
        q_b = bh_fdr(p_b)

        pass_counts += ((r_b > 0) & (q_b < alpha)).astype(int)

    return pass_counts / n_boot


# main

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
                    default=str(REPO_ROOT / "data/processed/deconvolution/latest"
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
    ap.add_argument("--diff_threshold", type=float, default=0.0,
                    help="Anti-confounding threshold: min corr(expr, DCT_resid) after "
                         "regressing out TAL/CD from DCT (default: 0.0 = positive correlation only)")
    ap.add_argument("--min_r", type=float, default=0.30,
                    help="Minimum marginal Pearson r(gene, DCT) to keep in fallback mode "
                         "(default: 0.30 — eliminates weak correlations that survive BH "
                         "correction at n=80)")
    ap.add_argument("--fallback_delta", type=float, default=0.10,
                    help="Stricter anti-confounding delta used ONLY in fallback mode "
                         "(default: 0.10 — gene must track DCT ≥0.10 more than any "
                         "confounder). Primary OLS path still uses --diff_threshold.")
    ap.add_argument("--max_panel", type=int, default=200,
                    help="Maximum marker panel size. If more genes pass all filters, "
                         "keep the top --max_panel ranked by delta (default: 200).")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    print("=" * 60)
    print("Phase 1.5: Dataset-Derived DCT Marker Discovery")
    print("=" * 60)

    # 1. Load VST expression (genes  samples)
    print(f"\n1) Loading VST expression: {args.vst}")
    vst = pd.read_csv(args.vst, index_col=0)
    # Strip Ensembl versions
    vst.index = vst.index.str.replace(r"\.\d+$", "", regex=True)
    if vst.index.duplicated().any():
        # mean of duplicate Ensembl IDs (sum is wrong in VST space)
        vst = vst.groupby(vst.index).mean()
    print(f"   {vst.shape[0]} genes × {vst.shape[1]} samples")

    # 2. Load metadata
    print(f"\n2) Loading metadata: {args.meta}")
    meta = pd.read_csv(args.meta, sep="\t")
    # Find sample column
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)

    # 3. Load CLR proportions
    print(f"\n3) Loading CLR proportions: {args.clr}")
    clr = pd.read_csv(args.clr, index_col=0)
    print(f"   Segments: {list(clr.columns)}")

    # If no single 'DCT' column, synthesize from DCT subtypes.
    if "DCT" not in clr.columns:
        # Look for DCT subtypes to combine
        dct_parts = [c for c in clr.columns if c.startswith("DCT")]
        # CNT is part of the distal convoluted tubule continuum
        if "CNT" in clr.columns:
            dct_parts.append("CNT")
        if dct_parts:
            # CLR values are log-ratio transformed, so we need to: 1) back-transform to proportions, 2) sum, 3) re-CLR
            import warnings
            dct_exp = np.exp(clr[dct_parts].values)  # undo log
            dct_combined_exp = dct_exp.sum(axis=1)    # sum proportional parts
            # Re-log to get a combined CLR-like score
            clr["DCT"] = np.log(np.maximum(dct_combined_exp, 1e-12))
            print(f"   Synthesized 'DCT' column from {dct_parts}")
            print(f"   DCT CLR range: [{clr['DCT'].min():.3f}, {clr['DCT'].max():.3f}]")
            # Check if there's actual variance
            dct_std = clr["DCT"].std()
            if dct_std < 1e-6:
                print(f"   WARNING: Synthesized DCT has near-zero variance (std={dct_std:.2e}).")
                print(f"   This likely means deconvolution assigned ~0 to all DCT subtypes.")
                print(f"   Marker discovery will proceed but may find no markers.")
        else:
            raise ValueError(
                f"CLR file must contain 'DCT' column or DCT subtypes (DCT1/DCT2/CNT). "
                f"Found: {list(clr.columns)}"
            )

    # 4. Align samples
    print("\n4) Aligning samples across VST, metadata, CLR...")
    common = sorted(set(vst.columns) & set(meta.index) & set(clr.index))
    if len(common) < 10:
        raise ValueError(f"Only {len(common)} common samples found. Check sample IDs.")
    vst = vst[common]
    meta = meta.loc[common]
    clr = clr.loc[common]
    print(f"   Aligned: {len(common)} samples")

    #  5. Filter genes (CPM  1 in  20% of samples)
    var = vst.var(axis=1)
    keep = var > 1e-8
    vst = vst.loc[keep]
    print(f"\n5) After low-variance filter: {vst.shape[0]} genes")

    # 6. Build design matrix
    cell_cols = [c.strip() for c in args.cell_cols.split(",")]
    tech_cols = [c.strip() for c in args.tech_cols.split(",") if c.strip()] if args.tech_cols else None
    print(f"\n6) Building design matrix...")
    X, dct_col, design_col_names = build_design_matrix(meta, clr, cell_cols, tech_cols, args.drop_segment)

    # Check rank
    rank = np.linalg.matrix_rank(X)
    print(f"   Rank: {rank}/{X.shape[1]} (should be full rank)")
    if rank < X.shape[1]:
        print("   WARNING: Design matrix is rank-deficient. Results may be unreliable.")

    # 7. Vectorised OLS for all genes
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
    dct_vec = clr["DCT"].values
    
    # Compute marginal correlations for DCT
    print(f"\n   Computing marginal correlations for DCT...")
    marginal_r = np.array([_pearsonr(Y[:, g], dct_vec)[0] for g in range(Y.shape[1])])
    marginal_p = np.array([_pearsonr(Y[:, g], dct_vec)[1] for g in range(Y.shape[1])])
    results["pearson_r_vst_dct"] = marginal_r
    results["pearson_p_vst_dct"] = marginal_p
    
    # Strict Sign Consistency: β_DCT > 0 AND Corr(expr, DCT) > 0
    results["sign_consistent"] = (results["beta_dct"] > 0) & (results["pearson_r_vst_dct"] > 0)
    
    #  Anti-confounding via DELTA approach

    possible_tal = ["TAL", "TAL_LOH", "LOH"]
    possible_cd = ["CD", "CD_PC", "CD_IC"]
    
    tal_col = next((c for c in clr.columns if c in possible_tal), None)
    cd_col = next((c for c in clr.columns if c in possible_cd), None)
    
    # Compute per-gene correlations with confounders
    confound_names = []
    confound_r = {}
    if tal_col:
        confound_names.append(tal_col)
        tal_vec = clr[tal_col].values
        confound_r[tal_col] = np.array([_pearsonr(Y[:, g], tal_vec)[0] for g in range(Y.shape[1])])
    if cd_col:
        confound_names.append(cd_col)
        cd_vec_ = clr[cd_col].values
        confound_r[cd_col] = np.array([_pearsonr(Y[:, g], cd_vec_)[0] for g in range(Y.shape[1])])
    
    if confound_names:
        print(f"\n   Anti-confounding: delta approach (DCT corr vs max confounder corr)")
        print(f"   Confounders: {confound_names}")
        from scipy.stats import spearmanr as _spearmanr
        for cn in confound_names:
            rho, _ = _spearmanr(dct_vec, clr[cn].values)
            print(f"   Spearman(DCT, {cn}) = {rho:.3f}")
        
        # Max |confounder correlation| per gene
        confound_stack = np.column_stack([confound_r[c] for c in confound_names])
        max_confound_r = np.max(np.abs(confound_stack), axis=1)  # (G,)
        
        # Anti-confounded = gene tracks DCT MORE than it tracks any confounder
        delta = marginal_r - max_confound_r
        results["delta_confound"] = delta
        results["anti_confounded"] = delta > args.diff_threshold
        
        # Also store a "DCT_resid" correlation for reference (using partial correlation)
        results["pearson_r_dct_resid"] = delta  # use delta as proxy for DCT-specific signal
        
        # Diagnostic: top 20 by delta
        top_idx = np.argsort(delta)[-20:][::-1]
        print(f"\n   Top 20 genes by delta(corr_DCT - max_corr_confounder):")
        for idx in top_idx:
            conf_str = ", ".join([f"r_{cn}={confound_r[cn][idx]:+.4f}" for cn in confound_names])
            print(f"     {gene_names[idx]:25s}  r_DCT={marginal_r[idx]:+.4f}  {conf_str}  δ={delta[idx]:+.4f}")
        
        n_pass = results["anti_confounded"].sum()
        print(f"\n   Genes passing delta > {args.diff_threshold}: {n_pass} / {len(gene_names)}")
    else:
        print("   WARNING: No TAL/CD columns found. Skipping anti-confounding (all pass).")
        results["pearson_r_dct_resid"] = marginal_r
        results["delta_confound"] = 1.0
        results["anti_confounded"] = True
    
    # Also store marginal TAL/CD correlations for reference (reuse computed values)
    if tal_col and tal_col in confound_r:
        results["pearson_r_vst_tal"] = confound_r[tal_col]
    elif tal_col:
        tal_vec = clr[tal_col].values
        results["pearson_r_vst_tal"] = np.array([_pearsonr(Y[:, g], tal_vec)[0] for g in range(Y.shape[1])])
    else:
        results["pearson_r_vst_tal"] = -1.0
    if cd_col and cd_col in confound_r:
        results["pearson_r_vst_cd"] = confound_r[cd_col]
    elif cd_col:
        cd_vec_ = clr[cd_col].values
        results["pearson_r_vst_cd"] = np.array([_pearsonr(Y[:, g], cd_vec_)[0] for g in range(Y.shape[1])])
    else:
        results["pearson_r_vst_cd"] = -1.0
    
    results["sign_flip"] = (results["beta_dct"] > 0) & (results["pearson_r_vst_dct"] < 0)

    sig_pos = (results["beta_dct"] > 0) & (results["q_BH"] < args.q_threshold)
    print(f"   Significant positive β_DCT (q<{args.q_threshold}): {sig_pos.sum()} genes")

    # 8. DCT specificity filter
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

    # 9. Bootstrap stability
    print(f"\n9) Bootstrap stability ({args.boot_n} stratified resamples)...")

    # Determine whether OLS has enough power to be useful
    n_ols_sig = sig_pos.sum()

    if n_ols_sig > 0:
        # OLS is powered — use the strict OLS-based bootstrap
        print(f"   OLS found {n_ols_sig} significant genes → using OLS bootstrap")
        boot_freq = bootstrap_stability(
            X, Y, meta, cell_cols, dct_col,
            n_boot=args.boot_n, alpha=args.q_threshold, rng=rng
        )
        fallback_mode = False
    else:
        # OLS is under-powered (29 predictors / 80 samples) — fall back to
        # marginal-correlation bootstrap which has 1 predictor and full power.
        print(f"   OLS found 0 significant genes (model has {X.shape[1]} predictors "
              f"for {X.shape[0]} samples)")
        print(f"   → Falling back to marginal-correlation bootstrap")
        boot_freq = bootstrap_marginal_stability(
            Y, dct_vec, meta, cell_cols,
            n_boot=args.boot_n, alpha=args.q_threshold, rng=rng
        )
        fallback_mode = True

    results["bootstrap_freq"] = boot_freq
    stable = boot_freq >= args.boot_freq
    print(f"   Stable at ≥{args.boot_freq:.0%}: {stable.sum()} genes")

    #  10. Final panel

    if not fallback_mode:
        # PRIMARY path: OLS was powered
        # 1. Significant positive beta
        c1 = (results["beta_dct"] > 0) & (results["q_BH"] < args.q_threshold)
        # 2. Bootstrap stability
        c2 = results["bootstrap_freq"] >= args.boot_freq
        # 3. Sign consistency (Marginal R > 0)
        c3 = results["sign_consistent"]
        # 4. Anti-confounding (Marginal R_DCT > R_TAL + delta)
        c4 = results["anti_confounded"]
        
        panel_mask = c1 & c2 & c3 & c4
        
        print(f"\n10) Filtering stats (OLS primary path):")
        print(f"    Significant (q<{args.q_threshold}): {c1.sum()}")
        print(f"    Stable (freq>={args.boot_freq}): {c2.sum()}")
        print(f"    Sign Consistent (r>0): {c3.sum()}")
        print(f"    Anti-Confounded (δ>{args.diff_threshold}): {c4.sum()}")
    else:
        # FALLBACK path: OLS under-powered, use marginal correlation
        # BH-correct the marginal p-values
        marginal_q = bh_fdr(marginal_p)
        results["marginal_q_BH"] = marginal_q

        # f1. Marginal correlation: r > min_r AND q < threshold
        #     min_r floor eliminates weak correlations that survive BH at n=80
        f1 = (marginal_r >= args.min_r) & (marginal_q < args.q_threshold)
        # f2. Bootstrap stability (marginal-correlation bootstrap)
        f2 = results["bootstrap_freq"] >= args.boot_freq
        # f3. Anti-confounding delta - use the STRICTER fallback_delta, not diff_threshold
        f3 = results["delta_confound"] > args.fallback_delta
        # f4. OLS β_DCT > 0 (sign must agree, even if not significant)
        f4 = results["beta_dct"] > 0

        panel_mask = f1 & f2 & f3 & f4

        print(f"\n10) Filtering stats (FALLBACK: marginal correlation path):")
        print(f"    Marginal sig (r>={args.min_r}, q<{args.q_threshold}): {f1.sum()}")
        print(f"    Marginal-boot stable (freq>={args.boot_freq}): {f2.sum()}")
        print(f"    Anti-Confounded (δ>{args.fallback_delta}): {f3.sum()}")
        print(f"    OLS β_DCT > 0 (sign): {f4.sum()}")
        print(f"    Combined (before cap): {panel_mask.sum()}")
    
    if args.require_specific:
        panel_mask = panel_mask & is_specific
        spec_label = ", DCT-dominant"
    else:
        spec_label = ""

    # Panel size cap: keep top genes by delta_confound
    n_passing = panel_mask.sum()
    if n_passing > args.max_panel:
        # Rank by delta_confound descending; keep top max_panel
        passing_idx = results.index[panel_mask]
        ranked = results.loc[passing_idx, "delta_confound"].sort_values(ascending=False)
        keep_idx = ranked.index[:args.max_panel]
        panel_mask[:] = False
        panel_mask[keep_idx] = True
        print(f"\n    Panel cap: {n_passing} → {args.max_panel} (top by δ)")

    panel_genes = results.loc[panel_mask, "gene"].tolist()
    print(f"\n    FINAL DCT marker panel: {len(panel_genes)} genes")
    if fallback_mode:
        print(f"    (r≥{args.min_r}, q<{args.q_threshold}, boot≥{args.boot_freq}, "
              f"δ>{args.fallback_delta}, β>0, max={args.max_panel}{spec_label})")
    else:
        print(f"    (β>0, q<{args.q_threshold}, boot≥{args.boot_freq}, "
              f"sign-consistent, anti-conf{spec_label})")

    # Sort by beta_dct descending
    results = results.sort_values("beta_dct", ascending=False)

    # 11. Save outputs
    cols_order = ["gene", "beta_dct", "t_stat", "p_value", "q_BH",
                  "partial_r2", "dct_top3", "bootstrap_freq",
                  "pearson_r_vst_dct", "pearson_r_dct_resid",
                  "marginal_q_BH",
                  "delta_confound", "sign_consistent", "anti_confounded",
                  "pearson_r_vst_tal", "pearson_r_vst_cd",
                  "pearson_p_vst_dct", "sign_flip"]
    # Only include columns that exist
    cols_order = [c for c in cols_order if c in results.columns]
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

    # 12. Quick sanity report
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    print(f"  Total genes tested:       {len(gene_names)}")
    if fallback_mode:
        print(f"  Mode:                     FALLBACK (marginal correlation)")
        print(f"  Reason:                   OLS under-powered ({X.shape[1]} predictors / {X.shape[0]} samples)")
        print(f"  Marginal sig (r≥{args.min_r}, q<{args.q_threshold}): {f1.sum()}")
        print(f"  + Anti-confounded (δ>{args.fallback_delta}): {(f1 & f3).sum()}")
        print(f"  + Bootstrap stable:       {(f1 & f3 & f2).sum()}")
        print(f"  + OLS β>0 sign:           {(f1 & f2 & f3 & f4).sum()}")
        print(f"  After panel cap ({args.max_panel}): {panel_mask.sum()}")
    else:
        print(f"  Mode:                     PRIMARY (OLS)")
        print(f"  Significant (β>0, q<{args.q_threshold}): {sig_pos.sum()}")
        print(f"  + DCT top-3:              {sig_specific.sum()}")
        print(f"  + Bootstrap stable:       {panel_mask.sum()}")
    n_flips = results.loc[panel_mask, "sign_flip"].sum()
    print(f"  Sign flips (β>0, r<0):    {n_flips}/{panel_mask.sum()}")
    print(f"\n  Top 10 DCT markers by {'marginal r' if fallback_mode else 'β_DCT'}:")
    if fallback_mode:
        top = results.loc[panel_mask].sort_values("pearson_r_vst_dct", ascending=False).head(10)
    else:
        top = results.loc[panel_mask].head(10)
    for _, row in top.iterrows():
        print(f"    {row['gene']:25s}  β={row['beta_dct']:+.4f}  "
              f"q={row['q_BH']:.2e}  r_DCT={row.get('pearson_r_vst_dct', float('nan')):.4f}  "
              f"δ={row.get('delta_confound', float('nan')):+.4f}  "
              f"boot={row['bootstrap_freq']:.2f}")

    print(f"\n[OK] Phase 1.5 complete → {outdir}")


if __name__ == "__main__":
    main()
