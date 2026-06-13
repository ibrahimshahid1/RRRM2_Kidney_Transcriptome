#!/usr/bin/env python3
"""Phase 1.5b: Generalized Segment Marker Discovery"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from src.common import REPO_ROOT, find_sample_col, bh_fdr as _bh_fdr


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """BH-FDR with NaN/Inf guard."""
    p = np.asarray(p, dtype=float)
    p = np.where(np.isfinite(p), p, 1.0)
    return _bh_fdr(p)


# Core regression (same as discover_dct.py, parameterized)

def build_design_matrix(
    meta: pd.DataFrame,
    clr: pd.DataFrame,
    target_segment: str,
    cell_cols: list[str],
    drop_segment: str | None = None,
) -> tuple[np.ndarray, int, list[str]]:
    """Build design matrix with target segment as the coefficient of interest.

    Returns (X, target_col_idx, column_names).
    """
    n = len(meta)
    parts = [np.ones((n, 1))]
    col_names = ["intercept"]

    # Target segment CLR — the coefficient of interest
    target_col_idx = 1
    target_vals = clr[[target_segment]].values.astype(float)
    t_mean, t_std = np.nanmean(target_vals), np.nanstd(target_vals)
    if t_std > 1e-8:
        target_vals = (target_vals - t_mean) / t_std
    parts.append(target_vals)
    col_names.append(f"CLR_{target_segment}")

    # Other segments (drop one for compositionality)
    if drop_segment is None:
        # Auto-pick: drop the largest segment (typically PT)
        drop_segment = clr.drop(columns=[target_segment]).var().idxmax()

    other_cols = [c for c in clr.columns
                  if c != target_segment and c != drop_segment]
    if other_cols:
        other_vals = clr[other_cols].values.astype(float)
        o_mean = np.nanmean(other_vals, axis=0)
        o_std = np.nanstd(other_vals, axis=0)
        o_std[o_std < 1e-8] = 1.0
        other_vals = (other_vals - o_mean) / o_std
        parts.append(other_vals)
        col_names.extend([f"CLR_{c}" for c in other_cols])

    # Cell fixed effects
    cell_labels = meta[cell_cols].astype(str).agg("|".join, axis=1)
    cell_dummies = pd.get_dummies(cell_labels, drop_first=True, dtype=float)
    parts.append(cell_dummies.values)
    col_names.extend(cell_dummies.columns.tolist())

    X = np.hstack(parts).astype(np.float64)
    return X, target_col_idx, col_names


def ols_all_genes(X: np.ndarray, Y: np.ndarray, target_col: int) -> pd.DataFrame:
    """Vectorised OLS for all genes. Returns beta, t-stat, p-value, partial R²."""
    n, p = X.shape
    rank = np.linalg.matrix_rank(X)

    X_pinv = np.linalg.pinv(X)
    B = X_pinv @ Y
    residuals = Y - X @ B
    rss = np.sum(residuals ** 2, axis=0)
    df_resid = max(n - rank, 1)
    sigma2 = rss / df_resid

    XtX_inv = np.linalg.pinv(X.T @ X)
    var_beta = XtX_inv[target_col, target_col] * sigma2

    beta = B[target_col, :]
    se = np.sqrt(np.maximum(var_beta, 1e-30))
    t_stat = beta / se
    p_value = 2.0 * sp_stats.t.sf(np.abs(t_stat), df_resid)

    # Partial R²
    X_red = np.delete(X, target_col, axis=1)
    B_red = np.linalg.pinv(X_red) @ Y
    rss_red = np.sum((Y - X_red @ B_red) ** 2, axis=0)
    partial_r2 = np.clip(1.0 - rss / np.maximum(rss_red, 1e-30), 0, 1)

    return pd.DataFrame({
        "beta": beta,
        "t_stat": t_stat,
        "p_value": p_value,
        "partial_r2": partial_r2,
    })


def bootstrap_stability(
    X: np.ndarray, Y: np.ndarray, meta: pd.DataFrame,
    cell_cols: list[str], target_col: int,
    n_boot: int, alpha: float, rng: np.random.Generator,
) -> np.ndarray:
    """Stratified bootstrap. Returns (G,) frequency of passing β>0 & q<α."""
    G = Y.shape[1]
    pass_counts = np.zeros(G, dtype=int)

    cell_labels = meta[cell_cols].astype(str).agg("|".join, axis=1).values
    cell_indices = {c: np.where(cell_labels == c)[0]
                    for c in np.unique(cell_labels)}

    for _ in range(n_boot):
        boot_idx = np.concatenate([
            rng.choice(idx, size=len(idx), replace=True)
            for idx in cell_indices.values()
        ])
        X_b, Y_b = X[boot_idx], Y[boot_idx]
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
        var_b = XtX_inv_b[target_col, target_col] * sigma2_b
        beta_b = B_b[target_col, :]
        t_b = beta_b / np.sqrt(np.maximum(var_b, 1e-30))
        p_b = 2.0 * sp_stats.t.sf(np.abs(t_b), df_b)
        q_b = bh_fdr(p_b)
        pass_counts += ((beta_b > 0) & (q_b < alpha)).astype(int)

    return pass_counts / n_boot


def bootstrap_marginal(
    Y: np.ndarray, target_vec: np.ndarray,
    meta: pd.DataFrame, cell_cols: list[str],
    n_boot: int, alpha: float, rng: np.random.Generator,
) -> np.ndarray:
    """Fallback: marginal correlation bootstrap."""
    G = Y.shape[1]
    pass_counts = np.zeros(G, dtype=int)

    cell_labels = meta[cell_cols].astype(str).agg("|".join, axis=1).values
    cell_indices = {c: np.where(cell_labels == c)[0]
                    for c in np.unique(cell_labels)}

    for _ in range(n_boot):
        boot_idx = np.concatenate([
            rng.choice(idx, size=len(idx), replace=True)
            for idx in cell_indices.values()
        ])
        t_b = target_vec[boot_idx]
        Y_b = Y[boot_idx]
        n = len(boot_idx)

        t_c = t_b - t_b.mean()
        Y_c = Y_b - Y_b.mean(axis=0)
        t_ss = np.dot(t_c, t_c)
        Y_ss = np.sum(Y_c ** 2, axis=0)
        r_b = (t_c @ Y_c) / np.sqrt(np.maximum(t_ss * Y_ss, 1e-30))
        t_stat = r_b * np.sqrt((n - 2) / np.maximum(1 - r_b ** 2, 1e-30))
        p_b = 2.0 * sp_stats.t.sf(np.abs(t_stat), n - 2)
        q_b = bh_fdr(p_b)
        pass_counts += ((r_b > 0) & (q_b < alpha)).astype(int)

    return pass_counts / n_boot


# Single-segment discovery

def discover_segment_markers(
    vst: pd.DataFrame,
    meta: pd.DataFrame,
    clr: pd.DataFrame,
    target_segment: str,
    cell_cols: list[str],
    drop_segment: str | None = None,
    q_threshold: float = 0.05,
    boot_n: int = 200,
    boot_freq: float = 0.70,
    min_r: float = 0.25,
    fallback_delta: float = 0.05,
    max_panel: int = 200,
    seed: int = 42,
) -> tuple[pd.DataFrame, list[str]]:
    """Discover marker genes for a single segment.

    Returns (results_df, panel_genes).
    """
    from scipy.stats import pearsonr as _pearsonr

    rng = np.random.default_rng(seed)

    print(f"\n{'=' * 60}")
    print(f"Discovering markers for: {target_segment}")
    print(f"{'=' * 60}")

    if target_segment not in clr.columns:
        print(f"  WARNING: {target_segment} not in CLR columns. Skipping.")
        return pd.DataFrame(), []

    # Check variance
    seg_std = clr[target_segment].std()
    if seg_std < 1e-6:
        print(f"  WARNING: {target_segment} has near-zero variance (std={seg_std:.2e}). Skipping.")
        return pd.DataFrame(), []

    # Build design matrix
    X, target_col, _ = build_design_matrix(
        meta, clr, target_segment, cell_cols, drop_segment
    )
    rank = np.linalg.matrix_rank(X)
    print(f"  Design: {X.shape[0]} samples × {X.shape[1]} predictors (rank={rank})")

    # OLS
    Y = vst.values.T  # (samples, genes)
    gene_names = list(vst.index)
    G = len(gene_names)

    results = ols_all_genes(X, Y, target_col)
    results.index = gene_names
    results["gene"] = gene_names
    results["q_BH"] = bh_fdr(results["p_value"].values)

    # Marginal correlation
    target_vec = clr[target_segment].values
    marginal_r = np.array([_pearsonr(Y[:, g], target_vec)[0] for g in range(G)])
    marginal_p = np.array([_pearsonr(Y[:, g], target_vec)[1] for g in range(G)])
    results["pearson_r"] = marginal_r
    results["pearson_p"] = marginal_p
    results["sign_consistent"] = (results["beta"] > 0) & (marginal_r > 0)

    # Anti-confounding: compute correlation with all OTHER segments
    other_segs = [c for c in clr.columns if c != target_segment]
    confound_r = {}
    for seg in other_segs:
        confound_r[seg] = np.array([_pearsonr(Y[:, g], clr[seg].values)[0]
                                     for g in range(G)])

    if confound_r:
        confound_stack = np.column_stack(list(confound_r.values()))
        max_confound_r = np.max(np.abs(confound_stack), axis=1)
        delta = marginal_r - max_confound_r
        results["delta_confound"] = delta
        results["anti_confounded"] = delta > 0
    else:
        results["delta_confound"] = 1.0
        results["anti_confounded"] = True

    # Determine OLS vs fallback path
    sig_pos = (results["beta"] > 0) & (results["q_BH"] < q_threshold)
    n_ols_sig = sig_pos.sum()
    print(f"  OLS significant (β>0, q<{q_threshold}): {n_ols_sig}")

    # Bootstrap
    if n_ols_sig > 0:
        boot = bootstrap_stability(
            X, Y, meta, cell_cols, target_col,
            n_boot=boot_n, alpha=q_threshold, rng=rng
        )
        fallback = False
    else:
        print(f"  OLS underpowered → falling back to marginal correlation")
        boot = bootstrap_marginal(
            Y, target_vec, meta, cell_cols,
            n_boot=boot_n, alpha=q_threshold, rng=rng
        )
        fallback = True

    results["bootstrap_freq"] = boot

    # Panel selection
    if not fallback:
        panel_mask = (
            sig_pos &
            (boot >= boot_freq) &
            results["sign_consistent"] &
            results["anti_confounded"]
        )
    else:
        marginal_q = bh_fdr(marginal_p)
        results["marginal_q_BH"] = marginal_q
        panel_mask = (
            (marginal_r >= min_r) &
            (marginal_q < q_threshold) &
            (boot >= boot_freq) &
            (results["delta_confound"] > fallback_delta) &
            (results["beta"] > 0)
        )

    # Cap panel size
    n_passing = panel_mask.sum()
    if n_passing > max_panel:
        passing_genes = results.index[panel_mask]
        ranked = results.loc[passing_genes, "delta_confound"].sort_values(ascending=False)
        keep = ranked.index[:max_panel]
        panel_mask = results.index.isin(keep)
        print(f"  Panel capped: {n_passing} → {max_panel}")

    panel_genes = results.loc[panel_mask, "gene"].tolist()
    mode = "FALLBACK (marginal)" if fallback else "OLS"
    print(f"  {target_segment} panel: {len(panel_genes)} markers ({mode})")

    # Top 5
    if panel_genes:
        top = results.loc[panel_mask].sort_values(
            "pearson_r" if fallback else "beta", ascending=False
        ).head(5)
        for _, row in top.iterrows():
            print(f"    {row['gene']:20s}  β={row['beta']:+.4f}  r={row['pearson_r']:.4f}  "
                  f"δ={row['delta_confound']:+.4f}  boot={row['bootstrap_freq']:.2f}")

    results = results.sort_values("beta", ascending=False)
    return results, panel_genes


# Main: run for all segments

def main():
    ap = argparse.ArgumentParser(
        description="Discover marker genes for all kidney segments"
    )
    ap.add_argument("--vst",
                    default=str(REPO_ROOT / "data/processed/vst_normalized"
                                / "GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"),
                    help="Path to VST expression CSV (genes × samples)")
    ap.add_argument("--meta",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals"
                                / "meta_phase1.tsv.gz"),
                    help="Path to metadata TSV")
    ap.add_argument("--clr",
                    default=str(REPO_ROOT / "data/processed/deconvolution/latest"
                                / "music_segment_direct_proportions_CLR.csv"),
                    help="Path to CLR proportions CSV")
    ap.add_argument("--segments", nargs="*", default=None,
                    help="Specific segments to discover (default: all in CLR file)")
    ap.add_argument("--cell_cols", default="Age,Arm,EnvGroup")
    ap.add_argument("--drop_segment", default=None,
                    help="Segment to drop for compositionality (default: auto = highest variance)")
    ap.add_argument("--q_threshold", type=float, default=0.05)
    ap.add_argument("--boot_n", type=int, default=200)
    ap.add_argument("--boot_freq", type=float, default=0.70)
    ap.add_argument("--min_r", type=float, default=0.25)
    ap.add_argument("--fallback_delta", type=float, default=0.05)
    ap.add_argument("--max_panel", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/processed/segment_markers"),
                    help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 1.5b: Segment Marker Discovery (All Segments)")
    print("=" * 60)

    # Load data once
    print(f"\nLoading VST expression: {args.vst}")
    vst = pd.read_csv(args.vst, index_col=0)
    vst.index = vst.index.str.replace(r"\.\d+$", "", regex=True)
    if vst.index.duplicated().any():
        vst = vst.groupby(vst.index).mean()
    print(f"  {vst.shape[0]} genes × {vst.shape[1]} samples")

    print(f"\nLoading metadata: {args.meta}")
    meta = pd.read_csv(args.meta, sep="\t")
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)

    print(f"\nLoading CLR: {args.clr}")
    clr = pd.read_csv(args.clr, index_col=0)
    segments = list(clr.columns)
    print(f"  Segments available: {segments}")

    # Align
    common = sorted(set(vst.columns) & set(meta.index) & set(clr.index))
    vst = vst[common]
    meta = meta.loc[common]
    clr = clr.loc[common]
    print(f"  Aligned: {len(common)} samples")

    # Filter low-variance genes
    vst = vst.loc[vst.var(axis=1) > 1e-8]
    print(f"  Genes after variance filter: {vst.shape[0]}")

    cell_cols = [c.strip() for c in args.cell_cols.split(",")]

    # Which segments to run
    if args.segments:
        targets = [s for s in args.segments if s in segments]
        missing = [s for s in args.segments if s not in segments]
        if missing:
            print(f"\n  WARNING: segments not in CLR: {missing}")
    else:
        targets = segments

    print(f"\n  Will discover markers for {len(targets)} segments: {targets}")

    # Run for each segment
    summary_rows = []
    for seg in targets:
        results, panel = discover_segment_markers(
            vst, meta, clr, seg, cell_cols,
            drop_segment=args.drop_segment,
            q_threshold=args.q_threshold,
            boot_n=args.boot_n,
            boot_freq=args.boot_freq,
            min_r=args.min_r,
            fallback_delta=args.fallback_delta,
            max_panel=args.max_panel,
            seed=args.seed,
        )

        if results.empty:
            continue

        # Save scores
        scores_path = outdir / f"{seg}_marker_scores.tsv"
        cols = [c for c in ["gene", "beta", "t_stat", "p_value", "q_BH",
                            "partial_r2", "pearson_r", "pearson_p",
                            "delta_confound", "sign_consistent",
                            "anti_confounded", "bootstrap_freq",
                            "marginal_q_BH"] if c in results.columns]
        results[cols].to_csv(scores_path, sep="\t", index=False)

        # Save panel
        panel_path = outdir / f"{seg}_marker_panel.txt"
        with open(panel_path, "w") as f:
            f.write("\n".join(sorted(panel)) + "\n")

        summary_rows.append({
            "segment": seg,
            "n_markers": len(panel),
            "scores_file": scores_path.name,
            "panel_file": panel_path.name,
        })

        print(f"  Saved: {scores_path.name} ({len(results)} genes), "
              f"{panel_path.name} ({len(panel)} markers)")

    # Summary
    summary = pd.DataFrame(summary_rows)
    summary_path = outdir / "segment_marker_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    for _, row in summary.iterrows():
        print(f"  {row['segment']:15s}  {row['n_markers']:4d} markers")
    print(f"\n[OK] All segment markers saved to: {outdir}")


if __name__ == "__main__":
    main()
