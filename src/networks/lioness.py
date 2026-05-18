# src/networks/lioness.py
"""
Phase 2 Step A1: Compute LIONESS-style sample-specific edge weights on skeleton E.

Uses the residualized, non-cell-standardized Rtech for sample-specific edge
weights. The primary output is a raw LIONESS correlation contribution that is
rank-normalized edge-wise across samples before regression. Fisher z is applied
only to pooled/leave-one-out correlation estimates in the explicit
``z_contribution`` sensitivity mode; sample-specific LIONESS weights are not
treated as correlations and are not Fisher transformed in the default mode.

Usage:
    python -m src.networks.lioness --lioness-transform raw_ranknorm
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata

from src.common import find_sample_col


# Numerical stabilization constants
CLIP_R = 0.9995   # Clip correlations before atanh to prevent explosion from r≈±1
ZCAP = 20.0       # Final cap on LIONESS z to catch outlier pathology


def pearson_from_sums(N, Sx, Sy, Sxx, Syy, Sxy, eps=1e-12):
    """Compute Pearson correlation from precomputed sums (safe version)."""
    num = N * Sxy - Sx * Sy
    # Fix: sqrt can fail if N*Sxx - Sx^2 goes slightly negative due to rounding
    denx = N * Sxx - Sx * Sx
    deny = N * Syy - Sy * Sy
    den = np.sqrt(np.maximum(denx, 0.0) * np.maximum(deny, 0.0))
    r = np.where(den > eps, num / den, 0.0)
    return np.clip(r, -1.0, 1.0)


def fisher_z_from_r(r: np.ndarray) -> np.ndarray:
    """Convert true Pearson correlations to Fisher z with robust clipping."""
    r = np.clip(r, -CLIP_R, CLIP_R)
    return np.arctanh(r)


def rank_normalize_edges(weights: np.ndarray) -> np.ndarray:
    """Rank-normalize each edge across samples to normal scores.

    This is the default regression scale for LIONESS contributions because the
    sample-specific values are linear contributions, not bounded correlations.
    """
    weights = np.asarray(weights, dtype=np.float64)
    n, e = weights.shape
    if n < 2:
        return np.zeros_like(weights, dtype=np.float32)
    out = np.empty_like(weights, dtype=np.float64)
    for col in range(e):
        ranks = rankdata(weights[:, col], method="average")
        q = (ranks - 0.5) / n
        out[:, col] = norm.ppf(q)
    return out.astype(np.float32)


def robust_scale_edges(weights: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Robust-scale each edge across samples using median/MAD with SD fallback."""
    weights = np.asarray(weights, dtype=np.float64)
    med = np.nanmedian(weights, axis=0, keepdims=True)
    mad = np.nanmedian(np.abs(weights - med), axis=0, keepdims=True) * 1.4826
    sd = np.nanstd(weights, axis=0, keepdims=True)
    scale = np.where(mad > eps, mad, sd)
    out = np.divide(weights - med, scale + eps, out=np.zeros_like(weights), where=scale > eps)
    return out.astype(np.float32)


def normalize_lioness_weights(weights: np.ndarray, transform: str) -> np.ndarray:
    """Apply the requested edge-wise regression transform."""
    if transform == "raw":
        return np.asarray(weights, dtype=np.float32)
    if transform == "raw_ranknorm":
        return rank_normalize_edges(weights)
    if transform == "raw_robust":
        return robust_scale_edges(weights)
    raise ValueError(f"Unsupported normalization transform for raw weights: {transform}")


def compute_lioness_weights(
    X: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    transform: str = "raw_ranknorm",
) -> tuple[np.ndarray, dict[str, object]]:
    """Compute LIONESS edge weights for X (genes x samples).

    ``raw*`` modes apply the LIONESS linear formula on Pearson correlations and
    then optionally normalize sample-specific contributions. ``z_contribution``
    applies the linear formula after Fisher-z-transforming the pooled and LOO
    correlations; it is retained only as a sensitivity analysis.
    """
    G, N = X.shape
    E = len(ii)
    if N < 3:
        raise ValueError("Need at least 3 samples for LIONESS leave-one-out correlations")

    Sx = X.sum(axis=1)
    Sxx = (X**2).sum(axis=1)
    Xi = X[ii, :]
    Xj = X[jj, :]
    Sxy = (Xi * Xj).sum(axis=1)

    r_all = pearson_from_sums(N, Sx[ii], Sx[jj], Sxx[ii], Sxx[jj], Sxy)
    if transform == "z_contribution":
        base_all = fisher_z_from_r(r_all)
    elif transform in {"raw", "raw_ranknorm", "raw_robust"}:
        base_all = r_all
    else:
        raise ValueError(
            f"Unknown --lioness-transform={transform}. "
            "Use raw_ranknorm, raw_robust, raw, or z_contribution."
        )

    raw_out = np.empty((N, E), dtype=np.float64)
    N1 = N - 1
    for s in range(N):
        Sx_i = Sx[ii] - X[ii, s]
        Sx_j = Sx[jj] - X[jj, s]
        Sxx_i = Sxx[ii] - X[ii, s] ** 2
        Sxx_j = Sxx[jj] - X[jj, s] ** 2
        Sxy_loo = Sxy - (X[ii, s] * X[jj, s])

        r_loo = pearson_from_sums(N1, Sx_i, Sx_j, Sxx_i, Sxx_j, Sxy_loo)
        if transform == "z_contribution":
            loo = fisher_z_from_r(r_loo)
        else:
            loo = r_loo

        contribution = N * base_all - (N - 1) * loo
        if transform == "z_contribution":
            contribution = np.clip(contribution, -ZCAP, ZCAP)
        raw_out[s, :] = contribution

    if transform == "z_contribution":
        out = raw_out.astype(np.float32)
    else:
        out = normalize_lioness_weights(raw_out, transform)

    meta = {
        "transform": transform,
        "n_genes": int(G),
        "n_samples": int(N),
        "n_edges": int(E),
        "pooled_r_min": float(np.nanmin(r_all)),
        "pooled_r_max": float(np.nanmax(r_all)),
        "raw_lioness_min": float(np.nanmin(raw_out)),
        "raw_lioness_max": float(np.nanmax(raw_out)),
        "output_min": float(np.nanmin(out)),
        "output_max": float(np.nanmax(out)),
        "sample_specific_values_are_correlations": False,
        "fisher_z_used_on_sample_specific_weights": transform == "z_contribution",
    }
    return out, meta


def main():
    ap = argparse.ArgumentParser(description="Compute LIONESS sample-specific edge weights on skeleton E")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz",
                    help="Path to Rtech.tsv.gz (genes x samples)")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz",
                    help="Path to meta_phase1.tsv.gz")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2",
                    help="Phase 2 directory with skeleton_edges.tsv")
    ap.add_argument("--out", default="lioness_edges.npy",
                    help="Output filename for LIONESS weights")
    ap.add_argument("--lioness-transform", default="raw_ranknorm",
                    choices=["raw_ranknorm", "raw_robust", "raw", "z_contribution"],
                    help="Default raw_ranknorm avoids Fisher z on sample-specific LIONESS weights. "
                         "z_contribution is a sensitivity mode only.")
    args = ap.parse_args()

    p2 = Path(args.phase2_dir)

    print("=" * 60)
    print("Phase 2 Step A1: LIONESS Sample-Specific Weights on E")
    print("=" * 60)

    # Load skeleton
    print(f"\nLoading Phase 2 skeleton from {p2}")
    genes = (p2 / "phase2_genes.txt").read_text().splitlines()
    genes = [g for g in genes if g.strip()]
    edges = pd.read_csv(p2 / "skeleton_edges.tsv", sep="\t")
    print(f"  Genes: {len(genes)}")
    print(f"  Edges: {len(edges)}")

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    
    # Validate gene mapping (fail hard on bad edges)
    ii = edges["gene_i"].map(gene_to_idx)
    jj = edges["gene_j"].map(gene_to_idx)
    bad = ii.isna() | jj.isna()
    if bad.any():
        bad_edges = edges.loc[bad, ["gene_i", "gene_j"]].head(20)
        raise ValueError(f"Skeleton edges reference genes not in phase2_genes.txt. Examples:\n{bad_edges}")
    ii = ii.to_numpy(np.int32)
    jj = jj.to_numpy(np.int32)

    # Load ORIGINAL Rtech (not cell-standardized!)
    print(f"\nLoading Rtech: {args.rtech}")
    rtech = pd.read_csv(args.rtech, sep="\t", compression="gzip", index_col=0)
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")

    # Find sample column
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)

    # Align
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]
    print(f"  Aligned: {len(common)} samples")
    
    # Save sample order for downstream regression alignment
    (p2 / "lioness_samples.txt").write_text("\n".join(common) + "\n")
    print(f"  → Saved sample order: lioness_samples.txt")

    # Extract gene matrix - only Phase 2 genes
    missing = [g for g in genes if g not in rtech.index]
    if missing:
        print(f"  WARNING: {len(missing)} genes not in Rtech")
        print(f"  First 20 missing: {missing[:20]}")
        # Keep only genes that exist
        genes = [g for g in genes if g in rtech.index]
        # Rebuild edge index for remaining genes with validation
        gene_to_idx = {g: i for i, g in enumerate(genes)}
        mask = edges["gene_i"].isin(genes) & edges["gene_j"].isin(genes)
        edges = edges[mask].reset_index(drop=True)
        ii = edges["gene_i"].map(gene_to_idx)
        jj = edges["gene_j"].map(gene_to_idx)
        bad = ii.isna() | jj.isna()
        if bad.any():
            raise ValueError("Gene mismatch after filtering - this should not happen")
        ii = ii.to_numpy(np.int32)
        jj = jj.to_numpy(np.int32)
        # Save effective gene list for reproducibility
        (p2 / "phase2_genes_effective.txt").write_text("\n".join(genes) + "\n")
        print(f"  Saved effective gene list: phase2_genes_effective.txt")

    X = rtech.loc[genes].values.astype(np.float64)  # (G x N)
    G, N = X.shape
    E = len(ii)
    print(f"  Expression matrix: {G} genes × {N} samples")
    print(f"  Computing LIONESS for {E} edges...")
    print(f"  Transform: {args.lioness_transform}")
    if args.lioness_transform == "z_contribution":
        print("  NOTE: z_contribution is a sensitivity mode; outputs are transformed contributions, not correlations.")

    out, transform_meta = compute_lioness_weights(X, ii, jj, transform=args.lioness_transform)

    # Save outputs
    np.save(p2 / args.out, out)
    if args.out != "lioness_z_edges.npy":
        # Compatibility sentinel for older scripts: do not duplicate data unless requested.
        legacy_note = {
            "current_edge_weight_file": args.out,
            "legacy_lioness_z_edges": "not written by default; pass --out=lioness_z_edges.npy for legacy paths",
        }
        (p2 / "lioness_legacy_note.json").write_text(json.dumps(legacy_note, indent=2) + "\n")
    edges.to_csv(p2 / "edge_index.tsv", sep="\t", index=False)
    (p2 / "lioness_transform_metadata.json").write_text(json.dumps(transform_meta, indent=2) + "\n")

    print(f"\n{'=' * 60}")
    print("LIONESS computation complete")
    print(f"{'=' * 60}")
    print(f"\nOutput shape: {out.shape} (samples × edges)")
    print(f"  Output range: [{out.min():.4f}, {out.max():.4f}]")
    print(f"  Mean |weight|: {np.abs(out).mean():.4f}")
    print(f"\nOutputs in {p2}:")
    print(f"  - {args.out}")
    print(f"  - edge_index.tsv")
    print(f"  - lioness_transform_metadata.json")


if __name__ == "__main__":
    main()
