# scripts/compute_phase2_lioness_on_skeleton.py
"""
Phase 2 Step A1: Compute LIONESS-style sample-specific edge weights on skeleton E

Uses the ORIGINAL (non-cell-standardized) Rtech for correlation computation.
Computes LIONESS in Fisher-z space for regression-friendly output.

Key points:
- Uses original Rtech (cell-standardization was only for topology selection)
- Correlations are clipped aggressively before atanh (CLIP_R=0.9995) to prevent explosion
- LIONESS formula applied in z-space: z_s = N*z_all - (N-1)*z_loo
- Final LIONESS z is clipped at ±ZCAP to catch remaining outliers
- Output is bounded and regression-friendly

Usage:
    python scripts/compute_phase2_lioness_on_skeleton.py
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd


# ===============================================================
# Numerical stabilization constants
# ===============================================================
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
    """Convert Pearson r to Fisher z with robust clipping.
    
    Clips |r| to CLIP_R before atanh to prevent explosion when r≈±1.
    This is critical for LIONESS where LOO can produce extreme correlations.
    """
    r = np.clip(r, -CLIP_R, CLIP_R)
    return np.arctanh(r)


def main():
    ap = argparse.ArgumentParser(description="Compute LIONESS Fisher-z weights on skeleton E")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz",
                    help="Path to Rtech.tsv.gz (genes x samples)")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz",
                    help="Path to meta_phase1.tsv.gz")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2",
                    help="Phase 2 directory with skeleton_edges.tsv")
    ap.add_argument("--out", default="lioness_z_edges.npy",
                    help="Output filename for LIONESS weights")
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
    sample_col = None
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            sample_col = col
            break
    if sample_col is None:
        sample_col = meta.columns[0]
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
        print(f"  WARNING: {len(missing)} genes not in Rtech!")
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

    # Precompute gene sums for efficient Pearson
    Sx = X.sum(axis=1)              # (G,)
    Sxx = (X**2).sum(axis=1)        # (G,)

    # Edge cross-sums
    Xi = X[ii, :]                   # (E x N)
    Xj = X[jj, :]
    Sxy = (Xi * Xj).sum(axis=1)     # (E,)

    # Pooled correlations (all samples)
    print("\n  Computing pooled correlations...")
    r_all = pearson_from_sums(N, Sx[ii], Sx[jj], Sxx[ii], Sxx[jj], Sxy)
    z_all = fisher_z_from_r(r_all)  # Fisher z with CLIP_R protection
    
    print(f"  Pooled r range: [{r_all.min():.4f}, {r_all.max():.4f}]")
    print(f"  Pooled z range: [{z_all.min():.4f}, {z_all.max():.4f}]")
    print(f"  Using CLIP_R={CLIP_R}, ZCAP={ZCAP}")

    # LIONESS z output
    out = np.empty((N, E), dtype=np.float32)
    N1 = N - 1

    print(f"\n  Computing leave-one-out LIONESS for {N} samples...")
    for s in range(N):
        if s % 20 == 0:
            print(f"    Sample {s+1}/{N}...")
        
        # Leave-one-out sums (memory-efficient: use X[ii,s] instead of Xi[:,s])
        Sx_i = Sx[ii] - X[ii, s]
        Sx_j = Sx[jj] - X[jj, s]
        Sxx_i = Sxx[ii] - X[ii, s]**2
        Sxx_j = Sxx[jj] - X[jj, s]**2
        Sxy_loo = Sxy - (X[ii, s] * X[jj, s])

        r_loo = pearson_from_sums(N1, Sx_i, Sx_j, Sxx_i, Sxx_j, Sxy_loo)
        z_loo = fisher_z_from_r(r_loo)  # Fisher z with CLIP_R protection

        # LIONESS in Fisher space
        z_s = N * z_all - (N - 1) * z_loo
        
        # Final clip to catch any remaining pathological values
        z_s = np.clip(z_s, -ZCAP, ZCAP)
        out[s, :] = z_s.astype(np.float32)

    # Save outputs
    np.save(p2 / args.out, out)
    edges.to_csv(p2 / "edge_index.tsv", sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print("LIONESS computation complete!")
    print(f"{'=' * 60}")
    print(f"\nOutput shape: {out.shape} (samples × edges)")
    print(f"  Fisher-z range: [{out.min():.4f}, {out.max():.4f}]")
    print(f"  Mean |z|: {np.abs(out).mean():.4f}")
    print(f"\nOutputs in {p2}:")
    print(f"  - {args.out}")
    print(f"  - edge_index.tsv")


if __name__ == "__main__":
    main()
