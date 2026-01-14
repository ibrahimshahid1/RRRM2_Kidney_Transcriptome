# scripts/build_phase2_skeleton.py
"""
Phase 2 Step 2.1-2.3: Build Shared Sparse Skeleton E

Cell-standardize expression within each (Age × Arm × EnvGroup) cell,
then build a sparse partial correlation network using Ledoit-Wolf
shrinkage and top-k neighbors per gene.

The skeleton E is fixed for all downstream sample-specific weighting.

Usage:
    python scripts/build_phase2_skeleton.py --max_genes 1200 --topk 80
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf


def load_rtech(path: str) -> pd.DataFrame:
    """Load Rtech matrix (genes x samples) from gzipped TSV."""
    return pd.read_csv(path, sep="\t", compression="gzip", index_col=0)


def load_meta(path: str) -> pd.DataFrame:
    """Load metadata from gzipped TSV."""
    return pd.read_csv(path, sep="\t", compression="gzip")


def pick_genes(rtech: pd.DataFrame, max_genes: int) -> list[str]:
    """Select top genes by variance, excluding ERCC."""
    keep = ~rtech.index.str.upper().str.startswith("ERCC")
    r = rtech.loc[keep]
    v = r.var(axis=1)
    return v.sort_values(ascending=False).head(max_genes).index.tolist()


def cell_standardize(
    rtech_gxs: pd.DataFrame,
    meta: pd.DataFrame,
    cell_cols: list[str],
    eps: float = 1e-8,
    sd_floor: float = 1e-3,
) -> np.ndarray:
    """
    Standardize within each experimental cell (defined by cell_cols).
    
    Args:
        rtech_gxs: genes x samples DataFrame
        meta: metadata with index = sample IDs
        cell_cols: columns defining experimental cells
        eps: small constant for numerical stability
        sd_floor: minimum SD to avoid division issues
    
    Returns:
        Z: (n_samples x n_genes) cell-standardized matrix
    """
    samples = rtech_gxs.columns.tolist()
    meta_aligned = meta.loc[samples]
    cell_key = meta_aligned[cell_cols].astype(str).agg("|".join, axis=1)
    
    # Work in samples x genes for vectorized operations
    X = rtech_gxs.T.values.astype(np.float64)  # (N x G)
    Z = np.empty_like(X)
    
    for ck in cell_key.unique():
        idx = np.where(cell_key.values == ck)[0]
        Xc = X[idx, :]  # (n_cell x G)
        mu = Xc.mean(axis=0)
        # Fix: ddof=1 blows up for n=1 cells
        n_cell = Xc.shape[0]
        ddof = 1 if n_cell >= 2 else 0
        sd = Xc.std(axis=0, ddof=ddof)
        sd = np.where(np.isfinite(sd), sd, 0.0)  # NaN guard
        sd = np.maximum(sd, sd_floor)
        Z[idx, :] = (Xc - mu) / (sd + eps)
    
    return Z  # (N x G)


def partial_corr_from_precision(P: np.ndarray) -> np.ndarray:
    """Convert precision matrix to partial correlations."""
    d = np.sqrt(np.diag(P))
    pc = -P / np.outer(d, d)
    np.fill_diagonal(pc, 0.0)
    return pc


def topk_skeleton(pc: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Build skeleton by taking top-k neighbors per gene.
    
    The union of all top-k neighbors gives ~G*k edges (not G*k/2).
    """
    G = pc.shape[0]
    abs_pc = np.abs(pc)
    # Fix: guard against k > G-1
    k_eff = min(k, G - 1)
    
    edges = set()
    for i in range(G):
        # Fix: exclude self from consideration
        row = abs_pc[i].copy()
        row[i] = -np.inf
        # Use argpartition for O(n) instead of O(n log n)
        idx = np.argpartition(row, -k_eff)[-k_eff:]
        for j in idx:
            # Store as (min, max) for deduplication
            a, b = (i, j) if i < j else (j, i)
            edges.add((a, b))
    
    ii = np.fromiter((e[0] for e in edges), dtype=np.int32)
    jj = np.fromiter((e[1] for e in edges), dtype=np.int32)
    return ii, jj


def main():
    ap = argparse.ArgumentParser(description="Build Phase 2 shared skeleton E")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz",
                    help="Path to Rtech.tsv.gz (genes x samples)")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz",
                    help="Path to meta_phase1.tsv.gz")
    ap.add_argument("--outdir", default="data/processed/networks/phase2",
                    help="Output directory")
    ap.add_argument("--max_genes", type=int, default=1200,
                    help="Maximum genes for skeleton")
    ap.add_argument("--cell_cols", default="Age,Arm,EnvGroup",
                    help="Comma-separated columns defining experimental cells")
    ap.add_argument("--topk", type=int, default=80,
                    help="Top-k neighbors per gene (~G*k edges)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 2: Build Shared Sparse Skeleton E")
    print("=" * 60)

    # Load data
    print(f"\nLoading Rtech: {args.rtech}")
    rtech = load_rtech(args.rtech)
    print(f"  Shape: {rtech.shape[0]} genes × {rtech.shape[1]} samples")

    print(f"Loading metadata: {args.meta}")
    meta = load_meta(args.meta)

    # Find sample column
    sample_col = None
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            sample_col = col
            break
    if sample_col is None:
        sample_col = meta.columns[0]
        print(f"  Warning: using first column '{sample_col}' as sample ID")
    meta = meta.set_index(sample_col, drop=False)

    # Align
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]
    print(f"  Aligned: {len(common)} samples")

    # Gene selection
    genes = pick_genes(rtech, args.max_genes)
    (outdir / "phase2_genes.txt").write_text("\n".join(genes) + "\n")
    print(f"\nSelected {len(genes)} top-variance genes")
    print(f"  → Saved to {outdir / 'phase2_genes.txt'}")

    rtech_gxs = rtech.loc[genes]
    cell_cols = [c.strip() for c in args.cell_cols.split(",") if c.strip()]
    print(f"\nCell columns: {cell_cols}")

    # Step 2.1: Cell-standardize for topology selection
    print("\nStep 2.1: Cell-standardizing within each experimental cell...")
    n_cells = meta[cell_cols].astype(str).agg("|".join, axis=1).nunique()
    print(f"  Found {n_cells} experimental cells (n≈{len(common)//n_cells} per cell)")
    Z = cell_standardize(rtech_gxs, meta, cell_cols=cell_cols)
    print(f"  Cell-standardized matrix: {Z.shape[0]} samples × {Z.shape[1]} genes")

    # Step 2.3: Shrinkage partial correlation
    print("\nStep 2.3: Computing shrinkage partial correlations (Ledoit-Wolf)...")
    lw = LedoitWolf().fit(Z)
    cov = lw.covariance_
    print(f"  Shrinkage: {lw.shrinkage_:.4f}")
    
    # Invert covariance to get precision (with ridge fallback)
    try:
        prec = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        ridge = 1e-6 * np.eye(cov.shape[0])
        prec = np.linalg.inv(cov + ridge)
        print("  Warning: added tiny ridge for inversion")
    pc = partial_corr_from_precision(prec)
    
    # Build skeleton
    print(f"\nBuilding skeleton with top-k={args.topk} neighbors per gene...")
    ii, jj = topk_skeleton(pc, k=args.topk)
    
    # Save edge indices for downstream determinism
    np.save(outdir / "edge_i.npy", ii)
    np.save(outdir / "edge_j.npy", jj)
    
    # Save edge dataframe with partial correlation weights
    w = pc[ii, jj]
    edge_df = pd.DataFrame({
        "gene_i": [genes[i] for i in ii],
        "gene_j": [genes[j] for j in jj],
        "pcorr": w,
        "abs_pcorr": np.abs(w),
        "i": ii,
        "j": jj,
    })
    edge_df.to_csv(outdir / "skeleton_edges.tsv", sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print(f"Skeleton E built successfully!")
    print(f"{'=' * 60}")
    print(f"  Genes: {len(genes)}")
    print(f"  Edges: {len(edge_df)} (target: ~{len(genes)*args.topk//2} to ~{len(genes)*args.topk})")
    print(f"\nOutputs in {outdir}:")
    print(f"  - phase2_genes.txt")
    print(f"  - skeleton_edges.tsv (with pcorr weights)")
    print(f"  - edge_i.npy, edge_j.npy (indices)")


if __name__ == "__main__":
    main()
