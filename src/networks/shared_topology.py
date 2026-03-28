# src/networks/shared_topology.py
"""
Phase 2 Step 2.1-2.3: Build Shared Sparse Skeleton E

Cell-standardize expression within each (Age × Arm × EnvGroup) cell,
then build a sparse partial correlation network using Ledoit-Wolf
shrinkage and top-k neighbors per gene.

The skeleton E is fixed for all downstream sample-specific weighting.

Usage:
    python -m src.networks.shared_topology --max_genes 2500 --topk 80

With biotype filtering (recommended):
    python -m src.networks.shared_topology --max_genes 2500 --topk 80 \\
        --id_map data/processed/resources/id_map.tsv \\
        --biotype_filter protein_coding
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

from src.common import REPO_ROOT
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf


def load_rtech(path: str) -> pd.DataFrame:
    """Load Rtech matrix (genes x samples) from gzipped TSV."""
    return pd.read_csv(path, sep="\t", compression="gzip", index_col=0)


def load_meta(path: str) -> pd.DataFrame:
    """Load metadata from gzipped TSV."""
    return pd.read_csv(path, sep="\t", compression="gzip")


def load_biotype_map(id_map_path: str) -> pd.DataFrame:
    """Load Ensembl biotype annotations from id_map.tsv.

    Returns DataFrame with columns: ensembl_gene_id, mgi_symbol, biotype
    """
    df = pd.read_csv(id_map_path, sep="\t", comment="#")
    cols = ["ensembl_gene_id", "mgi_symbol", "biotype"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"id_map.tsv missing required column: {c}")
    return df[cols].copy()


# Regex patterns for noise genes to exclude regardless of biotype
_NOISE_SYMBOL_PATTERNS = [
    re.compile(r"^Gm\d+$"),          # predicted genes (Gm12345)
    re.compile(r"Rik$|Rik\d+$"),     # RIKEN cDNA clones
    re.compile(r"^\d+-[A-Z]"),       # numeric-prefix BAC clones
]


def _is_noise_symbol(symbol: str) -> bool:
    """Check if a gene symbol matches known noise patterns."""
    if not isinstance(symbol, str) or not symbol:
        return True  # unmapped genes are noise
    return any(p.search(symbol) for p in _NOISE_SYMBOL_PATTERNS)


def pick_genes(
    rtech: pd.DataFrame,
    max_genes: int,
    force_include: list[str] | None = None,
    biotype_map: pd.DataFrame | None = None,
    allowed_biotypes: list[str] | None = None,
    exclude_noise_symbols: bool = True,
) -> list[str]:
    """Select top genes by variance, with biotype filtering and force-include.

    Pipeline:
        1. Drop ERCC spike-ins
        2. (Optional) Filter to allowed biotypes via id_map.tsv
        3. (Optional) Exclude Gm-prefix and other noise symbols
        4. Force-include genes get priority slots
        5. HVGs fill remaining slots so total ≤ max_genes

    Args:
        rtech: Expression matrix (genes x samples), index = Ensembl IDs
        max_genes: Maximum total genes in panel
        force_include: Ensembl IDs to force-include (bypass biotype filter)
        biotype_map: DataFrame with ensembl_gene_id, mgi_symbol, biotype
        allowed_biotypes: List of allowed biotypes (e.g. ["protein_coding"])
        exclude_noise_symbols: If True, drop Gm\\d+, Rik, etc. from HVG pool
    """
    keep = ~rtech.index.str.upper().str.startswith("ERCC")
    r = rtech.loc[keep]
    n_before = len(r)

    # ── Biotype filter ───────────────────────────────────────────────────
    if biotype_map is not None and allowed_biotypes:
        bt = biotype_map.set_index("ensembl_gene_id")
        # Genes in allowed biotypes
        bt_pass = set(bt[bt["biotype"].isin(allowed_biotypes)].index)
        # Force-include genes bypass the biotype filter
        force_set = set(force_include) if force_include else set()
        # Keep genes that pass biotype OR are force-included
        bio_mask = r.index.isin(bt_pass | force_set)
        n_removed_biotype = (~bio_mask).sum()
        r = r.loc[bio_mask]
        print(f"  Biotype filter ({', '.join(allowed_biotypes)}):")
        print(f"    Removed {n_removed_biotype} non-{'/'.join(allowed_biotypes)} genes")
        print(f"    Remaining: {len(r)} / {n_before}")
    else:
        bt = None

    # ── Noise-symbol filter (Gm\d+, Rik, etc.) ──────────────────────────
    if exclude_noise_symbols and biotype_map is not None:
        if bt is None:
            bt = biotype_map.set_index("ensembl_gene_id")
        # Build set of noise Ensembl IDs (only from HVG candidates, not forced)
        force_set = set(force_include) if force_include else set()
        noise_ids = set()
        for eid in r.index:
            if eid in force_set:
                continue  # never exclude forced genes
            sym = bt.loc[eid, "mgi_symbol"] if eid in bt.index else ""
            if _is_noise_symbol(sym):
                noise_ids.add(eid)
        if noise_ids:
            n_noise = len(noise_ids)
            r = r.loc[~r.index.isin(noise_ids)]
            print(f"  Noise-symbol filter (Gm\\d+, Rik, unmapped):")
            print(f"    Removed {n_noise} noise-symbol genes")
            print(f"    Remaining: {len(r)}")

    # ── Force-include ────────────────────────────────────────────────────
    forced: set[str] = set()
    if force_include:
        present = set(r.index)
        forced = {g for g in force_include if g in present}
        missing = len(force_include) - len(forced)
        if missing:
            print(f"  Note: {missing}/{len(force_include)} forced genes not in Rtech")
        if len(forced) > max_genes:
            raise ValueError(
                f"Force-include genes ({len(forced)}) exceed max_genes ({max_genes}). "
                f"Increase --max_genes or reduce the marker panel."
            )

    # Fill remaining slots with top-variance genes (excluding already-forced)
    remaining_budget = max(max_genes - len(forced), 0)
    v = r.drop(index=list(forced), errors="ignore").var(axis=1)
    hvg = set(v.sort_values(ascending=False).head(remaining_budget).index.tolist())

    genes = list(forced | hvg)
    print(f"  Gene panel: {len(forced)} forced + {len(hvg)} HVG = {len(genes)} total")
    return genes


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
    ap.add_argument("--max_genes", type=int, default=2500,
                    help="Maximum genes for skeleton")
    ap.add_argument("--cell_cols", default="Age,Arm,EnvGroup",
                    help="Comma-separated columns defining experimental cells")
    ap.add_argument("--topk", type=int, default=80,
                    help="Top-k neighbors per gene (~G*k edges)")
    ap.add_argument("--force_include", default="",
                    help="Path to gene list file (one gene per line) to force-include")
    # ── Biotype filter arguments ─────────────────────────────────────────
    ap.add_argument("--id_map", default="data/processed/resources/id_map.tsv",
                    help="Path to id_map.tsv with biotype annotations "
                         "(from build_id_map.py)")
    ap.add_argument("--biotype_filter", default="protein_coding",
                    help="Comma-separated allowed biotypes. "
                         "Set to 'none' to disable. Default: protein_coding")
    ap.add_argument("--no_noise_symbol_filter", action="store_true",
                    help="Disable filtering of Gm-prefix / Rik / unmapped symbols")
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
    force_include = None
    if args.force_include and Path(args.force_include).exists():
        force_include = [
            line.strip() for line in
            Path(args.force_include).read_text().strip().split("\n")
            if line.strip()
        ]
        print(f"\nLoaded {len(force_include)} force-include genes from {args.force_include}")

    # Load biotype annotations if filtering is enabled
    biotype_map = None
    allowed_biotypes = None
    if args.biotype_filter.lower() != "none":
        id_map_path = Path(args.id_map)
        if not id_map_path.is_absolute():
            id_map_path = REPO_ROOT / id_map_path
        if id_map_path.exists():
            biotype_map = load_biotype_map(str(id_map_path))
            allowed_biotypes = [b.strip() for b in args.biotype_filter.split(",")]
            print(f"\nBiotype filter enabled: {allowed_biotypes}")
            print(f"  Loaded {len(biotype_map)} annotations from {id_map_path}")
        else:
            print(f"\n  WARNING: id_map not found at {id_map_path}, "
                  f"skipping biotype filter. Run build_id_map.py first.")

    genes = pick_genes(
        rtech,
        args.max_genes,
        force_include=force_include,
        biotype_map=biotype_map,
        allowed_biotypes=allowed_biotypes,
        exclude_noise_symbols=not args.no_noise_symbol_filter,
    )
    (outdir / "phase2_genes.txt").write_text("\n".join(genes) + "\n")
    print(f"\nSelected {len(genes)} genes for skeleton")
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
    print("Skeleton E built successfully")
    print(f"{'=' * 60}")
    print(f"  Genes: {len(genes)}")
    print(f"  Edges: {len(edge_df)} (target: ~{len(genes)*args.topk//2} to ~{len(genes)*args.topk})")
    print(f"\nOutputs in {outdir}:")
    print(f"  - phase2_genes.txt")
    print(f"  - skeleton_edges.tsv (with pcorr weights)")
    print(f"  - edge_i.npy, edge_j.npy (indices)")


if __name__ == "__main__":
    main()
