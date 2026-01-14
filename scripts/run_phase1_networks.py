# scripts/run_phase1_networks.py
"""
Phase 1 Network Analysis Pipeline

Runs LIONESS on Phase 1 residualized expression data and computes
group-level networks and rewiring metrics.

Usage:
    python scripts/run_phase1_networks.py
    python scripts/run_phase1_networks.py --max_genes 1000 --compare "EnvGroup:FLT-vs-GC"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd

from src.networks.lioness import lioness_correlation_edges, save_lioness_bundle
from src.statistics.rewiring_metrics import node_rewiring_from_delta_edges


def load_rtech_tsv(path: str) -> pd.DataFrame:
    """Load Rtech matrix (genes x samples) from gzipped TSV."""
    df = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)
    print(f"Loaded Rtech: {df.shape[0]} genes × {df.shape[1]} samples")
    return df


def pick_top_variable_genes(rtech: pd.DataFrame, max_genes: int) -> list:
    """Select top genes by variance across samples."""
    v = rtech.var(axis=1)
    genes = v.sort_values(ascending=False).head(max_genes).index.tolist()
    print(f"Selected {len(genes)} top-variance genes")
    return genes


def main():
    ap = argparse.ArgumentParser(description="Run LIONESS on Phase 1 data")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz",
                    help="Path to Rtech.tsv.gz (genes x samples)")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz",
                    help="Path to meta_phase1.tsv.gz")
    ap.add_argument("--outdir", default="data/processed/networks/phase1",
                    help="Output directory for networks")
    ap.add_argument("--max_genes", type=int, default=1200,
                    help="Maximum genes for network (800-1500 recommended)")
    ap.add_argument("--group_cols", default="Age,Arm,EnvGroup",
                    help="Comma-separated columns for grouping")
    ap.add_argument("--compare", default="EnvGroup:FLT-vs-GC",
                    help="Comparison spec: COL:A-vs-B")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load data
    print("=" * 60)
    print("Phase 1 LIONESS Network Pipeline")
    print("=" * 60)
    
    print(f"\nLoading Rtech: {args.rtech}")
    rtech = load_rtech_tsv(args.rtech)

    ercc_mask = rtech.index.astype(str).str.upper().str.startswith("ERCC")
    n_ercc = int(ercc_mask.sum())
    if n_ercc:
        print(f"Removing {n_ercc} ERCC genes from Rtech entirely")
        rtech = rtech.loc[~ercc_mask]

    print(f"Loading metadata: {args.meta}")
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    
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
    
    # Align meta to Rtech columns
    common_samples = [s for s in rtech.columns if s in meta.index]
    if len(common_samples) < len(rtech.columns):
        print(f"  Warning: {len(rtech.columns) - len(common_samples)} samples not in metadata")
    rtech = rtech[common_samples]
    meta = meta.loc[common_samples]
    
    print(f"  Aligned: {len(common_samples)} samples")
    assert list(meta.index) == list(rtech.columns), "Meta not aligned to Rtech"

    # Gene universe
    genes = pick_top_variable_genes(rtech, args.max_genes)
    (outdir / "phase1_genes.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    print(f"  → Saved gene list to {outdir / 'phase1_genes.txt'}")

    # Prepare expression matrix (samples x genes)
    X = rtech.loc[genes].T.values
    gene_names = genes
    sample_names = rtech.columns.tolist()
    
    print(f"\nExpression matrix: {X.shape[0]} samples × {X.shape[1]} genes")

    # LIONESS
    print("\n" + "-" * 40)
    print("Running LIONESS (correlation-based)...")
    print("-" * 40)
    edges, index = lioness_correlation_edges(X, gene_names=gene_names, dtype=np.float32)
    print(f"LIONESS complete: {edges.shape}")

    # Diagnostics
    G = len(gene_names)
    expected_edges = G*(G-1)//2
    print("Expected edges:", expected_edges, "Actual:", edges.shape[1])
    print("edges dtype:", edges.dtype)
    print("min/max:", np.nanmin(edges), np.nanmax(edges))
    print("row0==row1:", np.allclose(edges[0], edges[1]))

    # Save LIONESS bundle
    bundle_dir = outdir / "lioness"
    save_lioness_bundle(bundle_dir, edges, index)

    # Grouping
    print("\n" + "-" * 40)
    print("Computing group-mean networks...")
    print("-" * 40)
    
    group_cols = [c.strip() for c in args.group_cols.split(",") if c.strip()]
    missing_cols = [c for c in group_cols if c not in meta.columns]
    if missing_cols:
        print(f"  Warning: columns not found: {missing_cols}")
        group_cols = [c for c in group_cols if c in meta.columns]
    
    if group_cols:
        group_keys = meta[group_cols].astype(str).agg("|".join, axis=1)
        meta["_group_key"] = group_keys.values

        # Compute group means
        group_means = {}
        for g, sub_df in meta.groupby("_group_key"):
            idxs = [sample_names.index(s) for s in sub_df.index if s in sample_names]
            if len(idxs) > 0:
                group_means[g] = edges[idxs, :].mean(axis=0)
                print(f"  {g}: {len(idxs)} samples")

        # Save group means
        group_mean_path = outdir / "group_mean_edges.npz"
        np.savez_compressed(group_mean_path, **{g: v.astype(np.float32) for g, v in group_means.items()})
        
        groups_df = pd.DataFrame({"group": list(group_means.keys()), "n_samples": [
            len([s for s in meta[meta["_group_key"] == g].index]) for g in group_means.keys()
        ]})
        groups_df.to_csv(outdir / "groups_present.csv", index=False)
        print(f"  → Saved to {group_mean_path}")
    else:
        print("  No valid group columns - skipping group means")

    # Comparison: FLT vs GC
    print("\n" + "-" * 40)
    print("Computing rewiring scores...")
    print("-" * 40)
    
    comp = args.compare
    try:
        col, rest = comp.split(":")
        a_name, b_name = rest.split("-vs-")
    except ValueError:
        print(f"  Error: Invalid compare format '{comp}'. Expected 'COL:A-vs-B'")
        return

    if col not in meta.columns:
        print(f"  Error: Column '{col}' not in metadata")
        return

    # Build stratum keys (excluding comparison column)
    strata_cols = [c for c in group_cols if c != col]
    if strata_cols:
        strata = meta[strata_cols].astype(str).agg("|".join, axis=1)
    else:
        strata = pd.Series(["ALL"] * len(meta), index=meta.index)
    meta["_stratum"] = strata

    rewiring_tables = []
    for s, sub in meta.groupby("_stratum"):
        A = sub[sub[col].astype(str) == a_name]
        B = sub[sub[col].astype(str) == b_name]
        
        if len(A) < 2 or len(B) < 2:
            print(f"  Stratum {s}: skipped (A={len(A)}, B={len(B)} samples)")
            continue

        # Compute mean edges for each group
        A_idx = [sample_names.index(x) for x in A.index if x in sample_names]
        B_idx = [sample_names.index(x) for x in B.index if x in sample_names]
        
        meanA = edges[A_idx, :].mean(axis=0)
        meanB = edges[B_idx, :].mean(axis=0)
        delta = meanA - meanB

        # Node rewiring scores
        df = node_rewiring_from_delta_edges(delta, index)
        df.insert(0, "stratum", s)
        df.insert(1, "comparison", f"{a_name}-vs-{b_name}")
        rewiring_tables.append(df)
        print(f"  Stratum {s}: {a_name}={len(A_idx)}, {b_name}={len(B_idx)} → computed")

    if rewiring_tables:
        out_df = pd.concat(rewiring_tables, ignore_index=True)
        out_path = outdir / f"rewiring_{a_name}_vs_{b_name}.tsv"
        out_df.to_csv(out_path, sep="\t", index=False)
        print(f"\n  → Saved rewiring table: {out_path}")
        print(f"     {len(out_df)} gene × stratum rows")
        print("\n  Preview (first 10 rows):")
        print(out_df.head(10).to_string(index=False))
    else:
        print("\n  No rewiring tables produced (insufficient samples per group)")

    # -------------------------------------------------------------------------
    # Second Comparison: Age (YNG vs OLD) within each Arm × EnvGroup
    # -------------------------------------------------------------------------
    print("\n" + "-" * 40)
    print("Computing Age rewiring scores (YNG vs OLD)...")
    print("-" * 40)
    
    age_comparisons = []
    age_col = "Age"
    age_a, age_b = "YNG", "OLD"
    
    if age_col not in meta.columns:
        print(f"  Warning: '{age_col}' column not found in metadata, skipping age comparison")
    else:
        # Stratify by Arm × EnvGroup
        age_strata_cols = ["Arm", "EnvGroup"]
        age_strata_cols = [c for c in age_strata_cols if c in meta.columns]
        
        if age_strata_cols:
            age_strata = meta[age_strata_cols].astype(str).agg("|".join, axis=1)
        else:
            age_strata = pd.Series(["ALL"] * len(meta), index=meta.index)
        meta["_age_stratum"] = age_strata
        
        for s, sub in meta.groupby("_age_stratum"):
            A = sub[sub[age_col].astype(str) == age_a]
            B = sub[sub[age_col].astype(str) == age_b]
            
            if len(A) < 2 or len(B) < 2:
                print(f"  Stratum {s}: skipped ({age_a}={len(A)}, {age_b}={len(B)} samples)")
                continue
            
            A_idx = [sample_names.index(x) for x in A.index if x in sample_names]
            B_idx = [sample_names.index(x) for x in B.index if x in sample_names]
            
            meanA = edges[A_idx, :].mean(axis=0)
            meanB = edges[B_idx, :].mean(axis=0)
            delta = meanA - meanB
            
            df = node_rewiring_from_delta_edges(delta, index)
            df.insert(0, "stratum", s)
            df.insert(1, "comparison", f"{age_a}-vs-{age_b}")
            age_comparisons.append(df)
            print(f"  Stratum {s}: {age_a}={len(A_idx)}, {age_b}={len(B_idx)} → computed")
        
        if age_comparisons:
            age_df = pd.concat(age_comparisons, ignore_index=True)
            age_out_path = outdir / f"rewiring_{age_a}_vs_{age_b}.tsv"
            age_df.to_csv(age_out_path, sep="\t", index=False)
            print(f"\n  → Saved age rewiring table: {age_out_path}")
            print(f"     {len(age_df)} gene × stratum rows")
            print("\n  Preview (first 10 rows):")
            print(age_df.head(10).to_string(index=False))
        else:
            print("\n  No age rewiring tables produced (insufficient samples per group)")

    print("\n" + "=" * 60)
    print("Pipeline complete!")
    print("=" * 60)
    print(f"\nOutputs in: {outdir}")
    print(f"  - phase1_genes.txt: gene universe")
    print(f"  - lioness/: sample-specific networks")
    print(f"  - group_mean_edges.npz: group-averaged networks")
    print(f"  - rewiring_*.tsv: node rewiring scores (FLT-vs-GC and YNG-vs-OLD)")


if __name__ == "__main__":
    main()
