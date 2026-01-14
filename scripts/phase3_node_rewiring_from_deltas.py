# scripts/phase3_node_rewiring_from_deltas.py
"""
Phase 3.1: Convert *_delta_z.npy to per-gene rewiring scores

Simple, dependency-free node rewiring metrics:
  - rewiring_abs: sum of abs(delta) over incident edges
  - rewiring_signed: sum of delta over incident edges  
  - degree: number of incident edges

Usage:
    python scripts/phase3_node_rewiring_from_deltas.py
"""
from __future__ import annotations

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def node_rewiring_from_delta_edges(
    delta: np.ndarray, 
    edge_i: np.ndarray, 
    edge_j: np.ndarray, 
    genes: list[str]
) -> pd.DataFrame:
    """
    Simple, dependency-free node rewiring metrics:
      - rewiring_abs: sum of abs(delta) over incident edges
      - rewiring_signed: sum of delta over incident edges
      - degree: number of incident edges counted (each edge contributes to both endpoints)
    """
    if delta.ndim != 1:
        raise ValueError(f"delta must be 1D, got shape {delta.shape}")
    if edge_i.shape[0] != delta.shape[0] or edge_j.shape[0] != delta.shape[0]:
        raise ValueError(
            f"Edge count mismatch: delta={delta.shape[0]} edge_i={edge_i.shape[0]} edge_j={edge_j.shape[0]}"
        )

    G = len(genes)
    abs_sum = np.zeros(G, dtype=np.float64)
    signed_sum = np.zeros(G, dtype=np.float64)
    deg = np.zeros(G, dtype=np.int64)

    dabs = np.abs(delta).astype(np.float64)

    np.add.at(abs_sum, edge_i, dabs)
    np.add.at(abs_sum, edge_j, dabs)

    np.add.at(signed_sum, edge_i, delta.astype(np.float64))
    np.add.at(signed_sum, edge_j, delta.astype(np.float64))

    np.add.at(deg, edge_i, 1)
    np.add.at(deg, edge_j, 1)

    df = pd.DataFrame(
        {
            "gene": genes,
            "rewiring_abs": abs_sum,
            "rewiring_signed": signed_sum,
            "degree": deg,
        }
    ).sort_values("rewiring_abs", ascending=False, kind="mergesort").reset_index(drop=True)

    return df


def main():
    os.chdir(REPO_ROOT)

    ap = argparse.ArgumentParser(description="Phase 3.1: Convert *_delta_z.npy to per-gene rewiring scores")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2")
    ap.add_argument("--reg_dir", default="data/processed/networks/phase2/regression")
    ap.add_argument("--outdir", default="data/processed/networks/phase3/node_rewiring")
    ap.add_argument("--pattern", default="*_delta_z.npy")
    args = ap.parse_args()

    phase2_dir = Path(args.phase2_dir)
    reg_dir = Path(args.reg_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 3.1: Node Rewiring from Delta-z")
    print("=" * 60)

    genes = [g for g in (phase2_dir / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    edge_i = np.load(phase2_dir / "edge_i.npy")
    edge_j = np.load(phase2_dir / "edge_j.npy")

    print(f"\nLoaded {len(genes)} genes, {len(edge_i)} edges")

    files = sorted(reg_dir.glob(args.pattern))
    if not files:
        raise FileNotFoundError(f"No files matched {args.pattern} in {reg_dir}")

    print(f"Found {len(files)} delta files to process\n")

    for f in files:
        delta = np.load(f).astype(np.float64)
        df = node_rewiring_from_delta_edges(delta, edge_i, edge_j, genes)

        out = outdir / (f.stem.replace("_delta_z", "") + "_node_rewiring.tsv")
        df.to_csv(out, sep="\t", index=False)
        print(f"[OK] {f.name} -> {out.name} (rows={len(df)})")

    print(f"\nDone. Outputs in: {outdir}")


if __name__ == "__main__":
    main()
