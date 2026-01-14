# scripts/pick_anchors.py
"""
Select anchor genes for Procrustes alignment.

Anchors are genes with minimal network rewiring between conditions,
making them suitable reference points for embedding alignment.

Usage:
    python scripts/pick_anchors.py
    python scripts/pick_anchors.py --groupA "YNG|ISS-T|FLT" --groupB "OLD|ISS-T|FLT" --k 150
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Add project root to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np
import pandas as pd

from src.statistics.rewiring_metrics import node_rewiring_from_delta_edges

OUTDIR = Path("data/processed/networks/phase1")
LIONESS_DIR = OUTDIR / "lioness"


def load_index():
    """Load edge index arrays (triu_i, triu_j)."""
    triu_i = np.load(LIONESS_DIR / "triu_i.npy")
    triu_j = np.load(LIONESS_DIR / "triu_j.npy")
    return (triu_i, triu_j)


def load_genes():
    """Load gene list from LIONESS bundle."""
    return (LIONESS_DIR / "genes.txt").read_text().splitlines()


def pick_anchors(groupA: str, groupB: str, k: int = 150, outname: str | None = None):
    """
    Select k anchor genes with minimal rewiring between two groups.
    
    Parameters
    ----------
    groupA, groupB : str
        Group keys (e.g., "YNG|ISS-T|FLT")
    k : int
        Number of anchors to select
    outname : str, optional
        Output filename (default: auto-generated)
    """
    means = np.load(OUTDIR / "group_mean_edges.npz")
    index = load_index()
    genes = load_genes()

    if groupA not in means.files:
        raise KeyError(f"{groupA} not found in group_mean_edges.npz")
    if groupB not in means.files:
        raise KeyError(f"{groupB} not found in group_mean_edges.npz")

    delta = means[groupA] - means[groupB]
    
    # Build EdgeIndex-like object for the rewiring function
    class EdgeIndex:
        def __init__(self, i, j, names):
            self.triu_i = i
            self.triu_j = j
            self.gene_names = names
    
    edge_index = EdgeIndex(index[0], index[1], genes)
    df = node_rewiring_from_delta_edges(delta, edge_index)

    # Anchors = smallest absolute rewiring
    anchors = df.nsmallest(k, "rewiring_abs").copy()

    if outname is None:
        outname = f"anchors_{groupA}_vs_{groupB}_k{k}.tsv".replace("|", "_")

    outpath = OUTDIR / outname
    anchors.to_csv(outpath, sep="\t", index=False)
    print(f"Saved {k} anchors to: {outpath}")
    print("Preview (first 10):")
    print(anchors.head(10).to_string(index=False))
    
    return anchors


def main():
    ap = argparse.ArgumentParser(description="Pick anchor genes for Procrustes alignment")
    ap.add_argument("--groupA", default="YNG|ISS-T|FLT",
                    help="First group key")
    ap.add_argument("--groupB", default="OLD|ISS-T|FLT",
                    help="Second group key")
    ap.add_argument("--k", type=int, default=150,
                    help="Number of anchors to select")
    ap.add_argument("--outname", default=None,
                    help="Output filename (default: auto)")
    args = ap.parse_args()
    
    pick_anchors(args.groupA, args.groupB, k=args.k, outname=args.outname)


if __name__ == "__main__":
    main()
