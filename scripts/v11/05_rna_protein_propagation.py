#!/usr/bin/env python3
"""Run the v11 RNA-to-protein propagation extension."""

from __future__ import annotations

import argparse
from pathlib import Path

import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.v11.rna_protein_propagation import (  # noqa: E402
    DEFAULT_ANCHOR_DIR,
    DEFAULT_GENE_SETS,
    DEFAULT_RUN_ROOT,
    PropagationConfig,
    run_propagation,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--anchor-dir", type=Path, default=DEFAULT_ANCHOR_DIR)
    parser.add_argument("--gene-sets", type=Path, default=DEFAULT_GENE_SETS)
    parser.add_argument("--n-null", type=int, default=10000)
    parser.add_argument("--n-bootstrap", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260606)
    parser.add_argument("--peptide-filter", type=int, default=2)
    args = parser.parse_args()

    cfg = PropagationConfig(
        run_root=args.run_root,
        anchor_dir=args.anchor_dir,
        gene_sets_path=args.gene_sets,
        n_null=args.n_null,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        peptide_filter=args.peptide_filter,
    )
    outputs = run_propagation(cfg)
    print("[propagation] wrote:")
    for name, path in outputs.items():
        print(f"  {name}: {path}")


if __name__ == "__main__":
    main()
