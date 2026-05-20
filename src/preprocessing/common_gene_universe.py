"""Common gene-universe utilities for the Cross-OSDR framework.

The implementation only intersects already processed within-study matrices.
It does not pool raw expression across studies.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

from src.common import REPO_ROOT


def load_gene_index(path: str | Path) -> pd.Index:
    """Load a gene-by-sample matrix and return normalized gene IDs."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Matrix not found: {path}")
    suffixes = path.suffixes
    compression = "gzip" if ".gz" in suffixes else None
    sep = "\t" if ".tsv" in suffixes else ","
    first_col = pd.read_csv(path, sep=sep, compression=compression, usecols=[0]).iloc[:, 0]
    return pd.Index(first_col.astype(str).str.replace(r"\.\d+$", "", regex=True))


def intersect_gene_universe(matrix_paths: Iterable[str | Path]) -> list[str]:
    """Return sorted Ensembl IDs present in every matrix."""
    paths = [Path(p) for p in matrix_paths]
    if not paths:
        raise ValueError("At least one matrix path is required")
    universe: set[str] | None = None
    for path in paths:
        genes = set(load_gene_index(path))
        universe = genes if universe is None else universe & genes
    return sorted(universe or set())


def write_gene_universe(universe: list[str], outpath: str | Path, manifest: dict | None = None) -> None:
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text("\n".join(universe) + ("\n" if universe else ""))
    if manifest is not None:
        with open(outpath.with_suffix(outpath.suffix + ".manifest.json"), "w") as fh:
            json.dump(manifest, fh, indent=2)


def main() -> int:
    ap = argparse.ArgumentParser(description="Build common gene universe from processed matrices")
    ap.add_argument("--matrices", nargs="+", required=True)
    ap.add_argument("--out", default="data/processed/multi_study_gene_universe.tsv")
    args = ap.parse_args()
    paths = [REPO_ROOT / p for p in args.matrices]
    universe = intersect_gene_universe(paths)
    write_gene_universe(
        universe,
        REPO_ROOT / args.out,
        manifest={"matrices": args.matrices, "n_genes": len(universe)},
    )
    print(f"Wrote {len(universe)} common genes to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
