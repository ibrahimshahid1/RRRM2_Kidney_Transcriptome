#!/usr/bin/env python3
"""Derive the reconstructed gene map from the Ensembl-116 map (config/gene_map_reconstruction.yaml)."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[2]


def build(source: pd.DataFrame, excluded_ids) -> pd.DataFrame:
    """Keep symbolized rows not in ``excluded_ids``, as (ensembl_gene_id, mgi_symbol)."""
    symbols = source["mgi_symbol"].fillna("").astype(str).str.strip()
    keep = (symbols != "") & ~source["ensembl_gene_id"].isin(set(excluded_ids))
    return source.loc[keep, ["ensembl_gene_id", "mgi_symbol"]].reset_index(drop=True)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=REPO / "config/gene_map_reconstruction.yaml")
    parser.add_argument(
        "--no-verify", action="store_true", help="skip the source/output sha256 checks"
    )
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    source_path = REPO / cfg["source_map"]["path"]
    output_path = REPO / cfg["output_map"]["path"]
    if not args.no_verify and sha256(source_path) != cfg["source_map"]["sha256"]:
        sys.exit(f"source map {source_path} does not match the recorded sha256")

    source = pd.read_csv(source_path, sep="\t", comment="#", dtype=str)
    missing = set(cfg["excluded_ids"]) - set(source["ensembl_gene_id"])
    if missing:
        sys.exit(f"excluded IDs absent from the source map: {sorted(missing)}")

    output = build(source, cfg["excluded_ids"])
    output.to_csv(output_path, sep="\t", index=False)
    digest = sha256(output_path)
    print(f"wrote {output_path} ({len(output)} genes, sha256 {digest})")
    if not args.no_verify and digest != cfg["output_map"]["sha256"]:
        sys.exit("reconstructed map does not match the recorded sha256")


if __name__ == "__main__":
    main()
