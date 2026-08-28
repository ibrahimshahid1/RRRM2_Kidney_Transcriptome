#!/usr/bin/env python3
"""Inventory local distal-nephron reference inputs without downloading data."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def _checksum(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _row(
    dataset: str,
    role: str,
    path: Path | None,
    required: bool,
    expected_sha256: str | None = None,
    expected_md5: str | None = None,
    note: str = "",
    compute_hash: bool = False,
) -> dict[str, Any]:
    present = bool(path and path.exists())
    actual_sha256 = _checksum(path, "sha256") if present and compute_hash else None
    actual_md5 = _checksum(path, "md5") if present and compute_hash and expected_md5 else None
    hash_matches = (
        actual_sha256 == expected_sha256
        if actual_sha256 is not None and expected_sha256 is not None
        else actual_md5 == expected_md5
        if actual_md5 is not None and expected_md5 is not None
        else None
    )
    return {
        "dataset": dataset,
        "role": role,
        "required": required,
        "path": str(path) if path else "",
        "present": present,
        "bytes": path.stat().st_size if present else None,
        "expected_sha256": expected_sha256,
        "actual_sha256": actual_sha256,
        "expected_md5": expected_md5,
        "actual_md5": actual_md5,
        "hash_matches": hash_matches,
        "note": note,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument(
        "--hash",
        action="store_true",
        help="Compute SHA-256 for large reference files (off by default).",
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    cfg = yaml.safe_load(config_path.read_text()) or {}
    inputs = cfg["reference_builder"]["inputs"]
    rows: list[dict[str, Any]] = []
    for dataset, entries in inputs.items():
        for entry in entries:
            path_text = entry.get("path")
            path = Path(path_text) if path_text else None
            rows.append(
                _row(
                    dataset=dataset,
                    role=entry["role"],
                    path=path,
                    required=bool(entry.get("required", False)),
                    expected_sha256=entry.get("sha256"),
                    expected_md5=entry.get("md5"),
                    note=entry.get("note", ""),
                    compute_hash=args.hash,
                )
            )

    table = pd.DataFrame(rows)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    table.to_csv(out / "reference_input_inventory.tsv", sep="\t", index=False)
    missing = table[table["required"] & ~table["present"]]
    report = {
        "freeze_id": cfg.get("freeze_id"),
        "config": str(config_path),
        "downloads_attempted": False,
        "required_inputs_complete": bool(missing.empty),
        "missing_required_inputs": missing[
            ["dataset", "role", "path", "note"]
        ].to_dict(orient="records"),
        "present_optional_inputs": table[
            ~table["required"] & table["present"]
        ][["dataset", "role", "path"]].to_dict(orient="records"),
    }
    (out / "reference_input_inventory.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
