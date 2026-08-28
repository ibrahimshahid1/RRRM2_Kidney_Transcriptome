#!/usr/bin/env python3
"""Write a non-recursive provenance manifest for a completed v13 run."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_NAME = "provenance_manifest.json"
CODE_PATHS = (
    "src/v13/continuous_phospho_inference.py",
    "scripts/v13/run_continuous_phospho_inference.py",
    "src/subtype_reference/reference_builder.py",
    "scripts/subtype_reference/04_build_frozen_signatures.py",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(args: Sequence[str]) -> str | None:
    try:
        return subprocess.check_output(
            ["git", *args],
            cwd=REPO_ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).rstrip("\n")
    except (OSError, subprocess.CalledProcessError):
        return None


def build_provenance_manifest(run_dir: Path) -> dict:
    run_dir = run_dir.resolve()
    inference_manifest_path = run_dir / "manifest.json"
    if not inference_manifest_path.exists():
        raise FileNotFoundError(
            f"Completed inference manifest is missing: {inference_manifest_path}"
        )
    inference = json.loads(inference_manifest_path.read_text())
    if not (run_dir / "claim_tier.tsv").exists():
        raise FileNotFoundError(
            f"Run is incomplete; claim_tier.tsv is missing from {run_dir}"
        )

    output_hashes = {}
    for path in sorted(item for item in run_dir.rglob("*") if item.is_file()):
        if path.name == MANIFEST_NAME:
            continue
        output_hashes[str(path.relative_to(run_dir))] = sha256_file(path)

    code_hashes = {}
    for value in CODE_PATHS:
        path = REPO_ROOT / value
        if path.exists():
            code_hashes[value] = sha256_file(path)

    recorded_inputs = inference.get("inputs", {})
    input_hashes = {
        key: value
        for key, value in recorded_inputs.items()
        if key.endswith("_sha256")
    }
    broad_flags = (
        Path(recorded_inputs.get("frozen_gene_sets", "")).with_name(
            "broad_expression_flags.tsv"
        )
        if recorded_inputs.get("frozen_gene_sets")
        else None
    )
    if broad_flags is not None and not broad_flags.is_absolute():
        broad_flags = (REPO_ROOT / broad_flags).resolve()
    if broad_flags is not None and broad_flags.exists():
        input_hashes["broad_expression_flags_sha256"] = sha256_file(broad_flags)

    status = _git(["status", "--porcelain=v1", "--untracked-files=all"])
    return {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "analysis_id": inference.get("analysis_id"),
        "permutation": inference.get("permutation"),
        "source_inference_manifest_sha256": sha256_file(inference_manifest_path),
        "input_hashes": input_hashes,
        "code_hashes": code_hashes,
        "output_hashes_excluding_this_manifest": output_hashes,
        "git_commit": _git(["rev-parse", "HEAD"]),
        "git_status_porcelain": status.splitlines() if status else [],
        "git_status_sha256": hashlib.sha256(
            ((status or "") + ("\n" if status else "")).encode()
        ).hexdigest(),
        "condition_reporter_position_confounded": inference.get(
            "condition_reporter_position_confounded"
        ),
        "analysis_recomputed": False,
        "note": (
            "This manifest was written after inference and reporting. It hashes "
            "all run artifacts except itself to avoid a recursive self-hash."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-dir", required=True)
    args = parser.parse_args(argv)
    run_dir = Path(args.run_dir)
    manifest = build_provenance_manifest(run_dir)
    output = run_dir / MANIFEST_NAME
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
