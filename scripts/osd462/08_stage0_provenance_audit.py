#!/usr/bin/env python3
"""Generate the OSD-462 Stage 0 design and phosphosite provenance audit."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.multiomics.osd462_stage0 import (  # noqa: E402
    audit_ncc_spak_phosphosites,
    build_sample_design,
    build_stage0_qc,
    isolated_canonical_assay_features,
    sha256_file,
    source_files,
)


DEFAULT_OUT = REPO_ROOT / "data" / "results" / "run_20260728_osd462_stage0"


def _git_commit() -> str:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unavailable"


def _git_dirty_state() -> dict[str, object]:
    """Return deterministic porcelain lines plus their digest."""
    result = subprocess.run(
        ["git", "status", "--porcelain=v1", "--untracked-files=all"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return {
            "git_dirty": None,
            "git_status_porcelain": [],
            "git_status_porcelain_sha256": "unavailable",
        }
    lines = sorted(line for line in result.stdout.splitlines() if line)
    serialized = "\n".join(lines) + ("\n" if lines else "")
    return {
        "git_dirty": bool(lines),
        "git_status_porcelain": lines,
        "git_status_porcelain_sha256": hashlib.sha256(
            serialized.encode("utf-8")
        ).hexdigest(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Trace OSD-462 samples, reporter channels, raw acquisitions, and "
            "NCC/SPAK phosphosite residue identities from source files."
        )
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT,
        help=f"output directory (default: {DEFAULT_OUT.relative_to(REPO_ROOT)})",
    )
    args = parser.parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    design, raw_inventory = build_sample_design()
    site_audit = audit_ncc_spak_phosphosites()
    isolated_features = isolated_canonical_assay_features(site_audit)
    qc = build_stage0_qc(design, raw_inventory, site_audit)

    outputs = {
        "sample_design": out_dir / "osd462_sample_plex_channel_design.tsv",
        "raw_inventory": out_dir / "osd462_raw_file_inventory.tsv",
        "phosphosite_audit": out_dir / "osd462_ncc_spak_phosphosite_provenance.tsv",
        "isolated_canonical_features": (
            out_dir / "osd462_isolated_canonical_assay_features.tsv"
        ),
        "qc": out_dir / "osd462_stage0_qc.tsv",
    }
    design.to_csv(outputs["sample_design"], sep="\t", index=False)
    raw_inventory.to_csv(outputs["raw_inventory"], sep="\t", index=False)
    site_audit.to_csv(outputs["phosphosite_audit"], sep="\t", index=False)
    isolated_features.to_csv(
        outputs["isolated_canonical_features"], sep="\t", index=False
    )
    qc.to_csv(outputs["qc"], sep="\t", index=False)

    manifest = {
        "analysis": "OSD-462 Stage 0 sample and phosphosite provenance audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit_at_generation": _git_commit(),
        **_git_dirty_state(),
        "command": (
            "venv/bin/python scripts/osd462/08_stage0_provenance_audit.py "
            f"--out-dir {out_dir.relative_to(REPO_ROOT)}"
        ),
        "source_sha256": {
            str(path.relative_to(REPO_ROOT)): sha256_file(path)
            for path in source_files()
        },
        "output_sha256": {
            str(path.relative_to(REPO_ROOT)): sha256_file(path)
            for path in outputs.values()
        },
        "row_counts": {
            "sample_design": int(len(design)),
            "raw_inventory": int(len(raw_inventory)),
            "phosphosite_audit": int(len(site_audit)),
            "isolated_canonical_features": int(len(isolated_features)),
            "qc": int(len(qc)),
        },
        "qc_status_counts": {
            str(key): int(value)
            for key, value in qc["status"].value_counts().sort_index().items()
        },
        "warnings_are_expected_design_or_metadata_findings": True,
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print(f"Wrote Stage 0 audit to {out_dir}")
    print(qc[["metric", "status", "value"]].to_string(index=False))
    failures = int(qc["status"].eq("FAIL").sum())
    if failures:
        print(f"Stage 0 failed {failures} QC check(s).", file=sys.stderr)
        return 1
    print(
        f"Stage 0 completed with {int(qc['status'].eq('WARN').sum())} "
        "documented warning(s) and no failures."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
