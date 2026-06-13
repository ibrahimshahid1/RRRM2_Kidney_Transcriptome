#!/usr/bin/env python3
"""Convenience auditor for the Phase 3 stability gate (Guardrail A)."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import yaml

SCRIPT_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_REPO_ROOT))

from src.networks.stability_test import gates_from_config  # noqa: E402


def _latest_contrast_vectors_dir() -> Path | None:
    results_root = SCRIPT_REPO_ROOT / "data" / "results"
    if not results_root.exists():
        return None
    matches = sorted(results_root.glob("run_*/contrast_vectors"), reverse=True)
    for match in matches:
        if match.is_dir():
            return match
    return None


def resolve_outdir(arg: str | None) -> Path:
    if arg:
        path = Path(arg)
        return path if path.is_absolute() else SCRIPT_REPO_ROOT / path
    latest = _latest_contrast_vectors_dir()
    if latest is None:
        raise SystemExit(
            "Could not locate any data/results/run_*/contrast_vectors directory. "
            "Pass --outdir explicitly."
        )
    return latest


def load_config(path: Path) -> dict:
    if not path.is_absolute():
        path = SCRIPT_REPO_ROOT / path
    if not path.exists():
        raise SystemExit(f"Config file not found: {path}")
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def check_gates_match_config(decisions: list[dict], gates: dict[str, dict[str, float]]) -> list[str]:
    """Verify each recorded decision's required thresholds match the config."""
    discrepancies: list[str] = []
    for row in decisions:
        resolution = row.get("resolution")
        gate = gates.get(resolution)
        if gate is None:
            discrepancies.append(
                f"resolution={resolution} has a recorded decision but no configured gate"
            )
            continue
        if abs(float(row.get("median_required", -1)) - gate["median"]) > 1e-9:
            discrepancies.append(
                f"{resolution}: recorded median_required={row.get('median_required')} "
                f"does not match configured {gate['median']}"
            )
        if abs(float(row.get("lower_required", -1)) - gate["lower"]) > 1e-9:
            discrepancies.append(
                f"{resolution}: recorded lower_required={row.get('lower_required')} "
                f"does not match configured {gate['lower']}"
            )
    return discrepancies


def main() -> int:
    ap = argparse.ArgumentParser(description="Audit the Cross-OSDR stability gate artifacts.")
    ap.add_argument("--outdir", default="",
                    help="Path to data/results/<run>/contrast_vectors. Defaults to latest run.")
    ap.add_argument("--config", default="config/contrast_vector_framework.yaml",
                    help="Path to the contrast-vector framework config.")
    ap.add_argument("--report-filename", default="agc_stability_report.tsv")
    ap.add_argument("--decision-filename", default="agc_stability_decision.json")
    ap.add_argument("--strict", action="store_true",
                    help="Treat any failed gate (other than configured fallbacks) as an error.")
    args = ap.parse_args()

    outdir = resolve_outdir(args.outdir)
    print(f"[audit] Inspecting {outdir}")

    report_path = outdir / args.report_filename
    decision_path = outdir / args.decision_filename
    missing: list[str] = []
    if not report_path.exists():
        missing.append(str(report_path))
    if not decision_path.exists():
        missing.append(str(decision_path))
    if missing:
        print("[audit] Missing artifacts:")
        for m in missing:
            print(f"  - {m}")
        print("[audit] Run scripts/run_contrast_vector_framework.py --phases stability first.")
        return 2

    config = load_config(Path(args.config))
    gates = gates_from_config(config)

    with open(decision_path) as fh:
        decision = json.load(fh)
    report = pd.read_csv(report_path, sep="\t")

    decisions = decision.get("decisions", [])
    discrepancies = check_gates_match_config(decisions, gates)
    if discrepancies:
        print("[audit] Threshold discrepancies vs configured gates:")
        for d in discrepancies:
            print(f"  - {d}")

    print("[audit] Stability report:")
    cols = ["arm", "resolution", "cosine_median", "cosine_lower", "cosine_upper", "n_bootstrap"]
    cols = [c for c in cols if c in report.columns]
    print(report[cols].to_string(index=False))

    print("[audit] Decisions:")
    decision_df = pd.DataFrame(decisions)
    if not decision_df.empty:
        decision_cols = [
            "arm", "resolution", "passed", "role",
            "median_required", "cosine_median",
            "lower_required", "cosine_lower",
            "on_fail_action",
        ]
        decision_cols = [c for c in decision_cols if c in decision_df.columns]
        print(decision_df[decision_cols].to_string(index=False))

    n_fail = int((~decision_df["passed"].astype(bool)).sum()) if not decision_df.empty else 0
    n_pass = int(decision_df["passed"].astype(bool).sum()) if not decision_df.empty else 0
    print(f"[audit] {n_pass} pass / {n_fail} fail across {len(decision_df)} (arm, resolution) cells.")
    fallback = bool(decision.get("fallback_to_external_axis_only", False))
    if fallback:
        print("[audit] fallback_to_external_axis_only=TRUE — within-RRRM-2 projection is gated off.")
        print("[audit] Phase 4 requires --bypass-stability to proceed and is restricted to sensitivity analyses.")
    else:
        print("[audit] fallback_to_external_axis_only=FALSE — within-RRRM-2 projection is allowed for passing resolutions.")

    if discrepancies:
        return 1
    if args.strict and n_fail > 0 and not fallback:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
