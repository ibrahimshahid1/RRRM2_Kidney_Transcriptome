#!/usr/bin/env python3
"""
Self-audit: report which of the 14 remediation-guide fixes are detectable in
the current source tree and configuration. This is a static check (no pipeline
run required) and is intended to be a quick pre-submission sanity check.

A fix is reported `present` if all of its detection probes succeed, `partial`
if some probes succeed, `absent` if none. The script does NOT verify behavior
correctness; it verifies that the implementing artifacts exist with the
expected semantic markers.

Usage:
    python scripts/audit_fix_status.py
"""
from __future__ import annotations

from pathlib import Path
import argparse

REPO = Path(__file__).resolve().parents[1]


def grep_any(path: Path, patterns: list[str]) -> list[bool]:
    """Return True/False for each pattern, indicating whether it appears in the file."""
    if not path.exists():
        return [False] * len(patterns)
    try:
        text = path.read_text()
    except Exception:
        return [False] * len(patterns)
    return [p in text for p in patterns]


def all_or(probes: list[bool]) -> str:
    n = sum(probes)
    if n == len(probes) and n > 0:
        return "present"
    if n > 0:
        return "partial"
    return "absent"


# ---------------------------------------------------------------- fix probes
FIXES = [
    {
        "id": 1,
        "name": "Brown's-method / label-permutation edge p-value aggregation",
        "file": REPO / "src/statistics/full_regression.py",
        "probes": [
            "signed_t_sum_sqrt_degree",
            "empirical_signed_gene_pvalues",
            "p_empirical_signed",
        ],
    },
    {
        "id": 2,
        "name": "Hierarchical Benjamini-Bogomolov FDR",
        "file": REPO / "src/statistics/permutation_bootstrap.py",
        "probes": [
            "hierarchical-fdr",
            "load_hierarchical_gene_sets",
            "q_BB_two_stage",
            "family_selected_q05",
        ],
    },
    {
        "id": 3,
        "name": "RESULTS.md reconciliation (deferred to post-fix doc pass)",
        "file": REPO / "RESULTS.md",
        "probes": ["run_20260505"],
        "note": "deferred per remediation guide; status optional",
    },
    {
        "id": 4,
        "name": "External validation on OSD-102 / OSD-513 (+ context cohorts)",
        "file": REPO / "src/validation/external_replication.py",
        "probes": [
            "ALLOWED_STUDIES",
            "context_detected",
            "directional_claim",
        ],
    },
    {
        "id": 5,
        "name": "Silent-shifter null calibration",
        "file": REPO / "src/statistics/silent_shifters.py",
        "probes": [
            "expected_strict_under_independence",
            "hypergeom_p_overrepresentation",
            "observed_over_expected",
        ],
    },
    {
        "id": 6,
        "name": "Anchor pre-registration enforcement",
        "file": REPO / "src/networks/procrustes.py",
        "probes": [
            "load_configured_anchor_indices",
            "anchor_genes.yaml",
            "anchors_report.tsv",
        ],
    },
    {
        "id": 7,
        "name": "Config-code drift removal",
        "file": REPO / "tests/test_config_consistency.py",
        "probes": [
            "DEFAULT_TOPK",
        ],
    },
    {
        "id": 8,
        "name": "MuSiC reference-atlas sensitivity analysis",
        "file": REPO / "src/preprocessing/deconvolution_sensitivity.R",
        "probes": [
            "MuSiC",
        ],
    },
    {
        "id": 9,
        "name": "Permutation-statistic rename + appendix-tier full-pipeline run",
        "file": REPO / "src/statistics/permutation_bootstrap.py",
        "probes": [
            "perm_edge_sum",
            "edge_sum_node_rewiring",
        ],
    },
    {
        "id": 10,
        "name": "LIONESS pooling across cells",
        "file": REPO / "src/networks/lioness.py",
        "probes": [
            "pool",
        ],
    },
    {
        "id": 11,
        "name": "Alternative SSN methods (SSN/CSN/scLink)",
        "file": REPO / "src/networks/alternative_methods.py",
        "probes": [
            "ssn_delta_z",
            "compute_alternative_network",
        ],
    },
    {
        "id": 12,
        "name": "Alternative feature engineering for Phase 8",
        "file": REPO / "src/validation/enhanced_cv.py",
        "probes": [
            "sparse_edges",
            "edge_pca",
            "pathway_strength",
        ],
    },
    {
        "id": 13,
        "name": "Continuous classification targets",
        "file": REPO / "src/validation/continuous_target.py",
        "probes": [
            "Havcr1",
            "kidney_stress",
        ],
    },
    {
        "id": 14,
        "name": "Multi-study cohort pooling (LAR-Young)",
        "file": REPO / "src/validation/multi_study_pool.py",
        "probes": [
            "OSD-102",
            "OSD-771",
        ],
    },
]


def main() -> None:
    ap = argparse.ArgumentParser(description="Audit which remediation fixes are present in the source tree.")
    ap.add_argument("--verbose", action="store_true", help="Print per-probe details for each fix.")
    args = ap.parse_args()

    print("=" * 80)
    print(" RRRM-2 Kidney Transcriptome: Remediation Fix Status Audit")
    print("=" * 80)

    summary_counts = {"present": 0, "partial": 0, "absent": 0}
    rows = []
    for fix in FIXES:
        probes = grep_any(fix["file"], fix["probes"])
        status = all_or(probes)
        summary_counts[status] += 1
        rows.append((fix["id"], status, fix["name"], fix["file"], probes, fix.get("note", "")))

    width_id = 4
    width_status = 9
    width_name = 56
    print(f"\n{'ID':>{width_id}}  {'Status':<{width_status}}  {'Fix':<{width_name}}")
    print("-" * (width_id + 2 + width_status + 2 + width_name))
    for fid, status, name, fpath, probes, note in rows:
        line = f"{fid:>{width_id}}  {status:<{width_status}}  {name:<{width_name}}"
        if note:
            line += f"   ({note})"
        print(line)
        if args.verbose:
            try:
                rel = fpath.relative_to(REPO)
            except ValueError:
                rel = fpath
            print(f"     file: {rel}")
            for probe, hit in zip(FIXES[fid - 1]["probes"], probes):
                marker = "+" if hit else "-"
                print(f"       {marker} probe: {probe!r}")

    print()
    print(f"Summary: present={summary_counts['present']}  "
          f"partial={summary_counts['partial']}  "
          f"absent={summary_counts['absent']}  "
          f"total={len(FIXES)}")

    test_dir = REPO / "tests"
    if test_dir.exists():
        n_tests = len(list(test_dir.glob("test_*.py")))
        print(f"Test files in tests/: {n_tests}")
    print()


if __name__ == "__main__":
    main()
