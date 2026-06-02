#!/usr/bin/env python3
"""Compact mechanism-triangulation summary for v11."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


RUN_ROOT = Path("data/results/run_20260526_v11_dct1_phospho_mediation")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {"status": "missing", "path": str(path)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=RUN_ROOT)
    args = ap.parse_args()
    out_dir = args.run_root / "perturbation"
    out_dir.mkdir(parents=True, exist_ok=True)

    lowk = read_json(out_dir / "lowk_alignment_verdict.json")
    occ = read_json(args.run_root / "h2_occupancy/h2_occupancy_verdict.json")
    spatial = read_json(args.run_root / "spatial_reference/visium_dct_transport_verdict.json")
    pxd = read_json(args.run_root / "h2_pxd/h2_pxd001729_ddavp_antialignment_verdict.json")
    klhl = args.run_root / "h2_klhl3/h2_klhl3_cul3_interpretation.md"

    rows = [
        {
            "mechanism_branch": "lowk_dct_wnk_reference",
            "analysis": "GSE228367 low-K DCT-enriched pseudobulk anti-alignment",
            "status": lowk.get("status"),
            "headline_result": lowk.get("classification"),
            "paper_use": "primary_candidate" if lowk.get("promote_to_primary_result") else "mechanism_triage_or_negative",
            "boundary": lowk.get("boundary"),
        },
        {
            "mechanism_branch": "parent_protein_normalized_phosphorylation",
            "analysis": "OSD-462 parent-protein-normalized DCT-subtype-prior enrichment",
            "status": occ.get("status"),
            "headline_result": f"primary={occ.get('primary_record', {})}",
            "paper_use": "robustness_upgrade" if occ.get("supports_parent_normalized_dct_prior_enrichment") else "bounded_robustness",
            "boundary": occ.get("boundary"),
        },
        {
            "mechanism_branch": "injury_repair_spatial_context",
            "analysis": "GSE269622 IRI Visium DCT transport score in DCT-adjacent niches",
            "status": spatial.get("status"),
            "headline_result": (
                "late_iri_dct_adjacent_suppression"
                if spatial.get("late_iri_dct_adjacent_transport_suppression_present")
                else "late_iri_no_clear_dct_adjacent_suppression"
            ),
            "paper_use": "spatial_prediction_context",
            "boundary": spatial.get("boundary"),
        },
        {
            "mechanism_branch": "vasopressin_camp",
            "analysis": "PXD001729 mpkDCT dDAVP phosphoproteomics comparison",
            "status": pxd.get("status", "available" if pxd else "missing"),
            "headline_result": pxd.get("classification", pxd.get("verdict", "see h2_pxd outputs")),
            "paper_use": "secondary_negative_or_context",
            "boundary": "cultured DCT-lineage phosphoproteomics; not spaceflight tissue",
        },
        {
            "mechanism_branch": "wnk_turnover_klhl3_cul3",
            "analysis": "OSD-462 KLHL3/CUL3/WNK coverage check",
            "status": "available" if klhl.exists() else "missing",
            "headline_result": "unresolved_from_current_public_omics" if klhl.exists() else "missing",
            "paper_use": "future_experiment_rationale",
            "boundary": "no ubiquitinomics or DCT-isolated WNK protein turnover assay",
        },
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(out_dir / "perturbation_triage_summary.tsv", sep="\t", index=False)

    verdict = {
        "status": "complete",
            "best_current_upgrade": (
                "lowK_primary_candidate"
                if lowk.get("promote_to_primary_result")
                else "parent_protein_normalization_and_triage"
                if occ.get("supports_parent_normalized_dct_prior_enrichment")
                else "bounded_negative_mechanism_triage"
            ),
        "rows": rows,
    }
    (out_dir / "perturbation_triage_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(f"[perturbation-triage] complete: {out_dir}")


if __name__ == "__main__":
    main()
