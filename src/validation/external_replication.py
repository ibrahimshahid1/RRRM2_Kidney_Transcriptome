# src/validation/external_replication.py
"""
Protocol-guarded independent external cohort analysis.

OSD-102 is the primary LAR-Young-like replication partner. OSD-513 is secondary
and limited to sex-robustness/sex-stratification checks. OSD-163 and OSD-253 are
context-mapping cohorts for the biology-first remodeling panel; they are not
used as strict one-to-one replication cohorts for RRRM-2 gene claims. OSD-568 is
explicitly excluded from validation claims in this remediation pass.

This module does not require ComBat-seq. Each external cohort is analyzed
independently and compared against pre-registered direction, q-value, and
pathway criteria. Multi-study pooling is handled separately by
src.validation.multi_study_pool.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT

ALLOWED_STUDIES = {"OSD-102", "OSD-513", "OSD-163", "OSD-253"}
EXCLUDED_STUDIES = {"OSD-568"}
EXPECTED_ANALYSIS_BY_STUDY = {
    "OSD-102": "lar_young",
    "OSD-513": "sex_robustness",
    "OSD-163": "cross_strain_context",
    "OSD-253": "rr7_context",
}
REQUIRED_HYPOTHESIS_COLUMNS = {
    "analysis",
    "study",
    "hypothesis_type",
    "feature",
    "expected_direction",
    "discovery_effect",
    "q_threshold",
    "claim_boundary",
}
PLACEHOLDER_FEATURES = {
    "pre_registered_lar_young_pathways",
    "sex_consistency",
    "tbd",
    "todo",
}


def validate_study_scope(study: str, analysis: str) -> None:
    if study in EXCLUDED_STUDIES:
        raise ValueError("OSD-568 is excluded from future-validation claims and must not be used here.")
    if study not in ALLOWED_STUDIES:
        raise ValueError(f"{study} is not in the approved independent-replication scope: {sorted(ALLOWED_STUDIES)}")
    expected = EXPECTED_ANALYSIS_BY_STUDY[study]
    if analysis != expected:
        raise ValueError(f"{study} is restricted to analysis={expected}.")


def write_template_protocol(protocol_dir: Path) -> None:
    protocol_dir.mkdir(parents=True, exist_ok=True)
    protocol = {
        "version": "2026-05-09",
        "independent_replication": {
            "primary": {
                "study": "OSD-102",
                "scope": "LAR-Young-like flight/ground findings and pathways only",
                "no_combat_seq": True,
                "q_threshold": 0.10,
                "criteria": [
                    "same direction for pre-registered genes/pathways",
                    "gene/pathway q-value threshold fixed in protocol before loading OSD-102 results",
                    "pathway replication judged on fixed pathway list and direction",
                ],
            },
            "secondary": {
                "study": "OSD-513",
                "scope": "sex-stratification and sex-robustness only",
                "no_combat_seq": True,
                "q_threshold": 0.10,
            },
            "context_mapping": {
                "studies": ["OSD-163", "OSD-253"],
                "scope": "fixed lipid/ECM/tubular remodeling panel; non-directional context mapping, not strict replication",
                "no_combat_seq": True,
                "q_threshold": 0.10,
            },
            "excluded": ["OSD-568"],
        },
        "non_claims": [
            "ISS-T validation",
            "old-arm validation",
            "global cohort expansion",
            "classifier validation",
        ],
    }
    (protocol_dir / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    directional_rows = [
        {
            "analysis": "lar_young",
            "study": "OSD-102",
            "hypothesis_type": "pathway",
            "feature": "PPAR_signaling",
            "expected_direction": "positive_edge_sum_enrichment",
            "discovery_effect": 1.0,
            "q_threshold": 0.10,
            "claim_boundary": "primary pathway-level replication only",
        },
        {
            "analysis": "lar_young",
            "study": "OSD-102",
            "hypothesis_type": "pathway",
            "feature": "translation_machinery",
            "expected_direction": "positive_edge_sum_enrichment",
            "discovery_effect": 1.0,
            "q_threshold": 0.10,
            "claim_boundary": "primary pathway-level replication only",
        },
        {
            "analysis": "sex_robustness",
            "study": "OSD-513",
            "hypothesis_type": "pathway",
            "feature": "PPAR_signaling",
            "expected_direction": "same_direction_as_OSD102_or_RRRM2",
            "discovery_effect": 1.0,
            "q_threshold": 0.10,
            "claim_boundary": "secondary sex-robustness check only",
        },
        {
            "analysis": "sex_robustness",
            "study": "OSD-513",
            "hypothesis_type": "pathway",
            "feature": "translation_machinery",
            "expected_direction": "same_direction_as_OSD102_or_RRRM2",
            "discovery_effect": 1.0,
            "q_threshold": 0.10,
            "claim_boundary": "secondary sex-robustness check only",
        },
    ]
    context_features = [
        "PPAR_signaling",
        "cholesterol_biosynthesis",
        "ECM_remodeling",
        "EMT_fibrosis",
        "tubular_ion_transport",
        "TGF_beta_Wnt",
        "oxidative_stress",
        "translation_machinery",
    ]
    context_rows = []
    for study, analysis in [("OSD-163", "cross_strain_context"), ("OSD-253", "rr7_context")]:
        for feature in context_features:
            context_rows.append({
                "analysis": analysis,
                "study": study,
                "hypothesis_type": "pathway",
                "feature": feature,
                "expected_direction": "context_dependent_screen",
                "discovery_effect": 0.0,
                "q_threshold": 0.10,
                "claim_boundary": "cross-cohort context mapping only; not a directional replication claim",
            })
    hypothesis = pd.DataFrame(directional_rows + context_rows)
    hypothesis.to_csv(protocol_dir / "hypothesis_registry.tsv", sep="\t", index=False)


def validate_hypothesis_registry(hypotheses: pd.DataFrame) -> pd.DataFrame:
    """Require concrete pre-registered features, directions, effects, and thresholds."""
    missing = REQUIRED_HYPOTHESIS_COLUMNS - set(hypotheses.columns)
    if missing:
        raise ValueError(f"Hypothesis registry missing required columns: {sorted(missing)}")
    registry = hypotheses.copy()
    registry["feature"] = registry["feature"].astype(str).str.strip()
    bad_feature = registry["feature"].str.lower().isin(PLACEHOLDER_FEATURES)
    bad_feature |= registry["feature"].str.lower().str.startswith("pre_registered_")
    if bad_feature.any():
        bad = sorted(registry.loc[bad_feature, "feature"].unique())
        raise ValueError(f"Hypothesis registry contains placeholder features: {bad}")
    registry["discovery_effect"] = pd.to_numeric(registry["discovery_effect"], errors="coerce")
    registry["q_threshold"] = pd.to_numeric(registry["q_threshold"], errors="coerce")
    if registry["discovery_effect"].isna().any():
        raise ValueError("Hypothesis registry discovery_effect must be numeric for every row")
    if registry["q_threshold"].isna().any() or (registry["q_threshold"] <= 0).any() or (registry["q_threshold"] >= 1).any():
        raise ValueError("Hypothesis registry q_threshold must be numeric and between 0 and 1")
    if registry.duplicated(["analysis", "study", "feature"]).any():
        dup = registry.loc[registry.duplicated(["analysis", "study", "feature"], keep=False), "feature"].tolist()
        raise ValueError(f"Hypothesis registry has duplicate analysis/study/feature rows: {dup}")
    return registry


def require_protocol(protocol_dir: Path) -> tuple[dict, pd.DataFrame]:
    protocol_path = protocol_dir / "protocol.json"
    hypothesis_path = protocol_dir / "hypothesis_registry.tsv"
    if not protocol_path.exists() or not hypothesis_path.exists():
        raise FileNotFoundError(
            f"External replication requires a committed protocol at {protocol_dir}. "
            "Run with --write-template-protocol before loading external data."
        )
    protocol = json.loads(protocol_path.read_text())
    hypotheses = validate_hypothesis_registry(pd.read_csv(hypothesis_path, sep="\t"))
    return protocol, hypotheses


def registry_for_analysis(hypotheses: pd.DataFrame, study: str, analysis: str) -> pd.DataFrame:
    registry = hypotheses[
        (hypotheses["study"].astype(str) == study) &
        (hypotheses["analysis"].astype(str) == analysis)
    ].copy()
    if registry.empty:
        raise ValueError(f"No registered hypotheses for study={study}, analysis={analysis}")
    return registry


def validate_discovery_table_against_registry(discovery: pd.DataFrame, registry: pd.DataFrame) -> None:
    required = {"feature", "effect"}
    missing = required - set(discovery.columns)
    if missing:
        raise ValueError(f"Discovery table missing columns: {sorted(missing)}")
    merged = registry[["feature", "discovery_effect"]].merge(
        discovery[["feature", "effect"]],
        on="feature",
        how="left",
    )
    if merged["effect"].isna().any():
        missing_features = merged.loc[merged["effect"].isna(), "feature"].tolist()
        raise ValueError(f"Discovery table is missing registered features: {missing_features}")
    directional = np.sign(merged["discovery_effect"]) != 0
    if directional.any() and not (
        np.sign(merged.loc[directional, "discovery_effect"]) == np.sign(merged.loc[directional, "effect"])
    ).all():
        raise ValueError("Discovery table effect signs disagree with the committed hypothesis registry")


def compare_directional_replication(
    discovery: pd.DataFrame,
    external: pd.DataFrame,
    q_threshold: float | None = None,
) -> pd.DataFrame:
    """Compare pre-registered genes/pathways by direction and external q-value."""
    required = {"feature", "effect"}
    missing_discovery = required - set(discovery.columns)
    missing_external = (required | {"q_value"}) - set(external.columns)
    if missing_discovery:
        raise ValueError(f"Discovery table missing columns: {sorted(missing_discovery)}")
    if missing_external:
        raise ValueError(f"External table missing columns: {sorted(missing_external)}")
    if q_threshold is None and "q_threshold" not in discovery.columns:
        raise ValueError("q_threshold must be supplied either as an argument or discovery column")

    missing_features = sorted(set(discovery["feature"]) - set(external["feature"]))
    if missing_features:
        raise ValueError(f"External table is missing registered features: {missing_features}")

    discovery_cols = ["feature", "effect"]
    if "q_threshold" in discovery.columns:
        discovery_cols.append("q_threshold")
    merged = discovery[discovery_cols].rename(columns={"effect": "discovery_effect"}).merge(
        external[["feature", "effect", "q_value"]].rename(columns={"effect": "external_effect"}),
        on="feature",
        how="inner",
    )
    if q_threshold is not None:
        merged["q_threshold"] = float(q_threshold)
    merged["directional_claim"] = np.sign(merged["discovery_effect"]) != 0
    same_direction = np.sign(merged["discovery_effect"]) == np.sign(merged["external_effect"])
    merged["same_direction"] = same_direction.where(merged["directional_claim"], pd.NA)
    merged["external_pass_q"] = merged["q_value"] <= merged["q_threshold"]
    merged["context_detected"] = merged["external_pass_q"]
    merged["replicated"] = merged["directional_claim"] & merged["same_direction"].fillna(False) & merged["external_pass_q"]
    return merged


def main() -> None:
    ap = argparse.ArgumentParser(description="Protocol-guarded independent OSD replication")
    ap.add_argument("--study", required=True, choices=sorted(ALLOWED_STUDIES | EXCLUDED_STUDIES))
    ap.add_argument("--analysis", required=True,
                    choices=["lar_young", "sex_robustness", "cross_strain_context", "rr7_context"])
    ap.add_argument("--protocol_dir", default=str(REPO_ROOT / "docs/external_replication_protocol"))
    ap.add_argument("--write-template-protocol", action="store_true")
    ap.add_argument("--discovery_table", default="",
                    help="Pre-registered RRRM-2 features with columns feature,effect")
    ap.add_argument("--external_table", default="",
                    help="External cohort results with columns feature,effect,q_value")
    ap.add_argument("--q_threshold", type=float, default=None,
                    help="Optional override; must match the committed registry threshold if supplied.")
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/external_replication"))
    args = ap.parse_args()

    protocol_dir = Path(args.protocol_dir)
    if args.write_template_protocol:
        write_template_protocol(protocol_dir)
        print(f"[OK] Wrote protocol template to {protocol_dir}")

    validate_study_scope(args.study, args.analysis)
    protocol, hypotheses = require_protocol(protocol_dir)
    registry = registry_for_analysis(hypotheses, args.study, args.analysis)
    if args.q_threshold is not None and not np.allclose(registry["q_threshold"].values, args.q_threshold):
        raise ValueError("CLI q_threshold does not match committed hypothesis registry thresholds")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    discovery = registry.rename(columns={"discovery_effect": "effect"})[
        ["feature", "effect", "q_threshold", "expected_direction", "hypothesis_type", "claim_boundary"]
    ].copy()
    if args.discovery_table:
        supplied_discovery = pd.read_csv(args.discovery_table, sep=None, engine="python")
        validate_discovery_table_against_registry(supplied_discovery, registry)

    if args.external_table:
        external = pd.read_csv(args.external_table, sep=None, engine="python")
        result = compare_directional_replication(discovery, external, args.q_threshold)
        result["study"] = args.study
        result["analysis"] = args.analysis
        result.to_csv(outdir / f"{args.study}_{args.analysis}_replication.tsv", sep="\t", index=False)
        print(f"[OK] Wrote replication comparison to {outdir}")
    else:
        manifest = {
            "study": args.study,
            "analysis": args.analysis,
            "status": "protocol_validated_no_external_tables_loaded",
            "claim_boundary": "independent replication only; no pooled cohort claim",
            "protocol_version": protocol.get("version"),
        }
        (outdir / f"{args.study}_{args.analysis}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        print(f"[OK] Protocol guard passed. Manifest written to {outdir}")


if __name__ == "__main__":
    main()
