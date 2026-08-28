"""Tests for read-only v13 inference reporting."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import yaml

from src.v13.reporting import (
    build_claim_decision_summary,
    build_leading_parent_gene_matrix,
    build_primary_compartment_table,
    build_report,
    build_robustness_summary,
    load_inference_artifacts,
)


PRIMARY_SCORE = "median_negative_site_effect"


def _write_minimal_run(root: Path) -> Path:
    root.mkdir()
    config = {
        "set_test": {
            "primary_family": ["DCT2_CNT_transition", "ASDN"],
            "minimum_observable_genes": 2,
        },
        "canonical_exclusions": {"primary": [], "strict": ["Wnk1"]},
        "comparator_sets": ["proximal_tubule"],
    }
    config_path = root / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))

    membership = pd.DataFrame(
        [
            ("A", "ASDN", "transport"),
            ("B", "ASDN", "transport"),
            ("C", "ASDN", "sensing"),
            ("D", "DCT2_CNT_transition", ""),
            ("E", "DCT2_CNT_transition", ""),
            ("D", "strict_DCT2_peaked", ""),
            ("A", "proximal_tubule", ""),
            ("F", "proximal_tubule", ""),
        ],
        columns=["gene_symbol", "gene_set", "category"],
    )
    membership["definition_source"] = "test"
    membership["status"] = "defined"
    membership["final_for_testing"] = True
    membership["analysis_role"] = "test"
    membership.to_csv(
        root / "gene_set_membership_frozen_copy.tsv", sep="\t", index=False
    )

    gene_rows = []
    profile_genes = {
        "primary": ["A", "B", "C", "D", "F"],
        "exclude_multicompartment_broad_expression": ["A", "B", "D", "F"],
        "exclude_non_distal_or_broad_expression": ["A", "D", "F"],
        "parent_protein_subtracted": ["A", "D", "F"],
    }
    z_values = {"A": 2.0, "B": 1.0, "C": -0.5, "D": 0.25, "F": 0.0}
    for profile, genes in profile_genes.items():
        for gene in genes:
            gene_rows.append(
                {
                    "profile": profile,
                    "gene_symbol": gene,
                    "gene_score": PRIMARY_SCORE,
                    "observed_raw_score": z_values[gene] / 10,
                    "null_mean": 0,
                    "null_sd": 0.1,
                    "observed_gene_z": z_values[gene],
                    "empirical_p_greater": 0.1,
                    "n_null_valid": 31,
                    "n_sites": 2,
                    "median_log2_signal": 3,
                    "mean_missing_fraction": 0,
                    "fixed_null_universe_eligible": True,
                }
            )
    pd.DataFrame(gene_rows).to_csv(
        root / "parent_gene_null_calibration.tsv", sep="\t", index=False
    )

    set_rows = [
        {
            "profile": "primary",
            "exclusion": "primary",
            "gene_score": PRIMARY_SCORE,
            "set_test": "competitive",
            "gene_set": "ASDN",
            "set_role": "primary",
            "n_observable_genes": 3,
            "n_eligible_background_genes": 2,
            "observed_statistic": 0.7,
            "null_mean": 0,
            "null_sd": 0.2,
            "null_ci_low": -0.4,
            "null_ci_high": 0.4,
            "empirical_p_greater": 0.03,
            "n_null_valid": 31,
            "maxT_fwer": 0.04,
            "comparator_bh_q": float("nan"),
        },
        {
            "profile": "primary",
            "exclusion": "primary",
            "gene_score": PRIMARY_SCORE,
            "set_test": "competitive",
            "gene_set": "proximal_tubule",
            "set_role": "kidney_comparator",
            "n_observable_genes": 2,
            "n_eligible_background_genes": 3,
            "observed_statistic": 0.1,
            "null_mean": 0,
            "null_sd": 0.2,
            "null_ci_low": -0.4,
            "null_ci_high": 0.4,
            "empirical_p_greater": 0.4,
            "n_null_valid": 31,
            "maxT_fwer": float("nan"),
            "comparator_bh_q": 0.4,
        },
        {
            "profile": "exclude_multicompartment_broad_expression",
            "exclusion": "primary",
            "gene_score": PRIMARY_SCORE,
            "set_test": "competitive",
            "gene_set": "ASDN",
            "set_role": "primary",
            "n_observable_genes": 2,
            "n_eligible_background_genes": 2,
            "observed_statistic": 0.8,
            "null_mean": 0,
            "null_sd": 0.2,
            "null_ci_low": -0.4,
            "null_ci_high": 0.4,
            "empirical_p_greater": 0.03,
            "n_null_valid": 31,
            "maxT_fwer": 0.04,
            "comparator_bh_q": float("nan"),
        },
    ]
    pd.DataFrame(set_rows).to_csv(
        root / "set_level_permutation_inference.tsv", sep="\t", index=False
    )

    pd.DataFrame(
        [
            {
                "gene_set": "DCT2_CNT_transition",
                "claim_gate_status": "non_evaluable",
                "gate_minimum_observable_genes": "fail",
                "gate_positive_competitive_effect": "non_evaluable",
                "n_leave_one_gene_out_rows": 0,
                "gate_leave_one_gene_out_completeness": "non_evaluable",
                "gate_all_leave_one_gene_out_positive": "non_evaluable",
                "gate_centered_vs_uncentered_direction": "non_evaluable",
                "gate_scaled_vs_summed_direction": "non_evaluable",
                "gate_normalization_direction_concordance": "non_evaluable",
                "gate_multicompartment_broad_expression_exclusion_direction": "non_evaluable",
                "gate_no_unrelated_compartment_equal_or_stronger": "non_evaluable",
                "boundary": "Parent-protein annotation only.",
            },
            {
                "gene_set": "ASDN",
                "claim_gate_status": "fail",
                "gate_minimum_observable_genes": "pass",
                "gate_positive_competitive_effect": "pass",
                "n_leave_one_gene_out_rows": 3,
                "gate_leave_one_gene_out_completeness": "pass",
                "gate_all_leave_one_gene_out_positive": "pass",
                "gate_centered_vs_uncentered_direction": "pass",
                "gate_scaled_vs_summed_direction": "pass",
                "gate_normalization_direction_concordance": "pass",
                "gate_multicompartment_broad_expression_exclusion_direction": "pass",
                "gate_no_unrelated_compartment_equal_or_stronger": "fail",
                "boundary": "Parent-protein annotation only.",
            },
        ]
    ).to_csv(root / "claim_gates.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "claim_tier": "neither",
                "DCT2_CNT_gate_status": "non_evaluable",
                "ASDN_gate_status": "fail",
                "statistical_beyond_axis_support": "False",
                "statistical_dct2_title_eligible": "False",
                "statistical_asdn_claim_eligible": "False",
                "publication_promotion_status": "blocked_by_design_or_provenance",
                "publication_design_and_provenance_gate": "fail",
                "publication_blockers": "reporter-tag aliasing",
                "dct2_title_permitted": "False",
                "asdn_beyond_axis_claim_permitted": "False",
                "publication_beyond_axis_claim_permitted": "False",
                "n_stage0_qualified_anchor_phosphoforms": 0,
                "rationale": (
                    "DCT2 was coverage-non-evaluable and ASDN failed a gate."
                ),
            }
        ]
    ).to_csv(root / "claim_tier.tsv", sep="\t", index=False)
    manifest = {
        "analysis_id": "test",
        "config": str(config_path),
        "primary_gene_score": PRIMARY_SCORE,
        "condition_reporter_position_confounded": True,
        "profiles": [
            {"name": "primary"},
            {"name": "exclude_multicompartment_broad_expression"},
            {"name": "exclude_non_distal_or_broad_expression"},
            {"name": "parent_protein_subtracted"},
        ],
        "permutation": {"mode": "smoke", "n_assignments_run": 32},
        "inputs": {},
        "normalization_equivalence_audit": {
            "independent_robustness_evidence": False
        },
    }
    (root / "manifest.json").write_text(json.dumps(manifest))
    (root / "run_summary.json").write_text(json.dumps({"smoke": True}))
    return root


def test_reporting_marks_missing_primary_set_non_evaluable(tmp_path: Path):
    artifacts = load_inference_artifacts(_write_minimal_run(tmp_path / "run"))
    claim = build_claim_decision_summary(artifacts).set_index("gene_set")
    assert claim.loc["DCT2_CNT_transition", "n_observable_genes"] == 1
    assert claim.loc["DCT2_CNT_transition", "claim_gate_status"] == "non_evaluable"
    assert (
        claim.loc[
            "DCT2_CNT_transition", "statistical_evaluation_status"
        ]
        == "coverage_non_evaluable"
    )
    assert not bool(claim.loc["DCT2_CNT_transition", "claim_permitted"])
    assert (
        claim.loc["ASDN", "statistical_evaluation_status"]
        == "evaluated_gate_failure"
    )
    assert not bool(claim.loc["ASDN", "statistical_claim_eligible"])
    assert not bool(claim.loc["ASDN", "claim_permitted"])
    assert claim.loc["ASDN", "gate_scaled_vs_summed_direction"] == "pass"
    assert (
        claim.loc["ASDN", "gate_leave_one_gene_out_completeness"]
        == "pass"
    )
    assert (
        claim.loc["ASDN", "publication_promotion_status"]
        == "blocked_by_design_or_provenance"
    )

    forest = build_primary_compartment_table(artifacts).set_index("gene_set")
    assert forest.loc["DCT2_CNT_transition", "evaluation_status"] == "non_evaluable"
    assert forest.loc["ASDN", "multiplicity_adjustment"] == "maxT_FWER"
    assert forest.loc["proximal_tubule", "multiplicity_adjustment"] == "BH_q"


def test_leading_matrix_uses_run_internal_broad_profile(tmp_path: Path):
    artifacts = load_inference_artifacts(_write_minimal_run(tmp_path / "run"))
    leading = build_leading_parent_gene_matrix(artifacts).set_index("gene_symbol")
    assert bool(leading.loc["C", "broadly_expressed"])
    assert not bool(leading.loc["A", "broadly_expressed"])
    assert not bool(leading.loc["B", "broadly_expressed"])
    assert (
        leading.loc["C", "broad_flag_source"]
        == "inference_profile_difference"
    )
    assert pd.notna(leading.loc["A", "parent_protein_observed_gene_z"])
    assert pd.isna(leading.loc["B", "parent_protein_observed_gene_z"])


def test_full_report_writes_tables_markdown_and_plots(tmp_path: Path):
    run = _write_minimal_run(tmp_path / "run")
    outputs = build_report(run, tmp_path / "report")
    assert all(path.exists() and path.stat().st_size > 0 for path in outputs.values())
    assert outputs["compartment_plot_pdf"].suffix == ".pdf"
    assert outputs["leading_matrix_pdf"].suffix == ".pdf"
    robustness = build_robustness_summary(load_inference_artifacts(run))
    missing = robustness[
        robustness["gene_set"].eq("DCT2_CNT_transition")
        & robustness["profile"].eq("primary")
        & robustness["exclusion"].eq("primary")
    ]
    assert missing.iloc[0]["evaluation_status"] == "non_evaluable"
    report_manifest = json.loads(outputs["report_manifest"].read_text())
    assert report_manifest["analysis_recomputed"] is False
