"""Regression contracts for the corrected v13 inference and claim layer."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.v13.continuous_phospho_inference import (
    PRIMARY_SCORE,
    AnalysisProfile,
    PreparedPhosphoData,
    _claim_gate_table,
    _first_pass_calibration,
    _minimum_localization_score,
    _second_pass_inference,
    audit_manual_anchor_exclusion,
    build_gene_layout,
    enumerate_balanced_labels,
    matched_observability_tests,
    prepare_osd462_phosphosites,
)


def _sample_metadata() -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for plex in ("Samp1-5", "Samp6-10"):
        for index in range(4):
            condition = "FL" if index < 2 else "GC"
            rows.append(
                {
                    "column": f"{plex}|{condition}|{plex}-{index}",
                    "plex": plex,
                    "observed_condition": condition,
                    "sample_id": f"{plex}-{index}",
                }
            )
    return pd.DataFrame(rows)


def _claim_config(
    *,
    reporter_tag_confounded: bool = False,
    require_anchor: bool = False,
) -> dict:
    return {
        "contrast": {
            "condition_reporter_position_confounded": reporter_tag_confounded,
        },
        "set_test": {
            "minimum_observable_genes": 8,
            "primary_family": ["DCT2_CNT_transition", "ASDN"],
        },
        "comparator_sets": ["proximal_tubule"],
        "robustness": {
            "require_normalization_direction_concordance": True,
            "broad_expression_exclusion": True,
        },
        "publication_promotion_gate": {
            "qualified_canonical_anchor_required_for_extension_claim": (
                require_anchor
            ),
        },
        "claim_gates": {"alpha": 0.05},
    }


def _asdn_claim_inputs(
    *,
    primary_statistic: float = 1.0,
    max_t: float = 0.01,
    uncentered_statistic: float = 0.9,
    summed_statistic: float | None = 0.8,
    broad_statistic: float = 0.7,
    strict_statistic: float = 0.6,
    comparator_statistic: float = 0.1,
    loo_rows: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = {
        "gene_score": PRIMARY_SCORE,
        "set_test": "competitive",
        "set_role": "primary",
        "gene_set": "ASDN",
        "n_observable_genes": 8,
        "maxT_fwer": np.nan,
    }
    rows = [
        {
            **common,
            "profile": "primary",
            "exclusion": "primary",
            "observed_statistic": primary_statistic,
            "maxT_fwer": max_t,
        },
        {
            **common,
            "profile": "primary",
            "exclusion": "strict",
            "observed_statistic": strict_statistic,
        },
        {
            **common,
            "profile": "official_scaled_uncentered",
            "exclusion": "primary",
            "observed_statistic": uncentered_statistic,
        },
        {
            **common,
            "profile": "exclude_multicompartment_broad_expression",
            "exclusion": "primary",
            "observed_statistic": broad_statistic,
        },
        {
            "profile": "primary",
            "exclusion": "primary",
            "gene_score": PRIMARY_SCORE,
            "set_test": "competitive",
            "set_role": "kidney_comparator",
            "gene_set": "proximal_tubule",
            "n_observable_genes": 20,
            "observed_statistic": comparator_statistic,
            "maxT_fwer": np.nan,
        },
    ]
    if summed_statistic is not None:
        rows.append(
            {
                **common,
                "profile": "signal_to_noise_sum_centered",
                "exclusion": "primary",
                "observed_statistic": summed_statistic,
            }
        )
    loo = pd.DataFrame(
        [
            {
                "profile": "primary",
                "exclusion": "primary",
                "gene_score": PRIMARY_SCORE,
                "gene_set": "ASDN",
                "removed_gene": f"A{index}",
                "observed_competitive_statistic": 0.4,
            }
            for index in range(loo_rows)
        ]
    )
    return pd.DataFrame(rows), loo


def test_claim_gate_requires_both_uncentered_and_summed_normalization_directions():
    set_results, loo = _asdn_claim_inputs(summed_statistic=-0.2)
    gates, tier = _claim_gate_table(
        set_results,
        loo,
        _claim_config(),
        observable_counts={"DCT2_CNT_transition": 5, "ASDN": 8},
    )
    asdn = gates.set_index("gene_set").loc["ASDN"]
    assert asdn["gate_centered_vs_uncentered_direction"] == "pass"
    assert asdn["gate_scaled_vs_summed_direction"] == "fail"
    assert asdn["gate_normalization_direction_concordance"] == "fail"
    assert asdn["claim_gate_status"] == "fail"
    assert tier.iloc[0]["claim_tier"] == "neither"

    without_summed, loo = _asdn_claim_inputs(summed_statistic=None)
    gates, _ = _claim_gate_table(
        without_summed,
        loo,
        _claim_config(),
        observable_counts={"DCT2_CNT_transition": 5, "ASDN": 8},
    )
    asdn = gates.set_index("gene_set").loc["ASDN"]
    assert asdn["gate_scaled_vs_summed_direction"] == "non_evaluable"
    assert asdn["claim_gate_status"] == "non_evaluable"


def test_claim_gate_fails_when_leave_one_gene_out_rows_are_incomplete():
    set_results, loo = _asdn_claim_inputs(loo_rows=7)
    gates, tier = _claim_gate_table(
        set_results,
        loo,
        _claim_config(),
        observable_counts={"DCT2_CNT_transition": 5, "ASDN": 8},
    )
    asdn = gates.set_index("gene_set").loc["ASDN"]
    assert asdn["n_leave_one_gene_out_rows"] == 7
    assert asdn["gate_leave_one_gene_out_completeness"] == "fail"
    assert asdn["gate_all_leave_one_gene_out_positive"] == "fail"
    assert asdn["claim_gate_status"] == "fail"
    assert tier.iloc[0]["statistical_beyond_axis_support"] == np.bool_(False)


def test_dct2_non_evaluable_and_asdn_failure_supports_no_beyond_axis_claim():
    set_results, loo = _asdn_claim_inputs(max_t=0.2)
    gates, tier = _claim_gate_table(
        set_results,
        loo,
        _claim_config(),
        observable_counts={"DCT2_CNT_transition": 5, "ASDN": 8},
    )
    by_set = gates.set_index("gene_set")
    assert by_set.loc["DCT2_CNT_transition", "claim_gate_status"] == "non_evaluable"
    assert by_set.loc["ASDN", "claim_gate_status"] == "fail"
    decision = tier.iloc[0]
    assert decision["claim_tier"] == "neither"
    assert not bool(decision["statistical_beyond_axis_support"])
    assert not bool(decision["publication_beyond_axis_claim_permitted"])
    assert "no beyond-axis claim is supported" in decision["rationale"]
    assert "does not establish a null DCT2/CNT effect" in decision["rationale"]


def test_publication_promotion_is_blocked_by_tag_aliasing_and_zero_anchors():
    set_results, loo = _asdn_claim_inputs()
    _, tier = _claim_gate_table(
        set_results,
        loo,
        _claim_config(reporter_tag_confounded=True, require_anchor=True),
        observable_counts={"DCT2_CNT_transition": 5, "ASDN": 8},
        promotion_context={"n_stage0_qualified_anchor_phosphoforms": 0},
    )
    decision = tier.iloc[0]
    assert decision["claim_tier"] == "ASDN_only"
    assert bool(decision["statistical_beyond_axis_support"])
    assert decision["publication_promotion_status"] == (
        "blocked_by_design_or_provenance"
    )
    assert decision["publication_design_and_provenance_gate"] == "fail"
    assert "reporter-tag position blocks" in decision["publication_blockers"]
    assert "no isolated Stage-0-qualified" in decision["publication_blockers"]
    assert not bool(decision["asdn_beyond_axis_claim_permitted"])
    assert not bool(decision["publication_beyond_axis_claim_permitted"])


def test_max_t_uses_the_joint_null_across_multiple_primary_sets():
    samples = _sample_metadata()
    observed_flight = samples["observed_condition"].eq("FL").to_numpy()
    genes = (
        [f"A{index}" for index in range(3)]
        + [f"C{index}" for index in range(3)]
        + [f"B{index}" for index in range(10)]
    )
    rng = np.random.default_rng(0)
    values = rng.normal(0, 0.4, (len(genes), len(samples)))
    values[:3, observed_flight] -= 0.55
    values[3:6, observed_flight] -= 0.50
    metadata = pd.DataFrame(
        {
            "gene_symbol": genes,
            "mean_log2_signal_uncentered": 8.0,
            "missing_fraction": 0.0,
        }
    )
    profile = AnalysisProfile(
        name="primary",
        min_localization_score=19,
        include_composite=False,
        require_singly_modified=True,
        min_finite_total=8,
        min_finite_per_plex=4,
        deduplication="accession_modified_peptide_phosphoform",
    )
    data = PreparedPhosphoData(
        values=values,
        site_metadata=metadata,
        sample_metadata=samples,
        filter_audit=pd.DataFrame(),
        profile=profile,
    )
    design = enumerate_balanced_labels(samples, mode="exact")
    layout = build_gene_layout(metadata)
    calibration = _first_pass_calibration(
        data, design, layout, chunk_size=36, min_per_group=2
    )
    cfg = {
        "set_test": {
            "minimum_observable_genes": 3,
            "primary_family": ["A", "C"],
        },
        "comparator_sets": [],
        "robustness": {"leave_one_gene_out": False},
    }
    result, _, _ = _second_pass_inference(
        data,
        design,
        layout,
        calibration,
        {"A": set(genes[:3]), "C": set(genes[3:6])},
        {"primary": set()},
        cfg,
        chunk_size=36,
        min_per_group=2,
    )
    primary = result[
        result["gene_score"].eq(PRIMARY_SCORE)
        & result["set_test"].eq("competitive")
    ].set_index("gene_set")
    assert (primary["maxT_fwer"] >= primary["empirical_p_greater"]).all()
    assert primary.loc["C", "empirical_p_greater"] == pytest.approx(1 / 36)
    assert primary.loc["C", "maxT_fwer"] == pytest.approx(3 / 36)


def _write_multimodified_workbook(path: Path) -> None:
    from openpyxl import Workbook

    samples = _sample_metadata()
    scaled_headers = [
        f"{row.plex}~rq_{126 + index}_sn scaled"
        for index, row in samples.reset_index(drop=True).iterrows()
    ]
    sample_labels = (
        samples["observed_condition"]
        .str.cat(
            samples.groupby("observed_condition").cumcount().astype(str),
            sep="-",
        )
        .tolist()
    )
    workbook = Workbook()
    single = workbook.active
    single.title = "siteQuant_360"
    single_metadata = [
        "Protein Id",
        "gene_symbol",
        "Site Position",
        "Motif",
        "Max Score",
        "redundancy",
        "sequence",
        "Samp1-5~num_quant",
        "Samp6-10~num_quant",
    ]
    single.append([None] * (len(single_metadata) + len(samples)))
    single.append([None] * len(single_metadata) + sample_labels)
    single.append(single_metadata + scaled_headers)
    single.append(
        [
            "sp|P1|ONE_MOUSE",
            "One",
            10,
            "AAAAAASAAAAAA",
            25,
            "R",
            "R.AAAAAAS#AAAAAA.K",
            4,
            4,
        ]
        + [16.0] * len(samples)
    )
    single.append(
        [
            "sp|P2|MULTI_MOUSE",
            "Multi",
            20,
            "AAAAAATAAAAAA",
            30,
            "R",
            "R.AAA#AAT#AAAAAA.K",
            4,
            4,
        ]
        + [16.0] * len(samples)
    )

    composite = workbook.create_sheet("siteQuant_360_compositeSite")
    composite_metadata = [
        "proteinID",
        "geneSymbol",
        "sitePosStr",
        "motifPeptideStr",
        "maxScoreStr",
        "redundancyStr",
        "sequence",
        "Samp1-5~num_quant",
        "Samp6-10~num_quant",
    ]
    composite.append([None] * (len(composite_metadata) + len(samples)))
    composite.append([None] * len(composite_metadata) + sample_labels)
    composite.append(composite_metadata + scaled_headers)
    pass_row = [
        "sp|P3|COMPOSITE_MOUSE",
        "CompositePass",
        "40;50",
        "AAAAAASAAAAAA;AAAAAATAAAAAA",
        "25;20",
        "R;R",
        "R.AAA#AAT#AAAAAA.K",
        4,
        4,
    ]
    composite.append(pass_row + [16.0] * len(samples))
    composite.append(
        [
            *pass_row[:4],
            "24;21",
            *pass_row[5:],
        ]
        + [16.0] * len(samples)
    )
    composite.append(
        [
            "sp|P4|LOWCOMP_MOUSE",
            "CompositeFail",
            "60;70",
            "AAAAAASAAAAAA;AAAAAATAAAAAA",
            "25;18",
            "R;R",
            "R.AAA#AAT#AAAAAA.K",
            4,
            4,
        ]
        + [16.0] * len(samples)
    )
    workbook.save(path)


def test_deduplicated_multimodified_profile_uses_minimum_composite_localization(
    tmp_path: Path,
):
    assert _minimum_localization_score("25;18") == 18
    assert _minimum_localization_score("24, 21 | 19") == 19
    assert np.isnan(_minimum_localization_score("not-reported"))

    workbook = tmp_path / "multimodified.xlsx"
    _write_multimodified_workbook(workbook)
    profile = AnalysisProfile(
        name="deduplicated_multimodified",
        min_localization_score=19,
        include_composite=True,
        require_singly_modified=False,
        min_finite_total=8,
        min_finite_per_plex=4,
        deduplication="accession_modified_peptide_phosphoform",
        channel_center=False,
    )
    result = prepare_osd462_phosphosites(workbook, profile)
    assert set(result.site_metadata["gene_symbol"]) == {
        "One",
        "Multi",
        "CompositePass",
    }
    composite = result.site_metadata.set_index("gene_symbol").loc["CompositePass"]
    assert composite["localization_score"] == 20
    assert composite["n_source_rows_collapsed"] == 2
    audit = result.filter_audit.set_index("reason")["n_rows"]
    assert audit["below_localization_score"] == 1
    assert audit["not_singly_modified"] == 0


def test_matched_observability_honors_primary_only_scope():
    rows: list[dict[str, object]] = []
    for profile in ("primary", "official_scaled_uncentered"):
        for score in (PRIMARY_SCORE, "mean_negative_site_effect"):
            for index in range(12):
                rows.append(
                    {
                        "profile": profile,
                        "gene_score": score,
                        "gene_symbol": f"G{index}",
                        "fixed_null_universe_eligible": True,
                        "observed_gene_z": 1.0 - index / 10,
                        "n_sites": 1 + index % 4,
                        "median_log2_signal": 8.0 + index / 20,
                        "mean_missing_fraction": (index % 3) / 10,
                    }
                )
    config = {
        "set_test": {
            "minimum_observable_genes": 3,
            "primary_family": ["target"],
        },
        "comparator_sets": ["other"],
        "matched_observability": {
            "profiles": ["primary"],
            "gene_scores": [PRIMARY_SCORE],
            "exclusions": ["primary"],
            "gene_sets": ["target"],
        },
    }
    result = matched_observability_tests(
        pd.DataFrame(rows),
        {
            "target": {"G0", "G1", "G2"},
            "other": {"G3", "G4", "G5"},
        },
        {"primary": set(), "strict": {"G0"}},
        config,
        n_draws=25,
        seed=7,
    )
    assert len(result) == 1
    only = result.iloc[0]
    assert only["profile"] == "primary"
    assert only["gene_score"] == PRIMARY_SCORE
    assert only["exclusion"] == "primary"
    assert only["gene_set"] == "target"


def test_manual_anchor_audit_is_noop_when_none_qualify_and_fails_closed(
    tmp_path: Path,
):
    config = {
        "canonical_exclusions": {
            "strict": ["Slc12a3", "Stk39"],
            "strict_remove_manual_anchor_phosphoforms": True,
        }
    }
    provenance = tmp_path / "anchors.tsv"
    pd.DataFrame(
        [
            {
                "gene_symbol": "Slc12a3",
                "source_sheet": "siteQuant_360",
                "source_excel_row": 10,
                "isolated_canonical_qualified": False,
            }
        ]
    ).to_csv(provenance, sep="\t", index=False)
    summary, audit = audit_manual_anchor_exclusion(provenance, config)
    assert summary["status"] == "implemented_no_qualified_anchor_features"
    assert summary["n_stage0_qualified_anchor_phosphoforms"] == 0
    assert summary["feature_level_removal_numerically_redundant"]
    assert len(audit) == 1

    pd.DataFrame(
        [
            {
                "gene_symbol": "Nedd4l",
                "source_sheet": "siteQuant_360",
                "source_excel_row": 12,
                "isolated_canonical_qualified": True,
                "qualification_reason": "synthetic qualified feature",
            }
        ]
    ).to_csv(provenance, sep="\t", index=False)
    with pytest.raises(
        RuntimeError,
        match="parent genes not removed by the strict exclusion",
    ):
        audit_manual_anchor_exclusion(provenance, config)
