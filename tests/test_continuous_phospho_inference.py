"""Tests for the prospective continuous phosphoproteomic inference."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src.v13.continuous_phospho_inference import (
    PRIMARY_SCORE,
    AnalysisProfile,
    PreparedPhosphoData,
    _claim_gate_table,
    _deduplicate,
    _first_pass_calibration,
    _second_pass_inference,
    aggregate_parent_gene_scores,
    build_analysis_profiles,
    build_gene_layout,
    enumerate_balanced_labels,
    load_parent_protein_values,
    load_analysis_config,
    normalization_equivalence_audit,
    prepare_osd462_phosphosites,
    site_contrasts,
)


def test_unique_key_deduplication_fast_path_preserves_sorted_rows_and_values():
    meta = pd.DataFrame(
        {
            "canonical_site_key": ["P2|B", "P1|A", "P3|C"],
            "raw_feature_id": ["sheet:5", "sheet:4", "sheet:6"],
            "localization_score": [20.0, 21.0, 22.0],
            "mean_log2_signal_uncentered": [2.0, 1.0, 3.0],
        }
    )
    values = np.array([[20.0, 21.0], [10.0, 11.0], [30.0, 31.0]])
    observed_meta, observed_values = _deduplicate(meta, values)
    assert observed_meta["canonical_site_key"].tolist() == [
        "P1|A",
        "P2|B",
        "P3|C",
    ]
    assert observed_meta["n_source_rows_collapsed"].eq(1).all()
    assert observed_meta["source_rows"].tolist() == [
        "sheet:4",
        "sheet:5",
        "sheet:6",
    ]
    np.testing.assert_array_equal(
        observed_values,
        np.array([[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]]),
    )


def _sample_metadata(n_per_plex: int = 4) -> pd.DataFrame:
    rows = []
    for plex in ("Samp1-5", "Samp6-10"):
        for i in range(n_per_plex):
            condition = "FL" if i < n_per_plex // 2 else "GC"
            rows.append(
                {
                    "column": f"{plex}|{condition}|{plex}-{i}",
                    "plex": plex,
                    "observed_condition": condition,
                    "sample_id": f"{plex}-{i}",
                }
            )
    return pd.DataFrame(rows)


def test_exact_balanced_assignments_have_observed_first():
    samples = _sample_metadata()
    design = enumerate_balanced_labels(samples, mode="exact")
    # choose(4, 2)^2
    assert design.n_total_possible == 36
    assert design.labels.shape == (36, 8)
    assert np.array_equal(
        design.labels[0], samples["observed_condition"].eq("FL").to_numpy()
    )
    for _, indexes in samples.groupby("plex").groups.items():
        assert (design.labels[:, list(indexes)].sum(axis=1) == 2).all()
    assert np.unique(design.labels, axis=0).shape[0] == 36


def test_site_contrast_equal_weights_valid_plexes_despite_missingness():
    samples = _sample_metadata(n_per_plex=8)
    design = enumerate_balanced_labels(samples, mode="exact")
    # Plex 1 has effect +2 with 4/4 observations.  Plex 2 has effect 0
    # with 3/3 observations.  The frozen estimator is their equal mean (+1),
    # not the missingness-weighted pooled OLS coefficient.
    y = np.array(
        [
            [
                2,
                2,
                2,
                2,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                np.nan,
                0,
                0,
                0,
                np.nan,
            ]
        ],
        dtype=float,
    )
    effect = site_contrasts(y, samples, design.labels[:1], min_per_group=3)[0, 0]
    np.testing.assert_allclose(effect, 1.0, atol=1e-12)

    observed = design.labels[0].astype(float)
    plex2 = samples["plex"].eq("Samp6-10").astype(float).to_numpy()
    keep = np.isfinite(y[0])
    x = np.column_stack([np.ones(keep.sum()), observed[keep], plex2[keep]])
    pooled_ols = np.linalg.lstsq(x, y[0, keep], rcond=None)[0][1]
    assert not np.isclose(effect, pooled_ols)


def test_site_contrast_rejects_two_of_four_in_one_plex_despite_adequate_totals():
    samples = _sample_metadata(n_per_plex=6)
    design = enumerate_balanced_labels(samples, mode="exact")
    y = np.ones((1, 12), dtype=float)
    # Observed assignment has only 2 FL observations in plex 1, while pooled
    # totals remain 5 FL and 6 GC.  Strict within-plex 3/3 validity must fail.
    y[0, 2] = np.nan
    effect = site_contrasts(y, samples, design.labels[:1], min_per_group=3)[0, 0]
    assert np.isnan(effect)


def test_parent_gene_scores_are_signed_and_continuous():
    sites = pd.DataFrame({"gene_symbol": ["A", "A", "A", "B"]})
    layout = build_gene_layout(sites)
    # Site effects are FL-GC.  Negative values mean suppression in flight.
    effects = np.array([[-2.0, -1.0, 4.0, 0.5]])
    scores = aggregate_parent_gene_scores(effects, layout)
    index = {gene: i for i, gene in enumerate(layout.genes)}
    assert scores[PRIMARY_SCORE][0, index["A"]] == 1.0
    assert scores["mean_negative_site_effect"][0, index["A"]] == -1 / 3
    assert scores["one_sided_maxmean"][0, index["A"]] == 1.0
    assert scores[PRIMARY_SCORE][0, index["B"]] == -0.5


def _synthetic_prepared_data() -> PreparedPhosphoData:
    samples = _sample_metadata()
    genes = [f"T{i}" for i in range(4)] + [f"B{i}" for i in range(8)]
    observed_flight = samples["observed_condition"].eq("FL").to_numpy()
    # Every target site has a one-log2-unit flight decrease.  Background genes
    # have balanced idiosyncratic patterns that do not align as a set.
    values = []
    for i, gene in enumerate(genes):
        base = np.array([0.1, -0.1, 0.05, -0.05, 0.0, 0.1, -0.1, 0.0])
        if gene.startswith("T"):
            row = base - observed_flight.astype(float)
        else:
            row = np.roll(base, i)
        values.append(row)
    values = np.asarray(values)
    site_meta = pd.DataFrame(
        {
            "gene_symbol": genes,
            "mean_log2_signal_uncentered": np.linspace(8.0, 9.0, len(genes)),
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
    return PreparedPhosphoData(
        values=values,
        site_metadata=site_meta,
        sample_metadata=samples,
        filter_audit=pd.DataFrame(),
        profile=profile,
    )


def _minimal_config() -> dict:
    return {
        "set_test": {
            "minimum_observable_genes": 3,
            "primary_family": ["target"],
        },
        "comparator_sets": ["background_set"],
        "robustness": {"leave_one_gene_out": True},
    }


def test_full_pipeline_exact_permutation_detects_injected_target():
    data = _synthetic_prepared_data()
    design = enumerate_balanced_labels(data.sample_metadata, mode="exact")
    layout = build_gene_layout(data.site_metadata)
    calibration = _first_pass_calibration(
        data, design, layout, chunk_size=7, min_per_group=2
    )
    set_result, gene_result, loo = _second_pass_inference(
        data,
        design,
        layout,
        calibration,
        {
            "target": {"T0", "T1", "T2", "T3"},
            "background_set": {"B0", "B1", "B2", "B3"},
        },
        {"none": set()},
        _minimal_config(),
        chunk_size=7,
        min_per_group=2,
    )
    primary = set_result[
        set_result["gene_set"].eq("target")
        & set_result["gene_score"].eq(PRIMARY_SCORE)
        & set_result["set_test"].eq("competitive")
    ].iloc[0]
    assert primary["observed_statistic"] > 0
    assert primary["empirical_p_greater"] <= 0.05
    assert primary["maxT_fwer"] <= 0.05
    assert gene_result["observed_gene_z"].notna().all()
    target_loo = loo[
        loo["gene_set"].eq("target") & loo["gene_score"].eq(PRIMARY_SCORE)
    ]
    assert len(target_loo) == 4
    assert (target_loo["observed_competitive_statistic"] > 0).all()


def test_set_membership_excludes_genes_not_estimable_for_every_assignment():
    samples = _sample_metadata(n_per_plex=6)
    observed_flight = samples["observed_condition"].eq("FL").to_numpy()
    genes = ["T0", "T1", "T2", "U", "B0", "B1", "B2", "B3"]
    values = np.zeros((len(genes), len(samples)), dtype=float)
    values[:4, observed_flight] = -1.0
    background_pattern = np.linspace(-0.3, 0.3, len(samples))
    for i in range(4, len(genes)):
        values[i] = np.roll(background_pattern, i)
    # U is estimable under the observed labels (2/2 per plex) but not under
    # every balanced assignment, where its four observed samples can split 3/1.
    for plex_indexes in samples.groupby("plex").groups.values():
        idx = np.asarray(list(plex_indexes))
        values[3, idx[[2, 5]]] = np.nan
    meta = pd.DataFrame(
        {
            "gene_symbol": genes,
            "mean_log2_signal_uncentered": 8.0,
            "missing_fraction": np.isnan(values).mean(axis=1),
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
        site_metadata=meta,
        sample_metadata=samples,
        filter_audit=pd.DataFrame(),
        profile=profile,
    )
    design = enumerate_balanced_labels(samples, mode="exact")
    layout = build_gene_layout(meta)
    calibration = _first_pass_calibration(
        data, design, layout, chunk_size=50, min_per_group=2
    )
    result, genes_out, _ = _second_pass_inference(
        data,
        design,
        layout,
        calibration,
        {
            "target": {"T0", "T1", "T2", "U"},
            "background_set": {"B0", "B1", "B2", "B3"},
        },
        {"none": set()},
        _minimal_config(),
        chunk_size=50,
        min_per_group=2,
    )
    target = result[
        result["gene_set"].eq("target")
        & result["gene_score"].eq(PRIMARY_SCORE)
        & result["set_test"].eq("competitive")
    ].iloc[0]
    assert target["n_observable_genes"] == 3
    unstable = genes_out[
        genes_out["gene_symbol"].eq("U")
        & genes_out["gene_score"].eq(PRIMARY_SCORE)
    ].iloc[0]
    assert not bool(unstable["fixed_null_universe_eligible"])


def _write_synthetic_workbook(path: Path) -> None:
    from openpyxl import Workbook

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "siteQuant_360"
    metadata = [
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
    samples = _sample_metadata()
    scaled_headers = [
        f"{row.plex}~rq_{126 + i}_sn scaled"
        for i, row in samples.reset_index(drop=True).iterrows()
    ]
    sum_headers = [
        f"{row.plex}~rq_{126 + i}_sn_sum"
        for i, row in samples.reset_index(drop=True).iterrows()
    ]
    worksheet.append([None] * (len(metadata) + 2 * len(samples)))
    worksheet.append(
        [None] * len(metadata)
        + samples["observed_condition"].str.cat(
            samples.groupby("observed_condition").cumcount().astype(str),
            sep="-",
        ).tolist()
        * 2
    )
    worksheet.append(metadata + scaled_headers + sum_headers)
    features = [
        ("sp|P1|ONE_MOUSE", "One", 10, "AAAAAASAAAAAA", 20.0, "R", "R.AAAAAS#AAAAAA.K"),
        # Exact duplicate phosphoform: must collapse.
        ("sp|P1|ONE_MOUSE", "One", 10, "AAAAAASAAAAAA", 25.0, "R", "R.AAAAAS#AAAAAA.K"),
        # Multiply modified: must be excluded from the primary universe.
        ("sp|P2|MULTI_MOUSE", "Multi", 20, "AAAAAATAAAAAA", 30.0, "R", "R.AAA#AAT#AAAAAA.K"),
        # Localization score 15: primary excludes; score-13 sensitivity retains.
        ("sp|P3|LOW_MOUSE", "Low", 30, "AAAAAAYAAAAAA", 15.0, "R", "R.AAAAA#YAAAAAA.K"),
    ]
    for feature in features:
        values = [16.0] * len(samples)
        worksheet.append(list(feature) + [4, 4] + values + values)
    workbook.save(path)


def _write_synthetic_protein_workbook(path: Path) -> pd.DataFrame:
    from openpyxl import Workbook

    samples = _sample_metadata()
    samples["sample_id"] = (
        samples["observed_condition"]
        + "-"
        + samples.groupby("observed_condition").cumcount().astype(str)
    )
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "protein_quant_2721"
    metadata = [
        "Protein Id",
        "Gene Symbol",
        "Samp1-5 Peptides",
        "Samp6-10 Peptides",
    ]
    headers = [
        f"{row.plex}~rq_{126 + i}_sn scaled"
        for i, row in samples.reset_index(drop=True).iterrows()
    ]
    worksheet.append([None] * (len(metadata) + len(samples)))
    worksheet.append([None] * len(metadata) + samples["sample_id"].tolist())
    worksheet.append(metadata + headers)
    # Two rows for one symbol.  Plex 1 should favor log2=1, while plex 2
    # should favor log2=3.
    worksheet.append(["sp|P1|G_MOUSE", "G", 9, 1] + [2.0] * len(samples))
    worksheet.append(["sp|P2|G_MOUSE", "G", 1, 9] + [8.0] * len(samples))
    workbook.save(path)
    return samples


def test_workbook_filter_is_phosphoform_unique_and_localization_aware(tmp_path):
    workbook = tmp_path / "phospho.xlsx"
    _write_synthetic_workbook(workbook)
    primary = AnalysisProfile(
        name="primary",
        min_localization_score=19,
        include_composite=False,
        require_singly_modified=True,
        min_finite_total=8,
        min_finite_per_plex=4,
        deduplication="accession_modified_peptide_phosphoform",
        channel_center=False,
    )
    result = prepare_osd462_phosphosites(workbook, primary)
    assert result.site_metadata["gene_symbol"].tolist() == ["One"]
    assert result.site_metadata.iloc[0]["n_source_rows_collapsed"] == 2
    audit = result.filter_audit.set_index("reason")["n_rows"]
    assert audit["below_localization_score"] == 1
    assert audit["not_singly_modified"] == 1
    assert audit["duplicate_rows_collapsed"] == 1

    lower = AnalysisProfile(
        **{
            **primary.__dict__,
            "name": "localization_score_13",
            "min_localization_score": 13,
        }
    )
    sensitivity = prepare_osd462_phosphosites(workbook, lower)
    assert set(sensitivity.site_metadata["gene_symbol"]) == {"One", "Low"}


def test_normalization_equivalence_audit_labels_redundant_inputs(tmp_path):
    workbook = tmp_path / "phospho.xlsx"
    _write_synthetic_workbook(workbook)
    summary = normalization_equivalence_audit(
        workbook, ["siteQuant_360"], tmp_path, min_per_group=2
    )
    assert summary["max_scaled_to_summed_ratio_cv"] == 0
    assert summary["max_abs_uncentered_FL_minus_GC_contrast_difference"] == 0
    assert summary["independent_robustness_evidence"] is False
    assert (tmp_path / "normalization_contrast_equivalence_audit.tsv").exists()


def test_parent_protein_rollup_uses_plex_specific_peptide_weights(tmp_path):
    workbook = tmp_path / "protein.xlsx"
    samples = _write_synthetic_protein_workbook(workbook)
    values = load_parent_protein_values(
        workbook, samples, ["G"], channel_center=False
    )
    plex1 = samples["plex"].eq("Samp1-5").to_numpy()
    plex2 = samples["plex"].eq("Samp6-10").to_numpy()
    np.testing.assert_allclose(values[0, plex1], 1.2)
    np.testing.assert_allclose(values[0, plex2], 2.8)


def test_claim_tier_marks_missing_independent_dct2_set_non_evaluable():
    rows = []
    for profile in (
        "primary",
        "official_scaled_uncentered",
        "signal_to_noise_sum_centered",
        "exclude_multicompartment_broad_expression",
    ):
        for exclusion in ("primary", "strict"):
            for gene_set, statistic, role in (
                ("DCT2_CNT_transition", 1.0, "primary"),
                ("ASDN", 0.9, "primary"),
                ("proximal_tubule", 0.1, "kidney_comparator"),
            ):
                rows.append(
                    {
                        "profile": profile,
                        "exclusion": exclusion,
                        "gene_score": PRIMARY_SCORE,
                        "set_test": "competitive",
                        "gene_set": gene_set,
                        "set_role": role,
                        "n_observable_genes": 10,
                        "observed_statistic": statistic,
                        "maxT_fwer": 0.01 if role == "primary" else np.nan,
                    }
                )
    set_results = pd.DataFrame(rows)
    loo = pd.DataFrame(
        [
            {
                "profile": "primary",
                "exclusion": "primary",
                "gene_score": PRIMARY_SCORE,
                "gene_set": gene_set,
                "removed_gene": f"{gene_set}_{i}",
                "observed_competitive_statistic": 0.5,
            }
            for gene_set in ("DCT2_CNT_transition", "ASDN")
            for i in range(10)
        ]
    )
    cfg = {
        "set_test": {
            "minimum_observable_genes": 8,
            "primary_family": ["DCT2_CNT_transition", "ASDN"],
        },
        "comparator_sets": ["proximal_tubule"],
        "robustness": {
            "require_normalization_direction_concordance": True,
            "broad_expression_exclusion": True,
        },
        "claim_gates": {"alpha": 0.05},
    }
    gates, tier = _claim_gate_table(set_results, loo, cfg)
    dct2 = gates[gates["gene_set"].eq("DCT2_CNT_transition")].iloc[0]
    asdn = gates[gates["gene_set"].eq("ASDN")].iloc[0]
    assert dct2["gate_independent_reference_direction"] == "non_evaluable"
    assert dct2["claim_gate_status"] == "non_evaluable"
    assert asdn["claim_gate_status"] == "pass"
    assert tier.iloc[0]["claim_tier"] == "ASDN_only"
    assert not bool(tier.iloc[0]["dct2_title_permitted"])


def test_production_config_maps_all_frozen_profiles():
    config = Path("config/dct_asdn_phospho_reanalysis.yaml")
    cfg = load_analysis_config(config)
    profiles = {p.name: p for p in build_analysis_profiles(cfg)}
    assert profiles["primary"].min_localization_score == 19
    assert profiles["primary"].signal_source == "official_scaled"
    assert profiles["primary"].channel_center
    assert profiles["localization_score_13"].min_localization_score == 13
    assert profiles["site_position_collapse"].deduplication == "accession_site_position"
    assert profiles["deduplicated_multimodified"].include_composite
    assert not profiles["deduplicated_multimodified"].require_singly_modified
    assert (
        profiles["exclude_non_distal_or_broad_expression"]
        .broad_expression_flag_column
        == "non_distal_specific_or_broad"
    )
    assert (
        profiles["exclude_multicompartment_broad_expression"]
        .broad_expression_flag_column
        == "broadly_expressed"
    )
    assert not profiles["official_scaled_uncentered"].channel_center
    assert profiles["signal_to_noise_sum_centered"].signal_source == "signal_to_noise_sum"
    assert profiles["parent_protein_subtracted"].parent_protein_subtraction
