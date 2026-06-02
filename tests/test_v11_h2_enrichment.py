import numpy as np
import pandas as pd
import pytest
from scipy import stats

from src.v11.core_analysis import (
    FISHER_ALTERNATIVE,
    build_parent_gene_table,
    fisher_table,
    one_single_position_representative_site_per_gene,
    one_representative_site_per_gene,
)
from src.v11.h2_composition_aware_phospho import covariate_diagnostics, fixed_set_adjusted_effects
from src.v11.publication_figures import axis_has_visible_data, fig_h2_primary_enrichment, legend_overlaps_visible_data, plt, qa_axis


RUN_ROOT = "data/results/run_20260526_v11_dct1_phospho_mediation"


def test_v11_fisher_table_uses_directional_greater_alternative():
    df = pd.DataFrame(
        {
            "gene_symbol": list("abcdef"),
            "dct1_enrichment_score": np.arange(6, dtype=float),
            "dct1_top_decile": [True, True, False, False, False, False],
            "suppressed": [True, True, True, False, False, False],
        }
    )
    odds, p, _, arr = fisher_table(df, "dct1_top_decile", "suppressed")
    expected_odds, expected_p = stats.fisher_exact(arr, alternative="greater")
    assert FISHER_ALTERNATIVE == "greater"
    assert odds == expected_odds
    assert p == expected_p


def test_one_representative_site_per_gene_is_not_single_position_filter():
    df = pd.DataFrame(
        {
            "gene_symbol": ["A", "A", "B", "B"],
            "site_id": ["A:1;2", "A:3", "B:4", "B:5"],
            "phospho_p_value": [0.01, 0.20, 0.04, 0.04],
            "phospho_effect": [-0.2, -2.0, 0.5, -0.8],
        }
    )
    selected = one_representative_site_per_gene(df)
    by_gene = dict(zip(selected["gene_symbol"], selected["site_id"]))
    assert by_gene == {"A": "A:1;2", "B": "B:5"}


def test_single_position_representative_site_per_gene_excludes_composites():
    df = pd.DataFrame(
        {
            "gene_symbol": ["A", "A", "B", "B"],
            "site_id": ["A:1;2", "A:3", "B:4;5", "B:6"],
            "phospho_p_value": [0.001, 0.20, 0.01, 0.04],
            "phospho_effect": [-0.2, -2.0, -1.0, -0.8],
            "is_single_site": [False, True, False, True],
        }
    )
    selected = one_single_position_representative_site_per_gene(df)
    by_gene = dict(zip(selected["gene_symbol"], selected["site_id"]))
    assert by_gene == {"A": "A:3", "B": "B:6"}


def test_parent_gene_table_collapses_any_suppressed_site():
    df = pd.DataFrame(
        {
            "gene_symbol": ["A", "A", "B"],
            "site_id": ["A:1", "A:2", "B:1"],
            "is_suppressed_p05": [False, True, False],
            "is_suppressed_q10": [False, False, True],
            "phospho_p_value": [0.2, 0.01, 0.03],
            "phospho_effect": [0.5, -1.0, -0.3],
            "is_single_site": [True, True, True],
            "dct1_enrichment_score": [2.0, 2.0, -1.0],
            "dct1_top_quartile": [True, True, False],
            "dct1_top_decile": [True, True, False],
            "dct2_bottom_quartile": [False, False, True],
            "dct2_bottom_decile": [False, False, True],
            "dct1_core_fdr": [False, False, False],
            "dct2_core_fdr": [False, False, False],
            "flight_effect": [0.1, 0.1, -0.2],
            "n_peptides": [10, 10, 5],
            "abundance_log2": [3.0, 3.0, 2.0],
            "n_fl": [10, 10, 8],
            "n_gc": [10, 10, 7],
        }
    )
    parent = build_parent_gene_table(df)
    row_a = parent[parent["gene_symbol"].eq("A")].iloc[0]
    row_b = parent[parent["gene_symbol"].eq("B")].iloc[0]
    assert bool(row_a["any_suppressed_p05"]) is True
    assert row_a["n_quantified_phosphosites"] == 2
    assert bool(row_b["any_suppressed_q10"]) is True
    assert bool(row_b["any_effect_bottom_quartile"]) is False
    assert row_b["mean_missing_samples"] == 5


def test_primary_enrichment_figure_accepts_current_sensitivity_labels(tmp_path):
    h2 = tmp_path / "h2_enrichment"
    h2.mkdir()
    rows = []
    analyses = [
        "primary_p05",
        "exclude_anchor_genes",
        "exclude_ncc_sites",
        "composite_sites_excluded",
        "one_site_per_parent_gene",
        "single_position_one_site_per_parent_gene",
        "strict_q10",
    ]
    for analysis in analyses:
        for test in ["fisher_dct1_top_decile", "fisher_dct2_bottom_decile", "fisher_dct1_top_quartile"]:
            rows.append(
                {
                    "analysis": analysis,
                    "test": test,
                    "statistic": 1.4 if test.endswith("decile") else 1.1,
                    "table_suppressed_in_flag": 10,
                    "table_suppressed_not_flag": 90,
                    "table_background_in_flag": 50,
                    "table_background_not_flag": 850,
                }
            )
    pd.DataFrame(rows).to_csv(h2 / "h2_dct1_sensitivity_summary.tsv", sep="\t", index=False)
    out = tmp_path / "figures"
    assert fig_h2_primary_enrichment(tmp_path, out)
    assert (out / "v11_dct1_parent_gene_enrichment.png").exists()
    assert (out / "v11_distal_nephron_prior_enrichment.png").exists()


def test_fixed_set_adjusted_effects_tracks_m0_suppressed_membership():
    effects = pd.DataFrame(
        {
            "site_row_id": ["s1", "s2", "s3", "s1", "s2", "s3"],
            "model": ["M0_raw", "M0_raw", "M0_raw", "M4_full", "M4_full", "M4_full"],
            "adjusted_flight_effect": [-0.5, -0.2, 0.3, -0.4, 0.1, -0.3],
            "adjusted_flight_p": [0.01, 0.02, 0.50, 0.03, 0.20, 0.04],
        }
    )
    site_meta = pd.DataFrame(
        {
            "site_row_id": ["s1", "s2", "s3"],
            "dct1_top_decile": [True, False, True],
            "dct1_top_quartile": [True, True, True],
        }
    )
    out = fixed_set_adjusted_effects(effects, site_meta)
    m4_all = out[(out["model"].eq("M4_full")) & (out["fixed_set"].eq("all_m0_suppressed"))].iloc[0]
    assert m4_all["n_sites_from_m0_set"] == 2
    assert m4_all["fraction_adjusted_effect_negative"] == 0.5
    assert "not redefined" in m4_all["boundary"]


def test_covariate_diagnostics_reports_correlations_and_vif_terms():
    long = pd.DataFrame(
        {
            "sample_key": [f"s{i}" for i in range(8)],
            "flight": [0, 0, 0, 0, 1, 1, 1, 1],
            "plex2": [0, 1, 0, 1, 0, 1, 0, 1],
            "dct_identity_score_z": [-1.2, -0.4, 0.2, 1.0, -0.9, -0.1, 0.5, 1.3],
            "endothelial_score_z": [-1.1, -0.7, -0.2, 0.1, 0.2, 0.5, 0.9, 1.4],
            "stromal_score_z": [-1.0, -0.5, -0.3, 0.0, 0.3, 0.7, 1.0, 1.5],
            "composition_pc1_z": [1.1, 0.7, 0.2, -0.1, -0.2, -0.5, -0.9, -1.4],
        }
    )
    out = covariate_diagnostics(long)
    assert "flight__vs__endothelial_score_z" in set(out["term"])
    assert {"flight", "dct_identity_score_z", "endothelial_score_z", "stromal_score_z"} <= set(
        out.loc[out["diagnostic"].eq("vif_M4_sample_level"), "term"]
    )


def test_figure_qa_helpers_detect_data_and_legend_overlap():
    fig, ax = plt.subplots(figsize=(3, 3))
    ax.plot([0.5], [0.5], marker="o", label="point")
    assert axis_has_visible_data(ax)
    ax.legend(loc="center")
    assert legend_overlaps_visible_data(ax)
    with pytest.raises(RuntimeError, match="legend overlaps"):
        qa_axis(ax, "test panel")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3, 3))
    with pytest.raises(RuntimeError, match="no visible plotted data"):
        qa_axis(ax, "empty panel")
    plt.close(fig)


def test_rna_recurrence_supplements_are_family_aware():
    members = pd.read_csv(f"{RUN_ROOT}/baseline/rna_recurrence_gene_set_members.tsv", sep="\t")
    families = pd.read_csv(f"{RUN_ROOT}/baseline/osd462_rna_recurrence_leave_one_family.tsv", sep="\t")
    boot = pd.read_csv(f"{RUN_ROOT}/baseline/osd462_rna_recurrence_paired_pathway_bootstrap.tsv", sep="\t")
    assert {"pathway", "pathway_family", "gene_symbol"} <= set(members.columns)
    assert "matrix_remodeling" in set(members["pathway_family"])
    assert {"dropped_pathway_family", "cosine_after_drop", "delta_vs_full"} <= set(families.columns)
    assert families["dropped_pathway_family"].nunique() >= 3
    assert set(boot["bootstrap_unit"]) == {"paired_curated_pathway"}


def test_tmt_qc_summaries_cover_protein_and_phosphosite_layers():
    qc = pd.read_csv(f"{RUN_ROOT}/baseline/osd462_tmt_channel_qc.tsv", sep="\t")
    missing = pd.read_csv(f"{RUN_ROOT}/baseline/osd462_tmt_missingness_by_condition.tsv", sep="\t")
    assert {"protein", "phosphosite_single", "phosphosite_composite"} <= set(qc["layer"])
    assert {"summed_signal_to_noise", "scaled_signal_to_noise_row_sum_100"} <= set(qc["metric"])
    assert {"BL", "FL", "GC"} <= set(qc["condition"])
    assert {"median_channel_median", "median_missing_fraction"} <= set(missing.columns)


def _apply_filter(df: pd.DataFrame, spec: str) -> pd.DataFrame:
    out = df.copy()
    for token in spec.split(";"):
        col, value = token.split("=", 1)
        out = out[out[col].astype(str).eq(value)]
    return out


def test_v11_headline_numbers_match_generated_artifacts():
    fixture_path = "tests/fixtures/v11_headline_numbers.tsv"
    expected = pd.read_csv(fixture_path, sep="\t")
    for _, row in expected.iterrows():
        artifact = f"{RUN_ROOT}/{row['relative_path']}"
        try:
            df = pd.read_csv(artifact, sep="\t")
        except FileNotFoundError:
            pytest.fail(f"missing headline artifact: {artifact}")
        matched = _apply_filter(df, row["filters"])
        assert len(matched) == 1, f"{row['key']} matched {len(matched)} rows"
        observed = float(matched.iloc[0][row["value_column"]])
        assert observed == pytest.approx(float(row["expected"]), abs=float(row["abs_tolerance"]))
