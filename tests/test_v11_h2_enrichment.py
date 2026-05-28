import numpy as np
import pandas as pd
from scipy import stats

from src.v11.core_analysis import (
    FISHER_ALTERNATIVE,
    build_parent_gene_table,
    fisher_table,
    one_representative_site_per_gene,
)
from src.v11.publication_figures import fig_h2_primary_enrichment


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
        }
    )
    parent = build_parent_gene_table(df)
    row_a = parent[parent["gene_symbol"].eq("A")].iloc[0]
    row_b = parent[parent["gene_symbol"].eq("B")].iloc[0]
    assert bool(row_a["any_suppressed_p05"]) is True
    assert row_a["n_quantified_phosphosites"] == 2
    assert bool(row_b["any_suppressed_q10"]) is True


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
        "strict_q10",
    ]
    for analysis in analyses:
        for test in ["fisher_dct1_top_decile", "fisher_dct1_top_quartile"]:
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
