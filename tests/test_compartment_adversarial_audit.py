from __future__ import annotations

import numpy as np
import pandas as pd

from src.v13.compartment_adversarial_audit import (
    bh_fdr,
    build_marker_tiers,
    competitive_statistic,
)


def _tiny_config():
    return {
        "marker_tiers": {
            "target_mean_cpm_min": 1.0,
            "target_source_study_detection_fraction_min": 0.5,
            "all_enriched_log2_target_to_max_other_min": 1.0,
            "within_kidney_broad_non_target_count_min": 2,
            "within_kidney_broad_non_target_mean_cpm_min": 1.0,
            "within_kidney_broad_non_target_detection_fraction_min": 0.5,
            "high_specificity_log2_target_to_max_other_min": 2.0,
            "high_specificity_target_detection_fraction_min": 0.75,
        },
        "compartment_mapping": {"A": "a", "B": "b", "C": "c"},
        "structural_scaffold_control": {
            "gene_set_name": "scaffold__all",
            "source_terms": ["Actin"],
        },
        "set_test": {
            "primary_family": [
                f"{compartment}__{tier}"
                for compartment in ("a", "b", "c")
                for tier in (
                    "all_enriched",
                    "within_kidney_not_broad",
                    "high_specificity",
                    "broad_enriched",
                    "scaffold_excluded",
                )
            ]
            + ["scaffold__all"]
        },
    }


def test_marker_tiers_partition_breadth_and_remove_scaffold():
    values = {
        "restricted": (8.0, 0.5, 0.2),
        "broad": (8.0, 2.0, 2.0),
        "scaffold": (8.0, 0.5, 0.2),
    }
    rows = []
    for gene, triple in values.items():
        for compartment, cpm in zip(("A", "B", "C"), triple):
            rows.append(
                {
                    "gene_symbol": gene,
                    "compartment": compartment,
                    "mean_cpm": cpm,
                    "source_study_detection_fraction": 1.0,
                }
            )
    # Ensure each compartment emits every frozen tier in the tiny fixture.
    for compartment, prefix in (("B", "b_only"), ("C", "c_only")):
        for gene, cpm_a, cpm_b, cpm_c in (
            (prefix + "_restricted", 0.2, 8.0, 0.2),
            (prefix + "_broad", 2.0, 8.0, 2.0),
        ):
            if compartment == "C":
                if gene.endswith("_restricted"):
                    cpm_a, cpm_b, cpm_c = 0.2, 0.2, 8.0
                else:
                    cpm_a, cpm_b, cpm_c = 2.0, 2.0, 8.0
            for label, cpm in zip(("A", "B", "C"), (cpm_a, cpm_b, cpm_c)):
                rows.append(
                    {
                        "gene_symbol": gene,
                        "compartment": label,
                        "mean_cpm": cpm,
                        "source_study_detection_fraction": 1.0,
                    }
                )
    table = build_marker_tiers(
        pd.DataFrame(rows), {"Actin": ["SCAFFOLD"]}, _tiny_config()
    )
    a_all = set(table.loc[table["gene_set"].eq("a__all_enriched"), "gene_symbol"])
    a_narrow = set(
        table.loc[table["gene_set"].eq("a__within_kidney_not_broad"), "gene_symbol"]
    )
    a_broad = set(table.loc[table["gene_set"].eq("a__broad_enriched"), "gene_symbol"])
    a_no_scaffold = set(
        table.loc[table["gene_set"].eq("a__scaffold_excluded"), "gene_symbol"]
    )
    assert a_narrow | a_broad == a_all
    assert not (a_narrow & a_broad)
    assert "scaffold" in a_all
    assert "scaffold" not in a_no_scaffold


def test_bh_and_competitive_statistic_are_deterministic():
    q = bh_fdr([0.01, 0.04, np.nan, 0.03])
    assert np.allclose(q[[0, 1, 3]], [0.03, 0.04, 0.04])
    assert np.isnan(q[2])
    values = np.array([3.0, 1.0, -1.0, -3.0])
    membership = np.array([True, True, False, False])
    assert competitive_statistic(values, membership) == 4.0
