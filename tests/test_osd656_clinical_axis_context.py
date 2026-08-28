from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.clinical_axes.osd656_context import (
    build_osd656_context,
    frozen_axis_gene_catalog,
)
from src.v11.human_concordance import parse_osd656_submitted


REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
OSD656 = (
    REPO
    / "data/external/human_spaceflight/OSD-656"
    / "LSDS-64_Multiplex_urine.immune.AlamarPanel_SUBMITTED.xlsx"
)


def _tiny_config():
    return {
        "primary_family": {
            "axis_a": {
                "subdomains": {
                    "core": {"genes": {"LCN2": 1, "HAVCR1": 1}}
                },
                "secondary_markers": {"SHOULD_NOT_ENTER": 1},
            }
        }
    }


def test_catalog_uses_only_frozen_primary_subdomain_genes():
    catalog = frozen_axis_gene_catalog(_tiny_config())
    assert catalog["gene_symbol"].tolist() == ["LCN2", "HAVCR1"]
    assert "SHOULD_NOT_ENTER" not in set(catalog["gene_symbol"])


def test_pairing_averages_preflight_within_subject_and_keeps_recovery_days():
    rows = []
    for subject, pre_values, recoveries in (
        ("S1", [8.0, 10.0], {1: 12.0, 45: 8.0}),
        ("S2", [4.0, 6.0], {1: 6.0, 45: 7.0}),
    ):
        for day, value in zip((-10, -2), pre_values, strict=True):
            rows.append(
                {
                    "analyte": "LCN2",
                    "normalized_marker": "LCN2",
                    "subject": subject,
                    "timepoint": f"L{day}",
                    "phase": "preflight",
                    "mission_day": day,
                    "concentration_npq": value,
                    "percent_normalized_value": value,
                    "unit": "NPQ",
                    "assay_type": "test",
                }
            )
        for day, value in recoveries.items():
            rows.append(
                {
                    "analyte": "LCN2",
                    "normalized_marker": "LCN2",
                    "subject": subject,
                    "timepoint": f"R+{day}",
                    "phase": "recovery",
                    "mission_day": day,
                    "concentration_npq": value,
                    "percent_normalized_value": value,
                    "unit": "NPQ",
                    "assay_type": "test",
                }
            )
    coverage, paired, summary = build_osd656_context(pd.DataFrame(rows), _tiny_config())
    assert coverage.set_index("gene_symbol").loc["HAVCR1", "panel_status"] == "absent_from_panel"
    s1_r1 = paired[(paired["subject"] == "S1") & (paired["timepoint"] == "R+1")].iloc[0]
    assert s1_r1["subject_preflight_mean_npq"] == pytest.approx(9.0)
    assert s1_r1["delta_recovery_minus_subject_preflight_npq"] == pytest.approx(3.0)
    assert set(summary["timepoint"]) == {"R+1", "R+45"}
    assert not any(column.startswith("p_") for column in summary.columns)


@pytest.mark.skipif(not OSD656.exists(), reason="OSD-656 workbook is not staged")
def test_staged_osd656_overlap_and_exact_paired_descriptives():
    config = yaml.safe_load(CONFIG.read_text())
    long = parse_osd656_submitted(OSD656)
    coverage, paired, summary = build_osd656_context(long, config)
    measured = coverage.loc[
        coverage["panel_status"].eq("measured_overlap"), "gene_symbol"
    ]
    assert set(measured) == {"Lcn2", "Ccl2", "Egf", "Timp1"}
    assert len(coverage) == 34
    assert not coverage.loc[
        coverage["axis"].eq("glomerular_barrier_identity_loss"), "panel_status"
    ].eq("measured_overlap").any()
    assert coverage.loc[
        coverage["gene_symbol"].eq("Havcr1"), "panel_status"
    ].iloc[0] == "absent_from_panel"

    lcn2_r1 = summary[
        summary["gene_symbol"].eq("Lcn2") & summary["timepoint"].eq("R+1")
    ].iloc[0]
    assert lcn2_r1["n_paired_subjects"] == 3
    assert lcn2_r1["mean_delta_npq"] == pytest.approx(1.1238503, abs=1e-6)
    assert lcn2_r1["subject_pattern"] == "all_higher"

    timp1_r82 = summary[
        summary["gene_symbol"].eq("Timp1") & summary["timepoint"].eq("R+82")
    ].iloc[0]
    assert timp1_r82["n_paired_subjects"] == 4
    assert timp1_r82["subject_pattern"] == "all_lower"
    assert paired.groupby(["gene_symbol", "timepoint"])["subject"].nunique().loc[
        "Lcn2", "R+1"
    ] == 3
