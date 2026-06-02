"""Unit tests for Reach D human urine/fluid-axis concordance."""

import numpy as np
import pandas as pd

from src.v11.human_concordance import (
    build_axis_concordance,
    concordance_verdict,
    figure_evidence_table,
    observed_direction,
    osd656_relevance_category,
    parse_osd656_timepoint,
    parse_timepoint,
    summarize_osd656_prepost,
    summarize_table_analytes,
)


def test_parse_timepoint_classifies_twins_labels_and_qc():
    assert parse_timepoint("L-162")["phase"] == "preflight"
    assert parse_timepoint("FD237")["phase"] == "inflight"
    assert parse_timepoint("R+34/35 (pooled)")["phase"] == "recovery"
    flagged = parse_timepoint("L-71 (tube frozen before aliquotting)")
    assert flagged["phase"] == "preflight"
    assert flagged["qc_excluded"] is True


def test_parse_osd656_timepoint_classifies_i4_recovery_axis():
    assert parse_osd656_timepoint("L-44", "Preflight")["phase"] == "preflight"
    recovered = parse_osd656_timepoint("R+45", "R+45")
    assert recovered["phase"] == "recovery"
    assert recovered["day"] == 45


def test_observed_direction_handles_sign_and_flat():
    assert observed_direction(1.2) == "up"
    assert observed_direction(-0.4) == "down"
    assert observed_direction(0.0) == "flat"
    assert observed_direction(np.nan) == "flat"


def _long_fixture():
    rows = []
    specs = {
        "Urinary sodium mmol/day": [100, 110, 140, 150],
        "24 h volume mL": [2000, 1800, 1500, 1400],
        "Urinary magnesium mg/day": [20, 22, 55, 60],
        "Urinary potassium mmol/day": [70, 70, 72, 73],
    }
    phases = ["preflight", "preflight", "inflight", "inflight"]
    for label, vals in specs.items():
        for i, (phase, val) in enumerate(zip(phases, vals)):
            rows.append({
                "analyte_label": label,
                "phase": phase,
                "value": val,
                "qc_excluded": False,
                "timepoint": f"t{i}",
            })
    return pd.DataFrame(rows)


def test_summarize_table_analytes_is_sign_faithful():
    summary = summarize_table_analytes(_long_fixture()).set_index("analyte")
    assert summary.loc["urinary_sodium", "observed_direction"] == "up"
    assert summary.loc["urinary_sodium", "concordant"] is True
    assert summary.loc["urine_volume_24h", "observed_direction"] == "down"
    assert summary.loc["urine_volume_24h", "concordant"] is True
    assert summary.loc["urinary_magnesium", "observed_direction"] == "up"
    assert summary.loc["urinary_magnesium", "concordant"] is True
    assert summary.loc["urinary_potassium", "scored"] == False


def test_figure_evidence_scores_aqp2_and_agt_but_not_renr():
    fig = figure_evidence_table().set_index("analyte")
    assert fig.loc["AQP2", "concordant"] is True
    assert fig.loc["AGT", "concordant"] is True
    assert fig.loc["RENR", "scored"] == False


def test_osd656_prepost_summary_is_recovery_context_not_primary_validation():
    long = pd.DataFrame({
        "analyte": ["CCL2", "CCL2", "AQP2", "AQP2", "PanelOnly", "PanelOnly"],
        "normalized_marker": ["CCL2", "CCL2", "AQP2", "AQP2", "PANELONLY", "PANELONLY"],
        "subject": ["C001", "C001", "C001", "C001", "C001", "C001"],
        "phase": ["preflight", "recovery", "preflight", "recovery", "preflight", "recovery"],
        "mission_day": [-44, 1, -44, 1, -44, 1],
        "unit": ["NPQ", "NPQ", "NPQ", "NPQ", "NPQ", "NPQ"],
        "assay_type": ["Alamar panel"] * 6,
        "concentration_npq": [10.0, 12.0, 2.0, 3.0, 5.0, 4.0],
        "percent_normalized_value": [100.0, 120.0, 100.0, np.inf, 100.0, 80.0],
        "relevance_category": [
            osd656_relevance_category("CCL2"),
            osd656_relevance_category("CCL2"),
            osd656_relevance_category("AQP2"),
            osd656_relevance_category("AQP2"),
            osd656_relevance_category("PanelOnly"),
            osd656_relevance_category("PanelOnly"),
        ],
    })
    summary = summarize_osd656_prepost(long).set_index("analyte")
    assert summary.loc["CCL2", "relevance_category"] == "recovery_inflammation_matrix_context"
    assert summary.loc["CCL2", "context_use"] == "optional_recovery_inflammation_context"
    assert summary.loc["AQP2", "relevance_category"] == "direct_distal_raas_marker"
    assert summary.loc["AQP2", "context_use"] == "potential_direct_marker_if_present"
    assert np.isfinite(summary.loc["AQP2", "recovery_mean_percent_normalized"]) == False
    assert summary.loc["CCL2", "observed_direction_npq"] == "up"
    assert "not independent inflight AQP2 validation" in summary.loc["AQP2", "caveat"]


def test_verdict_counts_only_scored_rows():
    table_summary = summarize_table_analytes(_long_fixture())
    axis = build_axis_concordance(table_summary, figure_evidence_table())
    verdict = concordance_verdict(axis)
    assert verdict["n_scored_analytes"] == 5
    assert verdict["n_concordant"] == 5
    assert verdict["status"] == "directionally_concordant_all_scored"
