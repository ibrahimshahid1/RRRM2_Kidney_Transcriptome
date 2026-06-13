"""Unit tests for the v10 regulator-activity layer."""
import numpy as np
import pandas as pd

from src.multiomics.regulator_activity import (
    grade_axis,
    ksea,
    ksea_positive_control_passes,
    load_kinase_substrate_net,
    recurrence_class,
)


def _net():
    return pd.DataFrame({
        "kinase": ["K_down"] * 4 + ["K_flat"] * 4,
        "substrate_gene": ["G"] * 8,
        "substrate_site": ["1", "2", "3", "4", "10", "11", "12", "13"],
    })


def _sites():
    # background: 200 sites ~N(0,1); K_down substrates strongly negative;
    # K_flat substrates ~0
    rng = np.random.default_rng(0)
    bg = rng.normal(0, 1, 200)
    rows = [{"gene_symbol": "BG", "site_position": str(i), "phospho_effect": bg[i]}
            for i in range(200)]
    for site in ["1", "2", "3", "4"]:
        rows.append({"gene_symbol": "G", "site_position": site, "phospho_effect": -2.5})
    for site in ["10", "11", "12", "13"]:
        rows.append({"gene_symbol": "G", "site_position": site, "phospho_effect": 0.02})
    return pd.DataFrame(rows)


def test_ksea_detects_suppressed_kinase():
    out = ksea(_sites(), _net(), min_substrates=3)
    kd = out[out["kinase"].eq("K_down")].iloc[0]
    assert kd["z_score"] < 0
    assert kd["p_value"] < 0.01
    assert kd["direction"] == "inferred_activity_down"
    assert kd["n_substrates_quantified"] == 4


def test_ksea_flat_kinase_not_significant():
    out = ksea(_sites(), _net(), min_substrates=3)
    kf = out[out["kinase"].eq("K_flat")].iloc[0]
    assert abs(kf["z_score"]) < 2.0
    assert kf["p_value"] > 0.05


def test_ksea_insufficient_substrates_flagged():
    net = pd.DataFrame({
        "kinase": ["K_small", "K_small"],
        "substrate_gene": ["G", "G"],
        "substrate_site": ["1", "2"],
    })
    out = ksea(_sites(), net, min_substrates=3)
    assert out.iloc[0]["status"] == "insufficient_substrates"
    assert not np.isfinite(out.iloc[0]["z_score"])


def test_ksea_collapses_duplicate_sites():
    sites = _sites()
    # add a duplicate measurement of an existing substrate site
    dup = pd.DataFrame([{"gene_symbol": "G", "site_position": "1", "phospho_effect": -2.5}])
    out = ksea(pd.concat([sites, dup], ignore_index=True), _net(), min_substrates=3)
    kd = out[out["kinase"].eq("K_down")].iloc[0]
    assert kd["n_substrates_quantified"] == 4  # not 5


def test_positive_control_helper():
    out = ksea(_sites(), _net(), min_substrates=3)
    # rename K_down -> SPAK_OSR1 to mimic the control kinase
    out2 = out.copy()
    out2["kinase"] = out2["kinase"].replace({"K_down": "SPAK_OSR1"})
    assert ksea_positive_control_passes(out2, control_kinases=["SPAK_OSR1"])
    assert not ksea_positive_control_passes(
        out.assign(kinase=out["kinase"].replace({"K_flat": "SPAK_OSR1"})),
        control_kinases=["SPAK_OSR1"],
    )


def test_recurrence_class():
    assert recurrence_class({"a": -1.5, "b": -2.0, "c": -1.2}) == "recurrent_down"
    assert recurrence_class({"a": 1.5, "b": 2.0}) == "recurrent_up"
    assert recurrence_class({"a": 1.5, "b": -2.0}) == "mixed_direction"
    assert recurrence_class({"a": 0.1, "b": 0.2}) == "not_recurrent"
    assert recurrence_class({"a": 1.5}) == "insufficient_cohorts"


def test_grade_axis():
    assert grade_axis(rna_recurrence="recurrent_down", phospho_support=True,
                      protein_boundary_null=True, is_positive_control=True) == "activity_anchor"
    assert grade_axis(rna_recurrence="recurrent_up", phospho_support=True,
                      protein_boundary_null=True) == "candidate_upstream_organizer"
    assert grade_axis(rna_recurrence="recurrent_up", phospho_support=False,
                      protein_boundary_null=True) == "candidate_context_axis"
    assert grade_axis(rna_recurrence="not_recurrent", phospho_support=False,
                      protein_boundary_null=True) == "negative_boundary"


def test_load_kinase_substrate_net_core():
    net = load_kinase_substrate_net(
        "data/external/kinase_substrate/renal_kinase_substrate_core.tsv")
    assert {"kinase", "substrate_gene", "substrate_site"}.issubset(net.columns)
    assert set(net["kinase"]) == {"SPAK_OSR1", "WNK"}
