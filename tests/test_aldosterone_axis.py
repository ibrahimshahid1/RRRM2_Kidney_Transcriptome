"""Unit tests for Reach F -- the aldosterone / mineralocorticoid-axis test."""
import numpy as np
import pandas as pd

from src.v11.aldosterone_axis import (
    competitive_permutation,
    meta_axis,
    panel_signs,
    score_axis_in_cohort,
)


def _panel():
    return pd.DataFrame({
        "gene_symbol": ["Sgk1", "Scnn1a", "Slc12a3", "Hsd11b2", "Nr3c2"],
        "expected_aldosterone_direction": ["up", "up", "up", "context", "context"],
        "role": ["k"] * 5,
        "evidence": ["e"] * 5,
    })


def test_panel_signs_maps_up_to_plus_one_and_context_to_zero():
    s = panel_signs(_panel())
    assert s["Sgk1"] == 1.0 and s["Scnn1a"] == 1.0 and s["Slc12a3"] == 1.0
    assert s["Hsd11b2"] == 0.0 and s["Nr3c2"] == 0.0


def test_score_axis_positive_when_up_genes_rise():
    signs = panel_signs(_panel())
    gs = pd.Series({"Sgk1": 2.0, "Scnn1a": 2.0, "Slc12a3": 2.0, "Hsd11b2": -5.0, "Nr3c2": 9.0})
    res = score_axis_in_cohort(gs, signs)
    assert res["n_present"] == 3                 # context genes excluded
    assert np.isclose(res["effect"], 2.0)        # mean of direction-corrected up-genes


def test_score_axis_negative_under_suppression():
    signs = panel_signs(_panel())
    gs = pd.Series({"Sgk1": -1.5, "Scnn1a": -2.5, "Slc12a3": -2.0})
    res = score_axis_in_cohort(gs, signs)
    assert res["effect"] < 0
    assert np.isclose(res["effect"], -2.0)


def test_score_axis_ignores_missing_genes():
    signs = panel_signs(_panel())
    gs = pd.Series({"Sgk1": -2.0})               # only one panel gene present
    res = score_axis_in_cohort(gs, signs)
    assert res["n_present"] == 1
    assert np.isclose(res["effect"], -2.0)


def test_competitive_permutation_flags_coherent_suppression():
    rng = np.random.default_rng(0)
    universe = pd.Series(rng.normal(0, 1, size=2000),
                         index=[f"g{i}" for i in range(2000)])
    panel = ["g0", "g1", "g2", "g3", "g4"]
    universe.loc[panel] = [-3.0, -3.2, -2.8, -3.1, -2.9]   # strongly, coherently negative
    signs = np.ones(5)
    res = competitive_permutation(universe, panel, signs, n_perm=2000, seed=1)
    assert res["p_suppression"] < 0.01                      # observed far in the suppression tail
    assert res["p_two_sided"] < 0.05


def test_competitive_permutation_not_significant_for_random_panel():
    rng = np.random.default_rng(2)
    universe = pd.Series(rng.normal(0, 1, size=2000),
                         index=[f"g{i}" for i in range(2000)])
    panel = list(universe.index[:5])                        # ordinary genes, no shift
    res = competitive_permutation(universe, panel, np.ones(5), n_perm=2000, seed=3)
    assert res["p_two_sided"] > 0.05


def test_meta_axis_all_negative_is_recurrent_down_and_sign_consistent():
    res = meta_axis({"c1": -1.2, "c2": -0.8, "c3": -1.5, "c4": -0.9})
    assert res["meta_effect"] < 0
    assert res["meta_z"] < 0
    assert res["sign_consistent"] is True
    assert res["n_negative"] == 4
    assert res["recurrence_class"] == "recurrent_down"


def test_meta_axis_mixed_signs_not_consistent():
    res = meta_axis({"c1": 1.2, "c2": -1.1, "c3": 0.9, "c4": -1.3})
    assert res["sign_consistent"] is False
    assert res["recurrence_class"] in {"mixed_direction", "not_recurrent"}


def test_meta_axis_all_positive_is_recurrent_up():
    res = meta_axis({"c1": 1.2, "c2": 0.8, "c3": 1.5})
    assert res["meta_z"] > 0
    assert res["recurrence_class"] == "recurrent_up"
    assert res["n_positive"] == 3
