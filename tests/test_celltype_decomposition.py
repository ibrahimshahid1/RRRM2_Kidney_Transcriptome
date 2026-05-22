"""Unit tests for the cell-type marker-panel decomposition (memo Rec 4)."""
import numpy as np
import pandas as pd

from src.multiomics.celltype_panels import (
    KIDNEY_PANELS, decide_scenario, panel_flight_effect, per_sample_panel_scores)


def test_panels_are_disjoint_identity_vs_transport():
    # the identity/transport split must not collapse to the same gene list
    ident = set(KIDNEY_PANELS["dct_identity"])
    transp = set(KIDNEY_PANELS["dct_transport"])
    assert ident != transp
    assert "Slc12a3" in transp and "Slc12a3" not in ident
    assert "Pvalb" in ident


def test_panel_flight_effect_means_mapped_members():
    gs = pd.DataFrame({"gene": ["Slc12a3", "Wnk1", "Stk39", "Pecam1", "Kdr"],
                       "stat": [-2.0, -1.0, -3.0, 4.0, 2.0]})
    pe = panel_flight_effect(gs, KIDNEY_PANELS, min_genes=2).set_index("panel")
    assert pe.loc["dct_transport", "n_mapped"] == 3
    assert np.isclose(pe.loc["dct_transport", "mean_stat"], -2.0)
    assert np.isclose(pe.loc["endothelial", "mean_stat"], 3.0)
    # panel with too few mapped members -> NaN mean
    assert np.isnan(pe.loc["podocyte", "mean_stat"])


def test_per_sample_panel_scores_zscore_mean():
    vst = pd.DataFrame(
        {"s1": [10.0, 9.0, 1.0], "s2": [8.0, 7.0, 5.0], "s3": [6.0, 5.0, 9.0]},
        index=["ENSG_a", "ENSG_b", "ENSG_c"])
    sym_to_ens = {"slc12a3": {"ENSG_a"}, "wnk1": {"ENSG_b"}, "pecam1": {"ENSG_c"}}
    scores = per_sample_panel_scores(vst, {"dct_transport": ["Slc12a3", "Wnk1"]},
                                     sym_to_ens, min_genes=2)
    # row-z then mean across the 2 genes; each gene z has mean 0 across samples
    assert np.isclose(scores.loc["dct_transport"].mean(), 0.0, atol=1e-9)


def test_decide_scenario_transport_suppression():
    d = decide_scenario(-0.5, 0.0, 0.0, 0.0)
    assert d["scenario"] == "transport_program_suppression"


def test_decide_scenario_composition_shift():
    d = decide_scenario(-0.5, -0.5, 0.5, 0.2)
    assert d["scenario"] == "segment_loss_or_composition_shift"
