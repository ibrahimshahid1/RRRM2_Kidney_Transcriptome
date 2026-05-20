import numpy as np
import pandas as pd

from src.networks.contrast_vectors import (
    categorize_bootstrap_distribution,
    classify_beta,
    compute_decomposition,
    redirected_component,
    summarize_bootstrap_decomposition,
)


def test_known_geometry_amplification_and_redirected_component():
    a_gc = np.array([1.0, 0.0])
    a_flt = np.array([2.0, 0.0])
    dec = compute_decomposition(a_flt, a_gc)
    assert dec.beta == 2.0
    assert dec.cos == 1.0
    assert dec.rho == 0.0
    assert np.allclose(redirected_component(a_flt, a_gc), [0.0, 0.0])


def test_known_geometry_redirect_category():
    a_gc = np.array([1.0, 0.0])
    a_flt = np.array([0.0, 2.0])
    dec = compute_decomposition(a_flt, a_gc)
    assert dec.beta == 0.0
    assert dec.cos == 0.0
    assert dec.rho == 1.0
    assert classify_beta(dec.beta, dec.rho) == "redirect"


def test_bootstrap_summary_reports_category_fractions():
    point = compute_decomposition(np.array([2.0, 0.0]), np.array([1.0, 0.0]))
    boot = pd.DataFrame({
        "beta": [1.2, 1.4, 1.6, 0.5],
        "cos": [0.9, 0.8, 0.7, 0.6],
        "rho": [0.1, 0.1, 0.1, 0.2],
        "a_flt_norm": [2.0, 2.0, 2.0, 2.0],
        "a_gc_norm": [1.0, 1.0, 1.0, 1.0],
        "residual_norm": [0.0, 0.0, 0.0, 0.0],
    })
    summary = summarize_bootstrap_decomposition(point, boot)
    beta_row = summary[summary["statistic"] == "beta"].iloc[0]
    assert beta_row["frac_amplify"] == 0.75
    assert beta_row["interpretation"] == "amplify"

    cats = categorize_bootstrap_distribution(np.array([1.1, 0.5]), np.array([0.2, 0.2]))
    assert cats["headline"] is None
