"""Unit tests for Repair C -- the continuous DCT1<->DCT2 phospho gradient."""
import numpy as np
import pandas as pd

from src.v11.dct_continuous_gradient import (
    fit_linear_gradient,
    fit_logistic_gradient,
    fit_spline_gradient,
    spearman_gradient,
)


def _coord(n=400, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(0.0, 1.0, size=n), rng


def test_linear_gradient_recovers_planted_positive_slope():
    x, rng = _coord(seed=1)
    y = 0.8 * x + rng.normal(0, 0.5, size=x.size)
    res = fit_linear_gradient(x, y)
    assert res["model_status"] == "fit"
    assert res["slope"] > 0.5            # planted +0.8, sign-faithful
    assert res["p_value"] < 1e-6
    assert res["ci_low"] > 0             # CI excludes zero on the supportive side


def test_linear_gradient_flat_data_gives_slope_near_zero():
    x, rng = _coord(seed=2)
    y = rng.normal(0, 1.0, size=x.size)  # no dependence on x
    res = fit_linear_gradient(x, y)
    assert res["model_status"] == "fit"
    assert abs(res["slope"]) < 0.1
    assert res["p_value"] > 0.05


def test_linear_gradient_symmetric_ushape_has_zero_linear_trend():
    # Symmetric U-shape: strong curvature, but zero net *linear* slope.
    x = np.linspace(-2, 2, 400)
    y = x**2  # symmetric about 0 => OLS linear slope ~ 0
    res = fit_linear_gradient(x, y)
    assert res["model_status"] == "fit"
    assert abs(res["slope"]) < 1e-6
    assert res["p_value"] > 0.5


def test_linear_gradient_recovers_planted_negative_slope():
    # The biologically supportive direction for phospho-effect outcomes.
    x, rng = _coord(seed=7)
    y = -0.6 * x + rng.normal(0, 0.5, size=x.size)
    res = fit_linear_gradient(x, y)
    assert res["slope"] < -0.3
    assert res["ci_high"] < 0
    assert res["p_value"] < 1e-6


def test_covariate_adjustment_partials_out_a_confounder():
    x, rng = _coord(seed=3)
    cov = 0.9 * x + rng.normal(0, 0.4, size=x.size)   # correlated nuisance
    y = 1.0 * cov + rng.normal(0, 0.3, size=x.size)    # y driven by cov, not x
    adj = fit_linear_gradient(x, y, covariates=pd.DataFrame({"cov": cov}))
    assert adj["model_status"] == "fit"
    assert adj["covariate_adjusted"] is True
    assert abs(adj["slope"]) < 0.2                     # x effect vanishes after adjustment


def test_spearman_is_sign_faithful_and_robust_to_outliers():
    x, rng = _coord(seed=4)
    y = 0.7 * x + rng.normal(0, 0.5, size=x.size)
    y[:5] = 1e6                                         # gross outliers
    res = spearman_gradient(x, y)
    assert res["model_status"] == "fit"
    assert res["spearman_rho"] > 0.5
    assert res["p_value"] < 1e-6


def test_logistic_gradient_recovers_positive_log_odds():
    x, rng = _coord(n=800, seed=5)
    p = 1.0 / (1.0 + np.exp(-(1.2 * x)))               # P(suppressed) rises with x
    ind = rng.uniform(size=x.size) < p
    res = fit_logistic_gradient(x, ind)
    assert res["model_status"] == "fit"
    assert res["log_odds"] > 0.4
    assert res["odds_ratio"] > 1.0
    assert res["p_value"] < 1e-3


def test_spline_flags_nonmonotone_structure_linear_misses():
    # U-shape: linear slope ~ 0 but spline should detect curvature.
    x = np.linspace(-2, 2, 400)
    y = x**2 + np.random.default_rng(6).normal(0, 0.05, size=x.size)
    lin = fit_linear_gradient(x, y)
    spl = fit_spline_gradient(x, y)
    assert abs(lin["slope"]) < 0.05                    # linear contrast sees nothing
    assert spl["model_status"] == "fit"
    assert spl["p_value"] < 1e-6                        # spline catches the curve
    assert spl["aic_linear_minus_spline"] > 0          # spline is the better fit


def test_spline_no_false_curvature_on_pure_linear_data():
    x, rng = _coord(seed=8)
    y = 0.5 * x + rng.normal(0, 0.5, size=x.size)
    spl = fit_spline_gradient(x, y)
    assert spl["model_status"] == "fit"
    assert spl["p_value"] > 0.05                        # no spurious non-linearity


def test_estimators_refuse_insufficient_data():
    assert fit_linear_gradient([1, 2, 3], [1, 2, 3])["model_status"].startswith("not_fit")
    assert spearman_gradient([1, 2, 3], [1, 2, 3])["model_status"].startswith("not_fit")
