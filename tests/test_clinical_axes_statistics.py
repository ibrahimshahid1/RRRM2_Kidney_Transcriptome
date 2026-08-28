from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.clinical_axes.statistics import (
    _batch_reml_mkh,
    blocked_meta_permutation,
    combine_fixed_effects,
    genewise_z_scores,
    hedges_g,
    leave_one_gene_out_scores,
    leave_one_mission_out,
    random_effects_reml_mkh,
    score_signed_axis,
)


def test_genewise_z_scores_preserve_missingness_and_neutralize_constants():
    expression = pd.DataFrame(
        [[1.0, 2.0, 3.0], [4.0, 4.0, 4.0], [1.0, np.nan, 3.0]],
        index=["variable", "constant", "missing"],
        columns=["s1", "s2", "s3"],
    )
    z = genewise_z_scores(expression)
    assert z.loc["variable"].mean() == pytest.approx(0.0)
    assert z.loc["variable"].std(ddof=1) == pytest.approx(1.0)
    assert np.allclose(z.loc["constant"], 0.0)
    assert np.isnan(z.loc["missing", "s2"])


def test_signed_axis_applies_directions_and_equal_weights_subdomains():
    expression = pd.DataFrame(
        [
            [0.0, 2.0],
            [0.0, 2.0],
            [0.0, 2.0],
            [2.0, 0.0],
        ],
        index=["g1", "g2", "g3", "g4"],
        columns=["control_like", "engaged_like"],
    )
    directions = {"g1": 1, "g2": 1, "g3": 1, "g4": 1}
    ordinary = score_signed_axis(expression, directions).scores
    equal_domain = score_signed_axis(
        expression,
        directions,
        subdomains={"large": ["g1", "g2", "g3"], "small": ["g4"]},
    )
    assert ordinary["engaged_like"] > 0
    assert equal_domain.scores["engaged_like"] == pytest.approx(0.0, abs=1e-12)
    assert list(equal_domain.subdomain_scores.columns) == ["large", "small"]

    reversed_gene = score_signed_axis(
        expression.loc[["g1", "g4"]], {"g1": 1, "g4": -1}
    ).scores
    assert reversed_gene["engaged_like"] > reversed_gene["control_like"]


def test_median_axis_and_missing_gene_audit():
    expression = pd.DataFrame(
        [[0, 1, 2], [0, 2, 4], [4, 2, 0]],
        index=["a", "b", "c"],
        columns=["s1", "s2", "s3"],
    )
    result = score_signed_axis(
        expression, {"a": 1, "b": 1, "not_measured": -1}, method="median"
    )
    assert result.genes_used == ("a", "b")
    assert result.genes_missing == ("not_measured",)
    assert result.scores["s3"] > result.scores["s1"]
    with pytest.raises(ValueError, match="required genes are missing"):
        score_signed_axis(
            expression, {"a": 1, "not_measured": -1}, require_all_genes=True
        )


def test_hedges_g_matches_declared_small_sample_formula():
    treatment = np.array([2.0, 3.0, 4.0, 5.0])
    control = np.array([0.0, 1.0, 1.0, 2.0])
    result = hedges_g(treatment, control)
    df = 6
    pooled = (
        3 * np.var(treatment, ddof=1) + 3 * np.var(control, ddof=1)
    ) / df
    expected = (1 - 3 / (4 * df - 1)) * (
        np.mean(treatment) - np.mean(control)
    ) / np.sqrt(pooled)
    expected_variance = 8 / 16 + expected**2 / (2 * df)
    assert result.estimate == pytest.approx(expected)
    assert result.variance == pytest.approx(expected_variance)
    assert result.ci_low < result.estimate < result.ci_high


def test_fixed_effect_combines_strata_by_inverse_variance():
    result = combine_fixed_effects([0.2, 0.8], [0.25, 1.0])
    assert result.estimate == pytest.approx(0.32)
    assert result.variance == pytest.approx(0.2)
    assert result.weights.tolist() == pytest.approx([0.8, 0.2])
    assert result.k == 2


def test_reml_mkh_reports_heterogeneity_and_prediction_interval():
    homogeneous = random_effects_reml_mkh(
        [0.45, 0.50, 0.55, 0.48], [0.08, 0.08, 0.08, 0.08]
    )
    assert homogeneous.tau2 == pytest.approx(0.0, abs=1e-9)
    conventional_se = np.sqrt(1 / np.sum(1 / np.repeat(0.08, 4)))
    assert homogeneous.standard_error >= conventional_se
    assert homogeneous.ci_low < homogeneous.estimate < homogeneous.ci_high
    assert homogeneous.prediction_df == 2
    assert homogeneous.prediction_low < homogeneous.prediction_high

    heterogeneous = random_effects_reml_mkh(
        [-1.0, 0.2, 1.5, 2.0], [0.05, 0.05, 0.05, 0.05]
    )
    assert heterogeneous.tau2 > 0
    assert heterogeneous.i_squared > 50
    assert heterogeneous.prediction_low < heterogeneous.ci_low
    assert heterogeneous.prediction_high > heterogeneous.ci_high


def test_vectorized_reml_matches_scalar_engine():
    effects = np.array([[0.2, 0.5, 0.7, 0.4], [-0.8, 0.1, 1.1, 1.7]])
    variances = np.array([[0.1, 0.2, 0.15, 0.1], [0.08, 0.12, 0.09, 0.11]])
    estimate, tau2, se, t_value = _batch_reml_mkh(effects, variances)
    for row in range(len(effects)):
        scalar = random_effects_reml_mkh(effects[row], variances[row])
        assert estimate[row] == pytest.approx(scalar.estimate, abs=2e-6)
        assert tau2[row] == pytest.approx(scalar.tau2, abs=2e-6)
        assert se[row] == pytest.approx(scalar.standard_error, abs=2e-6)
        assert t_value[row] == pytest.approx(scalar.t, abs=2e-6)


def _three_mission_design() -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(91)
    rows = []
    scores = []
    for mission in ("M1", "M2", "M3"):
        for stratum in ("young", "old"):
            for group in ("control", "flight"):
                for animal in range(4):
                    sample = f"{mission}_{stratum}_{group}_{animal}"
                    rows.append((sample, mission, stratum, group))
                    strong = rng.normal(0, 0.25) + (2.5 if group == "flight" else 0.0)
                    null = rng.normal(0, 1.0)
                    scores.append((sample, strong, null))
    design = pd.DataFrame(
        rows, columns=["sample", "mission", "stratum", "group"]
    ).set_index("sample")
    score_frame = pd.DataFrame(
        scores, columns=["sample", "injury_axis", "null_axis"]
    ).set_index("sample")
    return score_frame, design


def test_blocked_meta_permutation_controls_one_frozen_axis_family():
    scores, design = _three_mission_design()
    result = blocked_meta_permutation(
        scores,
        design,
        n_permutations=399,
        seed=122,
        chunk_size=57,
    )
    assert result.null_t.shape == (399, 2)
    assert result.null_max_abs_t.shape == (399,)
    assert set(result.observed_meta.index) == {"injury_axis", "null_axis"}
    assert len(result.mission_effects) == 6
    assert np.all(result.max_t_fwer >= result.empirical_p_two_sided)
    assert result.max_t_fwer["injury_axis"] <= 0.05
    assert result.observed_meta.loc["injury_axis", "estimate"] > 0

    repeated = blocked_meta_permutation(
        scores,
        design,
        n_permutations=399,
        seed=122,
        chunk_size=57,
    )
    assert np.array_equal(result.null_t.to_numpy(), repeated.null_t.to_numpy())


def test_leave_one_out_helpers_preserve_labels_and_sample_order():
    loo_meta = leave_one_mission_out(
        [0.2, 0.4, 0.5, 0.7],
        [0.1, 0.1, 0.1, 0.1],
        ["M1", "M2", "M3", "M4"],
    )
    assert loo_meta["omitted_mission"].tolist() == ["M1", "M2", "M3", "M4"]
    assert np.all(loo_meta["k_missions"] == 3)

    expression = pd.DataFrame(
        [[0, 1, 2], [0, 2, 4], [2, 1, 0]],
        index=["a", "b", "c"],
        columns=["s1", "s2", "s3"],
    )
    loo_genes = leave_one_gene_out_scores(expression, {"a": 1, "b": 1, "c": -1})
    assert list(loo_genes.index) == ["s1", "s2", "s3"]
    assert list(loo_genes.columns) == ["a", "b", "c"]
    assert loo_genes.columns.name == "omitted_gene"
