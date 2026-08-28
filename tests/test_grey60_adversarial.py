from __future__ import annotations

import numpy as np
import pandas as pd

from src.grey60.adversarial import (
    familywise_p,
    hedges_g,
    max_t_permutation,
    mean_z_score,
    osd462_animal_id,
    random_effects_reml_hk,
    stratified_bootstrap_iss_effect,
)


def _balanced_metadata() -> pd.DataFrame:
    rows = []
    for arm in ("ISS-T", "LAR"):
        for age in ("YNG", "OLD"):
            for env in ("GC", "FLT"):
                for i in range(5):
                    rows.append((f"{arm}_{age}_{env}_{i}", arm, age, env))
    return pd.DataFrame(rows, columns=["sample", "Arm", "Age", "EnvGroup"]).set_index("sample")


def test_mean_z_score_has_expected_sample_order():
    expr = pd.DataFrame(
        [[1, 2, 3, 4], [2, 4, 6, 8]],
        index=["a", "b"],
        columns=["s1", "s2", "s3", "s4"],
    )
    score = mean_z_score(expr, ["a", "b"])
    assert list(score.index) == ["s1", "s2", "s3", "s4"]
    assert score.iloc[-1] > score.iloc[0]


def test_stratified_bootstrap_detects_positive_iss_effect():
    meta = _balanced_metadata()
    score = pd.Series(0.0, index=meta.index)
    score.loc[(meta["Arm"] == "ISS-T") & (meta["EnvGroup"] == "FLT")] = 2.0
    boot = stratified_bootstrap_iss_effect(score, meta, n_boot=100, seed=4)
    assert np.allclose(boot, 2.0)


def test_max_t_null_has_declared_shape_and_valid_fwer():
    meta = _balanced_metadata()
    rng = np.random.default_rng(3)
    responses = pd.DataFrame(
        rng.normal(size=(40, 20)),
        index=meta.index,
        columns=[f"M{i}" for i in range(20)],
    )
    null = max_t_permutation(responses, meta, n_perm=100, seed=7, chunk_size=25)
    assert null.shape == (100,)
    assert np.all(np.isfinite(null))
    p = familywise_p(null, 0.0)
    assert p == 1.0


def test_hedges_g_and_random_effects_positive():
    g = hedges_g([2, 3, 4, 5], [0, 1, 1, 2])
    assert g.estimate > 0
    assert g.variance > 0
    meta = random_effects_reml_hk([0.4, 0.6, 0.5], [0.1, 0.1, 0.1])
    assert 0.4 < meta.estimate < 0.6
    assert meta.ci_low_hk > 0
    assert np.isclose(meta.weights.sum(), 1.0)


def test_modified_hartung_knapp_does_not_shrink_uncertainty():
    standard = random_effects_reml_hk([0.4, 0.6, 0.5], [0.1, 0.1, 0.1])
    modified = random_effects_reml_hk(
        [0.4, 0.6, 0.5], [0.1, 0.1, 0.1], modified=True
    )
    assert modified.standard_error_hk >= standard.standard_error_hk
    assert modified.ci_low_hk <= standard.ci_low_hk
    assert modified.ci_high_hk >= standard.ci_high_hk


def test_osd462_animal_id_strips_preparation_and_technical_replicate():
    assert (
        osd462_animal_id("RR10_KDN_WT_FLT_F1_UPX_techrep1")
        == "RR10_KDN_WT_FLT_F1"
    )
    assert osd462_animal_id("RR10_KDN_WT_GC_G3_totRNA") == "RR10_KDN_WT_GC_G3"
