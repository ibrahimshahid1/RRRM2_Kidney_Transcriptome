"""Contract tests for K0 (podocyte-vs-structural-scaffold specificity)."""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml

from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions
from src.clinical_axes.statistics import random_effects_reml_mkh, score_signed_axis
from scripts.clinical_axes.run_podocyte_scaffold_specificity import (
    PODOCYTE_SET,
    STRUCTURAL_SET,
    _design,
    _frozen_genes,
    _MissionPerm,
    _regression_rows,
)

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"


def test_fwl_closed_form_flight_coefficient_matches_statsmodels():
    """The vectorized FWL coefficient/variance used in the permutation must equal
    the model-based statsmodels OLS estimate on the full design."""
    rng = np.random.default_rng(0)
    n = 20
    idx = [f"a{i}" for i in range(n)]
    metadata = pd.DataFrame(
        {
            "condition": ["FLT"] * 10 + ["GC"] * 10,
            "block": (["young"] * 5 + ["old"] * 5) * 2,
        },
        index=idx,
    )
    adjuster = pd.Series(rng.normal(size=n), index=idx)
    outcome = pd.Series(rng.normal(size=n), index=idx)

    design = _design(metadata, adjuster)
    fit = sm.OLS(outcome, design).fit()
    beta_sm = float(fit.params["flight"])
    se_sm = float(fit.bse["flight"])

    mp = _MissionPerm(metadata, outcome, adjuster)
    flight = (metadata["condition"] == "FLT").to_numpy(float)[None, :]
    beta, var = mp.coeff_var(flight)
    assert beta[0] == np.float64(beta_sm) or np.isclose(beta[0], beta_sm, atol=1e-9)
    assert np.isclose(np.sqrt(var[0]), se_sm, atol=1e-9)


def test_permuted_flights_preserve_block_counts():
    """Blocked permutation must keep each block's flight/control counts fixed."""
    idx = [f"a{i}" for i in range(19)]
    metadata = pd.DataFrame(
        {
            "condition": ["FLT"] * 10 + ["GC"] * 9,
            "block": (["d25"] * 5 + ["d75"] * 5) + (["d25"] * 5 + ["d75"] * 4),
        },
        index=idx,
    )
    outcome = pd.Series(np.arange(19, dtype=float), index=idx)
    mp = _MissionPerm(metadata, outcome, None)
    flights = mp.permuted_flights(np.random.default_rng(1), 500)
    # total flights preserved every draw
    assert np.all(flights.sum(axis=1) == 10)
    # per-block flight counts preserved
    for positions, nt in mp.blocks:
        assert np.all(flights[:, positions].sum(axis=1) == nt)


def test_k0_reproduces_frozen_podocyte_scaffold_attenuation():
    """Independent recompute of the decisive observed estimates: the podocyte
    flight coefficient attenuates ~0.198 -> ~0.112 under structural adjustment
    and its adjusted interval crosses zero."""
    config = yaml.safe_load(CONFIG.read_text())
    missions = load_primary_missions(config, REPO)
    tiers = pd.read_csv(TIERS, sep="\t")
    threshold = float(config["eligibility"]["cpm_threshold"])

    podocyte_genes = _frozen_genes(tiers, PODOCYTE_SET)
    structural_disjoint = [
        g for g in _frozen_genes(tiers, STRUCTURAL_SET) if g not in set(podocyte_genes)
    ]

    p0_est, p0_var, p1_est, p1_var = [], [], [], []
    for _, data in missions.items():
        eligible = set(cpm_eligible_genes(data, threshold))
        pod_used = [g for g in podocyte_genes if g in eligible and g in data.expression.index]
        struc_used = [g for g in structural_disjoint if g in eligible and g in data.expression.index]
        podocyte = score_signed_axis(data.expression, {g: 1 for g in pod_used}).scores
        structural = score_signed_axis(data.expression, {g: 1 for g in struc_used}).scores
        hc3_0 = next(
            r for r in _regression_rows("m", data.metadata, podocyte, None, "P0")
            if r["variance_type"] == "HC3"
        )
        hc3_1 = next(
            r for r in _regression_rows("m", data.metadata, podocyte, structural, "P1")
            if r["variance_type"] == "HC3"
        )
        p0_est.append(hc3_0["estimate"]); p0_var.append(hc3_0["variance"])
        p1_est.append(hc3_1["estimate"]); p1_var.append(hc3_1["variance"])

    unadjusted = random_effects_reml_mkh(np.array(p0_est), np.array(p0_var))
    adjusted = random_effects_reml_mkh(np.array(p1_est), np.array(p1_var))

    assert unadjusted.estimate == np.float64(unadjusted.estimate)
    assert abs(unadjusted.estimate - 0.198) < 0.02
    assert unadjusted.ci_low > 0  # unadjusted interval excludes zero
    assert abs(adjusted.estimate - 0.112) < 0.02
    assert adjusted.ci_low < 0 < adjusted.ci_high  # adjusted interval crosses zero
    assert adjusted.estimate < unadjusted.estimate  # attenuation under adjustment
