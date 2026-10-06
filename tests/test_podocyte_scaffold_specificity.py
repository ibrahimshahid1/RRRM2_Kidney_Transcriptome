"""Contract tests for K0 (podocyte-vs-structural-scaffold specificity)."""

from argparse import Namespace
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
import yaml

from src.clinical_axes.data import MissionData, cpm_eligible_genes, load_primary_missions
from src.clinical_axes.statistics import (
    _batch_reml_mkh,
    random_effects_reml_mkh,
    score_signed_axis,
)
import scripts.clinical_axes.run_podocyte_scaffold_specificity as k0
from scripts.clinical_axes.run_podocyte_scaffold_specificity import (
    FAMILY,
    FREEDMAN_LANE_SCHEME,
    JOINT_LABEL_SCHEME,
    LEGACY_SCHEME,
    PERMUTATION_SCHEMES,
    PODOCYTE_SET,
    STRUCTURAL_SET,
    _design,
    _frozen_genes,
    _MissionPerm,
    _regression_rows,
    family_permutation_null,
    observed_family_t,
    permutation_table,
)

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"

# Synthetic mission layouts mirroring the five terminal missions:
# (block, n_flight, n_control) per block.
LAYOUTS = [
    [("all", 6, 6)],
    [("all", 6, 6)],
    [("day25", 5, 5), ("day75", 5, 4)],
    [("all", 10, 10)],
    [("YNG", 5, 5), ("OLD", 5, 5)],
]


def _k0_real_inputs_present() -> bool:
    if not (CONFIG.exists() and TIERS.exists()):
        return False
    config = yaml.safe_load(CONFIG.read_text())
    paths = [REPO / config["gene_mapping"]["path"]]
    for spec in config["primary_missions"].values():
        paths.extend(REPO / str(spec[key]) for key in ("vst", "counts") if key in spec)
    return all(path.exists() for path in paths)


def _metadata(layout, prefix="a") -> pd.DataFrame:
    condition, block = [], []
    for name, n_flight, n_control in layout:
        condition += ["FLT"] * n_flight + ["GC"] * n_control
        block += [name] * (n_flight + n_control)
    index = [f"{prefix}{i}" for i in range(len(condition))]
    return pd.DataFrame({"condition": condition, "block": block}, index=index)


def _synthetic_family(rng, *, flight_on_adjuster=1.0, adjuster_slope=0.8, effect=0.0):
    """P1/S1/D-like perm objects on synthetic missions with a flight-shifted
    adjuster (S0 > 0) and no direct flight effect on the outcome (H0 for P1)."""
    objects = {model: [] for model in FAMILY}
    for layout in LAYOUTS:
        md = _metadata(layout)
        flight = (md["condition"] == "FLT").to_numpy(float)
        block = pd.factorize(md["block"])[0].astype(float)
        structural = flight_on_adjuster * flight + 0.5 * block + rng.normal(size=len(md))
        podocyte = effect * flight + adjuster_slope * structural + 0.3 * block + rng.normal(size=len(md))
        pod = pd.Series(podocyte, index=md.index)
        struc = pd.Series(structural, index=md.index)
        objects["P1_podocyte_adjusted"].append(_MissionPerm(md, pod, struc))
        objects["S1_structural_adjusted"].append(_MissionPerm(md, struc, pod))
        objects["D_direct_difference"].append(_MissionPerm(md, pod - struc, None))
    return objects


# --- pre-correction reference implementation (verbatim logic of HEAD 26c1890) ---

def _old_permuted_flights(obj, rng, batch):
    out = np.zeros((batch, obj.n), dtype=float)
    for positions, nt in obj.blocks:
        m = len(positions)
        keys = rng.random((batch, m))
        chosen = np.argpartition(keys, nt - 1, axis=1)[:, :nt] if nt > 0 else np.empty((batch, 0), int)
        for b in range(batch):
            out[b, positions[chosen[b]]] = 1.0
    return out


def _old_null_and_pvalues(perm_objects, n_perm, seed, chunk):
    """Copy of the old run() permutation block (one RNG per (model, mission))."""
    family = list(perm_objects)
    n_missions = len(perm_objects[family[0]])
    rng_root = np.random.SeedSequence(seed)
    observed_t = {}
    for model in family:
        objs = perm_objects[model]
        betas = np.array([o.coeff_var(o.flight_obs[None, :])[0][0] for o in objs])
        varis = np.array([o.coeff_var(o.flight_obs[None, :])[1][0] for o in objs])
        observed_t[model] = random_effects_reml_mkh(betas, varis).t
    null_t = {model: np.empty(n_perm) for model in family}
    child_seeds = rng_root.spawn(len(family) * n_missions)
    rngs = {}
    si = 0
    for model in family:
        for mi in range(n_missions):
            rngs[(model, mi)] = np.random.default_rng(child_seeds[si])
            si += 1
    for model in family:
        objs = perm_objects[model]
        for start in range(0, n_perm, chunk):
            stop = min(start + chunk, n_perm)
            batch = stop - start
            betas = np.empty((len(objs), batch))
            varis = np.empty((len(objs), batch))
            for mi, o in enumerate(objs):
                flights = _old_permuted_flights(o, rngs[(model, mi)], batch)
                b, v = o.coeff_var(flights)
                betas[mi] = b
                varis[mi] = v
            _, _, _, t_vals = _batch_reml_mkh(betas.T, varis.T)
            null_t[model][start:stop] = t_vals
    null_max = np.max(np.abs(np.column_stack([null_t[m] for m in family])), axis=1)
    rows = []
    for model in family:
        obs = observed_t[model]
        emp = (1.0 + np.count_nonzero(np.abs(null_t[model]) >= abs(obs))) / (n_perm + 1.0)
        fwer = (1.0 + np.count_nonzero(null_max >= abs(obs))) / (n_perm + 1.0)
        rows.append(
            {
                "model": model,
                "observed_t_model_based": obs,
                "empirical_p_two_sided": emp,
                "max_t_fwer": fwer,
                "n_permutations": n_perm,
            }
        )
    return null_t, pd.DataFrame(rows)


# --- existing contracts --------------------------------------------------------

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


@pytest.mark.skipif(
    not _k0_real_inputs_present(),
    reason="K0 real-data inputs (config, id map, tiers, OSDR matrices) not present",
)
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


# --- WP-K0c: permutation correction ----------------------------------------------

@pytest.mark.parametrize("n_perm,chunk", [(5000, 2048), (1500, 700), (64, 2048)])
def test_legacy_scheme_reproduces_pre_correction_implementation(n_perm, chunk):
    """legacy_independent_streams is bit-for-bit the old code path (seed 0),
    including multi-chunk and partial-chunk draws."""
    objects = _synthetic_family(np.random.default_rng(11))
    old_null, old_table = _old_null_and_pvalues(objects, n_perm, 0, chunk)
    null = family_permutation_null(
        objects, n_permutations=n_perm, seed=0, chunk_size=chunk, schemes=(LEGACY_SCHEME,)
    )[LEGACY_SCHEME]
    for model in FAMILY:
        assert np.array_equal(null[model].to_numpy(), old_null[model])
    table = permutation_table(observed_family_t(objects), {LEGACY_SCHEME: null}, seed=0)
    pd.testing.assert_frame_equal(
        table[list(old_table.columns)].reset_index(drop=True), old_table, check_exact=True
    )


def test_vectorized_permuted_flights_equal_old_loop():
    objects = _synthetic_family(np.random.default_rng(3))
    for obj in objects["P1_podocyte_adjusted"]:
        new = obj.permuted_flights(np.random.default_rng(99), 333)
        old = _old_permuted_flights(obj, np.random.default_rng(99), 333)
        assert np.array_equal(new, old)


def test_identical_models_share_null_only_under_joint_schemes():
    """With two identical family models the joint schemes give identical null
    columns (same permutation per mission); legacy streams do not."""
    base = _synthetic_family(np.random.default_rng(5))["P1_podocyte_adjusted"]
    objects = {"A": base, "B": list(base)}
    null = family_permutation_null(objects, n_permutations=400, seed=0, chunk_size=128)
    for scheme in (JOINT_LABEL_SCHEME, FREEDMAN_LANE_SCHEME):
        assert np.array_equal(null[scheme]["A"].to_numpy(), null[scheme]["B"].to_numpy())
    assert not np.array_equal(
        null[LEGACY_SCHEME]["A"].to_numpy(), null[LEGACY_SCHEME]["B"].to_numpy()
    )
    # identical models -> max over the joint family equals the single-model |t|
    table = permutation_table(observed_family_t(objects), null, seed=0)
    joint = table[table["permutation_scheme"] != LEGACY_SCHEME]
    assert np.allclose(joint["max_t_fwer"], joint["empirical_p_two_sided"])


def test_freedman_lane_identity_permutation_reproduces_model_based_ols_t():
    """pi = identity in Freedman-Lane is the observed full-model OLS fit: same
    flight coefficient, model-based SE and t as statsmodels, for adjusted and
    unadjusted designs in every synthetic mission layout."""
    rng = np.random.default_rng(8)
    for layout in LAYOUTS:
        md = _metadata(layout)
        flight = (md["condition"] == "FLT").to_numpy(float)
        adjuster = pd.Series(rng.normal(size=len(md)) + flight, index=md.index)
        outcome = pd.Series(rng.normal(size=len(md)) + 0.7 * adjuster, index=md.index)
        for adj in (adjuster, None):
            obj = _MissionPerm(md, outcome, adj)
            fit = sm.OLS(outcome, _design(md, adj)).fit()
            beta, var = obj.freedman_lane_coeff_var(np.arange(obj.n)[None, :])
            assert np.isclose(beta[0], fit.params["flight"], rtol=0, atol=1e-10)
            assert np.isclose(np.sqrt(var[0]), fit.bse["flight"], rtol=0, atol=1e-10)
            assert np.isclose(beta[0] / np.sqrt(var[0]), fit.tvalues["flight"], rtol=0, atol=1e-8)
            # equals the permute-X closed form at the observed labels
            b_x, v_x = obj.coeff_var(obj.flight_obs[None, :])
            assert np.isclose(beta[0], b_x[0], rtol=0, atol=1e-12)
            assert np.isclose(var[0], v_x[0], rtol=0, atol=1e-12)
    # pooled: the observed family t (used for every scheme) is the identity draw
    objects = _synthetic_family(np.random.default_rng(21))
    observed = observed_family_t(objects)
    for model, objs in objects.items():
        identity = [o.freedman_lane_coeff_var(np.arange(o.n)[None, :]) for o in objs]
        betas = np.array([[b[0] for b, _ in identity]])
        varis = np.array([[v[0] for _, v in identity]])
        _, _, _, t_vals = _batch_reml_mkh(betas, varis)
        assert np.isclose(t_vals[0], observed.loc[model, "observed_t_model_based"], atol=1e-6)


def test_closed_forms_match_brute_force_refits_for_random_permutations():
    """Freedman-Lane: y* = Z g_hat + e[pi] refitted on the full design; permute-X:
    y refitted with the permuted flight column. Both closed forms must equal the
    statsmodels refit coefficient and model-based variance."""
    rng = np.random.default_rng(77)
    md = _metadata(LAYOUTS[4])
    flight = (md["condition"] == "FLT").to_numpy(float)
    adjuster = pd.Series(rng.normal(size=len(md)) + 1.2 * flight, index=md.index)
    outcome = pd.Series(rng.normal(size=len(md)) + 0.9 * adjuster, index=md.index)
    obj = _MissionPerm(md, outcome, adjuster)
    design = _design(md, adjuster)
    reduced = design.drop(columns=["flight"])
    fitted_z = sm.OLS(outcome, reduced).fit().fittedvalues.to_numpy()
    resid_z = outcome.to_numpy() - fitted_z
    perms = obj.permutation_indices(np.random.default_rng(3), 25)
    fl_beta, fl_var = obj.freedman_lane_coeff_var(perms)
    x_beta, x_var = obj.coeff_var(obj.flights_for_permutation(perms))
    for b, perm in enumerate(perms):
        y_star = fitted_z + resid_z[perm]
        fit = sm.OLS(y_star, design).fit()
        assert np.isclose(fl_beta[b], fit.params["flight"], rtol=0, atol=1e-10)
        assert np.isclose(fl_var[b], fit.bse["flight"] ** 2, rtol=1e-9, atol=0)
        permuted = design.copy()
        permuted["flight"] = obj.flights_for_permutation(perm[None, :])[0]
        fit_x = sm.OLS(outcome, permuted).fit()
        assert np.isclose(x_beta[b], fit_x.params["flight"], rtol=0, atol=1e-10)
        assert np.isclose(x_var[b], fit_x.bse["flight"] ** 2, rtol=1e-9, atol=0)


def test_permutation_indices_preserve_block_membership_and_match_label_draws():
    md = _metadata([("d25", 5, 5), ("d75", 5, 4), ("x", 3, 2)])
    rng = np.random.default_rng(0)
    obj = _MissionPerm(md, pd.Series(rng.normal(size=len(md)), index=md.index), None)
    perm = obj.permutation_indices(np.random.default_rng(42), 2000)
    for row in perm[:200]:
        assert sorted(row) == list(range(obj.n))  # bijection
    for positions, nt in obj.blocks:
        # every unit of a block is mapped inside the same block
        assert np.all(np.isin(perm[:, positions], positions))
    flights = obj.flights_for_permutation(perm)
    for positions, nt in obj.blocks:
        assert np.all(flights[:, positions].sum(axis=1) == nt)
    # same generator state -> same flight allocation as permuted_flights
    assert np.array_equal(flights, obj.permuted_flights(np.random.default_rng(42), 2000))
    # permutations are not degenerate: every unit of a block can receive a flight label
    for positions, nt in obj.blocks:
        assert np.all(flights[:, positions].mean(axis=0) > 0.3)


def test_freedman_lane_equals_joint_label_without_continuous_adjuster():
    """With Z = const + block only (model D), Freedman-Lane and permute-X are the
    same statistic for the shared permutation (residuals stay block-centred)."""
    objects = _synthetic_family(np.random.default_rng(13))
    null = family_permutation_null(
        objects, n_permutations=600, seed=4, chunk_size=256,
        schemes=(JOINT_LABEL_SCHEME, FREEDMAN_LANE_SCHEME),
    )
    np.testing.assert_allclose(
        null[FREEDMAN_LANE_SCHEME]["D_direct_difference"],
        null[JOINT_LABEL_SCHEME]["D_direct_difference"],
        rtol=1e-9, atol=1e-10,
    )
    # adjusted models genuinely differ between schemes
    assert not np.allclose(
        null[FREEDMAN_LANE_SCHEME]["P1_podocyte_adjusted"],
        null[JOINT_LABEL_SCHEME]["P1_podocyte_adjusted"],
    )


def test_joint_streams_reuse_legacy_p1_streams():
    """SeedSequence(seed).spawn(n_missions) equals the legacy P1 streams, so the
    joint_label P1 marginal null is bit-identical to legacy's; only the coupling
    across models (and hence the FWER) changes."""
    objects = _synthetic_family(np.random.default_rng(17))
    null = family_permutation_null(objects, n_permutations=3000, seed=0, chunk_size=2048)
    assert np.array_equal(
        null[LEGACY_SCHEME]["P1_podocyte_adjusted"].to_numpy(),
        null[JOINT_LABEL_SCHEME]["P1_podocyte_adjusted"].to_numpy(),
    )
    assert not np.array_equal(
        null[LEGACY_SCHEME]["S1_structural_adjusted"].to_numpy(),
        null[JOINT_LABEL_SCHEME]["S1_structural_adjusted"].to_numpy(),
    )


def test_permutation_table_has_three_schemes_and_declares_primary():
    objects = _synthetic_family(np.random.default_rng(23))
    null = family_permutation_null(objects, n_permutations=199, seed=0, chunk_size=64)
    table = permutation_table(observed_family_t(objects), null, seed=0)
    assert list(table["permutation_scheme"].unique()) == list(PERMUTATION_SCHEMES)
    assert len(table) == 3 * len(FAMILY)
    assert set(table.loc[table["primary_scheme"], "permutation_scheme"]) == {FREEDMAN_LANE_SCHEME}
    # observed statistic identical across schemes; FWER >= marginal p within scheme
    assert table.groupby("model")["observed_t_model_based"].nunique().eq(1).all()
    assert np.all(table["max_t_fwer"] >= table["empirical_p_two_sided"] - 1e-15)
    assert "observed_coef_score_units_model_based" in table.columns


def test_mission_perm_rejects_models_with_different_animals():
    objects = _synthetic_family(np.random.default_rng(1))
    other = _synthetic_family(np.random.default_rng(2))
    bad = {"A": objects["P1_podocyte_adjusted"], "B": other["P1_podocyte_adjusted"][::-1]}
    with pytest.raises(ValueError):
        family_permutation_null(bad, n_permutations=10, seed=0)


def test_freedman_lane_is_calibrated_with_flight_correlated_covariate():
    """Synthetic null: flight shifts the adjuster (S0 > 0) and the outcome tracks
    the adjuster, but flight has no direct effect. 300 simulations x 199
    permutations; FL rejection rate at 0.05 must lie in [0.02, 0.09]."""
    rng = np.random.default_rng(20261006)
    n_sim, n_perm = 300, 199
    rejections = 0
    for _ in range(n_sim):
        objects = {"P1": _synthetic_family(rng)["P1_podocyte_adjusted"]}
        null = family_permutation_null(
            objects, n_permutations=n_perm, seed=int(rng.integers(2**31)),
            chunk_size=n_perm, schemes=(FREEDMAN_LANE_SCHEME,),
        )
        table = permutation_table(observed_family_t(objects), null, seed=0)
        rejections += bool(table["empirical_p_two_sided"].iloc[0] <= 0.05)
    rate = rejections / n_sim
    assert 0.02 <= rate <= 0.09, rate


def test_freedman_lane_has_power_against_a_direct_flight_effect():
    """Guard against a trivially conservative implementation."""
    rng = np.random.default_rng(5)
    rejections = 0
    for _ in range(40):
        objects = {"P1": _synthetic_family(rng, effect=1.0)["P1_podocyte_adjusted"]}
        null = family_permutation_null(
            objects, n_permutations=199, seed=int(rng.integers(2**31)),
            chunk_size=199, schemes=(FREEDMAN_LANE_SCHEME,),
        )
        table = permutation_table(observed_family_t(objects), null, seed=0)
        rejections += bool(table["empirical_p_two_sided"].iloc[0] <= 0.05)
    assert rejections / 40 >= 0.5


# --- end-to-end synthetic run --------------------------------------------------------

def _synthetic_missions(rng, genes):
    missions = {}
    for mi, layout in enumerate(LAYOUTS):
        md = _metadata(layout, prefix=f"m{mi}_")
        md["source_sample"] = md.index
        n = len(md)
        flight = (md["condition"] == "FLT").to_numpy(float)
        latent = rng.normal(size=n) + 0.5 * flight
        expression = pd.DataFrame(
            rng.normal(size=(len(genes), n)) + latent[None, :] + 8.0,
            index=pd.Index(genes, name="gene"), columns=md.index,
        )
        counts = pd.DataFrame(
            rng.poisson(500, size=(len(genes), n)).astype(float),
            index=pd.Index(genes, name="gene"), columns=md.index,
        )
        data = MissionData(
            mission=f"SYN-{mi}", expression=expression, counts=counts,
            metadata=md[["condition", "block", "source_sample"]],
            qc=pd.DataFrame(index=md.index),
        )
        data.validate()
        missions[f"SYN-{mi}"] = data
    return missions


def _write_inputs(tmp_path):
    pod = [f"Pod{i}" for i in range(12)]
    struc = [f"Str{i}" for i in range(14)] + ["Pod0", "Pod1"]
    rows = [(PODOCYTE_SET, g) for g in pod] + [(STRUCTURAL_SET, g) for g in struc]
    tiers = pd.DataFrame(rows, columns=["gene_set", "gene_symbol"])
    tiers["final_for_testing"] = True
    tiers_path = tmp_path / "tiers.tsv"
    tiers.to_csv(tiers_path, sep="\t", index=False)
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG0\tPod0\n")
    alt_map = tmp_path / "alt_id_map.tsv"
    alt_map.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG0\tPod0\n")
    config = {
        "seed": 1,
        "eligibility": {"cpm_threshold": 0.1},
        "gene_mapping": {"path": str(id_map), "annotation_fallback": "unused.csv"},
    }
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump(config))
    genes = sorted(set(pod) | set(struc))
    return tiers_path, config_path, alt_map, genes


def test_run_end_to_end_on_synthetic_missions(tmp_path, monkeypatch):
    tiers_path, config_path, alt_map, genes = _write_inputs(tmp_path)
    missions = _synthetic_missions(np.random.default_rng(31), genes)
    seen = {}

    def fake_loader(config, root):
        seen["gene_map"] = config["gene_mapping"]["path"]
        return missions

    monkeypatch.setattr(k0, "load_primary_missions", fake_loader)
    config_before = config_path.read_text()
    out = tmp_path / "out"
    args = Namespace(
        config=config_path, tiers=tiers_path, results=out, permutations=150,
        seed=0, chunk_size=64, gene_map=alt_map,
    )
    k0.run(args)

    # --gene-map replaces the path in memory only
    assert seen["gene_map"] == str(alt_map.resolve())
    assert config_path.read_text() == config_before

    perm = pd.read_csv(out / "podocyte_scaffold_permutation.tsv", sep="\t")
    assert list(perm["permutation_scheme"].unique()) == list(PERMUTATION_SCHEMES)
    for scheme in PERMUTATION_SCHEMES:
        assert list(perm.loc[perm["permutation_scheme"] == scheme, "model"]) == FAMILY
    for column in ["model", "observed_t_model_based", "empirical_p_two_sided",
                   "max_t_fwer", "n_permutations"]:
        assert column in perm.columns  # legacy columns kept

    lomo = pd.read_csv(out / "podocyte_scaffold_leave_one_mission.tsv", sep="\t")
    assert set(lomo["model"]) == {
        "P1_podocyte_adjusted", "S1_structural_adjusted", "D_direct_difference"
    }
    assert (lomo.groupby("model").size() == len(LAYOUTS)).all()

    for name in ["podocyte_scaffold_meta.tsv", "podocyte_scaffold_mission_effects.tsv",
                 "podocyte_scaffold_leave_one_mission.tsv"]:
        table = pd.read_csv(out / name, sep="\t")
        assert np.array_equal(table["coef_score_units"], table["estimate"])
        assert not any(column.startswith("g") and column != "gene" for column in table.columns)

    manifest = json.loads((out / "podocyte_scaffold_manifest.json").read_text())
    assert set(manifest["permutation_schemes"]) == set(PERMUTATION_SCHEMES)
    assert manifest["primary_permutation_scheme"] == FREEDMAN_LANE_SCHEME
    assert all(rec["seed"] == 0 for rec in manifest["permutation_schemes"].values())
    assert manifest["permutation_schemes"][FREEDMAN_LANE_SCHEME]["primary"] is True
    assert manifest["gene_map"] == str(alt_map.resolve())
    assert manifest["gene_map_overridden"] is True
    assert "D_direct_difference" in manifest["leave_one_mission_models"]
