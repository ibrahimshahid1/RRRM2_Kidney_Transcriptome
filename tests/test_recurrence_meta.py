"""Unit tests for Repair B -- the cross-cohort recurrence meta-analysis."""
import numpy as np
import pandas as pd

from src.v11.recurrence_meta import (
    leave_one_cohort_out,
    meta_random_effects,
    per_cohort_effect,
    resolve_design,
    score_gene_set,
)


# #

def test_resolve_design_splits_flight_ground_and_excludes_viv_baseline_rerun():
    cols = [
        "RR23_KDN_FLT_F1", "RR23_KDN_FLT_F2",
        "RR23_KDN_GC_G1", "RR23_KDN_GC_G2",
        "RR23_KDN_VIV_V1",                       # excluded
        "Mmus_KDN_BSL_Rep1_B1",                  # excluded
        "RR7_KDN_GCrerun_25days_Blue_Rep1_GR01", # excluded
    ]
    d = resolve_design(cols)
    assert d["RR23_KDN_FLT_F1"] == "flight" and d["RR23_KDN_FLT_F2"] == "flight"
    assert d["RR23_KDN_GC_G1"] == "ground" and d["RR23_KDN_GC_G2"] == "ground"
    assert "RR23_KDN_VIV_V1" not in d
    assert "Mmus_KDN_BSL_Rep1_B1" not in d
    assert "RR7_KDN_GCrerun_25days_Blue_Rep1_GR01" not in d   # rerun batch never ground


# #

def _planted_vst(effects, n=6, noise=0.3, seed=0):
    """genes x samples VST with a planted flight-minus-ground shift per gene."""
    rng = np.random.default_rng(seed)
    fcols = [f"X_FLT_{i}" for i in range(n)]
    gcols = [f"X_GC_{i}" for i in range(n)]
    rows = {}
    for gene, eff in effects.items():
        base = 8.0
        f = base + eff + rng.normal(0, noise, n)
        g = base + rng.normal(0, noise, n)
        rows[gene] = np.concatenate([f, g])
    return pd.DataFrame.from_dict(rows, orient="index", columns=fcols + gcols)


def test_per_cohort_effect_is_sign_faithful():
    vst = _planted_vst({"up": +1.0, "down": -1.0, "flat": 0.0}, seed=1)
    res = per_cohort_effect(vst, resolve_design(vst.columns))
    assert res.loc["up", "effect"] > 0.5
    assert res.loc["down", "effect"] < -0.5
    assert abs(res.loc["flat", "effect"]) < 0.5   # flat clearly separates from the +/-1 signals
    assert (res["se"] > 0).all()
    assert (res["n_flight"] == 6).all() and (res["n_ground"] == 6).all()


def test_per_cohort_effect_refuses_too_few_replicates():
    vst = _planted_vst({"g": 1.0}, n=1)
    try:
        per_cohort_effect(vst, resolve_design(vst.columns))
        assert False, "expected ValueError on n=1 per arm"
    except ValueError:
        pass


# #

def _cohort_effect_table(effect, se, genes):
    return pd.DataFrame(
        {"effect": [effect] * len(genes), "se": [se] * len(genes)},
        index=genes,
    )


def test_meta_recovers_shared_effect_with_low_i2_when_homogeneous():
    genes = ["a", "b", "c"]
    # five cohorts, same true effect -0.8, identical SE -> homogeneous
    effects = {f"c{i}": _cohort_effect_table(-0.8, 0.2, genes) for i in range(5)}
    meta = meta_random_effects(effects)
    assert np.allclose(meta["pooled_effect"], -0.8, atol=1e-6)
    assert (meta["I2"] < 1.0).all()                  # essentially no heterogeneity
    assert (meta["ci_high"] < 0).all()               # CI excludes zero on the down side
    assert (meta["k_cohorts"] == 5).all()


def test_meta_flags_heterogeneity_with_high_i2():
    genes = ["g"]
    effects = {
        "c0": _cohort_effect_table(+2.0, 0.1, genes),
        "c1": _cohort_effect_table(-2.0, 0.1, genes),
        "c2": _cohort_effect_table(+2.0, 0.1, genes),
        "c3": _cohort_effect_table(-2.0, 0.1, genes),
    }
    meta = meta_random_effects(effects)
    assert meta.loc["g", "I2"] > 90.0                # cohorts strongly disagree
    assert meta.loc["g", "tau2"] > 0


def test_meta_fdr_separates_shared_signal_from_noise():
    rng = np.random.default_rng(3)
    real = [f"real{i}" for i in range(10)]
    noise = [f"noise{i}" for i in range(90)]
    genes = real + noise
    effects = {}
    for c in range(5):
        eff = np.concatenate([
            np.full(len(real), -1.0),                       # strong shared signal
            rng.normal(0, 0.05, len(noise)),                # ~null
        ])
        effects[f"c{c}"] = pd.DataFrame({"effect": eff, "se": [0.2] * len(genes)}, index=genes)
    meta = meta_random_effects(effects)
    assert (meta.loc[real, "fdr"] < 0.05).all()             # real signal passes
    assert (meta.loc[noise, "fdr"] > 0.05).mean() > 0.8     # most noise does not


def test_meta_requires_minimum_cohorts():
    genes = ["solo"]
    effects = {"only": _cohort_effect_table(-1.0, 0.2, genes)}
    meta = meta_random_effects(effects)
    assert meta.empty                                        # one cohort -> nothing pooled


# #

def _toy_meta(ens_effects):
    rows = []
    for ens, (eff, z, fdr, i2) in ens_effects.items():
        rows.append({
            "pooled_effect": eff, "se_pooled": 0.2, "z": z,
            "p_value": 0.01, "fdr": fdr, "I2": i2,
        })
    return pd.DataFrame(rows, index=pd.Index(list(ens_effects), name="gene_id"))


def test_score_gene_set_is_sign_faithful_and_counts_presence():
    meta = _toy_meta({
        "ENS1": (-1.0, -5.0, 0.001, 3.0),
        "ENS2": (-0.8, -4.0, 0.002, 5.0),
        "ENS3": (+0.1, 0.5, 0.9, 50.0),     # not in the set
    })
    s2e = {"slc12a3": {"ENS1"}, "stk39": {"ENS2"}}
    sc = score_gene_set(meta, ["SLC12A3", "STK39", "ABSENT"], s2e)
    assert sc["n_genes_present"] == 2
    assert sc["set_effect"] < 0                              # coherent suppression
    assert sc["stouffer_z"] < 0
    assert sc["median_fdr"] < 0.05


def test_leave_one_out_reports_baseline_and_each_drop():
    genes = ["ENS1", "ENS2"]
    effects = {f"c{i}": _cohort_effect_table(-0.9, 0.2, genes) for i in range(4)}
    s2e = {"slc12a3": {"ENS1"}, "stk39": {"ENS2"}}
    loo = leave_one_cohort_out(effects, {"transport": ["SLC12A3", "STK39"]}, s2e)
    assert set(loo["dropped_cohort"]) == {"__none__", "c0", "c1", "c2", "c3"}
    # set effect stays negative and stable across drops
    assert (loo["set_effect"] < 0).all()
    spread = loo["set_effect"].max() - loo["set_effect"].min()
    assert spread < 0.2
