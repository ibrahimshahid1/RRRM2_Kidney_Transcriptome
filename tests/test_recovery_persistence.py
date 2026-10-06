"""Contracts for G2 recovery/persistence (plan §3 WP-P).

Everything here runs on synthetic data except three skip-if-absent data contracts that
touch only already-known quantities (blinding): the ISS-T podocyte mission effect
(ISS-T is the primary terminal arm) and the LAR / OSD-513 four-clinical-axis effects.
No LAR or OSD-513 compartment, podocyte, structural or K0 quantity is computed.
"""

from __future__ import annotations

from pathlib import Path
import math
import subprocess

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
import yaml
from scipy import stats

from conftest import requires_paths
from src.clinical_axes import persistence as P
from src.clinical_axes.analysis import mission_effect_from_score
from src.clinical_axes.data import MissionData
from src.clinical_axes.statistics import blocked_meta_permutation
from scripts.clinical_axes import run_recovery_persistence as R

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"


def _configured_id_map() -> Path:
    """The baseline gene-ID map named by the config (data tests skip when it is absent)."""
    try:
        return REPO / yaml.safe_load(CONFIG.read_text())["gene_mapping"]["path"]
    except (OSError, KeyError, TypeError):
        return REPO / "data/processed/resources/id_map.tsv"


ID_MAP = _configured_id_map()
OSD771_VST = REPO / "data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
OSD513_VST = REPO / "data/external/osdr/OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"

# Plan §4.2 verbatim: the temporary preregistration every rule test parses.
PREREG_4_2 = """
analysis_id: clinical_axes_recovery_persistence_v1
parent_config: config/clinical_renal_axes_cross_mission.yaml
status: preregistered_secondary_analysis
lock_date: null
seed_offsets: {freedman_lane: 700, bootstrap: 701, power_simulation: 702, posterior: 703}
cohorts:
  primary_within_mission: {mission: OSD-771, arms: [ISS-T, LAR], groups: [FLT, GC], blocks: arm_x_age,
                           n_per_cell: 5, excluded_groups: [BSL, VIV]}
  descriptive_single_cohort: {mission: OSD-513, n_flight: 9, n_control: 9}
  not_done: [REML/mKH across recovery cohorts (k=2, df=1), confirmatory 7-cohort meta-regression]
targets:
  primary: podocyte__high_specificity
  key_secondary: [structural_disjoint, D_podocyte_minus_structural_disjoint]
  secondary: [P1_type_adjusted_podocyte]
  exploratory: all 49 evaluable compartment sets; 4 frozen axes as context
eligibility: genes CPM-eligible in BOTH arms (each arm's own FLT/GC rule) and in VST; set needs >=8
scoring: mean signed gene z; z over the 40-animal joint analysis set (never within blocks)
model: score ~ 0 + block(arm x age) + flight_ISS-T + flight_LAR; OLS; HC3; per-arm df 17; Satterthwaite for contrasts
contrasts: {interaction: theta_LAR - theta_ISS-T, margin: theta_LAR - 0.5*theta_ISS-T}
margin_sweep: {kappas: [0.0, 0.25, 0.5, 0.75, 1.0], role: secondary_descriptive}
inference:
  permutation: Freedman-Lane residual permutation within arm x age blocks, 100000, same permutation across sets
  single_arm: exact enumeration (63,504 per OSD-771 arm; 48,620 OSD-513), exact maxT across the 49-set family
  ratio: Fieller 90/95% (primary interval for rho); stratified bootstrap 10000 within 8 cells (secondary)
  bayesian_descriptive: "flat-prior normal draws: P(rho>=0.5), P(theta_LAR<=0)"  # quoted: unquoted §4.2 text is invalid YAML
classification:   # applied to primary and key secondary; exploratory labels are descriptive
  NOT_EVALUABLE_TECHNICAL: LAR primary coverage metric flagged by the frozen QC rule
  NOT_EVALUABLE_NO_TERMINAL_EFFECT: HC3 95% CI of theta_ISS-T includes 0, or sign differs from pooled terminal
  PERSISTS: FL one-sided p(margin>0) <= 0.05 AND HC3 90% CI lower(margin) > 0
  DECAYS: FL one-sided p(margin<0) <= 0.05 AND HC3 90% CI upper(margin) < 0  (sublabel REVERSED if theta_LAR 95% CI upper < 0)
  INCONCLUSIVE: otherwise, including FL/HC3 disagreement
multiplicity:
  primary: single target, no correction
  key_secondary: Holm across the 2 margin tests per direction; interaction maxT over {T1,T2,T3}
  exploratory: interaction maxT + BH over 49 sets; single-arm exact maxT; labels descriptive
confounds_stated: [recovery aliased with 32 vs 53-56 d flight, on-ISS vs on-Earth euthanasia/handling,
                   re-entry and landing; winner's-curse inflation of theta_ISS-T for the selected podocyte set
                   biases rho downward; OSD-513 is ~1 d post-landing with a 3 d ground-control offset]
power_analytic:   # SD units of the score; normal approximation; simulation in persistence_power.tsv
  se_theta_per_arm: 0.447
  se_interaction: 0.632
  mde80_interaction_two_sided: 1.77
  power_complete_reversal: {theta_0.689: 0.19, theta_1.174: 0.46}
  se_margin: 0.50
  power_PERSISTS_if_full_persistence: {theta_0.689: 0.17, theta_1.174: 0.32}
  power_DECAYS_if_complete_reversal: {theta_0.689: 0.17, theta_1.174: 0.32}
  theta_needed_for_80pct_margin_power: 2.49
  expected_modal_label: INCONCLUSIVE (approx 0.68-0.83)
  single_arm_LAR_mde80: 1.25    # power 0.34 / 0.75 if theta_LAR = 0.689 / 1.174
  osd513_mde80: 1.32
provenance_disclosure: >-
  Known before lock: LAR and OSD-513 4-axis effects. Not computed by anyone before lock:
  any LAR/OSD-513 compartment, podocyte, structural or K0 effect.
"""


def _rules() -> P.PersistenceRules:
    return P.parse_classification_rules(yaml.safe_load(PREREG_4_2), source="test")


# --------------------------------------------------------------------------- synthetic cohorts

GENES = [f"g{i:03d}" for i in range(260)]
SETS = {
    P.PODOCYTE_SET: GENES[0:40],
    P.STRUCTURAL_FULL_SET: GENES[35:100],  # 5 genes shared with the podocyte set
    "tubule__high_specificity": GENES[100:130],
    "immune__high_specificity": GENES[130:160],
    "endothelial__high_specificity": GENES[160:190],
    "small__high_specificity": GENES[190:199],  # 9 genes; 3 are dead in LAR only
    "tiny__undefined": GENES[199:205],  # < 8 defined genes: never in the family
}
AXES = {
    "axis_barrier": {
        "role": "primary_biological_axis",
        "subdomains": {"core": {"genes": {g: -1 for g in GENES[210:216]}, "minimum_present": 4}},
    },
    "axis_two_domains": {
        "role": "primary_biological_axis",
        "subdomains": {
            "a": {"genes": {g: 1 for g in GENES[220:224]}, "minimum_present": 2},
            "b": {"genes": {g: 1 for g in GENES[224:230]}, "minimum_present": 4},
        },
    },
}
CONFIG_SYNTH = {
    "seed": 20260811,
    "eligibility": {"cpm_threshold": 0.1},
    "technical_qc": {"primary_metric": "ratio_genebody_cov_3_to_5"},
    "primary_family": AXES,
}


def _mission(
    name: str,
    design: list[tuple[str, int, int]],
    rng: np.random.Generator,
    effects: dict[str, float],
    *,
    context: str = "terminal",
    prefix: str | None = None,
    dead: tuple[str, ...] = (),
    qc_shift: float = 0.0,
) -> MissionData:
    animals, condition, block = [], [], []
    for stratum, n_flt, n_gc in design:
        for label, count in (("FLT", n_flt), ("GC", n_gc)):
            for rep in range(count):
                animals.append(f"{prefix or name}_{stratum}_{label}_{rep}")
                condition.append(label)
                block.append(stratum)
    n = len(animals)
    meta = pd.DataFrame({"condition": condition, "block": block, "source_sample": animals}, index=animals)
    flight = np.array([c == "FLT" for c in condition], dtype=float)
    base = rng.normal(8.0, 1.0, size=(len(GENES), 1))
    animal_shift = rng.normal(0.0, 0.3, size=(1, n))  # shared animal factor, as in real VST data
    expr = base + animal_shift + rng.normal(0.0, 0.6, size=(len(GENES), n))
    for set_name, effect in effects.items():
        rows = [GENES.index(g) for g in SETS[set_name]]
        expr[rows] += effect * flight[None, :] * rng.uniform(0.5, 1.5, size=(len(rows), 1))
    counts = rng.poisson(300, size=(len(GENES), n)).astype(float)
    counts[[GENES.index(g) for g in dead]] = 0.0
    qc = pd.DataFrame(
        {
            "ratio_genebody_cov_3_to_5": rng.normal(1.0, 0.05, n) + qc_shift * flight,
            "rin": rng.normal(8.0, 0.3, n),
            "read_depth": rng.normal(5e7, 5e6, n),
            "uniquely_mapped_percent": rng.normal(85.0, 2.0, n),
        },
        index=animals,
    )
    data = MissionData(
        mission=name,
        expression=pd.DataFrame(expr, index=pd.Index(GENES, name="gene"), columns=animals),
        counts=pd.DataFrame(counts, index=pd.Index(GENES, name="gene"), columns=animals),
        metadata=meta,
        qc=qc,
        context=context,
    )
    data.validate()
    return data


def synthetic_cohorts(seed: int = 7, *, lar_qc_shift: float = 0.0, lar_effect: float = 0.5) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    age = [("YNG", 5, 5), ("OLD", 5, 5)]
    terminal = {P.PODOCYTE_SET: 0.9, P.STRUCTURAL_FULL_SET: 0.5, "tubule__high_specificity": -0.6}
    iss = _mission("OSD-771", age, rng, terminal, prefix="ISS")
    lar = _mission(
        "OSD-771-LAR",
        age,
        rng,
        {P.PODOCYTE_SET: lar_effect, P.STRUCTURAL_FULL_SET: 0.2},
        context="long_recovery",
        prefix="LAR",
        dead=tuple(SETS["small__high_specificity"][:3]),
        qc_shift=lar_qc_shift,
    )
    osd513 = _mission("OSD-513", [("all", 9, 9)], rng, {"immune__high_specificity": 0.4}, context="live_return")
    primary = {
        "OSD-102": _mission("OSD-102", [("all", 6, 6)], rng, terminal),
        "OSD-163": _mission("OSD-163", [("all", 6, 6)], rng, terminal),
        "OSD-253": _mission("OSD-253", [("day25", 5, 5), ("day75", 5, 4)], rng, terminal),
        "OSD-462": _mission("OSD-462", [("all", 10, 10)], rng, terminal),
        "OSD-771": iss,
    }
    P.check_persistence_cohorts(iss, lar, osd513)
    return {"iss": iss, "lar": lar, "osd513": osd513, "primary": primary}


def synthetic_tiers(path: Path) -> Path:
    rows = []
    for gene_set, genes in SETS.items():
        for gene in genes:
            rows.append(
                {
                    "gene_set": gene_set,
                    "gene_symbol": gene,
                    "final_for_testing": True,
                    "report_compartment": gene_set.split("__")[0],
                    "tier": gene_set.split("__")[-1],
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def synthetic_settings(**overrides) -> R.Settings:
    settings = R.prereg_settings(yaml.safe_load(PREREG_4_2), CONFIG_SYNTH, n_perm=199, n_boot=200, posterior_draws=20_000)
    return settings if not overrides else R.Settings(**{**settings.__dict__, **overrides})


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory):
    tiers = synthetic_tiers(tmp_path_factory.mktemp("tiers") / "tiers.tsv")
    cohorts = synthetic_cohorts()
    settings = synthetic_settings()
    defs = R.build_set_definitions(tiers, cohorts["primary"], AXES, settings.cpm_threshold, settings.min_genes)
    result = R.analyze(cohorts, defs, CONFIG_SYNTH, _rules(), settings, log=lambda _msg: None)
    return {"tiers": tiers, "cohorts": cohorts, "settings": settings, "defs": defs, "result": result}


# --------------------------------------------------------------------------- allocations / exact tests


def test_enumeration_gives_63504_unique_block_preserving_allocations_including_observed():
    meta = P.synthetic_meta40(5)
    iss = meta[meta["arm"] == "ISS-T"]
    observed = iss["condition"].eq("FLT").to_numpy()
    alloc = P.enumerate_block_allocations(iss["block"], {"ISS-T_YNG": 5, "ISS-T_OLD": 5}, observed=observed)
    assert alloc.exact and alloc.masks.shape == (63_504, 20) and alloc.n_total == 63_504
    assert len(np.unique(alloc.masks, axis=0)) == 63_504
    for block in ("ISS-T_YNG", "ISS-T_OLD"):
        assert np.all(alloc.masks[:, iss["block"].eq(block).to_numpy()].sum(axis=1) == 5)
    assert np.array_equal(alloc.masks[alloc.observed_index], observed)


def test_enumeration_osd513_and_monte_carlo_fallback():
    labels = ["all"] * 18
    observed = np.array([True] * 9 + [False] * 9)
    exact = P.enumerate_block_allocations(labels, {"all": 9}, observed=observed)
    assert exact.masks.shape == (48_620, 18) and exact.exact
    mc = P.enumerate_block_allocations(labels, {"all": 9}, observed=observed, max_exact=1000, n_monte_carlo=500, seed=3)
    assert not mc.exact and mc.masks.shape == (501, 18) and mc.observed_index == 0
    assert np.array_equal(mc.masks[0], observed) and np.all(mc.masks.sum(axis=1) == 9)


def test_arm_exact_test_observed_allocation_equals_mission_effect_from_score():
    cohorts = synthetic_cohorts()
    family = {name: P.single_set_spec(genes, role="t") for name, genes in SETS.items() if len(genes) >= 8}
    for data in (cohorts["iss"], cohorts["osd513"]):
        scores, _ = P.score_arm(data, family, 0.1)
        table = P.arm_exact_test(scores, data.metadata, {"all": list(scores.columns)}).set_index("set")
        n_alloc = int(table["n_allocations"].iloc[0])
        for name in scores.columns:
            summary, _ = mission_effect_from_score(scores[name], data.metadata)
            assert abs(table.loc[name, "estimate_g"] - summary["estimate"]) < 1e-10
            assert abs(table.loc[name, "variance"] - summary["variance"]) < 1e-10
            # the mirror allocation (labels swapped in every block) ties |z|, so p >= 2/N
            assert table.loc[name, "p_exact"] >= 2.0 / n_alloc - 1e-15
            assert table.loc[name, "max_t_fwer_exact"] >= table.loc[name, "p_exact"] - 1e-15


# --------------------------------------------------------------------------- HC3 / Freedman--Lane


def _joint_synthetic(rng, n_sets=3, theta_t=1.0, theta_l=0.5):
    meta = P.synthetic_meta40(5)
    design = P.design_matrix(meta)
    y = (
        design[list(P.BLOCK_COLUMNS)].to_numpy() @ rng.normal(0, 1, 4)
    )[:, None] + theta_t * design["F_T"].to_numpy()[:, None] + theta_l * design["F_L"].to_numpy()[:, None]
    hetero = np.where(meta["arm"].eq("LAR"), 1.7, 1.0)[:, None]
    y = y + rng.standard_normal((40, n_sets)) * hetero
    return meta, design, y


def test_ols_hc3_matches_statsmodels_and_arm_covariance_vanishes():
    meta, design, y = _joint_synthetic(np.random.default_rng(11))
    fit = P.ols_hc3(design, y)
    for s in range(y.shape[1]):
        reference = sm.OLS(y[:, s], design.to_numpy()).fit(cov_type="HC3")
        np.testing.assert_allclose(fit.beta[:, s], reference.params, atol=1e-10, rtol=0)
        np.testing.assert_allclose(fit.cov_hc3[s], reference.cov_params(), atol=1e-10, rtol=0)
    assert np.all(np.abs(fit.cov_hc3[:, 4, 5]) < 1e-12)
    assert P.arm_residual_df(design.to_numpy(), meta["arm"].eq("ISS-T").to_numpy()) == 17


@pytest.mark.parametrize("kappa", [1.0, 0.5])
@pytest.mark.parametrize("adjusted", [False, True])
def test_freedman_lane_identity_is_hc3_t_and_delta_is_reparameterized_contrast(kappa, adjusted):
    rng = np.random.default_rng(5)
    meta, design, y = _joint_synthetic(rng)
    blocks = design[list(P.BLOCK_COLUMNS)].to_numpy()
    f_t, f_l = design["F_T"].to_numpy(), design["F_L"].to_numpy()
    nuisance = blocks
    x = design.to_numpy()
    if adjusted:
        cov = rng.normal(size=40)
        iss = meta["arm"].eq("ISS-T").to_numpy(float)
        extra = np.column_stack([cov * iss, cov * (1 - iss)])
        nuisance = np.column_stack([blocks, extra])
        x = np.column_stack([x, extra])
    fit = P.ols_hc3(x, y)
    delta = fit.beta[5] - kappa * fit.beta[4]
    se = np.sqrt(fit.cov_hc3[:, 5, 5] + kappa**2 * fit.cov_hc3[:, 4, 4] - 2 * kappa * fit.cov_hc3[:, 4, 5])
    result = P.freedman_lane_contrast(y, nuisance, f_t, f_l, kappa, meta["block"], 99, seed=1)
    np.testing.assert_allclose(result.estimate, delta, atol=1e-12, rtol=0)
    np.testing.assert_allclose(result.t_obs, delta / se, atol=1e-10, rtol=0)
    np.testing.assert_allclose(result.standard_error, se, atol=1e-12, rtol=0)
    assert result.t_null.shape == (99, y.shape[1]) and np.all(np.isfinite(result.t_null))


def test_permutations_stay_within_blocks_and_are_shared_across_sets():
    meta = P.synthetic_meta40(5)
    positions = P.block_positions(meta["block"])
    perms = P.within_block_permutations(positions, 40, 500, np.random.default_rng(2))
    blocks = meta["block"].to_numpy()
    assert np.all(blocks[perms] == blocks[None, :])
    assert np.all(np.sort(perms, axis=1) == np.arange(40)[None, :])
    _, design, y = _joint_synthetic(np.random.default_rng(9), n_sets=3)
    blocks_x = design[list(P.BLOCK_COLUMNS)].to_numpy()
    args = (blocks_x, design["F_T"].to_numpy(), design["F_L"].to_numpy(), 1.0, meta["block"], 300, 44, 64)
    joint = P.freedman_lane_contrast(y, *args)
    alone = P.freedman_lane_contrast(y[:, [1]], *args)
    np.testing.assert_allclose(joint.t_null[:, 1], alone.t_null[:, 0], atol=1e-12)


def test_freedman_lane_max_t_matches_helper_and_bounds():
    meta, design, y = _joint_synthetic(np.random.default_rng(3), n_sets=4)
    blocks = design[list(P.BLOCK_COLUMNS)].to_numpy()
    res = P.freedman_lane(y, np.column_stack([blocks, design["F_T"] + design["F_L"]]), design["F_L"], meta["block"], 400, 8, family=[0, 1, 2, 3])
    np.testing.assert_allclose(res.p_max_t, P.max_t_pvalues(res.t_obs, res.t_null, [0, 1, 2, 3]))
    assert np.all(res.p_max_t >= res.p_two_sided - 1e-15)
    assert np.all((res.p_upper + res.p_lower) >= 1.0)  # ties on the observed value are counted twice


def _calibration_rate(theta_t, theta_l, kappa, side, n_sim=300, n_perm=199, seed=20260811):
    rng = np.random.default_rng(seed)
    meta = P.synthetic_meta40(5)
    design = P.design_matrix(meta)
    blocks = design[list(P.BLOCK_COLUMNS)].to_numpy()
    f_t, f_l = design["F_T"].to_numpy(), design["F_L"].to_numpy()
    mean = blocks @ np.array([0.2, -0.4, 0.1, 0.3]) + theta_t * f_t + theta_l * f_l
    seeds = np.random.SeedSequence(seed).spawn(n_sim)
    rejections = 0
    for child in seeds:
        y = mean + rng.standard_normal(40)
        res = P.freedman_lane_contrast(y, blocks, f_t, f_l, kappa, meta["block"], n_perm, child, keep_null=False)
        p = res.p_two_sided[0] if side == "two" else res.p_upper[0]
        rejections += p <= 0.05
    return rejections / n_sim


def test_calibration_interaction_two_sided_under_equal_effects():
    rate = _calibration_rate(1.0, 1.0, 1.0, "two")
    assert 0.02 <= rate <= 0.09, rate


def test_calibration_margin_one_sided_at_exact_half_persistence():
    rate = _calibration_rate(1.0, 0.5, 0.5, "upper", seed=20260812)
    assert 0.02 <= rate <= 0.09, rate


# --------------------------------------------------------------------------- ratio summaries


def test_fieller_hand_case_and_unbounded_case():
    tT, vT, tL, vL, df, level = 2.0, 0.04, 1.0, 0.04, 17, 0.95
    q = stats.t.ppf(0.975, df)
    a, b, c = tT**2 - q**2 * vT, -2 * tT * tL, tL**2 - q**2 * vL
    roots = sorted(((-b - math.sqrt(b * b - 4 * a * c)) / (2 * a), (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)))
    low, high, bounded = P.fieller_ratio(tT, vT, tL, vL, df, level)
    assert bounded and abs(low - roots[0]) < 1e-12 and abs(high - roots[1]) < 1e-12
    assert low < tL / tT < high
    for rho in (low, high):  # endpoints sit exactly on the Fieller boundary
        assert abs((tL - rho * tT) ** 2 - q**2 * (vL + rho**2 * vT)) < 1e-10
    # hand numbers with q = 1.96 (df -> infinity): (4 -/+ sqrt(16 - 4ac)) / 2a
    low_inf, high_inf, _ = P.fieller_ratio(2.0, 0.04, 1.0, 0.04, np.inf, 0.95)
    assert abs(low_inf - 0.2956) < 1e-3 and abs(high_inf - 0.7443) < 1e-3
    # theta_T not significant (tT^2 <= q^2 vT): no bounded interval
    assert P.fieller_ratio(0.3, 0.2, 0.5, 0.04, 17, 0.95)[2] is False
    assert all(np.isnan(P.fieller_ratio(0.3, 0.2, 0.5, 0.04, 17, 0.95)[:2]))


def test_bootstrap_is_exact_when_cells_are_constant_and_posterior_matches_normal_cdf():
    meta = P.synthetic_meta40(5)
    design = P.design_matrix(meta)
    y = 2.0 * design["F_T"] + 1.0 * design["F_L"] + design[list(P.BLOCK_COLUMNS)].to_numpy() @ np.arange(4.0)
    boot = P.stratified_bootstrap_ratio(pd.DataFrame({"s": y}), meta, B=50, seed=1)
    assert abs(boot.loc["s", "rho_boot95_low"] - 0.5) < 1e-12 and abs(boot.loc["s", "rho_boot95_high"] - 0.5) < 1e-12
    assert boot.loc["s", "boot_frac_thetaT_le0"] == 0 and not boot.loc["s", "bootstrap_unreliable"]
    post = P.posterior_summaries([1.0], [0.01], [0.3], [0.04], draws=200_000, seed=4)
    assert abs(post.loc[0, "post_p_thetaL_le_0"] - stats.norm.cdf(-0.3 / 0.2)) < 0.005
    expected = stats.norm.sf(0, loc=0.3 - 0.5, scale=math.sqrt(0.04 + 0.25 * 0.01))  # theta_T > 0 almost surely
    assert abs(post.loc[0, "post_p_rho_ge_0_5"] - expected) < 0.005


# --------------------------------------------------------------------------- classification


def _row(**kw):
    base = dict(
        technical_flag=0.0, theta_isst=1.2, se_isst=0.4, df_isst=17, theta_lar=1.0, se_lar=0.4, df_lar=17,
        margin=0.4, margin_se=0.2, margin_df=30.0, p_persist=0.01, p_decay=0.99, pooled_terminal_estimate=0.7,
    )
    base.update(kw)
    return base


Q90 = stats.t.ppf(0.95, 30.0)


@pytest.mark.parametrize(
    ("row", "label", "reason_part"),
    [
        (_row(), P.PERSISTS, "FL p(margin>0)"),
        (_row(technical_flag=1.0), P.NOT_EVALUABLE_TECHNICAL, "flagged"),
        (_row(technical_flag=1.0, theta_isst=0.0), P.NOT_EVALUABLE_TECHNICAL, "flagged"),
        (_row(theta_isst=0.5, se_isst=0.4), P.NOT_EVALUABLE_NO_TERMINAL_EFFECT, "includes 0"),
        (_row(theta_isst=-1.2), P.NOT_EVALUABLE_NO_TERMINAL_EFFECT, "sign differs"),
        (_row(theta_isst=1.2, pooled_terminal_estimate=-0.7), P.NOT_EVALUABLE_NO_TERMINAL_EFFECT, "sign differs"),
        (_row(margin=-0.5, p_persist=0.99, p_decay=0.01, theta_lar=0.1), P.DECAYS, "FL p(margin<0)"),
        (_row(margin=-0.9, p_persist=0.999, p_decay=0.001, theta_lar=-1.5, se_lar=0.3), P.DECAYS_REVERSED, "REVERSED"),
        # FL says persists, HC3 interval does not: disagreement -> INCONCLUSIVE
        (_row(margin=0.9 * Q90 * 0.2, p_persist=0.03), P.INCONCLUSIVE, "FL/HC3 disagreement"),
        # HC3 interval says persists, FL does not
        (_row(margin=0.4, p_persist=0.08), P.INCONCLUSIVE, "FL/HC3 disagreement"),
        (_row(margin=-0.4, p_persist=0.99, p_decay=0.07), P.INCONCLUSIVE, "FL/HC3 disagreement"),
        (_row(margin=0.05, p_persist=0.4, p_decay=0.6), P.INCONCLUSIVE, "neither margin test"),
        (_row(margin=np.nan, p_persist=np.nan, p_decay=np.nan), P.INCONCLUSIVE, "not available"),
        # boundary: FL exactly at alpha passes (<=); HC3 bound must be strictly > 0
        (_row(margin=1.0001 * Q90 * 0.2, p_persist=0.05), P.PERSISTS, "FL p"),
    ],
)
def test_classification_is_table_driven_across_all_labels(row, label, reason_part):
    got, reason = P.classify_persistence(row, _rules())
    assert got == label
    assert reason_part in reason


def test_rules_are_read_from_prereg_yaml_and_strict(tmp_path):
    path = tmp_path / "prereg.yaml"
    path.write_text(PREREG_4_2)
    rules = P.parse_classification_rules(R.load_prereg(path), source=str(path))
    assert (rules.terminal_ci_level, rules.persist_ci_level, rules.decay_ci_level) == (0.95, 0.90, 0.90)
    assert (rules.persist_alpha, rules.decay_alpha, rules.reversed_ci_level) == (0.05, 0.05, 0.95)
    assert (rules.margin_kappa, rules.interaction_kappa) == (0.5, 1.0)
    assert rules.require_terminal_sign_match
    # a changed YAML changes the classifier
    stricter = yaml.safe_load(PREREG_4_2)
    stricter["classification"]["PERSISTS"] = "FL one-sided p(margin>0) <= 0.01 AND HC3 95% CI lower(margin) > 0"
    tight = P.parse_classification_rules(stricter)
    assert P.classify_persistence(_row(p_persist=0.03), tight)[0] == P.INCONCLUSIVE
    # a structured mapping is accepted too
    stricter["classification"]["PERSISTS"] = {"fl_one_sided_alpha": 0.05, "hc3_ci_level": 0.9}
    assert P.parse_classification_rules(stricter).persist_alpha == 0.05
    broken = yaml.safe_load(PREREG_4_2)
    broken["classification"]["PERSISTS"] = "persistence is likely"
    with pytest.raises(P.RuleParseError):
        P.parse_classification_rules(broken)
    missing = yaml.safe_load(PREREG_4_2)
    del missing["classification"]["DECAYS"]
    with pytest.raises(P.RuleParseError):
        P.parse_classification_rules(missing)
    settings = R.prereg_settings(R.load_prereg(path), CONFIG_SYNTH)
    assert (settings.n_perm, settings.n_boot, settings.min_genes) == (100_000, 10_000, 8)
    assert settings.seed_for("freedman_lane") == 20260811 + 700


def test_vectorized_and_row_classifier_agree():
    rng = np.random.default_rng(0)
    rows = [
        _row(
            technical_flag=float(rng.random() < 0.1), theta_isst=rng.normal(1, 0.8), theta_lar=rng.normal(0.5, 1),
            margin=rng.normal(0, 0.4), p_persist=rng.random(), p_decay=rng.random(), margin_df=rng.uniform(10, 34),
        )
        for _ in range(200)
    ]
    rules = _rules()
    vector = P.classify_persistence_arrays({k: np.array([r[k] for r in rows]) for k in rows[0]}, rules)["label"]
    assert list(vector) == [P.classify_persistence(r, rules)[0] for r in rows]


def test_conservative_verdict_across_maps():
    assert P.conservative_verdict([P.PERSISTS, P.PERSISTS]) == (P.PERSISTS, True)
    assert P.conservative_verdict([P.PERSISTS, P.INCONCLUSIVE]) == (P.INCONCLUSIVE, False)
    assert P.conservative_verdict([P.PERSISTS, P.DECAYS]) == (P.INCONCLUSIVE, False)
    assert P.conservative_verdict([P.DECAYS, P.DECAYS_REVERSED]) == (P.DECAYS, False)
    assert P.conservative_verdict([P.INCONCLUSIVE, P.NOT_EVALUABLE_NO_TERMINAL_EFFECT]) == (P.NOT_EVALUABLE_NO_TERMINAL_EFFECT, False)


# --------------------------------------------------------------------------- pattern concordance


def test_pattern_concordance_copy_case_gives_r_one_and_p_one_over_n():
    cohorts = synthetic_cohorts()
    family = {name: P.single_set_spec(genes, role="t") for name, genes in SETS.items() if len(genes) >= 8}
    scores, _ = P.score_arm(cohorts["iss"], family, 0.1, add_difference=False)
    meta = cohorts["iss"].metadata
    reference = P.observed_arm_effects(scores, meta)["estimate_g"]
    result = P.pattern_concordance_exact(reference, scores.copy(), meta.copy())
    assert abs(result["pearson_r"] - 1.0) < 1e-12
    assert result["n_allocations"] == 63_504
    assert result["p_one_sided_exact"] == pytest.approx(1.0 / 63_504)


# --------------------------------------------------------------------------- power


def test_power_simulation_is_deterministic_given_seed():
    kwargs = dict(n_per_cell=5, n_sim=60, n_perm=49)
    first = P.simulate_power([(1.174, 1.0), (0.689, 0.0)], seed=123, **kwargs)
    second = P.simulate_power([(1.174, 1.0), (0.689, 0.0)], seed=123, **kwargs)
    pd.testing.assert_frame_equal(first, second)
    other = P.simulate_power([(1.174, 1.0), (0.689, 0.0)], seed=124, **kwargs)
    assert not first["power_interaction_hc3"].equals(other["power_interaction_hc3"])
    label_cols = [c for c in first.columns if c.startswith("p_label_") and c.split("p_label_")[1] in P.LABELS]
    np.testing.assert_allclose(first[label_cols].sum(axis=1), 1.0)
    analytic = P.analytic_power(0.689, 0.0)
    assert abs(analytic["analytic_se_theta"] - 0.447) < 1e-3 and abs(analytic["analytic_se_margin"] - 0.5) < 1e-3
    assert abs(analytic["analytic_power_interaction"] - 0.19) < 0.01
    assert abs(analytic["analytic_power_margin_decay"] - 0.17) < 0.01


def test_power_only_mode_loads_no_expression_data(tmp_path, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("power-only mode must not load expression data")

    for name in ("load_persistence_cohorts", "load_primary_missions", "load_moderator_missions"):
        monkeypatch.setattr(P, name, forbidden)
    monkeypatch.setattr(R, "POWER_THETAS", (1.174,))
    monkeypatch.setattr(R, "POWER_RHOS", (1.0, 0.0))
    code = R.main(
        ["--power-only", "--prereg", str(tmp_path / "absent.yaml"), "--results", str(tmp_path),
         "--power-sims", "40", "--power-permutations", "19"]
    )
    assert code == 0
    table = pd.read_csv(tmp_path / R.OUTPUT_SUBDIR / "persistence_power.tsv", sep="\t")
    assert list(table["rho"]) == [1.0, 0.0]
    manifest = (tmp_path / R.OUTPUT_SUBDIR / "persistence_power_manifest.json").read_text()
    assert "builtin_plan_section_4.2_draft" in manifest


# --------------------------------------------------------------------------- preregistration guard / CLI


def _fake_git(status_out="", status_code=0, tracked_code=0):
    def runner(cmd, **_kwargs):
        sub = cmd[1]
        if sub == "status":
            return subprocess.CompletedProcess(cmd, status_code, status_out, "")
        if sub == "ls-files":
            return subprocess.CompletedProcess(cmd, tracked_code, "", "")
        return subprocess.CompletedProcess(cmd, 0, "abc123\n", "")

    return runner


def test_prereg_guard_refuses_missing_dirty_untracked_and_records_clean(tmp_path):
    path = tmp_path / "prereg.yaml"
    with pytest.raises(R.PreregistrationError, match="missing"):
        R.prereg_guard(path)
    path.write_text(PREREG_4_2)
    with pytest.raises(R.PreregistrationError, match="uncommitted"):
        R.prereg_guard(path, runner=_fake_git(status_out=" M prereg.yaml"))
    with pytest.raises(R.PreregistrationError, match="not tracked"):
        R.prereg_guard(path, runner=_fake_git(tracked_code=1))
    with pytest.raises(R.PreregistrationError):  # a real git call outside any repository
        R.prereg_guard(path)
    record = R.prereg_guard(path, runner=_fake_git())
    assert record["prereg_guard"] == "clean" and record["git_head"] == "abc123"
    assert len(record["prereg_sha256"]) == 64
    bypass = R.prereg_guard(path, allow_uncommitted=True, runner=_fake_git(status_out="?? prereg.yaml"))
    assert bypass["prereg_guard"].startswith("bypassed: uncommitted")


def test_full_run_refuses_without_committed_prereg(tmp_path):
    with pytest.raises(R.PreregistrationError):
        R.main(["--prereg", str(tmp_path / "absent.yaml"), "--results", str(tmp_path)])
    assert not (tmp_path / R.OUTPUT_SUBDIR).exists()


def test_gene_map_override_is_in_memory_only(tmp_path):
    config = yaml.safe_load(CONFIG.read_text())
    original = config["gene_mapping"]["path"]
    alt = tmp_path / "alt_map.tsv"
    alt.write_text("ensembl_gene_id\tmgi_symbol\n")
    updated = R.apply_gene_map_override(config, alt)
    assert updated["gene_mapping"]["path"] == str(alt.resolve())
    assert config["gene_mapping"]["path"] == original
    assert R.apply_gene_map_override(config, None) == config
    with pytest.raises(FileNotFoundError):
        R.apply_gene_map_override(config, tmp_path / "nope.tsv")


# --------------------------------------------------------------------------- synthetic end-to-end


def test_synthetic_end_to_end_tables_and_internal_contracts(synthetic_run, tmp_path):
    result = synthetic_run["result"]
    defs = synthetic_run["defs"]
    cohorts = synthetic_run["cohorts"]
    contrasts = result["contrasts"].set_index("set")
    # targets present, ordered first, scopes as preregistered
    assert list(result["contrasts"]["set"][:3]) == [P.PODOCYTE_SET, P.STRUCTURAL_DISJOINT, P.D_CONTRAST]
    assert contrasts.loc[P.PODOCYTE_SET, "family"] == "primary"
    assert contrasts.loc[P.D_CONTRAST, "family"] == "key_secondary"
    assert contrasts.loc["axis_barrier", "family"] == "context"
    # the LAR-dead small set is evaluable in ISS-T and the primaries but not in Layer 2 / LAR
    assert "small__high_specificity" in defs.compartment and "small__high_specificity" not in contrasts.index
    assert "tiny__undefined" not in defs.compartment
    coverage = result["gene_coverage"].set_index("set")
    assert not coverage.loc["small__high_specificity", "evaluable_lar"]
    assert coverage.loc["small__high_specificity", "evaluable_isst"]
    # every spec'd contrast column is written
    spec_columns = [
        "set", "family", "theta_isst", "se_isst", "theta_lar", "se_lar", "g_isst_joint", "g_lar_joint", "interaction",
        "interaction_se", "interaction_df", "interaction_ci95_low", "interaction_ci95_high", "p_interaction_hc3",
        "p_interaction_fl", "interaction_maxT_fwer_fl", "interaction_bh_q", "margin", "margin_se", "margin_df",
        "margin_ci90_low", "margin_ci90_high", "p_persist_fl", "p_decay_fl", "rho_hat", "rho_fieller90_low",
        "rho_fieller90_high", "rho_fieller95_low", "rho_fieller95_high", "rho_boot95_low", "rho_boot95_high",
        "boot_frac_thetaT_le0", "post_p_rho_ge_0_5", "post_p_thetaL_le_0", "isst_effect_established",
        "classification", "classification_reason",
    ]
    assert set(spec_columns) <= set(result["contrasts"].columns)
    assert set(contrasts["classification_base"]) <= set(P.LABELS)
    # D is the difference of the joint scores; orientation follows the pooled terminal sign
    assert contrasts.loc["tubule__high_specificity", "orientation"] == -1.0
    assert contrasts.loc[P.PODOCYTE_SET, "orientation"] == 1.0
    # Satterthwaite with equal per-arm df and HC3 per-arm df 17
    assert set(contrasts["df_isst"]) == {17} and set(contrasts["df_lar"]) == {17}
    # FL observed t equals the HC3 t of each contrast
    np.testing.assert_allclose(contrasts["margin_t_fl_identity"], contrasts["margin_t"], atol=1e-10)
    # key secondary Holm and max-T bounds
    key = contrasts.loc[[P.STRUCTURAL_DISJOINT, P.D_CONTRAST]]
    assert np.all(key["p_persist_fl_holm_key_secondary"] >= key["p_persist_fl"] - 1e-15)
    assert np.all(contrasts["interaction_maxT_fwer_fl"].dropna() >= contrasts.loc[contrasts["interaction_maxT_fwer_fl"].notna(), "p_interaction_fl"] - 1e-15)
    # Layer 1: P1 rows carry coef_score_units (never g) and FL identity
    arm = result["arm_effects"]
    p1 = arm[arm["set"] == P.P1_ADJUSTED]
    assert len(p1) == 3 and p1["estimate_g"].isna().all() and p1["fl_t_identity_matches_hc3"].all()
    assert set(arm.loc[arm["cohort"] == "OSD-771", "n_allocations"].dropna()) == {63_504}
    assert set(arm.loc[arm["cohort"] == "OSD-513", "n_allocations"].dropna()) == {48_620}
    # Layer-1 ISS-T podocyte g equals the compartment-context code path on the same inputs
    family = {P.PODOCYTE_SET: defs.compartment[P.PODOCYTE_SET]}
    from src.clinical_axes.analysis import combined_score_design

    scores, design, _, _ = combined_score_design(cohorts["primary"], family, cpm_threshold=0.1)
    path_effect = blocked_meta_permutation(scores, design, n_permutations=1, seed=0).mission_effects
    reference = path_effect[path_effect["mission"] == "OSD-771"].iloc[0]
    ours = arm[(arm["cohort"] == "OSD-771") & (arm["set"] == P.PODOCYTE_SET)].iloc[0]
    assert abs(ours["estimate_g"] - reference["estimate"]) < 1e-9
    # Layer 3 and 4 and context are produced
    assert result["adjusted"]["set"].tolist() == [P.P1_ADJUSTED]
    assert set(result["adjusted"]["df_isst"]) == {16}
    assert set(result["concordance"]["comparison"]) == {"ISS-T_per_arm_g_vs_LAR_per_arm_g", "OSD-513_g_vs_pooled_5_mission_terminal"}
    assert {"OSD-253_duration", "OSD-771_joint_age_stratum"} <= set(result["duration_context"]["context"])
    # writing every output file
    out = tmp_path / "out"
    R.write_outputs(out, result, {"mode": "test"})
    for name in (
        "persistence_cohort_manifest.tsv", "persistence_technical_qc.tsv", "persistence_gene_coverage.tsv",
        "persistence_arm_effects.tsv", "persistence_contrasts.tsv", "persistence_adjusted.tsv",
        "persistence_pattern_concordance.tsv", "persistence_duration_context.tsv", "persistence_scores.tsv",
        "persistence_null_t.tsv.gz", "persistence_manifest.json",
    ):
        assert (out / name).exists(), name
    null = pd.read_csv(out / "persistence_null_t.tsv.gz", sep="\t")
    assert len(null) == synthetic_run["settings"].n_perm and f"margin::{P.PODOCYTE_SET}" in null


def test_headline_equals_full_run_rows_and_idmap_sensitivity(synthetic_run, tmp_path):
    cohorts, defs, settings = synthetic_run["cohorts"], synthetic_run["defs"], synthetic_run["settings"]
    head = R.headline(cohorts, defs, CONFIG_SYNTH, _rules(), settings).set_index("set")
    full = synthetic_run["result"]["contrasts"].set_index("set")
    targets = [P.PODOCYTE_SET, P.STRUCTURAL_DISJOINT, P.D_CONTRAST]
    assert list(head.index) == targets
    for column in ("theta_isst", "margin", "p_persist_rule", "p_decay_rule", "interaction_maxT_fwer_fl", "rho_boot95_low", "post_p_rho_ge_0_5"):
        np.testing.assert_allclose(head.loc[targets, column].astype(float), full.loc[targets, column].astype(float), atol=1e-12)
    assert list(head["classification"]) == list(full.loc[head.index, "classification"])
    # the idmap mode with a fake loader: one map per synthetic world
    prereg = tmp_path / "prereg.yaml"
    prereg.write_text(PREREG_4_2)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({**CONFIG_SYNTH, "gene_mapping": {"path": "x", "annotation_fallback": "y"}}))
    maps = [tmp_path / "map_a.tsv", tmp_path / "map_b.tsv"]
    for path in maps:
        path.write_text("ensembl_gene_id\tmgi_symbol\n")
    worlds = {str(maps[0].resolve()): synthetic_cohorts(), str(maps[1].resolve()): synthetic_cohorts(seed=8, lar_effect=1.2)}
    args = R.build_parser().parse_args(
        ["--config", str(config_path), "--prereg", str(prereg), "--tiers", str(synthetic_run["tiers"]),
         "--results", str(tmp_path), "--allow-uncommitted-prereg", "--permutations", "99", "--bootstrap", "50",
         "--posterior-draws", "1000", "--idmap-sensitivity", str(maps[0]), str(maps[1])]
    )
    table = R.run_idmap_sensitivity(args, loader=lambda config, _root: worlds[config["gene_mapping"]["path"]])
    robust = table[table["output"] == f"{P.PODOCYTE_SET}::robust_to_gene_id_reconstruction"]
    assert len(robust) == 1
    labels = table[table["output"] == f"{P.PODOCYTE_SET}::classification"]["value"].tolist()
    assert bool(robust["value"].iloc[0]) == (len(set(labels)) == 1)
    assert (tmp_path / R.OUTPUT_SUBDIR / "idmap_sensitivity.tsv").exists()


def test_lar_technical_flag_makes_every_layer2_label_not_evaluable(synthetic_run):
    cohorts = synthetic_cohorts(lar_qc_shift=0.5)
    settings = synthetic_settings()
    result = R.analyze(cohorts, synthetic_run["defs"], CONFIG_SYNTH, _rules(), settings, log=lambda _m: None)
    assert result["lar_technical_flag"]
    assert set(result["contrasts"]["classification"]) == {P.NOT_EVALUABLE_TECHNICAL}
    assert np.isfinite(result["contrasts"]["margin"]).all()  # still written as descriptive


def test_full_cli_run_on_synthetic_world_with_bypassed_guard(synthetic_run, tmp_path):
    prereg = tmp_path / "prereg.yaml"
    prereg.write_text(PREREG_4_2)
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.safe_dump({**CONFIG_SYNTH, "gene_mapping": {"path": "x", "annotation_fallback": "y"}}))
    args = R.build_parser().parse_args(
        ["--config", str(config_path), "--prereg", str(prereg), "--tiers", str(synthetic_run["tiers"]),
         "--results", str(tmp_path), "--allow-uncommitted-prereg", "--permutations", "99", "--bootstrap", "50",
         "--posterior-draws", "1000"]
    )
    R.run(args, loader=lambda _config, _root: synthetic_run["cohorts"])
    manifest = yaml.safe_load((tmp_path / R.OUTPUT_SUBDIR / "persistence_manifest.json").read_text())
    assert manifest["prereg_guard"].startswith("bypassed")
    assert manifest["prereg_sha256"] == R.sha256(prereg)
    assert any("Freedman-Lane permutations 99" in d for d in manifest["deviations"])


# --------------------------------------------------------------------------- real-data contracts (known quantities only)


@pytest.fixture(scope="module")
def real_config():
    return yaml.safe_load(CONFIG.read_text())


@requires_paths(CONFIG, ID_MAP, OSD771_VST, OSD513_VST)
def test_real_cohorts_have_the_frozen_persistence_design(real_config):
    cohorts = P.load_persistence_cohorts(real_config, REPO)  # loading only; nothing is scored
    assert set(cohorts["primary"]) == {"OSD-102", "OSD-163", "OSD-253", "OSD-462", "OSD-771"}
    _, _, meta40 = P.joint_frame(cohorts["iss"], cohorts["lar"])
    assert meta40.groupby(["block", "condition"]).size().eq(5).all() and len(meta40) == 40


@requires_paths(CONFIG, ID_MAP, OSD771_VST, TIERS)
def test_real_isst_podocyte_g_equals_compartment_context_mission_effect(real_config):
    """ISS-T is the primary terminal arm (already known). No LAR/OSD-513 data is loaded."""
    from scripts.clinical_axes.run_compartment_context import load_family
    from src.clinical_axes.analysis import combined_score_design
    from src.clinical_axes.data import load_primary_missions

    threshold = float(real_config["eligibility"]["cpm_threshold"])
    primary = load_primary_missions(real_config, REPO)
    family, _ = load_family(TIERS)
    podocyte = {P.PODOCYTE_SET: family[P.PODOCYTE_SET]}
    scores, design, _, _ = combined_score_design(primary, podocyte, cpm_threshold=threshold)
    effects = blocked_meta_permutation(scores, design, n_permutations=1, seed=0).mission_effects
    reference = effects[effects["mission"] == "OSD-771"].iloc[0]
    arm_scores, _ = P.score_arm(primary["OSD-771"], podocyte, threshold)
    ours = P.arm_exact_test(arm_scores, primary["OSD-771"].metadata).set_index("set").loc[P.PODOCYTE_SET]
    assert abs(ours["estimate_g"] - reference["estimate"]) < 1e-9
    assert abs(ours["variance"] - reference["variance"]) < 1e-9
    assert int(ours["n_allocations"]) == 63_504


@requires_paths(CONFIG, ID_MAP, OSD771_VST, OSD513_VST)
def test_real_recovery_four_axis_effects_reproduce_moderator_code_path(real_config):
    """Four clinical axes only (known before lock); compartment sets are never scored here."""
    from scripts.clinical_axes.run_cross_mission import _moderator_effects
    from src.clinical_axes.data import load_moderator_missions

    threshold = float(real_config["eligibility"]["cpm_threshold"])
    axes = real_config["primary_family"]
    moderators = load_moderator_missions(real_config, REPO)
    reference, _, _ = _moderator_effects(moderators, axes, threshold)
    reference = reference.set_index(["mission", "axis"])
    for mission, data in moderators.items():
        scores, _ = P.score_arm(data, axes, threshold, add_difference=False)
        assert list(scores.columns) == list(axes)
        observed = P.observed_arm_effects(scores, data.metadata)
        for axis in axes:
            assert abs(observed.loc[axis, "estimate_g"] - reference.loc[(mission, axis), "estimate"]) < 1e-9
            assert abs(observed.loc[axis, "variance"] - reference.loc[(mission, axis), "variance"]) < 1e-9


# --------------------------------------------------------------------------- margin sweep


def test_margin_sweep_kappas_are_parsed_and_must_include_the_primary():
    prereg = yaml.safe_load(PREREG_4_2)
    assert P.parse_margin_sweep(prereg, 0.5) == (0.0, 0.25, 0.5, 0.75, 1.0)
    assert P.parse_margin_sweep({}, 0.5) == (0.5,)
    with pytest.raises(P.RuleParseError):
        P.parse_margin_sweep({"margin_sweep": {"kappas": [0.0, 1.0]}}, 0.5)
    with pytest.raises(P.RuleParseError):
        P.parse_margin_sweep({"margin_sweep": {"kappas": [0.5, 0.5]}}, 0.5)
    assert synthetic_settings().margin_kappas == (0.0, 0.25, 0.5, 0.75, 1.0)


def test_margin_sweep_endpoints_and_primary_agree_with_main_contrasts(synthetic_run):
    result = synthetic_run["result"]
    sweep = result["margin_sweep"]
    contrasts = result["contrasts"].set_index("set")
    assert list(sweep.columns[:12]) == R.SWEEP_COLUMNS
    assert set(sweep["set"]) == {P.PODOCYTE_SET, P.STRUCTURAL_DISJOINT, P.D_CONTRAST}
    by = {k: part.set_index("set") for k, part in sweep.groupby("kappa")}
    targets = list(by[0.5].index)
    # kappa = 1 is the interaction (same estimate, SE, df and the same FL permutations)
    np.testing.assert_allclose(by[1.0]["margin"], contrasts.loc[targets, "interaction"], atol=1e-12)
    np.testing.assert_allclose(by[1.0]["se"], contrasts.loc[targets, "interaction_se"], atol=1e-12)
    np.testing.assert_allclose(by[1.0]["p_persist_fl"], contrasts.loc[targets, "p_interaction_fl_upper"], atol=0)
    # kappa = 0 is theta_LAR ("any retention")
    np.testing.assert_allclose(by[0.0]["margin"], contrasts.loc[targets, "theta_lar"], atol=1e-12)
    np.testing.assert_allclose(by[0.0]["se"], contrasts.loc[targets, "se_lar"], atol=1e-12)
    # kappa = 0.5 reproduces the primary margin row and its classification exactly
    np.testing.assert_allclose(by[0.5]["margin"], contrasts.loc[targets, "margin"], atol=1e-12)
    np.testing.assert_allclose(by[0.5]["p_persist_fl"], contrasts.loc[targets, "p_persist_fl"], atol=0)
    np.testing.assert_allclose(by[0.5]["p_decay_fl"], contrasts.loc[targets, "p_decay_fl"], atol=0)
    assert list(by[0.5]["label_descriptive"]) == list(contrasts.loc[targets, "classification"])
    assert by[0.5]["is_primary"].all() and not by[0.25]["is_primary"].any()
    # margins are linear in kappa and the persist p-value cannot fall as kappa rises
    for name in targets:
        series = sweep[sweep["set"] == name].sort_values("kappa")
        expected = series["theta_lar"] - series["kappa"] * series["theta_isst"]
        np.testing.assert_allclose(series["margin"], expected, atol=1e-12)


def test_rho_lower90_fieller_inverts_the_hc3_margin_test():
    rng = np.random.default_rng(21)
    meta = P.synthetic_meta40(5)
    design = P.design_matrix(meta)
    y = 1.5 * design["F_T"] + 1.1 * design["F_L"] + rng.standard_normal(40) * 0.6
    fit = P.ols_hc3(design, y.to_numpy())
    tT, tL = fit.beta[4, 0], fit.beta[5, 0]
    vT, vL = fit.cov_hc3[0, 4, 4], fit.cov_hc3[0, 5, 5]
    low, _, bounded = P.fieller_ratio(tT, vT, tL, vL, 17, 0.90)
    assert bounded
    q = stats.t.ppf(0.95, 17)
    t_at = lambda k: (tL - k * tT) / math.sqrt(vL + k * k * vT)  # noqa: E731
    assert abs(t_at(low) - q) < 1e-9  # the bound is where the one-sided 5% margin test flips
    assert t_at(low - 1e-3) > q and t_at(low + 1e-3) < q


def test_margin_sweep_primary_kappa_and_targets_are_validated():
    prereg = yaml.safe_load(PREREG_4_2)
    prereg["margin_sweep"] = {"kappas": [0.0, 0.5, 1.0], "primary_kappa": 0.25}
    with pytest.raises(P.RuleParseError):
        P.parse_margin_sweep(prereg, 0.5)
    prereg["margin_sweep"] = {"kappas": [0.0, 0.5, 1.0], "primary_kappa": 0.5, "applies_to": [P.PODOCYTE_SET]}
    assert P.parse_margin_sweep(prereg, 0.5) == (0.0, 0.5, 1.0)
    assert P.margin_sweep_targets(prereg) == (P.PODOCYTE_SET,)
    prereg["margin_sweep"]["applies_to"] = ["tubule__high_specificity"]
    with pytest.raises(P.RuleParseError):
        P.margin_sweep_targets(prereg)


def test_validate_prereg_only_checks_gene_map_hashes_and_loads_no_data(tmp_path, monkeypatch):
    def forbidden(*_args, **_kwargs):
        raise AssertionError("validate-only mode must not load expression data")

    for name in ("load_persistence_cohorts", "load_primary_missions", "load_moderator_missions"):
        monkeypatch.setattr(P, name, forbidden)
    baseline, runner = tmp_path / "base.tsv", tmp_path / "runner.tsv"
    baseline.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG1\tA\n")
    runner.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG1\tB\n")
    prereg = yaml.safe_load(PREREG_4_2)
    prereg["gene_map"] = {"baseline": str(baseline), "baseline_sha256": R.sha256(baseline),
                          "runner_up": str(runner), "runner_up_sha256": R.sha256(runner)}
    path = tmp_path / "prereg.yaml"
    path.write_text(yaml.safe_dump(prereg))
    config = tmp_path / "config.yaml"
    config.write_text(yaml.safe_dump({**CONFIG_SYNTH, "gene_mapping": {"path": str(baseline), "annotation_fallback": "y"}}))
    args = ["--validate-prereg-only", "--prereg", str(path), "--config", str(config), "--results", str(tmp_path)]
    with pytest.raises(R.PreregistrationError):  # untracked file, no bypass
        R.main(args)
    summary = R.run_validate_prereg(R.build_parser().parse_args(args + ["--allow-uncommitted-prereg"]))
    assert summary["analysis_gene_map"]["role"] == "baseline" and summary["analysis_gene_map"]["sha256_matches"]
    assert summary["declared_gene_maps"]["runner_up"]["sha256_matches"] and summary["deviations"] == []
    assert summary["margin_sweep"]["kappas"] == [0.0, 0.25, 0.5, 0.75, 1.0]
    runner_record = R.gene_map_record(runner, prereg)
    assert runner_record["role"] == "runner_up"
    prereg["gene_map"]["runner_up_sha256"] = "0" * 64
    tampered = R.gene_map_record(runner, prereg)
    assert tampered["sha256_matches"] is False
    with pytest.raises(R.PreregistrationError, match="differs from the preregistered"):
        R.check_gene_map(tampered, allow=False)
    assert R.check_gene_map(R.gene_map_record(tmp_path / "other.tsv", prereg), allow=False)  # deviation, not refusal
    assert not list(tmp_path.glob(f"{R.OUTPUT_SUBDIR}/*"))  # nothing written


@requires_paths(R.DEFAULT_PREREG)
def test_real_preregistration_file_parses_without_data():
    """Parses the committed (or pending) prereg exactly as a run would; reads no expression data."""
    prereg = R.load_prereg(R.DEFAULT_PREREG)
    rules = P.parse_classification_rules(prereg, source=str(R.DEFAULT_PREREG))
    assert (rules.margin_kappa, rules.interaction_kappa) == (0.5, 1.0)
    assert (rules.terminal_ci_level, rules.persist_ci_level, rules.decay_ci_level, rules.reversed_ci_level) == (0.95, 0.9, 0.9, 0.95)
    settings = R.prereg_settings(prereg, yaml.safe_load(CONFIG.read_text()))
    assert (settings.n_perm, settings.n_boot, settings.min_genes) == (100_000, 10_000, 8)
    assert settings.primary_margin_kappa == 0.5 and 0.5 in settings.margin_kappas
    assert set(settings.sweep_sets) <= set(P.KEY_TARGETS)
