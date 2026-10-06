"""WP-S2: podocyte sensitivity primitives and driver.

Everything runs on synthetic missions except one map-agnostic real-data contract (item 2 equals
the strict-matching full_common_observable target computed by that script's own code path),
which is skipped when the inputs are absent. No new real-data sensitivity is computed here.
"""

from __future__ import annotations

import math
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
import yaml

from scripts.clinical_axes import run_podocyte_sensitivities as driver
from scripts.clinical_axes.run_podocyte_scaffold_specificity import (
    PODOCYTE_SET,
    STRUCTURAL_SET,
    _design,
    _frozen_genes,
    _regression_rows,
)
from scripts.clinical_axes.run_sensitivities import (
    _meta_for_family,
    gene_contributions,
    secondary_qc_covariate_sensitivity,
)
from src.clinical_axes.analysis import combined_score_design, descriptive_gene_effects
from src.clinical_axes.data import MissionData, cpm_eligible_genes
from src.clinical_axes.sensitivity import (
    FAIL,
    K0_MODEL_SPECS,
    NOT_EVALUABLE,
    PASS,
    PreregistrationError,
    absolute_shares,
    ancova_standardized_effect,
    classify_overall,
    common_observable_family,
    direction_retained,
    evaluate_gate,
    evaluate_gene_influence,
    k0_observed,
    leave_top_k_curve,
    parse_gate_rule,
    residualized_mission_effect,
    resolve_top_ks,
    retention_ratio,
    signed_contributions,
    single_set_family,
    standardized_qc_covariate,
    variant_meta,
)
from src.clinical_axes.statistics import random_effects_reml_mkh, score_signed_axis


REPO = Path(__file__).resolve().parents[1]
THRESHOLD = 0.1

POD = [f"Pod{i:02d}" for i in range(24)]
STRUC = [f"Str{i:02d}" for i in range(30)] + POD[:3]
OTHER_SETS = {
    "tubule__high_specificity": [f"Tub{i:02d}" for i in range(12)],
    "immune__high_specificity": [f"Imm{i:02d}" for i in range(10)],
    "tiny__high_specificity": [f"Tny{i:02d}" for i in range(5)],
}
BACKGROUND = [f"Bg{i:03d}" for i in range(40)]
ALL_GENES = list(
    dict.fromkeys(POD + STRUC + [g for genes in OTHER_SETS.values() for g in genes] + BACKGROUND)
)

# Verbatim structure of plan section 4.1 (the real file is written by WP-PR).
PREREG_TEXT = """
analysis_id: clinical_axes_podocyte_sensitivities_v1
parent_config: config/clinical_renal_axes_cross_mission.yaml
status: post_hoc_sensitivities_preregistered_before_computation
lock_date: null            # set at commit; runs refuse a dirty file
seed_base_from_parent: true
seed_offsets: {median_family: 101, common_family: 102, median_common_family: 103,
               lomo_family: [601, 602, 603, 604, 605], atlas_loso: 800, second_atlas: 900}
reference_values: read from the Phase-0 reference run manifest (never typed constants)
targets:
  primary: podocyte__high_specificity
  descriptive_comparators: [broad_structural_scaffold_control__all, structural_disjoint]
  k0_models: [P0, P1, S0, S1, D]      # HC3 per mission -> REML/mKH; observed only
retention_ratio_min: 0.50             # parent yaml LOGO rule; Grey60 gates D/F precedent
gate_rule_default: "estimate has the sign of the reference AND estimate/reference >= 0.50"
items:
  1_signed_median:    {gate: "median pooled estimate > 0", descriptive: [ci, max_t_fwer, rank_abs_t]}
  2_common_intersection: {gate: "pooled estimate > 0",
                          contract: "equals strict-matching full_common_observable target within 1e-6"}
  3_leave_one_mission: {gate: "all 5 LOMO estimates > 0", descriptive: [lomo_ci, lomo_max_t_fwer_fixed_49]}
  4_gene_influence:
    gates: {logo_all_direction: true, logo_min_retention: 0.50, max_single_gene_signed_share: 0.25}
    secondary_bounds_v13_precedent: {top1_positive_share_max: 0.10, top10_positive_share_max: 0.50}
    descriptive: [leave_top_k_curve, k50, k_ci0, n_eff, absolute_share]
  5_osd462_preparation: {variants: [mRNA, UPX], gate: gate_rule_default}
  6_osd253_rerun_control: {variant: rerun_white, gate: gate_rule_default}
  7_osd163_mapping_rate: {primary_method: label_blind_residualization, secondary: ancova_descriptive,
                          gate: gate_rule_default}
  8_atlas_loso: {gate: "every QC-passing rebuild with >=8 podocyte HS genes: gate_rule_default",
                 thresholds: unchanged_frozen, descriptive: [jaccard, rank_abs_t, max_t_fwer, k0_p1]}
  9_second_atlas:
    reference: tabula_muris_senis_droplet_kidney
    evaluability: ">=3 mice with >=25 podocytes; else age-group units; else NOT_EVALUABLE"
    gates: ["TMS podocyte HS estimate: gate_rule_default", "consensus set estimate > 0"]
overall:
  ROBUST_DIRECTION: items 1-7 all pass
  FRAGILE: any of items 1-7 fails (name the items)
  DEFINITION_ROBUST: item 8 passes (and item 9 if run)
k0_rule: >-
  No sensitivity or permutation-scheme change can upgrade K0's verdict (HC3 P1 CI crosses 0 ->
  bounded observation). If P1 <= 0 in any variant, report "podocyte-leaning tail not robust to <variant>".
  K0 primary permutation going forward: freedman_lane_joint; legacy 0.061 retained for audit.
multiplicity: sensitivities are not new discovery families; only family-wide maxT within each variant is reported
provenance_disclosure: >-
  Known before lock: synthetic test copy.
"""


# --------------------------------------------------------------------------------------
# Synthetic fixtures
# --------------------------------------------------------------------------------------


def make_mission(
    name: str,
    seed: int,
    *,
    blocks: tuple[str, ...] = ("all",),
    n_flight: int = 5,
    n_control: int = 5,
    pod_shift: float = 0.9,
    struc_shift: float = 0.5,
    shifted_pod: list[str] | None = None,
    ineligible: tuple[str, ...] = (),
    missing: tuple[str, ...] = (),
) -> MissionData:
    rng = np.random.default_rng(seed)
    animals, conditions, labels = [], [], []
    for block in blocks:
        for condition, count in (("FLT", n_flight), ("GC", n_control)):
            for index in range(count):
                animals.append(f"{name}_{block}_{condition}_{index}")
                conditions.append(condition)
                labels.append(block)
    n = len(animals)
    flight = np.array([c == "FLT" for c in conditions], dtype=float)
    drift = rng.normal(size=n)
    shifted = set(POD if shifted_pod is None else shifted_pod)
    values = rng.normal(size=(len(ALL_GENES), n))
    for row, gene in enumerate(ALL_GENES):
        loading = 0.5 + rng.uniform()
        if gene in STRUC or gene in POD:
            values[row] += 0.6 * drift
        if gene in shifted:
            values[row] += pod_shift * loading * flight
        elif gene in STRUC:
            values[row] += struc_shift * loading * flight
    expression = pd.DataFrame(values + 8.0, index=pd.Index(ALL_GENES, name="gene"), columns=animals)
    expression = expression.drop(index=list(missing))
    counts = pd.DataFrame(
        rng.poisson(200, size=(len(ALL_GENES), n)).astype(float) + 1.0,
        index=pd.Index(ALL_GENES, name="gene"),
        columns=animals,
    )
    counts.loc[list(ineligible)] = 0.0
    metadata = pd.DataFrame(
        {"condition": conditions, "block": labels, "source_sample": animals},
        index=pd.Index(animals, name="animal"),
    )
    qc = pd.DataFrame(
        {"uniquely_mapped_percent": 80.0 + rng.normal(scale=3.0, size=n) + 1.5 * flight},
        index=animals,
    )
    data = MissionData(
        mission=name, expression=expression, counts=counts, metadata=metadata, qc=qc
    )
    data.validate()
    return data


def synthetic_missions(**overrides) -> dict[str, MissionData]:
    return {
        "OSD-102": make_mission("OSD-102", 1, **overrides),
        "OSD-163": make_mission("OSD-163", 2, ineligible=("Pod05",), **overrides),
        "OSD-253": make_mission(
            "OSD-253", 3, blocks=("day25", "day75"), n_control=4, **overrides
        ),
        "OSD-462": make_mission("OSD-462", 4, missing=("Pod07",), **overrides),
        "OSD-771": make_mission("OSD-771", 5, blocks=("YNG", "OLD"), **overrides),
    }


@pytest.fixture(scope="module")
def missions() -> dict[str, MissionData]:
    return synthetic_missions()


def write_tiers(path: Path) -> Path:
    rows = []
    sets = {PODOCYTE_SET: POD, STRUCTURAL_SET: STRUC, **OTHER_SETS}
    for gene_set, genes in sets.items():
        compartment = gene_set.split("__")[0]
        for gene in genes:
            rows.append(
                {
                    "gene_symbol": gene,
                    "gene_set": gene_set,
                    "report_compartment": compartment,
                    "tier": gene_set.split("__")[1],
                    "final_for_testing": True,
                }
            )
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def prereg_raw(**changes) -> dict:
    raw = yaml.safe_load(PREREG_TEXT)
    for dotted, value in changes.items():
        cursor = raw
        keys = dotted.split(".")
        for key in keys[:-1]:
            cursor = cursor[key]
        if value is None:
            cursor.pop(keys[-1], None)
        else:
            cursor[keys[-1]] = value
    return raw


def synthetic_inputs(tmp_path: Path, missions, *, failing_osd253: bool = False):
    tiers = write_tiers(tmp_path / "tiers.tsv")

    def broken():
        raise FileNotFoundError("synthetic rerun-control file absent")

    alternates = {
        "OSD462_mRNA": ("OSD-462", lambda: make_mission("OSD-462", 14)),
        "OSD462_UPX": ("OSD-462", lambda: make_mission("OSD-462", 24)),
        "OSD253_rerun_white": (
            "OSD-253",
            broken
            if failing_osd253
            else (lambda: make_mission("OSD-253", 13, blocks=("day25", "day75"))),
        ),
    }
    return driver.PipelineInputs(
        missions=dict(missions),
        alternates=alternates,
        tiers_path=tiers,
        threshold=THRESHOLD,
        seed_base=20260811,
    )


# --------------------------------------------------------------------------------------
# Families, ratios, variant meta
# --------------------------------------------------------------------------------------


def test_single_set_family_has_frozen_shape_and_dedupes():
    family = single_set_family(["A", "B", "A", "C"], "x", minimum=3)
    assert family == {
        "x": {
            "role": "podocyte_sensitivity_target",
            "subdomains": {
                "atlas_markers": {"genes": {"A": 1, "B": 1, "C": 1}, "minimum_present": 3}
            },
        }
    }
    with pytest.raises(ValueError):
        single_set_family(["A", "B"], "x", minimum=3)


def test_retention_ratio_with_zero_reference_is_na():
    assert math.isnan(retention_ratio(0.3, 0.0))
    assert math.isnan(retention_ratio(0.3, float("nan")))
    assert math.isnan(retention_ratio(float("nan"), 0.6))
    assert retention_ratio(0.3, 0.6) == pytest.approx(0.5)
    assert retention_ratio(-0.3, 0.6) == pytest.approx(-0.5)
    assert direction_retained(0.3, 0.0) is None
    assert direction_retained(-0.3, 0.6) is False
    assert direction_retained(-0.3, -0.6) is True
    rule = parse_gate_rule(yaml.safe_load(PREREG_TEXT)["gate_rule_default"])
    result = evaluate_gate(rule, estimate=0.3, reference=0.0)
    assert result["verdict"] == NOT_EVALUABLE  # never PASS on an undefined ratio


def test_variant_meta_counts_positive_missions(missions):
    family = single_set_family(POD, PODOCYTE_SET)
    result = variant_meta(missions, family, threshold=THRESHOLD)
    meta, effects, _ = _meta_for_family(missions, family, THRESHOLD)
    assert result.meta["estimate"].iloc[0] == meta["estimate"].iloc[0]
    assert result.meta["n_positive_missions"].iloc[0] == int((effects["estimate"] > 0).sum())
    assert result.meta["n_negative_missions"].iloc[0] == int((effects["estimate"] < 0).sum())


def test_common_observable_family_requires_eligibility_and_expression(missions):
    family = single_set_family(POD, PODOCYTE_SET)
    common = common_observable_family(missions, family, THRESHOLD)
    genes = list(common[PODOCYTE_SET]["subdomains"]["atlas_markers"]["genes"])
    assert "Pod05" not in genes  # CPM-ineligible in OSD-163
    assert "Pod07" not in genes  # absent from the OSD-462 expression matrix
    assert len(genes) == len(POD) - 2
    assert common[PODOCYTE_SET]["subdomains"]["atlas_markers"]["minimum_present"] == len(genes)
    with pytest.raises(ValueError):
        common_observable_family(missions, family, THRESHOLD, minimum=len(POD))


# --------------------------------------------------------------------------------------
# K0 observed
# --------------------------------------------------------------------------------------


def _manual_k0(missions, pod, struc, extra=None):
    rows = []
    for mission, data in missions.items():
        eligible = set(cpm_eligible_genes(data, THRESHOLD))
        pod_used = [g for g in pod if g in eligible and g in data.expression.index]
        struc_used = [g for g in struc if g in eligible and g in data.expression.index]
        podocyte = score_signed_axis(data.expression, {g: 1 for g in pod_used}).scores
        structural = score_signed_axis(data.expression, {g: 1 for g in struc_used}).scores
        scored = {"podocyte": podocyte, "structural": structural, "diff": podocyte - structural}
        for label, outcome, adjuster in K0_MODEL_SPECS:
            if extra and mission in extra:
                design = _design(data.metadata, scored[adjuster] if adjuster else None)
                for column in extra[mission]:
                    design[column] = extra[mission].loc[data.metadata.index, column].to_numpy()
                fit = sm.OLS(scored[outcome], design).fit().get_robustcov_results(cov_type="HC3")
                index = list(design.columns).index("flight")
                hc3 = {"estimate": fit.params[index], "variance": fit.bse[index] ** 2}
            else:
                hc3 = next(
                    row
                    for row in _regression_rows(
                        mission,
                        data.metadata,
                        scored[outcome],
                        scored[adjuster] if adjuster else None,
                        label,
                    )
                    if row["variance_type"] == "HC3"
                )
            rows.append({"mission": mission, "model": label, **hc3})
    return pd.DataFrame(rows)


def test_k0_observed_matches_regression_rows_hc3(missions):
    disjoint = [g for g in STRUC if g not in set(POD)]
    result = k0_observed(missions, POD, disjoint, threshold=THRESHOLD)
    manual = _manual_k0(missions, POD, disjoint)
    merged = result.mission_effects.merge(manual, on=["mission", "model"])
    assert len(merged) == len(missions) * len(K0_MODEL_SPECS)
    assert np.allclose(merged["coef_score_units"], merged["estimate"], rtol=0, atol=1e-12)
    assert np.allclose(merged["variance_x"], merged["variance_y"], rtol=0, atol=1e-12)
    meta = result.meta.set_index("model")
    for label, _, _ in K0_MODEL_SPECS:
        sub = manual[manual["model"] == label]
        fit = random_effects_reml_mkh(sub["estimate"].to_numpy(), sub["variance"].to_numpy())
        assert meta.loc[label, "coef_score_units"] == pytest.approx(fit.estimate, abs=1e-12)
        assert meta.loc[label, "ci_low_mkh"] == pytest.approx(fit.ci_low, abs=1e-12)
        assert meta.loc[label, "ci_high_mkh"] == pytest.approx(fit.ci_high, abs=1e-12)
        assert meta.loc[label, "p_mkh"] == pytest.approx(fit.p, abs=1e-12)
        assert meta.loc[label, "i_squared"] == pytest.approx(fit.i_squared, abs=1e-9)
    assert "g" not in result.meta.columns and "estimate" not in result.meta.columns


def test_k0_observed_extra_covariate_enters_only_its_mission(missions):
    disjoint = [g for g in STRUC if g not in set(POD)]
    z = standardized_qc_covariate(missions["OSD-163"], "uniquely_mapped_percent")
    extra = {"OSD-163": pd.DataFrame({"z_mapping": z})}
    adjusted = k0_observed(missions, POD, disjoint, threshold=THRESHOLD, extra_covariates=extra)
    plain = k0_observed(missions, POD, disjoint, threshold=THRESHOLD)
    manual = _manual_k0(missions, POD, disjoint, extra=extra)
    merged = adjusted.mission_effects.merge(manual, on=["mission", "model"])
    assert np.allclose(merged["coef_score_units"], merged["estimate"], atol=1e-12)
    other = adjusted.mission_effects["mission"] != "OSD-163"
    assert np.array_equal(
        adjusted.mission_effects.loc[other, "coef_score_units"].to_numpy(),
        plain.mission_effects.loc[other, "coef_score_units"].to_numpy(),
    )
    assert set(adjusted.mission_effects.loc[~other, "extra_covariates"]) == {"z_mapping"}
    with pytest.raises(ValueError):
        k0_observed(missions, POD, disjoint, threshold=THRESHOLD, extra_covariates={"OSD-999": extra["OSD-163"]})
    with pytest.raises(ValueError):
        k0_observed(
            missions,
            POD,
            disjoint,
            threshold=THRESHOLD,
            extra_covariates={"OSD-163": pd.DataFrame({"flight": z})},
        )


# --------------------------------------------------------------------------------------
# Signed contributions
# --------------------------------------------------------------------------------------


def _details(missions, family, method="mean"):
    return combined_score_design(missions, family, cpm_threshold=THRESHOLD, method=method)[3]


@pytest.mark.parametrize(
    "family",
    [
        single_set_family(POD, PODOCYTE_SET),
        {
            "two_domain": {
                "role": "test",
                "subdomains": {
                    "a": {"genes": {g: 1 for g in POD[:10]}, "minimum_present": 6},
                    "b": {
                        "genes": {**{g: 1 for g in STRUC[:12]}, "Str20": -1},
                        "minimum_present": 8,
                    },
                },
            }
        },
    ],
    ids=["single_set", "two_subdomains_with_negative_sign"],
)
def test_signed_contributions_sum_exactly_to_pooled_estimate(missions, family):
    result = signed_contributions(missions, family, _details(missions, family))
    meta, effects, _ = _meta_for_family(missions, family, THRESHOLD)
    axis = next(iter(family))
    pooled = float(meta["estimate"].iloc[0])
    genes = result.genes[result.genes["target"] == axis]
    assert abs(genes["signed_contribution"].sum() - pooled) < 1e-10
    assert genes["signed_share"].sum() == pytest.approx(1.0, abs=1e-9)
    for mission in missions:
        mission_effect = float(effects.loc[effects["mission"] == mission, "estimate"].iloc[0])
        assert abs(genes[f"contribution_{mission}"].sum() - mission_effect) < 1e-10
    concentration = result.concentration.iloc[0]
    assert abs(concentration["decomposition_residual"]) < 1e-10
    assert (result.strata["residual"].abs() < 1e-10).all()
    positive = genes.loc[genes["signed_contribution"] > 0, "signed_contribution"]
    assert concentration["n_eff_positive"] == pytest.approx(
        positive.sum() ** 2 / (positive**2).sum()
    )
    assert concentration["top1_positive_share"] == pytest.approx(positive.max() / positive.sum())
    assert concentration["max_signed_share"] == pytest.approx(genes["signed_share"].max())
    # Genes never scored in a mission carry NaN there, not a fake zero.
    if axis == PODOCYTE_SET:
        pod05 = genes.set_index("gene").loc["Pod05"]
        assert np.isnan(pod05["contribution_OSD-163"]) and pod05["n_missions"] == 4


def test_signed_contributions_refuse_median_scores(missions):
    family = single_set_family(POD, PODOCYTE_SET)
    with pytest.raises(ValueError, match="mean method"):
        signed_contributions(missions, family, _details(missions, family, method="median"))


def test_absolute_shares_match_existing_gene_contributions(missions, tmp_path):
    family = single_set_family(POD, PODOCYTE_SET)
    effects = descriptive_gene_effects(missions, _details(missions, family))
    path = tmp_path / "gene_effects.tsv"
    effects.to_csv(path, sep="\t", index=False)
    existing = gene_contributions(path).set_index("gene")
    ours = absolute_shares(effects).set_index("gene")
    assert np.allclose(
        ours.loc[existing.index, "absolute_share"],
        existing["absolute_contribution_share"],
        atol=1e-12,
    )
    assert np.allclose(ours.loc[existing.index, "pooled_signed_g"], existing["pooled_signed_g"])


# --------------------------------------------------------------------------------------
# Leave-top-k curve
# --------------------------------------------------------------------------------------


def test_leave_top_k_curve_nests_removals_and_shrinks_the_set(missions):
    family = single_set_family(POD, PODOCYTE_SET)
    contributions = signed_contributions(missions, family, _details(missions, family)).genes
    series = contributions.set_index("gene")["signed_contribution"]
    result = leave_top_k_curve(
        missions, POD, series, threshold=THRESHOLD, ks=(1, 2, 5, "25%", 3), name=PODOCYTE_SET
    )
    curve = result.curve
    assert list(curve["k_label"])[0] == "0"
    assert curve["n_removed"].is_monotonic_increasing
    assert curve["n_genes_defined_remaining"].is_monotonic_decreasing
    removed = [set(filter(None, text.split("|"))) for text in curve["removed_genes"]]
    for smaller, larger in zip(removed, removed[1:]):
        assert smaller <= larger
    ranked = series[series > 0].sort_values(ascending=False)
    assert curve.loc[curve["n_removed"] == 2, "removed_genes"].iloc[0] == "|".join(ranked.index[:2])
    baseline = _meta_for_family(missions, family, THRESHOLD)[0]["estimate"].iloc[0]
    assert curve["estimate"].iloc[0] == pytest.approx(baseline, abs=1e-12)
    assert curve.loc[curve["k_label"] == "25%", "k_requested"].iloc[0] == math.ceil(0.25 * len(series))
    assert set(curve["status"]) == {"descriptive_influence_diagnostic"}
    assert resolve_top_ks(["10%"], 157) == [("10%", 16)]
    assert resolve_top_ks(["25%"], 157) == [("25%", 40)]


def test_leave_top_k_identifies_concentrated_signal():
    concentrated = synthetic_missions(shifted_pod=POD[:3], pod_shift=4.0, struc_shift=0.0)
    family = single_set_family(POD, PODOCYTE_SET)
    series = (
        signed_contributions(concentrated, family, _details(concentrated, family))
        .genes.set_index("gene")["signed_contribution"]
    )
    result = leave_top_k_curve(concentrated, POD, series, threshold=THRESHOLD, ks=(1, 2, 3, 5, 1000))
    assert result.k50 is not None and result.k50 <= 5
    assert result.k_ci0 is not None
    last = result.curve.iloc[-1]
    assert bool(last["capped_at_n_positive"]) and last["n_removed"] == (series > 0).sum()


# --------------------------------------------------------------------------------------
# Gate evaluation (rules read from the preregistration)
# --------------------------------------------------------------------------------------

_DEFAULT = "gate_rule_default"
GATE_CASES = [
    ("pooled estimate > 0", {"estimate": 0.4}, PASS),
    ("pooled estimate > 0", {"estimate": 0.0}, FAIL),
    ("pooled estimate > 0", {"estimate": -0.1}, FAIL),
    ("pooled estimate > 0", {"estimate": float("nan")}, NOT_EVALUABLE),
    ("median pooled estimate > 0", {"estimate": 0.2}, PASS),
    ("median pooled estimate > 0", {"estimate": -0.2}, FAIL),
    ("all 5 LOMO estimates > 0", {"lomo_estimates": [0.1, 0.2, 0.3, 0.4, 0.5]}, PASS),
    ("all 5 LOMO estimates > 0", {"lomo_estimates": [0.1, -0.2, 0.3, 0.4, 0.5]}, FAIL),
    ("all 5 LOMO estimates > 0", {"lomo_estimates": [0.1, 0.2, 0.3, 0.4]}, NOT_EVALUABLE),
    (_DEFAULT, {"estimate": 0.4, "reference": 0.6}, PASS),
    (_DEFAULT, {"estimate": 0.3, "reference": 0.6}, PASS),  # ratio exactly 0.50 (>=)
    (_DEFAULT, {"estimate": 0.29, "reference": 0.6}, FAIL),
    (_DEFAULT, {"estimate": -0.4, "reference": 0.6}, FAIL),
    (_DEFAULT, {"estimate": 1.2, "reference": 0.6}, PASS),
    (_DEFAULT, {"estimate": -0.5, "reference": -0.6}, PASS),
    (_DEFAULT, {"estimate": 0.3, "reference": 0.0}, NOT_EVALUABLE),
    (_DEFAULT, {"estimate": float("nan"), "reference": 0.6}, NOT_EVALUABLE),
]


@pytest.mark.parametrize("rule_text, inputs, expected", GATE_CASES)
def test_gate_evaluator_table(rule_text, inputs, expected):
    raw = yaml.safe_load(PREREG_TEXT)
    rule = parse_gate_rule(rule_text, aliases={_DEFAULT: raw[_DEFAULT]})
    assert evaluate_gate(rule, **inputs)["verdict"] == expected


def test_preregistration_rules_are_read_from_yaml():
    prereg = driver.parse_preregistration(yaml.safe_load(PREREG_TEXT))
    default = " ".join(yaml.safe_load(PREREG_TEXT)[_DEFAULT].split())
    assert prereg.rules[5].text == prereg.rules[6].text == prereg.rules[7].text == default
    assert prereg.rules[3].clauses[0].count == 5
    assert prereg.required_items == (1, 2, 3, 4, 5, 6, 7)
    assert (prereg.robust_label, prereg.fragile_label) == ("ROBUST_DIRECTION", "FRAGILE")
    assert prereg.k0_rule_model == "P1_podocyte_adjusted"
    assert prereg.k0_note_template == "podocyte-leaning tail not robust to <variant>"
    assert prereg.osd462_variants == ("mRNA", "UPX") and prereg.osd253_variant == "rerun_white"
    # A changed threshold in the YAML changes the evaluation (nothing is hard-coded).
    stricter = driver.parse_preregistration(
        prereg_raw(**{
            "retention_ratio_min": 0.8,
            _DEFAULT: "estimate has the sign of the reference AND estimate/reference >= 0.80",
        })
    )
    assert evaluate_gate(prereg.rules[5], estimate=0.4, reference=0.6)["verdict"] == PASS
    assert evaluate_gate(stricter.rules[5], estimate=0.4, reference=0.6)["verdict"] == FAIL


@pytest.mark.parametrize(
    "changes",
    [
        {_DEFAULT: "estimate is broadly similar to the reference"},
        {_DEFAULT: "estimate has the sign of the reference AND estimate/reference >= 0.60"},
        {"items.3_leave_one_mission.gate": "pooled estimate > 0"},
        {"items.5_osd462_preparation.variants": ["polyA"]},
        {"items.7_osd163_mapping_rate.primary_method": "ancova"},
        {"items.4_gene_influence.gates": {"logo_all_direction": True, "unknown_gate": 1}},
        {"items.6_osd253_rerun_control": None},
        {"overall": {"ROBUST_DIRECTION": "items 1-7 all pass"}},
        {"targets.k0_models": ["P1", "Q9"]},
    ],
)
def test_preregistration_parser_refuses_inconsistent_or_unknown_rules(changes):
    with pytest.raises(PreregistrationError):
        driver.parse_preregistration(prereg_raw(**changes))


def _logo(directions, ratios):
    return pd.DataFrame({"direction_retained": directions, "retention_ratio": ratios})


GENE_INFLUENCE_CASES = [
    (_logo([True] * 3, [0.9, 0.8, 0.7]), {"max_signed_share": 0.2, "top1_positive_share": 0.5}, PASS),
    (_logo([True, False, True], [0.9, -0.1, 0.7]), {"max_signed_share": 0.2}, FAIL),
    (_logo([True] * 3, [0.9, 0.45, 0.7]), {"max_signed_share": 0.2}, FAIL),
    (_logo([True] * 3, [0.9, 0.8, 0.7]), {"max_signed_share": 0.30}, FAIL),
    (_logo([True] * 3, [0.9, 0.8, 0.7]), {"max_signed_share": float("nan")}, NOT_EVALUABLE),
    (_logo([], []), {"max_signed_share": 0.2}, NOT_EVALUABLE),
]


@pytest.mark.parametrize("logo, concentration, expected", GENE_INFLUENCE_CASES)
def test_gene_influence_gate_table(logo, concentration, expected):
    spec = yaml.safe_load(PREREG_TEXT)["items"]["4_gene_influence"]
    result = evaluate_gene_influence(
        spec["gates"],
        spec["secondary_bounds_v13_precedent"],
        logo=logo,
        concentration=concentration,
    )
    assert result["verdict"] == expected
    # The v13 precedent bounds are reported, never gating (a 0.5 top-1 share exceeds 0.10).
    if "top1_positive_share" in concentration:
        bound = next(b for b in result["secondary_bounds"] if b["bound"] == "top1_positive_share_max")
        assert bound["within_bound"] is False and expected == PASS


@pytest.mark.parametrize(
    "verdicts, expected, failed",
    [
        ({n: PASS for n in range(1, 8)}, "ROBUST_DIRECTION", []),
        ({**{n: PASS for n in range(1, 8)}, 3: FAIL}, "FRAGILE", [3]),
        ({**{n: PASS for n in range(1, 8)}, 3: FAIL, 6: NOT_EVALUABLE}, "FRAGILE", [3]),
        ({**{n: PASS for n in range(1, 8)}, 6: NOT_EVALUABLE}, "INCOMPLETE", []),
        ({n: PASS for n in range(1, 7)}, "INCOMPLETE", []),
    ],
)
def test_overall_classification(verdicts, expected, failed):
    result = classify_overall(
        verdicts, range(1, 8), robust_label="ROBUST_DIRECTION", fragile_label="FRAGILE"
    )
    assert result["classification"] == expected
    assert result["failed_items"] == failed


# --------------------------------------------------------------------------------------
# Item-7 helpers
# --------------------------------------------------------------------------------------


def test_residualization_matches_frozen_contract_and_ancova_formula(missions):
    family = single_set_family(POD, PODOCYTE_SET)
    contract = secondary_qc_covariate_sensitivity(missions, family, THRESHOLD).iloc[0]
    data = missions["OSD-163"]
    score = _details(missions, family)["OSD-163"][PODOCYTE_SET].scores
    z = standardized_qc_covariate(data, "uniquely_mapped_percent")
    adjusted = residualized_mission_effect(score, data.metadata, z)
    assert adjusted["estimate"] == pytest.approx(contract["osd163_adjusted_estimate"], abs=1e-12)

    ancova = ancova_standardized_effect(score, data.metadata, z)
    flight = (data.metadata["condition"] == "FLT").astype(float).to_numpy()
    matrix = np.column_stack([np.ones(len(flight)), flight, z.to_numpy()])
    beta = np.linalg.lstsq(matrix, score.to_numpy(), rcond=None)[0][1]
    s_t = score[flight == 1].to_numpy()
    s_c = score[flight == 0].to_numpy()
    nt, nc = len(s_t), len(s_c)
    sd = math.sqrt(((nt - 1) * s_t.var(ddof=1) + (nc - 1) * s_c.var(ddof=1)) / (nt + nc - 2))
    df = nt + nc - 3
    g = (1 - 3 / (4 * df - 1)) * beta / sd
    assert ancova["df"] == df
    assert ancova["estimate"] == pytest.approx(g, abs=1e-12)
    assert ancova["variance"] == pytest.approx((nt + nc) / (nt * nc) + g**2 / (2 * df), abs=1e-12)


# --------------------------------------------------------------------------------------
# Preregistration guard
# --------------------------------------------------------------------------------------

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")


def _git(repo: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-c", "user.email=test@example.invalid", "-c", "user.name=test", *args],
        cwd=repo,
        check=True,
        capture_output=True,
    )


def test_guard_refuses_missing_preregistration_even_when_allowed(tmp_path):
    for allow in (False, True):
        with pytest.raises(PreregistrationError, match="does not exist"):
            driver.preregistration_guard(
                tmp_path / "absent.yaml", repo=tmp_path, allow_uncommitted=allow
            )
    assert driver.main(["--prereg", str(tmp_path / "absent.yaml")]) == 2


@needs_git
def test_guard_requires_a_committed_clean_file(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    prereg = repo / "prereg.yaml"
    prereg.write_text(PREREG_TEXT)

    with pytest.raises(PreregistrationError, match="not committed"):
        driver.preregistration_guard(prereg, repo=repo)  # untracked
    allowed = driver.preregistration_guard(prereg, repo=repo, allow_uncommitted=True)
    assert allowed["committed"] is False and allowed["allow_uncommitted_prereg"] is True

    _git(repo, "add", "prereg.yaml")
    with pytest.raises(PreregistrationError):
        driver.preregistration_guard(prereg, repo=repo)  # staged, not committed
    _git(repo, "commit", "-q", "-m", "lock")
    record = driver.preregistration_guard(prereg, repo=repo)
    assert record["committed"] is True and record["tracked"] is True
    assert len(record["git_head"]) == 40 and record["prereg_last_commit"] == record["git_head"]
    assert record["sha256"] == driver.sha256(prereg)

    prereg.write_text(PREREG_TEXT + "\n# edited after lock\n")
    with pytest.raises(PreregistrationError, match="not committed"):
        driver.preregistration_guard(prereg, repo=repo)
    parsed, dirty = driver.load_preregistration(prereg, repo=repo, allow_uncommitted=True)
    assert dirty["committed"] is False and parsed.analysis_id.endswith("_v1")


def test_validate_prereg_only_reads_no_data(tmp_path, monkeypatch, capsys):
    prereg = tmp_path / "prereg.yaml"
    prereg.write_text(PREREG_TEXT)

    def forbidden(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("validation must not load mission data")

    monkeypatch.setattr(driver, "build_inputs", forbidden)
    assert driver.main(["--prereg", str(prereg), "--validate-prereg-only"]) == 2  # untracked
    code = driver.main(
        ["--prereg", str(prereg), "--validate-prereg-only", "--allow-uncommitted-prereg"]
    )
    assert code == 0
    assert "3_leave_one_mission" in capsys.readouterr().out


def test_gene_map_override_is_in_memory_only(tmp_path):
    config = yaml.safe_load((REPO / "config/clinical_renal_axes_cross_mission.yaml").read_text())
    original = config["gene_mapping"]["path"]
    gene_map = tmp_path / "map.tsv"
    gene_map.write_text("ensembl_gene_id\tmgi_symbol\n")
    replaced = driver.config_with_gene_map(config, gene_map)
    assert replaced["gene_mapping"]["path"] == str(gene_map.resolve())
    assert config["gene_mapping"]["path"] == original
    assert driver.config_with_gene_map(config, None) == config
    with pytest.raises(FileNotFoundError):
        driver.config_with_gene_map(config, tmp_path / "absent.tsv")


# --------------------------------------------------------------------------------------
# Driver end to end on synthetic missions
# --------------------------------------------------------------------------------------


def test_pipeline_end_to_end_outputs_and_gates(missions, tmp_path):
    prereg = driver.parse_preregistration(yaml.safe_load(PREREG_TEXT))
    inputs = synthetic_inputs(tmp_path, missions)
    options = driver.PipelineOptions(lomo_permutations=19, chunk_size=8, k0_leave_one_gene=True)
    result = driver.run_pipeline(inputs, prereg, options)
    outdir = tmp_path / "out"
    driver.write_outputs(result, outdir, {"analysis": "synthetic"})

    for filename in list(driver.OUTPUT_FILES.values()) + [driver.GATES_FILE, driver.MANIFEST_FILE]:
        assert (outdir / filename).is_file(), filename
    summary = pd.read_csv(outdir / "podocyte_sensitivity_summary.tsv", sep="\t")
    assert list(summary.columns[: len(driver.SUMMARY_COLUMNS)]) == driver.SUMMARY_COLUMNS
    missions_table = pd.read_csv(outdir / "podocyte_sensitivity_mission_effects.tsv", sep="\t")
    assert list(missions_table.columns) == driver.MISSION_COLUMNS
    lomo = pd.read_csv(outdir / "podocyte_leave_one_mission.tsv", sep="\t")
    assert list(lomo.columns[: len(driver.LOMO_COLUMNS)]) == driver.LOMO_COLUMNS
    k0 = pd.read_csv(outdir / "podocyte_k0_sensitivities.tsv", sep="\t")
    assert list(k0.columns[: len(driver.K0_COLUMNS)]) == driver.K0_COLUMNS
    assert "g" not in k0.columns
    contributions = pd.read_csv(outdir / "podocyte_gene_contributions.tsv", sep="\t")
    for column in ["gene", "n_missions", "pooled_signed_g", "absolute_share", "signed_contribution", "signed_share"]:
        assert column in contributions.columns
    assert all(f"contribution_{m}" in contributions.columns for m in missions)

    # Step 1 reference equals the frozen _meta_for_family on the same inputs.
    family = single_set_family(POD, PODOCYTE_SET)
    reference = _meta_for_family(missions, family, THRESHOLD)[0]["estimate"].iloc[0]
    row = summary[(summary["item"] == "0_reference") & (summary["target"] == PODOCYTE_SET)].iloc[0]
    assert row["estimate"] == pytest.approx(reference, abs=1e-12)

    # Item 2 uses the common family; the comparators are never gated.
    common = summary[summary["item"] == "2_common_intersection"].set_index("target")
    expected = variant_meta(
        missions, common_observable_family(missions, family, THRESHOLD), threshold=THRESHOLD
    ).meta["estimate"].iloc[0]
    assert common.loc[PODOCYTE_SET, "estimate"] == pytest.approx(expected, abs=1e-12)
    comparators = summary[summary["target"] != PODOCYTE_SET]
    assert comparators["gate_pass"].isna().all()

    # LOMO family p is attached to family members only (structural_disjoint is not one).
    by_target = lomo.groupby("target")["lomo_max_t_fwer"]
    assert by_target.apply(lambda s: s.notna().all())[PODOCYTE_SET]
    assert by_target.apply(lambda s: s.isna().all())["structural_disjoint"]
    seeds = result.info["lomo_family_permutation"]["seeds_by_omitted_mission"]
    assert seeds == {m: 20260811 + o for m, o in zip(missions, [601, 602, 603, 604, 605])}

    gates = result.gates
    assert set(gates["items"]) == {prereg.item_keys[n] for n in range(1, 8)}
    for key, item in gates["items"].items():
        assert item["verdict"] in {PASS, FAIL, NOT_EVALUABLE}, key
    assert gates["overall"]["classification"] in {"ROBUST_DIRECTION", "FRAGILE", "INCOMPLETE"}
    assert gates["items_not_run_by_this_driver"] == ["8_atlas_loso", "9_second_atlas"]
    item1 = summary[(summary["item"] == "1_signed_median") & (summary["target"] == PODOCYTE_SET)]
    assert bool(item1["gate_pass"].iloc[0]) == (gates["items"]["1_signed_median"]["verdict"] == PASS)
    # The synthetic podocyte shift is strong and diffuse, so the direction gates pass.
    for key in ("1_signed_median", "2_common_intersection", "3_leave_one_mission"):
        assert gates["items"][key]["verdict"] == PASS
    assert len(result.tables["k0_leave_one_gene"]) == len(POD)  # every gene is used somewhere


def test_pipeline_marks_missing_alternate_not_evaluable(missions, tmp_path):
    prereg = driver.parse_preregistration(yaml.safe_load(PREREG_TEXT))
    inputs = synthetic_inputs(tmp_path, missions, failing_osd253=True)
    options = driver.PipelineOptions(structural_gene_influence=False, top_ks=(1, 2))
    result = driver.run_pipeline(inputs, prereg, options)
    item6 = result.gates["items"]["6_osd253_rerun_control"]
    assert item6["verdict"] == NOT_EVALUABLE
    assert result.gates["overall"]["classification"] != "ROBUST_DIRECTION"
    k0 = result.tables["k0"]
    rows = k0[k0["variant"] == "OSD253_rerun_white"]
    assert rows["coef_score_units"].isna().all()
    assert rows["notes"].str.contains("not evaluable").all()


def test_reference_run_mismatch_is_refused(missions, tmp_path):
    prereg = driver.parse_preregistration(yaml.safe_load(PREREG_TEXT))
    inputs = synthetic_inputs(tmp_path, missions)
    reference = tmp_path / "reference"
    reference.mkdir()
    pd.DataFrame(
        {"axis": [PODOCYTE_SET], "estimate": [9.0], "ci_low_mkh": [8.0], "ci_high_mkh": [10.0]}
    ).to_csv(reference / "compartment_context_meta_results.tsv", sep="\t", index=False)
    options = driver.PipelineOptions(reference_run=reference, structural_gene_influence=False)
    with pytest.raises(RuntimeError, match="differs from"):
        driver.run_pipeline(inputs, prereg, options)


def test_gene_maps_mode_writes_idmap_sensitivity(missions, tmp_path, monkeypatch):
    prereg_path = tmp_path / "prereg.yaml"
    prereg_path.write_text(PREREG_TEXT)
    tiers = write_tiers(tmp_path / "tiers.tsv")
    winner, runner_up = tmp_path / "winner.tsv", tmp_path / "runner_up.tsv"
    winner.write_text("ensembl_gene_id\tmgi_symbol\n")
    runner_up.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG0\tX\n")
    alternative = synthetic_missions(pod_shift=-0.6)  # stands in for a map that flips results

    def fake_inputs(config, tiers_path, prereg):
        chosen = alternative if config["gene_mapping"]["path"] == str(runner_up.resolve()) else missions
        inputs = synthetic_inputs(tmp_path, chosen)
        inputs.tiers_path = Path(tiers_path)
        return inputs

    monkeypatch.setattr(driver, "build_inputs", fake_inputs)
    code = driver.main(
        [
            "--prereg", str(prereg_path),
            "--allow-uncommitted-prereg",
            "--tiers", str(tiers),
            "--results", str(tmp_path / "run"),
            "--gene-maps", str(winner), str(runner_up),
            "--skip-structural-gene-influence",
        ]
    )
    assert code == 0
    root = tmp_path / "run" / driver.OUTPUT_SUBDIR
    assert (root / driver.GATES_FILE).is_file()
    assert (root / "idmap_runner_up" / driver.GATES_FILE).is_file()
    table = pd.read_csv(root / "idmap_sensitivity.tsv", sep="\t")
    assert set(table["map"]) == {"winner", "runner_up"}
    assert {"output", "map", "value", "verdict", "verdict_flip", "conservative_verdict"} <= set(table.columns)
    primary = table[table["output"] == f"{PODOCYTE_SET}__primary_estimate"].set_index("map")
    assert primary.loc["winner", "verdict"] != primary.loc["runner_up", "verdict"]
    assert primary["verdict_flip"].all()
    assert set(primary["conservative_verdict"]) == {"non_positive"}
    manifest = yaml.safe_load((root / driver.MANIFEST_FILE).read_text())
    assert manifest["gene_map"]["sha256"] == driver.sha256(winner)
    assert manifest["preregistration"]["sha256"] == driver.sha256(prereg_path)
    assert manifest["preregistration"]["committed"] is False


def test_idmap_table_takes_conservative_verdict():
    rows = []
    for label, item3, overall in (("winner", PASS, "ROBUST_DIRECTION"), ("runner_up", FAIL, "FRAGILE")):
        base = {"map": label, "map_path": label, "map_sha256": None, "ci_low": np.nan, "ci_high": np.nan}
        rows += [
            {**base, "output": "gate__3_leave_one_mission", "value": 0.1, "verdict": item3},
            {**base, "output": "gate__2_common_intersection", "value": 0.5, "verdict": PASS},
            {**base, "output": "overall_classification", "value": np.nan, "verdict": overall},
        ]
    table, summary = driver.idmap_sensitivity_table(
        rows, robust_label="ROBUST_DIRECTION", fragile_label="FRAGILE"
    )
    assert summary["flipped_outputs"] == ["gate__3_leave_one_mission", "overall_classification"]
    assert summary["robust_to_gene_id_reconstruction"] is False
    assert summary["conservative_overall"] == "FRAGILE"
    flags = table.drop_duplicates("output").set_index("output")
    assert flags.loc["gate__3_leave_one_mission", "conservative_verdict"] == FAIL
    assert not flags.loc["gate__2_common_intersection", "verdict_flip"]


# --------------------------------------------------------------------------------------
# Real-data contract (map-agnostic self-consistency; skipped without inputs)
# --------------------------------------------------------------------------------------

CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
REAL_INPUTS = [
    TIERS,
    REPO / "data/processed/resources/id_map.tsv",
    REPO / "data/raw/metadata/s_OSD-771.txt",
    REPO / "data/external/osdr/OSD-102",
]


@pytest.mark.skipif(
    not all(path.exists() for path in REAL_INPUTS), reason="real mission/tier inputs absent"
)
def test_item2_podocyte_common_intersection_equals_strict_matching_target():
    """Item 2 contract: the podocyte common-intersection estimate is K4's
    full_common_observable target, computed through the strict-matching script's own path."""

    from scripts.clinical_axes.run_matched_panel_null import combined_gene_z
    from scripts.clinical_axes.run_strict_podocyte_matching_audit import (
        _observed_meta,
        _reference_covariates,
    )
    from src.clinical_axes.data import load_primary_missions

    config = yaml.safe_load(CONFIG.read_text())
    missions = load_primary_missions(config, REPO)
    threshold = float(config["eligibility"]["cpm_threshold"])
    tiers = pd.read_csv(TIERS, sep="\t")
    family = single_set_family(_frozen_genes(tiers, PODOCYTE_SET), PODOCYTE_SET)
    common_family = common_observable_family(missions, family, threshold)
    ours = variant_meta(missions, common_family, threshold=threshold).meta.iloc[0]

    common = set.intersection(*(cpm_eligible_genes(d, threshold) for d in missions.values()))
    common &= set.intersection(*(set(d.expression.index) for d in missions.values()))
    _, target, candidates = _reference_covariates(missions, common, TIERS)
    gene_z, design = combined_gene_z(missions, sorted(set(target + candidates)))
    strict, _ = _observed_meta(
        pd.DataFrame({"full_common_observable__target": gene_z[target].mean(axis=1)}),
        design,
    )
    strict = strict.iloc[0]
    assert sorted(common_family[PODOCYTE_SET]["subdomains"]["atlas_markers"]["genes"]) == target
    for column in ("estimate", "ci_low_mkh", "ci_high_mkh", "p_mkh"):
        assert abs(float(ours[column]) - float(strict[column])) < 1e-6, column
