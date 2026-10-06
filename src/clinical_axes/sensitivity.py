"""Observed-data sensitivity primitives for the frozen podocyte set, its structural comparators and K0.

Scoring, Hedges g, within-mission fixed-effect strata and REML/mKH pooling are imported from the
frozen pipeline and never re-implemented, so a sensitivity differs from the primary analysis only
by its declared perturbation. Nothing here chooses genes, samples or thresholds from flight labels.

K0 quantities are raw score-unit OLS flight coefficients (``coef_score_units``). They are not
Hedges g and are never labelled ``g``.

Gate rules are read from the preregistration YAML through a deliberately small grammar. Rule text
outside that grammar is refused, never guessed, so a reworded preregistration fails loudly instead
of being evaluated against a rule nobody wrote.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import math
import operator
import re
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm

# Frozen stable contracts (plan section 0). Importing them, rather than copying, keeps the
# sensitivity arithmetic identical to the published pipeline by construction.
from scripts.clinical_axes.run_podocyte_scaffold_specificity import _design
from scripts.clinical_axes.run_sensitivities import (
    _meta_for_family,
    common_intersection_family,
)

from .analysis import (
    axis_directions,
    descriptive_gene_effects,
    mission_effect_from_score,
)
from .data import MissionData, cpm_eligible_genes
from .statistics import (
    EPS,
    AxisScoreResult,
    combine_fixed_effects,
    random_effects_reml_mkh,
    score_signed_axis,
)


TARGET_ROLE = "podocyte_sensitivity_target"
MINIMUM_GENES = 8
RATIO_TOLERANCE = 1e-12
DECOMPOSITION_TOLERANCE = 1e-10
LEAVE_TOP_K_RETENTION = 0.5  # definition of k50, not a gate
DEFAULT_TOP_KS: tuple[int | str, ...] = (1, 2, 5, 10, 20, 40, "10%", "25%")

# (label, outcome, adjuster) exactly as the K0 script names its models; "diff" is
# podocyte minus structural, and adjuster None means unadjusted.
K0_MODEL_SPECS: tuple[tuple[str, str, str | None], ...] = (
    ("P0_podocyte_unadjusted", "podocyte", None),
    ("P1_podocyte_adjusted", "podocyte", "structural"),
    ("S0_structural_unadjusted", "structural", None),
    ("S1_structural_adjusted", "structural", "podocyte"),
    ("D_direct_difference", "diff", None),
)
K0_MODEL_CODES: dict[str, str] = {
    label.split("_", 1)[0]: label for label, _, _ in K0_MODEL_SPECS
}

PASS, FAIL, NOT_EVALUABLE = "PASS", "FAIL", "NOT_EVALUABLE"


class PreregistrationError(RuntimeError):
    """The preregistration file is missing, uncommitted, inconsistent or unparseable."""


# --------------------------------------------------------------------------------------
# Families and ratios
# --------------------------------------------------------------------------------------


def rank_positive_contributors(contributions: pd.Series) -> pd.Series:
    """Positive contributions, largest first; exact ties broken by gene name (deterministic)."""

    values = contributions.astype(float)
    positive = [(str(gene), float(value)) for gene, value in values.items() if value > 0]
    positive.sort(key=lambda item: (-item[1], item[0]))
    return pd.Series(
        [value for _, value in positive], index=[gene for gene, _ in positive], dtype=float
    )


def single_set_family(
    genes: Sequence[str],
    name: str,
    minimum: int = MINIMUM_GENES,
    *,
    role: str = TARGET_ROLE,
) -> dict[str, dict[str, Any]]:
    """One equal-weight, positively signed atlas-marker set in the frozen family shape."""

    unique = list(dict.fromkeys(str(gene) for gene in genes))
    if len(unique) < minimum:
        raise ValueError(f"{name} defines {len(unique)} genes; need {minimum}")
    return {
        str(name): {
            "role": role,
            "subdomains": {
                "atlas_markers": {
                    "genes": {gene: 1 for gene in unique},
                    "minimum_present": int(minimum),
                }
            },
        }
    }


def family_genes(family: Mapping[str, Any], axis: str) -> list[str]:
    """Defined genes of one axis in definition order."""

    return list(axis_directions(family[axis]))


def retention_ratio(
    estimate: float, reference: float, *, tolerance: float = RATIO_TOLERANCE
) -> float:
    """Signed estimate/reference; NaN when the reference is zero or either value is undefined."""

    try:
        estimate = float(estimate)
        reference = float(reference)
    except (TypeError, ValueError):
        return float("nan")
    if not (math.isfinite(estimate) and math.isfinite(reference)):
        return float("nan")
    if abs(reference) <= tolerance:
        return float("nan")
    return estimate / reference


def direction_retained(
    estimate: float, reference: float, *, tolerance: float = RATIO_TOLERANCE
) -> bool | None:
    """Whether the estimate keeps the reference sign; None when the reference has no sign."""

    try:
        estimate = float(estimate)
        reference = float(reference)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(estimate) and math.isfinite(reference)):
        return None
    if abs(reference) <= tolerance:
        return None
    return bool(np.sign(estimate) == np.sign(reference))


def ci_excludes_zero(ci_low: float, ci_high: float) -> bool | None:
    if not (np.isfinite(ci_low) and np.isfinite(ci_high)):
        return None
    return bool(ci_low > 0.0 or ci_high < 0.0)


def common_observable_genes(
    missions: Mapping[str, MissionData], threshold: float
) -> set[str]:
    """Genes CPM-eligible in every mission and present in every expression matrix (K4's rule)."""

    eligible = set.intersection(
        *(set(cpm_eligible_genes(data, threshold)) for data in missions.values())
    )
    observed = set.intersection(
        *(set(map(str, data.expression.index)) for data in missions.values())
    )
    return eligible & observed


def common_observable_family(
    missions: Mapping[str, MissionData],
    family: Mapping[str, Any],
    threshold: float,
    *,
    minimum: int = MINIMUM_GENES,
) -> dict[str, Any]:
    """Item 2: the frozen ``common_intersection_family`` after restricting to genes observed in
    every expression matrix, so its gene universe equals the strict-matching (K4) universe."""

    observed = set.intersection(
        *(set(map(str, data.expression.index)) for data in missions.values())
    )
    restricted = deepcopy(dict(family))
    for axis_spec in restricted.values():
        for subdomain in axis_spec["subdomains"].values():
            subdomain["genes"] = {
                gene: sign
                for gene, sign in subdomain["genes"].items()
                if str(gene) in observed
            }
    common = common_intersection_family(missions, restricted, threshold)
    for axis, axis_spec in common.items():
        for name, subdomain in axis_spec["subdomains"].items():
            if len(subdomain["genes"]) < minimum:
                raise ValueError(
                    f"{axis}/{name}: {len(subdomain['genes'])} common genes; need {minimum}"
                )
    return common


@dataclass(frozen=True)
class VariantMeta:
    """Observed REML/mKH synthesis of one family under one input variant."""

    meta: pd.DataFrame
    mission_effects: pd.DataFrame
    coverage: pd.DataFrame


def variant_meta(
    missions_variant: Mapping[str, MissionData],
    family: Mapping[str, Any],
    *,
    threshold: float,
    method: str = "mean",
) -> VariantMeta:
    """``run_sensitivities._meta_for_family`` plus per-axis counts of positive/negative missions."""

    meta, mission_effects, coverage = _meta_for_family(
        missions_variant, family, threshold, method=method
    )
    grouped = mission_effects.groupby("axis", sort=False)["estimate"]
    positive = grouped.apply(lambda values: int((values > 0).sum()))
    negative = grouped.apply(lambda values: int((values < 0).sum()))
    meta = meta.copy()
    meta["n_positive_missions"] = meta["axis"].map(positive).astype(int)
    meta["n_negative_missions"] = meta["axis"].map(negative).astype(int)
    return VariantMeta(meta=meta, mission_effects=mission_effects, coverage=coverage)


def genes_used_by_mission(coverage: pd.DataFrame, axis: str) -> dict[str, list[str]]:
    """Genes that passed observability per mission, from an ``eligible_axis_spec`` audit."""

    out: dict[str, list[str]] = {}
    for mission, sub in coverage[coverage["axis"] == axis].groupby("mission", sort=False):
        genes: list[str] = []
        for value in sub["genes_used"].fillna(""):
            genes.extend(gene for gene in str(value).split("|") if gene)
        out[str(mission)] = genes
    return out


# --------------------------------------------------------------------------------------
# K0 observed models
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class K0Observed:
    """Observed K0 models: HC3 flight coefficient per mission pooled by REML/mKH."""

    meta: pd.DataFrame
    mission_effects: pd.DataFrame
    coverage: pd.DataFrame


def k0_mission_fit(
    metadata: pd.DataFrame,
    outcome: pd.Series,
    adjuster: pd.Series | None,
    extra: pd.DataFrame | None = None,
) -> dict[str, float]:
    """HC3 flight coefficient for const + flight + block dummies + adjuster + extra columns.

    Without ``extra`` the design and fit are exactly K0's ``_regression_rows`` HC3 row.
    """

    design = _design(metadata, adjuster)
    extra_columns: list[str] = []
    if extra is not None and not extra.empty:
        missing = [animal for animal in metadata.index if animal not in extra.index]
        if missing:
            raise ValueError(f"extra covariates miss animals {missing[:5]}")
        clash = set(map(str, extra.columns)) & set(map(str, design.columns))
        if clash:
            raise ValueError(f"extra covariate names collide with the design: {sorted(clash)}")
        aligned = extra.loc[metadata.index].astype(float)
        if not np.all(np.isfinite(aligned.to_numpy())):
            raise ValueError("extra covariates must be finite")
        for column in aligned.columns:
            design[str(column)] = aligned[column].to_numpy()
            extra_columns.append(str(column))
    matrix = design.to_numpy(dtype=float)
    if np.linalg.matrix_rank(matrix) < matrix.shape[1]:
        raise ValueError("K0 design is rank deficient")
    fit = sm.OLS(outcome, design).fit()
    robust = fit.get_robustcov_results(cov_type="HC3")
    index = list(design.columns).index("flight")
    estimate = float(robust.params[index])
    standard_error = float(robust.bse[index])
    if not np.isfinite(standard_error) or standard_error <= 0.0:
        raise ValueError("HC3 variance is undefined (leverage-one observation?)")
    return {
        "coef_score_units": estimate,
        "standard_error": standard_error,
        "variance": standard_error**2,
        "p": float(robust.pvalues[index]),
        "n": int(len(outcome)),
        "n_parameters": int(matrix.shape[1]),
        "extra_covariates": "|".join(extra_columns),
    }


def k0_program_scores(
    data: MissionData,
    pod_genes: Sequence[str],
    struc_genes: Sequence[str],
    *,
    threshold: float,
    method: str = "mean",
    minimum: int = MINIMUM_GENES,
) -> tuple[dict[str, pd.Series], list[str], list[str]]:
    """K0's per-mission podocyte, structural and difference scores from eligible genes."""

    eligible = set(cpm_eligible_genes(data, threshold))
    pod_used = [g for g in pod_genes if g in eligible and g in data.expression.index]
    struc_used = [g for g in struc_genes if g in eligible and g in data.expression.index]
    if len(pod_used) < minimum or len(struc_used) < minimum:
        raise RuntimeError(
            f"{data.mission}: podocyte={len(pod_used)}, structural={len(struc_used)}"
        )
    podocyte = score_signed_axis(
        data.expression, {g: 1 for g in pod_used}, method=method
    ).scores
    structural = score_signed_axis(
        data.expression, {g: 1 for g in struc_used}, method=method
    ).scores
    if podocyte.isna().any() or structural.isna().any():
        raise ValueError(f"{data.mission}: nonfinite K0 program score")
    diff = (podocyte - structural).rename("signed_axis_score")
    return (
        {"podocyte": podocyte, "structural": structural, "diff": diff},
        pod_used,
        struc_used,
    )


def k0_observed(
    missions: Mapping[str, MissionData],
    pod_genes: Sequence[str],
    struc_disjoint_genes: Sequence[str],
    *,
    threshold: float,
    method: str = "mean",
    extra_covariates: Mapping[str, pd.DataFrame] | None = None,
    models: Sequence[tuple[str, str, str | None]] = K0_MODEL_SPECS,
    minimum: int = MINIMUM_GENES,
) -> K0Observed:
    """Observed K0 models, HC3 per mission then REML/mKH; no permutation family.

    ``extra_covariates`` maps a mission to animal-indexed columns appended to every model's
    design in that mission. With no extras and ``method="mean"`` the per-mission rows equal
    K0's ``_regression_rows`` HC3 rows.
    """

    extras = dict(extra_covariates or {})
    unknown = set(extras) - set(missions)
    if unknown:
        raise ValueError(f"extra covariates for unknown missions: {sorted(unknown)}")
    rows: list[dict[str, Any]] = []
    coverage: list[dict[str, Any]] = []
    for mission, data in missions.items():
        scores, pod_used, struc_used = k0_program_scores(
            data,
            pod_genes,
            struc_disjoint_genes,
            threshold=threshold,
            method=method,
            minimum=minimum,
        )
        coverage.append(
            {
                "mission": mission,
                "n_podocyte_genes": len(pod_used),
                "n_structural_disjoint_genes": len(struc_used),
                "podocyte_structural_score_pearson_r": float(
                    scores["podocyte"].corr(scores["structural"])
                ),
            }
        )
        for label, outcome_name, adjuster_name in models:
            fit = k0_mission_fit(
                data.metadata,
                scores[outcome_name],
                scores[adjuster_name] if adjuster_name else None,
                extras.get(mission),
            )
            rows.append(
                {
                    "mission": mission,
                    "model": label,
                    "variance_type": "HC3",
                    "scoring_method": method,
                    **fit,
                }
            )
    effects = pd.DataFrame(rows)
    order = list(missions)
    meta_rows = []
    for label, _, _ in models:
        sub = effects[effects["model"] == label].set_index("mission").loc[order]
        fit = random_effects_reml_mkh(
            sub["coef_score_units"].to_numpy(), sub["variance"].to_numpy()
        )
        meta_rows.append(
            {
                "model": label,
                "coef_score_units": fit.estimate,
                "standard_error_mkh": fit.standard_error,
                "t_mkh": fit.t,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh": fit.p,
                "tau2": fit.tau2,
                "i_squared": fit.i_squared,
                "maximum_weight": float(np.max(fit.weights)),
                "k_missions": fit.k,
                "n_positive_missions": int((sub["coef_score_units"] > 0).sum()),
                "scoring_method": method,
            }
        )
    return K0Observed(
        meta=pd.DataFrame(meta_rows),
        mission_effects=effects,
        coverage=pd.DataFrame(coverage),
    )


# --------------------------------------------------------------------------------------
# Exact signed contribution decomposition
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class SignedContributions:
    """Exact additive split of each pooled estimate into per-gene signed contributions."""

    genes: pd.DataFrame
    concentration: pd.DataFrame
    strata: pd.DataFrame


def _linear_gene_weights(
    axis_spec: Mapping[str, Any], result: AxisScoreResult
) -> pd.Series:
    """Weights w_g with score = sum_g w_g z_g for the equal-weight mean score."""

    used = [str(gene) for gene in result.signed_gene_z.index]
    if result.subdomain_scores is None:
        return pd.Series(1.0 / len(used), index=used, dtype=float)
    used_set = set(used)
    subdomains = axis_spec["subdomains"]
    n_domains = len(subdomains)
    weights: dict[str, float] = {}
    for name, subdomain in subdomains.items():
        members = [str(g) for g in subdomain["genes"] if str(g) in used_set]
        if not members:
            raise ValueError(f"subdomain {name!r} has no scored genes")
        for gene in members:
            if gene in weights:
                raise ValueError(f"gene {gene} occurs in two subdomains")
            weights[gene] = 1.0 / (n_domains * len(members))
    if set(weights) != used_set:
        raise ValueError("scored genes are not partitioned by the axis subdomains")
    return pd.Series(weights, dtype=float).loc[used]


def _mission_contributions(
    axis_spec: Mapping[str, Any],
    data: MissionData,
    result: AxisScoreResult,
) -> tuple[pd.Series, dict[str, float], list[dict[str, Any]]]:
    z = result.signed_gene_z.astype(float)
    if not np.all(np.isfinite(z.to_numpy())):
        raise ValueError(f"{data.mission}: signed gene z must be finite for decomposition")
    weights = _linear_gene_weights(axis_spec, result)
    score = result.scores.astype(float)
    reconstructed = z.mul(weights, axis=0).sum(axis=0)
    scale = max(1.0, float(np.max(np.abs(score.to_numpy()))))
    if float(np.max(np.abs(reconstructed - score))) > 1e-9 * scale:
        raise ValueError(
            "the axis score is not the equal-weight mean of signed gene z; "
            "the exact decomposition is defined for the mean method only"
        )
    summary, strata = mission_effect_from_score(score, data.metadata)
    fixed = combine_fixed_effects(strata["estimate"], strata["variance"])
    contribution = pd.Series(0.0, index=weights.index)
    audit: list[dict[str, Any]] = []
    blocks = list(data.metadata.groupby("block", sort=False))
    if len(blocks) != len(strata):
        raise RuntimeError("stratum bookkeeping mismatch")
    for position, (block, sub) in enumerate(blocks):
        flight = sub.index[sub["condition"] == "FLT"]
        control = sub.index[sub["condition"] == "GC"]
        nt, nc = len(flight), len(control)
        df = nt + nc - 2
        pooled_variance = (
            (nt - 1) * np.var(score.loc[flight].to_numpy(), ddof=1)
            + (nc - 1) * np.var(score.loc[control].to_numpy(), ddof=1)
        ) / df
        if pooled_variance <= EPS:
            raise ValueError(f"{data.mission}/{block}: zero pooled score variance")
        correction = 1.0 - 3.0 / (4.0 * df - 1.0)
        difference = z.loc[:, flight].mean(axis=1) - z.loc[:, control].mean(axis=1)
        stratum_part = correction * weights * difference / np.sqrt(pooled_variance)
        stratum_g = float(strata["estimate"].iloc[position])
        residual = float(stratum_part.sum()) - stratum_g
        if abs(residual) > DECOMPOSITION_TOLERANCE:
            raise RuntimeError(
                f"{data.mission}/{block}: stratum decomposition residual {residual:.3e}"
            )
        weight = float(fixed.weights[position])
        contribution = contribution + weight * stratum_part
        audit.append(
            {
                "mission": data.mission,
                "stratum": str(block),
                "hedges_g": stratum_g,
                "sum_gene_contributions": float(stratum_part.sum()),
                "fixed_effect_weight": weight,
                "residual": residual,
            }
        )
    return contribution, summary, audit


def absolute_shares(gene_effects: pd.DataFrame) -> pd.DataFrame:
    """``run_sensitivities.gene_contributions`` definition (|pooled gene g| / sum |pooled gene g|),
    extended so genes observed in fewer than two missions get NaN instead of an error."""

    rows = []
    for (axis, gene), sub in gene_effects.groupby(["axis", "gene"], sort=False):
        if sub["mission"].nunique() < 2:
            rows.append(
                {
                    "axis": axis,
                    "gene": gene,
                    "pooled_signed_g": np.nan,
                    "ci_low_mkh": np.nan,
                    "ci_high_mkh": np.nan,
                    "p_mkh_descriptive": np.nan,
                }
            )
            continue
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        rows.append(
            {
                "axis": axis,
                "gene": gene,
                "pooled_signed_g": fit.estimate,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh_descriptive": fit.p,
            }
        )
    table = pd.DataFrame(rows)
    table["absolute_share"] = table.groupby("axis")["pooled_signed_g"].transform(
        lambda values: values.abs() / values.abs().sum()
    )
    return table


def signed_contributions(
    missions: Mapping[str, MissionData],
    family: Mapping[str, Any],
    details: Mapping[str, Mapping[str, AxisScoreResult]],
) -> SignedContributions:
    """Exact additive decomposition of each pooled estimate into signed gene contributions.

    With S = sum_g w_g z_g (w_g = 1/m for one subdomain), the stratum effect
    g_s = J_s * Delta_s / SD_s splits exactly as c_gs = J_s * w_g * d_gs / SD_s, holding the score's
    pooled SD_s fixed (d_gs = flight minus ground mean of z_g). The mission effect uses the frozen
    fixed-effect stratum weights and the pooled estimate the fitted REML weights, both held fixed,
    so sum_g c_g equals the pooled estimate (checked to 1e-10). ``details`` is the fourth output
    of ``combined_score_design`` for the same missions and family.
    """

    mission_order = list(missions)
    gene_effects = descriptive_gene_effects(missions, details)
    shares = absolute_shares(gene_effects).set_index(["axis", "gene"])
    gene_rows: list[dict[str, Any]] = []
    concentration_rows: list[dict[str, Any]] = []
    strata_rows: list[dict[str, Any]] = []
    for axis, axis_spec in family.items():
        per_mission: dict[str, pd.Series] = {}
        estimates: list[float] = []
        variances: list[float] = []
        for mission in mission_order:
            contribution, summary, audit = _mission_contributions(
                axis_spec, missions[mission], details[mission][axis]
            )
            residual = float(contribution.sum()) - float(summary["estimate"])
            if abs(residual) > DECOMPOSITION_TOLERANCE:
                raise RuntimeError(f"{mission}/{axis}: mission residual {residual:.3e}")
            per_mission[mission] = contribution
            estimates.append(float(summary["estimate"]))
            variances.append(float(summary["variance"]))
            strata_rows.extend({"axis": axis, **row} for row in audit)
        fit = random_effects_reml_mkh(estimates, variances)
        defined = [g for g in axis_directions(axis_spec)]
        frame = pd.DataFrame(per_mission).reindex(
            [g for g in defined if any(g in s.index for s in per_mission.values())]
        )[mission_order]
        total = frame.fillna(0.0).mul(fit.weights, axis=1).sum(axis=1)
        pooled = float(fit.estimate)
        residual = float(total.sum()) - pooled
        if abs(residual) > DECOMPOSITION_TOLERANCE:
            raise RuntimeError(f"{axis}: pooled decomposition residual {residual:.3e}")
        has_sign = abs(pooled) > RATIO_TOLERANCE
        signed_share = total / pooled if has_sign else total * np.nan
        n_missions = frame.notna().sum(axis=1)
        ranks = total.rank(ascending=False, method="first")
        for gene in frame.index:
            key = (axis, gene)
            row = {
                "target": axis,
                "gene": gene,
                "n_missions": int(n_missions[gene]),
                "pooled_signed_g": float(shares.loc[key, "pooled_signed_g"])
                if key in shares.index
                else np.nan,
                "absolute_share": float(shares.loc[key, "absolute_share"])
                if key in shares.index
                else np.nan,
                "signed_contribution": float(total[gene]),
                "signed_share": float(signed_share[gene]),
                "rank_signed_contribution": int(ranks[gene]),
            }
            for mission in mission_order:
                row[f"contribution_{mission}"] = float(frame.loc[gene, mission])
            gene_rows.append(row)

        positive = rank_positive_contributors(total)
        negative = total[total < 0]
        sum_positive = float(positive.sum())
        concentration = {
            "target": axis,
            "pooled_estimate": pooled,
            "sum_signed_contributions": float(total.sum()),
            "decomposition_residual": residual,
            "n_genes_contributing": int(len(total)),
            "n_positive_contributors": int(len(positive)),
            "n_negative_contributors": int(len(negative)),
            "sum_positive_contributions": sum_positive,
            "sum_negative_contributions": float(negative.sum()),
            "n_eff_positive": (sum_positive**2 / float((positive**2).sum()))
            if len(positive)
            else np.nan,
            "top1_positive_gene": str(positive.index[0]) if len(positive) else "",
            "top1_positive_share": float(positive.iloc[0] / sum_positive)
            if len(positive)
            else np.nan,
            "top10_positive_share": float(positive.iloc[:10].sum() / sum_positive)
            if len(positive)
            else np.nan,
            "max_signed_share": float(signed_share.max()) if has_sign else np.nan,
            "max_signed_share_gene": str(signed_share.idxmax()) if has_sign else "",
            "min_signed_share": float(signed_share.min()) if has_sign else np.nan,
            "max_absolute_share": np.nan,
            "max_absolute_share_gene": "",
        }
        axis_shares = shares.xs(axis, level="axis")["absolute_share"].dropna()
        if len(axis_shares):
            concentration["max_absolute_share"] = float(axis_shares.max())
            concentration["max_absolute_share_gene"] = str(axis_shares.idxmax())
        for mission, weight, estimate in zip(mission_order, fit.weights, estimates):
            concentration[f"random_effect_weight_{mission}"] = float(weight)
            concentration[f"mission_estimate_{mission}"] = float(estimate)
        concentration_rows.append(concentration)
    return SignedContributions(
        genes=pd.DataFrame(gene_rows),
        concentration=pd.DataFrame(concentration_rows),
        strata=pd.DataFrame(strata_rows),
    )


# --------------------------------------------------------------------------------------
# Leave-top-k influence curve (descriptive)
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class LeaveTopK:
    """Descriptive influence curve after dropping the top-k positive contributors."""

    curve: pd.DataFrame
    k50: int | None
    k_ci0: int | None


def resolve_top_ks(ks: Sequence[int | str], n_genes: int) -> list[tuple[str, int]]:
    """Resolve integer and percentage k; a percentage is ceil(fraction x genes contributing)."""

    resolved: list[tuple[str, int]] = []
    for k in ks:
        text = str(k).strip()
        if text.endswith("%"):
            fraction = float(text[:-1]) / 100.0
            if not 0.0 < fraction <= 1.0:
                raise ValueError(f"percentage k out of range: {k!r}")
            value = int(math.ceil(fraction * n_genes - 1e-9))
        else:
            value = int(text)
        if value < 1:
            raise ValueError(f"k must be at least one: {k!r}")
        resolved.append((text, value))
    return sorted(resolved, key=lambda item: (item[1], item[0]))


def leave_top_k_curve(
    missions: Mapping[str, MissionData],
    genes: Sequence[str],
    contributions: pd.Series,
    *,
    threshold: float,
    ks: Sequence[int | str] = DEFAULT_TOP_KS,
    name: str = "target",
    minimum: int = MINIMUM_GENES,
) -> LeaveTopK:
    """Drop the top-k positive signed contributors, rescore and re-pool (descriptive only).

    Removal sets are nested in k by construction (one fixed ranking, ties broken by gene name).
    ``k50`` is the first k whose retention falls below 0.5 and ``k_ci0`` the first k (including
    the full set, k = 0) whose interval includes zero.
    """

    contributions = contributions.astype(float)
    ranked = list(rank_positive_contributors(contributions).index)
    defined = list(dict.fromkeys(str(gene) for gene in genes))
    rows: list[dict[str, Any]] = []
    cache: dict[int, dict[str, Any]] = {}

    def _evaluate(n_removed: int) -> dict[str, Any]:
        if n_removed in cache:
            return cache[n_removed]
        removed = set(ranked[:n_removed])
        remaining = [gene for gene in defined if gene not in removed]
        out: dict[str, Any] = {"n_genes_defined_remaining": len(remaining)}
        try:
            result = variant_meta(
                missions,
                single_set_family(remaining, name, minimum),
                threshold=threshold,
            )
        except ValueError as error:
            out.update({"evaluable": False, "reason": str(error)})
        else:
            meta = result.meta.iloc[0]
            out.update(
                {
                    "evaluable": True,
                    "reason": "",
                    "estimate": float(meta["estimate"]),
                    "ci_low_mkh": float(meta["ci_low_mkh"]),
                    "ci_high_mkh": float(meta["ci_high_mkh"]),
                    "p_mkh": float(meta["p_mkh"]),
                    "tau2": float(meta["tau2"]),
                    "i_squared": float(meta["i_squared"]),
                    "n_positive_missions": int(meta["n_positive_missions"]),
                }
            )
        cache[n_removed] = out
        return out

    baseline = _evaluate(0)
    if not baseline.get("evaluable", False):
        raise ValueError(f"{name}: full set is not evaluable ({baseline['reason']})")
    reference = baseline["estimate"]
    plan = [("0", 0)] + resolve_top_ks(ks, len(contributions))
    for label, k in plan:
        n_removed = min(k, len(ranked))
        result = _evaluate(n_removed)
        estimate = result.get("estimate", np.nan)
        low = result.get("ci_low_mkh", np.nan)
        high = result.get("ci_high_mkh", np.nan)
        rows.append(
            {
                "target": name,
                "k_label": label,
                "k_requested": k,
                "n_removed": n_removed,
                "capped_at_n_positive": bool(k > len(ranked)),
                "removed_genes": "|".join(ranked[:n_removed]),
                **{
                    key: result.get(key, np.nan)
                    for key in (
                        "n_genes_defined_remaining",
                        "evaluable",
                        "reason",
                        "estimate",
                        "ci_low_mkh",
                        "ci_high_mkh",
                        "p_mkh",
                        "tau2",
                        "i_squared",
                        "n_positive_missions",
                    )
                },
                "reference_estimate": reference,
                "retention_ratio": retention_ratio(estimate, reference),
                "ci_includes_zero": bool(low <= 0.0 <= high)
                if np.isfinite(low) and np.isfinite(high)
                else None,
                "status": "descriptive_influence_diagnostic",
            }
        )
    curve = pd.DataFrame(rows)
    evaluable = curve[curve["evaluable"].astype(bool)].sort_values(
        ["n_removed", "k_requested"], kind="mergesort"
    )
    below = evaluable[
        (evaluable["n_removed"] > 0)
        & (evaluable["retention_ratio"] < LEAVE_TOP_K_RETENTION)
    ]
    includes = evaluable[[value is True for value in evaluable["ci_includes_zero"]]]
    k50 = int(below["n_removed"].iloc[0]) if len(below) else None
    k_ci0 = int(includes["n_removed"].iloc[0]) if len(includes) else None
    return LeaveTopK(curve=curve, k50=k50, k_ci0=k_ci0)


# --------------------------------------------------------------------------------------
# Leave-one-gene annotation
# --------------------------------------------------------------------------------------


def annotate_leave_one_gene(
    logo: pd.DataFrame,
    reference: Mapping[str, float],
    used_by_mission: Mapping[str, Mapping[str, Sequence[str]]],
) -> pd.DataFrame:
    """Add signed retention and gene-usage columns to ``run_sensitivities.leave_one_gene``.

    Genes never eligible in any mission are reported but not counted by the item-4 gate,
    because omitting them cannot change any estimate.
    """

    out = logo.copy().rename(columns={"axis": "target"})
    out["reference_estimate"] = out["target"].map(reference)
    out["retention_ratio"] = [
        retention_ratio(estimate, ref) if bool(ok) else np.nan
        for estimate, ref, ok in zip(
            out.get("estimate", pd.Series(np.nan, index=out.index)),
            out["reference_estimate"],
            out["evaluable"],
        )
    ]
    counts = []
    for target, gene in zip(out["target"], out["omitted_gene"]):
        used = used_by_mission.get(target, {})
        counts.append(sum(gene in set(genes) for genes in used.values()))
    out["n_missions_used"] = counts
    out["gene_used_in_any_mission"] = out["n_missions_used"] > 0
    out["counted_in_gate"] = out["gene_used_in_any_mission"] & out["evaluable"].astype(bool)
    return out


# --------------------------------------------------------------------------------------
# Item 7 helpers (OSD-163 mapping rate)
# --------------------------------------------------------------------------------------


def standardized_qc_covariate(data: MissionData, metric: str) -> pd.Series:
    """Label-blind z of one per-animal QC metric (ddof 1), in metadata order."""

    if metric not in data.qc:
        raise RuntimeError(f"{data.mission}: QC metric {metric} unavailable")
    values = data.qc[metric].reindex(data.metadata.index).astype(float)
    if values.isna().any() or values.std(ddof=1) <= 0:
        raise RuntimeError(f"{data.mission}: QC metric {metric} unusable")
    return (values - values.mean()) / values.std(ddof=1)


def residualized_mission_effect(
    score: pd.Series, metadata: pd.DataFrame, covariate: pd.Series
) -> dict[str, Any]:
    """Hedges g after label-blind OLS residualization of the score on [1, covariate]."""

    covariate = covariate.loc[metadata.index].astype(float)
    outcome = score.loc[metadata.index].astype(float)
    matrix = np.column_stack([np.ones(len(covariate)), covariate.to_numpy()])
    coefficients = np.linalg.lstsq(matrix, outcome.to_numpy(), rcond=None)[0]
    residual = pd.Series(outcome.to_numpy() - matrix @ coefficients, index=outcome.index)
    summary, _ = mission_effect_from_score(residual, metadata)
    return summary


def ancova_standardized_effect(
    score: pd.Series, metadata: pd.DataFrame, covariate: pd.Series
) -> dict[str, float]:
    """Covariate-adjusted standardized flight effect (descriptive, item 7 secondary).

    OLS of score on const + flight + block dummies + covariate. The flight coefficient is divided
    by the unadjusted pooled within-(block x group) SD of the score and multiplied by
    J = 1 - 3/(4 df - 1); its variance is (nt + nc)/(nt nc) + g^2/(2 df) with
    df = n - rank(design) (n - 3 for a single block).
    """

    covariate = covariate.loc[metadata.index].astype(float)
    outcome = score.loc[metadata.index].astype(float)
    design = pd.DataFrame(
        {"flight": (metadata["condition"] == "FLT").astype(float)}, index=metadata.index
    )
    design = design.join(pd.get_dummies(metadata["block"], drop_first=True, dtype=float))
    design["covariate"] = covariate.to_numpy()
    design = sm.add_constant(design, has_constant="add")
    matrix = design.to_numpy(dtype=float)
    rank = int(np.linalg.matrix_rank(matrix))
    if rank < matrix.shape[1]:
        raise ValueError("ANCOVA design is rank deficient")
    fit = sm.OLS(outcome, design).fit()
    beta = float(fit.params["flight"])
    n = len(outcome)
    df = n - rank
    if df < 2:
        raise ValueError("ANCOVA leaves fewer than two residual degrees of freedom")
    numerator = 0.0
    denominator = 0
    for _, sub in metadata.groupby(["block", "condition"], sort=False):
        values = outcome.loc[sub.index].to_numpy()
        if len(values) >= 2:
            numerator += (len(values) - 1) * float(np.var(values, ddof=1))
            denominator += len(values) - 1
    if denominator <= 0 or numerator <= EPS:
        raise ValueError("pooled within-group SD is undefined")
    sd = math.sqrt(numerator / denominator)
    nt = int((metadata["condition"] == "FLT").sum())
    nc = int((metadata["condition"] == "GC").sum())
    correction = 1.0 - 3.0 / (4.0 * df - 1.0)
    estimate = correction * beta / sd
    variance = (nt + nc) / (nt * nc) + estimate**2 / (2.0 * df)
    return {
        "estimate": float(estimate),
        "variance": float(variance),
        "beta_flight_score_units": beta,
        "pooled_within_sd": sd,
        "df": float(df),
        "n_flight": float(nt),
        "n_control": float(nc),
    }


# --------------------------------------------------------------------------------------
# Gate rules read from the preregistration
# --------------------------------------------------------------------------------------

_OPERATORS = {">": operator.gt, ">=": operator.ge, "<": operator.lt, "<=": operator.le}
_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)"
_OP = r"(>=|<=|>|<)"
_CLAUSE_PATTERNS = (
    (
        "threshold",
        re.compile(
            rf"(?:(?:mean|median|common)\s+)?(?:pooled\s+)?estimate\s*{_OP}\s*({_NUMBER})",
            re.IGNORECASE,
        ),
    ),
    (
        "all",
        re.compile(
            rf"all\s+(\d+)\s+lomo\s+estimates\s*{_OP}\s*({_NUMBER})", re.IGNORECASE
        ),
    ),
    (
        "sign",
        re.compile(
            r"estimate\s+has\s+the\s+(?:same\s+)?sign\s+(?:of|as)\s+the\s+reference",
            re.IGNORECASE,
        ),
    ),
    (
        "ratio",
        re.compile(rf"estimate\s*/\s*reference\s*{_OP}\s*({_NUMBER})", re.IGNORECASE),
    ),
)


@dataclass(frozen=True)
class GateClause:
    kind: str
    text: str
    op: str | None = None
    threshold: float | None = None
    count: int | None = None


@dataclass(frozen=True)
class GateRule:
    text: str
    clauses: tuple[GateClause, ...]


def parse_gate_rule(
    text: str, *, aliases: Mapping[str, str] | None = None
) -> GateRule:
    """Parse 'clause AND clause ...'; an alias (e.g. gate_rule_default) resolves first."""

    raw = str(text).strip()
    if aliases and raw in aliases:
        raw = str(aliases[raw]).strip()
    normalized = " ".join(raw.split())
    if not normalized:
        raise PreregistrationError("empty gate rule")
    clauses = []
    for part in re.split(r"\s+AND\s+", normalized):
        part = part.strip().strip(".")
        for kind, pattern in _CLAUSE_PATTERNS:
            match = pattern.fullmatch(part)
            if match is None:
                continue
            if kind == "threshold":
                clauses.append(
                    GateClause(kind, part, op=match.group(1), threshold=float(match.group(2)))
                )
            elif kind == "all":
                clauses.append(
                    GateClause(
                        kind,
                        part,
                        op=match.group(2),
                        threshold=float(match.group(3)),
                        count=int(match.group(1)),
                    )
                )
            elif kind == "sign":
                clauses.append(GateClause(kind, part))
            else:
                clauses.append(
                    GateClause(kind, part, op=match.group(1), threshold=float(match.group(2)))
                )
            break
        else:
            raise PreregistrationError(f"unrecognised gate clause {part!r} in {raw!r}")
    return GateRule(text=normalized, clauses=tuple(clauses))


def verdict_label(passed: bool | None) -> str:
    if passed is None:
        return NOT_EVALUABLE
    return PASS if passed else FAIL


def combine_passes(values: Sequence[bool | None]) -> bool | None:
    """FAIL if any clause fails, PASS if all pass, otherwise not evaluable."""

    values = list(values)
    if any(value is False for value in values):
        return False
    if values and all(value is True for value in values):
        return True
    return None


def _finite(value: Any) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def evaluate_gate(
    rule: GateRule,
    *,
    estimate: float = float("nan"),
    reference: float = float("nan"),
    lomo_estimates: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Evaluate a parsed rule; undefined inputs give NOT_EVALUABLE, never PASS."""

    clause_rows = []
    for clause in rule.clauses:
        value: Any = None
        passed: bool | None
        detail = ""
        if clause.kind == "threshold":
            value = float(estimate) if _finite(estimate) else None
            passed = (
                None
                if value is None
                else bool(_OPERATORS[clause.op](value, clause.threshold))
            )
        elif clause.kind == "all":
            values = list(lomo_estimates) if lomo_estimates is not None else []
            value = [float(v) for v in values]
            if len(values) != clause.count:
                passed = None
                detail = f"expected {clause.count} estimates, found {len(values)}"
            elif not all(_finite(v) for v in values):
                passed = None
                detail = "undefined estimate"
            else:
                passed = bool(
                    all(_OPERATORS[clause.op](float(v), clause.threshold) for v in values)
                )
        elif clause.kind == "sign":
            passed = direction_retained(estimate, reference)
            value = passed
            if passed is None:
                detail = "reference sign undefined"
        elif clause.kind == "ratio":
            ratio = retention_ratio(estimate, reference)
            value = ratio if math.isfinite(ratio) else None
            passed = (
                None if value is None else bool(_OPERATORS[clause.op](ratio, clause.threshold))
            )
            if value is None:
                detail = "retention ratio undefined (reference zero or missing)"
        else:  # pragma: no cover - parse_gate_rule only emits the kinds above
            raise PreregistrationError(f"unknown clause kind {clause.kind}")
        clause_rows.append(
            {
                "clause": clause.text,
                "kind": clause.kind,
                "value": value,
                "threshold": clause.threshold,
                "passed": passed,
                "detail": detail,
            }
        )
    passed = combine_passes([row["passed"] for row in clause_rows])
    return {
        "rule": rule.text,
        "passed": passed,
        "verdict": verdict_label(passed),
        "clauses": clause_rows,
    }


ITEM4_GATE_KEYS = ("logo_all_direction", "logo_min_retention", "max_single_gene_signed_share")
ITEM4_SECONDARY_KEYS = ("top1_positive_share_max", "top10_positive_share_max")


def evaluate_gene_influence(
    gates: Mapping[str, Any],
    secondary: Mapping[str, Any] | None,
    *,
    logo: pd.DataFrame,
    concentration: Mapping[str, Any],
) -> dict[str, Any]:
    """Item 4: leave-one-gene and signed-share gates; secondary bounds are reported only.

    ``logo`` holds the counted leave-one-gene rows (``direction_retained``, ``retention_ratio``);
    ``concentration`` one row of ``SignedContributions.concentration``.
    """

    unknown = set(gates) - set(ITEM4_GATE_KEYS)
    if unknown:
        raise PreregistrationError(f"unknown item-4 gates: {sorted(unknown)}")
    unknown = set(secondary or {}) - set(ITEM4_SECONDARY_KEYS)
    if unknown:
        raise PreregistrationError(f"unknown item-4 secondary bounds: {sorted(unknown)}")
    clauses = []
    if "logo_all_direction" in gates:
        required = bool(gates["logo_all_direction"])
        directions = list(logo["direction_retained"]) if len(logo) else []
        if not required:
            passed: bool | None = True
            detail = "not required by the preregistration"
        elif not directions or any(value is None or pd.isna(value) for value in directions):
            passed, detail = None, "leave-one-gene direction undefined"
        else:
            passed, detail = bool(all(bool(v) for v in directions)), ""
        clauses.append(
            {
                "clause": "logo_all_direction",
                "value": int(sum(bool(v) for v in directions if v is not None and not pd.isna(v))),
                "threshold": len(directions),
                "passed": passed,
                "detail": detail,
            }
        )
    if "logo_min_retention" in gates:
        threshold = float(gates["logo_min_retention"])
        ratios = pd.to_numeric(logo["retention_ratio"], errors="coerce") if len(logo) else pd.Series(dtype=float)
        if not len(ratios) or ratios.isna().any():
            passed, value, detail = None, None, "retention ratio undefined"
        else:
            value = float(ratios.min())
            passed, detail = bool(value >= threshold), ""
        clauses.append(
            {
                "clause": "logo_min_retention",
                "value": value,
                "threshold": threshold,
                "passed": passed,
                "detail": detail,
            }
        )
    if "max_single_gene_signed_share" in gates:
        threshold = float(gates["max_single_gene_signed_share"])
        value = concentration.get("max_signed_share", np.nan)
        passed = bool(float(value) <= threshold) if _finite(value) else None
        clauses.append(
            {
                "clause": "max_single_gene_signed_share",
                "value": float(value) if _finite(value) else None,
                "threshold": threshold,
                "passed": passed,
                "detail": str(concentration.get("max_signed_share_gene", "")),
            }
        )
    secondary_rows = []
    for key, column in (
        ("top1_positive_share_max", "top1_positive_share"),
        ("top10_positive_share_max", "top10_positive_share"),
    ):
        if secondary and key in secondary:
            threshold = float(secondary[key])
            value = concentration.get(column, np.nan)
            secondary_rows.append(
                {
                    "bound": key,
                    "value": float(value) if _finite(value) else None,
                    "threshold": threshold,
                    "within_bound": bool(float(value) <= threshold) if _finite(value) else None,
                    "role": "secondary_v13_precedent_reported_not_gating",
                }
            )
    passed = combine_passes([row["passed"] for row in clauses])
    return {
        "passed": passed,
        "verdict": verdict_label(passed),
        "clauses": clauses,
        "secondary_bounds": secondary_rows,
    }


def classify_overall(
    item_verdicts: Mapping[int, str],
    required_items: Sequence[int],
    *,
    robust_label: str,
    fragile_label: str,
    incomplete_label: str = "INCOMPLETE",
) -> dict[str, Any]:
    """Robust only if every required item passes; fragile if any fails (named)."""

    failed = [item for item in required_items if item_verdicts.get(item) == FAIL]
    open_items = [
        item
        for item in required_items
        if item_verdicts.get(item) not in (PASS, FAIL)
    ]
    if failed:
        label = fragile_label
    elif not open_items:
        label = robust_label
    else:
        label = incomplete_label
    return {
        "classification": label,
        "failed_items": failed,
        "not_evaluable_items": open_items,
        "required_items": list(required_items),
    }


def estimate_verdict(estimate: float, ci_low: float, ci_high: float) -> str:
    """Coarse headline label used for gene-ID-map reconstruction comparisons."""

    if not (_finite(estimate) and _finite(ci_low) and _finite(ci_high)):
        return NOT_EVALUABLE
    if estimate <= 0:
        return "non_positive"
    if ci_low > 0:
        return "positive_ci_excludes_zero"
    return "positive_ci_includes_zero"


# Most conservative first.
VERDICT_ORDER: dict[str, int] = {
    "non_positive": 0,
    FAIL: 0,
    NOT_EVALUABLE: 1,
    "INCOMPLETE": 1,
    "positive_ci_includes_zero": 2,
    PASS: 3,
    "positive_ci_excludes_zero": 3,
}


def conservative_verdict(verdicts: Sequence[str], *, order: Mapping[str, int] | None = None) -> str:
    """The least favourable verdict; unknown labels rank as most conservative."""

    order = dict(VERDICT_ORDER if order is None else order)
    return min(verdicts, key=lambda value: (order.get(value, -1), str(value)))
