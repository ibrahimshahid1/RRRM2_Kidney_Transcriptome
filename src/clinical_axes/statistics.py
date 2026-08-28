"""Statistical primitives for cross-mission renal tissue-axis analyses.

This module deliberately contains no repository paths, cohort names, gene
panels, or biological interpretation.  It provides the label-blind scoring
and inferential operations used by the clinical-axis runner:

* signed gene-wise z scores, with optional equal-weight subdomains;
* Hedges' g within an exchangeability stratum;
* inverse-variance fixed-effect pooling of strata within a mission;
* REML random-effects pooling across missions with modified
  Hartung--Knapp uncertainty and a prediction interval; and
* whole-pipeline blocked label permutations with max-|T| family-wise error
  control across a frozen set of axes.

Positive scores and effects always mean movement in the direction encoded by
the supplied gene signs.  The caller, not this module, owns that biological
direction convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import optimize, stats


EPS = np.finfo(float).eps * 100.0


@dataclass(frozen=True)
class AxisScoreResult:
    """Per-animal axis scores and an audit of the genes contributing to them."""

    scores: pd.Series
    signed_gene_z: pd.DataFrame
    subdomain_scores: pd.DataFrame | None
    genes_requested: tuple[str, ...]
    genes_used: tuple[str, ...]
    genes_missing: tuple[str, ...]


@dataclass(frozen=True)
class HedgesGResult:
    """Bias-corrected standardized mean difference (treatment minus control)."""

    estimate: float
    variance: float
    standard_error: float
    ci_low: float
    ci_high: float
    n_treatment: int
    n_control: int


@dataclass(frozen=True)
class FixedEffectResult:
    """Inverse-variance fixed-effect synthesis, typically across strata."""

    estimate: float
    variance: float
    standard_error: float
    ci_low: float
    ci_high: float
    z: float
    p: float
    q: float
    q_p: float
    weights: np.ndarray
    k: int


@dataclass(frozen=True)
class RandomEffectsResult:
    """REML synthesis with modified Hartung--Knapp inference."""

    estimate: float
    tau2: float
    standard_error: float
    t: float
    df: int
    ci_low: float
    ci_high: float
    p: float
    prediction_low: float
    prediction_high: float
    prediction_df: int | None
    q: float
    q_p: float
    i_squared: float
    weights: np.ndarray
    k: int


@dataclass(frozen=True)
class MetaPermutationResult:
    """Observed cross-mission results and their blocked permutation null."""

    observed_meta: pd.DataFrame
    mission_effects: pd.DataFrame
    null_t: pd.DataFrame
    null_max_abs_t: np.ndarray
    empirical_p_two_sided: pd.Series
    max_t_fwer: pd.Series
    n_permutations: int
    seed: int


def _validate_expression(expression: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(expression, pd.DataFrame) or expression.empty:
        raise ValueError("expression must be a non-empty genes-by-samples DataFrame")
    if expression.index.has_duplicates:
        raise ValueError("expression gene index contains duplicates")
    if expression.columns.has_duplicates:
        raise ValueError("expression sample columns contain duplicates")
    return expression.astype(float)


def genewise_z_scores(
    expression: pd.DataFrame,
    genes: Sequence[str] | None = None,
    *,
    ddof: int = 1,
) -> pd.DataFrame:
    """Z-standardize each gene across the supplied analysis samples.

    Missing observations remain missing and are not converted to a neutral
    value.  A gene with at least two finite observations but zero variance is
    assigned zero for its finite observations.  A gene with fewer than two
    finite values cannot be standardized and is left entirely missing.
    """

    x = _validate_expression(expression)
    if genes is not None:
        requested = list(dict.fromkeys(genes))
        present = [gene for gene in requested if gene in x.index]
        if not present:
            raise ValueError("none of the requested genes is observable")
        x = x.loc[present]
    if ddof < 0:
        raise ValueError("ddof must be non-negative")

    values = x.to_numpy(dtype=float)
    finite = np.isfinite(values)
    n = finite.sum(axis=1)
    means = np.divide(
        np.where(finite, values, 0.0).sum(axis=1),
        n,
        out=np.full(len(x), np.nan, dtype=float),
        where=n > 0,
    )
    centered = values - means[:, None]
    centered[~finite] = np.nan
    ss = np.nansum(centered**2, axis=1)
    denom = n - ddof
    sd = np.sqrt(
        np.divide(
            ss,
            denom,
            out=np.full(len(x), np.nan, dtype=float),
            where=denom > 0,
        )
    )

    z = np.full_like(values, np.nan, dtype=float)
    variable = np.isfinite(sd) & (sd > EPS)
    z[variable] = centered[variable] / sd[variable, None]
    constant = (n > ddof) & np.isfinite(sd) & (sd <= EPS)
    z[constant] = np.where(finite[constant], 0.0, np.nan)
    return pd.DataFrame(z, index=x.index, columns=x.columns)


def _aggregate_rows(values: pd.DataFrame, method: str) -> pd.Series:
    if method == "mean":
        return values.mean(axis=0, skipna=True)
    if method == "median":
        return values.median(axis=0, skipna=True)
    raise ValueError("method must be 'mean' or 'median'")


def _validated_signs(directions: Mapping[str, float] | pd.Series) -> pd.Series:
    signs = pd.Series(dict(directions), dtype=float)
    if signs.empty:
        raise ValueError("directions must contain at least one gene")
    if signs.index.has_duplicates:
        raise ValueError("directions contains duplicate genes")
    if np.any(~np.isfinite(signs.to_numpy())) or np.any(
        ~np.isin(signs.to_numpy(), (-1.0, 1.0))
    ):
        raise ValueError("every gene direction must be exactly -1 or +1")
    return signs


def _validate_subdomains(
    subdomains: Mapping[str, Sequence[str]],
    requested_genes: Sequence[str],
) -> dict[str, tuple[str, ...]]:
    if not subdomains:
        raise ValueError("subdomains cannot be empty")
    normalized: dict[str, tuple[str, ...]] = {}
    seen: set[str] = set()
    for name, members in subdomains.items():
        if not str(name):
            raise ValueError("subdomain names must be non-empty")
        genes = tuple(dict.fromkeys(members))
        if not genes:
            raise ValueError(f"subdomain {name!r} contains no genes")
        overlap = seen.intersection(genes)
        if overlap:
            raise ValueError(
                "genes may occur in only one equal-weight subdomain: "
                + ", ".join(sorted(overlap))
            )
        seen.update(genes)
        normalized[str(name)] = genes

    requested = set(requested_genes)
    if seen != requested:
        missing = sorted(requested - seen)
        extra = sorted(seen - requested)
        detail = []
        if missing:
            detail.append("unassigned=" + ",".join(missing))
        if extra:
            detail.append("unknown=" + ",".join(extra))
        raise ValueError("subdomains must partition directions (" + "; ".join(detail) + ")")
    return normalized


def score_signed_axis(
    expression: pd.DataFrame,
    directions: Mapping[str, float] | pd.Series,
    *,
    method: str = "mean",
    subdomains: Mapping[str, Sequence[str]] | None = None,
    min_genes_per_subdomain: int = 1,
    require_all_genes: bool = False,
) -> AxisScoreResult:
    """Compute a signed, label-blind tissue-axis score for every sample.

    Genes are z-standardized separately across the samples in ``expression``
    before their frozen +/- direction is applied.  Without ``subdomains``, all
    observable genes receive equal weight.  With ``subdomains``, genes are
    first aggregated inside each subdomain and the subdomain scores are then
    averaged, ensuring (for example) that a large ECM list cannot outweigh a
    smaller maladaptive-repair list simply because it contains more genes.

    ``method='median'`` changes the within-domain gene aggregation; equal
    subdomain weighting remains an arithmetic mean across domains.
    """

    if min_genes_per_subdomain < 1:
        raise ValueError("min_genes_per_subdomain must be at least one")
    x = _validate_expression(expression)
    signs = _validated_signs(directions)
    requested = tuple(str(gene) for gene in signs.index)
    present = tuple(gene for gene in requested if gene in x.index)
    missing = tuple(gene for gene in requested if gene not in x.index)
    if require_all_genes and missing:
        raise ValueError("required genes are missing: " + ", ".join(missing))
    if not present:
        raise ValueError("none of the axis genes is observable")

    z = genewise_z_scores(x, present)
    signed = z.mul(signs.loc[list(present)], axis=0)

    domain_table: pd.DataFrame | None = None
    if subdomains is None:
        if len(present) < min_genes_per_subdomain:
            raise ValueError(
                f"only {len(present)} genes observable; "
                f"need {min_genes_per_subdomain}"
            )
        scores = _aggregate_rows(signed, method)
        enough = signed.notna().sum(axis=0) >= min_genes_per_subdomain
        scores = scores.where(enough)
    else:
        domains = _validate_subdomains(subdomains, requested)
        domain_scores: dict[str, pd.Series] = {}
        for name, genes in domains.items():
            observed = [gene for gene in genes if gene in present]
            if len(observed) < min_genes_per_subdomain:
                raise ValueError(
                    f"subdomain {name!r} has {len(observed)} observable genes; "
                    f"need {min_genes_per_subdomain}"
                )
            domain_score = _aggregate_rows(signed.loc[observed], method)
            enough = signed.loc[observed].notna().sum(axis=0) >= min_genes_per_subdomain
            domain_scores[name] = domain_score.where(enough)
        domain_table = pd.DataFrame(domain_scores, index=x.columns)
        # Requiring every subdomain to be finite preserves equal weighting.
        scores = domain_table.mean(axis=1, skipna=False)

    scores.name = "signed_axis_score"
    return AxisScoreResult(
        scores=scores,
        signed_gene_z=signed,
        subdomain_scores=domain_table,
        genes_requested=requested,
        genes_used=present,
        genes_missing=missing,
    )


def hedges_g(
    treatment: Sequence[float],
    control: Sequence[float],
    *,
    alpha: float = 0.05,
) -> HedgesGResult:
    """Return Hedges' g and its conventional sampling variance.

    The small-sample correction is ``J = 1 - 3/(4*df - 1)`` and the sampling
    variance is ``(n1+n0)/(n1*n0) + g^2/(2*df)``.  Non-finite observations are
    rejected rather than silently removed because exchangeability-block sample
    counts are part of the frozen design.
    """

    xt = np.asarray(treatment, dtype=float)
    xc = np.asarray(control, dtype=float)
    if xt.ndim != 1 or xc.ndim != 1:
        raise ValueError("treatment and control must be one-dimensional")
    if len(xt) < 2 or len(xc) < 2:
        raise ValueError("Hedges' g requires at least two observations per group")
    if np.any(~np.isfinite(xt)) or np.any(~np.isfinite(xc)):
        raise ValueError("Hedges' g inputs must be finite")
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie between zero and one")

    nt, nc = len(xt), len(xc)
    df = nt + nc - 2
    pooled_variance = (
        (nt - 1) * np.var(xt, ddof=1) + (nc - 1) * np.var(xc, ddof=1)
    ) / df
    difference = float(np.mean(xt) - np.mean(xc))
    if pooled_variance <= EPS:
        if abs(difference) > EPS:
            raise ValueError("standardized effect is undefined with zero pooled variance")
        estimate = 0.0
    else:
        correction = 1.0 - 3.0 / (4.0 * df - 1.0)
        estimate = float(correction * difference / np.sqrt(pooled_variance))
    variance = float((nt + nc) / (nt * nc) + estimate**2 / (2.0 * df))
    se = float(np.sqrt(variance))
    critical = float(stats.norm.ppf(1.0 - alpha / 2.0))
    return HedgesGResult(
        estimate=estimate,
        variance=variance,
        standard_error=se,
        ci_low=float(estimate - critical * se),
        ci_high=float(estimate + critical * se),
        n_treatment=nt,
        n_control=nc,
    )


def _validated_effects(
    effects: Sequence[float], variances: Sequence[float], *, min_k: int
) -> tuple[np.ndarray, np.ndarray]:
    y = np.asarray(effects, dtype=float)
    v = np.asarray(variances, dtype=float)
    if y.ndim != 1 or v.ndim != 1 or len(y) != len(v) or len(y) < min_k:
        raise ValueError(f"need at least {min_k} matched effect/variance pairs")
    if np.any(~np.isfinite(y)) or np.any(~np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("effects and strictly positive variances must be finite")
    return y, v


def combine_fixed_effects(
    effects: Sequence[float],
    variances: Sequence[float],
    *,
    alpha: float = 0.05,
) -> FixedEffectResult:
    """Combine exchangeability strata by inverse-variance fixed effects."""

    y, v = _validated_effects(effects, variances, min_k=1)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie between zero and one")
    w = 1.0 / v
    normalized = w / np.sum(w)
    estimate = float(np.sum(normalized * y))
    variance = float(1.0 / np.sum(w))
    se = float(np.sqrt(variance))
    z = float(estimate / se)
    critical = float(stats.norm.ppf(1.0 - alpha / 2.0))
    q = float(np.sum(w * (y - estimate) ** 2))
    q_df = len(y) - 1
    q_p = float(stats.chi2.sf(q, q_df)) if q_df > 0 else np.nan
    return FixedEffectResult(
        estimate=estimate,
        variance=variance,
        standard_error=se,
        ci_low=float(estimate - critical * se),
        ci_high=float(estimate + critical * se),
        z=z,
        p=float(2.0 * stats.norm.sf(abs(z))),
        q=q,
        q_p=q_p,
        weights=normalized,
        k=len(y),
    )


def _reml_objective(tau2: float, y: np.ndarray, v: np.ndarray) -> float:
    w = 1.0 / (v + tau2)
    estimate = float(np.sum(w * y) / np.sum(w))
    return float(
        0.5
        * (
            np.sum(np.log(v + tau2))
            + np.log(np.sum(w))
            + np.sum(w * (y - estimate) ** 2)
        )
    )


def _estimate_tau2_reml(y: np.ndarray, v: np.ndarray) -> float:
    """Profile restricted-likelihood estimate, including the zero boundary."""

    empirical = float(np.var(y, ddof=1)) if len(y) > 1 else 0.0
    upper = max(empirical * 16.0, float(np.max(v)) * 16.0, 1.0)
    baseline = _reml_objective(0.0, y, v)
    # Grow the search interval if its right edge is still descending.
    for _ in range(20):
        midpoint = upper / 2.0
        if _reml_objective(upper, y, v) >= _reml_objective(midpoint, y, v):
            break
        upper *= 4.0
    fit = optimize.minimize_scalar(
        _reml_objective,
        args=(y, v),
        bounds=(0.0, upper),
        method="bounded",
        options={"xatol": 1e-12},
    )
    if not fit.success:
        raise RuntimeError("REML tau-squared optimization failed")
    candidate = float(max(fit.x, 0.0))
    # Bounded minimizers do not evaluate the boundary exactly.
    if baseline <= _reml_objective(candidate, y, v) + 1e-12:
        return 0.0
    return candidate


def random_effects_reml_mkh(
    effects: Sequence[float],
    variances: Sequence[float],
    *,
    alpha: float = 0.05,
) -> RandomEffectsResult:
    """REML random-effects meta-analysis with modified HK uncertainty.

    The Hartung--Knapp scale is floored at one (the modified HK rule), so a
    very homogeneous small set of missions cannot yield a standard error below
    the conventional random-effects standard error.  The prediction interval
    uses ``t_(k-2) * sqrt(tau^2 + SE_mKH^2)`` and is reported only for at least
    three missions.
    """

    y, v = _validated_effects(effects, variances, min_k=2)
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie between zero and one")
    k = len(y)
    df = k - 1
    tau2 = _estimate_tau2_reml(y, v)
    w = 1.0 / (v + tau2)
    normalized = w / np.sum(w)
    estimate = float(np.sum(normalized * y))
    q_random = float(np.sum(w * (y - estimate) ** 2))
    hk_scale = max(1.0, q_random / df)
    se = float(np.sqrt(hk_scale / np.sum(w)))
    t_value = float(estimate / se)
    critical = float(stats.t.ppf(1.0 - alpha / 2.0, df))

    fixed_w = 1.0 / v
    fixed_estimate = float(np.sum(fixed_w * y) / np.sum(fixed_w))
    q = float(np.sum(fixed_w * (y - fixed_estimate) ** 2))
    q_p = float(stats.chi2.sf(q, df))
    i_squared = float(max(0.0, (q - df) / q) * 100.0) if q > EPS else 0.0

    prediction_df: int | None = k - 2 if k >= 3 else None
    if prediction_df is None:
        prediction_low = prediction_high = np.nan
    else:
        pred_critical = float(stats.t.ppf(1.0 - alpha / 2.0, prediction_df))
        pred_se = float(np.sqrt(tau2 + se**2))
        prediction_low = float(estimate - pred_critical * pred_se)
        prediction_high = float(estimate + pred_critical * pred_se)

    return RandomEffectsResult(
        estimate=estimate,
        tau2=tau2,
        standard_error=se,
        t=t_value,
        df=df,
        ci_low=float(estimate - critical * se),
        ci_high=float(estimate + critical * se),
        p=float(2.0 * stats.t.sf(abs(t_value), df)),
        prediction_low=prediction_low,
        prediction_high=prediction_high,
        prediction_df=prediction_df,
        q=q,
        q_p=q_p,
        i_squared=i_squared,
        weights=normalized,
        k=k,
    )


def _batch_reml_score(tau2: np.ndarray, y: np.ndarray, v: np.ndarray) -> np.ndarray:
    w = 1.0 / (v + tau2[:, None])
    sum_w = np.sum(w, axis=1)
    estimate = np.sum(w * y, axis=1) / sum_w
    residual = y - estimate[:, None]
    sum_w2 = np.sum(w**2, axis=1)
    return np.sum((w**2) * (residual**2), axis=1) - (
        sum_w - sum_w2 / sum_w
    )


def _batch_tau2_reml(y: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Vectorized solution of the intercept-only REML score equation."""

    n_rows = y.shape[0]
    zero = np.zeros(n_rows, dtype=float)
    score_zero = _batch_reml_score(zero, y, v)
    active = score_zero > 1e-12
    if not np.any(active):
        return zero

    empirical = np.var(y, axis=1, ddof=1)
    upper = np.maximum.reduce(
        [empirical * 4.0, np.max(v, axis=1) * 4.0, np.full(n_rows, 1e-6)]
    )
    upper[~active] = 0.0
    for _ in range(60):
        score_upper = _batch_reml_score(np.maximum(upper, EPS), y, v)
        unbracketed = active & (score_upper > 0.0)
        if not np.any(unbracketed):
            break
        upper[unbracketed] *= 2.0
    else:  # pragma: no cover - defensive guard for pathological inputs
        raise RuntimeError("could not bracket batched REML tau-squared roots")

    lower = np.zeros(n_rows, dtype=float)
    for _ in range(55):
        midpoint = (lower + upper) / 2.0
        score_mid = _batch_reml_score(np.maximum(midpoint, EPS), y, v)
        move_lower = active & (score_mid > 0.0)
        lower[move_lower] = midpoint[move_lower]
        upper[active & ~move_lower] = midpoint[active & ~move_lower]
    tau2 = (lower + upper) / 2.0
    tau2[~active] = 0.0
    return tau2


def _batch_reml_mkh(
    effects: np.ndarray, variances: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized REML/mKH estimates for rows of small meta-analyses."""

    if effects.shape != variances.shape or effects.ndim != 2 or effects.shape[1] < 2:
        raise ValueError("batch effects/variances must be matched rows with k >= 2")
    if (
        np.any(~np.isfinite(effects))
        or np.any(~np.isfinite(variances))
        or np.any(variances <= 0)
    ):
        raise ValueError("batch effects and positive variances must be finite")
    tau2 = _batch_tau2_reml(effects, variances)
    w = 1.0 / (variances + tau2[:, None])
    sum_w = np.sum(w, axis=1)
    estimate = np.sum(w * effects, axis=1) / sum_w
    q_random = np.sum(w * (effects - estimate[:, None]) ** 2, axis=1)
    scale = np.maximum(1.0, q_random / (effects.shape[1] - 1))
    standard_error = np.sqrt(scale / sum_w)
    t_value = estimate / standard_error
    return estimate, tau2, standard_error, t_value


@dataclass(frozen=True)
class _PermutationBlock:
    mission: str
    stratum: str
    values: np.ndarray
    observed_treatment: np.ndarray
    n_treatment: int
    n_control: int


def _prepare_permutation_blocks(
    scores: pd.DataFrame,
    design: pd.DataFrame,
    *,
    mission_col: str,
    stratum_col: str | None,
    group_col: str,
    treatment_label: str,
    control_label: str | None,
) -> tuple[list[_PermutationBlock], list[str], list[str]]:
    if not isinstance(scores, pd.DataFrame) or scores.empty:
        raise ValueError("scores must be a non-empty samples-by-axes DataFrame")
    if scores.index.has_duplicates or scores.columns.has_duplicates:
        raise ValueError("score sample and axis labels must be unique")
    if design.index.has_duplicates:
        raise ValueError("design sample index contains duplicates")
    if set(scores.index) != set(design.index):
        raise ValueError("scores and design must contain exactly the same samples")
    required = {mission_col, group_col}
    if stratum_col is not None:
        required.add(stratum_col)
    missing_columns = required - set(design.columns)
    if missing_columns:
        raise ValueError("design is missing columns: " + ", ".join(sorted(missing_columns)))
    values = scores.astype(float)
    if np.any(~np.isfinite(values.to_numpy())):
        raise ValueError("all axis scores must be finite before permutation")
    aligned = design.loc[values.index].copy()
    if stratum_col is None:
        aligned["__single_stratum__"] = "all"
        stratum_col = "__single_stratum__"

    groups = set(aligned[group_col].astype(str))
    if treatment_label not in groups:
        raise ValueError(f"treatment label {treatment_label!r} is absent")
    if control_label is None:
        alternatives = groups - {treatment_label}
        if len(alternatives) != 1:
            raise ValueError("control_label is required when the design has more than two groups")
        control_label = next(iter(alternatives))
    if groups != {treatment_label, control_label}:
        raise ValueError("design must contain exactly the declared treatment and control labels")

    missions = [str(value) for value in pd.unique(aligned[mission_col])]
    axes = [str(value) for value in values.columns]
    blocks: list[_PermutationBlock] = []
    for mission in missions:
        mission_mask = aligned[mission_col].astype(str).eq(mission)
        strata = pd.unique(aligned.loc[mission_mask, stratum_col])
        for stratum_value in strata:
            mask = mission_mask & aligned[stratum_col].eq(stratum_value)
            idx = aligned.index[mask]
            labels = aligned.loc[idx, group_col].astype(str).to_numpy()
            treatment = labels == treatment_label
            nt = int(np.sum(treatment))
            nc = int(len(treatment) - nt)
            if nt < 2 or nc < 2:
                raise ValueError(
                    f"block {mission}/{stratum_value} needs at least two samples per group"
                )
            blocks.append(
                _PermutationBlock(
                    mission=mission,
                    stratum=str(stratum_value),
                    values=values.loc[idx].to_numpy(dtype=float),
                    observed_treatment=treatment,
                    n_treatment=nt,
                    n_control=nc,
                )
            )
    return blocks, missions, axes


def _block_hedges_g_batch(
    values: np.ndarray, treatment_masks: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Hedges g for a batch of treatment allocations in one block."""

    mask = treatment_masks.astype(float)
    nt = treatment_masks.sum(axis=1).astype(float)
    nc = values.shape[0] - nt
    treatment_sum = mask @ values
    treatment_ss = mask @ (values**2)
    total_sum = np.sum(values, axis=0)[None, :]
    total_ss = np.sum(values**2, axis=0)[None, :]
    control_sum = total_sum - treatment_sum
    control_ss = total_ss - treatment_ss
    treatment_mean = treatment_sum / nt[:, None]
    control_mean = control_sum / nc[:, None]
    treatment_var = (treatment_ss - treatment_sum**2 / nt[:, None]) / (nt[:, None] - 1.0)
    control_var = (control_ss - control_sum**2 / nc[:, None]) / (nc[:, None] - 1.0)
    df = nt + nc - 2.0
    pooled = (
        (nt[:, None] - 1.0) * treatment_var
        + (nc[:, None] - 1.0) * control_var
    ) / df[:, None]
    if np.any(pooled <= EPS):
        raise ValueError("a label allocation produced zero within-group pooled variance")
    correction = 1.0 - 3.0 / (4.0 * df - 1.0)
    g = correction[:, None] * (treatment_mean - control_mean) / np.sqrt(pooled)
    variance = (nt + nc)[:, None] / (nt * nc)[:, None] + g**2 / (2.0 * df[:, None])
    return g, variance


def _mission_effect_batch(
    blocks: Sequence[_PermutationBlock],
    missions: Sequence[str],
    treatment_masks: Sequence[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    block_results = [
        _block_hedges_g_batch(block.values, mask)
        for block, mask in zip(blocks, treatment_masks, strict=True)
    ]
    batch_size, n_axes = block_results[0][0].shape
    effects = np.empty((batch_size, n_axes, len(missions)), dtype=float)
    variances = np.empty_like(effects)
    for mission_index, mission in enumerate(missions):
        indices = [index for index, block in enumerate(blocks) if block.mission == mission]
        stratum_effects = np.stack([block_results[index][0] for index in indices], axis=2)
        stratum_variances = np.stack([block_results[index][1] for index in indices], axis=2)
        weights = 1.0 / stratum_variances
        effects[:, :, mission_index] = np.sum(weights * stratum_effects, axis=2) / np.sum(
            weights, axis=2
        )
        variances[:, :, mission_index] = 1.0 / np.sum(weights, axis=2)
    return effects, variances


def _observed_mission_effect_table(
    effects: np.ndarray,
    variances: np.ndarray,
    axes: Sequence[str],
    missions: Sequence[str],
    blocks: Sequence[_PermutationBlock],
) -> pd.DataFrame:
    rows = []
    for axis_index, axis in enumerate(axes):
        for mission_index, mission in enumerate(missions):
            rows.append(
                {
                    "axis": axis,
                    "mission": mission,
                    "estimate": float(effects[0, axis_index, mission_index]),
                    "variance": float(variances[0, axis_index, mission_index]),
                    "standard_error": float(
                        np.sqrt(variances[0, axis_index, mission_index])
                    ),
                    "n_strata": sum(block.mission == mission for block in blocks),
                }
            )
    return pd.DataFrame(rows)


def blocked_meta_permutation(
    scores: pd.DataFrame,
    design: pd.DataFrame,
    *,
    mission_col: str = "mission",
    stratum_col: str | None = "stratum",
    group_col: str = "group",
    treatment_label: str = "flight",
    control_label: str | None = "control",
    n_permutations: int = 10_000,
    seed: int = 0,
    chunk_size: int = 1_024,
) -> MetaPermutationResult:
    """Run a whole-pipeline, exchangeability-blocked max-|T| permutation.

    Treatment labels are shuffled independently inside each mission/stratum
    block while retaining its observed treatment count.  Every allocation is
    transformed into stratum-level Hedges g, fixed-effect mission estimates,
    and a cross-mission REML/mKH t statistic.  The maximum absolute t statistic
    across the supplied axes controls the single frozen axis family.

    The plus-one Monte Carlo correction is used for both unadjusted and max-T
    p-values.  Gene standardization and scoring must be completed before this
    function; because they are label-blind, they need not be recalculated under
    each permutation.
    """

    if n_permutations < 1:
        raise ValueError("n_permutations must be positive")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    blocks, missions, axes = _prepare_permutation_blocks(
        scores,
        design,
        mission_col=mission_col,
        stratum_col=stratum_col,
        group_col=group_col,
        treatment_label=treatment_label,
        control_label=control_label,
    )
    if len(missions) < 2:
        raise ValueError("cross-mission permutation inference requires at least two missions")

    observed_masks = [block.observed_treatment[None, :] for block in blocks]
    observed_effects, observed_variances = _mission_effect_batch(
        blocks, missions, observed_masks
    )
    mission_table = _observed_mission_effect_table(
        observed_effects, observed_variances, axes, missions, blocks
    )

    observed_rows = []
    observed_t = np.empty(len(axes), dtype=float)
    for axis_index, axis in enumerate(axes):
        result = random_effects_reml_mkh(
            observed_effects[0, axis_index], observed_variances[0, axis_index]
        )
        observed_t[axis_index] = result.t
        observed_rows.append(
            {
                "axis": axis,
                "estimate": result.estimate,
                "tau2": result.tau2,
                "standard_error_mkh": result.standard_error,
                "t_mkh": result.t,
                "df": result.df,
                "ci_low_mkh": result.ci_low,
                "ci_high_mkh": result.ci_high,
                "p_mkh": result.p,
                "prediction_low": result.prediction_low,
                "prediction_high": result.prediction_high,
                "prediction_df": result.prediction_df,
                "q": result.q,
                "q_p": result.q_p,
                "i_squared": result.i_squared,
                "k_missions": result.k,
            }
        )

    null_t = np.empty((n_permutations, len(axes)), dtype=float)
    seed_sequence = np.random.SeedSequence(seed)
    block_rngs = [
        np.random.default_rng(child)
        for child in seed_sequence.spawn(len(blocks))
    ]
    for start in range(0, n_permutations, chunk_size):
        stop = min(start + chunk_size, n_permutations)
        batch = stop - start
        masks: list[np.ndarray] = []
        for block, rng in zip(blocks, block_rngs, strict=True):
            random_keys = rng.random((batch, block.values.shape[0]))
            chosen = np.argpartition(
                random_keys, block.n_treatment - 1, axis=1
            )[:, : block.n_treatment]
            treatment_mask = np.zeros_like(random_keys, dtype=bool)
            treatment_mask[np.arange(batch)[:, None], chosen] = True
            masks.append(treatment_mask)
        perm_effects, perm_variances = _mission_effect_batch(blocks, missions, masks)
        flattened_effects = perm_effects.reshape(batch * len(axes), len(missions))
        flattened_variances = perm_variances.reshape(batch * len(axes), len(missions))
        _, _, _, t_values = _batch_reml_mkh(flattened_effects, flattened_variances)
        null_t[start:stop] = t_values.reshape(batch, len(axes))

    null_max = np.max(np.abs(null_t), axis=1)
    empirical = np.array(
        [
            (1.0 + np.count_nonzero(np.abs(null_t[:, index]) >= abs(observed_t[index])))
            / (n_permutations + 1.0)
            for index in range(len(axes))
        ]
    )
    fwer = np.array(
        [
            (1.0 + np.count_nonzero(null_max >= abs(observed_t[index])))
            / (n_permutations + 1.0)
            for index in range(len(axes))
        ]
    )
    observed_meta = pd.DataFrame(observed_rows).set_index("axis")
    observed_meta["empirical_p_two_sided"] = empirical
    observed_meta["max_t_fwer"] = fwer
    return MetaPermutationResult(
        observed_meta=observed_meta,
        mission_effects=mission_table,
        null_t=pd.DataFrame(null_t, columns=axes),
        null_max_abs_t=null_max,
        empirical_p_two_sided=pd.Series(empirical, index=axes, name="empirical_p_two_sided"),
        max_t_fwer=pd.Series(fwer, index=axes, name="max_t_fwer"),
        n_permutations=n_permutations,
        seed=seed,
    )


def leave_one_mission_out(
    effects: Sequence[float],
    variances: Sequence[float],
    mission_labels: Sequence[str],
    *,
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Repeat REML/mKH synthesis after omitting each mission in turn."""

    y, v = _validated_effects(effects, variances, min_k=3)
    labels = list(mission_labels)
    if len(labels) != len(y) or len(set(labels)) != len(labels):
        raise ValueError("mission_labels must be unique and match the effects")
    rows = []
    for index, label in enumerate(labels):
        keep = np.arange(len(y)) != index
        result = random_effects_reml_mkh(y[keep], v[keep], alpha=alpha)
        rows.append(
            {
                "omitted_mission": label,
                "estimate": result.estimate,
                "tau2": result.tau2,
                "standard_error_mkh": result.standard_error,
                "ci_low_mkh": result.ci_low,
                "ci_high_mkh": result.ci_high,
                "p_mkh": result.p,
                "i_squared": result.i_squared,
                "k_missions": result.k,
            }
        )
    return pd.DataFrame(rows)


def leave_one_gene_out_scores(
    expression: pd.DataFrame,
    directions: Mapping[str, float] | pd.Series,
    *,
    method: str = "mean",
    subdomains: Mapping[str, Sequence[str]] | None = None,
    min_genes_per_subdomain: int = 1,
    require_all_genes: bool = False,
) -> pd.DataFrame:
    """Return one signed axis-score vector for every observable gene omission."""

    signs = _validated_signs(directions)
    observable = [gene for gene in signs.index if gene in expression.index]
    if len(observable) < 2:
        raise ValueError("leave-one-gene-out scoring requires at least two observable genes")
    columns: dict[str, pd.Series] = {}
    for omitted in observable:
        reduced_signs = signs.drop(index=omitted)
        reduced_domains = None
        if subdomains is not None:
            reduced_domains = {
                name: [gene for gene in genes if gene != omitted]
                for name, genes in subdomains.items()
            }
            if any(len(genes) == 0 for genes in reduced_domains.values()):
                raise ValueError(
                    f"omitting {omitted!r} empties an equal-weight subdomain"
                )
        columns[str(omitted)] = score_signed_axis(
            expression,
            reduced_signs,
            method=method,
            subdomains=reduced_domains,
            min_genes_per_subdomain=min_genes_per_subdomain,
            require_all_genes=require_all_genes,
        ).scores
    result = pd.DataFrame(columns, index=expression.columns)
    result.columns.name = "omitted_gene"
    return result
