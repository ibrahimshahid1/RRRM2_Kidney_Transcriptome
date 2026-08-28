"""Core statistics for the frozen Grey60 adversarial reanalysis.

The functions here contain no repository path logic.  The runner is
responsible for loading the frozen inputs and writing manifests.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import re
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import optimize, stats


EPS = 1e-12


@dataclass(frozen=True)
class HedgesG:
    estimate: float
    variance: float
    standard_error: float
    ci_low: float
    ci_high: float
    n_flight: int
    n_control: int


@dataclass(frozen=True)
class RandomEffects:
    estimate: float
    tau2: float
    standard_error_hk: float
    ci_low_hk: float
    ci_high_hk: float
    p_hk: float
    i_squared: float
    q: float
    weights: np.ndarray


def osd462_animal_id(sample: str) -> str:
    """Remove the RNA-preparation and optional technical-replicate suffix."""
    return re.sub(r"_(UPX|mRNA|totRNA)(?:_techrep\d+)?$", "", sample)


def zscore_rows(values: np.ndarray) -> np.ndarray:
    """Z-score genes (rows) across samples, with constant rows set to zero."""
    x = np.asarray(values, dtype=float)
    mu = np.nanmean(x, axis=1, keepdims=True)
    sd = np.nanstd(x, axis=1, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd < EPS)] = 1.0
    z = (x - mu) / sd
    return np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)


def mean_z_score(expression: pd.DataFrame, genes: Sequence[str]) -> pd.Series:
    """Return an unweighted mean gene-z score for the requested genes."""
    present = [g for g in genes if g in expression.index]
    if not present:
        raise ValueError("No requested genes were present in the expression matrix")
    z = zscore_rows(expression.loc[present].to_numpy(dtype=float))
    return pd.Series(z.mean(axis=0), index=expression.columns, name="mean_z")


def median_z_score(expression: pd.DataFrame, genes: Sequence[str]) -> pd.Series:
    """Return a median gene-z score for the requested genes."""
    present = [g for g in genes if g in expression.index]
    if not present:
        raise ValueError("No requested genes were present in the expression matrix")
    z = zscore_rows(expression.loc[present].to_numpy(dtype=float))
    return pd.Series(np.median(z, axis=0), index=expression.columns, name="median_z")


def weighted_mean_z_score(
    expression: pd.DataFrame,
    genes: Sequence[str],
    weights: pd.Series,
) -> pd.Series:
    """Return a fixed-weight mean gene-z score.

    Weights retain their frozen sign and are normalized by the sum of their
    absolute values, so an arbitrary scale change cannot alter the score.
    """
    present = [g for g in genes if g in expression.index and g in weights.index]
    if not present:
        raise ValueError("No weighted genes were present")
    w = weights.loc[present].to_numpy(dtype=float)
    denom = float(np.abs(w).sum())
    if denom < EPS:
        raise ValueError("Frozen weights sum to zero")
    z = zscore_rows(expression.loc[present].to_numpy(dtype=float))
    return pd.Series((w[:, None] * z).sum(axis=0) / denom, index=expression.columns)


def pooled_iss_effect(
    scores: pd.Series,
    metadata: pd.DataFrame,
    *,
    arm: str = "ISS-T",
    flight_label: str = "FLT",
    control_label: str = "GC",
) -> float:
    """Average the age-specific FLT-minus-GC mean differences."""
    diffs: list[float] = []
    for age in ("YNG", "OLD"):
        idx_f = metadata.index[
            (metadata["Arm"] == arm)
            & (metadata["Age"] == age)
            & (metadata["EnvGroup"] == flight_label)
        ]
        idx_c = metadata.index[
            (metadata["Arm"] == arm)
            & (metadata["Age"] == age)
            & (metadata["EnvGroup"] == control_label)
        ]
        if len(idx_f) == 0 or len(idx_c) == 0:
            raise ValueError(f"Missing cell for {arm}/{age}")
        diffs.append(float(scores.loc[idx_f].mean() - scores.loc[idx_c].mean()))
    return float(np.mean(diffs))


def age_specific_effects(
    scores: pd.Series,
    metadata: pd.DataFrame,
    *,
    arm: str = "ISS-T",
) -> dict[str, float]:
    return {
        age: float(
            scores.loc[
                metadata.index[
                    (metadata["Arm"] == arm)
                    & (metadata["Age"] == age)
                    & (metadata["EnvGroup"] == "FLT")
                ]
            ].mean()
            - scores.loc[
                metadata.index[
                    (metadata["Arm"] == arm)
                    & (metadata["Age"] == age)
                    & (metadata["EnvGroup"] == "GC")
                ]
            ].mean()
        )
        for age in ("YNG", "OLD")
    }


def stratified_bootstrap_iss_effect(
    scores: pd.Series,
    metadata: pd.DataFrame,
    *,
    n_boot: int,
    seed: int,
    arm: str = "ISS-T",
) -> np.ndarray:
    """Bootstrap animals inside Age x condition cells and pool over age."""
    rng = np.random.default_rng(seed)
    cell_values: dict[tuple[str, str], np.ndarray] = {}
    for age in ("YNG", "OLD"):
        for env in ("FLT", "GC"):
            idx = metadata.index[
                (metadata["Arm"] == arm)
                & (metadata["Age"] == age)
                & (metadata["EnvGroup"] == env)
            ]
            vals = scores.loc[idx].to_numpy(dtype=float)
            if len(vals) < 2:
                raise ValueError(f"Too few animals in {arm}/{age}/{env}")
            cell_values[(age, env)] = vals

    out = np.empty(n_boot, dtype=float)
    for start in range(0, n_boot, 2048):
        stop = min(start + 2048, n_boot)
        b = stop - start
        age_diffs = []
        for age in ("YNG", "OLD"):
            vf = cell_values[(age, "FLT")]
            vc = cell_values[(age, "GC")]
            sf = vf[rng.integers(0, len(vf), size=(b, len(vf)))].mean(axis=1)
            sc = vc[rng.integers(0, len(vc), size=(b, len(vc)))].mean(axis=1)
            age_diffs.append(sf - sc)
        out[start:stop] = (age_diffs[0] + age_diffs[1]) / 2.0
    return out


def _contrast_matrix() -> tuple[list[str], np.ndarray]:
    """Contrasts over eight cell means.

    Cell order:
    ISS-Y GC/F, ISS-O GC/F, LAR-Y GC/F, LAR-O GC/F.
    """
    names = [
        "Flight",
        "Age",
        "Arm",
        "Flight:Age",
        "Flight:Arm",
        "Age:Arm",
        "Flight:Age:Arm",
        "ISS-T_YNG",
        "ISS-T_OLD",
        "LAR_YNG",
        "LAR_OLD",
    ]
    c = np.array(
        [
            [-1, 1, 0, 0, 0, 0, 0, 0],
            [-1, 0, 1, 0, 0, 0, 0, 0],
            [-1, 0, 0, 0, 1, 0, 0, 0],
            [1, -1, -1, 1, 0, 0, 0, 0],
            [1, -1, 0, 0, -1, 1, 0, 0],
            [1, 0, -1, 0, -1, 0, 1, 0],
            [-1, 1, 1, -1, 1, -1, -1, 1],
            [-1, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, -1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, -1, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, -1, 1],
        ],
        dtype=float,
    )
    return names, c


def _stratum_indices(metadata: pd.DataFrame) -> list[np.ndarray]:
    strata = []
    for arm, age in (("ISS-T", "YNG"), ("ISS-T", "OLD"), ("LAR", "YNG"), ("LAR", "OLD")):
        idx = np.flatnonzero(
            (metadata["Arm"].to_numpy() == arm)
            & (metadata["Age"].to_numpy() == age)
        )
        if len(idx) != 10:
            raise ValueError(f"Expected 10 samples in {arm}/{age}; found {len(idx)}")
        strata.append(idx)
    return strata


def observed_saturated_t(
    responses: pd.DataFrame,
    metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Compute the 11 saturated-model t statistics for each response."""
    y = responses.loc[metadata.index].to_numpy(dtype=float)
    strata = _stratum_indices(metadata)
    cell_means = np.empty((8, y.shape[1]), dtype=float)
    sse = np.zeros(y.shape[1], dtype=float)
    for j, idx in enumerate(strata):
        labels = (metadata.iloc[idx]["EnvGroup"].to_numpy() == "FLT")
        if labels.sum() != 5:
            raise ValueError("Each stratum must contain 5 FLT and 5 GC animals")
        yf = y[idx][labels]
        yc = y[idx][~labels]
        cell_means[2 * j] = yc.mean(axis=0)
        cell_means[2 * j + 1] = yf.mean(axis=0)
        sse += ((yf - yf.mean(axis=0)) ** 2).sum(axis=0)
        sse += ((yc - yc.mean(axis=0)) ** 2).sum(axis=0)
    mse = sse / 32.0
    names, contrasts = _contrast_matrix()
    estimates = contrasts @ cell_means
    norm = (contrasts**2).sum(axis=1) / 5.0
    se = np.sqrt(norm[:, None] * mse[None, :])
    t = np.divide(estimates, se, out=np.zeros_like(estimates), where=se > EPS)
    rows = []
    for i, term in enumerate(names):
        for j, response in enumerate(responses.columns):
            rows.append(
                {
                    "response": response,
                    "term": term,
                    "estimate": estimates[i, j],
                    "standard_error": se[i, j],
                    "t": t[i, j],
                    "p_t": 2 * stats.t.sf(abs(t[i, j]), 32),
                }
            )
    return pd.DataFrame(rows)


def max_t_permutation(
    module_responses: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    n_perm: int,
    seed: int,
    chunk_size: int = 2048,
) -> np.ndarray:
    """Generate a blocked max-|t| null across 11 contrasts x modules.

    Labels are permuted independently inside the four 5/5 Age x Arm strata.
    The implementation uses saturated cell means and a pooled 32-df residual
    variance, exactly matching the historical factorial OLS model.
    """
    if module_responses.shape[1] != 20:
        raise ValueError(
            f"Selection family requires 20 non-grey modules; got {module_responses.shape[1]}"
        )
    y = module_responses.loc[metadata.index].to_numpy(dtype=float)
    strata = _stratum_indices(metadata)
    combos = np.array(list(combinations(range(10), 5)), dtype=int)
    rng = np.random.default_rng(seed)
    _, contrasts = _contrast_matrix()
    contrast_norm = (contrasts**2).sum(axis=1) / 5.0

    sst = np.zeros(y.shape[1], dtype=float)
    for idx in strata:
        ys = y[idx]
        sst += ((ys - ys.mean(axis=0)) ** 2).sum(axis=0)

    out = np.empty(n_perm, dtype=float)
    for start in range(0, n_perm, chunk_size):
        stop = min(start + chunk_size, n_perm)
        b = stop - start
        cell_means = np.empty((b, 8, y.shape[1]), dtype=float)
        sse = np.broadcast_to(sst, (b, y.shape[1])).copy()
        for j, idx in enumerate(strata):
            ys = y[idx]
            choices = combos[rng.integers(0, len(combos), size=b)]
            labels = np.zeros((b, 10), dtype=float)
            labels[np.arange(b)[:, None], choices] = 1.0
            fmean = labels @ ys / 5.0
            cmean = (1.0 - labels) @ ys / 5.0
            diff = fmean - cmean
            sse -= 2.5 * diff**2
            cell_means[:, 2 * j, :] = cmean
            cell_means[:, 2 * j + 1, :] = fmean
        mse = np.maximum(sse / 32.0, EPS)
        est = np.einsum("ce,bem->bcm", contrasts, cell_means)
        se = np.sqrt(mse[:, None, :] * contrast_norm[None, :, None])
        t = est / se
        out[start:stop] = np.max(np.abs(t), axis=(1, 2))
    return out


def familywise_p(null_max_abs_t: np.ndarray, observed_t: float) -> float:
    return float(
        (1 + np.count_nonzero(null_max_abs_t >= abs(observed_t)))
        / (len(null_max_abs_t) + 1)
    )


def hedges_g(flight: Iterable[float], control: Iterable[float], alpha: float = 0.05) -> HedgesG:
    """Bias-corrected standardized mean difference and conventional variance."""
    xf = np.asarray(list(flight), dtype=float)
    xc = np.asarray(list(control), dtype=float)
    nf, nc = len(xf), len(xc)
    if nf < 2 or nc < 2:
        raise ValueError("Hedges g requires at least two observations per group")
    df = nf + nc - 2
    pooled_var = ((nf - 1) * xf.var(ddof=1) + (nc - 1) * xc.var(ddof=1)) / df
    if pooled_var < EPS:
        estimate = 0.0
    else:
        d = (xf.mean() - xc.mean()) / np.sqrt(pooled_var)
        j = 1.0 - 3.0 / (4.0 * df - 1.0)
        estimate = j * d
    variance = (nf + nc) / (nf * nc) + estimate**2 / (2.0 * df)
    se = float(np.sqrt(variance))
    z = stats.norm.ppf(1 - alpha / 2)
    return HedgesG(
        estimate=float(estimate),
        variance=float(variance),
        standard_error=se,
        ci_low=float(estimate - z * se),
        ci_high=float(estimate + z * se),
        n_flight=nf,
        n_control=nc,
    )


def random_effects_reml_hk(
    effects: Sequence[float],
    variances: Sequence[float],
    alpha: float = 0.05,
    *,
    modified: bool = False,
) -> RandomEffects:
    """Random-effects REML synthesis with Hartung-Knapp uncertainty.

    When ``modified`` is true, the Hartung-Knapp scale factor is floored at
    one. This prevents unusually homogeneous small meta-analyses from
    producing intervals narrower than the conventional random-effects
    interval.
    """
    y = np.asarray(effects, dtype=float)
    v = np.asarray(variances, dtype=float)
    if len(y) < 2 or len(y) != len(v):
        raise ValueError("Need at least two effect/variance pairs")
    if np.any(~np.isfinite(y)) or np.any(~np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("Effects and positive variances must be finite")

    def reml_nll(tau2: float) -> float:
        w = 1.0 / (v + tau2)
        mu = float(np.sum(w * y) / np.sum(w))
        return 0.5 * (
            np.sum(np.log(v + tau2))
            + np.log(np.sum(w))
            + np.sum(w * (y - mu) ** 2)
        )

    upper = max(float(np.var(y, ddof=1) * 10.0), 1.0)
    fit = optimize.minimize_scalar(reml_nll, bounds=(0.0, upper), method="bounded")
    tau2 = float(max(fit.x, 0.0))
    w = 1.0 / (v + tau2)
    mu = float(np.sum(w * y) / np.sum(w))
    q_re = float(np.sum(w * (y - mu) ** 2))
    df = len(y) - 1
    hk_scale = q_re / df
    if modified:
        hk_scale = max(1.0, hk_scale)
    se_hk = float(np.sqrt(max(hk_scale, 0.0) / np.sum(w)))
    if se_hk < EPS:
        ci_low = ci_high = mu
        p = 0.0 if abs(mu) > EPS else 1.0
    else:
        crit = float(stats.t.ppf(1 - alpha / 2, df))
        ci_low, ci_high = mu - crit * se_hk, mu + crit * se_hk
        p = float(2 * stats.t.sf(abs(mu / se_hk), df))

    wf = 1.0 / v
    muf = float(np.sum(wf * y) / np.sum(wf))
    q = float(np.sum(wf * (y - muf) ** 2))
    i2 = float(max(0.0, (q - df) / q) * 100.0) if q > EPS else 0.0
    return RandomEffects(
        estimate=mu,
        tau2=tau2,
        standard_error_hk=se_hk,
        ci_low_hk=float(ci_low),
        ci_high_hk=float(ci_high),
        p_hk=p,
        i_squared=i2,
        q=q,
        weights=w / w.sum(),
    )


def label_permutation_p(
    flight: Sequence[float],
    control: Sequence[float],
    *,
    n_perm: int,
    seed: int,
) -> float:
    """Two-sided Monte Carlo label-permutation p value for a mean difference."""
    xf = np.asarray(flight, dtype=float)
    xc = np.asarray(control, dtype=float)
    values = np.concatenate([xf, xc])
    nf = len(xf)
    observed = float(xf.mean() - xc.mean())
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_perm):
        order = rng.permutation(len(values))
        diff = values[order[:nf]].mean() - values[order[nf:]].mean()
        exceed += int(abs(diff) >= abs(observed))
    return float((exceed + 1) / (n_perm + 1))
