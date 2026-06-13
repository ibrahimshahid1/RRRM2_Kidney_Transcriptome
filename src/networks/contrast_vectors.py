"""Contrast-vector decomposition core (Cross-OSDR Network-Contrast Framework)."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

EPS = 1e-12



@dataclass
class Decomposition:
    """One decomposition of A_FLT onto A_GC at a single resolution."""

    beta: float
    cos: float
    rho: float
    a_flt_norm: float
    a_gc_norm: float
    residual_norm: float
    weighted: bool

    def as_dict(self) -> dict:
        return asdict(self)


def _dot(u: np.ndarray, v: np.ndarray, w: np.ndarray | None = None) -> float:
    if w is None:
        return float(np.dot(u, v))
    return float(np.dot(u * w, v))


def _norm(u: np.ndarray, w: np.ndarray | None = None) -> float:
    if w is None:
        return float(np.linalg.norm(u))
    return float(np.sqrt(max(_dot(u, u, w), 0.0)))


def compute_decomposition(
    a_flt: np.ndarray,
    a_gc: np.ndarray,
    weights: np.ndarray | None = None,
) -> Decomposition:
    """Compute beta, cos, rho, ||R|| for one (A_FLT, A_GC) pair.

    Parameters
    ----------
    a_flt, a_gc : np.ndarray
        Flight and ground-control aging vectors. Must share length.
    weights : np.ndarray or None
        Optional non-negative per-feature precision weights (Guardrail C).
    """
    a_flt = np.asarray(a_flt, dtype=float).ravel()
    a_gc = np.asarray(a_gc, dtype=float).ravel()
    if a_flt.shape != a_gc.shape:
        raise ValueError(
            f"A_FLT and A_GC must share length; got {a_flt.shape} vs {a_gc.shape}"
        )
    if weights is not None:
        weights = np.asarray(weights, dtype=float).ravel()
        if weights.shape != a_flt.shape:
            raise ValueError("weights must share length with A_FLT/A_GC")
        if np.any(weights < 0):
            raise ValueError("weights must be non-negative")

    a_gc_norm_sq = _dot(a_gc, a_gc, weights)
    a_gc_norm = float(np.sqrt(max(a_gc_norm_sq, 0.0)))
    a_flt_norm = _norm(a_flt, weights)
    dot_fg = _dot(a_flt, a_gc, weights)

    if a_gc_norm_sq <= EPS:
        beta = float("nan")
        cos = float("nan")
        residual_vec = a_flt.copy()
    else:
        beta = dot_fg / a_gc_norm_sq
        denom = a_flt_norm * a_gc_norm
        cos = (dot_fg / denom) if denom > EPS else float("nan")
        residual_vec = a_flt - beta * a_gc

    residual_norm = _norm(residual_vec, weights)
    rho = (residual_norm / a_flt_norm) if a_flt_norm > EPS else float("nan")

    return Decomposition(
        beta=beta,
        cos=cos,
        rho=rho,
        a_flt_norm=a_flt_norm,
        a_gc_norm=a_gc_norm,
        residual_norm=residual_norm,
        weighted=weights is not None,
    )


def redirected_component(
    a_flt: np.ndarray,
    a_gc: np.ndarray,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    """Return the residual vector R = A_FLT - beta * A_GC.

    Note: beta is computed in the (possibly weighted) inner product space, but
    the residual is returned in the original coordinate basis so callers can
    interpret per-feature R values directly.
    """
    dec = compute_decomposition(a_flt, a_gc, weights=weights)
    if not np.isfinite(dec.beta):
        return np.asarray(a_flt, dtype=float).ravel().copy()
    return np.asarray(a_flt, dtype=float).ravel() - dec.beta * np.asarray(a_gc, dtype=float).ravel()


def cosine(u: np.ndarray, v: np.ndarray, weights: np.ndarray | None = None) -> float:
    """Weighted or unweighted cosine of two vectors."""
    u = np.asarray(u, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    if u.shape != v.shape:
        raise ValueError(f"shape mismatch: {u.shape} vs {v.shape}")
    nu = _norm(u, weights)
    nv = _norm(v, weights)
    if nu <= EPS or nv <= EPS:
        return float("nan")
    return _dot(u, v, weights) / (nu * nv)



def classify_beta(beta: float, rho: float | None = None,
                  amplify_min: float = 1.0,
                  redirect_abs_max: float = 0.10,
                  redirect_rho_min: float = 0.50) -> str:
    """Map a single (beta, rho) point to an interpretation category."""
    if not np.isfinite(beta):
        return "undefined"
    if abs(beta) < redirect_abs_max and rho is not None and np.isfinite(rho) and rho >= redirect_rho_min:
        return "redirect"
    if beta < 0:
        return "reverse"
    if beta > amplify_min:
        return "amplify"
    if 0.0 < beta < amplify_min:
        return "dampen"
    return "boundary"


def categorize_bootstrap_distribution(
    betas: np.ndarray,
    rhos: np.ndarray | None,
    headline_fraction: float = 0.70,
    amplify_min: float = 1.0,
    redirect_abs_max: float = 0.10,
    redirect_rho_min: float = 0.50,
) -> dict:
    """Return per-category fractions and the headline assignment if any.

    The headline category is assigned only when its mass passes
    ``headline_fraction`` of the bootstrap distribution (§4.4).
    """
    betas = np.asarray(betas, dtype=float).ravel()
    rhos = np.asarray(rhos, dtype=float).ravel() if rhos is not None else None
    if rhos is not None and rhos.shape != betas.shape:
        raise ValueError("betas and rhos must share shape when both provided")

    cats = {"amplify": 0, "dampen": 0, "reverse": 0, "redirect": 0,
            "boundary": 0, "undefined": 0}
    for i, b in enumerate(betas):
        r = float(rhos[i]) if rhos is not None else None
        cats[classify_beta(b, r, amplify_min=amplify_min,
                           redirect_abs_max=redirect_abs_max,
                           redirect_rho_min=redirect_rho_min)] += 1
    total = max(len(betas), 1)
    fractions = {k: v / total for k, v in cats.items()}

    headline = None
    for k, frac in fractions.items():
        if k in ("boundary", "undefined"):
            continue
        if frac >= headline_fraction:
            headline = k
            break

    return {
        "fractions": fractions,
        "headline": headline,
        "headline_fraction_required": headline_fraction,
        "n_bootstrap": int(total),
    }


def interpretation_label(summary: Mapping[str, object]) -> str:
    """Return the manuscript-safe bootstrap category label (§4.4).

    A single headline category is used only when the configured fraction
    threshold is met. Otherwise the label explicitly reports the category mix.
    """
    headline = summary.get("headline")
    fractions = summary.get("fractions", {})
    if headline:
        return str(headline)
    parts = []
    for key in ("amplify", "dampen", "reverse", "redirect"):
        frac = float(fractions.get(key, 0.0)) if isinstance(fractions, Mapping) else 0.0
        parts.append(f"{key}={100.0 * frac:.1f}%")
    return "mixed (" + ", ".join(parts) + ")"



def stratified_bootstrap_indices(
    strata: Sequence,
    rng: np.random.Generator,
) -> np.ndarray:
    """Resample with replacement within each unique stratum value.

    Returns an array of indices into the original sample table with the same
    length as ``strata``.
    """
    strata = pd.Series(strata).reset_index(drop=True)
    indices = np.empty(len(strata), dtype=np.int64)
    for value in pd.unique(strata):
        cell_idx = np.where(strata.values == value)[0]
        if cell_idx.size == 0:
            continue
        draws = rng.integers(low=0, high=cell_idx.size, size=cell_idx.size)
        indices[cell_idx] = cell_idx[draws]
    return indices


def stratified_permutation_indices(
    labels: Sequence,
    strata: Sequence,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute labels within each stratum and return the permuted label array."""
    labels = np.asarray(labels)
    strata = pd.Series(strata).reset_index(drop=True)
    out = labels.copy()
    for value in pd.unique(strata):
        idx = np.where(strata.values == value)[0]
        if idx.size <= 1:
            continue
        perm = rng.permutation(idx.size)
        out[idx] = labels[idx[perm]]
    return out



def build_aging_vectors(
    feature_matrix: np.ndarray,
    age: Sequence,
    env: Sequence,
    arm_mask: Sequence[bool] | None = None,
    old_label: str = "OLD",
    young_label: str = "YNG",
    gc_label: str = "GC",
    flt_label: str = "FLT",
) -> tuple[np.ndarray, np.ndarray]:
    """Return (A_GC, A_FLT) = (Old-Young) in GC and FLT subsets.

    ``feature_matrix`` is samples x features. Aggregation is the within-cell
    mean. The optional ``arm_mask`` restricts to a single arm before computing
    the contrast.
    """
    feature_matrix = np.asarray(feature_matrix, dtype=float)
    age = np.asarray(age)
    env = np.asarray(env)
    n = feature_matrix.shape[0]
    if age.shape[0] != n or env.shape[0] != n:
        raise ValueError("age/env length must equal samples in feature_matrix")
    keep = np.ones(n, dtype=bool) if arm_mask is None else np.asarray(arm_mask, dtype=bool)
    if keep.shape[0] != n:
        raise ValueError("arm_mask length must equal samples in feature_matrix")

    def cell_mean(env_value: str, age_value: str) -> np.ndarray:
        mask = keep & (env == env_value) & (age == age_value)
        if not mask.any():
            raise ValueError(
                f"No samples for env={env_value}, age={age_value} after arm mask"
            )
        return feature_matrix[mask].mean(axis=0)

    a_gc = cell_mean(gc_label, old_label) - cell_mean(gc_label, young_label)
    a_flt = cell_mean(flt_label, old_label) - cell_mean(flt_label, young_label)
    return a_gc, a_flt


def precision_weights(
    bootstrap_vectors: np.ndarray,
    floor_percentile: float = 10.0,
) -> np.ndarray:
    """Variance-based per-feature precision weights from bootstrap replicates.

    bootstrap_vectors is (B, K). Returns length-K weight vector w_k = 1 / (var_k + eps)
    where eps is the ``floor_percentile`` of the variance vector.
    """
    bootstrap_vectors = np.asarray(bootstrap_vectors, dtype=float)
    if bootstrap_vectors.ndim != 2:
        raise ValueError("bootstrap_vectors must be 2-D (B, K)")
    var_k = bootstrap_vectors.var(axis=0, ddof=1)
    if var_k.size == 0:
        return np.array([])
    eps = float(np.percentile(var_k, floor_percentile)) if var_k.size > 0 else EPS
    eps = max(eps, EPS)
    return 1.0 / (var_k + eps)


def percentile_ci(values: np.ndarray, alpha: float = 0.05) -> tuple[float, float, float]:
    """Return (median, low, high) percentile bounds for a 1-D array."""
    values = np.asarray(values, dtype=float).ravel()
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (float("nan"), float("nan"), float("nan"))
    low_pct = 100.0 * (alpha / 2.0)
    high_pct = 100.0 * (1.0 - alpha / 2.0)
    return (
        float(np.median(finite)),
        float(np.percentile(finite, low_pct)),
        float(np.percentile(finite, high_pct)),
    )


def summarize_bootstrap_decomposition(
    point: Decomposition,
    bootstrap: pd.DataFrame,
    permutation: pd.DataFrame | None = None,
    alpha: float = 0.05,
    category_config: Mapping[str, float] | None = None,
) -> pd.DataFrame:
    """Summarize a bootstrap decomposition distribution for output artifacts.

    The returned table matches the pre-registered artifact shape in
    agents_instruction.md §2.4: one row per statistic with point estimate,
    bootstrap median and percentile interval, interpretation-category fractions,
    and optional empirical permutation p-values.
    """
    cfg = dict(category_config or {})
    category_summary = categorize_bootstrap_distribution(
        bootstrap.get("beta", pd.Series(dtype=float)).to_numpy(dtype=float),
        bootstrap.get("rho", pd.Series(dtype=float)).to_numpy(dtype=float)
        if "rho" in bootstrap.columns else None,
        headline_fraction=float(cfg.get("headline_assignment_fraction", 0.70)),
        amplify_min=float(cfg.get("amplify_min_beta", 1.0)),
        redirect_abs_max=float(cfg.get("redirect_abs_beta_max", 0.10)),
        redirect_rho_min=float(cfg.get("redirect_rho_min", 0.50)),
    )
    fractions = category_summary["fractions"]
    point_values = point.as_dict()

    rows: list[dict[str, object]] = []
    for stat in ("beta", "cos", "rho", "a_flt_norm", "a_gc_norm", "residual_norm"):
        med, lo, hi = percentile_ci(
            bootstrap[stat].to_numpy(dtype=float) if stat in bootstrap.columns else np.array([]),
            alpha=alpha,
        )
        empirical_p = np.nan
        if permutation is not None and stat in permutation.columns:
            empirical_p = empirical_pvalue(float(point_values.get(stat, np.nan)), permutation[stat])
        rows.append({
            "statistic": stat,
            "weighted": bool(point.weighted),
            "point_estimate": float(point_values.get(stat, np.nan)),
            "median": med,
            "ci_low": lo,
            "ci_high": hi,
            "empirical_p": empirical_p,
            "frac_amplify": fractions["amplify"],
            "frac_dampen": fractions["dampen"],
            "frac_reverse": fractions["reverse"],
            "frac_redirect": fractions["redirect"],
            "frac_boundary": fractions["boundary"],
            "frac_undefined": fractions["undefined"],
            "interpretation": interpretation_label(category_summary),
            "n_bootstrap": int(category_summary["n_bootstrap"]),
        })
    return pd.DataFrame(rows)


def bootstrap_decomposition(
    vector_builder: Callable[[np.ndarray], tuple[np.ndarray, np.ndarray]],
    strata: Sequence,
    n_iterations: int = 2000,
    weights: np.ndarray | None = None,
    rng: np.random.Generator | None = None,
) -> pd.DataFrame:
    """Bootstrap the (beta, cos, rho) statistics.

    ``vector_builder(indices)`` must return ``(A_GC, A_FLT)`` recomputed on the
    resampled rows. ``strata`` is a per-sample stratum (e.g. concatenation of
    Age, Arm, EnvGroup labels) governing the stratified resample.
    """
    if rng is None:
        rng = np.random.default_rng()
    rows: list[dict] = []
    for b in range(int(n_iterations)):
        idx = stratified_bootstrap_indices(strata, rng)
        try:
            a_gc, a_flt = vector_builder(idx)
        except ValueError:
            rows.append({"iteration": b, "beta": np.nan, "cos": np.nan, "rho": np.nan})
            continue
        dec = compute_decomposition(a_flt, a_gc, weights=weights)
        rows.append({
            "iteration": b,
            "beta": dec.beta,
            "cos": dec.cos,
            "rho": dec.rho,
            "a_flt_norm": dec.a_flt_norm,
            "a_gc_norm": dec.a_gc_norm,
            "residual_norm": dec.residual_norm,
        })
    return pd.DataFrame(rows)


def permutation_decomposition(
    feature_matrix: np.ndarray,
    age: Sequence,
    env: Sequence,
    arm_mask: Sequence[bool] | None = None,
    weights: np.ndarray | None = None,
    n_iterations: int = 5000,
    rng: np.random.Generator | None = None,
    strata: Sequence | None = None,
    old_label: str = "OLD",
    young_label: str = "YNG",
    gc_label: str = "GC",
    flt_label: str = "FLT",
) -> pd.DataFrame:
    """Permute Old/Young labels within supplied strata and recompute stats.

    ``strata`` should exclude the permuted label itself. If omitted, the legacy
    behavior stratifies by EnvGroup plus the supplied arm mask.
    """
    if rng is None:
        rng = np.random.default_rng()
    age = np.asarray(age)
    env = np.asarray(env)
    arm_mask = (np.ones(age.shape[0], dtype=bool)
                if arm_mask is None else np.asarray(arm_mask, dtype=bool))
    if strata is None:
        strata = pd.Series([f"{e}|{int(m)}" for e, m in zip(env, arm_mask)])
    else:
        strata = pd.Series(strata).reset_index(drop=True)
        if len(strata) != len(age):
            raise ValueError("strata length must equal number of samples")

    rows: list[dict] = []
    for k in range(int(n_iterations)):
        permuted_age = stratified_permutation_indices(age, strata, rng)
        try:
            a_gc, a_flt = build_aging_vectors(
                feature_matrix, permuted_age, env, arm_mask,
                old_label=old_label, young_label=young_label,
                gc_label=gc_label, flt_label=flt_label,
            )
        except ValueError:
            rows.append({"iteration": k, "beta": np.nan, "cos": np.nan, "rho": np.nan})
            continue
        dec = compute_decomposition(a_flt, a_gc, weights=weights)
        rows.append({
            "iteration": k,
            "beta": dec.beta,
            "cos": dec.cos,
            "rho": dec.rho,
        })
    return pd.DataFrame(rows)


def empirical_pvalue(observed: float, null_distribution: Iterable[float],
                     two_sided: bool = True) -> float:
    """Empirical p-value of an observed statistic against a null distribution."""
    null_arr = np.asarray(list(null_distribution), dtype=float)
    null_arr = null_arr[np.isfinite(null_arr)]
    if null_arr.size == 0 or not np.isfinite(observed):
        return float("nan")
    if two_sided:
        extreme = np.sum(np.abs(null_arr) >= abs(observed))
    else:
        extreme = np.sum(null_arr >= observed)
    return float((extreme + 1) / (null_arr.size + 1))
