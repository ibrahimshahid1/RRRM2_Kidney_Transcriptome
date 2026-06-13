"""Cross-OSDR contrast-vector alignment utilities."""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Callable, Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import norm

from src.networks.contrast_vectors import cosine, stratified_bootstrap_indices, percentile_ci


@dataclass
class AlignmentResult:
    """One cosine-alignment result with recurrence verdict fields."""

    contrast: str
    point_estimate: float
    median: float
    ci_low: float
    ci_high: float
    recurrence_pass: bool
    n_bootstrap: int

    def as_dict(self) -> dict:
        return asdict(self)


@dataclass
class PermutationAlignmentResult:
    """Label-permutation null for one cosine-alignment result."""

    contrast: str
    point_estimate: float
    null_median: float
    null_ci_low: float
    null_ci_high: float
    p_greater: float
    p_less: float
    p_two_sided_abs: float
    p_directional: float
    n_permutation: int

    def as_dict(self) -> dict:
        return asdict(self)


def align_series(
    external_vector: pd.Series,
    reference_vector: pd.Series,
    min_features: int = 3,
) -> tuple[np.ndarray, np.ndarray, pd.Index]:
    """Align two feature-indexed vectors and return numeric arrays."""
    ext = external_vector.copy()
    ref = reference_vector.copy()
    ext.index = ext.index.astype(str)
    ref.index = ref.index.astype(str)
    common = ext.index.intersection(ref.index)
    if len(common) < min_features:
        raise ValueError(
            f"Need at least {min_features} shared features for cross-OSDR alignment; found {len(common)}"
        )
    return ext.loc[common].to_numpy(dtype=float), ref.loc[common].to_numpy(dtype=float), common


def cosine_alignment(
    external_vector: pd.Series,
    reference_vector: pd.Series,
    weights: pd.Series | None = None,
    min_features: int = 3,
) -> float:
    """Compute cosine(F_external, reference) on shared features."""
    ext, ref, common = align_series(external_vector, reference_vector, min_features=min_features)
    w = None
    if weights is not None:
        weights = weights.copy()
        weights.index = weights.index.astype(str)
        w = weights.reindex(common).to_numpy(dtype=float)
        if np.isnan(w).any():
            raise ValueError("weights must cover all shared alignment features")
    return float(cosine(ext, ref, weights=w))


def bootstrap_cosine_alignment(
    external_vector_builder: Callable[[np.ndarray | None], pd.Series],
    reference_vector: pd.Series,
    strata: Iterable,
    n_iterations: int,
    rng: np.random.Generator,
    contrast: str,
    thresholds: Mapping[str, object] | None = None,
) -> tuple[pd.DataFrame, AlignmentResult]:
    """Bootstrap external-study sample resampling for cosine alignment."""
    point = cosine_alignment(external_vector_builder(None), reference_vector)
    rows: list[dict[str, float | int]] = []
    strata_series = pd.Series(list(strata)).reset_index(drop=True)
    for b in range(int(n_iterations)):
        idx = stratified_bootstrap_indices(strata_series, rng)
        try:
            c = cosine_alignment(external_vector_builder(idx), reference_vector)
        except ValueError:
            c = np.nan
        rows.append({"iteration": b, "cos": c})
    boot = pd.DataFrame(rows)
    med, lo, hi = percentile_ci(boot["cos"].to_numpy(dtype=float))
    result = AlignmentResult(
        contrast=contrast,
        point_estimate=point,
        median=med,
        ci_low=lo,
        ci_high=hi,
        recurrence_pass=recurrence_verdict(point, lo, thresholds=thresholds),
        n_bootstrap=int(boot["cos"].notna().sum()),
    )
    return boot, result


def permutation_cosine_alignment(
    external_vector_builder: Callable[[np.ndarray | None], pd.Series],
    reference_vector: pd.Series,
    labels: Iterable,
    n_iterations: int,
    rng: np.random.Generator,
    contrast: str,
) -> tuple[pd.DataFrame, PermutationAlignmentResult]:
    """Build a FLT/GC label-permutation null for cosine alignment.

    ``external_vector_builder(None)`` must return the observed external
    flight vector. ``external_vector_builder(permuted_labels)`` must rebuild
    the external vector after assigning the provided permuted FLT/GC labels to
    the same samples. The reference vector is held fixed because the question
    is whether the external FLT/GC contrast aligns with the RRRM-2 direction
    more than expected under exchangeable external labels.
    """
    point = cosine_alignment(external_vector_builder(None), reference_vector)
    label_values = np.asarray(list(labels), dtype=object)
    rows: list[dict[str, float | int]] = []
    for k in range(int(n_iterations)):
        permuted = rng.permutation(label_values)
        try:
            c = cosine_alignment(external_vector_builder(permuted), reference_vector)
        except ValueError:
            c = np.nan
        rows.append({"iteration": k, "cos": c})

    perm = pd.DataFrame(rows)
    finite = perm["cos"].to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0 or not np.isfinite(point):
        med = lo = hi = p_greater = p_less = p_abs = p_directional = np.nan
    else:
        med, lo, hi = percentile_ci(finite)
        p_greater = float((np.sum(finite >= point) + 1.0) / (finite.size + 1.0))
        p_less = float((np.sum(finite <= point) + 1.0) / (finite.size + 1.0))
        p_abs = float((np.sum(np.abs(finite) >= abs(point)) + 1.0) / (finite.size + 1.0))
        p_directional = p_greater if point >= 0 else p_less

    result = PermutationAlignmentResult(
        contrast=contrast,
        point_estimate=point,
        null_median=float(med),
        null_ci_low=float(lo),
        null_ci_high=float(hi),
        p_greater=float(p_greater),
        p_less=float(p_less),
        p_two_sided_abs=float(p_abs),
        p_directional=float(p_directional),
        n_permutation=int(finite.size),
    )
    return perm, result


def recurrence_verdict(
    point_estimate: float,
    ci_low: float,
    same_sign: bool = True,
    thresholds: Mapping[str, object] | None = None,
) -> bool:
    """Apply the pre-registered cross-OSDR recurrence rule (§4.3)."""
    thresholds = dict(thresholds or {})
    point_min = float(thresholds.get("cosine_point_estimate_min", 0.20))
    ci_min = float(thresholds.get("cosine_ci_lower_must_exceed", 0.0))
    require_same_sign = bool(thresholds.get("require_same_sign_as_rrrm2", True))
    return bool(
        np.isfinite(point_estimate)
        and np.isfinite(ci_low)
        and point_estimate > point_min
        and ci_low > ci_min
        and (same_sign or not require_same_sign)
    )


def signed_stouffer_z(z_scores: Iterable[float], weights: Iterable[float] | None = None) -> dict[str, float]:
    """Combine signed study-level Z scores with optional precision weights."""
    z_raw = np.asarray(list(z_scores), dtype=float)
    mask = np.isfinite(z_raw)
    z = z_raw[mask]
    if z.size == 0:
        return {"z": np.nan, "p_two_sided": np.nan, "n": 0}
    if weights is None:
        w = np.ones_like(z)
    else:
        w_raw = np.asarray(list(weights), dtype=float)
        if w_raw.shape != z_raw.shape:
            raise ValueError("weights must have the same length as z_scores")
        w = w_raw[mask]
        finite_w = np.isfinite(w)
        z = z[finite_w]
        w = w[finite_w]
        if z.size == 0:
            return {"z": np.nan, "p_two_sided": np.nan, "n": 0}
    combined = float(np.sum(w * z) / np.sqrt(np.sum(w * w)))
    p = float(2.0 * norm.sf(abs(combined)))
    return {"z": combined, "p_two_sided": p, "n": int(z.size)}
