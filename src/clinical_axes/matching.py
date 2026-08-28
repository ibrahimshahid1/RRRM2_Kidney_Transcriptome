"""Flight-label-blind matching helpers for compartment-marker audits."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


def robust_standardize(table: pd.DataFrame) -> pd.DataFrame:
    """Center columns by their median and scale by the normal-consistent MAD."""

    if not isinstance(table, pd.DataFrame) or table.empty:
        raise ValueError("table must be a non-empty DataFrame")
    values = table.astype(float)
    if not np.isfinite(values.to_numpy()).all():
        raise ValueError("matching covariates must all be finite")
    median = values.median(axis=0)
    scale = (values - median).abs().median(axis=0) * 1.4826
    scale = scale.mask(scale <= 1e-12, values.std(axis=0, ddof=1))
    scale = scale.mask(scale <= 1e-12, 1.0)
    return (values - median) / scale


def euclidean_distance_matrix(
    targets: np.ndarray, candidates: np.ndarray
) -> np.ndarray:
    """Return target-by-candidate Euclidean distances."""

    target_values = np.asarray(targets, dtype=float)
    candidate_values = np.asarray(candidates, dtype=float)
    if target_values.ndim != 2 or candidate_values.ndim != 2:
        raise ValueError("targets and candidates must be two-dimensional")
    if target_values.shape[1] != candidate_values.shape[1]:
        raise ValueError("targets and candidates must have the same columns")
    if len(candidate_values) < len(target_values):
        raise ValueError("one-to-one matching needs at least as many candidates")
    if not np.isfinite(target_values).all() or not np.isfinite(candidate_values).all():
        raise ValueError("matching inputs must all be finite")
    return np.sqrt(
        np.sum(
            (target_values[:, None, :] - candidate_values[None, :, :]) ** 2,
            axis=2,
        )
    )


def draw_balanced_unique_panels(
    distance: np.ndarray,
    target_covariates: np.ndarray,
    candidate_covariates: np.ndarray,
    *,
    pool_size: int,
    n_draws: int,
    seed: int,
    balance_caliper: float,
    max_attempt_multiplier: int = 100,
) -> tuple[np.ndarray, int]:
    """Draw unique matched panels and enforce aggregate covariate balance.

    Each target draws from its ``pool_size`` closest candidates. Candidate
    genes cannot repeat within a panel. A draw is retained only when the
    absolute difference between every target and candidate panel-mean
    covariate is no larger than ``balance_caliper``. Inputs should already be
    standardized; no outcome or treatment label is used here.
    """

    distances = np.asarray(distance, dtype=float)
    target = np.asarray(target_covariates, dtype=float)
    candidate = np.asarray(candidate_covariates, dtype=float)
    if distances.shape != (len(target), len(candidate)):
        raise ValueError("distance shape does not match the covariate matrices")
    if not 1 <= pool_size <= len(candidate):
        raise ValueError("pool_size is outside the candidate range")
    if n_draws < 1:
        raise ValueError("n_draws must be positive")
    if balance_caliper <= 0:
        raise ValueError("balance_caliper must be positive")

    pools = np.argsort(distances, axis=1)[:, :pool_size]
    # Match the hardest targets first to avoid greedy depletion of their pools.
    kth_distance = distances[np.arange(len(target)), pools[:, -1]]
    order = np.argsort(-kth_distance, kind="stable")
    target_mean = target.mean(axis=0)
    rng = np.random.default_rng(seed)
    accepted: list[np.ndarray] = []
    attempts = 0
    maximum = n_draws * max_attempt_multiplier
    while len(accepted) < n_draws and attempts < maximum:
        attempts += 1
        selected = np.full(len(target), -1, dtype=int)
        used = np.zeros(len(candidate), dtype=bool)
        valid = True
        for target_index in order:
            available = pools[target_index][~used[pools[target_index]]]
            if len(available) == 0:
                valid = False
                break
            candidate_index = int(available[rng.integers(len(available))])
            selected[target_index] = candidate_index
            used[candidate_index] = True
        if not valid:
            continue
        mean_difference = candidate[selected].mean(axis=0) - target_mean
        if np.max(np.abs(mean_difference)) <= balance_caliper:
            accepted.append(selected)
    if len(accepted) != n_draws:
        raise RuntimeError(
            f"only {len(accepted)} balanced panels accepted in {attempts} attempts"
        )
    return np.stack(accepted), attempts


def optimal_unique_assignment(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the minimum-total-distance one-to-one assignment."""

    distances = np.asarray(distance, dtype=float)
    if distances.ndim != 2 or distances.shape[1] < distances.shape[0]:
        raise ValueError("distance must have at least as many columns as rows")
    if not np.isfinite(distances).all():
        raise ValueError("distance must be finite")
    rows, columns = linear_sum_assignment(distances)
    if len(rows) != distances.shape[0]:
        raise RuntimeError("a complete one-to-one assignment was not found")
    return rows, columns


def cross_group_nearest_distances(
    candidate_covariates: np.ndarray, groups: np.ndarray
) -> np.ndarray:
    """Nearest-neighbour distances restricted to a different marker group."""

    values = np.asarray(candidate_covariates, dtype=float)
    labels = np.asarray(groups, dtype=object)
    if values.ndim != 2 or labels.ndim != 1 or len(values) != len(labels):
        raise ValueError("candidate covariates and groups are misaligned")
    distances = np.sqrt(
        np.sum((values[:, None, :] - values[None, :, :]) ** 2, axis=2)
    )
    invalid = labels[:, None] == labels[None, :]
    distances[invalid] = np.inf
    nearest = distances.min(axis=1)
    if not np.isfinite(nearest).all():
        raise ValueError("every candidate needs a neighbour in another group")
    return nearest
