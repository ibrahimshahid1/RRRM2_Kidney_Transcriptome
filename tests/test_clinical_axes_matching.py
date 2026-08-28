from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.clinical_axes.matching import (
    cross_group_nearest_distances,
    draw_balanced_unique_panels,
    euclidean_distance_matrix,
    optimal_unique_assignment,
    robust_standardize,
)


def test_robust_standardize_is_finite_and_median_centered():
    table = pd.DataFrame(
        {"a": [1.0, 2.0, 3.0, 40.0], "b": [2.0, 4.0, 6.0, 8.0]}
    )
    standardized = robust_standardize(table)
    assert np.isfinite(standardized.to_numpy()).all()
    assert standardized.median().to_numpy() == pytest.approx([0.0, 0.0])


def test_unique_panel_draws_respect_pool_and_balance_constraints():
    target = np.array([[0.0], [1.0], [2.0]])
    candidate = np.array([[-0.1], [0.1], [0.9], [1.1], [1.9], [2.1]])
    distance = euclidean_distance_matrix(target, candidate)
    draws, attempts = draw_balanced_unique_panels(
        distance,
        target,
        candidate,
        pool_size=2,
        n_draws=30,
        seed=19,
        balance_caliper=0.15,
    )
    assert draws.shape == (30, 3)
    assert attempts >= len(draws)
    assert all(len(np.unique(draw)) == 3 for draw in draws)
    differences = candidate[draws].mean(axis=1) - target.mean(axis=0)
    assert np.max(np.abs(differences)) <= 0.15


def test_optimal_assignment_is_unique_and_cross_group_distances_exclude_group():
    target = np.array([[0.0], [10.0]])
    candidate = np.array([[0.1], [9.9], [30.0]])
    distance = euclidean_distance_matrix(target, candidate)
    rows, columns = optimal_unique_assignment(distance)
    assert rows.tolist() == [0, 1]
    assert columns.tolist() == [0, 1]

    nearest = cross_group_nearest_distances(
        np.array([[0.0], [0.2], [2.0], [2.3]]),
        np.array(["a", "a", "b", "b"]),
    )
    assert nearest.tolist() == pytest.approx([2.0, 1.8, 1.8, 2.1])
