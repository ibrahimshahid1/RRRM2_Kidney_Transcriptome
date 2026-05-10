import numpy as np

from src.validation.enhanced_cv import build_network_features_fold, run_classifier


def test_sparse_edge_selection_uses_training_variance_only():
    rng = np.random.default_rng(42)
    train_w = rng.normal(size=(6, 120))
    train_w[:, 0] = 0.0
    test_w = rng.normal(size=(1, 120))
    # Test-only huge value on edge 0 must not affect feature selection.
    test_w[0, 0] = 1_000_000
    edge_i = np.arange(120) % 10
    edge_j = (np.arange(120) + 1) % 10
    x_train, x_test = build_network_features_fold(
        "sparse_edges",
        train_w,
        test_w,
        edge_i.astype(np.int32),
        edge_j.astype(np.int32),
        n_genes=10,
        pathway_masks={},
    )
    assert x_train.shape[1] == 100
    assert 1_000_000 not in set(x_test.ravel())


def test_run_classifier_returns_one_prediction_per_test_sample():
    x_train = np.array([
        [-2.0, -1.8],
        [-1.7, -2.1],
        [-1.9, -1.6],
        [1.8, 2.0],
        [2.2, 1.7],
        [1.6, 2.1],
    ])
    y_train = np.array([0, 0, 0, 1, 1, 1])
    x_test = np.array([
        [-1.5, -1.7],
        [1.4, 1.8],
        [2.4, 2.0],
    ])

    pred, prob = run_classifier(x_train, y_train, x_test, "LogisticRegression")

    assert pred.shape == (3,)
    assert prob.shape == (3,)
    assert np.isfinite(prob).all()
