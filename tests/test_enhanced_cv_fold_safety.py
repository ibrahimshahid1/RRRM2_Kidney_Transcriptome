import numpy as np

from src.validation.enhanced_cv import build_network_features_fold


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
