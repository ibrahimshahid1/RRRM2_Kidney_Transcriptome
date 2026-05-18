import numpy as np
import pandas as pd
import pytest

from src.networks.edge_regression import validate_covariate_policy
from src.networks.lioness import compute_lioness_weights
from src.networks.shared_topology import union_skeleton_with_priors
from src.validation.cross_validation import transform_fold_edge_weights


def test_default_lioness_does_not_fisher_z_sample_specific_weights_outside_correlation_bounds():
    X = np.array([
        [10.1, -0.2, 0.7, 0.1, -0.5],
        [-9.9, 0.3, -0.8, 0.2, 0.6],
        [0.5, -0.1, 0.2, -0.4, 0.1],
    ])
    edge_i = np.array([0], dtype=np.int32)
    edge_j = np.array([1], dtype=np.int32)

    raw, meta = compute_lioness_weights(X, edge_i, edge_j, transform="raw")
    assert raw.min() < -1.0
    assert meta["sample_specific_values_are_correlations"] is False

    ranknorm, meta_rank = compute_lioness_weights(X, edge_i, edge_j, transform="raw_ranknorm")
    assert np.isfinite(ranknorm).all()
    assert meta_rank["fisher_z_used_on_sample_specific_weights"] is False


def test_z_contribution_is_explicit_sensitivity_mode():
    X = np.array([
        [10.1, -0.2, 0.7, 0.1, -0.5],
        [-9.9, 0.3, -0.8, 0.2, 0.6],
        [0.5, -0.1, 0.2, -0.4, 0.1],
    ])
    z, meta = compute_lioness_weights(
        X,
        np.array([0], dtype=np.int32),
        np.array([1], dtype=np.int32),
        transform="z_contribution",
    )
    assert np.isfinite(z).all()
    assert meta["fisher_z_used_on_sample_specific_weights"] is True


def test_fold_transform_uses_training_reference_for_test_values():
    train = np.array([[0.0, 1.0], [1.0, 2.0], [2.0, 3.0], [3.0, 4.0]])
    test = np.array([[1000.0, -1000.0]])
    train_z, test_z = transform_fold_edge_weights(train, test, transform="raw_ranknorm")
    assert np.isfinite(train_z).all()
    assert np.isfinite(test_z).all()
    assert test_z[0, 0] > train_z[:, 0].max()
    assert test_z[0, 1] < train_z[:, 1].min()


def test_residualized_edge_regression_rejects_nuisance_covariates():
    with pytest.raises(ValueError, match="residualized"):
        validate_covariate_policy("residualized", ["LibraryBatch"])
    validate_covariate_policy("non_residualized", ["LibraryBatch"])


def test_skeleton_union_labels_prior_and_data_driven_edges():
    genes = ["g1", "g2", "g3"]
    pc = np.array([
        [0.0, 0.4, -0.2],
        [0.4, 0.0, 0.1],
        [-0.2, 0.1, 0.0],
    ])
    prior = pd.DataFrame([
        {
            "gene_i": "g2",
            "gene_j": "g3",
            "edge_origin": "external_prior",
            "edge_source_detail": "kidney_prior",
            "is_fixed_prior": True,
        },
        {
            "gene_i": "g1",
            "gene_j": "g2",
            "edge_origin": "preregistered_pathway",
            "edge_source_detail": "dct_ncc_wnk",
            "is_fixed_prior": True,
        },
    ])

    edge_df, edge_i, edge_j = union_skeleton_with_priors(
        genes,
        pc,
        np.array([0], dtype=np.int32),
        np.array([1], dtype=np.int32),
        prior,
    )

    assert len(edge_df) == 2
    assert edge_i.tolist() == [0, 1]
    assert edge_j.tolist() == [1, 2]
    g1_g2 = edge_df[(edge_df["gene_i"] == "g1") & (edge_df["gene_j"] == "g2")].iloc[0]
    assert g1_g2["is_fixed_prior"] is True or bool(g1_g2["is_fixed_prior"])
    assert "osd771_data" in g1_g2["edge_origin"]
    assert "preregistered_pathway" in g1_g2["edge_origin"]
