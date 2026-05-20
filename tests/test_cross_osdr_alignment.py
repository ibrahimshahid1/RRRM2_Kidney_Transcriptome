import numpy as np
import pandas as pd

from src.networks.cross_osdr_projection import (
    cosine_alignment,
    permutation_cosine_alignment,
    recurrence_verdict,
    signed_stouffer_z,
)


def test_external_cosine_alignment_on_shared_features():
    external = pd.Series({"g1": 1.0, "g2": 1.0, "g3": 99.0})
    reference = pd.Series({"g1": 1.0, "g2": 1.0})
    assert np.isclose(cosine_alignment(external, reference, min_features=2), 1.0)


def test_recurrence_thresholds_are_pre_registered():
    assert recurrence_verdict(0.25, 0.01, same_sign=True)
    assert not recurrence_verdict(0.19, 0.01, same_sign=True)
    assert not recurrence_verdict(0.25, -0.01, same_sign=True)
    assert not recurrence_verdict(0.25, 0.01, same_sign=False)


def test_signed_stouffer_z_combines_finite_scores():
    result = signed_stouffer_z([1.0, np.nan, 2.0], weights=[1.0, 1.0, 1.0])
    assert result["n"] == 2
    assert result["z"] > 2.0
    assert 0.0 < result["p_two_sided"] < 1.0


def test_permutation_cosine_alignment_reports_tail_pvalues():
    feature_matrix = pd.DataFrame(
        {
            "p1": [2.0, 2.2, -1.0, -1.1],
            "p2": [1.0, 1.1, -0.5, -0.4],
            "p3": [1.6, 1.7, -0.8, -0.7],
        },
        index=["s1", "s2", "s3", "s4"],
    )
    samples = pd.DataFrame({"sample": feature_matrix.index})
    labels = np.array(["FLT", "FLT", "GC", "GC"], dtype=object)
    reference = pd.Series({"p1": 1.0, "p2": 1.0, "p3": 1.0})

    def builder(permuted_labels):
        rows = samples.copy()
        rows["condition"] = labels if permuted_labels is None else permuted_labels
        flt = rows.loc[rows["condition"].eq("FLT"), "sample"]
        gc = rows.loc[rows["condition"].eq("GC"), "sample"]
        return feature_matrix.loc[flt].mean(axis=0) - feature_matrix.loc[gc].mean(axis=0)

    perm, summary = permutation_cosine_alignment(
        builder,
        reference,
        labels=labels,
        n_iterations=50,
        rng=np.random.default_rng(1),
        contrast="toy",
    )
    assert len(perm) == 50
    assert summary.point_estimate > 0.95
    assert 0.0 < summary.p_greater <= 1.0
    assert 0.0 < summary.p_two_sided_abs <= 1.0
