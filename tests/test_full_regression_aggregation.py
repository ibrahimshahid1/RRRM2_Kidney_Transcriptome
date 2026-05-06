import numpy as np
import pandas as pd

from src.statistics.full_regression import (
    build_factorial_design,
    empirical_signed_gene_pvalues,
    fit_edge_regression,
    signed_incident_t_aggregate,
)


def test_signed_empirical_aggregation_has_no_zero_saturation_under_null():
    rng = np.random.default_rng(123)
    rows = []
    for age in ["YNG", "OLD"]:
        for arm in ["LAR", "ISS-T"]:
            for env in ["GC", "FLT"]:
                for rep in range(3):
                    rows.append({"Age": age, "Arm": arm, "EnvGroup": env, "sample": f"{age}_{arm}_{env}_{rep}"})
    meta = pd.DataFrame(rows)
    W = rng.normal(0, 1, size=(len(meta), 6))
    edge_i = np.array([0, 0, 1, 1, 2, 2], dtype=np.int32)
    edge_j = np.array([1, 2, 2, 3, 3, 0], dtype=np.int32)

    X, terms = build_factorial_design(meta)
    edge_results = fit_edge_regression(W, X, terms)
    agg = signed_incident_t_aggregate(edge_results, edge_i, edge_j, 4, "flight")
    pvals = empirical_signed_gene_pvalues(
        W,
        meta,
        edge_i,
        edge_j,
        4,
        ["flight"],
        {"flight": agg["signed_t_sum_sqrt_degree"]},
        n_perm=25,
        seed=7,
    )["flight"]
    assert np.nanmin(pvals) >= 1 / 26
    assert np.all(pvals[np.isfinite(pvals)] > 0)
    assert np.nanmedian(pvals) > 0.05
