import numpy as np

from src.networks.stability_test import (
    apply_stability_gate,
    estimate_agc_stability,
    gates_from_config,
)


def test_stable_vector_passes_gene_gate():
    rng = np.random.default_rng(7)
    age = np.array(["YNG"] * 8 + ["OLD"] * 8)
    x = np.zeros((16, 2))
    x[:8, 0] = rng.normal(0.0, 0.01, 8)
    x[8:, 0] = rng.normal(1.0, 0.01, 8)

    def builder(idx):
        rows = np.arange(len(age)) if idx is None else idx
        xb = x[rows]
        ab = age[rows]
        return xb[ab == "OLD"].mean(axis=0) - xb[ab == "YNG"].mean(axis=0)

    report = estimate_agc_stability(
        builder,
        strata=age,
        arm="ISS-T",
        resolution="gene",
        n_iterations=100,
        rng=np.random.default_rng(11),
    )
    decision = apply_stability_gate(report)
    assert report.cosine_median > 0.99
    assert decision.passed


def test_zero_vector_fails_module_gate():
    age = np.array(["YNG"] * 4 + ["OLD"] * 4)

    def builder(idx):
        return np.zeros(3)

    report = estimate_agc_stability(
        builder,
        strata=age,
        arm="LAR",
        resolution="module",
        n_iterations=10,
        rng=np.random.default_rng(3),
    )
    decision = apply_stability_gate(report)
    assert not decision.passed


def test_gates_from_yaml_config_shape():
    gates = gates_from_config({
        "resolutions": {
            "module": {"stability_median": 0.6, "stability_lower": 0.4},
            "gene": {"stability_median": 0.4, "stability_lower": 0.2},
        }
    })
    assert gates["module"] == {"median": 0.6, "lower": 0.4}
    assert gates["gene"] == {"median": 0.4, "lower": 0.2}
