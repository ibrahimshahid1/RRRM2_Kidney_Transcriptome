from pathlib import Path

import numpy as np
import pandas as pd

from src.networks.lar_reversal import (
    STATE_COMPONENT_FEATURES,
    STATE_FEATURES,
    classify_vector_relationship,
    cosine_similarity,
    interaction_table,
    projection_beta,
    reversal_summary_for_features,
    score_symbol_panels,
    vector_stats,
)


def test_state_vector_features_exclude_derived_scalar():
    assert "matrix_minus_dct" not in STATE_FEATURES
    assert "matrix_minus_dct" in STATE_COMPONENT_FEATURES


def test_cosine_and_projection_detect_exact_reversal():
    iss = pd.Series({"a": 1.0, "b": 2.0, "c": -1.0})
    lar = -iss
    assert np.isclose(cosine_similarity(lar, iss), -1.0)
    assert np.isclose(projection_beta(lar, iss), -1.0)
    row = vector_stats("synthetic", "pooled", iss, lar).__dict__
    row["bootstrap_ci_high"] = -0.9
    assert classify_vector_relationship(row) == "model_B_reversal_candidate"


def test_classification_detects_attenuation():
    iss = pd.Series({"a": 3.0, "b": -2.0, "c": 1.0})
    lar = pd.Series({"a": 0.05, "b": -0.02, "c": 0.01})
    row = vector_stats("synthetic", "pooled", iss, lar).__dict__
    row["bootstrap_ci_high"] = 0.5
    assert classify_vector_relationship(row) == "model_A_attenuation_candidate"


def test_reversal_summary_bootstrap_on_synthetic_rows():
    rows = []
    for arm, sign in [("ISS-T", 1.0), ("LAR", -1.0)]:
        for age in ["YNG", "OLD"]:
            for condition, shift in [("GC", 0.0), ("FLT", sign)]:
                for rep in range(4):
                    rows.append({
                        "Sample Name": f"{arm}_{age}_{condition}_{rep}",
                        "study": "RRRM-2",
                        "scenario": "primary",
                        "Arm": arm,
                        "Age": age,
                        "condition": condition,
                        "f1": shift + 0.01 * rep,
                        "f2": 2 * shift + 0.01 * rep,
                    })
    df = pd.DataFrame(rows)
    summary, external = reversal_summary_for_features(
        df,
        ["f1", "f2"],
        "synthetic",
        age_scope="pooled",
        n_bootstrap=20,
        n_permutation=20,
        rng=np.random.default_rng(1),
        include_osd513=False,
    )
    assert external.empty
    assert summary.loc[0, "cos_lar_iss"] < -0.99
    assert summary.loc[0, "beta_lar_to_iss"] < -0.99


def test_lar_opposes_requires_stable_lar_nonzero_effect():
    rows = []
    for arm, flight_shift in [("ISS-T", 2.0), ("LAR", -0.05)]:
        for age in ["YNG", "OLD"]:
            for condition in ["GC", "FLT"]:
                for rep in range(6):
                    rows.append({
                        "Sample Name": f"{arm}_{age}_{condition}_{rep}",
                        "study": "RRRM-2",
                        "Arm": arm,
                        "Age": age,
                        "condition": condition,
                        "score": (flight_shift if condition == "FLT" else 0.0) + (0.2 * (-1) ** rep if arm == "LAR" else 0.01 * rep),
                    })
    out = interaction_table(
        pd.DataFrame(rows),
        ["score"],
        feature_family="synthetic",
        n_bootstrap=100,
        n_permutation=100,
        rng=np.random.default_rng(7),
    )
    assert out.loc[0, "iss_t_effect"] > 0
    assert out.loc[0, "lar_effect"] < 0
    assert not bool(out.loc[0, "lar_opposes_iss"])


def test_score_symbol_panels_handles_missing_genes(tmp_path: Path):
    expr = pd.DataFrame(
        {
            "s1": [1.0, 2.0],
            "s2": [2.0, 4.0],
            "s3": [3.0, 6.0],
        },
        index=["g1", "g2"],
    )
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("ensembl_gene_id\tmgi_symbol\ng1\tPer2\ng2\tRbm3\n")
    scores, coverage = score_symbol_panels(
        expr,
        str(id_map),
        panels={"clock_test": ("Per2", "MissingGene")},
    )
    assert "clock_test" in scores.columns
    assert "Rbm3_expression" in scores.columns
    assert coverage.loc[coverage["query_symbol"].eq("MissingGene"), "n_resolved_in_expression"].iloc[0] == 0
