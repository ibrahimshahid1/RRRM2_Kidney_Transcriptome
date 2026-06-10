import numpy as np
import pandas as pd
import pytest

from src.v11.matched_null import ols_slope, prepare_matched_pool, run_matched_null
from src.v11.rna_protein_propagation import classify_gene, direction


RUN_ROOT = "data/results/run_20260606_v11_layer_specificity"


def test_ols_slope_recovers_intercept_inclusive_slope():
    x = np.array([-2.0, -1.0, 0.0, 1.0, 2.0])
    y = 1.5 + 2.0 * x
    assert ols_slope(x, y) == pytest.approx(2.0)
    assert np.isnan(ols_slope([1.0, 1.0, 1.0], [0.0, 1.0, 2.0]))


def test_gene_layer_classification_prioritizes_protein_then_phospho():
    row = pd.Series(
        {
            "osd462_rna_effect": 0.20,
            "protein_flight_effect": 0.06,
            "phospho_parent_effect": -0.10,
        }
    )
    ternary, detail, protein_dir, phospho_dir = classify_gene(row)
    assert ternary == "rna_to_protein"
    assert detail == "rna_to_protein"
    assert protein_dir == "up"
    assert phospho_dir == "down"

    row["protein_flight_effect"] = -0.06
    row["phospho_parent_effect"] = 0.10
    ternary, detail, _, _ = classify_gene(row)
    assert ternary == "rna_to_phospho"
    assert detail == "rna_to_phospho_with_protein_discordance"

    row["phospho_parent_effect"] = -0.10
    ternary, detail, _, _ = classify_gene(row)
    assert ternary == "rna_only"
    assert detail == "rna_protein_discordant"


def test_direction_thresholds_flat_and_unobserved_values():
    assert direction(0.039, 0.04) == "flat"
    assert direction(0.041, 0.04) == "up"
    assert direction(np.nan, 0.04) == "unobserved"


def test_permuted_labels_return_non_significant_matched_null():
    rng = np.random.default_rng(0)
    n = 240
    x = np.tile(np.linspace(-1, 1, 24), 10) + rng.normal(0, 0.01, n)
    y = rng.permutation(x) + rng.normal(0, 0.25, n)
    df = pd.DataFrame(
        {
            "ENSEMBL": [f"g{i}" for i in range(n)],
            "osd462_rna_effect": x,
            "protein_flight_effect": y,
            "n_peptides": np.tile(np.arange(2, 26), 10),
            "abundance_log2": np.tile(np.linspace(1, 5, 24), 10),
        }
    )
    pool = prepare_matched_pool(df, ["osd462_rna_effect", "protein_flight_effect"])
    target_mask = np.zeros(len(pool), dtype=bool)
    target_mask[::8] = True
    result = run_matched_null(
        pool,
        target_mask,
        lambda part: ols_slope(part["osd462_rna_effect"], part["protein_flight_effect"]),
        "protein_slope",
        n_null=1000,
        rng=np.random.default_rng(100),
    )
    assert result.p_greater > 0.05
    assert result.p_two_sided > 0.05


def _apply_filter(df: pd.DataFrame, spec: str) -> pd.DataFrame:
    out = df.copy()
    for token in spec.split(";"):
        col, value = token.split("=", 1)
        out = out[out[col].astype(str).eq(value)]
    return out


def test_v11_layer_specificity_fixture_numbers_match_generated_artifacts():
    expected = pd.read_csv("tests/fixtures/v11_layer_specificity_numbers.tsv", sep="\t")
    for _, row in expected.iterrows():
        artifact = f"{RUN_ROOT}/{row['relative_path']}"
        try:
            df = pd.read_csv(artifact, sep="\t")
        except FileNotFoundError:
            pytest.fail(f"missing layer-specificity artifact: {artifact}")
        matched = _apply_filter(df, row["filters"])
        assert len(matched) == 1, f"{row['key']} matched {len(matched)} rows"
        observed = float(matched.iloc[0][row["value_column"]])
        assert observed == pytest.approx(float(row["expected"]), abs=float(row["abs_tolerance"]))


def test_layer_assignments_keep_ecm_inversion_and_dct_phospho_candidate():
    summary = pd.read_csv(f"{RUN_ROOT}/propagation/rna_protein_propagation_summary.tsv", sep="\t")
    labels = dict(zip(summary["pathway"], summary["layer_assignment"]))
    assert labels["ecm_organization"] == "protein_inverted_calibrated"
    assert labels["dct_ncc_wnk_transport"] == "RNA_to_phospho_candidate"
