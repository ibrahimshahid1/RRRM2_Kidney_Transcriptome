"""v11 Module 3 unit tests."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.v11.observability_audit import (
    assign_observability_strata,
    collapse_observability_to_gene,
    detectability_gradient,
    high_coverage_subset,
    merge_observability_into_pool,
    ncc_site_observability,
    propagation_with_strata,
)


def _toy_by_row():
    return pd.DataFrame({
        "gene_symbol": ["A", "A", "B", "C", "D", "E"],
        "n_channels_used": [30, 28, 30, 15, 5, 30],
        "n_peptides": [10, 5, 8, 4, 3, 12],
        "abundance_log2": [20.0, 19.5, 21.0, 18.0, 16.5, 22.0],
    })


def test_collapse_observability_weights_by_peptides():
    out = collapse_observability_to_gene(_toy_by_row())
    a = out[out["gene_symbol"].eq("A")].iloc[0]
    # n_channels_used weighted by peptides: (30*10 + 28*5) / 15 = 29.333
    assert a["n_channels_used_mean"] == pytest.approx((30 * 10 + 28 * 5) / 15, rel=1e-9)
    # missing_fraction = 1 - n_channels_used_mean / n_channels_total (=30)
    assert a["missing_fraction"] == pytest.approx(1 - (30 * 10 + 28 * 5) / 15 / 30, rel=1e-9)


def test_assign_observability_strata_falls_back_when_all_missing_equal():
    pool = pd.DataFrame({
        "abundance_log2": np.linspace(15, 22, 30),
        "n_peptides": np.tile(np.arange(2, 8), 5),
        "missing_fraction": np.zeros(30),
    })
    strata = assign_observability_strata(pool, n_missing_bins=3)
    # Falls back to 5×4 → at most 20 strata, no missing-fraction suffix in keys.
    assert "x0" not in strata.iloc[0]


def test_assign_observability_strata_extends_keys_when_missing_varies():
    pool = pd.DataFrame({
        "abundance_log2": np.linspace(15, 22, 30),
        "n_peptides": np.tile(np.arange(2, 8), 5),
        "missing_fraction": np.linspace(0, 0.5, 30),
    })
    strata = assign_observability_strata(pool, n_missing_bins=3)
    assert any("x0" in s for s in strata)


def test_detectability_gradient_separates_quantified_vs_not():
    rng = np.random.default_rng(0)
    rna = pd.DataFrame({
        "ENSEMBL": [f"g{i}" for i in range(200)],
        "rrrm2_iss_t_rna_effect": rng.normal(0, 1, 200),
    })
    # Top-half by |effect| are all protein-quantified.
    top_idx = rna["rrrm2_iss_t_rna_effect"].abs().rank() > 100
    protein_pool = pd.DataFrame({"ENSEMBL": rna.loc[top_idx, "ENSEMBL"].values})
    grad = detectability_gradient(rna, protein_pool, n_bins=5)
    # The largest-|RNA-effect| bin should be ~100% quantified.
    top_frac = grad.iloc[-1]["fraction_protein_quantified"]
    assert top_frac >= 0.95


def test_high_coverage_subset_filters_by_both_criteria():
    pool = pd.DataFrame({
        "gene_symbol": list("ABCDE"),
        "n_peptides": [2, 5, 10, 1, 3],
        "missing_fraction": [0.0, 0.05, 0.30, 0.1, 0.18],
    })
    sub = high_coverage_subset(pool, min_peptides=3, max_missing_fraction=0.2)
    assert set(sub["gene_symbol"]) == {"B", "E"}


def test_ncc_site_observability_finds_regulatory_and_control_sites():
    sites = pd.DataFrame({
        "gene_symbol": ["Slc12a3", "Slc12a3", "Slc12a3", "Slc12a3",
                        "Slc12a3", "Slc12a3", "Stk39", "Stk39",
                        "Other", "Other"],
        "site_position": ["53", "58", "65", "68", "96", "120",
                          "382", "383", "10", "11"],
        "site_kind": ["S"] * 10,
        "phospho_effect": [-0.8, -0.6, -0.7, -0.5, 0.1, 0.0, -1.0, -0.9, 0.5, 0.4],
        "phospho_se": [0.1] * 10,
        "phospho_p_value": [0.01, 0.02, 0.01, 0.05, 0.7, 0.9, 0.001, 0.002, 0.5, 0.6],
        "n_fl": [5, 5, 5, 5, 4, 3, 5, 5, 5, 5],
        "n_gc": [5, 5, 5, 5, 4, 3, 5, 5, 5, 5],
    })
    out = ncc_site_observability(sites)
    assert (out["role"] == "regulatory").sum() == 6
    assert (out["role"] == "non-regulatory control").sum() == 2
    # Regulatory sites should have higher coverage_percentile (n=10 each) than
    # the two NCC nonreg sites we included that have lower n_fl/n_gc.
    reg_med = out.loc[out["role"].eq("regulatory"), "coverage_percentile"].median()
    ctrl_med = out.loc[out["role"].eq("non-regulatory control"), "coverage_percentile"].median()
    assert reg_med >= ctrl_med


def test_propagation_with_strata_returns_expected_columns():
    rng = np.random.default_rng(1)
    n = 80
    rna = rng.normal(0, 1, n)
    pool = pd.DataFrame({
        "ENSEMBL": [f"g{i}" for i in range(n)],
        "gene_symbol": [f"G{i}" for i in range(n)],
        "gene_upper": [f"G{i}" for i in range(n)],
        "osd462_rna_effect": rna,
        "protein_flight_effect": rna + rng.normal(0, 0.2, n),  # strong concordance
        "n_peptides": rng.integers(2, 10, n),
        "abundance_log2": rng.normal(20, 1, n),
        "match_stratum": "a0_p0",
    })
    gene_sets = {"signal": [f"G{i}" for i in range(20)]}
    out = propagation_with_strata(pool, gene_sets, stratum_col="match_stratum",
                                  n_null=500, seed=0)
    assert "protein_slope" in out.columns
    assert "signed_mean_protein_by_rna" in out.columns
    row = out.iloc[0]
    assert row["protein_slope"] > 0.5  # strong concordance recovered
