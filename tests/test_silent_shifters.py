from pathlib import Path

import pandas as pd
import pytest

from src.statistics.silent_shifters import (
    build_silent_shifter_tables,
    de_path_for_contrast,
    load_de_table,
)


def test_de_is_required_and_strict_silent_counts_use_de(tmp_path: Path):
    contrast = "ISS_T_YNG_FLT_minus_GC"
    expected = de_path_for_contrast(contrast, tmp_path)
    assert expected.name == "ISS_T_YNG_FLT_vs_GC_gene_DE.tsv"
    with pytest.raises(FileNotFoundError):
        load_de_table(expected)

    rw = pd.DataFrame({
        "gene": [f"g{i}" for i in range(10)],
        "rewiring_mean": list(range(10)),
        "log2FC_shrunken": [0.0, 2.0, 0.1, 0.2, 3.0, 0.1, 0.1, 0.1, 0.1, 5.0],
        "lfc_ci_low": [-0.1, 1.8, -0.2, -0.1, 2.8, -0.2, -0.2, -0.2, -0.2, 4.8],
        "lfc_ci_high": [0.1, 2.2, 0.2, 0.25, 3.2, 0.2, 0.2, 0.2, 0.2, 5.2],
        "FDR": [0.5, 0.01, 0.5, 0.5, 0.01, 0.5, 0.5, 0.5, 0.5, 0.001],
    })
    annotated, candidates, strict, summary = build_silent_shifter_tables(rw, top_quantile=0.8)
    assert len(candidates) == 2
    assert len(strict) == 1
    assert summary["n_candidate_rewired"] == 2
    assert summary["n_strict_silent_shifters"] == 1
    assert "expected_strict_under_independence" in summary
    assert "lfc_ci_inside_fraction" in annotated.columns
    assert "bounded_lfc" in annotated.columns


def test_de_table_requires_ci_for_bounded_expression_claim(tmp_path: Path):
    de_path = tmp_path / "bad_gene_DE.tsv"
    pd.DataFrame({
        "gene": ["g1"],
        "log2FC": [0.1],
        "FDR": [0.8],
    }).to_csv(de_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="lfc_ci"):
        load_de_table(de_path)
