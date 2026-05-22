"""Unit tests for the phenotype-anchoring matched-score layer."""
import numpy as np
import pandas as pd

from src.multiomics.phenotype_anchor import (
    channel_to_animal,
    compare_scores,
    condition_adjusted_residuals,
    per_sample_score,
    rna_sample_animal,
)


def test_channel_to_animal():
    assert channel_to_animal("FL-01") == ("Space Flight", "1")
    assert channel_to_animal("GC-12") == ("Ground Control", "12")
    assert channel_to_animal("BL-09") == ("Basal", "9")
    assert channel_to_animal("junk") is None


def test_rna_sample_animal():
    assert rna_sample_animal("RR10_KDN_WT_FLT_F1_mRNA_techrep1") == ("Space Flight", "1")
    assert rna_sample_animal("RR10_KDN_WT_GC_G12_mRNA") == ("Ground Control", "12")
    assert rna_sample_animal("RR10_KDN_WT_BSL_B9_mRNA") == ("Basal", "9")


def test_channel_and_rna_map_to_same_animal_key():
    # the matched-design assumption: FL-01 and FLT_F1 are the same animal
    assert channel_to_animal("FL-01") == rna_sample_animal("RR10_KDN_WT_FLT_F1_mRNA")
    assert channel_to_animal("GC-20") == rna_sample_animal("RR10_KDN_WT_GC_G20_mRNA")


def test_condition_adjusted_residuals_removes_group_mean():
    score = pd.Series({"a": 10.0, "b": 12.0, "c": 0.0, "d": 2.0})
    cond = pd.Series({"a": "FL", "b": "FL", "c": "GC", "d": "GC"})
    res = condition_adjusted_residuals(score, cond)
    # within each condition the residuals are centered at zero
    assert abs(res[["a", "b"]].mean()) < 1e-9
    assert abs(res[["c", "d"]].mean()) < 1e-9


def test_per_sample_score_skips_missing_features():
    mat = pd.DataFrame(
        {"s1": [1.0, 2.0, 3.0], "s2": [2.0, 4.0, 6.0], "s3": [3.0, 1.0, 2.0]},
        index=["g1", "g2", "g3"],
    )
    sc = per_sample_score(mat, ["g1", "g2", "missing"])
    assert list(sc.index) == ["s1", "s2", "s3"]
    assert np.isfinite(sc).all()


def test_compare_scores_group_concordant_underpowered():
    # flight low on both scores, ground high; no within-condition structure
    rng = np.random.default_rng(0)
    rows = {}
    cond = {}
    for i in range(8):
        rows[f"FL{i}"] = (-1.0 + rng.normal(0, 0.05), -1.0 + rng.normal(0, 0.05))
        cond[f"FL{i}"] = "Space Flight"
        rows[f"GC{i}"] = (1.0 + rng.normal(0, 0.05), 1.0 + rng.normal(0, 0.05))
        cond[f"GC{i}"] = "Ground Control"
    rna = pd.Series({k: v[0] for k, v in rows.items()})
    pho = pd.Series({k: v[1] for k, v in rows.items()})
    res = compare_scores(rna, pho, pd.Series(cond))
    assert res.n_animals == 16
    assert np.sign(res.group_rna_flt_minus_gc) == np.sign(res.group_phospho_flt_minus_gc)
    assert res.interpretation in {
        "group_level_concordant_animal_level_underpowered",
        "animal_level_link_supported",
    }


def test_compare_scores_no_link():
    rng = np.random.default_rng(1)
    cond = {f"FL{i}": "Space Flight" for i in range(8)}
    cond.update({f"GC{i}": "Ground Control" for i in range(8)})
    rna = pd.Series({k: rng.normal() for k in cond})
    # phospho high in flight, low in ground -> opposite group direction to noise rna
    pho = pd.Series({k: (1.0 if k.startswith("FL") else -1.0) for k in cond})
    res = compare_scores(rna, pho, pd.Series(cond))
    assert res.n_animals == 16
