"""Unit tests for Reach E LINCS/CMap appendix screen."""

import h5py
import numpy as np
import pandas as pd

from src.v11.cmap_screen import (
    build_query_genes,
    connectivity_score,
    merge_signature_metadata,
    mouse_to_human_symbol,
    score_level5_matrix,
    summarize_moa,
    summarize_perturbagens,
)


def _gene_info():
    return pd.DataFrame({
        "human_symbol": ["UPA", "UPB", "DNA", "DNB", "OTHER"],
        "human_entrez_id": ["1", "2", "3", "4", "5"],
        "is_landmark": [True, False, True, False, False],
        "is_bing": [True, True, True, True, False],
        "lincs_space": ["landmark", "bing", "landmark", "bing", "aig"],
    })


def test_mouse_to_human_symbol_prefers_ortholog_table_then_uppercase():
    assert mouse_to_human_symbol("Bmal1", {"bmal1": "ARNTL"}) == ("ARNTL", "ortholog_table")
    assert mouse_to_human_symbol("Slc12a3", {}) == ("SLC12A3", "uppercase_symbol_assumption")


def test_build_query_genes_selects_top_up_and_down_bing_genes():
    meta = pd.DataFrame({
        "gene_symbol": ["Upa", "Upb", "Other", "Dna", "Dnb"],
        "pooled_effect": [2.0, 1.5, 3.0, -2.1, -1.9],
        "z": [5.0, 4.0, 9.0, -6.0, -3.0],
        "fdr": [0.01, 0.02, 0.001, 0.01, 0.03],
    })
    q = build_query_genes(meta, _gene_info(), n_per_direction=1)
    included = q[q["included_in_query"]].copy()
    assert set(included["human_gene_symbol"]) == {"UPA", "DNA"}
    assert (included["human_gene_symbol"] == "OTHER").sum() == 0
    assert included.groupby("direction").size().to_dict() == {"down": 1, "up": 1}


def test_connectivity_score_positive_for_mimic_negative_for_reversal():
    values = np.array([
        [2.0, 1.0, -1.0, -2.0],   # up high, down low -> mimic
        [-2.0, -1.0, 1.0, 2.0],   # reversed
    ])
    scores = connectivity_score(values, np.array([0, 1]), np.array([2, 3]))
    assert scores[0] > 0
    assert scores[1] < 0


def test_score_level5_matrix_reads_signature_x_gene_layout(tmp_path):
    gctx = tmp_path / "toy.gctx"
    with h5py.File(gctx, "w") as h5:
        h5.create_dataset("0/DATA/0/matrix", data=np.array([
            [2.0, 1.0, -1.0, -2.0],
            [-2.0, -1.0, 1.0, 2.0],
        ], dtype="float32"))
        h5.create_dataset("0/META/ROW/id", data=np.array([b"1", b"2", b"3", b"4"]))
        h5.create_dataset("0/META/COL/id", data=np.array([b"sig_mimic", b"sig_reverse"]))
    query = pd.DataFrame({
        "included_in_query": [True, True, True, True],
        "direction": ["up", "up", "down", "down"],
        "human_entrez_id": ["1", "2", "3", "4"],
    })
    scored = score_level5_matrix(gctx, query, chunk_size=1).set_index("sig_id")
    assert scored.loc["sig_mimic", "connectivity_score"] > 0
    assert scored.loc["sig_reverse", "connectivity_score"] < 0


def test_metadata_merge_and_summaries():
    scores = pd.DataFrame({"sig_id": ["s1", "s2"], "connectivity_score": [-2.0, 1.0]})
    sig = pd.DataFrame({
        "sig_id": ["s1", "s2"],
        "pert_id": ["p1", "p2"],
        "pert_iname": ["drug_a", "drug_b"],
        "pert_type": ["trt_cp", "trt_cp"],
    })
    pert = pd.DataFrame({"pert_id": ["p1", "p2"], "moa": ["kinase inhibitor", "agonist"]})
    merged = merge_signature_metadata(scores, sig, pert)
    assert set(merged["connectivity_class"]) == {"reversal", "mimic"}
    ps = summarize_perturbagens(merged)
    assert ps.iloc[0]["pert_iname"] == "drug_a"
    moa = summarize_moa(merged)
    assert "kinase inhibitor" in set(moa["moa"])
