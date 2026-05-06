from pathlib import Path

import numpy as np
import pandas as pd

from src.statistics.permutation_bootstrap import (
    edge_sum_node_rewiring_abs,
    hierarchical_fdr_tables,
    load_candidate_mask,
    load_hierarchical_gene_sets,
)


def test_edge_sum_statistic_and_candidate_mask_are_prespecified(tmp_path):
    genes = ["ENS1", "ENS2", "ENS3"]
    delta = np.array([1.0, -2.0])
    rew = edge_sum_node_rewiring_abs(delta, np.array([0, 1]), np.array([1, 2]), 3)
    assert rew.tolist() == [1.0, 3.0, 2.0]

    cand = tmp_path / "candidates.txt"
    cand.write_text("ENS1\nENS3\n")
    mask, loaded = load_candidate_mask(cand, genes)
    assert loaded == {"ENS1", "ENS3"}
    assert mask.tolist() == [True, False, True]


def test_hierarchical_gene_sets_and_overlap_conservative_columns(tmp_path):
    genes = ["ENS1", "ENS2", "ENS3", "ENS4"]
    id_map = tmp_path / "id_map.tsv"
    pd.DataFrame({
        "ensembl_gene_id": genes,
        "mgi_symbol": ["A", "B", "C", "D"],
    }).to_csv(id_map, sep="\t", index=False)
    gene_sets = tmp_path / "sets.yaml"
    gene_sets.write_text(
        "set_one:\n  genes: [A, B, C]\n"
        "set_two:\n  genes: [B, C, D]\n"
    )
    masks = load_hierarchical_gene_sets(gene_sets, id_map, genes, min_size=3)
    assert sorted(masks) == ["curated::set_one", "curated::set_two"]

    rew_obs = np.array([4.0, 3.0, 2.0, 1.0])
    null = np.vstack([rew_obs[::-1], np.ones(4), np.arange(4)])
    p_gene = np.array([0.01, 0.2, 0.3, 0.4])
    fam, gene_df = hierarchical_fdr_tables(rew_obs, null, p_gene, genes, masks)
    assert "q_family_BY_overlap_conservative" in fam.columns
    assert "q_BY_all_family_gene_rows_overlap_conservative" in gene_df.columns


def test_hierarchical_gene_q_values_only_confirm_selected_families():
    genes = ["ENS1", "ENS2", "ENS3", "ENS4", "ENS5", "ENS6"]
    masks = {
        "selected_family": np.array([True, True, True, False, False, False]),
        "unselected_family": np.array([False, False, False, True, True, True]),
    }
    rew_obs = np.array([10.0, 9.0, 8.0, 0.5, 0.4, 0.3])
    null = np.zeros((30, 6))
    p_gene = np.array([0.001, 0.01, 0.02, 0.001, 0.01, 0.02])

    fam, gene_df = hierarchical_fdr_tables(
        rew_obs,
        null,
        p_gene,
        genes,
        masks,
        alpha=0.10,
    )

    assert set(fam.loc[fam["q_family_BH"] < 0.10, "gene_set"]) == {"selected_family"}
    selected_rows = gene_df[gene_df["gene_set"] == "selected_family"]
    unselected_rows = gene_df[gene_df["gene_set"] == "unselected_family"]
    assert selected_rows["confirmatory_gene_tested"].all()
    assert selected_rows["q_BB_two_stage"].notna().all()
    assert not unselected_rows["confirmatory_gene_tested"].any()
    assert unselected_rows["q_BB_two_stage"].isna().all()
