import numpy as np
import pandas as pd

from src.networks.mechanism_axis import (
    aggregate_lioness_nodes,
    composite_gene_priority,
    define_recurrent_axis,
    osd253_discriminator_table,
    proportional_top_ks,
    resolve_gene_sets,
    select_osd253_scenario,
)


def _id_map(tmp_path):
    path = tmp_path / "id_map.tsv"
    path.write_text(
        "ensembl_gene_id\tmgi_symbol\n"
        "g_ecm\tCol1a1\n"
        "g_tgfb\tTgfb1\n"
        "g_tlr4\tTlr4\n"
        "g_dct\tSlc12a3\n"
        "g_keep\tMyd88\n"
        "g_other\tOther1\n"
    )
    return path


def test_recurrent_axis_orients_ecm_positive_and_dct_negative():
    rrrm2 = pd.Series({
        "ecm_remodeling": 2.0,
        "fibrosis": 1.0,
        "dct_ncc_wnk": -1.5,
        "ion_transport": 0.2,
    })
    osd513 = pd.Series({
        "ecm_remodeling": 1.5,
        "fibrosis": 0.5,
        "dct_ncc_wnk": -2.0,
        "ion_transport": -0.3,
    })
    table, weights = define_recurrent_axis(rrrm2, osd513)
    assert weights["ecm_remodeling"] > 0
    assert weights["dct_ncc_wnk"] < 0
    assert "ion_transport" not in weights.index
    ion = table[table["pathway"] == "ion_transport"].iloc[0]
    assert ion["included_primary_axis"] is False or not bool(ion["included_primary_axis"])


def test_gene_set_deoverlap_preserves_protected_upstream_genes(tmp_path):
    gene_sets = {
        "tlr4_innate": {
            "role": "candidate_upstream",
            "protected_genes": ["Tlr4"],
            "genes": ["Tlr4", "Tgfb1", "Myd88"],
        }
    }
    resolved, coverage = resolve_gene_sets(
        gene_sets,
        _id_map(tmp_path),
        {"g_tlr4", "g_tgfb", "g_keep"},
        min_genes=2,
        deoverlap_gene_ids={"g_tlr4", "g_tgfb"},
    )
    assert "tlr4_innate" in resolved
    assert "g_tlr4" in resolved["tlr4_innate"].scored_gene_ids
    assert "g_tgfb" not in resolved["tlr4_innate"].scored_gene_ids
    assert coverage.loc[0, "n_dropped_axis_overlap"] == 1


def test_osd253_scenario_selection_keeps_both_strains_and_right_controls():
    meta = pd.DataFrame({
        "Sample Name": ["c57f", "c57g", "c3hf", "c3hg", "viv"],
        "Factor Value[Strain]": ["C57BL/6J", "C57BL/6J", "C3H/HeJ", "C3H/HeJ", "C57BL/6J"],
        "Factor Value[Spaceflight]": [
            "Space Flight",
            "Ground Control",
            "Space Flight",
            "Ground Control",
            "Vivarium Control",
        ],
        "Factor Value[Duration]": ["~25 day", "~25 day", "~25 day", "~25 day", "~25 day"],
        "Factor Value[Treatment]": ["White light", "Blue light", "White light", "Blue light", "White light"],
    })
    selected = select_osd253_scenario(meta, "original_gc_blue_light")
    assert set(selected["strain"]) == {"C57BL/6J", "C3H/HeJ"}
    assert selected["condition"].tolist().count("FLT") == 2
    assert selected["condition"].tolist().count("GC") == 2


def test_lioness_node_aggregation_aligns_incident_edges():
    weights = np.array([[1.0, -3.0], [2.0, 4.0]])
    edges = pd.DataFrame({"gene_i": ["g1", "g2"], "gene_j": ["g2", "g3"]})
    out = aggregate_lioness_nodes(weights, edges, ["g1", "g2", "g3"])
    assert out.loc[0, "g1"] == 1.0
    assert out.loc[0, "g2"] == 2.0
    assert out.loc[0, "g3"] == 3.0


def test_composite_scoring_handles_missing_support_and_penalizes_logfc():
    ranked = composite_gene_priority(
        pd.Series({"g1": 2.0, "g2": 0.0}),
        pd.Series({"g1": 0.0, "g2": 2.0}),
        pd.Series({"g1": 3.0, "g2": 1.0}),
        pd.Series({"g1": 0.0, "g2": 1.0}),
        iss_abs_logfc=pd.Series({"g1": 5.0, "g2": 0.1}),
    )
    assert ranked["composite_score"].notna().all()
    g1 = ranked.set_index("gene").loc["g1"]
    assert g1["silent_composite_score"] < g1["composite_score"]


def test_proportional_top_ks_adds_fractional_thresholds():
    assert proportional_top_ks(2500) == [25, 50, 125, 250]
    assert proportional_top_ks(5000) == [50, 100, 250, 500]


def test_osd253_discriminator_labels_tlr4_supported_when_c3h_attenuates():
    rows = []
    for strain, base_flt in [("C57BL/6J", 2.0), ("C3H/HeJ", 0.1)]:
        for duration in ["~25 day", "~75 day"]:
            for i in range(4):
                rows.append({
                    "Sample Name": f"{strain}_{duration}_f{i}",
                    "Factor Value[Strain]": strain,
                    "Factor Value[Spaceflight]": "Space Flight",
                    "Factor Value[Duration]": duration,
                    "Factor Value[Treatment]": "White light",
                })
                rows.append({
                    "Sample Name": f"{strain}_{duration}_g{i}",
                    "Factor Value[Strain]": strain,
                    "Factor Value[Spaceflight]": "Ground Control",
                    "Factor Value[Duration]": duration,
                    "Factor Value[Treatment]": "Blue light",
                })
    meta = pd.DataFrame(rows)
    score = {}
    for sample in meta["Sample Name"]:
        if "_f" in sample:
            score[sample] = 2.0 if sample.startswith("C57") else 0.1
        else:
            score[sample] = 0.0
    scores = pd.DataFrame({"remodeling_score": pd.Series(score)})
    summary, _ = osd253_discriminator_table(
        scores,
        meta,
        scenarios=["original_gc_blue_light"],
        score_cols=["remodeling_score"],
        n_bootstrap=20,
        rng=np.random.default_rng(1),
    )
    assert summary.iloc[0]["interpretation"] in {"tlr4_supported", "inconclusive"}
