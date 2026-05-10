import numpy as np
import pandas as pd

from src.validation.osd_external_validation import (
    feature_gene_ids,
    gsea_enrichment_score,
    infer_flt_gc_samples,
    pathway_permutation_table,
)


def test_external_sample_inference_and_feature_sets():
    cols = ["A_FLT_1", "A_FLT_2", "A_GC_1", "A_GC_2", "A_VIV_1", "A_C3H_FLT_1"]
    flt, gc = infer_flt_gc_samples(cols, r"_FLT_", r"_GC_", exclude_pattern=r"C3H")
    assert flt == ["A_FLT_1", "A_FLT_2"]
    assert gc == ["A_GC_1", "A_GC_2"]

    id_map = pd.DataFrame({
        "ensembl_gene_id": ["ENS1", "ENS2", "ENS3", "ENS4", "ENS5", "ENS6", "ENS7", "ENS8"],
        "mgi_symbol": ["Ppara", "Acox1", "Cpt1a", "Rpl3", "Rps4x", "Mrpl1", "Col1a1", "Slc12a3"],
    })
    assert feature_gene_ids("PPAR_signaling", id_map, set(id_map["ensembl_gene_id"])) == ["ENS1", "ENS2", "ENS3"]
    assert feature_gene_ids("translation_machinery", id_map, set(id_map["ensembl_gene_id"])) == ["ENS4", "ENS5", "ENS6"]
    assert feature_gene_ids("ECM_remodeling", id_map, set(id_map["ensembl_gene_id"])) == ["ENS7"]
    assert feature_gene_ids("tubular_ion_transport", id_map, set(id_map["ensembl_gene_id"])) == ["ENS8"]


def test_external_pathway_permutation_table_is_finite():
    genes = ["ENS1", "ENS2", "ENS3", "ENS4", "ENS5", "ENS6"]
    samples = ["S_FLT_1", "S_FLT_2", "S_FLT_3", "S_GC_1", "S_GC_2", "S_GC_3"]
    values = np.array([
        [6, 6, 7, 2, 2, 2],
        [5, 6, 6, 2, 1, 2],
        [7, 7, 6, 3, 3, 2],
        [4, 4, 4, 4, 4, 4],
        [4, 5, 4, 4, 4, 4],
        [5, 4, 4, 4, 5, 4],
    ], dtype=float)
    vst = pd.DataFrame(values, index=genes, columns=samples)
    id_map = pd.DataFrame({
        "ensembl_gene_id": genes,
        "mgi_symbol": ["Ppara", "Acox1", "Cpt1a", "Rpl3", "Rps4x", "Mrpl1"],
    })
    registry = pd.DataFrame({
        "feature": ["PPAR_signaling", "translation_machinery"],
        "discovery_effect": [1.0, 1.0],
    })
    result = pathway_permutation_table(
        vst,
        samples[:3],
        samples[3:],
        registry,
        id_map,
        k_perm=25,
        seed=0,
    )
    assert set(result["feature"]) == {"PPAR_signaling", "translation_machinery"}
    assert result["p_value"].between(1 / 26, 1).all()
    assert result["q_value"].between(0, 1).all()
    assert (result["pathway_method"] == "gsea").all()


def test_gsea_enrichment_score_tracks_rank_direction():
    scores = np.array([3.0, 2.0, -1.0, -2.0])
    top_hits = np.array([True, True, False, False])
    bottom_hits = np.array([False, False, True, True])

    assert gsea_enrichment_score(scores, top_hits) > 0
    assert gsea_enrichment_score(scores, bottom_hits) < 0
