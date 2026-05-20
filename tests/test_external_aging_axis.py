import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from src.networks.external_aging_axis import (
    build_tms_kidney_female_aging_axis,
    build_rrrm2_flight_vector,
    bootstrap_tms_aging_axis,
    load_tms_donor_axis_source,
    load_aging_axis,
    load_or_build_aging_axis,
    project_vector_onto_axis,
)


def test_external_axis_loads_and_projects(tmp_path):
    axis_path = tmp_path / "axis.tsv"
    pd.DataFrame({
        "ensembl_gene_id": ["g1", "g2"],
        "old_vs_young_logfc": [1.0, 0.0],
    }).to_csv(axis_path, sep="\t", index=False)
    axis = load_aging_axis(axis_path, "ensembl_gene_id", "old_vs_young_logfc")
    dec, per_feature = project_vector_onto_axis(
        pd.Series({"g1": 2.0, "g2": 0.0}),
        axis,
        min_features=2,
    )
    assert dec.beta == 2.0
    assert dec.cos == 1.0
    assert per_feature["feature"].tolist() == ["g1", "g2"]


def test_rrrm2_age_pooled_flight_vector():
    features = pd.DataFrame({
        "g1": [0.0, 1.0, 3.0, 4.0],
        "g2": [0.0, 0.0, 1.0, 1.0],
    }, index=["gc_y", "gc_o", "flt_y", "flt_o"])
    meta = pd.DataFrame({
        "Sample Name": ["gc_y", "gc_o", "flt_y", "flt_o"],
        "Arm": ["ISS-T"] * 4,
        "EnvGroup": ["GC", "GC", "FLT", "FLT"],
        "Age": ["YNG", "OLD", "YNG", "OLD"],
    })
    flight = build_rrrm2_flight_vector(features, meta, arm="ISS-T")
    assert np.allclose(flight.loc[["g1", "g2"]], [3.0, 1.0])


def test_builds_tms_axis_from_local_h5ad_shape(tmp_path):
    h5ad_path = tmp_path / "tiny_tms.h5ad"
    axis_path = tmp_path / "axis.tsv"
    adata = ad.AnnData(
        X=sparse.csr_matrix([
            [100.0, 10.0],
            [120.0, 10.0],
            [500.0, 10.0],
            [520.0, 10.0],
        ]),
        obs=pd.DataFrame({
            "age": ["3m", "3m", "18m", "18m"],
            "donor_id": ["young1", "young2", "old1", "old2"],
            "sex": ["female"] * 4,
        }, index=["c1", "c2", "c3", "c4"]),
        var=pd.DataFrame(index=["g1", "g2"]),
    )
    adata.write_h5ad(h5ad_path)

    axis_df = build_tms_kidney_female_aging_axis(h5ad_path, axis_path)
    assert axis_path.exists()
    donor_expr, donor_meta = load_tms_donor_axis_source(axis_path)
    boot_axis = bootstrap_tms_aging_axis(donor_expr, donor_meta, np.random.default_rng(1))
    assert {"ensembl_gene_id", "old_vs_young_logfc", "se", "p_value", "fdr"} <= set(axis_df.columns)
    assert axis_df.loc[axis_df["ensembl_gene_id"] == "g1", "old_vs_young_logfc"].item() > 0
    assert boot_axis.index.tolist() == ["g1", "g2"]

    axis = load_or_build_aging_axis(axis_path, "ensembl_gene_id", "old_vs_young_logfc")
    assert axis.index.tolist() == ["g1", "g2"]
