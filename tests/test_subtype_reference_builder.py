"""Tests for flight-blind distal-nephron reference construction."""

from pathlib import Path

import pandas as pd
import pytest
import yaml

from src.subtype_reference.reference_builder import (
    _assert_reference_only,
    _broad_expression_table,
    build_signatures,
    load_config,
)


CONFIG = Path("config/dct_subtype_reference_freeze_v1.yaml")


def _write(frame: pd.DataFrame, path: Path) -> Path:
    frame.to_csv(path, sep="\t", index=False)
    return path


def test_reference_config_is_flight_blind_and_asdn_matches_production():
    cfg = load_config(CONFIG)
    categories = cfg["reference_builder"]["curated_sets"]["ASDN"]
    reference = {gene for genes in categories.values() for gene in genes}
    production = yaml.safe_load(
        Path("config/dct_asdn_phospho_reanalysis.yaml").read_text()
    )
    assert reference == set(production["asdn_gene_set"]["genes"])
    assert len(reference) == 29


def test_target_absent_unrelated_specific_gene_is_not_broad():
    cfg = load_config(CONFIG)
    atlas = pd.DataFrame(
        {
            "gene_symbol": ["UnrelatedSpecific", "UnrelatedSpecific"],
            "compartment": ["DCT_general", "proximal_tubule"],
            "mean_cpm": [0.0, 50.0],
            "source_study_detection_fraction": [1.0, 1.0],
        }
    )
    result = _broad_expression_table(atlas, cfg).set_index("gene_symbol")
    assert not bool(result.loc["UnrelatedSpecific", "target_expressed"])
    assert not bool(result.loc["UnrelatedSpecific", "broadly_expressed"])
    assert not bool(
        result.loc["UnrelatedSpecific", "non_distal_specific_or_broad"]
    )


def test_true_breadth_is_separate_from_single_compartment_nonspecificity():
    cfg = load_config(CONFIG)
    rows = []
    for gene, unrelated_values in {
        "SingleOffTarget": [30.0, 0.1, 0.1, 0.1],
        "MultiCompartment": [10.0, 10.0, 10.0, 10.0],
    }.items():
        rows.append(
            {
                "gene_symbol": gene,
                "compartment": "DCT_general",
                "mean_cpm": 20.0,
                "source_study_detection_fraction": 1.0,
            }
        )
        for compartment, cpm in zip(
            ("proximal_tubule", "endothelial", "immune", "fibroblast"),
            unrelated_values,
        ):
            rows.append(
                {
                    "gene_symbol": gene,
                    "compartment": compartment,
                    "mean_cpm": cpm,
                    "source_study_detection_fraction": 1.0,
                }
            )
    result = _broad_expression_table(pd.DataFrame(rows), cfg).set_index(
        "gene_symbol"
    )
    assert not bool(result.loc["SingleOffTarget", "broadly_expressed"])
    assert bool(
        result.loc["SingleOffTarget", "non_distal_specific_or_broad"]
    )
    assert bool(result.loc["MultiCompartment", "broadly_expressed"])
    assert bool(
        result.loc["MultiCompartment", "non_distal_specific_or_broad"]
    )


def test_builder_constructs_distinct_reference_sets(tmp_path):
    genes = ["Dct1Gene", "Dct2Peak", "Transition", "BroadDct2"]
    discovery = pd.DataFrame(
        {
            "gene_symbol": genes,
            "log2_fc_dct2_vs_dct1": [-2.0, 2.0, 1.5, 2.0],
            "fdr": [0.001] * 4,
            "pct_detected_dct1": [0.8, 0.2, 0.2, 0.2],
            "pct_detected_dct2": [0.2, 0.8, 0.8, 0.8],
            "n_consistent_pairs": [3] * 4,
            "n_pairs": [3] * 4,
        }
    )
    fine = pd.DataFrame(
        {
            "gene_symbol": genes,
            "log2_fc_dct2_vs_dct1": [-1.5, 1.5, 1.2, 1.5],
            "fdr_dct2_vs_dct1": [0.01] * 4,
            "log2_fc_dct2_vs_cnt": [-1.0, 0.8, -0.2, 0.8],
            "fdr_dct2_vs_cnt": [0.01] * 4,
            "log2_fc_dct1_vs_cnt": [1.5, -1.0, -1.0, -1.0],
            "fdr_dct1_vs_cnt": [0.01] * 4,
            "pct_detected_dct1": [0.8, 0.2, 0.2, 0.2],
            "pct_detected_dct2": [0.2, 0.8, 0.8, 0.8],
            "pct_detected_cnt": [0.2, 0.3, 0.8, 0.3],
            "n_consistent_dct2_vs_dct1": [3] * 4,
            "n_pairs": [3] * 4,
        }
    )
    segment = pd.DataFrame(
        {
            "gene_symbol": genes,
            "log2_fc_cnt_vs_dct": [-1.5, -1.5, -0.1, -0.1],
            "ci_low_cnt_vs_dct": [-2.0, -2.0, -0.4, -0.4],
            "p_noninferiority_cnt_vs_dct": [0.5, 0.5, 0.001, 0.001],
            "fdr_noninferiority_cnt_vs_dct": [0.5, 0.5, 0.01, 0.01],
            "pct_replicates_detected_cnt": [1.0] * 4,
            "n_consistent_cnt_retained": [0, 0, 3, 3],
            "n_pairs": [3] * 4,
        }
    )
    atlas_rows = []
    for gene in genes:
        atlas_rows.extend(
            [
                {
                    "gene_symbol": gene,
                    "compartment": "DCT_general",
                    "mean_cpm": 20.0,
                    "source_study_detection_fraction": 1.0,
                },
                {
                    "gene_symbol": gene,
                    "compartment": "proximal_tubule",
                    "mean_cpm": 0.1 if gene != "BroadDct2" else 30.0,
                    "source_study_detection_fraction": 1.0,
                },
                {
                    "gene_symbol": gene,
                    "compartment": "endothelial",
                    "mean_cpm": 0.1,
                    "source_study_detection_fraction": 1.0,
                },
                {
                    "gene_symbol": gene,
                    "compartment": "immune",
                    "mean_cpm": 0.1,
                    "source_study_detection_fraction": 1.0,
                },
                {
                    "gene_symbol": gene,
                    "compartment": "fibroblast",
                    "mean_cpm": 0.1,
                    "source_study_detection_fraction": 1.0,
                },
            ]
        )
    atlas = pd.DataFrame(atlas_rows)

    out = tmp_path / "out"
    gate = build_signatures(
        CONFIG,
        out,
        _write(discovery, tmp_path / "discovery.tsv"),
        _write(fine, tmp_path / "fine.tsv"),
        _write(segment, tmp_path / "segment.tsv"),
        _write(atlas, tmp_path / "atlas.tsv"),
    )
    membership = pd.read_csv(out / "data_driven_gene_set_membership.tsv", sep="\t")
    observed = set(zip(membership["gene_symbol"], membership["gene_set"]))
    assert ("Dct1Gene", "DCT1_core") in observed
    assert ("Dct2Peak", "strict_DCT2_peaked") in observed
    assert ("Transition", "DCT2_CNT_transition") in observed
    assert not any(gene == "BroadDct2" for gene, _ in observed)
    flags = pd.read_csv(out / "broad_expression_flags.tsv", sep="\t").set_index(
        "gene_symbol"
    )
    assert not bool(flags.loc["BroadDct2", "broadly_expressed"])
    assert bool(flags.loc["BroadDct2", "non_distal_specific_or_broad"])
    assert gate["transition_set_final"] is True
    assert gate["fine_subtype_sets_final"] is True


def test_transition_floor_and_external_signature_independence(tmp_path):
    discovery = pd.DataFrame(
        {
            "gene_symbol": ["Retained", "WeakRetention"],
            "log2_fc_dct2_vs_dct1": [1.5, 1.5],
            "fdr": [0.001, 0.001],
            "pct_detected_dct1": [0.2, 0.2],
            "pct_detected_dct2": [0.8, 0.8],
            "n_consistent_pairs": [3, 3],
            "n_pairs": [3, 3],
        }
    )
    # ExternalOnly is deliberately absent from discovery. Its membership must
    # therefore be independent of GSE228367 discovery membership.
    validation = pd.DataFrame(
        {
            "gene_symbol": ["Retained", "WeakRetention", "ExternalOnly"],
            "log2_fc_dct2_vs_dct1": [1.0, 1.0, 1.0],
            "fdr_dct2_vs_dct1": [0.01, 0.01, 0.01],
            "log2_fc_dct2_vs_cnt": [0.2, 0.2, 0.2],
            "fdr_dct2_vs_cnt": [0.5, 0.5, 0.5],
            "log2_fc_dct1_vs_cnt": [-0.8, -0.8, -0.8],
            "fdr_dct1_vs_cnt": [0.01, 0.01, 0.01],
            "pct_detected_dct1": [0.2, 0.2, 0.2],
            "pct_detected_dct2": [0.8, 0.8, 0.8],
            "pct_detected_cnt": [0.8, 0.8, 0.8],
            "n_consistent_dct2_vs_dct1": [3, 3, 3],
            "n_pairs": [3, 3, 3],
        }
    )
    segment = pd.DataFrame(
        {
            "gene_symbol": ["Retained", "WeakRetention", "ExternalOnly"],
            "log2_fc_cnt_vs_dct": [-0.1, -0.75, -0.1],
            "ci_low_cnt_vs_dct": [-0.4, -0.9, -0.4],
            "p_noninferiority_cnt_vs_dct": [0.001, 0.001, 0.001],
            "fdr_noninferiority_cnt_vs_dct": [0.01, 0.01, 0.01],
            "pct_replicates_detected_cnt": [1.0, 1.0, 1.0],
            "n_consistent_cnt_retained": [3, 3, 3],
            "n_pairs": [3, 3, 3],
        }
    )
    atlas_rows = []
    for gene in ["Retained", "WeakRetention", "ExternalOnly"]:
        for compartment, cpm in (
            ("DCT_general", 20.0),
            ("proximal_tubule", 0.1),
            ("endothelial", 0.1),
            ("immune", 0.1),
            ("fibroblast", 0.1),
        ):
            atlas_rows.append(
                {
                    "gene_symbol": gene,
                    "compartment": compartment,
                    "mean_cpm": cpm,
                    "source_study_detection_fraction": 1.0,
                }
            )
    out = tmp_path / "out"
    gate = build_signatures(
        CONFIG,
        out,
        _write(discovery, tmp_path / "discovery.tsv"),
        _write(validation, tmp_path / "validation.tsv"),
        _write(segment, tmp_path / "segment.tsv"),
        _write(pd.DataFrame(atlas_rows), tmp_path / "atlas.tsv"),
    )
    membership = pd.read_csv(out / "frozen_gene_sets.tsv", sep="\t")
    transition = set(
        membership.loc[
            membership["gene_set"].eq("DCT2_CNT_transition"), "gene_symbol"
        ]
    )
    external = set(
        membership.loc[
            membership["gene_set"].eq("DCT2_CNT_external_validation"),
            "gene_symbol",
        ]
    )
    assert "Retained" in transition
    assert "WeakRetention" not in transition
    assert "ExternalOnly" in external
    external_rows = membership[
        membership["gene_set"].eq("DCT2_CNT_external_validation")
    ]
    assert external_rows["analysis_role"].eq("external_validation").all()
    assert gate["external_validation_set_final"] is True


def test_missing_reference_data_emits_curated_sets_but_no_data_driven_claim(tmp_path):
    gate = build_signatures(CONFIG, tmp_path)
    frozen = pd.read_csv(tmp_path / "frozen_gene_sets.tsv", sep="\t")
    assert set(frozen.loc[frozen["gene_set"] == "ASDN", "gene_symbol"])
    assert not (tmp_path / "data_driven_gene_set_membership.tsv").read_text().strip().endswith(
        "DCT2_CNT_transition"
    )
    assert gate["transition_set_final"] is False
    assert gate["fine_subtype_sets_final"] is False


def test_outcome_leaking_reference_input_is_refused(tmp_path):
    path = tmp_path / "flight_phospho.tsv"
    path.write_text("gene_symbol\tvalue\nX\t1\n")
    with pytest.raises(ValueError, match="outcome-leaking"):
        _assert_reference_only(path, pd.DataFrame({"gene_symbol": ["X"]}))
