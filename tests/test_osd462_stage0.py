"""Tests for the OSD-462 Stage 0 provenance audit."""
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from src.multiomics.osd462_stage0 import (
    PHOSPHO_WORKBOOK,
    PROTEIN_WORKBOOK,
    audit_ncc_spak_phosphosites,
    build_sample_design,
    build_stage0_qc,
    derive_site_components,
    isolated_canonical_assay_features,
    map_sequence_phosphoform,
    motif_centre_residue,
)
from src.multiomics.phenotype_anchor import (
    NCC_REGULATORY_SITES,
    OSD462_COMODIFIED_CANONICAL_INDEX_FEATURES,
    OSD462_ISOLATED_CANONICAL_ASSAY_FEATURES,
    RENAL_AXIS_LITERATURE_CANONICAL_SITES,
)


SOURCE_WORKBOOKS_PRESENT = PROTEIN_WORKBOOK.exists() and PHOSPHO_WORKBOOK.exists()


def test_motif_centre_residue_and_component_mapping():
    assert motif_centre_residue("STLYMRTFGYNTI") == "T"
    assert motif_centre_residue("IDVVPAYEHYANS") == "Y"
    assert derive_site_components(
        "53;65", "STLYMRTFGYNTI;IDVVPAYEHYANS"
    ) == [
        (53, "T", "STLYMRTFGYNTI"),
        (65, "Y", "IDVVPAYEHYANS"),
    ]


def test_component_mapping_rejects_position_motif_mismatch():
    with pytest.raises(ValueError, match="position/motif component mismatch"):
        derive_site_components("53;65", "STLYMRTFGYNTI")


def test_peptide_phosphoform_mapping_uses_motif_anchor():
    ncc = map_sequence_phosphoform(
        "53",
        "STLYMRTFGYNTI",
        "R.T#FGYNTIDVVPAY#EHYANSALPGEPR.K",
    )
    assert ncc["sequence_phosphoform_mapping_status"] == "mapped"
    assert ncc["sequence_phosphoform_site_labels"] == "T53;Y65"
    spak = map_sequence_phosphoform(
        "383",
        "RRVPGSSGHLHKT",
        "R.VPGS#S#GHLHK.T",
    )
    assert spak["sequence_phosphoform_site_labels"] == "S382;S383"


def test_literature_sites_are_separate_from_qualified_assay_features():
    literature = set(RENAL_AXIS_LITERATURE_CANONICAL_SITES)
    assert literature == {
        ("Slc12a3", "53"),
        ("Slc12a3", "58"),
        ("Slc12a3", "71"),
        ("Stk39", "243"),
        ("Stk39", "383"),
    }
    assert OSD462_ISOLATED_CANONICAL_ASSAY_FEATURES == ()
    assert NCC_REGULATORY_SITES == ()
    assert set(OSD462_COMODIFIED_CANONICAL_INDEX_FEATURES) == {
        ("Slc12a3", "53"),
        ("Stk39", "383"),
    }


def test_curated_renal_net_is_residue_explicit_and_strict():
    path = Path("data/external/kinase_substrate/renal_kinase_substrate_core.tsv")
    net = pd.read_csv(path, sep="\t", dtype={"substrate_site": str})
    assert {
        "kinase",
        "substrate_gene",
        "substrate_site",
        "substrate_residue",
        "annotation_status",
        "reference_scope",
        "evidence_id",
    }.issubset(net.columns)
    assert set(net["annotation_status"]) == {"strict_canonical"}
    assert set(net["reference_scope"]) == {
        "literature_definition_not_assay_qualification"
    }
    assert set(net["substrate_residue"]) <= {"S", "T"}
    ncc = net[net["substrate_gene"].eq("Slc12a3")]
    assert set(zip(ncc["substrate_site"], ncc["substrate_residue"])) == {
        ("53", "T"),
        ("58", "T"),
        ("71", "S"),
    }
    spak = net[net["substrate_gene"].eq("Stk39")]
    assert set(zip(spak["substrate_site"], spak["substrate_residue"])) == {
        ("243", "T"),
        ("383", "S"),
    }
    assert not net["substrate_site"].isin(["65", "68", "382", "339"]).any()


@pytest.mark.skipif(
    not SOURCE_WORKBOOKS_PRESENT,
    reason="OSD-462 source workbooks are not available",
)
def test_source_sample_design_and_raw_inventory():
    design, raw = build_sample_design()
    assert len(design) == 60
    assert len(raw) == 96
    assert set(design.groupby("modality").size()) == {30}
    assert set(
        design.groupby(["modality", "plex", "condition_code"]).size()
    ) == {5}
    assert design["reporter_tag_matches_isa"].all()
    assert design["workbook_sample_matches_isa_source"].all()
    assert design["condition_perfectly_confounded_with_reporter_tag"].all()
    assert not design["reporter_assignment_swapped_between_plexes"].any()
    assert design["raw_files_shared_across_plex"].all()
    assert not design["raw_file_sample_specific"].any()
    assert set(raw.groupby(["modality", "plex"]).size()) == {24}
    assert not raw.duplicated(["modality", "raw_acquisition_id"]).any()

    duplicated_aliases = raw[
        raw.duplicated(["modality", "plex", "raw_alias"], keep=False)
    ]
    observed = set(
        zip(
            duplicated_aliases["modality"],
            duplicated_aliases["plex"],
            duplicated_aliases["raw_alias"],
        )
    )
    assert observed == {
        ("protein", "Samp6-10", "tc-885_11"),
        ("protein", "Samp6-10", "tc-885_12"),
    }


@pytest.mark.skipif(
    not SOURCE_WORKBOOKS_PRESENT,
    reason="OSD-462 source workbooks are not available",
)
def test_source_workbook_resolves_ncc_y65_y68_as_noncanonical():
    audit = audit_ncc_spak_phosphosites()
    assert len(audit) == 21

    def row(gene: str, position: str, kind: str = "single_site_rollup"):
        match = audit[
            audit["gene_symbol"].eq(gene)
            & audit["workbook_site_position"].eq(position)
            & audit["site_feature_kind"].eq(kind)
        ]
        assert len(match) == 1
        return match.iloc[0]

    assert row("Slc12a3", "53")["residue_resolved_site_labels"] == "T53"
    assert (
        row("Slc12a3", "53")["sequence_phosphoform_site_labels"]
        == "T53;Y65"
    )
    assert (
        row("Slc12a3", "53")["canonical_claim_class"]
        == "canonical_site_indexed_but_comodified_context_only"
    )
    assert bool(row("Slc12a3", "53")["all_components_strict_canonical"])
    assert not bool(row("Slc12a3", "53")["isolated_canonical_assay_feature"])
    assert bool(
        row("Slc12a3", "53")[
            "sequence_has_additional_modifications_beyond_feature_components"
        ]
    )
    assert row("Slc12a3", "65")["residue_resolved_site_labels"] == "Y65"
    assert row("Slc12a3", "68")["residue_resolved_site_labels"] == "Y68"
    assert not bool(row("Slc12a3", "65")["any_component_strict_canonical"])
    assert not bool(row("Slc12a3", "68")["any_component_strict_canonical"])

    y_composite = row("Slc12a3", "65;68", "composite_site_rollup")
    assert y_composite["residue_resolved_site_labels"] == "Y65;Y68"
    assert y_composite["canonical_claim_class"] == "not_in_strict_canonical_set"
    mixed = row("Slc12a3", "53;65", "composite_site_rollup")
    assert (
        mixed["canonical_claim_class"]
        == "mixed_composite_not_attributable_to_canonical_site"
    )
    assert row("Stk39", "382")["residue_resolved_site_labels"] == "S382"
    assert not bool(row("Stk39", "382")["all_components_strict_canonical"])
    assert row("Stk39", "383")["residue_resolved_site_labels"] == "S383"
    assert (
        row("Stk39", "383")["sequence_phosphoform_site_labels"]
        == "S382;S383"
    )
    assert (
        row("Stk39", "383")["canonical_claim_class"]
        == "canonical_site_indexed_but_comodified_context_only"
    )
    assert bool(row("Stk39", "383")["all_components_strict_canonical"])
    assert not bool(row("Stk39", "383")["isolated_canonical_assay_feature"])
    assert bool(
        row("Stk39", "383")[
            "sequence_has_additional_modifications_beyond_feature_components"
        ]
    )
    assert isolated_canonical_assay_features(audit).empty


@pytest.mark.skipif(
    not SOURCE_WORKBOOKS_PRESENT,
    reason="OSD-462 source workbooks are not available",
)
def test_stage0_qc_has_no_failures_and_expected_warnings():
    design, raw = build_sample_design()
    audit = audit_ncc_spak_phosphosites()
    isolated = isolated_canonical_assay_features(audit)
    qc_table = build_stage0_qc(design, raw, audit)
    qc = qc_table.set_index("metric")
    assert not qc["status"].eq("FAIL").any()
    assert qc.loc["raw_alias_uniqueness", "status"] == "WARN"
    assert qc.loc["condition_reporter_tag_assignment", "status"] == "WARN"
    assert qc.loc["detailed_protocol_TMTpro_evidence", "status"] == "PASS"
    assert qc.loc["legacy_investigation_iTRAQ_term", "status"] == "WARN"
    assert qc.loc["canonical_feature_sequence_co_modification", "status"] == "WARN"
    assert qc.loc["isolated_canonical_assay_coverage", "status"] == "WARN"
    assert qc.loc["NCC_Y65_Y68_annotation", "status"] == "PASS"
    assert qc.loc["peptide_phosphoform_absolute_site_mapping", "status"] == "PASS"

    # Exact source-derived output contract. These hashes intentionally fail if
    # workbook parsing, row order, fields, or QC semantics drift.
    outputs = {
        "sample_design": design,
        "raw_inventory": raw,
        "phosphosite_audit": audit,
        "isolated_canonical_features": isolated,
        "qc": qc_table,
    }
    assert {name: len(frame) for name, frame in outputs.items()} == {
        "sample_design": 60,
        "raw_inventory": 96,
        "phosphosite_audit": 21,
        "isolated_canonical_features": 0,
        "qc": 17,
    }
    hashes = {
        name: hashlib.sha256(
            frame.to_csv(sep="\t", index=False).encode("utf-8")
        ).hexdigest()
        for name, frame in outputs.items()
    }
    assert hashes == {
        "sample_design": "8997d9aba9db054c4456cd86ce7fbb64fe852cd80e5110b694e4266e1abbd536",
        "raw_inventory": "d119a5b787340c7740c533ac5a72e2a29fbb91d72b037bdce6e3f1849e29fc37",
        "phosphosite_audit": "65eacf725c4b7649258a3a26099584a399589fd785c3b67b2add726ebbc64b45",
        "isolated_canonical_features": "2843fb8a9796ac74e169758537216f389760da0fa662b9f01cf5b0e7535dc83c",
        "qc": "e0fbd13c79526b2d3b308b0070c6bd91b19f967a88df2e490972da24b10bf84c",
    }
