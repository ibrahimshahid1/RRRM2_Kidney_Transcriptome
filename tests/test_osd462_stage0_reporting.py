from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    REPO_ROOT / "scripts" / "osd462" / "09_stage0_manuscript_reporting.py"
)
STAGE0_DIR = (
    REPO_ROOT / "data" / "results" / "run_20260728_osd462_stage0"
)


def _load_reporting_module():
    spec = importlib.util.spec_from_file_location("stage0_reporting", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_stage0_reporting_tables_preserve_frozen_contract():
    reporting = _load_reporting_module()
    design = pd.read_csv(
        STAGE0_DIR / reporting.DESIGN_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    provenance = pd.read_csv(
        STAGE0_DIR / reporting.PROVENANCE_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )

    inclusion = reporting.build_inclusion_table(design)
    assert len(inclusion) == 30
    assert inclusion["primary_FL_vs_GC_included"].sum() == 20
    assert inclusion["condition_code"].value_counts().to_dict() == {
        "BL": 10,
        "FL": 10,
        "GC": 10,
    }
    assert inclusion["condition_tag_block_aliased"].all()
    assert not inclusion["cross_plex_label_swap"].any()

    compact = reporting.build_compact_phosphoform_table(provenance)
    assert len(compact) == len(provenance) == 21
    assert not compact["isolated_canonical_qualified"].eq("Yes").any()

    ncc_t53 = compact.loc[
        compact["gene_symbol"].eq("Slc12a3")
        & compact["feature_type"].eq("single_site_rollup")
        & compact["residue_resolved_feature"].eq("T53")
    ]
    assert len(ncc_t53) == 1
    assert ncc_t53.iloc[0]["observed_peptide_phosphoform"] == "T53;Y65"

    spak_s383 = compact.loc[
        compact["gene_symbol"].eq("Stk39")
        & compact["feature_type"].eq("single_site_rollup")
        & compact["residue_resolved_feature"].eq("S383")
    ]
    assert len(spak_s383) == 1
    assert spak_s383.iloc[0]["observed_peptide_phosphoform"] == "S382;S383"


def test_reporter_map_renders_all_publication_formats(tmp_path):
    reporting = _load_reporting_module()
    design = pd.read_csv(
        STAGE0_DIR / reporting.DESIGN_FILE,
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    inclusion = reporting.build_inclusion_table(design)
    outputs = reporting.plot_reporter_tag_map(inclusion, tmp_path, dpi=120)
    assert {path.suffix for path in outputs} == {".png", ".pdf", ".svg"}
    assert all(path.exists() and path.stat().st_size > 1_000 for path in outputs)
