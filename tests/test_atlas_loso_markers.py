"""Tests for item 8: leave-one-atlas-source-study-out marker-tier reconstruction.

Everything except the clearly marked data-contract tests uses synthetic data.  The
data-contract tests reproduce already-published, flight-blind atlas artefacts and are
skipped when the files are absent; no test reads any mission effect on real data.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
import yaml

from scripts.clinical_axes import run_atlas_loso_markers as rl
from src.clinical_axes import atlas_loso as al
from src.clinical_axes.data import MissionData, cpm_eligible_genes
from src.clinical_axes.statistics import random_effects_reml_mkh
from src.v13.compartment_adversarial_audit import build_marker_tiers, load_config

REPO = Path(__file__).resolve().parents[1]
ATLAS_DIR = REPO / "data/processed/v13_dct_reference/mouse_kidney_atlas"
COUNTS = ATLAS_DIR / "atlas_pseudobulk_counts.tsv.gz"
SAMPLE_META = ATLAS_DIR / "atlas_pseudobulk_sample_metadata.tsv"
EXPRESSION = ATLAS_DIR / "atlas_compartment_expression.tsv.gz"
QC_JSON = ATLAS_DIR / "atlas_pseudobulk_qc.json"
FROZEN_TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
KEGG = REPO / "data/processed/gene_sets/KEGG_2019_Mouse.json"
AUDIT_CONFIG = REPO / "config/v13_compartment_adversarial_audit.yaml"
FREEZE_CONFIG = REPO / "config/dct_subtype_reference_freeze_v1.yaml"
SCRIPT_03 = REPO / "scripts/subtype_reference/03_atlas_pseudobulk.py"

real_atlas = pytest.mark.skipif(
    not all(p.exists() for p in (COUNTS, SAMPLE_META, EXPRESSION, FROZEN_TIERS, KEGG)),
    reason="atlas pseudobulk / frozen tiers / KEGG not present",
)


# --- synthetic atlas -----------------------------------------------------------------

STUDIES = ["S1", "S2", "S3", "S4", "S5"]
LAYOUT = {
    "podocyte": ["S1", "S2", "S3", "S4"],
    "endothelial": ["S1", "S2", "S3", "S4", "S5"],
    "DCT2_CNT_atlas": ["S5"],
}
MARKERS = {
    "podocyte": [f"Pa{i}" for i in range(10)],
    "endothelial": [f"Ea{i}" for i in range(10)],
    "DCT2_CNT_atlas": [f"Da{i}" for i in range(10)],
}
BACKGROUND = [f"Bg{i}" for i in range(20)]
# Pedge is podocyte-specific in S1-S3 only: detection fraction 3/4 = 0.75 in the full
# atlas (high-specificity), 2/3 after omitting S1, S2 or S3 (drops out), 3/3 after S4.
PEDGE = "Pedge"
SYN_GENES = [g for v in MARKERS.values() for g in v] + [PEDGE] + BACKGROUND
POD_HS = MARKERS["podocyte"] + [PEDGE]
STRUCTURAL_TERMS = {"T1": ["pa0"] + [f"bg{i}" for i in range(10)]}
MAPPING = {
    "podocyte": "podocyte",
    "endothelial": "endothelial",
    "DCT2_CNT_atlas": "DCT2_CNT_context",
}
TIER_NAMES = [
    "all_enriched",
    "within_kidney_not_broad",
    "high_specificity",
    "broad_enriched",
    "scaffold_excluded",
]
CONTROL = "broad_structural_scaffold_control__all"


def synthetic_pseudobulk() -> tuple[pd.DataFrame, pd.DataFrame]:
    columns: dict[str, pd.Series] = {}
    meta_rows = []
    for study in STUDIES:
        for compartment, contributors in LAYOUT.items():
            if study not in contributors:
                continue
            counts = pd.Series(0, index=SYN_GENES, dtype=np.int64)
            counts[MARKERS[compartment]] = 200
            counts[BACKGROUND] = 50
            if compartment == "podocyte" and study in ("S1", "S2", "S3"):
                counts[PEDGE] = 200
            sample = f"{study}::{compartment}"
            columns[sample] = counts
            meta_rows.append((sample, 100, study, compartment, True))
    # an excluded (<25 cells) pseudobulk, absent from the count matrix like in script 03
    meta_rows.append(("S5::podocyte", 10, "S5", "podocyte", False))
    counts = pd.DataFrame(columns)
    counts.index.name = "gene_symbol"
    meta = pd.DataFrame(
        meta_rows,
        columns=["sample_id", "n_cells", "source_study_id", "compartment", "included"],
    ).sort_values("sample_id", ignore_index=True)
    counts = counts[[s for s in meta["sample_id"] if s in counts.columns]]
    return counts, meta


def synthetic_audit_cfg() -> dict:
    return {
        "marker_tiers": {
            "target_mean_cpm_min": 1.0,
            "target_source_study_detection_fraction_min": 0.5,
            "all_enriched_log2_target_to_max_other_min": 1.0,
            "within_kidney_broad_non_target_count_min": 2,
            "within_kidney_broad_non_target_mean_cpm_min": 1.0,
            "within_kidney_broad_non_target_detection_fraction_min": 0.5,
            "high_specificity_log2_target_to_max_other_min": 2.0,
            "high_specificity_target_detection_fraction_min": 0.75,
        },
        "compartment_mapping": dict(MAPPING),
        "structural_scaffold_control": {
            "gene_set_name": CONTROL,
            "source_terms": ["T1"],
        },
        "set_test": {
            "primary_family": [
                f"{report}__{tier}" for report in MAPPING.values() for tier in TIER_NAMES
            ]
            + [CONTROL],
            "minimum_observable_genes": 8,
        },
    }


def synthetic_qc_cfg(minimum: int = 3) -> dict:
    return {
        "minimum_source_studies_per_compartment": minimum,
        "required_compartments": ["podocyte", "endothelial"],
    }


# --- hand calculation, exclusion semantics ---------------------------------------------


def test_expression_matches_hand_calculation():
    counts = pd.DataFrame(
        {
            "S1::A": [10, 90, 900],
            "S2::A": [0, 50, 50],
            "S1::B": [1, 1, 0],
            "S2::B": [3, 0, 1],
        },
        index=pd.Index(["g1", "g2", "g3"], name="gene_symbol"),
    )
    meta = pd.DataFrame(
        {
            "sample_id": ["S1::A", "S2::A", "S3::A", "S1::B", "S2::B"],
            "n_cells": [30, 40, 5, 50, 60],
            "source_study_id": ["S1", "S2", "S3", "S1", "S2"],
            "compartment": ["A", "A", "A", "B", "B"],
            "included": [True, True, False, True, True],
        }
    )
    out = al.compartment_expression_from_pseudobulk(counts, meta, min_cpm=1.0)
    assert list(out.columns) == list(al.EXPRESSION_COLUMNS)
    a = out[out["compartment"] == "A"].set_index("gene_symbol")
    b = out[out["compartment"] == "B"].set_index("gene_symbol")
    # CPM of S1::A = (1e4, 9e4, 9e5); S2::A = (0, 5e5, 5e5)
    np.testing.assert_allclose(a["mean_cpm"], [5_000.0, 295_000.0, 700_000.0])
    np.testing.assert_allclose(a["median_cpm"], [5_000.0, 295_000.0, 700_000.0])
    np.testing.assert_allclose(a["source_study_detection_fraction"], [0.5, 1.0, 1.0])
    assert (a["n_source_studies"] == 2).all()
    # S1::B = (5e5, 5e5, 0); S2::B = (7.5e5, 0, 2.5e5)
    np.testing.assert_allclose(b["mean_cpm"], [625_000.0, 250_000.0, 125_000.0])
    np.testing.assert_allclose(b["source_study_detection_fraction"], [1.0, 0.5, 0.5])
    assert (out["reference_label"] == "mouse_kidney_atlas").all()
    # the excluded (included=False) S3::A pseudobulk never contributes
    assert list(out["compartment"].unique()) == ["A", "B"]


def test_median_uses_all_pseudobulks_of_the_compartment():
    counts = pd.DataFrame(
        {"S1::A": [1, 1], "S2::A": [1, 3], "S3::A": [1, 9]},
        index=pd.Index(["g1", "g2"], name="gene_symbol"),
    )
    meta = pd.DataFrame(
        {
            "sample_id": list(counts.columns),
            "source_study_id": ["S1", "S2", "S3"],
            "compartment": ["A"] * 3,
        }
    )
    out = al.compartment_expression_from_pseudobulk(counts, meta).set_index("gene_symbol")
    # g2 CPM: 5e5, 7.5e5, 9e5 -> median 7.5e5, mean 7.1666...e5
    assert out.loc["g2", "median_cpm"] == pytest.approx(750_000.0)
    assert out.loc["g2", "mean_cpm"] == pytest.approx((500_000 + 750_000 + 900_000) / 3)
    assert out.loc["g1", "source_study_detection_fraction"] == 1.0


def test_exclusion_removes_exactly_that_studys_columns():
    counts, meta = synthetic_pseudobulk()
    full = al.compartment_expression_from_pseudobulk(counts, meta)
    dropped = al.compartment_expression_from_pseudobulk(counts, meta, exclude_sources=["S2"])
    keep = meta[meta["source_study_id"] != "S2"]
    manual = al.compartment_expression_from_pseudobulk(
        counts.drop(columns=[c for c in counts.columns if c.startswith("S2::")]), keep
    )
    pd.testing.assert_frame_equal(dropped, manual)
    n_full = full.groupby("compartment")["n_source_studies"].first()
    n_dropped = dropped.groupby("compartment")["n_source_studies"].first()
    assert (n_full - n_dropped).to_dict() == {"endothelial": 1, "podocyte": 1, "DCT2_CNT_atlas": 0}
    # excluding nothing is the identity on an explicit empty list too
    pd.testing.assert_frame_equal(
        full, al.compartment_expression_from_pseudobulk(counts, meta, exclude_sources=())
    )


def test_excluding_the_only_contributor_removes_the_compartment_and_unknown_raises():
    counts, meta = synthetic_pseudobulk()
    out = al.compartment_expression_from_pseudobulk(counts, meta, exclude_sources=["S5"])
    assert "DCT2_CNT_atlas" not in set(out["compartment"])
    with pytest.raises(ValueError, match="unknown source studies"):
        al.compartment_expression_from_pseudobulk(counts, meta, exclude_sources=["NOPE"])
    with pytest.raises(ValueError, match="no pseudobulks remain"):
        al.compartment_expression_from_pseudobulk(counts, meta, exclude_sources=STUDIES)


def test_missing_counts_column_is_rejected():
    counts, meta = synthetic_pseudobulk()
    with pytest.raises(ValueError, match="lack samples"):
        al.compartment_expression_from_pseudobulk(counts.drop(columns=["S1::podocyte"]), meta)


def test_exclusion_equals_script_03_summary_for_the_same_synthetic_input(tmp_path):
    """Excluding nothing reproduces 03_atlas_pseudobulk.py's compartment summary."""

    ad = pytest.importorskip("anndata")
    import scipy.sparse as sp

    rng = np.random.default_rng(7)
    genes = [f"G{i}" for i in range(40)]
    rows, labels = [], []
    plan = {
        "S1": {"PT": 30, "PODO": 30, "ENDO": 30, "weird": 5},
        "S2": {"PT": 40, "PODO": 26, "ENDO": 30},
        "S3": {"PT": 35, "PODO": 30, "ENDO": 33},
        "S4": {"PT": 60, "PODO": 10, "ENDO": 28},
    }
    base = {"PT": rng.gamma(2.0, 3.0, 40), "PODO": rng.gamma(2.0, 3.0, 40), "ENDO": rng.gamma(2.0, 3.0, 40)}
    for study, cell_types in plan.items():
        for cell_type, n in cell_types.items():
            mean = base.get(cell_type, np.full(40, 2.0))
            rows.append(rng.poisson(mean, size=(n, 40)))
            labels += [(study, cell_type)] * n
    matrix = np.vstack(rows).astype(np.float32)
    obs = pd.DataFrame(labels, columns=["Origin", "Celltype_finest"])
    obs["Celltype"] = obs["Celltype_finest"]
    obs.index = [f"c{i}" for i in range(len(obs))]
    var = pd.DataFrame(index=genes)
    adata = ad.AnnData(X=sp.csr_matrix(matrix), obs=obs, var=var)
    adata.raw = adata
    h5 = tmp_path / "synthetic.h5ad"
    adata.write_h5ad(h5)
    config = {
        "reference_builder": {
            "whole_kidney_compartment_aliases": {
                "proximal_tubule": ["^pt$"],
                "podocyte": ["^podo$"],
                "endothelial": ["^endo$"],
            },
            "thresholds": {"min_cells_per_pseudobulk": 25},
            "broad_expression": {"min_cpm": 1.0},
            "atlas_qc": {
                "minimum_mapping_fraction": 0.5,
                "minimum_source_studies_per_compartment": 3,
                "required_compartments": ["proximal_tubule", "endothelial"],
            },
        }
    }
    config_path = tmp_path / "freeze.yaml"
    config_path.write_text(yaml.safe_dump(config))
    out = tmp_path / "pseudobulk"
    subprocess.run(
        [
            sys.executable, str(SCRIPT_03), "--config", str(config_path), "--input", str(h5),
            "--output-dir", str(out), "--reference-label", "mouse_kidney_atlas",
            "--source-study-column", "Origin", "--cell-type-column", "Celltype_finest",
            "--fallback-cell-type-column", "Celltype",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    counts, meta = al.load_pseudobulk(
        out / "atlas_pseudobulk_counts.tsv.gz", out / "atlas_pseudobulk_sample_metadata.tsv"
    )
    assert not meta["included"].all()  # S4 podocyte (10 cells) and S2 podocyte (26 -> kept)
    ours = al.compartment_expression_from_pseudobulk(counts, meta, min_cpm=1.0)
    theirs = pd.read_csv(
        out / "atlas_compartment_expression.tsv.gz", sep="\t", float_precision="round_trip"
    )
    assert list(ours.columns) == list(theirs.columns)
    for column in ("gene_symbol", "compartment", "reference_label"):
        assert (ours[column].to_numpy() == theirs[column].to_numpy()).all()
    for column in ("mean_cpm", "median_cpm", "source_study_detection_fraction", "n_source_studies"):
        np.testing.assert_allclose(ours[column], theirs[column], rtol=0, atol=1e-9)
    # and the omitted-study table equals dropping that study before script 03's arithmetic
    omitted = al.compartment_expression_from_pseudobulk(counts, meta, exclude_sources=["S3"])
    assert (omitted.groupby("compartment")["n_source_studies"].first() <= 3).all()


# --- QC rule ------------------------------------------------------------------------------


def test_atlas_qc_flags_compartments_below_minimum():
    _, meta = synthetic_pseudobulk()
    ok = al.atlas_qc(meta, synthetic_qc_cfg(3))
    assert ok["qc_pass"] and ok["studies_by_compartment"]["podocyte"] == 4
    # S5::podocyte is included=False and must not count
    assert ok["studies_by_compartment"]["endothelial"] == 5
    strict = al.atlas_qc(meta, synthetic_qc_cfg(4), exclude_sources=["S1"])
    assert not strict["qc_pass"]
    assert strict["failed_compartments"] == {"podocyte": 3}
    gone = al.atlas_qc(meta, {"minimum_source_studies_per_compartment": 1,
                              "required_compartments": ["DCT2_CNT_atlas"]}, exclude_sources=["S5"])
    assert gone["failed_compartments"] == {"DCT2_CNT_atlas": 0}


# --- tier rebuild ----------------------------------------------------------------------------


def _baseline_atlas():
    counts, meta = synthetic_pseudobulk()
    return counts, meta, al.compartment_expression_from_pseudobulk(counts, meta)


def test_strict_builder_aborts_on_empty_set_but_tolerant_builder_reports_it():
    _, _, atlas = _baseline_atlas()
    cfg = synthetic_audit_cfg()
    with pytest.raises(ValueError, match="disagrees with frozen family"):
        build_marker_tiers(atlas, STRUCTURAL_TERMS, cfg)
    tiers, info = al.build_tiers_tolerant(atlas, STRUCTURAL_TERMS, cfg)
    assert info["status"] == "OK"
    assert info["missing_sets"] == sorted(f"{r}__broad_enriched" for r in MAPPING.values())
    sets = al.set_membership(tiers)
    assert set(sets[al.PODOCYTE_SET]) == set(POD_HS)
    assert set(sets[CONTROL]) == {"Pa0"} | {f"Bg{i}" for i in range(10)}
    assert not (set(sets) & set(info["missing_sets"]))


def test_tolerant_builder_does_not_swallow_other_errors():
    _, _, atlas = _baseline_atlas()
    with pytest.raises(ValueError, match="lacks"):
        al.build_tiers_tolerant(atlas.drop(columns=["mean_cpm"]), STRUCTURAL_TERMS, synthetic_audit_cfg())


def test_absent_compartment_makes_the_rebuild_not_evaluable():
    counts, meta, _ = _baseline_atlas()
    rebuild = al.rebuild_tiers(
        counts, meta, STRUCTURAL_TERMS, synthetic_audit_cfg(), synthetic_qc_cfg(),
        omitted="S5", index=4,
    )
    assert rebuild.status == al.STATUS_NOT_EVALUABLE
    assert rebuild.tiers is None
    assert "DCT2_CNT_atlas" in rebuild.reason
    assert rebuild.podocyte_contributing is False
    row = al.tier_summary_row(rebuild, {}, hs_detection_fraction=0.75)
    assert row["tiers_evaluable"] is False
    assert row["n_podocyte_hs"] == 0


def test_leave_one_study_out_tiers_hand_checked_jaccard_and_contributors():
    counts, meta, atlas = _baseline_atlas()
    cfg, qc = synthetic_audit_cfg(), synthetic_qc_cfg(3)
    base, _ = al.build_tiers_tolerant(atlas, STRUCTURAL_TERMS, cfg)
    full = al.set_membership(base)
    rebuilds = al.leave_one_study_out_tiers(counts, meta, STRUCTURAL_TERMS, cfg, qc)
    assert [r.omitted_origin for r in rebuilds] == STUDIES
    assert [r.index for r in rebuilds] == list(range(5))
    rows = {
        r.omitted_origin: al.tier_summary_row(r, full, hs_detection_fraction=0.75)
        for r in rebuilds
    }
    contributing = {o: rows[o]["podocyte_contributing"] for o in STUDIES}
    assert contributing == {"S1": "yes", "S2": "yes", "S3": "yes", "S4": "yes", "S5": "no"}
    assert rows["S5"]["tiers_evaluable"] is False  # DCT2/CNT only exists in S5
    for omitted in ("S1", "S2", "S3"):
        row = rows[omitted]
        assert row["n_full_podocyte_hs"] == 11
        assert row["n_podocyte_hs"] == 10  # Pedge detection 2/3 < 0.75
        assert row["n_retained_from_full"] == 10
        assert row["jaccard_vs_full"] == pytest.approx(10 / 11)
        assert row["n_source_studies_podocyte"] == 3
        assert row["effective_podocyte_hs_rule"] == "3-of-3"
        assert row["qc_pass"] is True
    s4 = rows["S4"]
    assert s4["n_podocyte_hs"] == 11 and s4["jaccard_vs_full"] == 1.0  # 3/3 detection keeps Pedge
    assert s4["effective_podocyte_hs_rule"] == "3-of-3"


def test_failing_qc_is_flagged_but_the_rebuild_is_still_computed():
    counts, meta, _ = _baseline_atlas()
    rebuild = al.rebuild_tiers(
        counts, meta, STRUCTURAL_TERMS, synthetic_audit_cfg(), synthetic_qc_cfg(4),
        omitted="S1", index=0,
    )
    assert rebuild.qc["qc_pass"] is False
    assert rebuild.tiers is not None and rebuild.status == "OK"
    row = al.tier_summary_row(rebuild, {}, hs_detection_fraction=0.75)
    assert row["qc_pass"] is False and row["qc_failed_compartments"] == "podocyte:3"
    assert row["tiers_evaluable"] is True


def test_podocyte_set_below_eight_genes_is_not_evaluable_not_a_failure():
    counts, meta, atlas = _baseline_atlas()
    cfg = synthetic_audit_cfg()
    tiers, _ = al.build_tiers_tolerant(atlas, STRUCTURAL_TERMS, cfg)
    small = tiers[~((tiers["gene_set"] == al.PODOCYTE_SET) & tiers["gene_symbol"].isin(["Pa0", "Pa1", "Pa2", "Pa3"]))]
    rebuild = al.TierRebuild(
        omitted_origin="S1", index=0, podocyte_contributing=True, qc=al.atlas_qc(meta, synthetic_qc_cfg()),
        status="OK", reason="", tiers=small,
    )
    row = al.tier_summary_row(rebuild, al.set_membership(tiers), hs_detection_fraction=0.75)
    assert row["n_podocyte_hs"] == 7 and row["tiers_evaluable"] is False
    assert "< 8" in row["tier_status_reason"]


def test_set_helpers():
    assert al.jaccard({"a", "b"}, {"b", "c"}) == pytest.approx(1 / 3)
    assert math.isnan(al.jaccard([], []))
    full = {"x": frozenset("abc"), "y": frozenset("de")}
    rebuilt = {"x": frozenset("abd")}
    table = al.compare_sets(full, rebuilt, "S1").set_index("gene_set")
    assert table.loc["x", ["n_full", "n_rebuilt", "n_shared", "n_gained", "n_lost"]].tolist() == [3, 3, 2, 1, 1]
    assert table.loc["y", "n_rebuilt"] == 0 and table.loc["y", "jaccard_vs_full"] == 0.0
    left = pd.DataFrame({"gene_set": ["s", "s"], "gene_symbol": ["a", "b"]})
    right = pd.DataFrame({"gene_set": ["s", "s"], "gene_symbol": ["b", "a"]})
    assert al.tiers_equal_setwise(left, right)["identical"]
    other = pd.DataFrame({"gene_set": ["s", "t"], "gene_symbol": ["a", "b"]})
    cmp = al.tiers_equal_setwise(left, other)
    assert not cmp["identical"] and cmp["only_in_right"] == ["t"] and cmp["differing_sets"] == ["s"]


@pytest.mark.parametrize(
    "n, fraction, expected",
    [(4, 0.75, "3-of-4"), (3, 0.75, "3-of-3"), (8, 0.75, "6-of-8"), (7, 0.75, "6-of-7"),
     (3, 0.5, "2-of-3"), (4, 0.5, "2-of-4"), (0, 0.75, "0-of-0")],
)
def test_effective_detection_rule_table(n, fraction, expected):
    assert al.effective_detection_rule(n, fraction) == expected


@pytest.mark.parametrize(
    "estimate, reference, ratio, direction, passed",
    [
        (0.4, 0.8, 0.5, True, True),
        (0.39, 0.8, 0.4875, True, False),
        (0.8, 0.8, 1.0, True, True),
        (-0.1, 0.8, -0.125, False, False),
        (0.0, 0.5, 0.0, False, False),
        (-0.5, -0.8, 0.625, True, True),
        (0.2, -0.8, -0.25, False, False),
    ],
)
def test_retention_gate_table(estimate, reference, ratio, direction, passed):
    out = al.retention_gate(estimate, reference)
    assert out["retention_ratio"] == pytest.approx(ratio)
    assert out["direction_retained"] is direction
    assert out["gate_pass"] is passed


@pytest.mark.parametrize("estimate, reference", [(0.5, 0.0), (float("nan"), 0.5), (0.5, float("nan")), (None, 0.5)])
def test_retention_gate_is_na_without_a_usable_reference(estimate, reference):
    out = al.retention_gate(estimate, reference)
    assert math.isnan(out["retention_ratio"])
    assert out["gate_pass"] is pd.NA and out["direction_retained"] is pd.NA


def test_retention_threshold_is_a_parameter():
    assert al.retention_gate(0.3, 1.0, retention_min=0.25)["gate_pass"] is True
    assert al.retention_gate(0.3, 1.0, retention_min=0.5)["gate_pass"] is False


def test_seeds_follow_the_plan_offsets():
    assert al.SEED_OFFSET_REFERENCE == 100 and al.SEED_OFFSET_LOSO == 800
    assert [al.loso_seed(20260811, j) for j in (0, 7)] == [20261611, 20261618]


# --- preregistration guard ----------------------------------------------------------------


def _git(repo: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-c", "user.email=t@example.org", "-c", "user.name=t", *args],
        cwd=repo, check=True, capture_output=True,
    )


@pytest.fixture
def tmp_repo(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init", "-q")
    return repo


def test_guard_refuses_missing_untracked_dirty_and_ignored_prereg(tmp_repo):
    prereg = tmp_repo / "config" / "prereg.yaml"
    with pytest.raises(al.PreregistrationError, match="missing"):
        al.check_preregistration(prereg, repo=tmp_repo)
    prereg.parent.mkdir()
    prereg.write_text("analysis_id: t\nretention_ratio_min: 0.5\n")
    with pytest.raises(al.PreregistrationError, match="not committed"):
        al.check_preregistration(prereg, repo=tmp_repo)  # untracked
    _git(tmp_repo, "add", "config/prereg.yaml")
    _git(tmp_repo, "commit", "-q", "-m", "prereg")
    record = al.check_preregistration(prereg, repo=tmp_repo)
    assert record["uncommitted_changes"] is False
    assert record["sha256"] == al.sha256_file(prereg)
    assert record["git_head"] == al.git_head(tmp_repo) and record["git_head"]
    assert record["path"] == str(Path("config") / "prereg.yaml")
    prereg.write_text("analysis_id: t\nretention_ratio_min: 0.4\n")
    with pytest.raises(al.PreregistrationError, match="not committed"):
        al.check_preregistration(prereg, repo=tmp_repo)  # modified after commit
    waived = al.check_preregistration(prereg, repo=tmp_repo, allow_uncommitted=True)
    assert waived["uncommitted_changes"] is True and waived["allow_uncommitted_prereg"] is True
    # missing is never waived
    with pytest.raises(al.PreregistrationError, match="missing"):
        al.check_preregistration(tmp_repo / "nope.yaml", repo=tmp_repo, allow_uncommitted=True)
    # ignored files are untracked: refused
    (tmp_repo / ".gitignore").write_text("ignored.yaml\n")
    ignored = tmp_repo / "ignored.yaml"
    ignored.write_text("analysis_id: i\n")
    with pytest.raises(al.PreregistrationError, match="not tracked"):
        al.check_preregistration(ignored, repo=tmp_repo)


def test_guard_refuses_files_outside_the_repository(tmp_repo, tmp_path):
    outside = tmp_path / "outside.yaml"
    outside.write_text("analysis_id: o\n")
    with pytest.raises(al.PreregistrationError, match="outside the repository"):
        al.check_preregistration(outside, repo=tmp_repo)
    assert al.check_preregistration(outside, repo=tmp_repo, allow_uncommitted=True)["uncommitted_changes"]


# --- end-to-end on synthetic inputs ---------------------------------------------------------


def make_mission(name, rng, *, n_per_group=5, blocks=("all",), shift=1.0, dead_genes=()):
    animals, conditions, block_of = [], [], []
    for block in blocks:
        for condition in ("FLT", "GC"):
            for i in range(n_per_group):
                animals.append(f"{name}_{block}_{condition}{i}")
                conditions.append(condition)
                block_of.append(block)
    n = len(animals)
    expression = pd.DataFrame(rng.normal(8.0, 1.0, (len(SYN_GENES), n)), index=SYN_GENES, columns=animals)
    flight = [a for a, c in zip(animals, conditions) if c == "FLT"]
    expression.loc[POD_HS, flight] += shift
    counts = pd.DataFrame(rng.poisson(100, (len(SYN_GENES), n)), index=SYN_GENES, columns=animals)
    for gene in dead_genes:
        counts.loc[gene, :] = 0
    metadata = pd.DataFrame(
        {"condition": conditions, "block": block_of, "source_sample": animals}, index=animals
    )
    data = MissionData(
        mission=name, expression=expression, counts=counts, metadata=metadata,
        qc=pd.DataFrame(index=animals),
    )
    data.validate()
    return data


def synthetic_missions(seed=11, **kwargs):
    rng = np.random.default_rng(seed)
    return {
        "M1": make_mission("M1", rng, **kwargs),
        "M2": make_mission("M2", rng, **kwargs),
        "M3": make_mission("M3", rng, blocks=("A", "B"), n_per_group=4, **kwargs),
    }


def _write_workspace(tmp_path: Path, *, minimum_studies: int = 3, parent_seed: int = 123):
    counts, meta = synthetic_pseudobulk()
    counts.to_csv(tmp_path / "counts.tsv.gz", sep="\t")
    meta.to_csv(tmp_path / "meta.tsv", sep="\t", index=False)
    kegg = tmp_path / "kegg.json"
    kegg.write_text(json.dumps(STRUCTURAL_TERMS))
    audit = synthetic_audit_cfg()
    audit["input"] = {"structural_scaffold_reference": str(kegg)}
    (tmp_path / "audit.yaml").write_text(yaml.safe_dump(audit))
    (tmp_path / "freeze.yaml").write_text(
        yaml.safe_dump(
            {"reference_builder": {"atlas_qc": synthetic_qc_cfg(minimum_studies),
                                   "broad_expression": {"min_cpm": 1.0}}}
        )
    )
    (tmp_path / "analysis.yaml").write_text(
        yaml.safe_dump(
            {"seed": parent_seed, "eligibility": {"cpm_threshold": 0.1},
             "gene_mapping": {"path": "data/processed/resources/id_map.tsv",
                              "annotation_fallback": "x.csv"}}
        )
    )
    atlas = al.compartment_expression_from_pseudobulk(counts, meta)
    atlas.to_csv(tmp_path / "expression.tsv.gz", sep="\t", index=False, compression="gzip")
    tiers, _ = al.build_tiers_tolerant(atlas, STRUCTURAL_TERMS, audit)
    tiers.to_csv(tmp_path / "frozen_tiers.tsv", sep="\t", index=False)
    (tmp_path / "prereg.yaml").write_text("analysis_id: synthetic_prereg\nretention_ratio_min: 0.5\n")
    return [
        "--audit-config", str(tmp_path / "audit.yaml"),
        "--freeze-config", str(tmp_path / "freeze.yaml"),
        "--analysis-config", str(tmp_path / "analysis.yaml"),
        "--pseudobulk-counts", str(tmp_path / "counts.tsv.gz"),
        "--pseudobulk-meta", str(tmp_path / "meta.tsv"),
        "--expression", str(tmp_path / "expression.tsv.gz"),
        "--frozen-tiers", str(tmp_path / "frozen_tiers.tsv"),
        "--prereg", str(tmp_path / "prereg.yaml"),
        "--results", str(tmp_path / "out"),
    ]


def test_skip_scoring_reads_no_mission_data_and_needs_no_prereg(tmp_path, monkeypatch):
    base = _write_workspace(tmp_path)

    def _forbidden(*args, **kwargs):  # pragma: no cover - must not be reached
        raise AssertionError("mission data must not be loaded with --skip-scoring")

    monkeypatch.setattr(rl, "load_primary_missions", _forbidden)
    monkeypatch.setattr(rl, "evaluate_rebuild", _forbidden)
    (tmp_path / "prereg.yaml").unlink()  # the guard does not apply to --skip-scoring
    assert rl.main(base + ["--skip-scoring"]) == 0
    out = tmp_path / "out"
    summary = pd.read_csv(out / "atlas_loso_summary.tsv", sep="\t")
    assert list(summary.columns) == rl.SUMMARY_COLUMNS
    assert summary["omitted_origin"].tolist() == STUDIES
    assert (summary["status"].isin([al.STATUS_NOT_SCORED, al.STATUS_NOT_EVALUABLE])).all()
    assert summary.set_index("omitted_origin").loc["S5", "status"] == al.STATUS_NOT_EVALUABLE
    assert summary["estimate"].isna().all() and summary["gate_pass"].isna().all()
    assert not summary["gate_included"].any()
    assert summary.set_index("omitted_origin").loc["S1", "n_retained_from_full"] == 10
    for study in STUDIES[:-1]:
        assert (out / f"atlas_loso_tiers_{study}.tsv").is_file()
    assert not (out / "atlas_loso_tiers_S5.tsv").exists()
    assert not (out / "atlas_loso_family_meta.tsv.gz").exists()
    manifest = json.loads((out / "atlas_loso_manifest.json").read_text())
    assert manifest["skip_scoring"] is True and manifest["preregistration"]["applicable"] is False
    assert manifest["baseline_contract"]["tiers_identical_to_frozen"] is True
    assert manifest["baseline_contract"]["expression_max_abs_diff"] == 0.0
    assert manifest["thresholds"]["status"] == "unchanged_frozen"
    assert "3-of-4" in manifest["effective_stringency_disclosure"]
    sets = pd.read_csv(out / "atlas_loso_sets.tsv.gz", sep="\t")
    assert set(sets["omitted_origin"]) == set(STUDIES[:-1])
    assert {"omitted_origin", "gene_set", "gene_symbol", "in_full_set"} <= set(sets.columns)
    comparison = pd.read_csv(out / "atlas_loso_set_comparison.tsv", sep="\t")
    row = comparison[(comparison["omitted_origin"] == "S1") & (comparison["gene_set"] == al.PODOCYTE_SET)].iloc[0]
    assert (row["n_full"], row["n_rebuilt"], row["n_lost"], row["n_gained"]) == (11, 10, 1, 0)


def test_baseline_mismatch_aborts_unless_allowed(tmp_path):
    base = _write_workspace(tmp_path)
    tiers = pd.read_csv(tmp_path / "frozen_tiers.tsv", sep="\t")
    tiers[tiers["gene_symbol"] != "Pa3"].to_csv(tmp_path / "frozen_tiers.tsv", sep="\t", index=False)
    with pytest.raises(RuntimeError, match="does not reproduce the frozen tiers"):
        rl.main(base + ["--skip-scoring"])
    assert rl.main(base + ["--skip-scoring", "--allow-baseline-mismatch"]) == 0


def test_scoring_refuses_without_committed_prereg(tmp_path, monkeypatch, capsys):
    base = _write_workspace(tmp_path)
    monkeypatch.setattr(rl, "load_primary_missions", lambda *a, **k: pytest.fail("loaded flight data"))
    (tmp_path / "prereg.yaml").unlink()
    assert rl.main(base + ["--permutations", "50"]) == 2
    assert "REFUSED" in capsys.readouterr().err
    assert not (tmp_path / "out" / "atlas_loso_summary.tsv").exists()
    (tmp_path / "prereg.yaml").write_text("analysis_id: x\n")  # exists but outside any commit
    assert rl.main(base + ["--permutations", "50"]) == 2
    assert not (tmp_path / "out" / "atlas_loso_summary.tsv").exists()


def test_default_prereg_path_is_the_g1_preregistration():
    assert rl.DEFAULT_PREREG.name == "clinical_axes_podocyte_sensitivities.yaml"


def test_scoring_end_to_end_on_synthetic_missions(tmp_path):
    base = _write_workspace(tmp_path, parent_seed=123)
    missions = synthetic_missions()
    args = rl.build_parser().parse_args(
        base + ["--allow-uncommitted-prereg", "--permutations", "200", "--chunk-size", "64"]
    )
    out = rl.run(args, missions=missions)
    summary = out["summary"].set_index("omitted_origin")
    manifest = json.loads((tmp_path / "out" / "atlas_loso_manifest.json").read_text())

    # seeds are parent + 800 + j (alphabetical), reference uses parent + 100
    assert summary["seed"].tolist() == [123 + 800 + j for j in range(5)]
    assert manifest["reference_full_atlas"]["seed"] == 223
    assert manifest["preregistration"]["uncommitted_changes"] is True
    assert manifest["preregistration"]["analysis_id"] == "synthetic_prereg"
    assert manifest["preregistration"]["sha256"] == al.sha256_file(tmp_path / "prereg.yaml")

    # S5 removes DCT2/CNT altogether -> NOT_EVALUABLE, excluded from the gate, not a failure
    assert summary.loc["S5", "status"] == al.STATUS_NOT_EVALUABLE
    assert not summary.loc["S5", "gate_included"] and pd.isna(summary.loc["S5", "gate_pass"])
    evaluable = summary.drop(index="S5")
    assert (evaluable["status"] == al.STATUS_EVALUABLE).all()
    assert evaluable["gate_included"].all() and evaluable["qc_pass"].all()
    assert (evaluable["n_observable_all_missions"] == evaluable["n_podocyte_hs"]).all()
    assert evaluable["n_family_sets"].tolist() == [13] * 4
    reference = manifest["reference_full_atlas"]
    assert reference["estimate"] > 0 and reference["n_family_sets"] == 13
    assert reference["n_observable_all_missions"] == 11
    assert (evaluable["reference_estimate"] == reference["estimate"]).all()
    np.testing.assert_allclose(
        evaluable["retention_ratio"], evaluable["estimate"] / reference["estimate"]
    )
    assert evaluable["direction_retained"].all() and evaluable["gate_pass"].all()
    assert (evaluable["rank_by_abs_t"] >= 1).all()
    assert ((evaluable["max_t_fwer"] > 0) & (evaluable["max_t_fwer"] <= 1)).all()
    assert (evaluable["k0_status"] == "ok").all()
    assert "k0_p1_coef_score_units" in summary.columns and not any(
        c.startswith("k0_p1_g") for c in summary.columns
    )
    assert evaluable["k0_p1_positive"].all()
    gate = manifest["gate"]
    assert gate["overall"]["n_included"] == 4 and gate["overall"]["all_pass"] is True
    assert gate["by_podocyte_contributing"]["podocyte_contributing_yes"]["n_included"] == 4
    assert gate["by_podocyte_contributing"]["podocyte_contributing_no"]["n_included"] == 0
    assert gate["n_by_status"] == {al.STATUS_EVALUABLE: 4, al.STATUS_NOT_EVALUABLE: 1}
    meta = pd.read_csv(tmp_path / "out" / "atlas_loso_family_meta.tsv.gz", sep="\t")
    assert set(meta["omitted_origin"]) == {rl.REFERENCE_LABEL} | set(STUDIES[:-1])
    assert meta[meta["omitted_origin"] == "S1"]["axis"].is_unique
    pod = meta[(meta["axis"] == al.PODOCYTE_SET) & (meta["omitted_origin"] == "S1")].iloc[0]
    assert pod["estimate"] == pytest.approx(summary.loc["S1", "estimate"])

    # deterministic: same seeds, same numbers
    again = rl.run(
        rl.build_parser().parse_args(
            base[:-1] + [str(tmp_path / "out2"), "--allow-uncommitted-prereg",
                         "--permutations", "200", "--chunk-size", "64"]
        ),
        missions=missions,
    )["summary"].set_index("omitted_origin")
    pd.testing.assert_series_equal(
        summary["estimate"], again["estimate"], check_names=False
    )
    pd.testing.assert_series_equal(
        summary["max_t_fwer"], again["max_t_fwer"], check_names=False
    )


def test_reference_crosscheck_against_a_compartment_context_meta_file(tmp_path):
    base = _write_workspace(tmp_path)
    missions = synthetic_missions()
    args = rl.build_parser().parse_args(
        base + ["--allow-uncommitted-prereg", "--permutations", "100", "--chunk-size", "50"]
    )
    rl.run(args, missions=missions)
    meta = pd.read_csv(tmp_path / "out" / "atlas_loso_family_meta.tsv.gz", sep="\t")
    ref = meta[meta["omitted_origin"] == rl.REFERENCE_LABEL].drop(columns="omitted_origin")
    path = tmp_path / "compartment_context_meta_results.tsv"
    ref.to_csv(path, sep="\t", index=False)
    args = rl.build_parser().parse_args(
        base[:-1] + [str(tmp_path / "out_x"), "--allow-uncommitted-prereg", "--permutations", "100",
                     "--chunk-size", "50", "--reference-meta", str(path)]
    )
    rl.run(args, missions=missions)
    check = json.loads((tmp_path / "out_x" / "atlas_loso_manifest.json").read_text())[
        "reference_full_atlas"]["crosscheck"]
    assert check["estimate_matches"] and check["abs_diff_estimate"] < 1e-9
    assert check["abs_diff_max_t_fwer"] < 1e-12
    shifted = ref.copy()
    shifted.loc[shifted["axis"] == al.PODOCYTE_SET, "estimate"] += 0.1
    shifted.to_csv(path, sep="\t", index=False)
    rl.run(args, missions=missions)
    mismatch = json.loads((tmp_path / "out_x" / "atlas_loso_manifest.json").read_text())[
        "reference_full_atlas"]["crosscheck"]
    assert not mismatch["estimate_matches"]


def test_failing_qc_rebuilds_are_computed_but_excluded_from_the_gate(tmp_path):
    base = _write_workspace(tmp_path, minimum_studies=4)
    args = rl.build_parser().parse_args(
        base + ["--allow-uncommitted-prereg", "--permutations", "100", "--chunk-size", "50"]
    )
    summary = rl.run(args, missions=synthetic_missions())["summary"].set_index("omitted_origin")
    for study in ("S1", "S2", "S3", "S4"):
        row = summary.loc[study]
        assert not row["qc_pass"] and row["qc_failed_compartments"] == "podocyte:3"
        assert row["status"] == al.STATUS_EVALUABLE and not pd.isna(row["estimate"])
        assert not row["gate_included"] and pd.isna(row["gate_pass"])
    manifest = json.loads((tmp_path / "out" / "atlas_loso_manifest.json").read_text())
    assert manifest["gate"]["overall"] == {"n_included": 0, "n_pass": 0, "all_pass": None}


def test_rebuild_with_too_few_observable_podocyte_genes_is_not_evaluable(tmp_path):
    """3 dead genes: the 11-gene reference (8 observable) and S4 (11 genes) score, S1-S3 (10 genes) do not."""

    base = _write_workspace(tmp_path)
    missions = synthetic_missions(dead_genes=["Pa0", "Pa1", "Pa2"])
    args = rl.build_parser().parse_args(
        base + ["--allow-uncommitted-prereg", "--permutations", "50"]
    )
    summary = rl.run(args, missions=missions)["summary"].set_index("omitted_origin")
    for study in ("S1", "S2", "S3"):
        row = summary.loc[study]
        assert row["status"] == al.STATUS_NOT_EVALUABLE and "observable" in row["status_reason"]
        assert row["n_podocyte_hs"] == 10 and not row["gate_included"] and pd.isna(row["gate_pass"])
        assert pd.isna(row["estimate"]) and pd.isna(row["k0_p1_coef_score_units"])
    assert summary.loc["S4", "status"] == al.STATUS_EVALUABLE
    assert summary.loc["S4", "n_observable_all_missions"] == 8
    manifest = json.loads((tmp_path / "out" / "atlas_loso_manifest.json").read_text())
    assert manifest["gate"]["overall"]["n_included"] == 1
    assert manifest["gate"]["not_evaluable_is_not_a_failure"] is True


def test_scoring_aborts_when_the_full_atlas_reference_is_not_evaluable(tmp_path):
    base = _write_workspace(tmp_path)
    missions = synthetic_missions(dead_genes=["Pa0", "Pa1", "Pa2", "Pa3"])
    args = rl.build_parser().parse_args(
        base + ["--allow-uncommitted-prereg", "--permutations", "50"]
    )
    with pytest.raises(RuntimeError, match="reference is not evaluable"):
        rl.run(args, missions=missions)


def test_evaluate_rebuild_reports_unobservable_podocyte_rebuild(tmp_path):
    base = _write_workspace(tmp_path)
    tiers = pd.read_csv(tmp_path / "frozen_tiers.tsv", sep="\t")
    missions = synthetic_missions(dead_genes=["Pa0", "Pa1", "Pa2", "Pa3"])
    out = rl.evaluate_rebuild(
        tmp_path / "frozen_tiers.tsv", missions, threshold=0.1, n_permutations=20, seed=1, chunk_size=20
    )
    assert out["status"] == al.STATUS_NOT_EVALUABLE and "observable" in out["reason"]
    # < 8 defined genes: not even loaded as a family member
    small = tiers[~((tiers["gene_set"] == al.PODOCYTE_SET) & tiers["gene_symbol"].str.startswith("Pa"))]
    small.to_csv(tmp_path / "small.tsv", sep="\t", index=False)
    out = rl.evaluate_rebuild(
        tmp_path / "small.tsv", synthetic_missions(), threshold=0.1, n_permutations=20, seed=1, chunk_size=20
    )
    assert out["status"] == al.STATUS_NOT_EVALUABLE and "defined with 1 genes" in out["reason"]


def test_gene_map_override_is_in_memory_only(tmp_path, monkeypatch):
    base = _write_workspace(tmp_path)
    gene_map = tmp_path / "alt_map.tsv"
    gene_map.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG1\tPa0\n")
    captured = {}

    def fake_loader(config, root):
        captured["path"] = config["gene_mapping"]["path"]
        return synthetic_missions()

    monkeypatch.setattr(rl, "load_primary_missions", fake_loader)
    before = (tmp_path / "analysis.yaml").read_bytes()
    args = rl.build_parser().parse_args(
        base + ["--allow-uncommitted-prereg", "--permutations", "50", "--gene-map", str(gene_map)]
    )
    rl.run(args)
    assert captured["path"] == str(gene_map.resolve())
    assert (tmp_path / "analysis.yaml").read_bytes() == before  # config file untouched
    manifest = json.loads((tmp_path / "out" / "atlas_loso_manifest.json").read_text())
    assert manifest["gene_map"]["override"] is True
    assert manifest["gene_map"]["sha256"] == al.sha256_file(gene_map)
    # without the flag the config's own path is used
    captured.clear()
    args = rl.build_parser().parse_args(
        base[:-1] + [str(tmp_path / "out_default"), "--allow-uncommitted-prereg", "--permutations", "50"]
    )
    rl.run(args)
    assert captured["path"] == "data/processed/resources/id_map.tsv"


def test_gene_map_override_requires_an_existing_file_and_copies_the_config(tmp_path):
    config = {"gene_mapping": {"path": "data/processed/resources/id_map.tsv"}}
    assert rl.apply_gene_map_override(config, None) is config
    assert config["gene_mapping"]["path"] == "data/processed/resources/id_map.tsv"
    with pytest.raises(FileNotFoundError, match="does not exist"):
        rl.apply_gene_map_override(config, tmp_path / "nope.tsv")
    alt = tmp_path / "m.tsv"
    alt.write_text("ensembl_gene_id\tmgi_symbol\n")
    assert rl.apply_gene_map_override(config, alt)["gene_mapping"]["path"] == str(alt.resolve())


def test_k0_p1_matches_independent_recompute():
    missions = synthetic_missions(seed=5)
    eligible = {m: set(cpm_eligible_genes(d, 0.1)) for m, d in missions.items()}
    structural = ["Pa0"] + [f"Bg{i}" for i in range(10)]
    got = rl.k0_p1(missions, POD_HS, structural, eligible)
    assert got["k0_status"] == "ok"
    disjoint = [g for g in structural if g not in set(POD_HS)]
    assert disjoint == [f"Bg{i}" for i in range(10)]  # overlap removed from the comparator

    def score(expression, genes):
        x = expression.loc[genes]
        z = x.sub(x.mean(axis=1), axis=0).div(x.std(axis=1, ddof=1), axis=0)
        return z.mean(axis=0)

    estimates, variances = [], []
    for data in missions.values():
        podocyte = score(data.expression, POD_HS)
        comparator = score(data.expression, disjoint)
        design = pd.DataFrame({"flight": (data.metadata["condition"] == "FLT").astype(float)})
        design = design.join(pd.get_dummies(data.metadata["block"], drop_first=True, dtype=float))
        design["adjuster"] = comparator
        fit = sm.OLS(podocyte, sm.add_constant(design, has_constant="add")).fit(cov_type="HC3")
        estimates.append(float(fit.params["flight"]))
        variances.append(float(fit.bse["flight"]) ** 2)
    pooled = random_effects_reml_mkh(np.array(estimates), np.array(variances))
    assert got["k0_p1_coef_score_units"] == pytest.approx(pooled.estimate, abs=1e-9)
    assert got["k0_p1_ci_low"] == pytest.approx(pooled.ci_low, abs=1e-9)
    assert got["k0_p1_ci_high"] == pytest.approx(pooled.ci_high, abs=1e-9)
    assert got["k0_p1_positive"] is True


def test_k0_p1_not_evaluable_with_too_few_comparator_genes():
    missions = synthetic_missions(seed=5)
    eligible = {m: set(cpm_eligible_genes(d, 0.1)) for m, d in missions.items()}
    out = rl.k0_p1(missions, POD_HS, ["Pa0"] + [f"Bg{i}" for i in range(7)], eligible)
    assert out["k0_status"].startswith("NOT_EVALUABLE") and "k0_p1_coef_score_units" not in out


def test_family_scoring_matches_direct_blocked_permutation():
    """evaluate_rebuild is the compartment-context code path on the rebuilt family."""

    from scripts.clinical_axes.run_compartment_context import load_family
    from src.clinical_axes.analysis import combined_score_design
    from src.clinical_axes.statistics import blocked_meta_permutation
    import tempfile

    counts, meta, atlas = _baseline_atlas()
    tiers, _ = al.build_tiers_tolerant(atlas, STRUCTURAL_TERMS, synthetic_audit_cfg())
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "t.tsv"
        tiers.to_csv(path, sep="\t", index=False)
        missions = synthetic_missions(seed=3)
        got = rl.evaluate_rebuild(path, missions, threshold=0.1, n_permutations=120, seed=77, chunk_size=40)
        family, _ = load_family(path)
    scores, design, _, _ = combined_score_design(missions, family, cpm_threshold=0.1)
    direct = blocked_meta_permutation(scores, design, n_permutations=120, seed=77, chunk_size=40)
    row = got["podocyte_row"]
    assert got["n_family_sets"] == len(family)
    assert row["estimate"] == pytest.approx(direct.observed_meta.loc[al.PODOCYTE_SET, "estimate"], abs=1e-12)
    assert row["max_t_fwer"] == pytest.approx(direct.max_t_fwer[al.PODOCYTE_SET], abs=1e-12)
    ranks = direct.observed_meta["t_mkh"].abs().rank(ascending=False, method="min")
    assert row["rank_by_abs_t"] == int(ranks[al.PODOCYTE_SET])


# --- real-data contracts (flight-blind, already-published atlas artefacts) -----------------


@real_atlas
def test_real_excluding_nothing_reproduces_the_stored_compartment_expression():
    counts, meta = al.load_pseudobulk(COUNTS, SAMPLE_META)
    ours = al.compartment_expression_from_pseudobulk(counts, meta, min_cpm=1.0)
    stored = pd.read_csv(EXPRESSION, sep="\t", float_precision="round_trip")
    assert list(ours.columns) == list(stored.columns)
    assert len(ours) == len(stored)
    for column in ("gene_symbol", "compartment", "reference_label"):
        assert (ours[column].to_numpy() == stored[column].to_numpy()).all()
    for column in ("mean_cpm", "median_cpm", "source_study_detection_fraction", "n_source_studies"):
        assert np.abs(ours[column].to_numpy() - stored[column].to_numpy()).max() <= 1e-9


@real_atlas
def test_real_build_marker_tiers_on_rebuilt_expression_reproduces_frozen_tiers_setwise():
    counts, meta = al.load_pseudobulk(COUNTS, SAMPLE_META)
    audit = load_config(AUDIT_CONFIG)
    atlas = al.compartment_expression_from_pseudobulk(counts, meta, min_cpm=1.0)
    tiers = build_marker_tiers(atlas, json.loads(KEGG.read_text()), audit)
    frozen = pd.read_csv(FROZEN_TIERS, sep="\t")
    comparison = al.tiers_equal_setwise(tiers, frozen)
    assert comparison["identical"], comparison
    assert set(al.set_membership(tiers)) == set(audit["set_test"]["primary_family"])
    assert len(al.set_membership(frozen)[al.PODOCYTE_SET]) == len(
        al.set_membership(tiers)[al.PODOCYTE_SET]
    )
    # the tolerant builder is the strict builder when no set is empty
    tolerant, info = al.build_tiers_tolerant(atlas, json.loads(KEGG.read_text()), audit)
    assert info["missing_sets"] == [] and al.tiers_equal_setwise(tolerant, frozen)["identical"]


@real_atlas
def test_real_study_counts_match_the_pseudobulk_qc_record():
    meta = pd.read_csv(SAMPLE_META, sep="\t")
    recorded = json.loads(QC_JSON.read_text())["source_studies_by_compartment"]
    assert al.studies_by_compartment(meta) == recorded
    freeze = yaml.safe_load(FREEZE_CONFIG.read_text())["reference_builder"]["atlas_qc"]
    assert al.atlas_qc(meta, freeze)["qc_pass"]


@real_atlas
def test_real_leave_one_study_out_rebuild_structure():
    counts, meta = al.load_pseudobulk(COUNTS, SAMPLE_META)
    audit = load_config(AUDIT_CONFIG)
    freeze = yaml.safe_load(FREEZE_CONFIG.read_text())["reference_builder"]
    full_n = al.studies_by_compartment(meta)
    rebuilds = al.leave_one_study_out_tiers(
        counts, meta, json.loads(KEGG.read_text()), audit, freeze["atlas_qc"],
        min_cpm=float(freeze["broad_expression"]["min_cpm"]),
    )
    origins = sorted(meta["source_study_id"].unique())
    assert [r.omitted_origin for r in rebuilds] == origins
    included = meta[meta["included"]]
    for rebuild in rebuilds:
        contributed = set(included.loc[included["source_study_id"] == rebuild.omitted_origin, "compartment"])
        expected = {c: n - (c in contributed) for c, n in full_n.items()}
        assert rebuild.qc["studies_by_compartment"] == {c: n for c, n in expected.items() if n > 0}
        assert rebuild.podocyte_contributing == ("podocyte" in contributed)
        if rebuild.tiers is not None:
            assert rebuild.status == "OK"
            assert set(al.set_membership(rebuild.tiers)) <= set(audit["set_test"]["primary_family"])


@real_atlas
def test_real_skip_scoring_cli_writes_a_blind_tier_table(tmp_path, monkeypatch):
    def _forbidden(*args, **kwargs):  # pragma: no cover
        raise AssertionError("mission data must not be loaded with --skip-scoring")

    monkeypatch.setattr(rl, "load_primary_missions", _forbidden)
    assert rl.main(["--skip-scoring", "--results", str(tmp_path)]) == 0
    summary = pd.read_csv(tmp_path / "atlas_loso_summary.tsv", sep="\t")
    meta = pd.read_csv(SAMPLE_META, sep="\t")
    assert len(summary) == meta["source_study_id"].nunique()
    expected_yes = int(
        meta[meta["included"] & (meta["compartment"] == "podocyte")]["source_study_id"].nunique()
    )
    assert (summary["podocyte_contributing"] == "yes").sum() == expected_yes
    manifest = json.loads((tmp_path / "atlas_loso_manifest.json").read_text())
    assert manifest["baseline_contract"]["tiers_identical_to_frozen"] is True
    assert manifest["baseline_contract"]["expression_max_abs_diff"] <= 1e-9
    assert summary["gate_pass"].isna().all() and summary["estimate"].isna().all()


def test_json_default_serializes_preregistration_dates():
    prereg = yaml.safe_load((REPO / "config/clinical_axes_podocyte_sensitivities.yaml").read_text())
    encoded = json.loads(json.dumps(prereg, default=rl._json_default))
    assert encoded["lock_date"] == "2026-10-06"
