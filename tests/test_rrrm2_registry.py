"""Registry invariants: every revision parses, resolves, and points at real entries."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from conftest import missing_paths, requires_paths

from rrrm2 import panels as panels_mod
from rrrm2 import runner
from rrrm2.paths import REPO_ROOT
from rrrm2.registry import (
    RegistryError,
    VALID_RUNNER,
    VALID_STATUS,
    load_registry,
    load_revision,
)

REGISTRY = load_registry()


def test_registry_is_not_empty():
    assert REGISTRY, "config/revisions/ contains no revision files"


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_revision_metadata_is_complete(rev_id):
    rev = REGISTRY[rev_id]
    assert rev.status in VALID_STATUS
    assert rev.title.strip()
    assert rev.question.strip(), f"{rev_id} must state its scientific question"
    assert rev.stages, f"{rev_id} has no stages"


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_stage_entries_exist_on_disk(rev_id):
    for stage in REGISTRY[rev_id].stages:
        assert stage.entry_path.exists(), f"{rev_id}/{stage.id}: missing {stage.entry}"
        assert stage.runner in VALID_RUNNER


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_spec_config_and_docs_exist(rev_id):
    rev = REGISTRY[rev_id]
    if rev.spec_config:
        assert (REPO_ROOT / rev.spec_config).exists(), rev.spec_config
    for doc in rev.docs:
        assert (REPO_ROOT / doc).exists(), f"{rev_id}: missing doc {doc}"


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_stage_order_is_resolvable_and_respects_dependencies(rev_id):
    rev = REGISTRY[rev_id]
    ordered = rev.ordered_stages(with_optional=True)
    assert len(ordered) == len(rev.stages)
    seen: set[str] = set()
    for stage in ordered:
        for need in stage.needs:
            assert need in seen, f"{rev_id}/{stage.id} runs before its dependency {need}"
        seen.add(stage.id)


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_stage_placeholders_all_resolve(rev_id):
    rev = REGISTRY[rev_id]
    subs = runner.substitutions(rev, Path("/tmp/run"), "RUNID")
    for stage in rev.stages:
        stage.resolved_args(subs)  # raises RegistryError on an unknown placeholder


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_results_root_is_namespaced_by_revision(rev_id):
    rev = REGISTRY[rev_id]
    assert rev.results_root == f"data/results/{rev_id}", (
        "each revision must write into its own results folder"
    )


def test_upstream_of_includes_transitive_dependencies():
    rev = REGISTRY["v13-phospho"]
    assert [s.id for s in rev.upstream_of("provenance")] == ["infer", "report", "provenance"]


def test_id_must_match_filename(tmp_path):
    bad = tmp_path / "alpha.yaml"
    bad.write_text(yaml.safe_dump({
        "id": "beta", "status": "complete", "results_root": "data/results/beta",
        "stages": [],
    }))
    with pytest.raises(RegistryError, match="does not match the filename"):
        load_revision(bad)


def test_unknown_stage_dependency_is_rejected(tmp_path):
    bad = tmp_path / "alpha.yaml"
    bad.write_text(yaml.safe_dump({
        "id": "alpha", "status": "complete", "results_root": "data/results/alpha",
        "stages": [{"id": "a", "title": "A", "runner": "python",
                    "entry": "scripts/audit_fix_status.py", "needs": ["ghost"]}],
    }))
    with pytest.raises(RegistryError, match="unknown stage"):
        load_revision(bad)


def test_bad_status_is_rejected(tmp_path):
    bad = tmp_path / "alpha.yaml"
    bad.write_text(yaml.safe_dump({
        "id": "alpha", "status": "probably-fine",
        "results_root": "data/results/alpha", "stages": [],
    }))
    with pytest.raises(RegistryError, match="status"):
        load_revision(bad)


def test_panel_registry_parses_and_ids_are_unique():
    registry = panels_mod.load_panels()
    for panel in registry.values():
        assert panel.genes, f"panel {panel.id} has no genes"
        assert panel.species in {"mouse", "human"}
        assert len(panel.digest()) == 64


# --- stage-level ``requires`` ----------------------------------------------------

REAL_ENTRY = "scripts/audit_fix_status.py"  # any script that exists in the repo
PRESENT_PATH = "config/revisions/_schema.yaml"
ABSENT_PATH = "data/definitely/not/here/input.tsv"


def _write_revision(tmp_path: Path, stages: list[dict], **top) -> Path:
    path = tmp_path / "alpha.yaml"
    path.write_text(yaml.safe_dump({
        "id": "alpha", "status": "complete", "results_root": "data/results/alpha",
        "stages": stages, **top,
    }))
    return path


def _stage(stage_id: str, **extra) -> dict:
    return {"id": stage_id, "title": stage_id.upper(), "runner": "python",
            "entry": REAL_ENTRY, **extra}


def test_stage_requires_parses_to_a_tuple(tmp_path):
    rev = load_revision(_write_revision(tmp_path, [
        _stage("a", requires=[PRESENT_PATH, ABSENT_PATH]),
        _stage("b"),
    ]))
    assert rev.stage("a").requires == (PRESENT_PATH, ABSENT_PATH)
    assert rev.stage("b").requires == ()
    assert rev.stage("a").missing_requires() == [ABSENT_PATH]


@pytest.mark.parametrize("bad", ["{run_dir}/x.tsv", "/etc/passwd"])
def test_stage_requires_must_be_literal_repo_relative_paths(tmp_path, bad):
    path = _write_revision(tmp_path, [_stage("a", requires=[bad])])
    with pytest.raises(RegistryError, match="requires"):
        load_revision(path)


def test_stage_requires_must_be_a_list(tmp_path):
    path = _write_revision(tmp_path, [_stage("a", requires=PRESENT_PATH)])
    with pytest.raises(RegistryError, match="list of repo-relative paths"):
        load_revision(path)


def test_preflight_flags_a_missing_stage_requirement(tmp_path):
    rev = load_revision(_write_revision(tmp_path, [
        _stage("has-input", requires=[PRESENT_PATH]),
        _stage("lacks-input", requires=[PRESENT_PATH, ABSENT_PATH]),
    ]))
    check = runner.preflight(rev)
    assert not check.ok
    assert check.problems == [f"stage lacks-input: missing required input {ABSENT_PATH}"]
    assert list(check.blocked_stages) == ["lacks-input"]
    assert check.degraded_stages == {}


def test_preflight_only_inspects_the_selected_stages(tmp_path):
    rev = load_revision(_write_revision(tmp_path, [
        _stage("ok-stage"), _stage("lacks-input", requires=[ABSENT_PATH]),
    ]))
    assert runner.preflight(rev, stages=[rev.stage("ok-stage")]).ok
    assert not runner.preflight(rev, stages=[rev.stage("lacks-input")]).ok


def test_preflight_optional_stage_is_a_warning_for_doctor_only(tmp_path):
    rev = load_revision(_write_revision(tmp_path, [
        _stage("required-ok"),
        _stage("extra", optional=True, requires=[ABSENT_PATH]),
    ]))
    strict = runner.preflight(rev)
    assert not strict.ok and strict.blocked_stages.keys() == {"extra"}

    lenient = runner.preflight(rev, optional_as_warnings=True)
    assert lenient.ok and not lenient.problems and not lenient.blocked_stages
    assert lenient.degraded_stages.keys() == {"extra"}
    assert any(f"stage extra: missing required input {ABSENT_PATH}" in w
               for w in lenient.warnings)


def test_preflight_still_reports_revision_level_data_and_entry(tmp_path):
    rev = load_revision(_write_revision(
        tmp_path,
        [{**_stage("a"), "entry": "scripts/no_such_entry.py"}],
        requires={"data": [ABSENT_PATH]},
    ))
    check = runner.preflight(rev)
    assert not check.ok
    assert any("missing declared data input" in p for p in check.problems)
    assert "stage a: entry not found: scripts/no_such_entry.py" in check.problems
    assert "a" in check.blocked_stages


# --- no two stages may write the same files --------------------------------------

def _output_collisions(rev) -> dict[tuple[str, str], list[str]]:
    """Group stages by (entry script, resolved output directory); return duplicates.

    Stages that pass none of --results/--outdir/--output-dir write wherever the
    script's own default points. That is not a statement about collisions, so
    they are skipped rather than compared (several legacy stages call one entry
    with different positional modes).
    """
    subs = runner.substitutions(rev, Path("/tmp/run"), "RUNID")
    groups: dict[tuple[str, str], list[str]] = {}
    for stage in rev.stages:
        out = stage.output_dir(subs)
        if out is None:
            continue
        groups.setdefault((stage.entry, str(Path(out))), []).append(stage.id)
    return {key: ids for key, ids in groups.items() if len(ids) > 1}


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_no_two_stages_share_entry_script_and_output_dir(rev_id):
    collisions = _output_collisions(REGISTRY[rev_id])
    assert not collisions, (
        f"{rev_id}: stages run the same script into the same directory and will "
        f"overwrite each other: {collisions}"
    )


def test_output_collision_check_catches_the_overwrite_bug(tmp_path):
    """The pre-fix clinical-axes registry: podocyte-disjoint reused {run_dir}."""
    same_dir = ["--config", "{spec_config}", "--results", "{run_dir}"]
    rev = load_revision(_write_revision(tmp_path, [
        _stage("compartment-context", args=same_dir),
        _stage("podocyte-disjoint", args=[*same_dir, "--add-podocyte-disjoint-variants"]),
    ]))
    assert _output_collisions(rev) == {
        (REAL_ENTRY, "/tmp/run"): ["compartment-context", "podocyte-disjoint"]
    }


def test_output_collision_check_accepts_equals_form_and_distinct_subfolders(tmp_path):
    rev = load_revision(_write_revision(tmp_path, [
        _stage("a", args=["--outdir={run_dir}"]),
        _stage("b", args=["--outdir", "{run_dir}/b"]),
        _stage("c", args=["--output-dir={run_dir}/c"]),
        _stage("no-flag-1", args=["--mode", "x"]),
        _stage("no-flag-2", args=["--mode", "y"]),
    ]))
    assert _output_collisions(rev) == {}
    rev_dup = load_revision(_write_revision(tmp_path, [
        _stage("a", args=["--outdir={run_dir}"]),
        _stage("b", args=["--output-dir", "{run_dir}"]),
    ]))
    assert list(_output_collisions(rev_dup).values()) == [["a", "b"]]


# --- clinical-axes registry hygiene (WP-H) ---------------------------------------

TIERS = "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
ATLAS_EXPR = "data/processed/v13_dct_reference/mouse_kidney_atlas/atlas_compartment_expression.tsv.gz"
OSD771_META = "data/raw/metadata/s_OSD-771.txt"


def _flag_value(stage, flag: str) -> str:
    args = list(stage.args)
    return args[args.index(flag) + 1]


def test_podocyte_disjoint_writes_to_its_own_subfolder():
    ca = REGISTRY["clinical-axes"]
    assert _flag_value(ca.stage("compartment-context"), "--results") == "{run_dir}"
    assert _flag_value(ca.stage("podocyte-disjoint"), "--results") == (
        "{run_dir}/disjoint_podocyte_family"
    )
    assert "--add-podocyte-disjoint-variants" in ca.stage("podocyte-disjoint").args


def test_c3h_sensitivity_stage_is_registered_with_its_own_folder():
    ca = REGISTRY["clinical-axes"]
    st = ca.stage("podocyte-disjoint-c3h")
    assert st.entry == "scripts/clinical_axes/run_compartment_context.py"
    assert st.needs == () and not st.optional
    assert _flag_value(st, "--results") == "{run_dir}/disjoint_podocyte_family_c3h"
    assert _flag_value(st, "--osd253-strain") == "C3H/HeJ"
    assert "--add-podocyte-disjoint-variants" in st.args
    assert _flag_value(st, "--permutations") == "100000"


def test_podxl_nid1_stage_is_registered():
    st = REGISTRY["clinical-axes"].stage("podxl-nid1")
    assert st.entry == "scripts/clinical_axes/run_podocyte_podxl_nid1_sensitivity.py"
    assert st.needs == () and not st.optional
    assert _flag_value(st, "--config") == "{spec_config}"
    assert _flag_value(st, "--results") == "{run_dir}/podxl_nid1_forced_sensitivity"
    assert _flag_value(st, "--permutations") == "100000"


def test_figures_need_every_input_they_read():
    figures = REGISTRY["clinical-axes"].stage("figures")
    assert set(figures.needs) == {
        "cross-mission", "compartment-context", "sensitivities",
        "matched-null", "barrier-adjustment",
    }


def test_k0_has_no_stage_dependency_and_pins_seed_zero():
    k0 = REGISTRY["clinical-axes"].stage("podocyte-scaffold-specificity")
    assert k0.needs == ()
    assert _flag_value(k0, "--seed") == "0"


@pytest.mark.parametrize("stage_id", ["atlas-pseudobulk", "marker-tiers"])
def test_reference_build_stages_are_optional_and_say_they_write_to_data_processed(stage_id):
    ca = REGISTRY["clinical-axes"]
    st = ca.stage(stage_id)
    assert st.optional
    assert "data/processed" in st.title
    assert stage_id not in {n for s in ca.stages for n in s.needs}, (
        "no stage may silently depend on a stage that rewrites data/processed"
    )


def test_reference_build_stages_are_the_atlas_and_tier_commands():
    ca = REGISTRY["clinical-axes"]
    atlas = ca.stage("atlas-pseudobulk")
    assert atlas.entry == "scripts/subtype_reference/03_atlas_pseudobulk.py"
    assert _flag_value(atlas, "--input").endswith("mouse_kidney_atlas/mka.h5ad")
    assert _flag_value(atlas, "--output-dir") == (
        "data/processed/v13_dct_reference/mouse_kidney_atlas"
    )
    assert _flag_value(atlas, "--source-study-column") == "Origin"
    assert _flag_value(atlas, "--cell-type-column") == "Celltype_finest"
    assert _flag_value(atlas, "--fallback-cell-type-column") == "Celltype"
    tiers = ca.stage("marker-tiers")
    assert tiers.entry == "scripts/v13/compartment_adversarial_audit.py"
    assert list(tiers.args) == [
        "--config", "config/v13_compartment_adversarial_audit.yaml", "prepare",
    ]


def test_reference_builds_run_before_every_consumer_when_optional_stages_are_included():
    ca = REGISTRY["clinical-axes"]
    order = [s.id for s in ca.ordered_stages(with_optional=True)]
    reference_builds = ("gene-map-reconstruction", "atlas-pseudobulk", "marker-tiers")
    builds = sorted(order.index(stage_id) for stage_id in reference_builds)
    assert builds == [0, 1, 2], "rebuilt inputs must exist before any stage reads them"


def test_clinical_axes_stage_inputs_are_declared():
    ca = REGISTRY["clinical-axes"]
    assert "data/processed/resources/id_map.tsv" in ca.requires_data
    tiers_stages = {
        "compartment-context", "podocyte-disjoint", "podocyte-disjoint-c3h",
        "podxl-nid1", "gene-coherence", "strict-matching", "barrier-adjustment",
        "podocyte-scaffold-specificity",
    }
    for stage_id in tiers_stages:
        assert TIERS in ca.stage(stage_id).requires, stage_id
    for stage_id in ("matched-null", "strict-matching"):
        assert ATLAS_EXPR in ca.stage(stage_id).requires, stage_id
    osd771_stages = tiers_stages | {"cross-mission", "sensitivities", "matched-null"}
    for stage_id in osd771_stages:
        assert OSD771_META in ca.stage(stage_id).requires, stage_id


@pytest.mark.parametrize("rev_id", sorted(REGISTRY))
def test_stage_requires_are_literal_repo_relative_paths(rev_id):
    for stage in REGISTRY[rev_id].stages:
        for path in stage.requires:
            assert "{" not in path and not Path(path).is_absolute(), (rev_id, stage.id, path)


def test_v13_compartment_audit_runs_postprocess_and_declares_its_legacy_inputs():
    st = REGISTRY["v13-compartment"].stage("audit")
    args = list(st.args)
    assert args[:3] == ["--config", "{spec_config}", "postprocess"], (
        "the audit script needs a subcommand right after --config"
    )
    assert _flag_value(st, "--exact-dir") == "data/archive/results/run_20260802_v13_compartment_exact"
    assert "data/archive/results/run_20260802_v13_compartment_exact" in st.requires
    assert "data/results/run_20260729_v13_continuous_phospho_exact_final" in st.requires
    assert any("grey60_adversarial" in r for r in st.requires)
    assert any("layer_block_shift" in r for r in st.requires)
    assert any("osd462_slc12a3_stk39_phosphoform_provenance" in r for r in st.requires)


# --- conftest.requires_paths ------------------------------------------------------

def test_requires_paths_marker_skips_only_when_something_is_missing():
    present = requires_paths(PRESENT_PATH, REPO_ROOT / "README.md")
    absent = requires_paths(PRESENT_PATH, ABSENT_PATH)
    assert present.name == "skipif" and present.args == (False,)
    assert absent.name == "skipif" and absent.args == (True,)
    assert ABSENT_PATH in absent.kwargs["reason"]
    assert requires_paths().args == (False,)
    assert missing_paths(PRESENT_PATH, [ABSENT_PATH]) == [ABSENT_PATH]
