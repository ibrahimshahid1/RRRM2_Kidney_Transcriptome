"""CLI behaviour: listing, showing, dry-run planning, and the retired-revision guard."""

from __future__ import annotations

from dataclasses import replace
import json

import pytest
import yaml

from rrrm2 import cli
from rrrm2.cli import main
from rrrm2.registry import load_revision


def test_list_runs_clean(capsys):
    assert main(["list"]) == 0
    assert "registered revisions" in capsys.readouterr().out


def test_list_json_is_parseable(capsys):
    assert main(["list", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert {r["id"] for r in payload} >= {"v13-phospho", "clinical-axes"}


def test_show_reports_question_and_conclusion(capsys):
    assert main(["show", "v13-phospho"]) == 0
    out = capsys.readouterr().out
    assert "QUESTION" in out and "CONCLUSION" in out and "STAGES" in out


def test_show_unknown_revision_lists_valid_ids(capsys):
    assert main(["show", "no-such-revision"]) == 2
    assert "known revisions" in capsys.readouterr().err


def test_dry_run_plans_without_writing(tmp_path, capsys):
    assert main(["run", "v13-phospho", "--dry-run", "--run-dir", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "[dry-run]" in out
    assert "run_continuous_phospho_inference.py" in out
    assert not list(tmp_path.iterdir()), "dry-run must not create artifacts"


def test_stage_selection_pulls_in_dependencies(tmp_path, capsys):
    main(["run", "v13-phospho", "--stage", "provenance", "--dry-run",
          "--run-dir", str(tmp_path)])
    out = capsys.readouterr().out
    assert "infer, report, provenance" in out


def test_only_flag_runs_a_single_stage(tmp_path, capsys):
    main(["run", "v13-phospho", "--stage", "provenance", "--only", "--dry-run",
          "--run-dir", str(tmp_path)])
    assert "stages    provenance\n" in capsys.readouterr().out


def test_retired_revision_is_refused_without_override(tmp_path, capsys):
    code = main(["run", "network-rewiring", "--run-dir", str(tmp_path)])
    assert code == 3
    assert "refusing to run" in capsys.readouterr().err


def test_optional_stages_are_excluded_by_default(tmp_path, capsys):
    main(["run", "v13-phospho", "--dry-run", "--run-dir", str(tmp_path)])
    default = capsys.readouterr().out
    main(["run", "v13-phospho", "--dry-run", "--with-optional",
          "--run-dir", str(tmp_path)])
    expanded = capsys.readouterr().out
    assert "intensity-confound" not in default
    assert "intensity-confound" in expanded


def test_skip_r_drops_rscript_stages(tmp_path, capsys):
    main(["run", "subtype-reference", "--dry-run", "--skip-r",
          "--run-dir", str(tmp_path)])
    out = capsys.readouterr().out
    assert "Rscript" not in out
    assert "inventory" in out


def test_doctor_reports_every_revision(capsys):
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "revision runnability" in out and "v13-phospho" in out


def test_panels_list_does_not_crash_when_empty_or_populated(capsys):
    assert main(["panels", "list"]) == 0


def test_podocyte_disjoint_dry_run_targets_its_own_subfolder(tmp_path, capsys):
    assert main(["run", "clinical-axes", "--stage", "podocyte-disjoint", "--only",
                 "--dry-run", "--run-dir", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "stages    podocyte-disjoint\n" in out
    assert f"--results {tmp_path}/disjoint_podocyte_family " in out
    assert "--add-podocyte-disjoint-variants" in out
    assert not list(tmp_path.iterdir()), "dry-run must not create artifacts"


def test_default_clinical_axes_dry_run_plans_the_new_sensitivity_stages(tmp_path, capsys):
    assert main(["run", "clinical-axes", "--dry-run", "--run-dir", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert f"--results {tmp_path}/disjoint_podocyte_family_c3h " in out
    assert "--osd253-strain C3H/HeJ" in out
    assert f"--results {tmp_path}/podxl_nid1_forced_sensitivity " in out
    assert "--permutations 100000 --seed 0" in out  # K0 pins its seed
    assert "atlas-pseudobulk" not in out, "reference builds are optional"


def test_reference_builds_are_only_planned_with_optional(tmp_path, capsys):
    main(["run", "clinical-axes", "--dry-run", "--with-optional", "--run-dir", str(tmp_path)])
    out = capsys.readouterr().out
    assert "03_atlas_pseudobulk.py" in out
    assert "compartment_adversarial_audit.py --config config/v13_compartment_adversarial_audit.yaml prepare" in out
    assert out.index("atlas-pseudobulk") < out.index("cross-mission")


def test_show_json_and_text_report_stage_requires(capsys):
    assert main(["show", "clinical-axes", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    stage = next(s for s in payload["stages"] if s["id"] == "compartment-context")
    assert "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv" in stage["requires"]
    assert main(["show", "clinical-axes"]) == 0
    assert "requires data/raw/metadata/s_OSD-771.txt" in capsys.readouterr().out


# --- doctor / preflight with stage-level requires (temporary registry) -----------

ABSENT = "data/definitely/not/here/input.tsv"
PRESENT = "config/revisions/_schema.yaml"


@pytest.fixture
def temp_registry(tmp_path, monkeypatch):
    """Point the CLI at one temporary revision with a blocked and an optional stage."""
    path = tmp_path / "alpha.yaml"
    stage = {"runner": "python", "entry": "scripts/audit_fix_status.py"}
    path.write_text(yaml.safe_dump({
        "id": "alpha", "title": "Alpha", "status": "complete",
        "results_root": "data/results/alpha", "question": "q?",
        "stages": [
            {**stage, "id": "fine", "title": "fine", "requires": [PRESENT]},
            {**stage, "id": "blocked-one", "title": "blocked", "requires": [ABSENT]},
            {**stage, "id": "extra", "title": "extra", "optional": True,
             "requires": [ABSENT + ".extra"]},
        ],
    }))
    rev = load_revision(path)
    monkeypatch.setattr(cli, "load_registry", lambda: {"alpha": rev})
    monkeypatch.setattr(cli, "get_revision", lambda rev_id: rev)
    return rev


def test_doctor_exits_zero_and_names_blocked_stages(temp_registry, capsys):
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "alpha" in out and "blocked" in out
    assert f"stage blocked-one: missing required input {ABSENT}" in out
    assert "blocked stages: blocked-one" in out
    assert "0/1 revisions runnable" in out
    # the optional stage is a warning, never a blocked stage
    assert f"stage extra: missing required input {ABSENT}.extra" in out
    assert "blocked stages: blocked-one, extra" not in out


def test_doctor_treats_a_revision_with_only_optional_gaps_as_runnable(
        temp_registry, monkeypatch, capsys):
    only_optional = replace(
        temp_registry,
        stages=tuple(s for s in temp_registry.stages if s.id != "blocked-one"),
    )
    monkeypatch.setattr(cli, "load_registry", lambda: {"alpha": only_optional})
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "alpha                ok" in out
    assert "blocked stages" not in out
    assert "optional stage, skipped by default" in out


def test_real_doctor_lists_every_revision_and_still_exits_zero(capsys):
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    for rev_id in ("clinical-axes", "v13-compartment", "v13-phospho"):
        assert rev_id in out
    assert "revisions runnable here" in out


def test_run_refuses_a_stage_with_a_missing_required_input(temp_registry, tmp_path, capsys):
    code = main(["run", "alpha", "--stage", "blocked-one", "--only",
                 "--run-dir", str(tmp_path / "out")])
    assert code == 2
    err = capsys.readouterr().err
    assert f"stage blocked-one: missing required input {ABSENT}" in err


def test_dry_run_still_prints_the_plan_when_an_input_is_missing(temp_registry, tmp_path, capsys):
    code = main(["run", "alpha", "--stage", "blocked-one", "--only", "--dry-run",
                 "--run-dir", str(tmp_path / "out")])
    captured = capsys.readouterr()
    assert code == 0
    assert "preflight failed" in captured.err
    assert "[dry-run]" in captured.out


def test_selecting_an_optional_stage_explicitly_makes_its_inputs_required(
        temp_registry, tmp_path, capsys):
    code = main(["run", "alpha", "--stage", "extra", "--only",
                 "--run-dir", str(tmp_path / "out")])
    assert code == 2
    assert f"stage extra: missing required input {ABSENT}.extra" in capsys.readouterr().err
