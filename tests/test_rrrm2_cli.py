"""CLI behaviour: listing, showing, dry-run planning, and the retired-revision guard."""

from __future__ import annotations

import json

import pytest

from rrrm2.cli import main


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
