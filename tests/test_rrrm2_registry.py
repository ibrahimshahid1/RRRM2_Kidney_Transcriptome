"""Registry invariants: every revision parses, resolves, and points at real entries."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

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
