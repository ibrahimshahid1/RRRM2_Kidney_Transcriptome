"""Tests for completed-run provenance finalization."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "v13"
    / "finalize_inference_provenance.py"
)
SPEC = importlib.util.spec_from_file_location("v13_provenance", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_provenance_hashes_completed_outputs_without_self_hash(tmp_path: Path):
    (tmp_path / "manifest.json").write_text(
        json.dumps(
            {
                "analysis_id": "test",
                "inputs": {},
                "permutation": {"mode": "exact", "n_assignments_run": 36},
                "condition_reporter_position_confounded": True,
            }
        )
        + "\n"
    )
    (tmp_path / "claim_tier.tsv").write_text("claim_tier\nneither\n")
    (tmp_path / "result.tsv").write_text("x\n1\n")
    manifest = MODULE.build_provenance_manifest(tmp_path)
    assert manifest["analysis_id"] == "test"
    assert manifest["condition_reporter_position_confounded"] is True
    assert "manifest.json" in manifest["output_hashes_excluding_this_manifest"]
    assert "result.tsv" in manifest["output_hashes_excluding_this_manifest"]
    assert MODULE.MANIFEST_NAME not in manifest[
        "output_hashes_excluding_this_manifest"
    ]


def test_provenance_refuses_incomplete_run(tmp_path: Path):
    (tmp_path / "manifest.json").write_text('{"analysis_id":"test"}\n')
    try:
        MODULE.build_provenance_manifest(tmp_path)
    except FileNotFoundError as exc:
        assert "claim_tier.tsv" in str(exc)
    else:
        raise AssertionError("Incomplete run should be refused")
