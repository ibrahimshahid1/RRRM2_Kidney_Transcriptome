"""The committed frozen-resource bundle restores byte-identical files."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
BUNDLE = REPO / "resources/frozen/2026-10-06"


def test_manifest_lists_existing_packed_files():
    rows = (BUNDLE / "MANIFEST.tsv").read_text().splitlines()[1:]
    assert len(rows) == 12
    for row in rows:
        packed, target, compression, digest = row.split("\t")
        assert (BUNDLE / packed).exists(), packed
        assert target.startswith("data/processed/")
        assert compression in {"gzip", "none"}
        assert len(digest) == 64


def test_restore_into_empty_root_verifies_every_sha256(tmp_path):
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts/restore_frozen_resources.py"), "--root", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "12/12 resources verified" in result.stdout
    assert (tmp_path / "data/processed/resources/id_map.tsv").exists()
