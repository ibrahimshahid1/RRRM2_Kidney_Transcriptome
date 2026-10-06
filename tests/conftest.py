"""Shared pytest helpers.

``requires_paths`` is the one sanctioned way to skip a data-dependent test when
a fresh clone lacks gitignored inputs (see plan section 0: data-dependent tests
must be skipped, never failed, when the data is absent). Import it with::

    from conftest import requires_paths

    @requires_paths("data/processed/resources/id_map.tsv")
    def test_something_that_reads_real_data(): ...
"""

from __future__ import annotations

from pathlib import Path
import sys
from typing import Iterable

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]

# Existing tests import ``src.*`` and ``scripts.*`` relative to the repository
# root. ``python -m pytest`` puts the root on sys.path, a bare ``pytest`` does not.
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def missing_paths(*paths: str | Path) -> list[str]:
    """Return the entries of ``paths`` that do not exist.

    Relative paths are resolved against the repository root, so a test behaves
    the same no matter which directory pytest was launched from. Entries are
    reported as given, which keeps the skip reason readable.
    """
    missing: list[str] = []
    for path in _flatten(paths):
        candidate = Path(path)
        resolved = candidate if candidate.is_absolute() else REPO_ROOT / candidate
        if not resolved.exists():
            missing.append(str(path))
    return missing


def requires_paths(*paths: str | Path) -> pytest.MarkDecorator:
    """``pytest.mark.skipif`` that skips when any of ``paths`` is absent.

    With no ``paths`` the marker never skips. The reason names the missing
    files so a skipped test explains what a contributor has to download.
    """
    missing = missing_paths(*paths)
    reason = "missing data input(s): " + ", ".join(missing) if missing else "inputs present"
    return pytest.mark.skipif(bool(missing), reason=reason)


def _flatten(items: Iterable) -> Iterable[str | Path]:
    for item in items:
        if isinstance(item, (str, Path)):
            yield item
        else:
            yield from _flatten(item)
