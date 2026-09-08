"""Repository-root discovery and the canonical directory layout."""

from __future__ import annotations

import os
from pathlib import Path

_MARKERS = ("config/revisions", "src", "scripts")


def find_repo_root(start: Path | None = None) -> Path:
    """Walk upward from ``start`` to the directory holding the registry."""
    env = os.environ.get("RRRM2_REPO_ROOT")
    if env:
        return Path(env).resolve()
    here = (start or Path(__file__)).resolve()
    for candidate in (here, *here.parents):
        if all((candidate / marker).exists() for marker in _MARKERS):
            return candidate
    raise FileNotFoundError(
        "cannot locate the repository root; run from inside the repository "
        "or set RRRM2_REPO_ROOT"
    )


REPO_ROOT = find_repo_root()
REVISIONS_DIR = REPO_ROOT / "config" / "revisions"
PANELS_DIR = REPO_ROOT / "config" / "panels"
RESULTS_DIR = REPO_ROOT / "data" / "results"


def relative(path: Path) -> str:
    """Render ``path`` relative to the repository root when it is inside it."""
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)
