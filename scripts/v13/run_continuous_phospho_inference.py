#!/usr/bin/env python3
"""CLI wrapper for the prospectively frozen v13 phosphoproteomic analysis."""

from __future__ import annotations

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.v13.continuous_phospho_inference import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
