#!/usr/bin/env python3
"""CLI wrapper for the frozen v13 compartment adversarial audit."""

from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.v13.compartment_adversarial_audit import main  # noqa: E402


if __name__ == "__main__":
    raise SystemExit(main())
