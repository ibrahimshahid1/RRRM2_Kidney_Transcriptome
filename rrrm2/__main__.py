"""Allow `python -m rrrm2` without installing the package."""

from __future__ import annotations

import sys

from .cli import main

if __name__ == "__main__":
    sys.exit(main())
