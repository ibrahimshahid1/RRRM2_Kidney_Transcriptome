"""
RRRM-2 Enrichment Module

Phase 7: Biological Grounding

This module contains:
- biological_grounding: Gene set enrichment + pathway analysis
"""

from pathlib import Path

# Module directory
MODULE_DIR = Path(__file__).parent

# Expose submodules
try:
    from . import biological_grounding
except ImportError:
    pass
