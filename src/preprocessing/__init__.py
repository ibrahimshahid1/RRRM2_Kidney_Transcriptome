"""
RRRM-2 Preprocessing Module

Phase 0: Pre-processing, Deconvolution, QC
Phase 1: Global Residualization (VST + SVA)

This module contains:
- data_alignment: Metadata-to-counts alignment
- export_counts: Raw count export utilities

Note: R scripts (residualization.R, deconvolution.R) are called via Rscript
and are stored alongside Python modules for organizational purposes.
"""

from pathlib import Path

# Module directory for locating R scripts
MODULE_DIR = Path(__file__).parent

# Expose Python modules when they exist
try:
    from . import data_alignment
except ImportError:
    pass

try:
    from . import export_counts
except ImportError:
    pass
