"""
RRRM-2 Visualization Module

Publication-ready plots and network diagnostics.

This module contains:
- publication_plots: Generate publication-ready plots for all phases
- network_diagnostics: Network skeleton visualization and diagnostics
"""

from pathlib import Path

# Module directory
MODULE_DIR = Path(__file__).parent

# Expose submodules
try:
    from . import publication_plots
except ImportError:
    pass

try:
    from . import network_diagnostics
except ImportError:
    pass
