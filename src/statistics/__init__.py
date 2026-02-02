"""
RRRM-2 Statistics Module

Phase 5: Rewiring Metrics + Silent Shifters
Phase 6: Uncertainty/Null Models (Permutation + Bootstrap)

This module contains:
- silent_shifters: Silent shifter identification (high rewiring, low DE)
- interaction_metrics: Interaction persistence analysis
- permutation_bootstrap: Permutation/bootstrap uncertainty estimation
- full_regression: Full edge regression (all 80 samples)

Note: R script (differential_expression.R) is called via Rscript.
"""

from pathlib import Path

# Module directory
MODULE_DIR = Path(__file__).parent

# Expose submodules
try:
    from . import silent_shifters
except ImportError:
    pass

try:
    from . import interaction_metrics
except ImportError:
    pass

try:
    from . import permutation_bootstrap
except ImportError:
    pass

try:
    from . import full_regression
except ImportError:
    pass
