"""
RRRM-2 Validation Module

Leakage-safe cross-validation and sample-level feature extraction.

This module contains:
- cross_validation: Leakage-safe CV framework
- sample_features: Mouse-level feature extraction for LIONESS networks
"""

from pathlib import Path

# Module directory
MODULE_DIR = Path(__file__).parent

# Expose submodules
try:
    from . import cross_validation
except ImportError:
    pass

try:
    from . import sample_features
except ImportError:
    pass
