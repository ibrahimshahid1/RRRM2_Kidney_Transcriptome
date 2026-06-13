"""RRRM-2 Validation Module"""

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
