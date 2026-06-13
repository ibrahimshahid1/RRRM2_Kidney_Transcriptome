"""RRRM-2 Statistics Module"""

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
