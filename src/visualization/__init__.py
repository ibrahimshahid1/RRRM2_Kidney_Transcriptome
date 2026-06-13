"""RRRM-2 Visualization Module"""

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
