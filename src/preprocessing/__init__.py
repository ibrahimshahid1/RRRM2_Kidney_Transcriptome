"""RRRM-2 Preprocessing Module"""

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
