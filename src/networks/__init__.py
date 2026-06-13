"""RRRM-2 Networks Module"""

from pathlib import Path

# Module directory
MODULE_DIR = Path(__file__).parent

# Expose submodules
try:
    from . import shared_topology
except ImportError:
    pass

try:
    from . import lioness
except ImportError:
    pass

try:
    from . import edge_regression
except ImportError:
    pass

try:
    from . import embeddings
except ImportError:
    pass

try:
    from . import procrustes
except ImportError:
    pass
