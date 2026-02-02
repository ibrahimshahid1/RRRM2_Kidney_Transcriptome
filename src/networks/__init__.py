"""
RRRM-2 Networks Module

Phase 2: Shared Topology Construction + LIONESS + Edge Regression
Phase 3: node2vec Embeddings + Procrustes Alignment

This module contains:
- shared_topology: Cell-standardized shared skeleton construction
- lioness: LIONESS sample-specific edge weights
- edge_regression: Edge-wise regression with full factorial design
- embeddings: PecanPy node2vec embeddings
- procrustes: Procrustes alignment + rewiring metrics
"""

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
