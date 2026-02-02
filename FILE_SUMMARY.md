# Repository File Summary

## Reorganized Module Structure (src/)

### src/preprocessing/
Phase 0-1: Deconvolution, VST Normalization, Global Residualization
| File | Description |
|------|-------------|
| `__init__.py` | Module exports |
| `deconvolution.R` | MuSiC cell-type deconvolution with Tabula Muris Senis reference |
| `residualization.R` | DESeq2 VST + SVA confounder regression |
| `data_alignment.py` | Metadata-to-counts alignment utilities |
| `export_counts.py` | Raw count export utilities |
| `export_phase1.R` | R-to-Python data export |

### src/networks/
Phase 2-3: Shared Topology + LIONESS + Edge Regression + node2vec + Procrustes
| File | Description |
|------|-------------|
| `__init__.py` | Module exports |
| `shared_topology.py` | Cell-standardized skeleton (Ledoit-Wolf shrinkage + top-k) |
| `lioness.py` | LIONESS sample-specific edge weights |
| `edge_regression.py` | Edge-wise regression + limma EB + predicted networks |
| `embeddings.py` | PecanPy node2vec embeddings (multi-seed) |
| `procrustes.py` | Procrustes alignment + rewiring quantification |

### src/statistics/
Phase 5-6: Silent Shifters, Uncertainty Estimation, Full Regression
| File | Description |
|------|-------------|
| `__init__.py` | Module exports |
| `silent_shifters.py` | Silent shifter identification (high rewiring + low DE) |
| `interaction_metrics.py` | Interaction persistence analysis across contrasts |
| `permutation_bootstrap.py` | Permutation tests + bootstrap CIs (n=2000) |
| `full_regression.py` | Full factorial edge regression (n=80) |
| `differential_expression.R` | Gene-level DE tables (limma) |

### src/enrichment/
Phase 7: Biological Grounding
| File | Description |
|------|-------------|
| `__init__.py` | Module exports |
| `biological_grounding.py` | Gene set enrichment (GO/KEGG/Reactome) + module analysis |

### src/visualization/
Publication-ready plots and diagnostics
| File | Description | Source |
|------|-------------|--------|
| `__init__.py` | Module exports | NEW |
| `publication_plots.py` | All-phase publication plots | `plot_output.py` |
| `network_diagnostics.py` | Skeleton visualization | `plot_skeleton_diagnostics.py` |

### src/validation/
Leakage-safe cross-validation (Section 5 of methodology)
| File | Description |
|------|-------------|
| `__init__.py` | Module exports |
| `cross_validation.py` | 5-fold stratified CV with fold-wise feature engineering |
| `sample_features.py` | Mouse-level feature extraction (pathway strength, node strength, PCA) |

### src/utils.py
Common utility functions (unchanged)

---

## Scripts Remaining in scripts/

### Pipeline Orchestration
- **`src/run_all_phases.py`** - ⭐ Main pipeline orchestrator (runs phases 0-7)
- `run_full_pipeline.py` - Legacy placeholder (not currently used)

### Phase-Specific Scripts
- `DESeq2.R` - VST normalization source
- `run_deconvolution.R` - Full deconvolution pipeline
- `phase2_edge_regression.py` - Edge-wise regression implementation
- `phase3_node2vec_embedding.py` - PecanPy node2vec implementation
- `phase3_procrustes_rewiring.py` - Procrustes alignment + rewiring

### Utility Scripts
- `pick_anchors.py`, `make_consensus_anchors.py` - Procrustes anchor selection
- `plot_output.py`, `plot_skeleton_diagnostics.py` - Visualization utilities

---

## Research Phase → Module Mapping

| Phase | Description | Module | Key Files |
|-------|-------------|--------|-----------|
| 0 | Deconvolution | `src/preprocessing/` | `deconvolution.R` |
| 1 | VST + Residualization | `src/preprocessing/` | `residualization.R`, `data_alignment.py` |
| 2 | Network Construction | `src/networks/` | `shared_topology.py`, `lioness.py`, `edge_regression.py` |
| 3 | Embeddings + Alignment | `src/networks/` | `embeddings.py`, `procrustes.py` |
| 5 | Silent Shifters | `src/statistics/` | `silent_shifters.py`, `interaction_metrics.py` |
| 6 | Uncertainty + Full Reg | `src/statistics/` | `permutation_bootstrap.py`, `full_regression.py` |
| 7 | Biological Grounding | `src/enrichment/` | `biological_grounding.py` |
| - | Validation | `src/validation/` | `cross_validation.py`, `sample_features.py` |

> **Note:** Phase 4 (edge regression) is now integrated into Phase 2. The pipeline runs phases {0, 1, 2, 3, 5, 6, 7}.

---

## Usage Examples

```bash
# Run full pipeline
python src/run_all_phases.py

# Run specific phases
python src/run_all_phases.py --phases 2 3 5

# Run without R dependencies
python src/run_all_phases.py --skip-r

# Run as modules (from repo root)
python -m src.networks.shared_topology --max_genes 1200 --topk 80
python -m src.networks.lioness
python -m src.networks.edge_regression

# Import in Python
from src.networks import shared_topology, lioness
from src.statistics import silent_shifters
from src.validation import cross_validation
```
