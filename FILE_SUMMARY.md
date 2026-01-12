# Repository File Summary

## Created Template/Scaffold Files

### Core Modules (src/)

#### Preprocessing
- ✅ `src/preprocessing/residualization.py` - Global confounder regression with SVA
- ✅ `src/preprocessing/deconvolution.py` - Cell-type proportion estimation
- ⚠️  `src/preprocessing/normalization.py` - Existing, needs update for VST
- ⚠️  `src/preprocessing/qc.py` - Existing, needs update for variance partitioning

#### Networks
- ✅ `src/networks/shared_topology.py` - Cell-standardized edge selection (KEY INNOVATION)
- ✅ `src/networks/lioness.py` - LIONESS sample-specific networks
- ✅ `src/networks/edge_regression.py` - Edge-wise factorial modeling
- ✅ `src/networks/embeddings.py` - node2vec + Procrustes alignment
- ⚠️  `src/networks/graph_builder.py` - Exists, may need updates for signed edges
- ⚠️  `src/networks/metrics.py` - Exists, centrality/module detection

#### Statistics
- ✅ `src/statistics/rewiring_metrics.py` - Cosine distance & silent shifters
- ⚠️  `src/statistics/bootstrap.py` - Exists, needs fixed-topology implementation
- ⚠️  `src/statistics/permutation.py` - Exists, needs leakage-safe updates

#### Validation (NEW MODULE)
- ✅ `src/validation/cross_validation.py` - Leakage-safe CV framework
- ✅ `src/validation/sample_features.py` - Mouse-level feature extraction
- ✅ `src/validation/__init__.py`

#### Utils
- ✅ `src/utils.py` - Updated with all starter code functions

### Configuration Files (config/)
- ✅ `config/metadata_design.yaml` - Full n=80 factorial specification
- ✅ `config/anchor_genes.yaml` - Pre-registered Procrustes anchors
- ✅ `config/gene_sets.yaml` - DCT/NCC-WNK + spaceflight pathways
- ✅ `config/hyperparameters.yaml` - Complete pipeline parameters

### Scripts
- ✅ `scripts/run_full_pipeline.py` - Master orchestration (with placeholders)
- ⚠️  Individual phase scripts (phase0-7) - Not yet created, use master script

### Documentation
- ✅ `METHODOLOGY.md` - Complete detailed methodology
- ✅ `README.md` - Updated with new research framework

## Status Legend
- ✅ Newly created template file with full scaffold
- ⚠️ Existing file that may need updates
- ❌ Not yet created

## Next Steps for User

1. **Test Imports**
   ```bash
   cd c:\Users\Ibrah\Documents\GitHub\RRRM2_Kidney_Transcriptome
   python -c "from src.preprocessing import residualization, deconvolution; print('✓ Imports OK')"
   ```

2. **Load Configuration**
   ```bash
   python -c "import yaml; print(yaml.safe_load(open('config/metadata_design.yaml')))"
   ```

3. **Add Actual Data**
   - Place OSD-771 raw counts in `data/raw/`
   - Update placeholder data loaders in phase scripts

4. **Implement Production Features** (Priority Order)
   - Partial correlation / graphical lasso in `shared_topology.py`
   - Fast LIONESS optimization in `lioness.py`
   - Full empirical Bayes in `edge_regression.py`
   -  PecanPy integration in `embeddings.py`

5. **Run Tests**
   - Create synthetic test data (small n=16, 100 genes)
   - Validate pipeline end-to-end on toy data

## Key Design Decisions

1. **Fixed Topology E**: Same edges across all conditions, different weights
2. **Cell-Standardization**: Within Age×Arm×Group cells before pooling
3. **Fisher z**: For approximate normality in edge regression
4. **Leakage-Safe CV**: All feature engineering within training folds
