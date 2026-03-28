# 🚀 Age-Dependent Network Rewiring in Mouse Kidney After Spaceflight

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![NASA GeneLab](https://img.shields.io/badge/Data-NASA%20OSD--771-red.svg)](https://genelab.nasa.gov/data/study?acc=OSD-771)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

> **A node2vec-Driven Transcriptomic Framework with Cell-Type Deconvolution, Shared-Topology Network Construction, and Leakage-Safe Validation**

Computational pipeline for analyzing transcriptomic network rewiring in mouse kidney tissue from the NASA RRRM-2 spaceflight mission (OSD-771).

---

## 🎯 Overview

This pipeline implements a rigorous graph embedding approach to identify **"silent shifters"** — genes whose network context changes substantially despite minimal differential expression. We analyze the full factorial design (n=80): **Age × Arm × Environment Group** (2 × 2 × 4).

### Key Innovations

1. **Cell-Standardized Shared Topology**: Prevent Simpson's paradox by standardizing within experimental cells before pooling for edge selection
2. **LIONESS Sample-Specific Networks**: Estimate per-mouse edge weights on a fixed skeleton
3. **Edge-Wise Regression**: Model edge weights over full factorial design with empirical Bayes
4. **node2vec Embeddings**: 128-dimensional graph embeddings with Procrustes alignment
5. **Leakage-Safe Validation**: Fold-wise feature engineering for honest predictive assessment

---

## 📊 Dataset (NASA GeneLab OSD-771)

**Rodent Research Reference Mission 2 (RRRM-2)**  
Female C57BL/6NTac mice, whole kidney RNA-seq

### Full Factorial Design (n=80)
| Factor | Levels |
|--------|--------|
| **Age** | Young (16 weeks), Old (34 weeks) |
| **Arm** | ISS-T (terminal on ISS), LAR (live animal return) |
| **Environment Group** | FLT (flight), HGC (hardware control), VIV (vivarium), BSL (basal) |

5 mice per age within each Age × Arm × Group cell.

---

## 🔬 Pipeline Phases

### Phase 0: Preprocessing & Deconvolution
- VST normalization (DESeq2)
- Cell-type deconvolution using **MuSiC** with Tabula Muris Senis kidney atlas reference
- Nephron segment proportion estimation (DCT, PT, CD, Glom, Immune, Stroma)
- QC and outlier detection

### Phase 1: Global Residualization
- SVA for surrogate variables
- Regression: `Y ~ batch + SVs + CLR(cell_props)`
- Output: `R_tech` (biology preserved), `R_all` (fully residualized)

### Phase 2: Network Construction (Shared Topology + LIONESS + Edge Regression)
- **Cell-standardization**: `(R - μ_cell) / σ_cell` per gene within each Age×Arm×Group cell
- Pool all samples → build fixed skeleton **E** using Ledoit-Wolf shrinkage + top-k partial correlations
- **LIONESS** sample-specific edge weights on fixed **E**: `w_e(s) = N·w_e(all) - (N-1)·w_e(-s)`
- Fisher z-transform: `z = atanh(r)`
- **Edge-wise regression**: Model each edge `z_e ~ Age + Arm + Group + interactions + covariates`
- Empirical Bayes variance moderation (limma-style) → predicted contrast networks

### Phase 3: node2vec Embeddings & Procrustes Alignment
- Multi-seed **PecanPy** node2vec (10 seeds, d=128, p=0.25, q=4)
- Orthogonal Procrustes alignment using pre-registered anchor genes
- Consensus embedding across seeds

### Phase 5: Silent Shifters & Interaction Metrics
- **Silent shifters**: high rewiring (top 10%, FDR<0.1) + low DE (|log₂FC|<0.3, FDR>0.2)
- Interaction persistence analysis across contrasts

### Phase 6: Uncertainty Estimation & Full Regression
- Permutation tests (n=2000, stratified within Age/Arm)
- Bootstrap CIs (n=2000) for rewiring metrics
- Full factorial edge regression (n=80)
- Gene-level differential expression (limma)

### Phase 7: Biological Grounding
- Gene set enrichment analysis (GO, KEGG, Reactome)
- Pathway-level rewiring quantification
- Module analysis via k-means on embeddings

### Phase 8: Leakage-Safe Predictive Validation
- Stratified K-fold CV (FLT vs GC, stratified by Age×Arm)
- Fold-wise skeleton E built on training data only
- LIONESS features extracted per fold; no test data leakage
- Classification performance (accuracy, AUC) with LogisticRegression and RandomForest

---

## 📦 Installation

### Option 1: Conda (Recommended)
```bash
git clone https://github.com/ibrahimshahid1/RRRM2_Kidney_Transcriptome.git
cd RRRM2_Kidney_Transcriptome

conda env create -f environment.yml
conda activate rrrm2_kidney

pip install -e .
```

### Option 2: pip
```bash
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

pip install -r requirements.txt
pip install -e .
```

---

## 📥 Data Download

Due to file size limitations, input data files are not included in this repository. Download the following before running the pipeline:

### 1. NASA GeneLab OSD-771 (Bulk RNA-seq)

Download from [NASA GeneLab OSD-771](https://genelab.nasa.gov/data/study?acc=OSD-771):

| File | Destination |
|------|-------------|
| `GLDS-771_rna_seq_ERCCnorm_counts.csv` | `data/raw/counts/` |
| `GLDS-771_metadata_OSD-771-samples.csv` | `data/raw/metadata/` |
| `GLDS-771_annotations.csv` | `data/raw/metadata/` |

### 2. Single-Cell Reference Atlas (for Deconvolution)

The MuSiC deconvolution requires a Tabula Muris Senis kidney reference atlas. Download from [CZI CellxGene](https://cellxgene.cziscience.com/):

```bash
# Search for: "Tabula Muris Senis" + "kidney" + "female"
# Download the .h5ad file and place in:
mkdir -p data/external/single_cell_atlases/
# → data/external/single_cell_atlases/tms_kidney_female_ALLDATASETS_counts_innerGenes.h5ad
```

**Required files** (place in `data/external/single_cell_atlases/`):
- `tms_kidney_female_ALLDATASETS_counts_innerGenes.h5ad` (~200MB)
- `tms_kidney_female_ALLDATASETS_obs.csv` (cell metadata)
- `tms_kidney_female_ALLDATASETS_var.csv` (gene metadata)

### 3. Verify Data Structure

After downloading, your data directory should look like:
```
data/
├── raw/
│   ├── counts/
│   │   └── GLDS-771_rna_seq_ERCCnorm_counts.csv
│   └── metadata/
│       ├── GLDS-771_metadata_OSD-771-samples.csv
│       └── GLDS-771_annotations.csv
└── external/
    └── single_cell_atlases/
        ├── tms_kidney_female_ALLDATASETS_counts_innerGenes.h5ad
        ├── tms_kidney_female_ALLDATASETS_obs.csv
        └── tms_kidney_female_ALLDATASETS_var.csv
```

## 🚀 Quick Start

### Run Full Pipeline
```bash
# Run all phases
python src/run_all_phases.py

# Run specific phases
python src/run_all_phases.py --phases 2 3 5

# Skip R-dependent steps (if R not available)
python src/run_all_phases.py --skip-r

# Dry run (show commands without executing)
python src/run_all_phases.py --dry-run
```

### Configuration
Edit `config/hyperparameters.yaml` to adjust:
- node2vec parameters (dimensions=128, p=0.25, q=4.0)
- Network topology method (top-k=50 neighbors, Ledoit-Wolf shrinkage)
- Statistical thresholds (FDR, percentiles, bootstrap/permutation iterations)

### Key Gene Sets
Pre-configured in `config/gene_sets.yaml`:
- **DCT/NCC-WNK pathway**: WNK1, WNK4, STK39, SLC12A3, KCNJ10
- **Positive controls**: ECM remodeling, oxidative stress, calcium handling

---

## 📁 Repository Structure

```
RRRM2_Kidney_Transcriptome/
├── config/
│   ├── metadata_design.yaml     # Full factorial design specification
│   ├── anchor_genes.yaml        # Pre-registered Procrustes alignment anchors
│   ├── gene_sets.yaml           # DCT/NCC-WNK & pathway genes
│   └── hyperparameters.yaml     # Pipeline parameters
├── src/
│   ├── run_all_phases.py        # ⭐ Main pipeline orchestrator
│   ├── preprocessing/
│   │   ├── deconvolution.R      # MuSiC cell-type deconvolution
│   │   ├── residualization.R    # DESeq2 VST + SVA residualization
│   │   ├── data_alignment.py    # Metadata-to-counts alignment
│   │   └── export_phase1.R      # R-to-Python export utilities
│   ├── networks/
│   │   ├── shared_topology.py   # Cell-standardized skeleton (Ledoit-Wolf)
│   │   ├── lioness.py           # Sample-specific edge weights
│   │   ├── edge_regression.py   # Edge-wise modeling + limma EB
│   │   ├── embeddings.py        # PecanPy node2vec embeddings
│   │   └── procrustes.py        # Procrustes alignment + rewiring
│   ├── statistics/
│   │   ├── silent_shifters.py   # Silent shifter identification
│   │   ├── interaction_metrics.py # Interaction persistence
│   │   ├── permutation_bootstrap.py # Uncertainty estimation
│   │   └── full_regression.py   # Full factorial regression (n=80)
│   ├── enrichment/
│   │   └── biological_grounding.py # GO/KEGG/Reactome enrichment
│   └── validation/
│       ├── cross_validation.py  # Leakage-safe CV framework
│       └── sample_features.py   # Mouse-level feature extraction
├── scripts/                     # Standalone utility scripts
├── data/results/                # Pipeline outputs by phase
├── METHODOLOGY.md               # Detailed methods
└── README.md
```

---

## 🔍 Key Outputs (Per Run)

Each pipeline run produces versioned outputs in `data/results/<run_id>/`:
- **Rewiring scores**: per-gene cosine-distance shifts for each contrast
- **Permutation p-values**: gene-level significance with BH-FDR
- **Silent shifters**: genes with high rewiring but low differential expression
- **Cross-validation results**: leakage-safe accuracy and AUC for FLT vs GC classification
- **Enrichment tables**: pathway enrichment for rewired gene sets
- **Publication figures**: volcano plots, rewiring distributions, QQ plots

---

## 📖 Citation

If you use this pipeline, please cite:

```bibtex
@misc{rrrm2_kidney_rewiring_2026,
  title={Age-Dependent Network Rewiring in Mouse Kidney After Spaceflight},
  author={Shahid, Ibrahim},
  year={2026},
  note={NASA GeneLab OSD-771 analysis pipeline}
}
```

---

## 📚 References

1. Kuijjer et al. (2019) LIONESS. *iScience*.
2. Grover & Leskovec (2016) node2vec. *KDD*.
3. NASA GeneLab Consortium. RRRM-2 Kidney Dataset (OSD-771).
4. Smyth (2004) limma empirical Bayes. *Stat Appl Genet Mol Biol*.

---

## 📧 Contact

**Ibrahim Shahid**  
GitHub: [@ibrahimshahid1](https://github.com/ibrahimshahid1)

---

## 📄 License

MIT License - see [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgments

- NASA GeneLab for open data access
- RRRM-2 research team
- Single-cell kidney atlas contributors (Tabula Muris Senis, Park et al.)