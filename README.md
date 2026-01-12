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
- Cell-type deconvolution using murine kidney atlases
- Nephron segment proportion estimation (DCT, PT, CD, Glom)
- QC and outlier detection

### Phase 1: Global Residualization
- SVA for surrogate variables
- Regression: `Y ~ batch + SVs + CLR(cell_props)`
- Output: `R_tech` (biology preserved), `R_all` (fully residualized)

### Phase 2: Shared Topology Construction
- **Cell-standardization**: `(R - μ_cell) / σ_cell` per gene within each Age×Arm×Group cell
- Pool all samples → build fixed edge list **E** (top-k or partial correlation)
- Prevents edges driven by between-cell mean shifts (Simpson's paradox)

### Phase 3: LIONESS Sample-Specific Networks
- Compute per-sample edge weights on fixed **E**
- `w_e(s) = N·w_e(all) - (N-1)·w_e(-s)`
- Fisher z-transform: `z = atanh(r)`

### Phase 4: Edge-Wise Regression
- Model each edge: `z_e ~ Age + Arm + Group + interactions + covariates`
- Empirical Bayes variance moderation (limma-style)
- Generate predicted contrast networks

### Phase 5: node2vec Embeddings & Alignment
- Multi-seed node2vec (10 seeds, d=128, p=0.25, q=4)
- Orthogonal Procrustes alignment using pre-registered anchors
- Consensus embedding across seeds

### Phase 6: Rewiring Quantification & Statistics
- Cosine distance rewiring: `Δ = 1 - cos(v₁, v₂)`
- **Silent shifters**: high Δ (top 10%, FDR<0.1) + low DE (|log₂FC|<0.3)
- Bootstrap CIs + permutation tests (n=2000 each)

### Phase 7: Leakage-Safe Cross-Validation
- 5-fold CV with fold-wise feature engineering
- Sample-level features: pathway strength, node strength, PCA
- RandomForest classification (Environment Group)

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

## 🚀 Quick Start

### Run Full Pipeline
```bash
python scripts/run_full_pipeline.py \
    --config config/hyperparameters.yaml \
    --data data/raw \
    --output results
```

### Configuration
Edit `config/hyperparameters.yaml` to adjust:
- node2vec parameters (dimensions, p, q)
- Network topology method (top-k, partial correlation)
- Statistical thresholds (FDR, percentiles)

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
│   ├── anchor_genes.yaml        # Pre-registered alignment anchors
│   ├── gene_sets.yaml           # DCT/NCC-WNK & pathway genes
│   └── hyperparameters.yaml     # Pipeline parameters
├── src/
│   ├── preprocessing/
│   │   ├── deconvolution.py     # Cell-type proportion estimation
│   │   ├── residualization.py   # Global confounder regression
│   │   ├── normalization.py     # VST transformation
│   │   └── qc.py                # Quality control
│   ├── networks/
│   │   ├── shared_topology.py   # Cell-standardized edge selection
│   │   ├── lioness.py           # Sample-specific networks
│   │   ├── edge_regression.py   # Edge-wise modeling
│   │   └── embeddings.py        # node2vec + Procrustes
│   ├── statistics/
│   │   ├── rewiring_metrics.py  # Cosine distance & silent shifters
│   │   ├── bootstrap.py         # Bootstrap CIs
│   │   └── permutation.py       # Permutation tests
│   └── validation/
│       ├── cross_validation.py  # Leakage-safe CV
│       └── sample_features.py   # Mouse-level feature extraction
├── scripts/
│   └── run_full_pipeline.py     # Master orchestration
├── METHODOLOGY.md               # Detailed methods
└── README.md
```

---

## 🔍 Key Results (Example)

### Silent Shifters (High Rewiring + Low DE)
Genes with substantial network context changes but minimal expression changes:
- WNK4 (NCC-WNK pathway kinase)
- KCNJ10 (potassium channel)
- Calcium handling regulators

### Validation Performance
Cross-validated classification (5-fold):
- **Accuracy**: ~75-85% (Environment Group prediction)
- **AUC**: ~0.80-0.90

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