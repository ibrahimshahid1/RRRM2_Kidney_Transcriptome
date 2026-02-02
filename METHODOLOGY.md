# OSD-771 RRRM-2 Kidney Rewiring Analysis: Detailed Methodology

This document provides the complete methodological framework for the age-dependent network rewiring analysis.

## Overview

This pipeline analyzes the NASA GeneLab OSD-771 (Rodent Research Reference Mission 2) kidney transcriptome dataset using a novel combination of:
- Cell-type deconvolution (**MuSiC** with Tabula Muris Senis kidney reference)
- Cell-standardized shared topology construction (prevents Simpson's paradox)
- LIONESS sample-specific network inference
- Edge-wise regression over full factorial design
- **PecanPy** node2vec graph embeddings with Procrustes alignment
- "Silent shifter" detection (high rewiring, low differential expression)
- Permutation and bootstrap uncertainty estimation
- Biological grounding via gene set enrichment

> **Note:** In the current implementation, Phase 4 (edge regression) is integrated into Phase 2. The pipeline runs phases {0, 1, 2, 3, 5, 6, 7}.

---

## Full Factorial Design (n=80)

**Primary design for all causal inference:**
- **Age**: Young (16 weeks) vs Old (34 weeks)
- **Arm**: ISS-T (terminal on ISS) vs LAR (live animal return)
- **Environment Group**: FLT (flight), HGC (hardware ground control), VIV (vivarium control), BSL (basal pre-flight)

Each cell (Age × Arm × Group) contains 5 mice per age level.

**Key contrasts:**
1. ISS-T young flight effect: `(Young, ISS-T, FLT) - (Young, ISS-T, HGC)`
2. ISS-T old flight effect: `(Old, ISS-T, FLT) - (Old, ISS-T, HGC)`
3. Age-dependent rewiring: `|Δ_old - Δ_young|`
4. LAR comparisons (persistence/recovery)

---

## Phase 0: Preprocessing & Deconvolution

### Expression Normalization
- **Input**: Raw counts from NASA GeneLab
- **Method**: DESeq2 VST (variance-stabilizing transformation)
- **Filtering**: Remove genes with CPM < 1 in > 80% of samples

### Cell-Type Deconvolution
- **Reference**: Tabula Muris Senis kidney atlas
- **Method**: **MuSiC** (Multi-Subject Single Cell Deconvolution)
- **Output**: Nephron segment proportions (DCT, PT, CD, Glom, Immune, Stroma)
- **Transformation**: Centered log-ratio (CLR) to avoid collinearity

### Quality Control
- PCA/UMAP outlier detection
- Variance partitioning (technical vs biological sources)
- Pre-registered exclusion criteria (> 3 SD in robust PC space)

---

## Phase 1: Global Residualization

**Goal**: Remove technical and compositional confounding while preserving biological signal.

### Model
```
Y_g ~ batch + lane + shipment + SVs + CLR(DCT) + CLR(PT) + CLR(CD) + CLR(Glom)
```

### Two Output Matrices
1. **R_tech**: Residuals with biology preserved (used for networks)
2. **R_all**: Residuals with Age/Arm/Group also removed (QC only)

### Surrogate Variable Analysis
- Capture hidden confounding via PCA on residuals
- Retain SVs explaining up to 50% variance (max 10 SVs)

---

## Phase 2: Network Construction (Shared Topology + LIONESS + Edge Regression)

**Critical Innovation**: Prevent Simpson's paradox in edge selection.

### Within-Cell Standardization
For each gene g and each experimental cell c (Age × Arm × Group):
```
R*_ig = (R_ig - μ_cg) / σ_cg
```

Pool all R* across cells → use for topology selection.

### Edge Selection Methods
- **Method**: Ledoit-Wolf shrinkage → precision matrix → partial correlations
- **Top-k neighbors**: k=50-80 per gene (union of all top-k gives skeleton **E**)
- **Output**: Fixed edge list **E** (50k-200k edges on 1k-3k gene panel)

### LIONESS Sample-Specific Networks
**LIONESS formula**:
```
w_e(s) = N · w_e(all) - (N-1) · w_e(-s)
```

Where:
- `w_e(all)`: pooled correlation for edge e
- `w_e(-s)`: leave-one-out correlation excluding sample s
- `N`: number of samples

**Fisher z-Transform**:
```
z = atanh(r)
```
Approximates normal distribution for regression modeling.

**Output**:
- **W_samp**: Sample-specific edge weights (n_samples × |E|)
- **Z_samp**: Fisher z-transformed weights

### Edge-Wise Regression
**Model (per edge e)**:
```
z_e ~ Age + Arm + Group + Age:Group + Arm:Group + Age:Arm + Age:Arm:Group + batch + SVs + cell_props
```

**Empirical Bayes Variance Moderation**:
- Borrow information across edges (limma-style)  
- Moderated variance: `σ²_mod = (d₀·s₀² + df·σ²) / (d₀ + df)`
- Improves stability with high-dimensional edge data

**Predicted Contrast Networks**:
For target condition profile P, predict each edge weight:
```
ẑ_e(P) = X(P) · β̂_e
```
Assemble weighted graph using E and predicted weights.

---

## Phase 3: node2vec Embeddings & Procrustes Alignment

### Multi-Seed PecanPy node2vec
- **Implementation**: PecanPy (fast Python node2vec)
- **Dimensions**: 128
- **Parameters**: p=0.25 (return), q=4.0 (in-out)
- **Random seeds**: 10 (for stability assessment)
- **Walk length**: 80, **Walks per node**: 200

### Procrustes Alignment
Using pre-registered anchor genes:
1. Align target embedding to reference via orthogonal rotation R
2. Minimize `||Y_anchor - X_anchor @ R||_F`
3. Apply R to all nodes

### Consensus Embedding
Average across 10 aligned embeddings.

### Cosine Distance Rewiring
```
Δ(P₁, P₂; g) = 1 - cos(v_P₁(g), v_P₂(g))
```

---

## Phase 5: Silent Shifters & Interaction Metrics

### Silent Shifter Criteria
- **High rewiring**: Δ in top 10% AND FDR < 0.1
- **Low DE**: |log₂FC| < 0.3 AND DE FDR > 0.2

### Interaction Persistence
- Compare rewiring across contrasts (ISS-T vs LAR, Young vs Old)
- Identify genes with consistent vs variable network changes

---

## Phase 6: Uncertainty Estimation & Full Regression

### Statistical Inference
- **Permutation** (n=2000): stratified within Age/Arm
- **Bootstrap** (n=2000): CIs for Δ, topology fixed
- **FDR correction**: Benjamini-Hochberg + Westfall-Young for top hits

### Full Factorial Regression (n=80)
- Edge-wise regression on all 80 samples
- Used for final contrast network prediction

### Gene-Level Differential Expression
- limma for per-gene DE analysis
- Used for silent shifter classification

---

## Phase 7: Biological Grounding

### Gene Set Enrichment
- GO biological process, molecular function, cellular component
- KEGG pathway analysis
- Reactome pathway analysis

### Pre-Registered Gene Sets
- **DCT/NCC-WNK pathway**: WNK1, WNK4, STK39, SLC12A3, KCNJ10, etc.
- **Positive controls**: ECM remodeling, oxidative stress, calcium handling, lipid metabolism

### Module Analysis
- k-means clustering on embeddings
- GO/KEGG/Reactome enrichment for high-Δ modules

---

## Validation: Leakage-Safe Cross-Validation

**Note**: This is implemented in `src/validation/` but not run as part of the main pipeline by default.

### Fold-Wise Operations (5-fold stratified CV)
**For each fold:**
1. Fit residualization on training only
2. Build skeleton E from training (with cell-standardization)
3. Select anchors from training
4. Compute LIONESS features per sample
5. Train classifier on training features
6. Evaluate on test features

### Sample-Level Features
- Pathway edge strength (mean weight in NCC-WNK subnetwork)
- Node strength for candidate genes
- Shifter-centered connectivity
- PCA features (10 components on stable edges)

### Classifier
RandomForest (100 trees) predicting Environment Group.

---

## Biological Interpretation

### Pre-Registered Gene Sets
- **DCT/NCC-WNK pathway**: WNK1, WNK4, STK39, SLC12A3, KCNJ10, etc.
- **Positive controls**: ECM remodeling, oxidative stress, calcium handling, lipid metabolism

### Module Analysis
- k-means clustering on embeddings
- GO/KEGG/Reactome enrichment for high-Δ modules

---

## References

1. Kuijjer et al. (2019) LIONESS. *iScience*.
2. Grover & Leskovec (2016) node2vec. *KDD*.
3. Smyth (2004) limma empirical Bayes. *Stat Appl Genet Mol Biol*.
4. Wang et al. (2019) MuSiC. *Nature Communications*.
5. NASA GeneLab OSD-771 dataset.
6. Tabula Muris Consortium (2020) Tabula Muris Senis. *Nature*.

---

**Document Version**: 2.0  
**Last Updated**: 2026-02-01
