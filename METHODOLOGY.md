# OSD-771 RRRM-2 Kidney Rewiring Analysis: Detailed Methodology

This document provides the complete methodological framework for the age-dependent network rewiring analysis.

## Overview

This pipeline analyzes the NASA GeneLab OSD-771 (Rodent Research Reference Mission 2) kidney transcriptome dataset using a novel combination of:
- Cell-type deconvolution
- Cell-standardized shared topology construction (prevents Simpson's paradox)
- LIONESS sample-specific network inference
- Edge-wise regression over full factorial design
- node2vec graph embeddings with Procrustes alignment
- "Silent shifter" detection (high rewiring, low differential expression)
- Leakage-safe cross-validation

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
- **Reference**: Murine kidney single-cell atlases (Tabula Muris Senis, Park et al. 2018)
- **Method**: Ridge regression or NNLS
- **Output**: Nephron segment proportions (DCT, PT, CD, Glom)
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

## Phase 2: Cell-Standardized Shared Topology

**Critical Innovation**: Prevent Simpson's paradox in edge selection.

### Within-Cell Standardization
For each gene g and each experimental cell c (Age × Arm × Group):
```
R*_ig = (R_ig - μ_cg) / σ_cg
```

Pool all R* across cells → use for topology selection.

### Edge Selection Methods
1. **Top-k neighbors** (default: k=50)
2. **Partial correlations** (Ledoit-Wolf shrinkage + precision matrix)
3. **Graphical lasso** (sparse precision estimation)

### Output
Fixed edge list **E** (50k-200k edges on 1k-3k gene panel).

---

## Phase 3: LIONESS Sample-Specific Networks

**LIONESS formula**:
```
w_e(s) = N · w_e(all) - (N-1) · w_e(-s)
```

Where:
- `w_e(all)`: pooled correlation for edge e
- `w_e(-s)`: leave-one-out correlation excluding sample s
- `N`: number of samples

### Fisher z-Transform
```
z = atanh(r)
```
Approximates normal distribution for regression modeling.

### Output
- **W_samp**: Sample-specific edge weights (n_samples × |E|)
- **Z_samp**: Fisher z-transformed weights

---

## Phase 4: Edge-Wise Regression

### Model (per edge e)
```
z_e ~ Age + Arm + Group + Age:Group + Arm:Group + Age:Arm + Age:Arm:Group + batch + SVs + cell_props
```

### Empirical Bayes Variance Moderation
- Borrow information across edges (limma-style)  
- Moderated variance: `σ²_mod = (d₀·s₀² + df·σ²) / (d₀ + df)`
- Improves stability with high-dimensional edge data

### Predicted Contrast Networks
For target condition profile P, predict each edge weight:
```
ẑ_e(P) = X(P) · β̂_e
```

Assemble weighted graph using E and predicted weights.

---

## Phase 5: node2vec Embeddings & Alignment

### Multi-Seed node2vec
- **Dimensions**: 128
- **Parameters**: p=0.25 (return), q=4.0 (in-out)
- **Random seeds**: 10 (for stability assessment)

### Procrustes Alignment
Using pre-registered anchor genes:
1. Align target embedding to reference via orthogonal rotation R
2. Minimize `||Y_anchor - X_anchor @ R||_F`
3. Apply R to all nodes

### Consensus Embedding
Average across 10 aligned embeddings.

---

## Phase 6: Rewiring Metrics & Statistics

### Cosine Distance Rewiring
```
Δ(P₁, P₂; g) = 1 - cos(v_P₁(g), v_P₂(g))
```

### Silent Shifter Criteria
- **High rewiring**: Δ in top 10% AND FDR < 0.1
- **Low DE**: |log₂FC| < 0.3 AND DE FDR > 0.2

### Statistical Inference
- **Bootstrap** (n=2000): CIs for Δ, topology fixed
- **Permutation** (n=2000): stratified within Age/Arm
- **FDR correction**: Benjamini-Hochberg + Westfall-Young for top hits

---

## Phase 7: Leakage-Safe Cross-Validation

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
4. NASA GeneLab OSD-771 dataset.

---

**Document Version**: 1.0  
**Last Updated**: 2026-01-07
