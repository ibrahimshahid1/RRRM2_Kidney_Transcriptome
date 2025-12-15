# 🚀 Spaceflight Kidney Network Analysis (OSD-771)

[![Python](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![NASA GeneLab](https://img.shields.io/badge/Data-NASA%20GeneLab-red.svg)](https://genelab.nasa.gov/data/study?acc=OSD-771)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-brightgreen.svg)]()

> **Uncovering "Silent Shifters":** A computational pipeline for analyzing transcriptomic network rewiring in mouse kidney tissue following spaceflight exposure.

---

##  Table of Contents
- [Overview](#-overview)
- [Scientific Context](#-scientific-context)
- [Dataset](#-dataset)
- [Methodology & Pipeline](#-methodology--pipeline)
- [Installation](#-installation)
- [Quick Start](#-quick-start)
- [Repository Structure](#-repository-structure)
- [Results & Interpretation](#-results--interpretation)
- [Citation](#-citation)
- [Contact](#-contact)

---

##  Overview

This repository contains the complete computational pipeline for the analysis of the **NASA OSD-771 (RRRM-2)** dataset. We apply graph embedding techniques (**node2vec**) to identify genes whose network context shifts substantially despite minimal differential expression.

### Key Innovation
Traditional differential expression (DE) analysis misses regulatory changes that occur without large fold-changes. By embedding co-expression networks into vector spaces and quantifying geometric shifts, we capture **topological rewiring** that reveals hidden mechanisms of spaceflight-induced renal dysfunction.

---

##  Scientific Context

* **The Problem:** Astronauts experience kidney stone rates 2-7× higher than pre-flight within one year post-mission. Microgravity deactivates the NCC/WNK ion transport hub and remodels the distal convoluted tubule (DCT).
* **The Gap:** Transcriptomic mechanisms underlying these changes remain unclear, specifically regarding genes that rewire networks without changing expression levels.
* **Our Approach:** We use **node2vec** embeddings to transform gene co-expression networks into 128-dimensional vectors, quantifying rewiring via cosine distance shifts ($\Delta$) across age and flight conditions.

---

## Dataset

**Source:** NASA GeneLab OSD-771 (Rodent Research Reference Mission 2)  
**Species:** Female C57BL/6NTac mice  
**Tissue:** Whole kidney  
**Design:** 2×2 factorial (Young/Old × Flight/Control)

| Label | Condition | Age | Flight Status | Samples |
| :--- | :--- | :--- | :--- | :--- |
| **YC** | Young Control | 16 weeks | Ground (GC+VIV) | 10 |
| **YF** | Young Flight | 16 weeks | **ISS 35-day mission** | 10 |
| **OC** | Old Control | 34 weeks | Ground (GC+VIV) | 10 |
| **OF** | Old Flight | 34 weeks | **ISS 35-day mission** | 10 |

---

##  Methodology & Pipeline

The analysis is divided into 10 distinct phases, ranging from raw data processing to multi-omics integration.

#### Phase 0: Preprocessing & QC
* VST normalization (DESeq2)
* TPM filtering (< 1 in >80% samples removed)
* Cell-type deconvolution (DCT proportion estimation)
* PCA/UMAP outlier detection

#### Phase 1: Network Construction
* Spearman correlation matrices (controlling for cell-type covariates)
* Fisher-Z transformation & Weighted, signed networks (top 1% edges)

#### Phase 2-3: Graph Embedding & Alignment
* **node2vec:** 10 random seeds ($p=0.25, q=4, d=128$)
* **Alignment:** Orthogonal Procrustes using housekeeping/ribosomal anchors
* **Consensus:** Average embeddings across seeds

#### Phase 4-5: Rewiring Quantification
* Compute $\Delta$ (cosine distance shift)
* Permutation testing (2000 iterations) & Bootstrap CIs
* FDR correction (Benjamini-Hochberg)

#### Phase 6-7: Interpretation & "Silent Shifters"
* **Silent Shifter Definition:** Genes with high $\Delta$ (top 10%, FDR < 0.1) but low Log₂FC (< 0.3).
* GSEA, Module detection (k-means), and Centrality dynamics.

#### Phase 8-9: Validation & Integration
* Random Forest classifiers (transfer learning).
* Triangulation with phosphoproteomics (kinase activity correlation).

---

##  Installation

### Option 1: Conda (Recommended)
```bash
# Clone repository
git clone [https://github.com/yourname/rrrm2-kidney-node2vec.git](https://github.com/yourname/rrrm2-kidney-node2vec.git)
cd rrrm2-kidney-node2vec

# Create and activate environment
conda env create -f environment.yml
conda activate rrrm2_kidney

# Install package in editable mode
pip install -e .
# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -e .