This repository contains the complete computational pipeline for analyzing transcriptomic network rewiring in mouse kidney tissue following spaceflight exposure. Using NASA's OSD-771 dataset (RRRM-2 mission), we apply graph embedding techniques (node2vec) to identify genes whose network context shifts substantially despite minimal differential expression—termed "silent shifters".
Key Innovation
Traditional differential expression analysis misses regulatory changes that occur without large fold-changes. By embedding co-expression networks into vector spaces and quantifying geometric shifts, we capture topological rewiring that reveals hidden mechanisms of spaceflight-induced renal dysfunction.

Scientific Context
Problem: Astronauts experience kidney stone rates 2-7× higher than pre-flight within one year post-mission. Microgravity deactivates the NCC/WNK ion transport hub and remodels the distal convoluted tubule (DCT).
Gap: Transcriptomic mechanisms underlying these changes remain unclear, especially genes that rewire networks without changing expression levels.
Approach: We use node2vec embeddings to transform gene co-expression networks into 128-dimensional vectors, then quantify rewiring via cosine distance shifts (Δ) across age and flight conditions.

Dataset
NASA GeneLab OSD-771 (Rodent Research Reference Mission 2)

Species: Female C57BL/6NTac mice
Tissue: Whole kidney
Design: 2×2 factorial (Young/Old × Flight/Control)
Samples: 80 bulk RNA-seq (10 per condition)
Genes: ~15,000 after filtering

ConditionAgeFlight StatusSamplesLabelYoung Control16 weeksGround (GC+VIV)10YCYoung Flight16 weeksISS 35-day mission10YFOld Control34 weeksGround (GC+VIV)10OCOld Flight34 weeksISS 35-day mission10OF

Pipeline Overview
Phase 0: Preprocessing & QC

VST normalization (DESeq2)
TPM filtering (< 1 in >80% samples removed)
Cell-type deconvolution (DCT proportion estimation)
PCA/UMAP outlier detection

Phase 1: Network Construction

Spearman correlation matrices (controlling for cell-type covariates)
Fisher-Z transformation
Weighted, signed networks (top 1% edges retained)

Phase 2-3: Graph Embedding & Alignment

node2vec with 10 random seeds (p=0.25, q=4, d=128)
Orthogonal Procrustes alignment using housekeeping/ribosomal anchors
Consensus embeddings (average across seeds)

Phase 4-5: Rewiring Quantification

Compute Δ = cosine distance shift between aligned embeddings
Permutation testing (2000 iterations, stratified by condition)
Bootstrap confidence intervals (1000 iterations)
FDR correction (Benjamini-Hochberg)

Phase 6: Biological Interpretation

Gene Set Enrichment Analysis (GSEA) over ranked Δ
GO/KEGG pathway enrichment on top 5% Δ genes
Module detection (k-means clustering in embedding space)
Centrality dynamics (betweenness, eigenvector)

Phase 7: Silent Shifter Identification
Definition: Genes with high Δ (top 10%, FDR < 0.1) but low |log₂FC| (< 0.3, DE FDR > 0.2)
Phase 8: Predictive Validation

Random Forest classifiers trained on embeddings
Cross-condition transfer (train on control, test on flight)
Accuracy, AUROC, confusion matrices

Phase 9: Multi-Omics Integration

Triangulation with phosphoproteomics (if available)
Correlation of Δ scores with kinase activity changes


Installation
Option 1: Conda (Recommended)
bash# Clone repository
git clone https://github.com/yourname/rrrm2-kidney-node2vec.git
cd rrrm2-kidney-node2vec

# Create conda environment
conda env create -f environment.yml
conda activate rrrm2_kidney

# Install package in editable mode
pip install -e .
Option 2: Pip + Virtual Environment
bash# Create virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
pip install -e .

Quick Start
1. Download Data
bashpython scripts/download_genelab.py
2. Run Analysis Pipeline
bash# Launch Jupyter
jupyter lab

# Execute notebooks sequentially:
# 00_data_download.ipynb          → 01_qc_and_normalization.ipynb
# 02_cell_type_deconvolution.ipynb → 03_batch_effect_analysis.ipynb
# 04_network_construction.ipynb    → 05_node2vec_embedding.ipynb
# 06_rewiring_analysis.ipynb       → 07_pathway_enrichment.ipynb
# 08_silent_shifters.ipynb         → 09_predictive_validation.ipynb
3. Generate Figures
bash# Run visualization dashboard
jupyter notebook notebooks/11_visualization_dashboard.ipynb

Repository Structure
RRRM2_Kidney_Project/
├── data/                    # Raw, processed, and external datasets
├── notebooks/               # Analysis notebooks (00-11)
├── src/                     # Reusable Python modules
│   ├── preprocessing/       # Normalization, QC, deconvolution
│   ├── networks/            # Graph construction, embeddings
│   ├── statistics/          # Permutation, bootstrap, FDR
│   ├── enrichment/          # GSEA, GO/KEGG enrichment
│   └── visualization/       # Plotting utilities
├── results/                 # Generated figures, tables, embeddings
├── scripts/                 # Standalone batch scripts
├── tests/                   # Unit tests
├── config/                  # Hyperparameters, gene sets
└── README.md

Key Results (Expected)

Ranked Δ Atlas: Genome-wide rewiring scores for 4 comparisons:

Age effect in controls (OC vs YC)
Age effect in flight (OF vs YF)
Flight effect in young mice (YF vs YC)
Flight effect in old mice (OF vs OC)


Silent Shifter Candidates: ~50-200 genes with high network rewiring but low differential expression, enriched for:

DCT/NCC-WNK pathway components
Calcium signaling
Ion transport regulation


Spaceflight-Age Interaction: Quantitative evidence whether spaceflight exacerbates or mitigates age-related network changes
Therapeutic Targets: Pharmacologically tractable modules for countermeasure development (e.g., thiazides, amiloride, CaSR modulators)


Citation
If you use this code or methodology, please cite:
bibtex@article{yourname2025kidney,
  title={Age-Dependent Network Rewiring in Mouse Kidney After Spaceflight: A node2vec-Driven Framework},
  author={Your Name},
  journal={TBD},
  year={2025},
  doi={TBD}
}
Dataset Citation:
NASA GeneLab. (2024). Rodent Research Reference Mission 2 (RRRM-2) Kidney Transcriptome. OSD-771. https://genelab.nasa.gov/data/study?acc=OSD-771

Dependencies
Core Analysis

Python 3.10+
pandas, numpy, scipy
scikit-learn, networkx
pecanpy (node2vec implementation)

Bioinformatics

pydeseq2 (normalization)
scanpy (single-cell reference)
gseapy, goatools (enrichment)

Visualization

matplotlib, seaborn, plotly

See requirements.txt for complete list.





Contact
Ibrahim Shahid
The University of North Carolina at Greensboro
Email: imshahid@uncg.edu
GitHub: @ibrahimshahid1