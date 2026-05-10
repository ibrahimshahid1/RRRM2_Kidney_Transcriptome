# Methodology

This repository analyzes RRRM-2 kidney bulk RNA-seq with a graph-rewiring pipeline. The remediated methodology separates exploratory prioritization from inferential claims.

## Phase 2: Shared Topology

The executed topology default is top-k = 80 neighbors per gene, matching `config/hyperparameters.yaml` and `src/networks/shared_topology.py`.

Phase 2 now merges the DCT marker force-include list with configured Procrustes anchor Ensembl IDs resolved from `config/anchor_genes.yaml` through `data/processed/resources/id_map.tsv`. This ensures anchors are present in the 2,500-gene panel when they exist in the expression universe.

## Phase 3: Embeddings And Procrustes

Phase 3 computes node2vec embeddings and Procrustes-aligned cosine-distance rewiring. Procrustes anchors are loaded from `config/anchor_genes.yaml`; they are not selected from observed low-rewiring genes. The pipeline fails below 20 mapped anchors and warns below 50.

Ribosomal/translation genes are excluded from the primary anchor list because translation biology is a result of interest, not a neutral alignment reference.

## Phase 5: Silent Shifters

Strict silent shifters are now defined as high rewiring plus DE-null:

`high_rewiring == true`, `|log2FC| < 0.3`, and `DE FDR > 0.2`.

The generator requires contrast-matched gene-level DE unless explicitly run in exploratory mode. Outputs distinguish candidate rewired genes, DE-supported rewired genes, strict silent shifters, supported strict silent shifters, and observed-vs-null calibration.

## Phase 6: Inference

The fast permutation/bootstrap module tests **edge-sum node rewiring**:

`sum(abs(delta_edge))` over edges incident on each gene.

These p-values are not direct inference for Phase 3 cosine-distance rewiring. Appendix-tier full-pipeline cosine-distance permutation is implemented separately for pre-selected top hits in `src/statistics/full_pipeline_permutation.py`.

Post-hoc focused permutation has been removed. Restricted testing is valid only through:

- `--candidate-genes`: pre-registered candidate-only BH.
- `--hierarchical-fdr`: two-stage hierarchical FDR over pre-specified gene sets, with conservative BY columns for overlapping families.

`src/statistics/full_regression.py` uses a signed incident-edge t-statistic calibrated by within-stratum label permutation. Unsigned Stouffer aggregation is no longer primary inference.

## Phase 8 And 8b

Phase 8 remains the leakage-safe network-vs-expression classifier baseline.

Phase 8b now computes all transforms inside folds: network pooling, skeleton construction, sample-specific network weights, feature selection, scaling, PCA, and classifier fitting. Supported pooling modes are `age_arm_envgroup`, `arm`, `age`, and `global`. Supported network methods include LIONESS and SSN; CSN/scLink wrappers report dependency limitations when unavailable. Feature sets include node strength, pathway strength, sparse edges, edge PCA, and combined features.

Classifier outputs are validation only if they beat expression baselines under the fold-safe rerun. The audited run is negative.

## External Replication, Context Mapping, And Pooling

Independent replication is separate from multi-study pooling.

OSD-102 is primary for LAR-Young-like findings only and does not require ComBat-seq in independent replication mode. OSD-513 is secondary for sex robustness. Multi-study OSD-102 plus OSD-771 LAR-Young pooling is implemented separately and requires ComBat-seq plus PCA gates before pooled-network claims.

OSD-163 and OSD-253 extend the external layer from strict validation to biology-first context mapping. They are analyzed independently against a frozen remodeling panel: PPAR/fatty-acid metabolism, cholesterol biosynthesis, ECM remodeling, EMT/fibrosis, tubular ion transport, TGF-beta/Wnt signaling, oxidative stress, and translation machinery. Registry rows for these cohorts use `discovery_effect = 0`, so a passing q-value is reported as `context_detected`, not as a directional replication claim.

The main runner exposes this as `--external-validation`, which appends protocol-gated OSD validation/context mapping to a pipeline run. `--external-validation-only` runs just the downloaded external cohorts. By default `--external-studies auto` runs every supported study folder present under `data/external/osdr`; pass `--external-studies OSD-102,OSD-513,OSD-163,OSD-253` to require the four-cohort panel.

The external pathway statistic defaults to preranked GSEA over all mapped genes. For directional OSD-102/513 rows, the enrichment score must agree with the registered discovery direction and pass the registered q threshold. For OSD-163/253 context rows, the GSEA statistic is two-sided through `abs(ES)` and is reported as `context_detected`, not as `replicated`. The legacy mean pathway t-statistic remains available through `--external-pathway-method mean_t` for sensitivity analysis.

This validation builds independent expression-level pathway feature tables from the GeneLab VST matrices and then applies the committed hypothesis registry; it does not pool OSD-771 with external cohorts.

References for data access: [NASA OSDR](https://www.nasa.gov/osdr/) and [OSDR Developer API](https://www.nasa.gov/reference/osdr-developer-api/).
