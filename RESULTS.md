# Results Status After Rigor Remediation

This document reflects the remediated code path. It intentionally narrows claims where the available evidence is weak.

## Current Audited Run

Latest audited run used here: `data/results/run_20260408_191759_2500g`.

Phase 6 full-domain edge-sum node-rewiring permutation results from that run:

| Contrast | Full-domain BH q < 0.05 | Minimum q |
|---|---:|---:|
| ISS-T Young FLT minus GC | 24 | 0.0208 |
| ISS-T Old FLT minus GC | 0 | 0.0833 |
| LAR Young FLT minus GC | 0 | 0.0833 |
| LAR Old FLT minus GC | 0 | 0.1250 |

These 24 ISS-T Young genes are **edge-sum node-rewiring candidates**, not confirmed per-gene discoveries. Phase 6 p-values test summed incident LIONESS edge changes. They are not direct p-values for the Phase 3 node2vec/Procrustes cosine-distance statistic.

The old focused-permutation narrative has been removed. Restricting BH to the observed top decile does not provide strict FDR control. Valid restricted modes are now limited to pre-registered `--candidate-genes` or `--hierarchical-fdr` over pre-specified gene sets.

## Gene Aggregation

`src/statistics/full_regression.py` no longer uses unsigned Stouffer aggregation as primary inference. Gene-level regression evidence is now a signed incident-edge t-statistic calibrated by within-Age by Arm label permutation. Legacy Stouffer output is available only as an explicitly named diagnostic.

## Silent Shifters

The previous strict silent-shifter files contained 250 genes per contrast because DE was not attached and all genes were treated as DE-null. In the audited run those files had blank DE columns; the 55-68 values were support counts inside those rewiring-only files, not true silent-shifter counts.

The remediated generator requires contrast-matched gene-level DE and writes separate outputs for candidate rewired genes, DE-supported rewired genes, strict DE-null silent shifters, supported strict silent shifters, and observed-vs-null calibration.

A smoke-test rerun for `ISS_T_YNG_FLT_minus_GC` with the available DE table produced 250 high-rewiring candidates, 208 strict DE-null silent shifters, 36 DE-supported high-rewiring genes, and 3 supported strict silent shifters. Full contrast-level counts require rerunning Phase 5 on the remediated code.

## Phase 8

The classifier result is negative in the audited run:

| Model | Accuracy | AUC |
|---|---:|---:|
| Random Forest | 0.400 | 0.525 |
| Logistic Regression | 0.375 | 0.359 |

Enhanced Phase 8b network advantage was negative in every stratum: -0.8, -0.8, -0.4, -0.3, and -0.1 pooled. The old Phase 8b implementation computed LIONESS/features/PCA before LOO splitting, so it was not fully leakage-safe. The remediated implementation computes network pools, skeletons, LIONESS/SSN weights, feature selection, scaling, PCA, and classifiers inside folds. Until a remediated rerun changes the numbers, classifier validation remains unsupported.

## External Validation Scope

Independent external replication is protocol-gated in `src/validation/external_replication.py` and `docs/external_replication_protocol/`.

OSD-102 is the primary partner only for LAR-Young-like findings and pathways. OSD-513 is secondary and limited to sex-stratification or sex-robustness checks. Neither validates ISS-T claims, old-arm claims, classifier validation, or global cohort expansion.

Multi-study OSD-102 plus OSD-771 LAR-Young pooling is a separate mode in `src/validation/multi_study_pool.py`; it requires ComBat-seq and PCA gates before any pooled-network claim is allowed.

OSDR source references: [NASA OSDR](https://www.nasa.gov/osdr/) and [OSDR Developer API](https://www.nasa.gov/reference/osdr-developer-api/).

## Defensible Claims

Defensible now:

- Pathway-level replication, if a remediated external run validates it.
- LAR-Young external validation using OSD-102, if the pre-registered checks pass.
- Sex-robustness checks using OSD-513, clearly marked as secondary.
- ISS-T Young full-domain edge-sum hits as prioritized candidates only.

Provisional or unsupported:

- Broad per-gene discoveries.
- ISS-T external validation.
- Old-arm external validation.
- Classifier validation.
- Global cohort expansion.
