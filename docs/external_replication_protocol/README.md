# External Replication Protocol

This protocol is committed before external cohort results are interpreted. The hypothesis registry contains concrete feature names, expected effect directions, discovery-effect signs, and fixed q-value thresholds; placeholder feature rows are rejected by the validation code.

Primary independent replication uses OSD-102 only for LAR-Young-like RRRM-2 findings and pathway signatures. Each cohort is analyzed independently; ComBat-seq is not required for this independent replication mode.

Secondary analysis uses OSD-513 only for sex-stratification and sex-robustness checks. It does not validate ISS-T claims, old-arm claims, classifier validation, or global cohort expansion.

OSD-163 and OSD-253 are context-mapping cohorts, not strict replication cohorts. They test a frozen biology-first remodeling panel covering PPAR/fatty-acid metabolism, cholesterol biosynthesis, ECM remodeling, EMT/fibrosis, tubular ion transport, TGF-beta/Wnt signaling, oxidative stress, and translation machinery. These rows use `discovery_effect = 0`, so passing the q threshold is reported as `context_detected`, not as a directional replication claim.

The preferred external pathway statistic is preranked GSEA over all mapped genes with sample-label permutation. Directional replication rows use the signed enrichment score; context-mapping rows use `abs(ES)` for two-sided pathway activity. The older mean pathway t-statistic is available only as a sensitivity analysis.

OSD-568 is excluded from validation claims in this repository.

Multi-study OSD-102 plus OSD-771 LAR-Young pooling is a separate analysis implemented in `src/validation/multi_study_pool.py`; that path requires pre-registered ComBat-seq and PCA checks before any pooled-network claim is allowed.
