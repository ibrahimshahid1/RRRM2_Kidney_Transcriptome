# RRRM2 Kidney Transcriptome

Rigor-remediated network rewiring analysis for RRRM-2 mouse kidney transcriptomics.

The repository now draws a hard line between candidates, calibrated inference, validation, and context mapping. Post-hoc focused permutation is removed; Procrustes anchors are configured and force-included; silent shifters require DE; Phase 8b is fold-safe; external replication and multi-study pooling are separate.

## Key Claims

Defensible after remediation:

- A biology-first remodeling panel spanning lipid/PPAR, cholesterol, ECM/fibrosis, tubular transport, TGF-beta/Wnt, oxidative stress, and translation programs across downloaded OSD cohorts.
- LAR-Young-like directional replication through OSD-102, if criteria pass.
- Sex-robustness checks through OSD-513 as secondary analysis.
- Context mapping through OSD-163 and OSD-253, reported as recurrence/detection rather than strict replication.
- ISS-T Young full-domain edge-sum hits as candidates only.

Still provisional:

- Broad per-gene discoveries.
- ISS-T or old-arm external validation.
- Classifier validation.
- Global cohort expansion.

## Important Modules

- `src/networks/shared_topology.py`: top-k 80 shared skeleton, DCT plus anchor force-inclusion.
- `src/networks/procrustes.py`: configured-anchor Procrustes alignment.
- `src/statistics/permutation_bootstrap.py`: edge-sum node-rewiring permutation/bootstrap, candidate-only BH, hierarchical FDR.
- `src/statistics/full_regression.py`: signed empirical gene aggregation.
- `src/statistics/silent_shifters.py`: DE-aware silent shifters and calibration.
- `src/statistics/full_pipeline_permutation.py`: appendix-tier cosine-distance full-pipeline permutation manifest/driver.
- `src/validation/enhanced_cv.py`: fold-safe Phase 8b with pooling modes and feature sets.
- `src/validation/external_replication.py`: independent OSD protocol guard for confirmatory and context-mapping cohorts.
- `src/validation/osd_external_validation.py`: builds protocol-bound OSD validation/context feature tables from downloaded VST matrices.
- `src/validation/multi_study_pool.py`: OSD-102 plus OSD-771 LAR-Young pooling guard with ComBat-seq/PCA requirements.
- `src/preprocessing/deconvolution_sensitivity.R`: MuSiC reference sensitivity.
- `docs/biology_first_framing.md`: current manuscript framing around lipid/ECM/tubular remodeling across cohorts.

## Quick Checks

```bash
venv/bin/python -m pytest tests
venv/bin/python -m py_compile src/statistics/permutation_bootstrap.py src/statistics/full_regression.py
```

Run the main remediated pipeline and append external validation:

```bash
venv/bin/python -m src.run_all_phases --hierarchical-fdr --external-validation
```

Run only downloaded supported OSD validation/context mapping:

```bash
venv/bin/python -m src.run_all_phases --external-validation-only
```

Run an explicit four-cohort external pass after downloading all folders:

```bash
venv/bin/python -m src.run_all_phases --external-validation-only --external-studies OSD-102,OSD-513,OSD-163,OSD-253
```

External validation now defaults to preranked GSEA over all mapped genes. The legacy pathway-mean statistic remains available for sensitivity checks:

```bash
venv/bin/python -m src.run_all_phases --external-validation-only --external-pathway-method mean_t
```

## External Data

External replication should store data under:

- `data/external/osdr/OSD-102/`
- `data/external/osdr/OSD-513/`
- `data/external/osdr/OSD-163/`
- `data/external/osdr/OSD-253/`

OSD-102/513 are directional validation cohorts. OSD-163/253 are context-mapping cohorts for lipid/ECM/tubular remodeling. Do not use excluded studies for validation claims. OSDR references: [NASA OSDR](https://www.nasa.gov/osdr/) and [OSDR Developer API](https://www.nasa.gov/reference/osdr-developer-api/).
