# RRRM2 Kidney Transcriptome

Cross-cohort and cross-omic analysis of mouse kidney responses to spaceflight, centered on RRRM-2/OSD-771 and the matched OSD-462 RNA, protein, and phosphoprotein dataset.

## Current Status

The maintained analysis and manuscript are v11. The project has moved beyond the original DCT network-rewiring proposal and now focuses on what the public bulk datasets support directly:

- A recurrent matrix/endothelial-high and DCT transport-low RNA response across mouse kidney spaceflight cohorts.
- RNA-to-protein decoupling in OSD-462, with the strongest distal-nephron signal appearing at the regulatory phosphosite layer.
- Suppression of NCC/SPAK/WNK regulatory phosphorylation while total NCC abundance remains approximately flat.
- Enrichment of flight-suppressed phosphosites among DCT1-high parent genes from the GSE228367 reference, strongest in the top-decile prior subset and persistent after conservative parent-protein and bulk-compartment adjustment.
- Layer-specific propagation and observability audits showing that the main RNA/protein/phosphosite mismatch is not explained by simple protein-detectability bias.

These are whole-kidney observational results. The repository does not claim DCT1 cell of origin, causal mediation, fibrosis, or a newly discovered NCC dephosphorylation phenotype.

## Repository Map

- `src/v11/`: current DCT-subtype-prior, phosphoproteomic, recurrence, robustness, and layer-specificity analyses.
- `scripts/v11/`: executable v11 runners, including RNA-to-protein propagation and observability audits.
- `src/multiomics/`: OSD-462 harmonization, cell-type panels, regulator activity, and phenotype comparisons.
- `src/networks/` and `src/statistics/`: supporting and historical network, recurrence, permutation, and regression analyses.
- `config/`: curated gene sets, priors, metadata design, and analysis parameters.
- `tests/`: unit, regression, and fixture-locked headline-number checks.
- `docs/v11_execution_results.md`: core v11 execution report and interpretation boundaries.
- `docs/v11_layer_specificity_execution_summary_2026-06-07.md`: latest propagation and observability extensions.
- `latex_paper/manuscript_v11.tex`: maintained manuscript source.

Earlier manuscript versions are historical. New changes to those iterations are ignored; v11 is the maintained manuscript artifact.

## Setup

The environment definition targets Python 3.10 and includes the required R/Bioconductor packages:

```bash
conda env create -f environment.yml
conda activate rrrm2_kidney
```

The analysis expects downloaded OSDR and external reference files under `data/external/`. Large source data and generated run directories are intentionally excluded from Git.

## Run The Current Analysis

Run the v11 stack through the main entry point:

```bash
python -m src.run_all_phases --v11-only --run-id RUN_ID
```

Useful options include `--v11-site-scope`, `--v11-skip-spatial`, `--v11-skip-visium`, `--v11-skip-xenium`, and `--v11-skip-figures`.

Run the layer-specificity extensions directly:

```bash
python scripts/v11/05_rna_protein_propagation.py
python scripts/v11/06_observability_audit.py
```

## Verification

Run the full test suite:

```bash
python -m pytest -q
```

The v11 headline and layer-specificity values are locked in:

- `tests/fixtures/v11_headline_numbers.tsv`
- `tests/fixtures/v11_layer_specificity_numbers.tsv`

## Manuscript

Build the current manuscript from `latex_paper/`:

```bash
latexmk -pdf manuscript_v11.tex
```

The manuscript-facing claim language and allowed interpretation boundaries are summarized in `docs/v11_reviewer_ready_revision_packet.md`.
