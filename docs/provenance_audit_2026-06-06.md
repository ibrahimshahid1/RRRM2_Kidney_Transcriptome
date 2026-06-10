# Data-provenance audit — hardcoded values masquerading as computed results

**Date:** 2026-06-06
**Scope:** `src/`, `scripts/`, `config/`, `data/external/`, `latex_paper/`
**Trigger:** review of `src/v11/human_concordance.py` flagged hand-typed marker
buckets, scoring rules, and figure-evidence directions presented as analysis output.

The goal of this pass is to separate three very different things that all look
alike in code: (1) **asserted results** — numbers/directions that feed a claim but
were typed in, not computed; (2) **buried curation** — defensible analyst choices
(predictions, gene sets) with no recorded provenance; and (3) **legitimate config**
— file paths, design constants, and already-sourced gene sets. Only (1) and (2) are
problems. Blanket "delete all hardcoded values" would destroy legitimate config and
preregistered predictions, so each item is classified individually.

---

## Tier 1 — Asserted results that feed a claim (highest priority)

These are numbers or directions that enter a reported result or figure but are not
computed from data.

### 1.1 `src/v11/human_concordance.py` — AQP2 figure direction — **FIXED & VERIFIED**
- **Was:** `FIGURE_EVIDENCE_ROWS` hardcoded AQP2 `observed_direction="up"`,
  `scored=True`, flowing into `concordance_verdict()` as a Bernoulli trial. The code
  itself admitted "no numeric digitization."
- **Fix:** Digitized Twins SM Fig. S4A from the on-disk PDF
  (`aau8650_garrett-bakelman_sm.pdf`, p.26, rendered ~432 dpi). Flight-twin (TW)
  in-flight vs pre-flight medians: **AQP2 up** (~0.43→0.62), **AGT down**
  (~30000→10000), **RENR down** (~0.40→0.27). Directions now live in
  `config/human_concordance_prereg.yaml` with provenance pointing at
  `data/external/human_spaceflight/figS4A_provenance/` (rendered panel crops +
  `figS4A_digitization.json` carrying the source PDF SHA-256).
- **Effect on results:** none — AQP2's direction was already `up`, so the verdict is
  byte-for-byte preserved (3/3 axes concordant, exact two-sided p=0.25, 4/4 scored
  analytes). RENR resolved from `requires_digitization`→`down` (still unscored
  context). All 7 `tests/test_human_concordance.py` pass.

### 1.2 `scripts/integration/plot_sixrec_integration.py:125` — hardcoded figure caption — **OPEN (recommend fix next)**
- **Issue:** the panel title hardcodes `"endo vs NCC phospho rho=-0.76, p<.001"`.
  The script reads `phospho_axis_summary.tsv` but **never computes** an endothelial-
  vs-NCC-phospho Spearman correlation, and no results file documents −0.76. This is
  a statistic baked into a figure label with no traceable computation.
- **Recommended fix:** either (a) compute the correlation from the cell-type
  deconvolution vs NCC-phospho data in-script and format the caption from the
  computed value (`f"rho={rho:.2f}, p={p:.1e}"`), or (b) if that correlation is not
  actually part of the analysis, remove the numeric claim from the caption. Do not
  leave a typed-in coefficient in a published figure.

---

## Tier 2 — Buried predictions / curation lacking provenance — **FIXED**

Defensible analyst choices that were correct in substance but undocumented and
editable-in-code. All moved to version-controlled, cited config.

### 2.1 `human_concordance.py` `TABLE_ANALYTE_SPECS` (prediction + scoring table)
- **Note:** this was *half* the critique — the observed directions are genuinely
  computed from Table S8; only the `expected_direction` and `scored/excluded` flags
  were analyst predictions. The real defect was that a de-facto **preregistration**
  sat in Python where it could be silently edited.
- **Fix:** moved to `config/human_concordance_prereg.yaml` with `prediction_basis`
  and `source` per analyte, under a documented prereg discipline. Code loads it via
  `load_prereg()`. Behavior unchanged.

### 2.2 `human_concordance.py` OSD-656 marker buckets + `CCL/CXCL/IL/MMP/TIMP` prefix heuristic
- **Note:** lower-stakes than it looks — these only label/sort an **optional context
  table** that never enters the sign test. But they were hand-typed biology with no
  source, plus a `startswith()` heuristic.
- **Fix:** moved to `config/human_urine_marker_panel.yaml`. The distal/RAAS category
  is now **composed from the already-sourced gene sets** in
  `config/mechanism_gene_sets.yaml` (`dct_ncc_wnk_transport`, `macrophage_inflammation`,
  `mmp_adam_proteolysis`) plus explicitly cited extra members. The prefix match is
  retained only as an **explicit, toggleable, citation-noted** nomenclature fallback
  (`family_prefix_fallback.enabled`), not buried logic. Context counts unchanged
  (203 markers, 92 recovery, 0 direct).

---

## Tier 3 — Gene-set duplication / single-source-of-truth violations — **OPEN (itemized)**

Lower stakes: these are duplication, not fabrication. The repo already has a correct
externalization pattern (`config/gene_sets.yaml`, `config/mechanism_gene_sets.yaml`,
`config/anchor_genes.yaml`, resolved through the Ensembl-backed
`src/enrichment/gene_set_loader.py`). The following files re-declare gene sets as
Python literals instead of loading from config. Recommended fix: load from config so
there is one sourced definition per program.

| File | Constant(s) | Note |
|---|---|---|
| `src/v11/core_analysis.py` | `ANCHOR_GENES` | DCT anchor; overlaps `config/gene_sets.yaml:dct_ncc_wnk` |
| `scripts/v11/02_run_v11_core_analysis.py` | `ANCHOR_GENES`, `TRANSPORT_TARGETS` | duplicates core_analysis + config |
| `src/v11/spatial_dct_transport_check.py` | `DCT_TRANSPORT_GENES` | duplicates config dct set |
| `scripts/osd462/02_phospho_axis.py` | `AXIS_GENES`, `GATE_GENES`, `NCC_REGULATORY_SITES` | axis genes duplicate config; sites should cite a phosphosite source |
| `scripts/osd462/01_protein_concordance.py` | `TARGET_SETS` | overlaps mechanism_gene_sets |
| `scripts/regulator_activity/run_phenotype_anchor.py` | `DCT_GENES` | duplicates config dct set |
| `src/networks/lar_reversal.py` | `CORE_CLOCK_GENES`, `PER_CRY_GENES`, `BMAL_CLOCK_GENES`, `DCT_WNK_GENES`, `ALDOSTERONE_ENAC_GENES`, `S1P_AXIS_GENES`, `RBM3_PRESERVATION_GENES` | clock genes are textbook but uncited; others duplicate config |
| `src/visualization/publication_figures.py` | `RETINOL_GENES` | hand-typed; should cite a GO/KEGG retinol set or move to config |
| `src/v11/spatial_reference_projection.py`, `scripts/v11/04_*` | `NICHE_MARKERS` | spatial niche markers; need a cited source |
| `src/v11/perturbation_gse228367_lowk.py`, `h2_occupancy_normalized_phospho.py`, `dct_continuous_gradient.py`, `aldosterone_axis.py` | `TARGET_GENES` etc. | confirm each cites or loads from config |

These should each `from src.enrichment.gene_set_loader import load_gene_sets` (or read
the relevant `config/*.yaml`) rather than hold a private copy. Textbook sets (core
clock) still want a one-line citation comment.

---

## Tier 4 — Non-functional stubs / scaffolding — **OPEN (low risk, label or remove)**

- `scripts/run_full_pipeline.py` — every phase is a `logger.info("[Placeholder] ...")`
  with no computation. Not imported anywhere and not referenced by the manuscript.
  The real analyses live in `scripts/v11/`, `scripts/osd462/`, etc. **Risk:** someone
  mistakes it for the production pipeline. Recommend renaming to
  `run_full_pipeline_skeleton.py` / moving to `docs/`, or deleting.
- `src/validation/multi_study_pool.py:67` — "Guarded placeholder for ComBat-seq."
  Confirm it raises or is clearly inert rather than silently returning unbatched data.

---

## Cleared — NOT problems (do not "fix")

- `src/validation/osd_external_validation.py` `STUDY_SPECS` — legitimate dataset
  config (VST filenames, ISA zips, FLT/GC regex). Correct as-is.
- `src/visualization/publication_figures.py` `CONTRAST_DE_NAMES` — DE output filenames.
- `src/markers/discover_markers.py`, `discover_dct.py` `p_value = 2.0 * sf(...)` —
  computed two-sided p-values, not literals.
- `src/validation/external_replication.py:182` — **good practice**: raises if the
  hypothesis registry contains placeholder features (a guard *against* fabrication).
- `config/*.yaml` gene sets — already sourced and resolver-backed; this is the model
  the Tier-3 offenders should follow.

---

## Verification performed for the fixed items

- `python -m py_compile src/v11/human_concordance.py` → OK.
- Baseline vs refactor run diff: all scored metrics identical (status, 3 axes, p=0.25,
  4 scored analytes); only RENR's unscored direction changed (the intended fix).
- `pytest tests/test_human_concordance.py` → 7 passed.
- AQP2 direction independently re-read from the rendered Fig. S4A panel (saved in
  `data/external/human_spaceflight/figS4A_provenance/`).

## Suggested order for the remaining work

1. **Tier 1.2** (figure caption `rho=-0.76`) — highest integrity priority.
2. **Tier 3** gene-set dedup — mechanical; do per-file with a test run each.
3. **Tier 4** stub labeling.
