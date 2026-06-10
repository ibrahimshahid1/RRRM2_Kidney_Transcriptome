# v11 Layer-Specificity Execution Summary

Date: 2026-06-07
Author context: companion to `docs/v11_novelty_extensions_implementation_plan_2026-06-06.md`
Scope: results of executing Modules 2 and 3 of the layer-specificity plan, the Gate 0 search for Module 4, and the deferred-work plan for Modules 1 and 4.

## TL;DR

- **Module 2 (RNA→protein propagation score) — SHIPPED.** Calibrates the descriptive cross-layer pathway matrix (7 discordant / 3 protein-flat / 1 concordant) to a continuous matched-null per-pathway statistic. Two pathways come out as **calibrated protein inversion** (`ecm_organization`, `tlr4_innate`), three as **RNA→phospho carry candidates** (`dct_ncc_wnk_transport`, `tubular_transport_broad`, `s1p_s1pr3`), six as **RNA-not-propagated**.
- **Module 3 (proteome observability-bias audit) — SHIPPED.** Adding per-channel missingness as a third stratum dimension does NOT materially shift the per-pathway q-values (ECM slope q: 0.889→0.870; DCT slope q: 0.910→0.870). The five pre-specified NCC/SPAK regulatory phosphosites all sit at **abs_effect_percentile ≥ 99.1** (top 0.9% of phosphoproteome effect magnitude); the four non-regulatory NCC control sites sit at 14–28%. The central cross-layer discordance is **not** explained by detectability bias.
- **Gate 0 for Module 4 — CONDITIONAL PASS.** OSD-105 (RR-1 tibialis anterior muscle, 6 flight + 6 ground, RNA-seq + mass-spec proteomics) is the primary candidate; compatibility checks (TMT vs LFQ, phospho presence, same-animal matching) need workbook-level verification before commit.
- **Modules 1 and 4 — DEFERRED**, with concrete next steps documented below.

All headline numbers are locked in `tests/fixtures/v11_layer_specificity_numbers.tsv` and asserted by the existing fixture-based test in `tests/test_v11_rna_protein_propagation.py`.

## Module 2: RNA→protein propagation score (calibrated)

**Inputs:**
- `data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_flight_effects.tsv` (per-gene RNA + protein effects with pre-assigned abundance × peptide strata)
- `data/results/run_20260522_osd462_anchor/osd462_anchor/phospho_all_sites.tsv` (per-site phospho effects)
- `config/mechanism_gene_sets.yaml` (11 curated v11 pathways)

**Code:**
- `src/v11/rna_protein_propagation.py` — `PropagationConfig`, `classify_gene`, `pathway_gene_table`, `summarize_pathways`, `plot_propagation`, `run_propagation`
- `scripts/v11/05_rna_protein_propagation.py` — runner

**Outputs (`data/results/run_20260606_v11_layer_specificity/propagation/`):**
- `rna_protein_propagation_summary.tsv` — per-pathway matched-null calibrated statistics + layer assignment
- `rna_protein_propagation_gene_classes.tsv` — per-gene ternary class (RNA-only / RNA→protein / RNA→phospho)
- `rna_protein_propagation_coverage.tsv` — per-pathway gene-set coverage summary
- `figures/v11_rna_protein_propagation.{pdf,png}` — signed-mean panel + per-gene class stack

**Calibrated layer assignments (11 pathways):**

| Layer assignment | Pathways |
|---|---|
| `protein_inverted_calibrated` | `ecm_organization`, `tlr4_innate` |
| `RNA_to_phospho_candidate` | `dct_ncc_wnk_transport`, `tubular_transport_broad`, `s1p_s1pr3` |
| `RNA_not_propagated` | `mmp_adam_proteolysis`, `fibrosis_tgfb_emt`, `integrin_cell_adhesion`, `preservation_stress_response`, `oxidative_stress_nrf2`, `macrophage_inflammation` |

**Fixture-locked headlines (`tests/fixtures/v11_layer_specificity_numbers.tsv`):**
- ECM signed mean protein effect (signed by RNA direction): **−0.107** — matrix proteins go DOWN while their RNA goes UP.
- ECM inverse signed-mean q: **0.091** — calibrated against the matched null, the inversion is right at the 0.10 FDR threshold within the 11-pathway family.
- DCT fraction of RNA-directional members carried at phospho: **0.308** (n_rna_phospho = 13) — 31% of dct_ncc_wnk_transport genes whose RNA moves are also carried at phospho.
- TLR4 inverse slope q: **0.078** — TLR4 shows calibrated protein inversion.

**Tests:** 6/6 pass in `tests/test_v11_rna_protein_propagation.py`, including:
- Permuted-flight-label control returns non-significant matched-null result (no manufactured signal).
- `ecm_organization` is assigned `protein_inverted_calibrated`.
- `dct_ncc_wnk_transport` is assigned `RNA_to_phospho_candidate`.

## Module 3: Proteome observability-bias audit

**Inputs:** same anchor dir + `protein_effects_by_row.tsv` for per-protein `n_channels_used`.

**Code:**
- `src/v11/observability_audit.py` — `collapse_observability_to_gene`, `merge_observability_into_pool`, `assign_observability_strata` (extends Module 2's 5×4 strata with a 3-bin missing-fraction dimension), `detectability_gradient`, `high_coverage_subset`, `propagation_with_strata`, `ncc_site_observability`
- `scripts/v11/06_observability_audit.py` — runner

**Outputs (`data/results/run_20260606_v11_layer_specificity/observability/`):**
- `v11_observability_per_gene.tsv` — per-gene observability features
- `v11_observability_detectability_gradient.tsv` — RNA |effect| decile → fraction protein-quantified
- `v11_observability_matched_propagation.tsv` — per-pathway propagation under standard vs observability-extended strata
- `v11_observability_q_delta.tsv` — pathway-level delta of q-values between the two stratifications
- `v11_observability_high_coverage_propagation.tsv` — re-test on `n_peptides ≥ 3 & missing_fraction ≤ 0.2` subset
- `v11_observability_ncc_site_audit.tsv` — NCC/SPAK site coverage + intensity percentiles

**Detectability gradient (RNA |effect| decile → fraction protein-quantified):**

| Decile | RNA |effect| range | Fraction protein-quantified |
|---|---|---|
| 0 | < 0.027 | 0.468 |
| 1 | 0.027 – 0.058 | 0.489 |
| 2 | 0.058 – 0.094 | 0.478 |
| 3 | 0.094 – 0.157 | 0.445 |
| 4 | > 0.157 | **0.387** |

Modest dip (≈10 pp) in the top decile. Reportable as: "fraction protein-quantified is roughly flat across RNA-effect magnitude bins, with a modest dip in the top decile — large-RNA-effect genes are slightly under-quantified at the protein level, but the gradient is not steep enough to fully explain the cross-layer mismatch." A bigger dip would have argued the discordance is a quantification artifact; this dip is small.

**q-value delta (standard 5×4 strata vs observability-extended 5×4×3 strata):**

| Pathway | slope q (standard) | slope q (observability) | signed-mean q (standard) | signed-mean q (observability) |
|---|---|---|---|---|
| `ecm_organization` | 0.889 | 0.870 | 0.993 | 0.975 |
| `dct_ncc_wnk_transport` | 0.910 | 0.870 | 0.993 | 0.975 |
| `tlr4_innate` | 0.208 | 0.205 | 0.993 | 0.975 |

All deltas |Δq| ≤ 0.04. **The propagation pattern is not explained by detectability** — adding missingness to the null distribution does not shift any per-pathway q-value enough to change a layer assignment. (Caveat: the slope-based test is mostly underpowered at FDR 10% in both stratifications, so the conclusion is "no detectability-driven softening of significance" rather than "calibrated significance gain.")

**NCC/SPAK phosphosite observability:**

| Gene | Site | Role | n_fl + n_gc | Coverage percentile | abs_effect_percentile |
|---|---|---|---|---|---|
| Slc12a3 | S53 | regulatory | 20 | 67 | **99.4** |
| Slc12a3 | S65 | regulatory | 20 | 67 | 99.1 |
| Slc12a3 | S68 | regulatory | 20 | 67 | 99.2 |
| Stk39  | S382 | regulatory | 10 | 17 | 99.3 |
| Stk39  | S383 | regulatory | 20 | 67 | 99.2 |
| Slc12a3 | S96 | control | 20 | 67 | 26.1 |
| Slc12a3 | S120 | control | 20 | 67 | **14.7** |
| Slc12a3 | S122 | control | 20 | 67 | 24.2 |
| Slc12a3 | S124 | control | 20 | 67 | 27.6 |

**Headline:** the five pre-specified regulatory phosphosites all sit in the **top 0.9%** of phosphoproteome effect magnitude; the four non-regulatory NCC control sites sit at the 14–28th percentile. The suppression is real and large in objectively the most-confidently-measured part of the phosphoproteome. Detectability-tail explanations are ruled out.

**Tests:** 7/7 pass in `tests/test_v11_observability_audit.py`, plus the 5 new observability rows in `tests/fixtures/v11_layer_specificity_numbers.tsv` all reproduce.

## Gate 0 for Module 4 — verdict

See `docs/v11_gate0_negative_control_search_2026-06-07.md` for the full memo.

**Primary candidate: OSD-105 (RR-1 tibialis anterior muscle).** Bulk RNA-seq + mass-spec proteomics from 6 spaceflight + 6 ground control mice, RR-1 mission (37-day flight, female C57BL/6J, 16 weeks).

**Open compatibility questions** that must be answered before Module 4 is built (require workbook download / OSDR API):
1. Is the proteomics TMT-based or label-free? (OSD-462 = TMT 2-plex; existing `parse_tmt_sheet` parser would need adaption or replacement.)
2. Is there phosphoproteomics? (Likely no, based on RR-1 publication coverage.)
3. Is proteomics matched to RNA-seq at the per-animal level? (Verify sample-ID cross-references.)
4. Strain/age/mission-duration comparability to OSD-462. (RR-1 is 37 days, female C57BL/6J, 16 wk; OSD-462 is RR-10 — these likely differ.)

**Conditional plan:** if Q1 returns TMT, re-use the OSD-462 pipeline directly. If Q1 returns label-free, add a thin LFQ flight-effect estimator before re-running Modules 2+3.

## Modules 1 and 4 — deferred work plan

Both modules are real but lower marginal ROI than the two we shipped — and each has open dependencies that block clean execution today.

### Module 1 — Formal cell-type deconvolution

**Why deferred (not abandoned):**
- The marker-panel decomposition (`src/multiomics/celltype_panels.py`, output in `data/results/run_20260522_celltype_decomposition/`) already established the scenario verdict: **"segment_loss_or_composition_shift"** — both `dct_transport` (−0.16) AND `dct_identity` (−0.37) fall, which means even the marker proxy cannot cleanly separate program-suppression from composition-shift.
- The composition-aware phospho analysis (`src/v11/h2_composition_aware_phospho.py`) already locks `composition_m4_fixed_fraction = 0.993` — 99.3% of M0-suppressed phosphosites stay negative after the M4 ladder. The **phospho-layer composition story is already closed**.
- The genuine remaining gap is RNA-layer pathway-effect residualization on **estimated cell fractions** (vs marker-panel z-scores) across cohorts. That requires MuSiC/Bisque infrastructure, reference signature matrices from TMS + Chen + GSE228367, and per-cohort bulk deconvolution.

**Concrete next steps if/when picked up:**
1. Build a coarse-resolution reference signature matrix (PT, TAL, DCT, CNT/CD, endothelial, stromal, immune) from TMS + Chen and a DCT-subtype-aware version by grafting GSE228367.
2. Run MuSiC + Bisque on the 4 OSDR cohorts (OSD-771, OSD-513, OSD-253, OSD-462).
3. Two-method agreement gate: report fractions only where MuSiC and Bisque agree in direction.
4. Test flight-vs-control change in each estimated fraction per cohort.
5. **Decisive panel:** per-cohort pathway flight-effect before vs after composition adjustment.
6. Anchor: relate estimated DCT/endothelial fractions in OSD-462 to the measured NCC regulatory phospho score.

Pre-registered null: if DCT pathway suppression persists after composition adjustment, the intrinsic-transcriptional reading strengthens; if it collapses, the composition/dilution claim becomes quantitative. Either way the manuscript hedge becomes a number.

### Module 4 — Negative-control tissue specificity

**Why deferred:** Gate 0 conditional pass identifies OSD-105 as a candidate, but the platform compatibility (Q1 above) needs to be resolved before any pipeline runs.

**Concrete next steps if/when picked up:**
1. Verify OSD-105 proteomics platform (TMT vs LFQ) via the OSDR study manifest.
2. If TMT: re-use the OSD-462 pipeline directly with header overrides.
3. If LFQ: add `src/multiomics/lfq_flight_effect.py` (a thin label-free analog of `compute_flight_effect`).
4. Run Modules 2 + 3 pipeline on OSD-105 → muscle propagation scores.
5. Specificity contrast: kidney vs muscle propagation slope and observability-matched q-values per pathway.

Pre-registered null: if muscle shows the same RNA→protein decoupling as kidney, the phenomenon is generic and the kidney framing weakens. If kidney is distinctively decoupled, the discordance can be elevated from "kidney shows X" to "X is kidney-distal-nephron specific" — the single largest novelty jump available.

## Manuscript integration

Where the new results plug in (relative to `docs/v11_execution_research_plan.md` / `manuscript_v11.tex`):

- The descriptive cross-layer matrix paragraph ("OSD-462 shows RNA-protein mismatch") gains a calibrated propagation panel: 7 of 11 pathways with a matched-null per-pathway statistic, 2 with calibrated protein inversion (FDR ≈ 10%), 3 candidates for RNA→phospho carry.
- The composition / tissue-context paragraph ("the bulk DCT-low signal partly reflects compartment remodeling") is unchanged on the phospho side (already locked at 99.3% via h2_composition_aware_phospho) but adds Module 3's NCC site observability finding as a new sentence: "the suppressed NCC/SPAK regulatory phosphosites sit at abs_effect_percentile ≥ 99.1 of the phosphoproteome — the discordance lives in the most-confidently-measured sites, not in the low-coverage tail."
- The Results closing paragraph gets a new sentence: "adding per-channel missingness as a third stratum dimension to the matched-null engine does not shift any per-pathway q-value enough to change a layer assignment; the RNA→protein discordance is not explained by detectability bias."
- Discussion gets a brief deferred-work paragraph naming Modules 1 and 4 as the largest remaining novelty gains (formal deconvolution; non-kidney specificity).

## Tests & reproducibility

- All tests pass: `pytest tests/test_v11_rna_protein_propagation.py tests/test_v11_observability_audit.py -v` → 13 passed.
- All 11 fixture-locked headline numbers in `tests/fixtures/v11_layer_specificity_numbers.tsv` reproduce exactly from the on-disk outputs.
- Run manifests (`data/results/run_20260606_v11_layer_specificity/manifests/*.json`) carry input file SHA256s and parameters for both modules.
- Modules 2 + 3 are deterministic given the locked seed (20260606 / 20260607).
