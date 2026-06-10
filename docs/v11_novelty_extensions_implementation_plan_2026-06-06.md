# v11 Novelty-Extension Implementation Plan

**Date:** 2026-06-06
**Author context:** companion to `manuscript_v11.tex`, `v11_novel_directions_2026_05_29.tex`, and `docs/v11_reaudit_2026-06-03.md`
**Status:** design only — no computation run yet

---

## 0. Guiding principle

The reaudit established the binding constraint: the phospho anchor re-analyzes the same OSD-462/RR-10 data Siew 2024 already published, so the only genuinely **Siew-independent** contribution v11 owns is the **cross-layer RNA→protein discordance and its composition interpretation**. Every module below is chosen to *harden, quantify, or generalize that one claim* rather than to chase a bigger result the data do not contain.

The four builds, in priority order:

1. **Formal cell-type deconvolution** — converts the central "tissue-context / composition" claim from marker-panel proxies to quantitative cell fractions (tests composition vs. intrinsic regulation directly).
2. **RNA→protein propagation score** — turns the descriptive cross-layer direction matrix into a calibrated per-pathway statistic with a null.
3. **Proteome observability-bias audit** — pre-empts the first reviewer objection ("the mismatch is just proteome detectability").
4. **Negative-control tissue specificity** — *data-gated*; if matched non-kidney multi-omic spaceflight data exist, shows the discordance is kidney/distal-nephron–specific rather than a generic spaceflight-omics artifact. Biggest novelty jump, lowest feasibility certainty.

Plus an assessment of the **DCT regulatory-network reconstruction (#10)** the user asked about (§6).

Each module is specified with a **pre-registered null interpretation** so a negative outcome is as reportable as a positive one — consistent with the project's reviewer-2 posture.

---

## 1. Shared infrastructure & conventions

- **New run directory:** `data/results/run_20260606_v11_layer_specificity/` with one subfolder per module: `deconvolution/`, `propagation/`, `observability/`, `negative_control_tissue/`, plus `figures/`, `logs/`, `manifests/`.
- **New code location:** `src/v11/` (extend the existing module set) + thin runner scripts in `scripts/v11/`.
- **Reuse, don't reinvent:** these modules already exist and should be imported, not duplicated —
  - `src/multiomics/osd462_anchor.py` — OSD-462 RNA/protein/phospho flight-effect harmonization (Layer 0–4).
  - `src/multiomics/celltype_panels.py` — `KIDNEY_PANELS`, `per_sample_panel_scores`, `panel_flight_effect`.
  - `scripts/celltype/run_celltype_decomposition.py` — marker-panel decomposition scaffold (the abandoned MuSiC attempt lives here).
  - `scripts/osd462/01_protein_concordance.py` — the abundance- & peptide-count-matched random-gene-set null machinery (10,000 draws). **This is the null engine reused by modules 2 and 3.**
  - `src/common.py` (`id_map_lookup`) — Ensembl/symbol identifier bridging.
- **Identifier discipline:** RNA = ENSMUSG; protein/phospho = gene symbol bridged through the project ID map. All cross-layer joins go through `id_map_lookup`.
- **Testing:** each module ships a `tests/test_v11_<module>.py` with at least (a) a fixture-locked headline number and (b) a null-control assertion (e.g. permuted labels return non-significant). Follows the pattern of `tests/test_v11_h2_enrichment.py` / `tests/fixtures/v11_headline_numbers.tsv`.
- **Run manifest:** every module writes input file SHA + row counts to `manifests/` (pattern of `v11_baseline_input_manifest.json`).

---

## 2. Module 1 — Formal cell-type deconvolution

### 2.1 Estimand / question
Is the recurrent matrix-high / DCT-low **bulk RNA** signature explained by **shifts in cell-type composition** (interstitial/endothelial expansion diluting DCT signal) versus **within-cell-state transcriptional regulation**? The manuscript currently asserts the composition reading from marker-panel z-scores and an endothelial↔NCC-phospho correlation; this replaces the proxy with estimated cell fractions and a composition-residualized re-test.

### 2.2 Inputs (all on disk)
- **Bulk targets:** OSD-771 (RRRM-2), OSD-513, OSD-253, OSD-462 VST/count matrices under `data/external/osdr/OSD-*/`.
- **Single-cell/-nucleus references:**
  - `data/external/single_cell_atlases/tms_kidney_female_ALLDATASETS_counts_innerGenes.h5ad` (Tabula Muris Senis kidney).
  - `data/external/single_cell_atlases/chen_atlas/GSE150338_Seurat_IntegrateData.RData` (broad kidney atlas).
  - `data/external/dct_reference/GSE228367` (DCT1/DCT2 resolution — already the manuscript's DCT prior).

### 2.3 Method
1. **Build reference signature matrices** at two resolutions: (a) coarse nephron compartments (PT, TAL, DCT, CNT/CD, endothelial, stromal/fibroblast, immune); (b) DCT-subtype-aware (DCT1, DCT2/CNT) by grafting GSE228367 onto the coarse atlas. Emit cell-type marker/signature matrices with detection-aware gene filtering.
2. **Run ≥2 deconvolution methods for robustness** (do not rely on one):
   - **MuSiC** (R; reuses the empty `deconv_music_out/` target and `install_r_packages.R`) — tree-guided, handles cross-subject variance.
   - **BisqueRNA** or **CIBERSORTx-style NNLS** as an orthogonal estimator.
   - Report fractions only where the two methods agree in direction; treat disagreement as an uncertainty flag.
3. **Composition flight effect:** per cohort, test flight-vs-control change in each estimated fraction (the quantitative version of "endothelial/stromal up, DCT down").
4. **Residualized pathway effects (the key panel):** regress each curated pathway score on estimated cell fractions; recompute the flight effect on residuals. The decisive comparison is pathway flight effect **before vs after** composition adjustment, per cohort.
5. **Anchor cross-check:** in OSD-462, relate estimated DCT/endothelial fractions to the measured NCC regulatory phospho score (does the bulk endothelial↔phospho correlation survive when endothelial is an estimated *fraction* rather than a marker z-score?).

### 2.4 Deliverable
- **Figure:** multi-cohort panel — estimated fractions by flight (left); pathway flight effect before/after composition adjustment (right).
- **Table:** per-cohort fraction shifts + residualized pathway effects with CIs.
- **Manuscript slot:** upgrades §"The bulk DCT-low signal partly reflects compartment remodeling" from proxy argument to quantitative result.

### 2.5 Risks / pre-registered null
- Deconvolution uncertainty is large in injured/remodeled tissue; references are female/aging-biased (TMS) and strain-mismatched. **Mitigation:** two-method agreement gate; report confidence bands, not point fractions.
- **Null interpretation:** if DCT pathway suppression *persists* after composition adjustment, that argues for an **intrinsic** transcriptional component (strengthens a DCT-cell-autonomous reading) — reportable and interesting. If it *collapses*, that quantitatively confirms the composition/dilution claim. Either way the manuscript's hedge becomes a number.

### 2.6 Effort: **Medium-High** (R+Python, reference construction is the bulk of the work; scaffold already exists).

---

## 3. Module 2 — RNA→protein propagation score

### 3.1 Estimand / question
For each pathway/gene family, **how much of the RNA flight effect propagates to protein abundance**, calibrated against a null? Replaces the descriptive `v11_cross_layer_heatmap` direction calls (↑/↓/≈) with a continuous, testable propagation statistic.

### 3.2 Inputs (on disk)
- OSD-462 matched RNA flight effects, gene-level protein flight effects, parent-gene phospho effects — all already produced by `src/multiomics/osd462_anchor.py` and emitted to the v11 run's `cross_layer/`.
- The 11 curated pathway panels (frozen member lists in the v11 run artifacts).

### 3.3 Method
1. **Define a per-gene propagation coefficient:** for genes with both RNA and protein effects, propagation = signed agreement scaled by magnitude (e.g. protein effect regressed on RNA effect within pathway; slope β is the propagation rate). Classify each gene into RNA-only, RNA→protein (concordant), or RNA→phospho (carried at phospho not abundance) — the ternary/Sankey structure from the table.
2. **Replication calibration (the "null model"):** the slope/agreement is meaningless without a baseline because genome-wide RNA–protein correlation is generically weak (your Spearman −0.034). Calibrate each pathway's propagation against the **abundance- and peptide-count-matched random-gene-set null** already implemented in `scripts/osd462/01_protein_concordance.py` (10,000 draws). A pathway "propagates" only if its propagation score exceeds the matched-null distribution.
3. **Layer assignment per pathway:** emit, for each of the 11 panels, the fraction of member genes in each propagation class with bootstrap CIs, oriented so the matrix/ECM inversion and the DCT→phospho carry are visible.

### 3.4 Deliverable
- **Figure:** ternary or Sankey of RNA-only / RNA→protein / RNA→phospho classes across pathways (directly the figure proposed in the table).
- **Table:** per-pathway propagation score vs matched-null p, with the matrix/ECM inversion and DCT/NCC-WNK→phospho rows as the anchors.
- **Manuscript slot:** replaces/augments the descriptive cross-layer matrix in §"OSD-462 shows RNA-protein mismatch"; makes the headline a calibrated quantity.

### 3.5 Risks / pre-registered null
- Proteome coverage is incomplete and dynamic range differs by layer → propagation undefined for unobserved genes. **Mitigation:** restrict to co-observed genes and report coverage per pathway (this also feeds Module 3).
- Phospho layer collapses to parent-gene mean, so RNA→phospho is parent-gene resolution, not site resolution — state explicitly.
- **Null interpretation:** the *expected* result is near-zero propagation for most pathways (that's the point). The reportable structure is *which* pathways deviate: matrix/ECM (inverted) and DCT/NCC-WNK (carried at phospho). If some pathway unexpectedly *does* propagate cleanly RNA→protein, that's a genuine new finding.

### 3.6 Effort: **Low-Medium** (mostly reuse of existing anchor + null engine; the new work is the scoring function and figure).

---

## 4. Module 3 — Proteome observability-bias audit

### 4.1 Estimand / question
Could the RNA→protein non-propagation be an artifact of **proteome detectability** — peptide coverage, missingness, abundance-dependent quantifiability — rather than real biological decoupling? This reviewer-proofs the central claim.

### 4.2 Inputs (on disk)
- OSD-462 proteomics/phospho **raw workbooks** with peptide counts and per-channel missingness (`GLDS-462_proteomics...WorkUp...xlsx`, `GLDS-462_phosphoproteomics...Pho_WorkUp_JM.xlsx`) under `data/external/osdr/OSD-462/` and `data/external/phosphoproteomics/`.
- The QC missingness summaries already emitted by the v11 baseline (`external_qc/`, TMT channel-centering QC).

### 4.3 Method
1. **Characterize observability:** for every RNA-detected gene, record whether it is protein-quantified, peptide count, mean intensity, and per-channel missing fraction. Quantify the detectability gradient (which transcripts never reach the proteome and why).
2. **Observability-matched null benchmark (core):** the observed RNA–protein mismatch statistic (e.g. the genome-wide Spearman and the per-pathway concordance) is recomputed on **observability-matched resamples** — draw control gene sets matched on peptide count + abundance + missingness, and ask whether the *matched* null already produces a mismatch as large as observed. If observed ≈ matched-null, the mismatch is partly detectability-driven; if observed exceeds it, the decoupling is real beyond observability.
3. **Sensitivity:** restrict the cross-layer and propagation analyses (Module 2) to high-coverage, low-missingness genes and confirm the mismatch persists. Re-test the matrix/ECM inversion specifically in the high-confidence subset.
4. **Phospho parallel:** confirm the suppressed NCC/SPAK sites are not in low-observability tails (they aren't — they're high-intensity quantified sites — but show it).

### 4.4 Deliverable
- **Figure:** QC panel — detectability gradient + observed mismatch vs observability-matched null distribution.
- **Table:** mismatch statistics in full vs high-confidence subsets; matched-null comparison.
- **Manuscript slot:** a Supplementary subsection + one sentence in Results closing the detectability objection; cite from the main mismatch paragraph.

### 4.5 Risks / pre-registered null
- Matching on observed covariates can't rule out unobserved detectability confounds — state as a bound, not proof.
- **Null interpretation:** if the matched null reproduces much of the mismatch, the central claim must be **softened** to "decoupling beyond what detectability explains is modest" — this is the honest, important outcome and exactly what the audit posture demands. If observed clearly exceeds the matched null, the claim is materially strengthened.

### 4.6 Effort: **Low-Medium** (reuses the matched-null engine; new work is the observability feature table).

---

## 5. Module 4 — Negative-control tissue specificity (DATA-GATED)

### 5.1 Estimand / question
Is RNA→protein non-propagation **specific to the kidney distal nephron**, or a **generic property of spaceflight multi-omics**? If specific, the discordance becomes a kidney biology finding rather than a methods artifact — the single largest novelty jump available.

### 5.2 Gate 0 — data availability (DO THIS FIRST; blocks the rest)
There is **no non-kidney matched multi-omic spaceflight data on disk** (the download script only covers kidney OSDR cohorts + OSD-656 urine). Before any analysis:
1. Search NASA OSDR / GeneLab for datasets with **matched RNA-seq + TMT proteomics (ideally + phospho)** from the **same animals** in a **non-kidney** tissue. Leading candidates:
   - **Other RR-10 / OSD-462 tissues** (RR-10 profiled multiple organs — if other tissues have matched RNA+protein, this is the cleanest same-mission control).
   - Other GeneLab multi-omic missions (e.g. liver, muscle, heart) with paired transcriptome+proteome.
2. Verify per candidate: same-animal matching, TMT proteomics present, Ensembl/symbol ID compatibility, flight-vs-control design. Tabulate verdicts (pattern of the `v11_novel_directions` Appendix A compatibility matrix).
3. **If ≥1 compatible non-kidney matched cohort exists → proceed.** If none → **stop and document the negative**: "matched non-kidney spaceflight multi-omics do not currently exist publicly," which itself is a reportable limitation and a wet-lab/data-generation recommendation. Do not force an underpowered or mismatched comparison.

### 5.3 Method (conditional on Gate 0 passing)
1. Run the **identical** propagation + observability pipeline (Modules 2–3) on each non-kidney matched cohort — same code, same nulls.
2. **Specificity contrast:** compare propagation scores and the mismatch magnitude kidney vs non-kidney. The claim is supported if the kidney distal-transport / matrix pathways show *stronger* RNA-protein decoupling (or the specific matrix/ECM inversion) than the non-kidney tissue baseline.
3. Frame strictly as **tissue/pathway specificity of the decoupling**, controlling for the generic weak RNA-protein correlation that exists everywhere.

### 5.4 Deliverable
- **Figure/Table:** comparative propagation/mismatch statistics across tissues; specificity contrast.
- **Manuscript slot:** new Discussion paragraph elevating the discordance from "kidney shows X" to "X is kidney-distal-nephron specific."

### 5.5 Risks / pre-registered null
- Tissue confounding (different proteome depth, different missions) likely; platform/depth differences can masquerade as specificity — **the observability audit (Module 3) must be run per tissue** to separate true specificity from coverage differences.
- **Null interpretation:** if non-kidney tissues show the *same* RNA-protein decoupling, the phenomenon is **generic** and the kidney framing weakens — important to know and report. If kidney is distinctively decoupled, large novelty gain.

### 5.6 Effort: **Medium** if data exists (pipeline is reused); the gate is the real cost.

---

## 6. On the DCT regulatory-network reconstruction (#10)

**Recommendation: do not build as a novelty module; optionally produce as a single summary figure.**

Rationale, grounded in the repo:
- The project already invested heavily in **data-driven network inference** — `src/networks/` (LIONESS, node2vec, WGCNA, edge-regression, Procrustes, contrast vectors, silent-shifter scoring) and the v2–v9 "network era" framing. The data falsified it twice: RRRM-2 aging vectors were bootstrap-unstable (3/4 arm settings), and network-prioritized genes did **not** enrich among OSD-462 protein/phospho-changing genes (pre-registered **H4 null**). v11 keeps it only as a *negative control* ("the method didn't manufacture the result").
- The #10 proposal is a *different object* — a **curated** mechanistic regulatory network (WNK-SPAK/NCC, KLHL3/CUL3, ion-transport, aldosterone nodes) with edges colored by RNA/protein/phospho evidence, built on OmniPath/kinase-substrate priors. That specific artifact was **not** built.
- **But** it would largely re-package evidence already in the paper: the WNK-SPAK-NCC axis is Fig 1 (TikZ), the per-layer evidence is the cross-layer heatmap, and KLHL3/CUL3 coverage is already inventoried (and found empty: Klhl3 = 0 phosphosites, Cul3 = 1). Your own audit warns such a network "may look too schematic unless quantitatively scored." It adds presentation, not inference, and risks reviving the framing the data already rejected.
- **If wanted at all:** scope it as one curated figure that overlays the existing OSD-462 layer effects onto a literature WNK-SPAK-NCC-KLHL3/CUL3 diagram (edges = measured RNA/protein/phospho direction), explicitly labeled as a summary schematic, not new analysis. Low effort, low risk, low novelty — a Discussion aid at most.

---

## 7. Sequencing, dependencies, deliverables

| Order | Module | Depends on | New code | Effort | Novelty ROI |
|-------|--------|-----------|----------|--------|-------------|
| 1 | RNA→protein propagation score | osd462_anchor, matched-null engine | `src/v11/rna_protein_propagation.py` | Low-Med | High |
| 2 | Observability-bias audit | raw peptide tables, matched-null engine | `src/v11/observability_audit.py` | Low-Med | High (defensive) |
| 3 | Formal deconvolution | sc references, celltype_panels | `src/v11/deconvolution.py` + R | Med-High | High |
| 4 | Negative-control tissue (Gate 0 first) | Modules 1–3 pipeline; external data | `src/v11/negative_control_tissue.py` | Med (if data) | Highest if feasible |
| opt | DCT network figure | existing cross-layer outputs | figure script only | Low | Low |

**Recommended execution order rationale:** Modules 1 and 2 are the cheapest (pure reuse of the anchor + null engine), directly harden the central claim, and produce the figures the propagation argument needs — do them first. Module 3 (deconvolution) is the highest-value single addition but the most build-heavy; do it next. Module 4's gate runs in parallel from the start (it's just a data search) but the analysis waits on the Modules 2–3 pipeline being stable, since it reruns the same code per tissue.

**Cross-module reuse:** the abundance/peptide/missingness **matched-null engine** (`scripts/osd462/01_protein_concordance.py`) is the shared dependency of Modules 2, 3-anchor-check, and 4 — refactor it once into `src/v11/matched_null.py` before starting so all three call the same tested implementation.

---

## 8. Per-module test & verification plan

For each module:
1. **Fixture lock:** headline numbers written to `tests/fixtures/v11_layer_specificity_numbers.tsv`; test asserts reproduction.
2. **Null control:** permuted flight labels must return non-significant (guards against pipeline-manufactured signal — the same discipline that caught the network framing).
3. **Coverage report:** every cross-layer join logs n genes in / n co-observed / n dropped, so attrition is auditable.
4. **Figure-vs-text check:** each emitted figure caption's numbers are pulled from the same TSV the test locks (no hand-typed values).
5. **Independent re-read:** after the run, a fresh audit pass (the project's recurring practice) confirms no figure over-claims relative to its table.

---

## 9. What this does and does not buy

**Does:** makes the cross-layer discordance — the one thing not owned by Siew — quantitative (propagation score), reviewer-proof (observability audit), mechanistically interpretable (deconvolution: composition vs intrinsic), and potentially *kidney-specific* (negative-control tissue). These convert the paper's soft spots into pre-registered tests whose null outcomes are reportable.

**Does not:** produce a result larger than what the data contain. Per the reaudit, a categorically bigger claim needs DCT-enriched / spatial / inflight phosphoproteomics or matched physiology — wet-lab, not public-data computation. These modules maximize the defensible novelty available *in silico*; they do not substitute for that experiment.
