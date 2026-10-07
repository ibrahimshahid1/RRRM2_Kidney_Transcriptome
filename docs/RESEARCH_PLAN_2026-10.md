# Research plan: cross-mission recurrence of renal transcriptomic responses to spaceflight

**Version:** 1.0 (2026-10-07). Supersedes the "next steps" sections of earlier decision documents.
**Status:** plan of record. Each analysis below is registered separately (§5) before it is run.
**Consolidates:**
- `docs/NEXT_DIRECTIONS_2026-10-07.md`
- `docs/strategy/2026-10-07/` (reader notes, draft, adversarial review, GSE295428 feasibility)
- the owner's external assessment (2026-10-07)
- the G1/G2 results (`docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md`)

> **Do not give this document to the G3 gate-writer.** It contains observed values. The gate-writer
> receives only the blind brief in Appendix A (§5).

---

## 1. Summary

**Central question.** *Which reported renal transcriptomic responses to spaceflight recur across
independent mouse missions, and which of them survive calibration for cohort-level (habitat,
handling, batch) variation?*

**What the existing data establish:**
- Across five terminal missions (83 animals), the four clinically anchored renal axes do not recur. The
  glomerular-barrier axis even moves opposite to the declared "identity loss".
- A low-heterogeneity structural/compartment band shifts upward. An externally defined 157-gene
  podocyte-associated program sits at its top (g = 0.689), and it is the only one of 49 sets whose
  interval excludes zero.
- That direction is robust to every preregistered perturbation (G1).
- It is not separable from the shared structural shift (K0).
- Persistence after live return cannot be identified (G2).

**What remains open.**

1. **Attribution.** Do these shifts survive when cohort-level variance is calibrated, and when
   alternative control groups are used?
2. **Specificity.** Is the podocyte program's position special relative to matched gene sets under that
   calibration?
3. **Literature recurrence.** Which specific published claims recur?
4. **Tissue specificity.** Is the structural shift kidney-specific, or systemic or handling-related
   (the same mice's heart)?

Every outcome yields a publishable answer to the central question (§7).

**Parked.** Recovery and podocyte-specific biology are not identifiable with available public data.
§8 lists exactly which data would reopen each.

## 2. What is established, refuted and retired

### 2.1 Established (frozen, preregistered or locked; reproduced 2026-10-06)

| Result | Numbers | Source |
|---|---|---|
| 4 clinical axes, 5 missions | barrier −0.716 (−1.361, −0.072) max-T 0.0018 (opposite to declared); tubular injury 0.570, FWER 0.026; fibrosis 0.311, 0.785; distal −0.154, 0.986 | clinical-axes revision |
| 49-set compartment family | podocyte HS 0.689 (0.042, 1.336), max-T 0.019, I² 0, mission effects −0.12 / 1.17 / 0.68 / 0.54 / 1.17; structural control 0.549 (−0.095, 1.193); ranks 2–5 within 0.02; 32/49 sets positive | compartment-context |
| K4 strict matching | 10,000 balanced same-tier 142-gene panels: none as extreme (p = 0.0001); target-minus-matched g 0.54 (−0.17, 1.25), blocked p 0.015 | `PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md` |
| K0 specificity | P0 0.198 → P1 0.112 (−0.031, 0.255) score units; S1 0.031; Freedman–Lane P1 max-T 0.042 (baseline map) / 0.056 (runner-up) → **not separable** | K0 stage |
| G1 sensitivities | ROBUST_DIRECTION (items 1–7); DEFINITION_ROBUST (8/8 atlas LOSO); no gene-map flips; leave-top-5 lower CI −0.044 | G1/G2 results doc |
| G2 recovery | podocyte INCONCLUSIVE (ISS-T g 1.17, LAR 0.14; interaction −0.28 (−0.62, 0.06)); structural interaction −0.35 (−0.63, −0.07), FL p 0.009, max-T 0.024 (label INCONCLUSIVE after Holm) | G1/G2 results doc |

### 2.2 Found after G1/G2 (post hoc, unregistered, so these are motivation, not results)

| Finding | Numbers | Consequence |
|---|---|---|
| Control choice | podocyte flight−vivarium meta 0.28 (−0.45, 1.01), k = 4 C57BL/6J; structural flight−vivarium 0.55 | Attribution to spaceflight is exposed → H1 |
| Cohort-level variance | control-vs-control contrasts imply about 0.36 excess between-group variance; calibrated interval about (−0.31, 1.68) | → H1 |
| OSD-771 batch structure | ERCC mix and dissection day are confounded with group. Flight−vivarium is the only mix-matched flight contrast (podocyte −0.17) | → H1, A2 |
| G2 "rebound" (r = −0.65) | reproduced by basal−ground (−0.65) and vivarium−ground (−0.48) without flight animals | **Refuted**: a cohort artifact |
| G2 exploratory interaction hits | ISS-T flight carries on-orbit handling signatures (immediate-early genes g 2.14, medulla panel 2.71, blood 1.73) | **Retired** as biology |
| Technical correction | GC/length/expression residualization keeps podocyte rank 1 (g 0.65–0.71) | → H4 |
| RIN | adjustment moves LAR podocyte from 0.14 to 0.17 | negligible; formalize in A4 |
| Sixth kidney cohort | OSD-457 (JAXA MHU-3) unscored; descriptive podocyte g 1.16 (wild-type 1.36) | → H7 |
| Eligibility rule | `cpm_eligible_genes` uses flight/ground membership: "arm-aware", not label-blind | → A4 |
| Prior art (Siew 2024, full text) | reports qualitative nephron structural remodeling (focal adhesion, tight junction, actin); PODXL/NID1 among the most frequently *down*-regulated by vote count | structural remodeling is not novel; podocyte direction conflicts → H5 |
| GSE295428 (RRRM-1 single-cell) | the "Kidney, LAR" files are spleen and duplicate the "Spleen, LAR" libraries (99.6–99.9% shared barcodes); kidney exists only for habitat controls (14–23 strict podocytes per animal) | **closed** (§8) |

### 2.3 Retired, do not reopen

- network rewiring
- contrast-vector framing
- LAR reversal and any "rebound"
- Grey60 rescue
- OSD-462 phospho mining (stop rule)
- DCT2/CNT and ASDN headlines
- the OSD-462 data note, the matched-library paper and the merged methods paper
- the existing Casaletto email draft
- podocyte as headline or title
- gene-list mechanism stories (glycocalyx, sialylation) without preregistered sub-modules

## 3. Hypotheses

Each hypothesis states its rationale, prediction, test and data status. **Seen** means the readers have
already computed a version post hoc, so the analysis is a *registered sensitivity analysis (post hoc in
origin)*. **Unseen** means no one has computed it, so a directional prediction can be registered
honestly.

| ID | Hypothesis | Test (analysis) | Data status |
|---|---|---|---|
| **H1** | The recurrent podocyte and structural shifts are spaceflight-associated, i.e. robust to cohort-level variance and control choice | A1, A2 | **Seen** |
| **H2** | Under that calibration, the podocyte program's shift exceeds matched gene sets (specificity of its position) | A3 (extends K4) | Partly seen (K4 done without calibration) |
| **H3** | A podocyte residual is separable from structural drift (K0 level 1) | Frozen: K0 says not separable. Nothing in this plan can upgrade it | Closed in-sample |
| **H4** | The direction is not an artifact of eligibility, normalization, scoring or RNA quality | A4 | Seen (partly) |
| **H5** | Specific published renal spaceflight claims recur across missions | A5 (ledger) | **Unseen** for most claims |
| **H6** | The structural shift is kidney-specific; the cohort artifact and handling signatures are animal-level, so they recur in another organ of the same mice | A6 (heart) | **Unseen** |
| **H7** | A live-return cohort sampled about 2 days after landing shows the same sign | A7 (OSD-457) | Seen descriptively; supplement only |
| **H8** | Podocyte and structural proteins shift in the same direction as the transcripts (same animals, second assay) | A8 (OSD-163, then OSD-102) | **Unseen** |
| **H9** | The bulk podocyte signal reflects (a) more podocytes, (b) more transcription per podocyte, or (c) structural/compositional change | A9 (morphometry, external) | Unseen; long lead |

**Recovery is not a hypothesis in this plan** (§8).

## 4. Datasets

### 4.1 In hand (on disk; inputs verified 2026-10-05/06)

| Dataset | Design | Role |
|---|---|---|
| OSD-102 (RR-1) | female C57BL/6J, 37 d, terminal; 6 flight / 6 ground | primary |
| OSD-163 (RR-3) | female BALB/c, 39–42 d, terminal; 6 / 6 | primary; mapping-rate flag |
| OSD-253 (RR-7) | C57BL/6J and C3H/HeJ, 25 and 75 d; 10 / 9 for B6; original and rerun controls | primary; duration and strain context |
| OSD-462 (RR-10) | B6129 female, 28–29 d; total RNA, mRNA and UPX libraries; 10 / 10 | primary; library-preparation sensitivity |
| OSD-771 (RRRM-2) | female B6NTac; ISS-T 53–56 d and LAR 32 d + 24 d; flight, ground, basal and vivarium groups, 5 per age block | primary (ISS-T); G2; control-choice |
| OSD-513 (RR-23) | male C57BL/6J, 37 d, euthanized about 1 d after landing; 9 / 9 | descriptive moderator |
| Mouse Kidney Atlas (Zenodo 17395591) | 141,401 cells, 8 source studies | frozen marker tiers |
| Tabula Muris Senis kidney (CELLxGENE) | 10x + Smart-seq2 | deferred second atlas |
| KEGG_2019_Mouse | Enrichr | structural control definition |
| Gene maps | reconstructed baseline and raw Ensembl 116 | `config/gene_map_reconstruction.yaml` |

### 4.2 To acquire (files verified on OSDR 2026-10-07)

| Dataset | Files | Size | Design | Analysis |
|---|---|---|---|---|
| **OSD-580** (RRRM-2 heart, same mice: 78/78 animal IDs shared with kidney) | `GLDS-573_rna_seq_RSEM_Unnormalized_Counts.csv`, `_Normalized_Counts.csv`, runsheet, SampleTable (**no GeneLab VST**) | 13–31 MB | ISS-T whole heart: flight 10 (on-orbit), ground 10, basal 10, vivarium 8. LAR: right ventricle, 10 per group | A6 |
| OSD-561 / 562 / 613 (RRRM-2 cerebellum / hippocampus / cerebrum) | GLDS-556 VST, counts and qc for 561; others similar | 5–12 MB each | flight vs ground, pooled controls; no basal or vivarium; metadata irregularities | A6 (secondary) |
| **OSD-457** (JAXA MHU-3, multi-tissue) | `GLDS-457_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv`, v2 runsheet, qc_metrics | 31 MB | male C57BL/6J, 31 d, live return, about 2 d post-landing after behavioural tests; wild-type 3 / 3, Nrf2-KO 3 / 3 (kidney 12 of 192 samples) | A7 |
| **OSD-163 proteomics** | `GLDS-163_proteomics_RR3_KDN_processed.zip` | 21 MB | same animals as OSD-163 RNA | A8 (first) |
| OSD-102 proteomics | `GLDS-102_proteomics_TMT_KDNL.processed.tar.gz` | 3.0 GB | same animals as OSD-102 RNA, TMT | A8 (only if OSD-163 passes) |
| Ensembl gene GC% and length | BioMart, release matching the gene map | small | — | A3, A4 |
| Region marker reference | Ransick 2019 (GSE129798: cortex, outer and inner medulla) | moderate | external, label-free | A4 sampling covariates |
| Prior-claim sources | Siew 2024 Supplementary Data (pathways, most-frequent genes); Suzuki 2022 (Kidney Int, MHU-3); other claims enumerated in A5 | small | — | A5 |

### 4.3 Watchlist and closed paths

| Item | Status | What would change it |
|---|---|---|
| GSE295428 (RRRM-1 single-cell) | **closed**: "Kidney, LAR" is a duplicate spleen upload | a corrected deposit of the live-return kidney libraries; report the duplicate to GEO and the authors |
| RRRM-1 bulk kidney mRNA | not found on OSDR 1–1150 or GEO (only miRNA, OSD-913) | an OSDR reply to the data request |
| Bion-M2 (2025) kidney | not public | a deposit |
| Siew GCR kidney (OSD-706..712) | not public | a deposit |
| GSE235042 (hindlimb-unloading kidney) | dropped (n = 3 per group, FPKM only, sex unknown) | — |
| Second atlas (GSE146912, GSE111107 + GSE129798) | deferred | a reviewer demand, or "podocyte" returning to the title |

## 5. Governance: making the new tests honest

1. **Blind gate-writing (A0).**
   - A fresh agent, or a person, who has not read §2.2, `reader_leads.md` or any scratchpad output
     writes the decision gates for A1–A4. They receive only Appendix A: designs, estimands and
     decision principles, with no observed values.
   - The gates are committed as `config/clinical_axes_g3_registration.yaml`.
   - The code is hash-locked through the rrrm2 CLI before execution.
2. **Directional predictions before download** for the unseen analyses (A5, A6, A8). They are
   committed before any data are fetched.
3. **Labels.**
   - A1–A4 are "registered sensitivity analyses (post hoc in origin)".
   - A5, A6 and A8 are "preregistered secondary analyses".
   - Nothing is "confirmatory".
4. **Provenance ledger in the paper:** targeted axes → 3/4 no-go → 49-set scan → podocyte → K1–K4 →
   G1/K0/G2 → reader post-hoc checks → G3 → unseen-data analyses.
5. **Freeze discipline.** The primary G1/G2 results are never rewritten. New analyses add rows; they
   never replace frozen rows.
6. **Gene maps.** Every verdict is computed under the baseline and the runner-up map. Flips are
   reported as "not robust to gene-ID reconstruction".

## 6. Analysis pathways

Effort assumes reuse of `src/clinical_axes/*` and the G1/G2 drivers. Compute for each is minutes, unless
noted.

### A0. Blind gate registration (H1–H4); week 1
- **Inputs:** Appendix A only.
- **Output:** `config/clinical_axes_g3_registration.yaml` (gates, seeds, families) and a hash-lock
  manifest.
- **Gate content the writer must decide:**
  - the τ²_g estimator;
  - the calibrated-interval rule;
  - the agreement rule across control types;
  - matched-null percentile thresholds;
  - the technical-sensitivity "robust" rule.

### A1. Cohort-level variance calibration and control choice (H1); weeks 1–2
- **Estimand.** For each of the 49 sets, plus the 4 axes and the structural control: the pooled
  flight−habitat-ground effect, with each mission's variance inflated by τ²_g, the between-cohort
  variance that animal replication cannot capture.
- **Estimating τ²_g.**
  - From every control-vs-control contrast available: basal−ground and vivarium−ground in OSD-771
    (both arms); vivarium−ground in OSD-253, 462 and 513; the original-vs-rerun ground controls in
    OSD-253.
  - Use a dependence-aware method: a multilevel model with mission and shared-group random effects,
    or one contrast per independent control pair.
  - Never use the raw count of overlapping contrasts (24 of 35 come from OSD-253).
- **Control-choice arms.** Flight−vivarium meta (k = 4) and flight−basal (k = 3), each reported beside
  flight−ground.
- **Method.** REML/mHK with the variance offset. The calibrated max-T uses the existing synchronized
  permutation, with the offset applied in the pooled t.
- **Default gate**, to be confirmed blind:
  - "spaceflight-associated (robust to cohort-level variance)" requires the calibrated 95% CI to
    exclude 0 **and** flight−vivarium to have the same sign with a ratio of at least 0.5 to
    flight−ground;
  - otherwise the label is "shift relative to habitat controls; not robust to control choice".
- **Reuse:** `statistics.random_effects_reml_mkh`, G1 driver loaders. Effort about 4–5 days.

### A2. OSD-771 batch-matched substitution (H1); week 2
- In OSD-771, group equals ERCC mix equals dissection day.
- Re-run the 5-mission meta with OSD-771's flight−ground replaced by flight−vivarium, the mix-matched
  contrast, and separately with OSD-771 dropped.
- Report it for all 49 sets. Effort about 1 day.

### A3. Calibrated, extended matched null (H2); week 2
- K4 already matched the podocyte panel against 10,000 same-tier panels on 7 label-blind variables.
- This analysis adds three things:
  - GC content and gene length as matching variables;
  - a second, genome-wide candidate pool (expressed non-marker genes), because K4's pool was 81/133
    endothelial;
  - A1's calibrated variance in the pooled statistic.
- **Report:** the podocyte percentile in both pools, and the structural control's percentile in the
  genome-wide pool.
- **Gate (blind):** for example, "position special" requires at least the 95th percentile in both pools
  under calibration.
- Effort about 2–3 days, reusing the K4 code.

### A4. Technical robustness bundle (H4); weeks 1–2
Each item is reported as robust or not robust (direction kept and retention ≥ 0.5). None of them can
upgrade a verdict.

1. **Label-invariant eligibility.** Replace the arm-aware rule with pooled-sample rules: CPM ≥ 0.1 in
   ≥ 25% of all mission samples, plus a stricter 50% version. Re-run the 49-set max-T.
2. **GC/length normalization.** CQN- or EDASeq-style correction using Ensembl GC and length.
3. **Orthogonal single-sample scoring.** A rank-based score (singscore-type) instead of mission z-means.
4. **RIN adjustment** for the OSD-771 LAR arm (descriptive; G2's label is unchanged).
5. **Sampling covariates.**
   - Region panels are derived *by rule* from GSE129798 (top-N region-specific genes, excluding
     target-set genes).
   - Slopes are fitted in control animals only.
   - Apply them only to sets not overlapping the panels.
   - Report robust or not robust; never "rank after adjustment".

Effort about 4–5 days.

### A5. Prior-claim recurrence ledger (H5); weeks 2–3, unseen
- **Enumerate the claims before scoring:**
  - Siew 2024's top kidney pathways and most-frequent gene products (PODXL, NID1, FN1, CTSB, APOH,
    PSAP; transporters NCC/SLC12A3 and others; focal adhesion, tight junction, actin);
  - Suzuki 2022's MHU-3 kidney claims;
  - DCT remodeling;
  - any other enumerable renal claims found by a scoped literature pass.
- Map each claim to a gene set by a frozen rule:
  - pathway → KEGG/GO term genes;
  - gene product → a single-gene contrast;
  - direction as claimed.
- Register the ledger and directional predictions.
- Score with the G1 machinery plus A1 calibration.
- **Verdict per claim:**
  - RECURS: calibrated CI excludes 0 in the claimed direction;
  - DOES NOT RECUR: the CI excludes 0 opposite to the claim, or the estimate is near 0 with an
    interval excluding the claimed magnitude;
  - INDETERMINATE: otherwise;
  - NOT EVALUABLE: insufficient genes.
- This yields the paper's denominator and turns the Siew PODXL conflict into a tested estimand
  difference. Effort about 1 week.

### A6. Same-mouse cross-tissue control (H6); weeks 2–3, unseen
- **Data:** OSD-580 heart. The ISS-T arm mirrors the kidney design (flight on-orbit, ground, basal,
  vivarium).
- **VST.** OSD-580 has no GeneLab VST. Compute DESeq2 `vst(blind = TRUE)` with GeneLab settings, and
  validate the procedure on a kidney mission where GeneLab VST exists (requires R/DESeq2; `DESeq2.R`
  is in the repo).
- **Within the ISS-T arm only** (LAR is right ventricle, a different region):
  1. **Structural program** flight−ground, and its control-vs-control contrasts (basal−ground,
     vivarium−ground).
  2. **Handling panels:** immediate-early genes, blood/hemoglobin, hypoxia. These are defined from
     external lists before download.
  3. **Cohort-artifact test.** Correlate the kidney and heart control-group deviations (basal−ground
     and vivarium−ground) across matched gene sets. Also correlate per-animal structural scores between
     kidney and heart (38 animals; within-group correlation).
- **Brain** (OSD-561/562/613) is secondary: flight−ground for handling panels only.
- **Registered predictions** go in before download.
- Effort about 1 week.

### A7. OSD-457 moderator (H7); week 3, supplement only
- Load the kidney samples (12 of 192) through the project loaders.
- Score the 49 sets, wild-type primary and Nrf2-KO separate, genotype-blocked and descriptive.
- No fitted days-since-return curve, and no main-text figure.
- Effort 2–3 days.

### A8. Protein layer (H8); week 4, unseen
- **Coverage gate first, before any effect is computed:** at least 10 podocyte HS proteins and at least
  50 structural-set proteins quantified in OSD-163.
- If it passes, score podocyte, structural and difference programs as protein z-means, flight−ground.
- OSD-102 (3 GB) only if OSD-163 passes.
- **Never combine p-values with the RNA results** (same animals). Report direction concordance only.
- Effort about 1 week.

### A9. Morphometry and outreach (H9); start week 1, gates nothing
- **Request:** WT1+ podocytes per glomerulus, glomerular density and tuft area, cortical sampling
  fraction, DCT density, and ECM/cytoskeletal staining on RR-10 and RRRM-2 sections.
- **Send to:** the PI, the Siew consortium and NASA LSDA/ALSDA biospecimens.
- **Also send:**
  - the GEO/author report of the GSE295428 duplicate;
  - an OSDR query about RRRM-1 bulk kidney mRNA.

## 7. Outcomes and what they mean

### 7.1 Per hypothesis

| Hypothesis | Outcome | Interpretation | Paper consequence |
|---|---|---|---|
| H1 | calibrated CI excludes 0 **and** vivarium agrees (≥ 0.5) | the shift is spaceflight-associated beyond cohort variation | strongest version of the recurrence claim |
| H1 | flight−ground holds but calibration or vivarium fails | a shift relative to habitat controls, not robust to cohort variation or control choice | the recurrence claim is narrowed; the calibration becomes a headline methods result |
| H1 | flight−vivarium opposes flight−ground | habitat or housing effect, not spaceflight | the reproducibility message dominates: control choice drives the signal |
| A2 | sign holds with OSD-771 mix-matched or dropped | not driven by OSD-771 batch structure | one robustness row |
| A2 | sign lost | OSD-771 batch structure drives the pooled effect | a major caveat; OSD-771 downgraded |
| H2 | ≥ 95th percentile in both pools under calibration | the podocyte program's position is special relative to matched genes | "podocyte-leaning" stays as a descriptive top row |
| H2 | fails in either pool | the podocyte position is a generic coherent-set shift | the paper talks about a structural band only |
| H4 | all items robust | direction is not an artifact of eligibility, normalization, scoring or RNA quality | robustness table |
| H4 | any item fails | name the item; that dependency becomes a stated limitation | narrower wording |
| H5 | per claim | RECURS / DOES NOT RECUR / INDETERMINATE / NOT EVALUABLE | the recurrence map (main figure) |
| H5 | PODXL/NID1 up vs Siew down | an estimand difference: animal-level set abundance vs dataset vote count | stated as a difference, not a correction |
| H6 | structural shift in kidney, not in heart | the drift is kidney-attributable | strengthens the kidney claim |
| H6 | structural shift in both | systemic or handling-related; handling panels decide which | reframe as an organism-level or handling component |
| H6 | ground-group deviation shared across organs | cohort-level (housing, handling, day) artifact | validates A1's calibration logic |
| H6 | deviation not shared | tissue-processing batch | calibration still needed; cause is library-level |
| H6 | handling panels elevated in flight heart | an on-orbit dissection/euthanasia signature | explains the G2 exploratory hits; a handling caveat for all on-orbit kidney data |
| H7 | same sign (wild-type) | consistent with the shift at about 2 d post-landing | supplementary row only |
| H7 | opposite or null | live-return heterogeneity | supplementary row only |
| H8 | protein agrees in sign | cross-layer consistency (same animals) | supporting panel |
| H8 | protein disagrees, or the coverage gate fails | RNA-level only, or not evaluable | a stated limitation |
| H9 | (a) / (b) / (c) | composition / per-cell transcription / structural remodeling | Paper II material |

### 7.2 Overall decision: what the paper says

| H1 | H2 | Central statement |
|---|---|---|
| pass | pass | "A podocyte-leaning structural transcriptional shift recurs across missions and survives cohort-level calibration; most targeted clinical axes do not recur." |
| pass | fail | "A broad structural/compartment band recurs and survives calibration; no single compartment program is distinguishable from matched sets." |
| fail | any | "Apparent recurrence of renal transcriptional shifts does not survive calibration for cohort-level variation; control choice and handling explain much of the reported signal." This is a methods and reproducibility result. |

The title, *"Which reported renal transcriptomic responses to spaceflight recur across independent mouse
missions?"*, fits every row. H6 modifies attribution in every row: kidney-specific vs systemic or
handling.

### 7.3 What would genuinely damage the paper

- An error in the frozen primary pipeline.
- A1 *and* A2 *and* A4 all reversing the direction of the shift.
- The ledger (A5) proving unscorable for most claims, which would leave no denominator.

Intervals crossing zero alone do not damage it; they narrow the language.

## 8. Parked questions and what would reopen them

| Question | Why parked | Data that would reopen it |
|---|---|---|
| Recovery after live return | RRRM-2 aliases recovery with shorter flight (32 vs 53–56 d) and euthanasia site; OSD-253 shows exposure-dose dependence | a terminal arm matched to 32 d; corrected RRRM-1 live-return kidney single-cell data (opposite aliasing: 22–24 d terminal vs 39.5 d + 2–4 d return); RRRM-1 bulk kidney mRNA |
| Podocyte separability, level 1 (statistical) | K0: P1 interval crosses 0; the map-dependent max-T | new missions; no further in-sample model tuning |
| Podocyte separability, levels 2–3 (external and cell-resolved) | no usable cell-resolved flight kidney data | corrected GSE295428 LAR kidney; spaceflight single-nucleus or spatial kidney; morphometry (A9) |
| Second atlas | the title is no longer podocyte-centric | a reviewer request; use glomerular references (GSE146912, GSE111107) plus whole-kidney GSE129798 |
| Duration dependence | the cross-study grid crosses strain, sex and lab | a matched-design deposit (Bion-M2, future RR missions) |

## 9. Manuscript plan

- **Title:** *Which reported renal transcriptomic responses to spaceflight recur across independent
  mouse missions?*
- **Claims, in order:**
  1. Recurrence map: ledger, 4 axes, compartment band.
  2. Robustness: G1, A3, A4.
  3. Attribution: A1, A2, A6.
  4. The interpretation boundary: K0 non-separability; bulk cannot resolve composition.
  5. Return-flight limitation: G2, parked.
- **Main figure:** the recurrence map.
  - Rows: ledger claims, 4 axes, structural and podocyte sets.
  - Columns: missions, a calibrated meta column and I².
  - Each row carries its verdict.
- **Named methods result.**
  - Animals replicate within cohorts, but each mission's exposure environment and matched control
    environment are each realized once (one cohort, one dissection day, one processing batch).
  - Cohort-level deviations are therefore estimated from control-vs-control contrasts across
    missions.
  - Do not write "n = 1".
- **Excluded from the paper:**
  - podocyte as hero;
  - recovery or rebound;
  - a days-since-return curve;
  - rank-after-adjustment;
  - mechanism stories;
  - un-rerun post-hoc numbers;
  - the ccRCC H&E project.
- **Sequence:**
  1. Hard freeze: tag `v1.0-manuscript` with configs, hashes, figure source tables, environment lock
     and CI.
  2. Internal hostile review (one quantitative and one biological reader).
  3. A compact package to Dr. Casaletto: skeleton, figures, repository, and two questions (is this a
     contribution, and what would a reviewer attack).
  4. Preprint, then a journal: npj Microgravity first; Communications Biology if A1 and A5 land
     cleanly; Life Sciences in Space Research as the alternative.
- **Drafting gate.** In weeks 1–2, draft only the Introduction, Methods and 4-axis Results. Write no
  structural or podocyte prose until A1–A3 and A5 resolve.

## 10. Timeline and decision points

| Week | Work | Decision point |
|---|---|---|
| 1 | A0 blind gates and hash-lock; A9 requests sent; A5 claim enumeration and predictions; A6 and A8 predictions; start A4 | gates locked before any A1–A4 run |
| 1–2 | A1, A2, A4; draft Introduction, Methods, 4-axis | **D1: H1 label** |
| 2 | A3 | **D2: H2 label → choose the §7.2 row** |
| 2–3 | A5 ledger; A6 heart | **D3: attribution (kidney vs systemic)** |
| 3 | A7 OSD-457 (supplement) | — |
| 4 | A8 OSD-163 coverage gate, then effects; OSD-102 if it passes | — |
| 4–6 | Results prose by §7.2 row; figures; supplement | — |
| 6 | Software freeze (environment lock, CI, Docker, README); `v1.0-manuscript` tag | hostile review → Casaletto → preprint (target December 2026) |

## 11. Compute and agent budget

- **Code and compute.** All analyses reuse existing loaders and drivers. The largest downloads are
  OSD-102 proteomics (3 GB, conditional) and the GSE129798 reference.
- **Agent use** (`docs/strategy/2026-10-07/README.md`):
  - at most 2–3 agents in flight;
  - commit every result before the next stage;
  - Sonnet for reading and implementation of well-specified packages; Opus for judging and statistical
    design;
  - compact handoffs;
  - the A0 gate-writer is a fresh agent given only Appendix A.

## 12. Risks

| Risk | Mitigation |
|---|---|
| H1 fails, so "spaceflight-associated" is lost | Expected and acceptable; §7.2 row 3 is a publishable methods result |
| Gate-writer contamination | Fresh agent; Appendix A only; hash-lock; disclosed in provenance |
| τ²_g estimate unstable (few independent control pairs) | Pre-specify the estimator and a sensitivity (upper-bound τ²); report both |
| The ledger is subjective (claim → gene-set mapping) | Frozen mapping rule; register before scoring; a second person audits the mapping |
| Heart VST computed differently from GeneLab | Validate on a kidney mission with GeneLab VST before use |
| Mediator/collider bias in sampling covariates | Rule-derived external panels, control-only slopes, robust/not-robust reporting only |
| Scope creep | §7.3 defines what matters; no new analyses unless they correct an error, answer a concrete reviewer objection, or form a separate project |

## 13. Housekeeping before the freeze

- Call CPM eligibility "arm-aware", not "label-blind" (`PROJECT_TECHNICAL_AUDIT_2026-08-25.md`
  line 595).
- Structural I² is 0.67%, not 67.3%.
- The structural control ranks 2nd, not 4th.
- "No sixth mission" is wrong (OSD-457).
- Record the K0, G1/G2, item-9 and this plan's decisions in `docs/owner_decisions.md`.
- Refresh the README (headline, test counts).
- Unify the Python and R environments and pin GitHub dependencies.
- Add CI.
- Retire the stale Docker entrypoint.
- Delete or rewrite the Casaletto email draft.

---

## Appendix A. Blind brief for the G3 gate-writer (contains no observed values)

You are writing the decision gates for registered sensitivity analyses of an existing cross-mission
kidney transcriptomics result. Do not read any file except this appendix and the code paths named
below. Do not look at results directories.

- **Design.**
  - Five terminal spaceflight missions; animals are the unit; flight vs habitat ground control is the
    designed contrast.
  - Some missions also have vivarium and/or basal control groups.
  - In each mission each group is one cohort, dissected and processed together, so cohort-level
    deviations are not estimable from animal replication.
  - Effects are Hedges g per mission (age/duration strata combined within mission), pooled by
    REML/modified Hartung–Knapp.
  - Max-T family-wise correction is applied over a 49-set family via synchronized blocked label
    permutation.
- **Estimands to gate:**
  1. **Cohort-variance-calibrated pooled effect.** Add τ²_g, estimated from control-vs-control
     contrasts, to each mission's variance. You choose the dependence-aware estimator and an
     upper-bound sensitivity.
  2. **Control-choice agreement.** Flight−vivarium and flight−basal pooled effects compared with
     flight−ground. You choose the sign and ratio rule.
  3. **Batch-matched substitution.** In the one mission where group equals processing batch, the
     flight contrast is replaced by its batch-matched alternative. You choose what counts as "holds".
  4. **Matched-null position.** The target set's percentile among 10,000 matched gene sets, in two
     candidate pools. You choose the threshold, and whether both pools are required.
  5. **Technical sensitivities:** label-invariant eligibility, GC/length normalization, rank-based
     single-sample scoring, RIN adjustment, rule-derived sampling covariates. You choose the
     "robust" rule (sign and retention).
- **Constraints.**
  - No gate may upgrade a frozen verdict.
  - Every gate is applied identically to all 49 sets, the 4 axes and the structural control.
  - Labels are "registered sensitivity (post hoc in origin)".
- **Output:** `config/clinical_axes_g3_registration.yaml` with gates, seeds (parent seed plus new
  offsets) and families. Return it for hash-locking before any code runs.

## Appendix B. How this plan handles the owner's external assessment (2026-10-07)

| Owner recommendation | Disposition in this plan |
|---|---|
| Inventory GSE295428 for podocytes before more bulk work | **Done 2026-10-07.** "Kidney, LAR" is a duplicate spleen upload, so the path is closed (§4.3, §8). Report the duplicate (A9). |
| Label-invariant (pooled-sample) eligibility sensitivity | Adopted: A4.1 |
| Independent second atlas | Deferred (§8). The title is no longer podocyte-centric. If run, use glomerular references plus GSE129798 rather than Tabula Muris Senis. |
| Orthogonal single-sample scoring | Adopted: A4.3 |
| RIN-adjusted LAR sensitivity | Adopted: A4.4 (post-hoc preview 0.14 → 0.17; G2 label unchanged) |
| Lock the computational environment and CI | Adopted: §13 and the week-6 freeze |
| Pursue RRRM-1 bulk kidney mRNA through OSDR | Adopted as an outreach query (A9); no analysis depends on a reply |
| GSE295428 animal-level pseudobulk validation | Closed until a corrected live-return kidney deposit (§8) |
| OSD-913 kidney miRNA, narrow regulatory question | Not scheduled: orthogonal but gives no cell localization; it carries the RRRM-1 dissection-method confound |
| Narrow podocyte-vs-structural proteomics (OSD-102/163/462) | Adopted for OSD-163, then OSD-102 (A8). **OSD-462 is excluded**: condition is perfectly aliased with TMT reporter-tag block (the reason the OSD-462 phospho line was closed). |
| OSD-513 as a heterogeneity/dissociation case | Kept as a descriptive moderator; never recovery evidence |
| Latent structural-axis or factor analysis | Optional, not scheduled; only as a new predefined question |
| Duration meta-regression | Parked (§8): the cross-study grid crosses strain, sex and lab |
| Matched random gene sets | Adopted as A3, **extending K4** (which already ran 10,000 same-tier matched panels) with GC/length, a genome-wide pool and calibrated variance |
| Decompose the endothelial/TAL/stromal structural components | Not scheduled |
| Hard freeze → hostile review → Casaletto → preprint → npj Microgravity | Adopted (§9) |
| Stop list: P1 tuning, OSD-513 as recovery, podocyte injury, redundant sensitivities, p-threshold robustness, silent changes to frozen results, merging the ccRCC project | Adopted (§2.3, §5, §9) |
| Title "Which reported renal transcriptomic responses … recur …" and careful cohort-replication wording | Adopted (§9) |
| §4 of the assessment: the negative ISS-T/LAR correlation as "suggestive of rebound" | **Superseded.** Basal−ground reproduces it without flight animals, so it is a cohort artifact (§2.2) |
| Go/no-go judged on direction, not on p-value thresholds | Adopted (§7.3), with one addition. Because 32/49 sets shift positive, direction alone cannot establish specificity, so H2 (A3) and H1 calibration carry that weight. |
