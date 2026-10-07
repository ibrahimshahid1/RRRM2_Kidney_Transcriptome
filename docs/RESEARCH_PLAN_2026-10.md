# Research plan: cross-mission recurrence of renal transcriptomic responses to spaceflight

**Version:** 1.1 (2026-10-07). This revision applies every required fix from an adversarial review
(`docs/strategy/2026-10-07/PLAN_REVIEW.md`). It supersedes the "next steps" sections of earlier
decision documents.

**Status:** plan of record. Each analysis is registered separately (§5) before it runs.

**Consolidates:**
- `docs/NEXT_DIRECTIONS_2026-10-07.md`;
- `docs/strategy/2026-10-07/` (reader notes, draft, two adversarial reviews, GSE295428 feasibility);
- the owner's external assessment (crosswalk in Appendix B);
- `docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md`.

> **Do not give this document to the G3 gate-writer or the G3 executor.** It contains observed values
> and previews. Each of them receives only Appendix A.

---

## 1. Summary

**Central question.** *Which reported renal transcriptomic responses to spaceflight recur across
independent mouse missions, and which survive calibration for cohort-level (habitat, handling, batch)
variation?*

**Established.** Across five terminal missions (83 animals):

- **The four clinically anchored renal axes do not recur.** The glomerular-barrier axis even moves
  opposite to the declared "identity loss".
- **A low-heterogeneity structural/compartment band shifts upward.** An externally defined 157-gene
  podocyte-associated program sits at its top (g = 0.689). It is the only one of 49 sets whose interval
  excludes zero.
- **Robustness and limits.**
  - Its direction survives every preregistered perturbation (G1).
  - It is not separable from the shared structural shift (K0).
  - Persistence after live return cannot be identified (G2).

**Open questions.**

1. **Attribution.** Does any shift survive calibration for cohort-level variance and alternative control
   groups?
2. **Specificity.** Is the podocyte program's position extreme relative to matched, equally coherent
   gene sets?
3. **Literature recurrence.** Which specific published claims recur, under an animal-level estimand, in
   the same missions the literature used?
4. **Tissue attribution.** Do the structural shift and the control-group deviations also appear in the
   heart of the same mice?

**Expected outcome, stated in advance.** Post-hoc previews (§2.2) imply that the podocyte set will
**not** pass the calibrated attribution test. The most likely paper is therefore the
methods/reproducibility version (§7.2, rows R2–R4). Every row answers the central question.

**Parked.** Recovery and podocyte-specific biology are not identifiable with public data. §8 lists what
would reopen each.

## 2. What is established, found post hoc, and retired

### 2.1 Established (frozen, preregistered or locked; reproduced 2026-10-06)

| Result | Numbers | Source |
|---|---|---|
| 4 clinical axes, 5 missions | barrier −0.716 (−1.361, −0.072), max-T 0.0018 (opposite to declared); tubular injury 0.570 (FWER 0.026); fibrosis 0.311 (0.785); distal −0.154 (0.986) | clinical-axes revision |
| 49-set compartment family | podocyte HS 0.689 (0.042, 1.336), max-T 0.019, I² 0; mission effects −0.12 / 1.17 / 0.68 / 0.54 / 1.17; structural control 0.549 (**−0.095**, 1.193); ranks 2–5 within 0.02; 32/49 sets positive | compartment-context |
| K4 strict matching | 10,000 balanced same-tier 142-gene panels: none as extreme (p = 0.0001); one-to-one target-minus-matched g 0.54 (−0.17, 1.25), blocked p 0.015 (81/133 of the *trimmed one-to-one matches* endothelial) | `PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md` |
| K0 specificity | P0 0.198 → P1 0.112 (−0.031, 0.255) score units; S1 0.031; Freedman–Lane P1 max-T 0.042 (baseline map) / 0.056 (runner-up) → **not separable** | K0 stage |
| G1 sensitivities | ROBUST_DIRECTION (items 1–7); DEFINITION_ROBUST (8/8 atlas LOSO); no gene-map flips; leave-top-5 lower CI −0.044 | G1/G2 results doc |
| G2 recovery | podocyte INCONCLUSIVE (ISS-T g 1.17, LAR 0.14; interaction −0.28 (−0.62, 0.06)); structural interaction −0.35 (−0.63, −0.07), FL p 0.009, max-T 0.024 (INCONCLUSIVE after Holm) | G1/G2 results doc |

### 2.2 Found after G1/G2 (post hoc, unregistered: motivation and previews, not results)

| Finding | Numbers | Consequence |
|---|---|---|
| Control choice | podocyte flight−vivarium meta 0.28 (−0.45, 1.01), k = 4, C57BL/6J; ratio to flight−ground 0.41. Structural flight−vivarium 0.55 | previews A1 → H1 |
| Cohort-level variance | control-vs-control contrasts imply about 0.36 excess between-group variance (an upper bound: housing and age effects are included, and contrasts overlap). The podocyte calibrated interval is about (−0.31, 1.68). At that τ², k = 5 needs g ≥ ~1.0 to exclude 0 | previews A1 |
| OSD-771 batch structure | **group = ERCC spike-in mix** (ISS-T: flight and vivarium Mix 1, ground and basal Mix 2; reversed in LAR). Dissection dates **partly overlap** (flight 16–19 Sep; ground 18–20 Sep; vivarium 18–20 Sep); basal animals were dissected 22–23 Jul. Podocyte flight−vivarium in ISS-T, the mix-matched contrast, is −0.17 | previews A2 |
| G2 "rebound" (r = −0.65) | reproduced by basal−ground (−0.65) and vivarium−ground (−0.48) without flight animals | **refuted**: a cohort artifact |
| G2 exploratory interaction hits | ISS-T flight carries on-orbit handling signatures (immediate-early genes g 2.14, medulla panel 2.71, blood 1.73) | **retired** as biology |
| Technical correction | GC/length/expression residualization keeps podocyte rank 1 (g 0.65–0.71) | previews A4.2 |
| RIN | adjustment moves LAR podocyte from 0.14 to 0.17 | previews A4.4 |
| Sixth kidney cohort | OSD-457 (JAXA MHU-3), scored descriptively by a reader's standalone pipeline (validated on OSD-513, r = 0.995): podocyte g 1.16 (wild-type 1.36) | previews A7 |
| Eligibility rule | `cpm_eligible_genes` selects on flight/ground membership: "arm-aware", not label-blind | A4.1 |
| Single genes already seen | per-mission effects for the 4 axes' genes, including Podxl, Nid1, Slc12a3 and Fn1, exist in `descriptive_signed_gene_effects.tsv`; K4 noted that neither PODXL nor NID1 alone was a significant cross-mission effect | A5 is partly seen |
| Prior art (Siew 2024, full text) | qualitative nephron structural remodeling (focal adhesion, tight junction, actin) by enrichment and vote count; PODXL/NID1 among the most frequently *down*-regulated products. **Siew's 11 mouse missions include RR-1, RR-3, RR-7, RR-10 and RR-23** (named in its acknowledgments), which are our OSD-102, 163, 253, 462 and 513 | H5 is re-analysis under a different estimand, **not** independent replication |
| GSE295428 (RRRM-1 single-cell) | "Kidney, LAR" is spleen and duplicates "Spleen, LAR" (99.6–99.9% shared barcodes); kidney exists only for habitat controls | **closed** (§8) |

### 2.3 Retired: do not reopen

- network rewiring
- contrast-vector framing
- LAR reversal and any "rebound"
- Grey60 rescue
- OSD-462 phospho mining (stop rule)
- DCT2/CNT and ASDN headlines
- the OSD-462 data note, the matched-library paper and the merged methods paper
- the existing Casaletto email draft
- podocyte as headline or title
- gene-list mechanism stories without preregistered sub-modules

## 3. Hypotheses and questions

**Data-status labels.**
- **Seen:** previewed post hoc, so the analysis is a registered sensitivity analysis (post hoc in
  origin).
- **Unseen:** no one has computed it, so a directional prediction can be registered honestly.

| ID | Statement (object) | Test | Data status |
|---|---|---|---|
| **H1** | A set's recurrent shift is spaceflight-associated: robust to cohort-level variance, control choice and OSD-771 batch structure. Evaluated **per set**. Reported for the podocyte set, the structural control, the 4 axes and all 49 sets | A1 + A2 | **Seen** (podocyte; expected fail) / partly seen (structural) / unseen (others calibrated) |
| **H2** | The podocyte set's position is extreme among matched, equally coherent gene sets | A3 (extends K4) | Partly seen (K4 without coherence matching or calibration) |
| **H3** | A podocyte residual is separable from structural drift | **Closed in-sample** (K0: not separable). A3, A4 and A8 cannot reopen it | Closed |
| **H4** | The podocyte and structural directions are not artifacts of eligibility, normalization, scoring or RNA quality | A4 (verdict restricted to podocyte and structural) | Partly seen |
| **H5** | Specific published renal spaceflight claims recur under an animal-level estimand in the same missions | A5 (ledger) | Unseen except the 4-axis single genes |
| **H6a** | The structural shift is absent in the heart of the same mice (kidney-attributable) | A6 | **Unseen** |
| **H6b** | Control-group deviations (basal−ground, vivarium−ground) are shared across organs of the same mice (organism-wide cohort effect) | A6 | **Unseen** |
| **H7** | A live-return cohort sampled about 2 days after landing shows the same sign | A7 (supplement) | Seen descriptively |
| **H8** | Within a mission, podocyte and structural protein shifts follow that mission's RNA direction | A8 (OSD-163, then OSD-102) | **Unseen** |
| **Q9** (outreach question; no test in this plan) | Is the bulk podocyte signal (a) more podocytes, (b) more transcription per podocyte, or (c) structural/compositional change? | A9 (morphometry request) | Long lead |

## 4. Datasets

### 4.1 In hand (inputs verified 2026-10-05/06)

| Dataset | Design | Role |
|---|---|---|
| OSD-102 (RR-1) | female C57BL/6J, 37 d, terminal; 6 flight / 6 ground | primary |
| OSD-163 (RR-3) | female BALB/c, 39–42 d, terminal; 6 / 6 | primary; mapping-rate flag |
| OSD-253 (RR-7) | C57BL/6J and C3H/HeJ, 25 and 75 d; 10 / 9 (B6); original and rerun ground controls | primary; **clean τ²_cohort pair**; duration and strain context |
| OSD-462 (RR-10) | B6129 female, 28–29 d; total RNA, mRNA and UPX libraries; 10 / 10; has vivarium controls | primary; library-preparation sensitivity |
| OSD-771 (RRRM-2) | female B6NTac; ISS-T 53–56 d and LAR 32 d + 24 d; flight, ground, basal and vivarium groups, 5 per age block | primary (ISS-T); G2; control choice |
| OSD-513 (RR-23) | male C57BL/6J, 37 d, euthanized about 1 d after landing; 9 / 9; has vivarium controls | descriptive moderator |
| Mouse Kidney Atlas (Zenodo 17395591); Tabula Muris Senis kidney; KEGG_2019_Mouse; both gene maps | — | frozen definitions |

### 4.2 To acquire (files verified on OSDR 2026-10-07)

| Dataset | Files (note the GLDS prefixes) | Size | Design | Analysis |
|---|---|---|---|---|
| **OSD-580** (RRRM-2 heart; 78/78 animal IDs shared with kidney: 38 ISS-T + 40 LAR) | `GLDS-573_rna_seq_RSEM_Unnormalized_Counts.csv`, `GLDS-573_rna_seq_Normalized_Counts.csv`, runsheet, SampleTable. **No GeneLab VST** | 13–31 MB | ISS-T whole heart: flight 10 (on orbit), ground 10, basal 10, vivarium 8. LAR: right ventricle, 10 per group | A6 |
| OSD-561 / 562 / 613 (RRRM-2 cerebellum / hippocampus / cerebrum) | `GLDS-556_*` (OSD-561) VST, counts, qc; the others similar | 5–12 MB each | flight vs ground, pooled controls, metadata irregularities | A6, descriptive only |
| **OSD-457** (JAXA MHU-3) | `GLDS-457_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv`, v2 runsheet, qc_metrics | 31 MB | male C57BL/6J, 31 d, live return, about 2 d after landing, after behavioural tests; wild-type 3 / 3, Nrf2-KO 3 / 3 (kidney is 12 of 192 samples) | A7 |
| **OSD-163 proteomics** | `GLDS-163_proteomics_RR3_KDN_processed.zip` | 21 MB | intended to be the same animals as the OSD-163 RNA (sample-ID linkage to be confirmed) | A8 (first) |
| OSD-102 proteomics | `GLDS-102_proteomics_TMT_KDNL.processed.tar.gz` | 3.0 GB | same, TMT | A8 (only if OSD-163 is informative) |
| Ensembl gene GC% and length | BioMart, release matching the gene map | small | — | A3, A4 |
| Region marker reference | GSE129798 (Ransick 2019: cortex, outer and inner medulla) | moderate | external, label-free | A4.5 |
| Prior-claim sources | Siew 2024 Supplementary Data; Suzuki 2022 (Kidney Int, MHU-3); other enumerable renal claims | small | — | A5 |

### 4.3 Watchlist and closed paths

| Item | Status | What would change it |
|---|---|---|
| GSE295428 (RRRM-1 single-cell) | **closed**: "Kidney, LAR" is a duplicate spleen upload | a corrected live-return kidney deposit (report the duplicate: A9) |
| RRRM-1 bulk kidney mRNA | not found on OSDR 1–1150 or GEO (miRNA only, OSD-913) | an OSDR reply |
| Bion-M2 (2025); Siew GCR kidney (OSD-706..712) | not public | deposits |
| GSE235042 (hindlimb-unloading kidney) | dropped (n = 3 per group, FPKM only, sex unknown) | — |
| Second atlas (GSE146912 / GSE111107 + GSE129798) | deferred | a reviewer request, or "podocyte" back in the title |

## 5. Governance: making the new tests honest

1. **Blind gate-writing (A0).**
   - A fresh agent or person who has not read §2.2, `reader_leads.md` or any scratchpad output writes
     every decision threshold for A1–A5 from Appendix A alone.
   - **The blind writer's YAML is adopted verbatim.** The owner may reject it only for internal
     inconsistency, in writing, before any run. A rejection triggers a new blind writer, not an edit.
   - The YAML is committed as `config/clinical_axes_g3_registration.yaml`, and the code is hash-locked
     through the rrrm2 CLI.
2. **Separate executor.**
   - A1–A4 are coded and executed by an agent that has not seen §2.2.
   - Output is opened only after the hash-locked run completes.
   - All runs are logged, failed ones included.
3. **Known-outcome disclosure (recorded before gates run).** From post-hoc previews, the podocyte
   set's A1/A2 result is **expected to fail** under any reasonable gate:
   - vivarium ratio about 0.41;
   - calibrated interval about (−0.31, 1.68);
   - OSD-771 mix-matched contrast −0.17;
   - the minimum detectable calibrated effect at k = 5 is about g ≥ 1.0.

   The structural control's uncalibrated interval already crosses 0, so it cannot pass a calibrated
   gate; its flight−vivarium preview is 0.55. A1/A2 are genuinely informative for the 4 axes
   (notably the barrier axis) and the remaining sets, whose calibrated values have not been previewed.
4. **Directional predictions before download** for the unseen analyses (A5, A6, A8) are committed
   before any data are fetched.
5. **Labels.**
   - A1–A4 are "registered sensitivity analyses (post hoc in origin)".
   - A5, A6 and A8 are "preregistered secondary analyses".
   - A7 is "descriptive".
   - Nothing is "confirmatory".
6. **Provenance ledger in the paper:** targeted axes → 3/4 no-go → 49-set scan → podocyte → K1–K4 →
   G1/K0/G2 → reader post-hoc previews → G3 → unseen-data analyses.
7. **Freeze discipline.** Primary G1/G2 results are never rewritten; new analyses add rows.
8. **Gene maps.** Every verdict in A1–A5 is computed under both maps. If the maps disagree, the more
   conservative label is reported, together with "not robust to gene-ID reconstruction".
9. **Scope control.** A9 replies, GEO corrections or new deposits that arrive before the freeze are
   logged, not analysed.

## 6. Analysis pathways

Thresholds below are written as `[blind]`: the gate-writer sets them (§5.1). Effort assumes reuse of
`src/clinical_axes/*` and the G1/G2 drivers.

### A0. Blind registration (week 1)

- **Inputs:** Appendix A only.
- **Output:** `config/clinical_axes_g3_registration.yaml`, containing the gates for A1–A5, seeds,
  families and multiplicity rules. Then the hash-lock manifest.

### A1. Cohort-variance calibration and control choice (H1; coding week 1, execution week 2)

**Variance model.** Total between-mission variance = **max(τ̂²_REML, τ̂²_cohort)**.
- τ̂²_cohort comes from a joint REML model. Control-vs-control contrasts enter it as extra observations
  that share the cohort variance component, with a random mission intercept and a covariance term for
  contrasts that share a control group.
- **The clean estimator** is same-condition replicate cohorts: OSD-253 original vs rerun ground
  controls.
- **Upper-bound sensitivity only:** vivarium−ground (which contains a housing effect) and basal−ground
  (which contains about 2 months of age and time), each with a fixed effect for contrast type.
- τ²_cohort is assumed exchangeable across missions and applied to missions without extra controls
  (OSD-102, OSD-163). This assumption is stated.
- **Common across sets** (a set-specific τ² is not identifiable from about 5 effective pairs).
  Set-specific values are a sensitivity.
- Report the χ²-based CI of τ̂²_cohort, its effective degrees of freedom, and the range of
  leave-one-control-contrast-out τ̂² (a calibration-validity check).

**Control-choice arms.**
- Flight−vivarium meta (k = 4) and flight−basal (k = 3, reported only).
- Shown beside flight−ground for every set.

**Inference.** REML/mHK with the calibrated variance. The calibrated max-T uses the existing
synchronized permutation, with the variance offset applied in the pooled t.

**H1 label per set** (thresholds `[blind]`):

| Uncalibrated CI excludes 0 | Calibrated CI excludes 0 | Flight−vivarium sign/ratio | A2 keeps sign | Label |
|---|---|---|---|---|
| — | yes | same sign, ratio ≥ r* | yes | **spaceflight-associated (calibrated)** |
| — | yes | fails | yes | calibrated; control-dependent |
| — | yes | any | no | calibrated; OSD-771-dependent |
| yes | no | any | any | habitat-relative shift; not robust to cohort variation |
| no | no | any | any | not established (direction only) |

**Overlay.** If the flight−vivarium CI excludes 0 *opposite* to flight−ground, the set is labelled
"habitat/housing effect". An opposite-signed point estimate alone does not qualify.

Effort about 4–5 days, plus τ² model validation.

### A2. OSD-771 batch-matched substitution (H1; week 2)

- Re-run the 5-mission meta with OSD-771's flight−ground replaced by flight−vivarium (the ERCC-mix-
  matched contrast), and separately with OSD-771 dropped.
- Add a dissection-date covariate sensitivity within OSD-771.
- Report it for all sets. "Keeps sign" feeds the H1 label. Effort about 1–2 days.

### A3. Calibrated, coherence-matched null (H2; week 2)

Extends K4:
- **Panel size** 142 (genes observable in all 5 missions, as in K4).
- **Matching variables:** K4's 7 variables, plus GC content, gene length and **mean pairwise
  inter-gene correlation in control animals** (coherence), so the null contains equally coherent sets.
- **Two pools:** K4's same-tier pool, and a genome-wide pool of expressed non-marker genes.
- **Ranking statistic:** calibrated t.
- **Report:** the podocyte percentile in both pools, and the structural control's percentile in the
  genome-wide pool.
- **Gates `[blind]`:** H2 pass for podocyte; whether the structural control's position counts as
  "band is special", or gets no gate.
- A3 cannot reopen H3. Effort about 3 days.

### A4. Technical robustness bundle (H4; coding week 1, execution week 2)

The verdict covers the **podocyte set and the structural control only**; other sets are descriptive.
"Retention" means the ratio of calibrated g to the reference. The rule is `[blind]`.

1. **Label-invariant eligibility:** pooled-sample rules (for example, CPM ≥ 0.1 in ≥ 25% of all mission
   samples, plus a stricter version), then re-run the 49-set max-T.
2. **GC/length normalization:** CQN- or EDASeq-style correction with Ensembl GC and length.
3. **Rank-based single-sample scoring** (singscore-type).
4. **RIN adjustment** for the OSD-771 LAR arm. Descriptive; it does not touch the primary ISS-T meta or
   G2's label.
5. **Sampling covariates.**
   - Region panels derived *by rule* from GSE129798 (top-N region-specific genes, excluding target-set
     genes).
   - Slopes fitted in control animals only.
   - Applied only to sets that do not overlap the panels.
   - Reported as robust or not robust; never "rank after adjustment".

None of these can upgrade a verdict. Effort about 4–5 days.

### A5. Prior-claim recurrence ledger (H5; weeks 2–3; unseen except the 4-axis single genes)

**Enumeration** happens before scoring:
- Siew 2024's top kidney pathways and most-frequent gene products;
- Suzuki 2022's MHU-3 kidney claims;
- DCT remodeling;
- other enumerable renal claims from a scoped literature pass.

Then:
- **Frozen mapping rule:** pathway → KEGG/GO term genes (a set score); gene product → a single-gene
  effect. Claimed direction recorded.
- **Two families:** pathway/set claims use the set τ²; single-gene claims use a gene-level τ². Holm is
  applied within each family, together with the expected number of false RECURS calls.
- **Verdicts** (δ and α `[blind]`):
  - **RECURS:** calibrated CI excludes 0 in the claimed direction.
  - **CONTRADICTED:** it excludes 0 opposite to the claim.
  - **ABSENT:** the calibrated CI lies entirely within (−δ, δ).
  - **INDETERMINATE:** otherwise.
  - **NOT EVALUABLE:** too few genes.
- **No interpretation is pre-written** for any claim, Siew's PODXL/NID1 included.
- **Mapping audit.** A second person audits the claim → gene-set mapping before scoring.
- **Framing.** Siew's missions overlap ours, so for those missions the ledger is a re-analysis under an
  animal-level, per-mission estimand, not independent replication. OSD-771 (RRRM-2) is the likely
  non-overlapping mission; confirm against Siew's mission table before claiming it.

Effort about 1 week.

### A6. Same-mouse cross-tissue analysis (H6a, H6b; weeks 2–3; unseen)

**Data.** OSD-580 heart, using the ISS-T arm only. It mirrors the kidney design (flight on orbit,
ground, basal, vivarium). The LAR arm is right ventricle, a different region, so it is excluded.

**VST.** DESeq2 `vst(blind = TRUE)` with GeneLab settings. **Validation gate:** on a kidney mission that
has a GeneLab VST, per-gene Pearson ≥ 0.99 against GeneLab *and* 49-set mission g identical within
±0.02. Without that, heart is not used.

**H6a.**
- Heart flight−ground for the structural control and the KEGG structural sets.
- Prediction `[registered]`.
- Outcomes:
  - structural CI excludes 0 in heart in the kidney direction → present;
  - CI within the equivalence margin → absent;
  - otherwise indeterminate.

**H6b.**
- **Estimand:** the correlation between kidney and heart *deviation vectors*, separately for
  basal−ground and vivarium−ground, across a registered panel of tissue-agnostic gene sets (KEGG
  structural and housekeeping-adjacent sets; no kidney compartment sets).
- **Null:** gene-set rotation, or permutation of animals within arm.

**Handling panels.**
- Immediate-early genes, blood/hemoglobin and hypoxia, defined from external lists before download.
- Flight−ground and the two control contrasts: 3 panels × 2 contrasts, Holm.

**Per-animal kidney–heart structural correlation** (38 ISS-T animals, within group). Labelled
*animal-level* covariation, not cohort-level.

**Brain** (OSD-561/562/613) is descriptive only. Effort about 1 week.

### A7. OSD-457 moderator (H7; week 3; supplement only)

- Kidney samples (12/192) are loaded through the project loaders, replacing the reader's standalone
  pipeline.
- Wild-type is primary; Nrf2-KO is reported separately, genotype-blocked and descriptive.
- No days-since-return curve and no main-text figure. Effort 2–3 days.

### A8. Protein layer (H8; week 4; unseen)

1. **Confirm the sample-ID linkage** of the OSD-163 proteomics to the OSD-163 RNA animals.
2. **Coverage gate, before any effect:** at least 10 podocyte HS and at least 50 structural-set proteins
   quantified in ≥ 80% of samples.
3. **Compare with the same mission's RNA effect** (OSD-163), not the pooled sign:
   - "agrees": the protein CI excludes 0 in that mission's RNA direction;
   - "disagrees": it excludes 0 in the opposite direction;
   - otherwise "uninformative".
4. The podocyte-minus-structural difference is **descriptive** and is not H3/K0 evidence.
5. OSD-102 (3 GB) only if OSD-163 is informative.

P-values are never combined with RNA. Effort about 1 week.

### A9. Outreach (Q9; start week 1; gates nothing)

- **Morphometry request:** WT1+ podocytes per glomerulus, glomerular density and tuft area, cortical
  sampling fraction, DCT density, and ECM/cytoskeletal staining on RR-10 and RRRM-2 sections. Send to
  the PI, the Siew consortium and NASA LSDA/ALSDA.
- **GEO/author report** of the GSE295428 duplicate upload.
- **OSDR query** about RRRM-1 bulk kidney mRNA.

## 7. Outcomes and what they mean

### 7.1 Per hypothesis (rows are exclusive within each hypothesis)

- **H1:** see the label table in A1. Every set receives exactly one label, plus the optional overlay.
- **H2 (podocyte):**
  - pass in both pools → the position is extreme among equally coherent matched sets;
  - pass in one pool only → pool-dependent (report both);
  - fail in both → a generic coherent-set shift.
- **H4 (podocyte, structural):**
  - all items robust → direction is not an artifact of eligibility, normalization, scoring or RNA
    quality;
  - any item not robust → name it; it becomes a stated dependency.
- **H5 (per claim):** RECURS / CONTRADICTED / ABSENT / INDETERMINATE / NOT EVALUABLE. These are the rows
  of the recurrence map. Claims are interpreted only after scoring.
- **H6a:**
  - present in heart → the shift is not kidney-specific (systemic or handling); handling panels and
    H6b inform which;
  - absent → kidney-attributable;
  - indeterminate → unresolved.
- **H6b:**
  - shared (correlation beyond the null) → an organism-wide cohort effect, supporting A1's calibration
    logic;
  - not shared → the cohort artifact is not organism-wide; cause unresolved (tissue-processing batch
    or tissue-specific housing biology).
- **Handling panels:** elevated in flight heart → an on-orbit dissection/euthanasia signature; a
  caveat for all on-orbit kidney data.
- **H7 (wild-type):**
  - same sign → consistent at about 2 d post-landing;
  - opposite or null → live-return heterogeneity.

  Supplement only.
- **H8:** agrees / disagrees / uninformative / not evaluable (coverage gate failed).

### 7.2 The paper's central statement

The rows are exhaustive. R1 requires the **podocyte set's** H1 label to be "spaceflight-associated
(calibrated)". R2–R4 are keyed on the podocyte set's H2 and on whether *any* set or axis attains a
calibrated label.

| Row | Condition | Central statement |
|---|---|---|
| R1 | podocyte H1 calibrated, H2 pass | "A podocyte-leaning structural shift recurs across missions and survives cohort-level calibration; most targeted clinical axes do not recur." |
| R2 | podocyte H1 not calibrated, H2 pass | "A recurrent habitat-relative shift whose top program is extreme among matched coherent gene sets, but no compartment program survives cohort-level calibration." |
| R3 | podocyte H1 not calibrated, H2 fail; ≥ 1 other set or axis calibrated | "Only [named sets/axes] survive cohort-level calibration; the apparent podocyte-led band is a generic habitat-relative shift." |
| R4 | podocyte H1 not calibrated, H2 fail; no set or axis calibrated | "No renal transcriptional response is established beyond direction once cohort-level variation is calibrated; control choice and handling explain much of the apparent recurrence." |

**Expected row:** R2 or R4 (from §5.3). The title, *"Which reported renal transcriptomic responses to
spaceflight recur across independent mouse missions?"*, fits every row. H5 supplies the denominator,
and H6a/H6b qualify attribution in every row.

### 7.3 What would genuinely damage the paper

- An error in the frozen primary pipeline.
- **Any two of A1, A2 and A4** reversing the sign of the podocyte set or the structural control. That
  would make even the direction claim unreliable.
- The ledger (A5) proving unscorable for most claims.

Intervals crossing zero alone narrow the language; they do not damage the paper.

## 8. Parked questions and what would reopen them

| Question | Why parked | Data that would reopen it |
|---|---|---|
| Recovery after live return | RRRM-2 aliases recovery with shorter flight (32 vs 53–56 d) and euthanasia site; OSD-253 shows exposure-dose dependence | a terminal arm matched to 32 d; corrected RRRM-1 live-return kidney single-cell data (opposite aliasing); RRRM-1 bulk kidney mRNA |
| Podocyte separability, statistical | K0: P1 interval crosses 0; map-dependent max-T | new missions only; no in-sample model tuning |
| Podocyte separability, external or cell-resolved | no usable cell-resolved flight kidney data | corrected GSE295428 LAR kidney; spaceflight single-nucleus or spatial kidney; morphometry (Q9) |
| Second atlas | the title is no longer podocyte-centric | a reviewer request (glomerular references + GSE129798) |
| Duration dependence | the cross-study grid crosses strain, sex and lab | a matched-design deposit |

## 9. Manuscript plan

**Title:** *Which reported renal transcriptomic responses to spaceflight recur across independent mouse
missions?*

**Claims:**
1. **Recurrence map:** ledger, 4 axes, compartment band; row R1–R4 statement.
2. **Robustness:** G1, A3, A4.
3. **Attribution:** A1, A2, A6.
4. **Interpretation boundary:** K0 non-separability; bulk cannot resolve composition.
5. **Return-flight limitation:** G2, parked.

**Main figure: the recurrence map.**
- Rows: ledger claims, 4 axes, structural and podocyte sets.
- Columns: missions, a calibrated meta column and I².
- Each row carries its verdict.

**Named methods result.**
- Animals replicate within cohorts, but each mission's exposure environment and matched control
  environment are each realized once.
- Condition and cohort are therefore inseparable within a mission. In OSD-771, group is also
  confounded with the ERCC spike-in batch.
- Cohort-level deviations are estimated from control-vs-control contrasts across missions. Do not write
  "n = 1".

**Excluded:**
- podocyte as hero;
- recovery or rebound;
- a days-since-return curve;
- rank-after-adjustment;
- mechanism stories;
- un-rerun post-hoc numbers;
- the ccRCC H&E project.

**Sequence:**
1. Hard freeze `v1.0-manuscript`: configs, hashes, figure source tables, environment lock, CI.
2. Internal hostile review (one quantitative and one biological reader).
3. A compact package to Dr. Casaletto: skeleton, figures, repository, and two questions (is this a
   contribution, and what would a reviewer attack first).
4. Preprint, then a journal: npj Microgravity first; Communications Biology if A1 and A5 land
   cleanly; Life Sciences in Space Research as the alternative.

**Drafting gate.** In weeks 1–2, draft only the Introduction, Methods and 4-axis Results. Write no
structural or podocyte prose until A1–A3 and A5 resolve.

## 10. Timeline (target, not a gate)

| Week | Work | Decision point |
|---|---|---|
| 1 | A0 blind gates and hash-lock; A9 sent; A5 claim enumeration, mapping audit and predictions; A6 and A8 predictions; code A1–A4 (no execution) | gates locked |
| 2 | Execute A1, A2, A3, A4 (separate executor); draft Introduction, Methods, 4-axis | **D1: H1 and H2 labels → §7.2 row** |
| 2–3 | A5 ledger; A6 heart (VST validation first) | **D2: attribution** |
| 3 | A7 (supplement) | — |
| 4 | A8: OSD-163 linkage and coverage gate, then effects; OSD-102 if informative | — |
| 4–6 | Results prose by §7.2 row; figures; supplement | — |
| 6–7 | Software freeze (environment, CI, Docker, README); `v1.0-manuscript`; one week of buffer for a gate rejection or a VST failure | hostile review → Casaletto → preprint (target December 2026) |

**Summed effort** is about 30–33 agent-days, plus R/VST setup, the ledger audit and the
environment/CI work. If time runs short, cut in this order:
1. brain (A6 secondary);
2. OSD-102 proteomics;
3. A7.

## 11. Compute and agent budget

- **Code and compute.** All analyses reuse existing loaders and drivers. The largest downloads are
  OSD-102 proteomics (3 GB, conditional) and GSE129798.
- **Agent use** (`docs/strategy/2026-10-07/README.md`):
  - at most 2–3 agents in flight;
  - commit every result before the next stage;
  - Sonnet for reading and well-specified implementation; Opus for judging and statistical design;
  - compact handoffs;
  - the gate-writer and the executor are fresh agents given only Appendix A.

## 12. Risks

| Risk | Mitigation |
|---|---|
| H1 fails for podocyte | Expected and disclosed (§5.3); rows R2–R4 are publishable |
| Gate-writer or executor contamination | Fresh agents; Appendix A only; verbatim adoption; hash-lock; full run log |
| τ̂²_cohort unstable (few clean pairs) | Common τ²; χ² CI; leave-one-contrast-out range; upper-bound sensitivity |
| Exchangeability assumption for missions without extra controls | Stated; sensitivity excluding OSD-102/163 |
| The ledger is subjective | Frozen rule; second-person audit; registered before scoring |
| Ledger read as replication | Stated as a re-analysis of the same missions under a different estimand |
| Heart VST differs from GeneLab | Numeric validation gate before use |
| Mediator/collider bias in sampling covariates | External rule-derived panels; control-only slopes; robust/not-robust reporting only |
| Scope creep | §5.9; §7.3; no new analyses unless they correct an error, answer a concrete reviewer objection, or form a separate project |

## 13. Housekeeping before the freeze

- Call CPM eligibility "arm-aware", not "label-blind" (`PROJECT_TECHNICAL_AUDIT_2026-08-25.md`
  line 595).
- Structural I² is 0.67%, not 67.3%.
- The structural control ranks 2nd, not 4th.
- "No sixth mission" is wrong (OSD-457).
- OSD-771 groups have overlapping dissection dates: do not write "one dissection day per group".
- Record the K0, G1/G2, item-9 and this plan's decisions in `docs/owner_decisions.md`.
- Refresh the README.
- Unify the Python and R environments and pin GitHub dependencies.
- Add CI.
- Retire the stale Docker entrypoint.
- Delete or rewrite the Casaletto email draft.

---

## Appendix A. Blind brief for the G3 gate-writer and executor (contains no observed values)

You are writing, or executing, registered sensitivity analyses of an existing cross-mission kidney
transcriptomics result. Read only this appendix and the code paths listed. Do not open any results
directory, the `docs/strategy/` folder, or any other document in `docs/`.

**Code paths you may read:**
- `src/clinical_axes/statistics.py`
- `src/clinical_axes/data.py`
- `src/clinical_axes/analysis.py`
- `src/clinical_axes/sensitivity.py`
- `scripts/clinical_axes/run_compartment_context.py`
- `scripts/clinical_axes/run_strict_podocyte_matching_audit.py`
- `scripts/clinical_axes/run_podocyte_sensitivities.py`
- `config/clinical_renal_axes_cross_mission.yaml`

**Design.**
- Five terminal spaceflight missions. The animal is the unit. Flight vs habitat ground control is the
  designed contrast.
- Some missions also have vivarium and/or basal control groups. One mission has two replicate
  ground-control cohorts (original and rerun).
- In each mission each group is one cohort, so cohort-level deviations are not estimable from animal
  replication. In one mission, group is also confounded with a processing (spike-in) batch.
- Effects are Hedges g per mission (strata combined within mission), pooled by REML/modified
  Hartung–Knapp.
- Max-T family-wise correction over a 49-set family via synchronized blocked label permutation.

**Decide and write, as YAML, each of the following.** Choose every threshold yourself; none is
suggested.

1. **Cohort variance.** The joint model for τ²_cohort, given that same-condition replicate cohorts are
   the clean estimator and other control contrasts include systematic differences. Specify the
   upper-bound sensitivity, how total variance combines with REML τ², whether τ² is common across sets,
   and the leave-one-contrast-out check.
2. **H1 label table.** The ratio threshold r* for alternative-control agreement; what "batch-matched
   substitution keeps sign" requires; the CI level.
3. **H2 matched-null gate.** The percentile threshold; whether both candidate pools are required; and
   whether a second reference set's percentile gets a gate.
4. **H4 robustness rule.** The sign and calibrated-retention thresholds per technical sensitivity.
5. **H5 ledger verdicts.** The CI level; the equivalence margin δ for ABSENT; and the multiplicity rule
   within the set-level and gene-level families.

**Constraints.**
- No gate may upgrade a frozen verdict.
- Each gate applies identically to every set and axis.
- Labels are "registered sensitivity (post hoc in origin)", or "preregistered secondary" for H5.
- Seeds are the parent seed plus new offsets.

**Output.** `config/clinical_axes_g3_registration.yaml`. It is adopted verbatim; it may be rejected
only for internal inconsistency, which triggers a new writer.

**Executor rule.** Run only the hash-locked code. Log every run. Open outputs only after completion.

## Appendix B. How this plan handles the owner's external assessment (2026-10-07)

| Owner recommendation | Disposition |
|---|---|
| Inventory GSE295428 for podocytes before more bulk work | **Done.** "Kidney, LAR" is a duplicate spleen upload, so the path is closed (§4.3, §8). Report the duplicate (A9). |
| Label-invariant (pooled-sample) eligibility | A4.1 |
| Independent second atlas | Deferred (§8); if run, glomerular references + GSE129798 rather than Tabula Muris Senis |
| Orthogonal single-sample scoring | A4.3 |
| RIN-adjusted LAR sensitivity | A4.4 (G2 label unchanged) |
| Lock the environment and add CI | §13, week 6–7 freeze |
| Pursue RRRM-1 bulk kidney mRNA | Outreach query (A9); nothing depends on a reply |
| GSE295428 animal-level pseudobulk validation | Closed until a corrected deposit (§8) |
| OSD-913 kidney miRNA | Not scheduled (no cell localization; RRRM-1 dissection-method confound) |
| Narrow proteomics comparison (OSD-102/163/462) | A8 for OSD-163, then OSD-102. **OSD-462 excluded**: condition is aliased with TMT reporter-tag block. |
| OSD-513 as a heterogeneity case | Descriptive moderator; never recovery evidence |
| Latent structural-axis or factor analysis | Optional, not scheduled |
| Duration meta-regression | Parked (§8) |
| Matched random gene sets | A3, extending K4 with coherence, GC/length, a genome-wide pool and calibrated variance |
| Decompose the structural components | Not scheduled |
| Freeze → hostile review → Casaletto → preprint → npj Microgravity | §9 |
| Stop list (P1 tuning, OSD-513 as recovery, podocyte injury, redundant sensitivities, p-threshold robustness, silent changes, ccRCC) | §2.3, §5, §9 |
| Title and careful cohort-replication wording | §9 |
| The negative ISS-T/LAR correlation as "suggestive of rebound" | **Superseded:** a cohort artifact (§2.2) |
| Go/no-go judged on direction, not p thresholds | Adopted (§7.3), with one addition: 32/49 sets shift positive, so direction alone cannot establish specificity. H1 calibration and H2 carry that weight. |
