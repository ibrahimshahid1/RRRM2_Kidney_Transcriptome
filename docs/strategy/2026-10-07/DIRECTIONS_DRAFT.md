# Next directions: coordinator draft (for adversarial review)

**Status:** draft. It goes to one adversarial judge before it is finalized.
**Inputs:**
- `reader_history.md`, `reader_leads.md`, `reader_external_data.md` and `reader_priorart.md`
- the G1/G2 results doc
- the K0 revision conclusion

## 1. Where the project actually stands

**Live workstreams.** Of about 28 workstreams, one is alive: the cross-mission clinical-axes analysis
(5 terminal missions, 83 animals). Its defensible content:

- **The 4 clinical axes do not recur.** Three of four are no-go. The glomerular barrier axis moves
  *opposite* to the declared "identity loss" (g = −0.716).
- **A structural band shifts upward, with podocyte HS at the top.**
  - Podocyte HS g = 0.689 (0.042, 1.336); structural control 0.549.
  - The direction and magnitude are robust (G1, 8/8 atlas rebuilds, both gene maps), but the
    family-wise margin is thin.
  - The podocyte shift is not separable from structural drift (K0).
- **Recovery is not established.** G2 is INCONCLUSIVE, as powered.

**What the 2026-10-06 readers added (post hoc, unregistered; must be preregistered before use):**

1. **The main remaining threat is control choice.**
   - Against vivarium controls the podocyte meta effect is 0.28 (−0.45, 1.01).
   - Control-vs-control contrasts show about 0.36 g² of excess group-level variance. A calibrated
     interval would cross zero (−0.31, 1.68).
   - Structural is stronger than podocyte against vivarium (0.55), so the podocyte-over-structural
     ordering is control-dependent.
2. **Technical and sampling correction may sharpen K0 rather than weaken it.**
   - Under GC/length/expression correction and prespecified-in-spirit sampling-marker covariates
     (cortex, medulla, blood, hilar), podocyte stays rank 1 (0.198 → 0.165 score units, CI above 0).
   - Under the same correction structural, endothelial and mesenchymal mostly collapse (structural
     0.549 → 0.187 g).
   - The panels were chosen after seeing the data, and the covariates may be flight-responsive.
3. **The G2 "rebound" (r = −0.65) is an artifact.** The same anti-correlation appears between the two
   ground-control groups without any flight animals. The ERCC spike-in mix is confounded with
   condition within each arm and reversed between arms, and each group was dissected on its own day.
4. **The exploratory G2 interaction hits are on-orbit handling and sampling signatures, not biology.**
   The ISS-T flight group carries them (immediate-early genes 2.14, medulla 2.71, blood 1.73), and they
   explain the endothelial, mesenchymal, TAL and fibrosis passes.
5. **An unscored sixth kidney cohort exists: OSD-457 (JAXA MHU-3).**
   - Male C57BL/6J, 31 d flight, sampled 2 days after return, 3+3 wild-type and 3+3 Nrf2-KO.
   - Descriptively the podocyte HS g is 1.16 (wild-type 1.36), rank 4 of 49, and its 49-set pattern
     correlates r = 0.61 with the pooled terminal pattern.
   - The 2026-08-25 audit's "no sixth mission exists" is wrong.
6. **No public terminal arm matches LAR's 32-day flight.** RRRM-1 kidney data are miRNA-only, or
   single-cell with labels reset to "unknown". Bion-M2 and Siew's GCR kidney data are not public.

## 2. Ranked directions

### Tier 1: decides the paper (start now, about 4–6 weeks)

**D1. G3: a preregistered threat package, run before any drafting is finalized.** These are the checks
most likely to change the headline:

- **(a) Control choice.** Flight vs vivarium and flight vs baseline, wherever those groups exist, plus
  an empirical group-variance calibration from all control-vs-control contrasts, reported as a
  calibrated interval.
- **(b) Technical normalization.** CQN/EDASeq GC-length normalization.
- **(c) Sampling composition.** Prespecified tissue-sampling marker covariates. Apply them only to sets
  whose markers are not in the panels, to avoid circularity.
- **(d) Group-level batch flags.** ERCC mix and dissection day.

Gates are frozen before running, and the provenance section must disclose that the readers' post-hoc
numbers were seen before the lock.

Outcomes:
- If (a) holds, the headline survives as "spaceflight-associated".
- If (a) fails, the paper becomes "recurrence across missions is not robust to control choice". That is
  still a reproducibility result.
- If (b) and (c) show podocyte surviving while structural collapses, the K0 framing is *sharpened* to a
  podocyte-leaning signal that survives technical and sampling correction. It is never upgraded to
  "separable" unless K0's own score-unit interval clears zero.

Effort about 2 weeks. Code reuse is high (the G1 driver and the loaders).

**D2. Add OSD-457 (MHU-3) as a frozen descriptive moderator cohort, plus a preregistered
days-since-return grid.**
- Ingest it through the project loaders, with GeneLab counts and the same eligibility and z rules.
- Grid: ISS-T 0 d; OSD-513 about 1 d; MHU-3 2 d; LAR 24 d. Flight lengths are 31–38 d for the
  live-return cohorts against LAR's 32 d, which partly decouples recovery from exposure length.
- This is descriptive only (n = 3+3 wild-type) and crosses strain, sex and lab.
- Effort about 1 week. It corrects the audit's "no sixth mission".

**D3. Write the cross-mission reproducibility paper (the critical path; start drafting in parallel
with D1).**
- **Working title:** "Most reported renal transcriptomic responses to spaceflight do not recur across
  independent mouse missions."
- **Results:** 4-axis non-recurrence; the structural band with a podocyte-leaning top; G1/K0/G2 (and
  G3) as the robustness table.
- **Methods contribution:** preregistered locks, adversarial audits, the reconstructed gene map, the
  Freedman–Lane correction, and a recovery design that states its own aliasing.
- **Target:** bioRxiv by December 2026, then a journal (npj Microgravity or Communications Biology
  tier).
- **Before submission:** tag the freeze, refresh the stale README, technical audit and owner_decisions
  (structural I² is 0.67%, not 67.3%), and retire the Casaletto email draft.

### Tier 2: cheap, high information (alongside or after Tier 1)

**D4. Same-animal cross-tissue control.** Score RRRM-2 heart (OSD-580), cerebellum (OSD-561),
hippocampus (OSD-562) and brain (OSD-613) bulk RNA from the same animals. These accessions were
verified on 2026-10-07; OSD-580/561/562 carry the euthanasia-location factor, i.e. the ISS-T vs LAR
split. Score the structural program and handling panels (immediate-early genes, blood). Three tests:
- If the structural drift appears in other organs, it is systemic or handling-related.
- If it is kidney-specific, kidney attribution is stronger.
- Whether the ground-control-group deviation behind the G2 "rebound" artifact (BSL − GC and VIV − GC)
  recurs across organs of the same mice. If it does, it is a group-level batch effect, not kidney
  biology.

Effort about 1 week.

**D5. Protein layer: OSD-102 and OSD-163 kidney proteomics** (never touched). Check NPHS1, NPHS2,
PODXL, SYNPO and the structural set at protein level.
- This is the same animals, a second assay, and p-values must not be combined.
- OSD-163 is 21 MB and OSD-102 is 3 GB.

Effort about 1 week.

**D6. Item 9, re-scoped.**
- Use podocyte-rich glomerular references (GSE146912 controls; GSE111107) paired with a whole-kidney
  reference (Ransick GSE129798) instead of Tabula Muris Senis, which has few podocytes and mislabels
  them.
- Mandatory only if "podocyte" stays in the title.

### Tier 3: long-lead and outreach (start the asks now; they gate nothing)

**D7. Morphometry request.** This is the only way to separate composition from per-cell expression.
- Ask the PI, the Siew consortium and NASA LSDA/ALSDA biospecimens for: WT1+ podocytes per glomerulus,
  glomerular density and tuft area, cortical sampling fraction, DCT density, and ECM/cytoskeletal
  staining.
- Target sections: RR-10 and RRRM-2.

**D8. Ground analog GSE235042** (hindlimb-unloading kidney, which reports glomerular widening). Asks
whether unloading alone moves the scores. Requantify from SRA; n = 3 per group; descriptive.

**D9. Watchlist.** Corrected labels for the RRRM-1 kidney single-cell data (GSE295428), Bion-M2 (2025)
kidney deposits, and Siew GCR kidney data (OSD-706..712).

## 3. Stop or retire

- The G2 "rebound" interpretation: it is a ground-control-group and ERCC artifact.
- The exploratory G2 interaction hits as biology: they are on-orbit handling and sampling signatures.
- Any recovery claim. No public matched 32-day terminal arm exists.
- Podocyte as headline or title. Gene-level storytelling (glycocalyx, sialylation) unless sub-modules
  are preregistered.
- Already-closed lines: OSD-462 phospho mining (stop rule), network rewiring, Grey60 rescue, LAR
  reversal framing, the OSD-462 data note, the matched-library paper and the merged methods paper.
- Sending the existing Casaletto email draft: it pitches retracted results.

## 4. Sequencing (about 6 weeks)

| Week | Work |
|---|---|
| 1 | D2 ingestion (OSD-457); write and lock the G3 preregistration (D1); send the D7 morphometry ask |
| 1–3 | Run D1; draft D3 (Introduction, Methods, the 4-axis Results) |
| 3–5 | D4 and D5; revise D3 Results around the G3 outcome |
| 6 | Finalize D3 figures and supplement; D6 only if "podocyte" is in the title |

## 5. Risks to the plan

- **G3 control choice could remove the "spaceflight" attribution.** That is the point of running it
  first. The paper survives as a reproducibility result either way.
- **G3 is preregistered after the readers' post-hoc numbers were seen.** That must be disclosed. The
  sampling panels need a justification independent of these data, such as published marker lists.
- **Every external addition is underpowered.** OSD-457, GSE235042 and the duration grid are
  descriptive.
- **Novelty risk for a "systemic structural/ECM drift" framing.** See `reader_priorart.md`.
