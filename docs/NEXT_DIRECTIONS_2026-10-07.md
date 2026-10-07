# Next directions after the G1/G2 results

> **Superseded for planning by `docs/RESEARCH_PLAN_2026-10.md`** (v1.1, 2026-10-07). Kept as the decision
> record that led to it.

**Date:** 2026-10-07
**Status:** recommendation for the project owner
**Basis:**
- five research readers' notes, the coordinator draft and one adversarial review, all in
  `docs/strategy/2026-10-07/`;
- the G1/G2 results (`docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md`);
- Siew et al. 2024, read in full text.

## The situation in one paragraph

Of about 28 workstreams, one is alive: the cross-mission clinical-axes study (5 terminal missions, 83
animals).

**What it shows:**
- The four clinical axes do not recur. The barrier axis even runs opposite to "identity loss".
- A low-heterogeneity structural band shifts upward. The 157-gene podocyte set is at the top
  (g = 0.689), and it is the only one of 49 sets whose interval clears zero.
- The direction is robust (G1). It is not separable from the structural drift (K0). Recovery is
  unresolved (G2).

**What the readers and reviewer added:**
- **The significance margin is thin.** Bonferroni over the 53 tests actually examined gives about
  p 0.05, and dropping the top 5 genes takes the lower bound below zero.
- **The attribution to spaceflight is exposed.** Each group is one cohort, and in OSD-771 group is also
  confounded with the ERCC spike-in batch (dissection dates partly overlap). Against vivarium controls,
  the podocyte effect is 0.28 (−0.45, 1.01).
- **The G2 "rebound" and the exploratory recovery hits are artifacts**, of group batch and on-orbit
  handling respectively.
- **Siew 2024 already reports qualitative nephron structural remodeling**, and its PODXL/NID1 vote
  count runs opposite in direction to our podocyte estimate.
- **The new contribution is therefore measurement:** a preregistered, calibrated estimate of which
  reported renal responses actually recur across independent missions.

## Recommended directions, in order

**1. Make the next tests non-theatrical (week 1).**
- **Blind gate-writing.** An agent or person who has not seen `reader_leads.md` writes the G3 gates and
  decision table. The config and code are hash-locked through the rrrm2 CLI before anything runs.
- **Labelling.** G3 is labelled "registered sensitivity analysis (post hoc in origin)", never
  "confirmatory".
- **Directional predictions** for the genuinely unseen data (steps 4 and 6) are written before download.
- **Morphometry requests go out now**: WT1+ podocytes per glomerulus, glomerular density and tuft area,
  cortical fraction, DCT density. Send them to the PI, the Siew consortium and NASA LSDA biospecimens.
  They gate nothing, but they are the only route to composition vs per-cell expression.

**2. G3 core: the tests most likely to change the headline (weeks 1–2).**
- **Group-variance-calibrated control choice, across all 49 sets.**
  - The primary contrast stays flight vs habitat ground control.
  - Each mission's variance is inflated by a group-level τ², estimated with a dependence-aware
    estimator (not the raw 35 overlapping contrasts).
  - Frozen label rule: "spaceflight-associated" requires the calibrated 95% CI to exclude 0 **and**
    flight-vs-vivarium to agree in sign with a ratio of at least 0.5. Otherwise the label is "shift
    relative to habitat controls; not robust to control choice".
- **OSD-771 batch-matched substitution.** Flight-vs-vivarium is the only within-batch-matched flight
  contrast there (podocyte −0.17). Report the 5-mission meta with it substituted, and with OSD-771
  dropped.
- **Matched-random-gene-set null.** Use about 10,000 sets matched on size, expression, GC and length,
  scored through the same pipeline, and report the podocyte percentile. This is the cleanest
  "specificity" test.
- **Compact technical sensitivities, reported as robust or not robust only:** GC/length normalization,
  and sampling covariates whose panels come from an external atlas by rule, with slopes fitted in
  control animals only.
- **Expected outcome, stated now:** the "spaceflight-associated" label probably fails the calibrated
  gate. The paper survives either way, because its claim becomes recurrence measured honestly.

**3. Prior-claim recurrence ledger (weeks 2–3).**
- Enumerate the field's specific renal transcriptomic claims before scoring: Siew 2024's pathways and
  most-frequent genes (PODXL, NID1, transporters, focal adhesion, tight junction, actin), Suzuki 2022
  (MHU-3), and DCT remodeling.
- Map each claim to a gene set by a frozen rule and score it with the G1 machinery and calibrated
  variance.
- Give each claim a verdict: recurs, does not recur, or not evaluable.
- This gives the paper its denominator and turns the Siew conflict into a tested estimand difference.

**4. Same-animal cross-tissue control (weeks 2–3; genuinely unseen data).** RRRM-2 heart (OSD-580),
cerebellum (OSD-561), hippocampus (OSD-562) and brain (OSD-613) come from the same mice, and the first
three carry the ISS-T/LAR split. Three questions:
- Does the structural drift appear outside the kidney (systemic or handling) or not (kidney-specific)?
- Do the handling panels (immediate-early genes, blood) move the same way?
- Does the ground-control-group deviation behind the fake "rebound" recur across organs (group batch)?

**5. Write the paper, gated (weeks 1–6).**
- **Working title:** "Which reported renal transcriptomic responses to spaceflight recur across
  independent mouse missions?"
- **Weeks 1–2:** draft only the Introduction, Methods and the 4-axis non-recurrence. Write no podocyte
  or structural prose until steps 2–4 resolve.
- **Main figure:** a recurrence map. Rows are the prior claims, the 4 axes, and the structural and
  podocyte sets. Columns are the missions plus a calibrated meta column and I². Each row carries its
  verdict.
- **Named methods result:** animals provide biological replication within each cohort, but each
  mission's exposure environment and matched control environment are each realized once (one cohort
  per group; in OSD-771 also one ERCC spike-in batch per group). Cohort-level handling, habitat and batch
  deviations therefore cannot be estimated from within-mission animal replication. The project
  estimates them from control-versus-control contrasts across missions. This calibration is the most
  transferable contribution to the field. Do not phrase it as "n = 1 vs 1": the animal-level contrast
  does have within-study replication.
- **Disclose the provenance chain:** targeted axes → compartment scan → podocyte → G1/K0/G2 → reader
  post-hoc checks → G3.
- **Venue:** npj Microgravity is realistic. Communications Biology if the ledger and the calibration
  land cleanly.
- **Keep out of it:**
  - podocyte as the headline;
  - recovery or "rebound";
  - a days-since-return curve;
  - rank-after-adjustment;
  - glycocalyx mechanism;
  - any un-rerun post-hoc number.

**6. Smaller additions (weeks 3–4).**
- **OSD-457 (JAXA MHU-3).** A supplementary moderator row only: wild-type primary, Nrf2-KO separate,
  2–3 days of work. It corrects the August audit's "no sixth mission".
- **Kidney proteomics.** Run OSD-163 (21 MB) first, behind a protein-coverage gate. Run OSD-102 (3 GB)
  only if OSD-163 passes. These are the same animals in a second assay, so p-values are never combined.

**7. Watch, defer or drop.**
- **Watch:** a corrected deposit of the RRRM-1 live-return kidney single-cell data, and Bion-M2 and
  Siew GCR kidney deposits. GSE295428's "Kidney, LAR" samples were tested on 2026-10-07 and are
  duplicated spleen libraries; see `docs/strategy/2026-10-07/GSE295428_FEASIBILITY.md`. Report the
  duplicate upload to GEO and the authors.
- **Defer:** the second atlas (item 9), unless a reviewer asks; G1 already has 8/8 rebuilds.
- **Drop:** GSE235042, the hindlimb-unloading kidney study (n = 3 per group, FPKM only, sex unknown).

## Stop doing

- Recovery or rebound interpretations. No public 32-day terminal arm exists, and the rebound is a batch
  artifact.
- The exploratory G2 hits as biology. They are on-orbit handling signatures.
- Podocyte as the title; mechanism stories built from gene lists; "rank 1 after adjustment" as a claim.
- Reopening closed lines: OSD-462 phospho mining, network rewiring, Grey60 rescue, the OSD-462 data
  note, the matched-library paper and the merged methods paper.
- The existing Casaletto email draft. It pitches retracted results.

## Housekeeping before submission

- Tag the freeze and commit the run manifests.
- Refresh the README, the technical audit and owner_decisions:
  - the CPM eligibility rule uses FLT/GC membership, so call it "arm-aware", not "label-blind"
    (`PROJECT_TECHNICAL_AUDIT_2026-08-25.md` line 595);
  - structural I² is 0.67%, not 67.3%;
  - the structural control ranks 2nd, not 4th;
  - "no sixth mission" is wrong;
  - record the K0, G1/G2 and item-9 decisions.

## How the agent work is run from now on

See `docs/strategy/2026-10-07/README.md`.
- Commit every agent result before the next stage.
- Keep at most 2–3 agents in flight.
- Use Sonnet for reading and Opus only for judging.
- Pass compact digests plus file paths downstream.
- Give each agent a hard scope and a 1-page handoff.

This round of strategy work used about 0.2M subagent tokens across 2 agents. The first attempt used
1.3M across 10 agents and lost 7 of them.
