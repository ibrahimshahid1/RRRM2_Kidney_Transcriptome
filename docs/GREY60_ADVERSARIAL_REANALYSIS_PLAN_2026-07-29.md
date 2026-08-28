# Grey60 adversarial go/no reanalysis

**Authorized:** 2026-07-29  
**Status:** locked before the main adversarial rerun  
**Configuration:** `config/grey60_adversarial_reanalysis.yaml`

## Scope and disclosure

This analysis asks whether the historical 48-gene Grey60 result can support a
standalone paper. It is a prospective lock for the tests below, not a
preregistration that predates Grey60's discovery. The Grey60 genes, the
OSD-771 stratum effects, the old external projections, the preservation tables,
and a preliminary OSD-462 library-discordance screen were already known when
this document was written.

The primary proposition is deliberately narrow:

> In OSD-771 bulk kidney, the frozen 48-gene Grey60 score shows an
> ISS-terminal-associated FLT-minus-GC shift across ages that is not explained
> by a few genes or mice, technical or bulk-composition covariates, or merely
> the larger GC-derived green module, and recurs directionally in compatible
> terminal kidney cohorts.

The analysis cannot establish causal microgravity activation, a cell-intrinsic
response, kidney specificity, TLR4 causality, or network rewiring. OSD-771 arm
is inseparable from mission duration, live return/recovery, euthanasia
location, preservation, and dissection context. The correct exposure label is
therefore **ISS-terminal-associated**, not simply **spaceflight-induced**.

## Frozen object and inputs

- Membership is the 48 Ensembl IDs in
  `data/results/run_20260505_remediated_2500g/wgcna/module_gene_lists/grey60.txt`.
- Membership SHA-256:
  `cfefe6ff93204aa3a9607dd8cf8456a7f07b8980bf1b8ef51d731f074445a685`.
- Primary internal expression is the technically and composition-residualized
  Rtech matrix from `run_20260408_191759_2500g`.
- The primary score is the unweighted mean of per-gene z scores across the 48
  frozen genes. Its orientation is fixed to agree with saved ME17.
- Sensitivity scores are median-z, fixed historical kME-weighted mean-z, and
  saved ME17.
- OSD-771 inference uses one biological kidney per animal. External studies
  use the first NASA-designated technical replicate.
- External biological units are animals and cohorts, never genes.

## Gate A — internal, selection-adjusted evidence

1. Fit `score ~ Flight * Age * Arm` in the 40 OSD-771 FLT/GC mice.
2. Permute FLT/GC labels within the four Age-by-Arm strata, preserving 5/5
   allocation, for 100,000 Monte Carlo assignments.
3. On every assignment, refit all 20 non-grey saved module eigengenes and
   evaluate all seven factorial terms plus four stratum-specific flight
   contrasts. The null maximum is therefore over 220 inspected tests.
4. Grey60's expected-direction family-wise maxT p value must be at most 0.05.
5. A 10,000-replicate stratified animal bootstrap must place the pooled ISS-T
   FLT-minus-GC 95% interval entirely above zero. Both age-specific ISS-T point
   estimates must be positive.

Failure retires Grey60 as a post-selection discovery.

## Gate B — flight-blind recovery and compactness

1. Rebuild signed bicor WGCNA on BSL+VIV only (40 mice), selecting variable
   genes inside that reference only. Match the resulting modules to Grey60 by
   overlap without using FLT/GC labels.
2. Flight-blind recovery requires at least 24/48 genes, overlap BH q at most
   0.05, Jaccard at least 0.30, and a positive projected ISS-T effect with
   blocked-permutation p at most 0.05.
3. Audit the existing GC-only network by overlap rather than module color.
   The known 20.57 preservation statistic belongs to a 439-gene GC module
   containing 40 Grey60 genes; it is not compact-Grey60 preservation. The
   combined-network Grey60 statistic is 7.80.
4. Compare the Grey60 subset inside that GC module against the remainder and
   10,000 reference-mean/variance-matched subsets of equal size. Grey60 must
   exceed the remainder and the 95th matched-subset percentile.
5. Run the 27-member WGCNA grid:
   max genes {2,500, 5,000, 10,000}, minimum module size {20, 30, 50}, and
   merge height {0.15, 0.25, 0.35}. Recovery/co-clustering must pass in at
   least 18/27 specifications and the projected ISS-T effect must be positive
   in at least 22/27.

Failure means there is no compact-module paper, even if a broad ECM program
shifts.

## Gate C — animal/gene influence and coherence

- Every leave-one-mouse-out and leave-one-gene-out pooled ISS-T effect remains
  positive.
- Mean-z, median-z, kME-weighted mean-z, and ME17 agree in sign.
- Removing the top 1, 5, and 10 kME genes preserves a positive point estimate.
- Removing Grey60 genes overlapping independent endothelial/fibroblast marker
  panels preserves a positive point estimate.
- At least 32/48 gene-wise ISS-T contrasts are positive, with exact sign-test
  p at most 0.05.
- No gene contributes more than 20% of the sum of absolute gene contrasts.

## Gate D — technical, composition, and baseline falsification

1. Recompute the score from local VST data.
2. Fit a parsimonious technical model using library batch and scaled log depth,
   rRNA contamination, and RNA QA; then add composition PC1-PC2 derived from
   the ten CLR proportions.
3. Use HC3 uncertainty plus blocked Freedman-Lane permutation.
4. Rtech and VST effects must agree in direction; the composition-adjusted
   ISS-T coefficient must be positive, retain at least 50% of the unadjusted
   magnitude, have VIF below 5, and keep its direction when each nuisance
   covariate is omitted.
5. In BSL+VIV samples, an analogous arm-related pattern must not reproduce a
   large baseline artifact: maxT p must exceed 0.05 and absolute Hedges g must
   be below 0.5.

## Gate E — compatible external recurrence

- Terminal-study family: OSD-102, OSD-163, OSD-253 C57BL/6J, and OSD-462 WT.
- OSD-253 is synthesized within study across duration. The same-run original
  GC is primary but has a light-spectrum alias; white-light rerun GC is a
  sensitivity and is never double-counted.
- OSD-462's primary cohort score averages the three preparation-specific
  standardized scores within each animal. Preparation-by-flight interaction
  and each preparation are reported; no preparation may be selected because
  of its observed sign.
- OSD-513 is live-return, is a recovery-context moderator, and cannot rescue a
  failed terminal analysis.
- Report Hedges g and sampling variance per study, random-effects REML with
  Hartung-Knapp uncertainty, I-squared, study weights, and leave-one-study-out
  synthesis. DerSimonian-Laird is sensitivity only.

External recurrence passes only if at least 3/4 terminal-study point estimates
are positive, the pooled 95% interval is entirely above zero, I-squared is
below 60%, every leave-one-study pooled point estimate is positive, and no
study contributes more than 50% of the weight.

Failure retires Grey60 as a recurrent biological module paper. A positive
OSD-513 result alone is explicitly insufficient.

## Gate F — generic ECM/injury specificity

Compare Grey60 against frozen endothelial, fibroblast, podocyte, nephron,
ECM-remodeling, fibrosis, angiogenesis, EMT, TGF-beta, hypoxia, and wound
healing signatures. The Grey60 context effect must remain positive and retain
at least 50% of its magnitude after independent endothelial/fibroblast score
adjustment, and it must exceed the 95th percentile of expression/variance
matched 48-gene sets.

Failure does not change Gate E, but restricts wording to a generic
vascular/interstitial ECM or composition-associated state. Kidney specificity
cannot be tested with the kidney-only local cohort collection.

## Overall decision

Grey60 is a **GO** only if Gates A-E all pass. Gate F determines whether the
word *distinct* is justified. There are no rescue databases or additional
post-hoc signatures after a failed gate.

If Grey60 fails, the repository-wide alternative-direction audit ranks
directions by identifiability, independence, novelty, and local data
availability rather than by nominal p value.

## Execution audit amendment and final status

The completed analysis produced a NO-GO decision. Before the final external
result was issued, code audit identified and corrected:

- OSD-462 technical-replicate suffix parsing;
- scenario-dependent OSD-253 score scaling, replaced by one fixed
  cross-scenario score reference;
- unmodified small-study Hartung-Knapp scaling, replaced by the conservative
  modified form for the primary meta-analysis;
- silent signature-coverage and sample-count assumptions, replaced by
  fail-closed assertions.

These are conservative implementation corrections, not new biological
endpoints. They did not alter the declared cohort family, frozen 48-gene
membership, primary recurrence question, or all-gates-required decision rule.
The final results and interpretation are in
`docs/GREY60_GO_NO_GO_REPORT_2026-07-29.md`.
