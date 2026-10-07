# Adversarial review of RESEARCH_PLAN_2026-10

Reviewer stance: biostatistician and space-biology editor. Object: `docs/RESEARCH_PLAN_2026-10.md` v1.0, including Appendices A and B.

## Verdict

**Ready with fixes.** The structure, freeze discipline and parking decisions are sound. Five defects would make results ambiguous or the blinding theatrical: the A1 variance model, the undefined H1 object and non-exclusive outcome rows, the known A1/A2 outcomes, the A3 null that ignores set coherence, and A5 being labelled "unseen". None of the fixes needs new data.

## Required fixes

**1. A1 variance model double-counts heterogeneity and conflates systematic and random effects.** Independent cohort deviations across missions *are* between-mission heterogeneity. Adding τ²_g on top of REML τ̂² double-counts it whenever τ̂² > 0. The control-vs-control contrasts are also not pure noise:
- vivarium−ground contains a systematic housing effect, the same one A1 then uses as a control-choice arm;
- basal−ground contains about 2 months of age and time (OSD-771 basal animals were euthanized 22–23 Jul 2019; ground animals 18–20 Sep 2019).

Replace A1 "Estimand/Method" with:
> "Primary: total between-mission variance = max(τ̂²_REML, τ̂²_cohort). τ̂²_cohort comes from a joint REML model in which control-vs-control contrasts enter as extra observations that share the cohort variance component, with a random mission intercept and a covariance term for contrasts sharing a control group. Same-condition replicate cohorts (OSD-253 original vs rerun ground) are the clean estimator. Vivarium−ground and basal−ground are used only in an upper-bound sensitivity, with a fixed effect for contrast type. τ²_cohort is assumed exchangeable across missions and applied to missions without extra controls (OSD-102, OSD-163); this assumption is stated. Report the χ²-based CI of τ̂²_cohort and its effective degrees of freedom."

Also state whether τ² is set-specific or common. Recommended: common across sets (a set-specific τ² from about 5 effective pairs is not identifiable), with set-specific values as a sensitivity.

**2. Disclose the known A1/A2 outcomes.** The plan's own §2.2 already shows that podocyte H1 fails under the default gate:
- calibrated interval about (−0.31, 1.68);
- flight−vivarium/flight−ground ratio 0.28/0.689 = 0.41, below 0.5;
- OSD-771 substitution −0.17.

Under the plan's τ² of about 0.36, the minimum calibrated effect that can exclude 0 at k = 5 with mHK is about g ≥ 1.0. Add to §5:
> "Known-outcome disclosure: for the podocyte set, A1 and A2 under any reasonable gate are expected to fail, from post-hoc previews. The registration records this expectation and the minimum detectable calibrated effect before the gates are run. A1/A2 are informative for the structural control and the 4 axes, whose calibrated values have not been previewed."

If structural previews exist, list them too.

**3. Define the H1 object and make §7.1 exclusive and exhaustive.** H1 names two objects (podocyte and structural), but §7.2 uses a single pass/fail. Also, the structural control's *uncalibrated* CI already crosses 0 (−0.095, 1.193), so none of the H1 rows applies to it.

Replace the H1 rows with:

| Uncalibrated CI excludes 0 | Calibrated CI excludes 0 | Vivarium sign/ratio | Label |
|---|---|---|---|
| — | yes | same sign, ratio ≥ r* | spaceflight-associated (calibrated) |
| — | yes | fails | calibrated; control-dependent |
| yes | no | any | habitat-relative shift; not robust to cohort variation |
| no | no | any | not established (direction only) |

Add a separate overlay: "flight−vivarium CI excludes 0 opposite to flight−ground → habitat/housing effect". A point estimate of opposite sign is not enough for that label. Flight−basal is reported only. State that §7.2's H1 refers to the **structural control** for rows 2–3 and to the **podocyte set** for row 1. Add §7.2 row 4: "H1 fail, H2 pass → the podocyte position is extreme among matched sets but not calibration-robust; report it descriptively."

**4. Tie A2 into the decision.** Add to A1's gate:
> "H1 'pass' additionally requires A2 (substitution and drop) to keep the sign."

§7.3's damage criterion ("A1 and A2 and A4 all reverse") is too lenient and contradicts the A2 row ("major caveat"). Replace it with "any two of A1, A2, A4 reverse the sign of the structural control."

**5. Correct the OSD-771 batch claim.** The euthanasia dates in `s_OSD-771.txt` overlap across ISS-T groups:
- flight on ISS 16–19 Sep;
- ground 18–20 Sep;
- vivarium 18–20 Sep.

"Group equals ERCC mix equals dissection day" is therefore not supported for dissection day. The ERCC-mix confound is confirmed (DIRECTIONS_REVIEW spot-check 3). Rewrite §2.2 and A2 as "group equals ERCC mix; dissection days overlap." Then add a dissection-day covariate sensitivity to A2.

**6. A3 null must contain coherent sets.** Random genome-wide panels of non-co-expressed genes have set-score SDs near zero, so any coherent set looks extreme. Add "mean pairwise inter-gene correlation in control animals" as a matching variable in both pools. Specify the ranking statistic as calibrated t, not g. If τ² is common across sets (fix 1), calibration cannot change a ranking by g. Specify the panel size: 142 observable genes, as in K4, or the full 157. Add a gate for the structural control's percentile (≥ 95th → "band is special"), or state that it has none.

**7. A5 is not unseen and needs an equivalence margin.**
- `descriptive_signed_gene_effects.tsv` already holds per-mission effects for axis genes, including the ledger genes (Podxl/Nid1/Slc12a3/Fn1 matches), and the K4 audit states "neither PODXL nor NID1 alone was a significant cross-mission effect."
- Relabel A5 "unseen except single genes in the 4 axes (listed)".
- Delete the pre-written H5 row "PODXL/NID1 up vs Siew down … estimand difference": it is a conclusion written before scoring.
- Replace the verdicts with:
  - RECURS: calibrated CI excludes 0 in the claimed direction;
  - CONTRADICTED: it excludes 0 in the opposite direction;
  - ABSENT: the calibrated CI lies inside (−δ, δ), with δ = 0.3 g registered;
  - INDETERMINATE: otherwise.
- Single-gene claims use a gene-level τ² (not the set τ²) and are a separate family.
- Add multiplicity: report the expected number of false RECURS calls at α = 0.05 across the ledger, or Holm within the ledger.
- State whether Siew 2024 used overlapping OSDR missions. If it did, H5 is re-analysis under a different estimand, not independent replication.

**8. Split H6 and give A6 a real estimand.** H6 joins two opposing predictions. Split it:
- H6a: the structural shift is absent in heart;
- H6b: the control-group deviations are shared across organs.

For H6b the estimand is the correlation, across a registered panel of tissue-agnostic sets (KEGG structural plus housekeeping-adjacent sets, not kidney compartment sets), between kidney and heart deviation vectors for basal−ground and vivarium−ground. The null comes from gene-set rotation or from permuting animals within arm.

The per-animal within-group correlation tests animal-level, not cohort-level, covariation; label it so. "Not shared → tissue-processing batch" is not the only alternative (tissue-specific housing biology is another). Replace that row with "not shared → cohort artifact not organism-wide; cause unresolved."

VST validation gate: on the kidney mission, require per-gene Pearson ≥ 0.99 against GeneLab VST *and* identical 49-set mission g to ±0.02, before heart use. Multiplicity family: 3 handling panels × 2 contrasts, Holm; brain is descriptive only.

**9. A8 concordance rule.** With one mission and n = 6/6, a sign match occurs half the time by chance. Compare the protein effect with the **same mission's** RNA effect: OSD-163 RNA podocyte is not the pooled sign, and the mission effects differ. Require "protein CI excludes 0 in the RNA direction" for "agrees" and "protein CI excludes 0 opposite" for "disagrees". Anything else is "uninformative".

Define coverage as "quantified in ≥ 80% of samples". Confirm the sample-ID linkage before calling the data the "same animals". The "difference program" must be labelled descriptive and cannot count as H3/K0 evidence.

**10. Appendix A defects.**
- It says "the code paths named below", but names none. List them, or say "no code".
- Delete the §6 "default gates" (ratio ≥ 0.5, ≥ 95th percentile, retention ≥ 0.5). They were written by unblinded authors and are anchored to known values (0.41 sits just below 0.5).
- Add to §5:
  > "The blind writer's YAML is adopted verbatim. The owner may reject it only for internal inconsistency, in writing, before any run; a rejection triggers a new blind writer, not an edit."
- Add: "A1–A4 are executed by an agent that has not seen §2.2. The output is opened only after the hash-locked run completes, and all runs are logged, including failed ones."

## Recommended improvements

- **H3 back doors.** State explicitly that A3 and A8 cannot reopen separability.
- **H9.** H9 has no test in this plan; relabel it "outreach question" rather than a hypothesis.
- **A4.** Restrict the H4 verdict to the podocyte set and the structural control; with 54 sets × 5 items × 2 maps, "any item fails" is near-certain. Define "retention" as calibrated-g ratio. Note that A4.4 (RIN, LAR) does not touch the primary ISS-T meta.
- **Gene-map rule (§5.6).** Apply it to A3 and A5 as well, and say which map's label is reported when they disagree (the more conservative one).
- **§4 consistency.** OSD-457 is listed "to acquire", yet §2.2 reports a descriptive g of 1.16; state where that came from. OSD-580 files carry a GLDS-573 prefix and OSD-561 files a GLDS-556 prefix. That is plausible on OSDR, but note it so a loader does not reject them. Dataset designs (RR-1/3/7/10, RRRM-2, MHU-3) look correct; the 78 heart animals (38 ISS-T + 40 LAR) are consistent with the kidney design.
- **Missing analysis.** Add a calibration-validity check, leave-one-control-contrast-out τ̂², and report its range.
- **Timeline.** Summed effort is about 30–33 agent-days, plus the R/VST setup, the second-person ledger audit, and environment/CI/Docker work. The plan has no buffer for a gate rejection. Week 1 says "start A4", but gates must lock before any A1–A4 run; move A4 execution to week 2, so only coding happens in week 1. Cut first if time runs short:
  1. brain (A6 secondary);
  2. OSD-102 proteomics;
  3. A7.

  Move "December preprint" to a target, not a gate.
- **Scope creep.** Add to §12: "A9 replies, GEO corrections or new deposits arriving before the freeze are logged, not analysed."

## Spot-checks performed

1. `src/clinical_axes/data.py::cpm_eligible_genes` uses the FLT and GC sample sets separately ("either contrasted arm"). The "arm-aware" description is **confirmed**.
2. K4 section of the audit: 142-gene panels, 10,000 panels, p = 0.0001, target-minus-matched 0.54 (−0.17, 1.25), blocked p 0.0148, 81/133 endothelial (among the one-to-one *trimmed matches*, not the random-panel pool). **Confirmed**, with that wording caveat. The same section states that PODXL/NID1 single-gene results were already seen.
3. `data/raw/metadata/s_OSD-771.txt`, group × euthanasia site × date: basal 22–23 Jul; ground 18–20 Sep; vivarium 18–20 Sep; flight on ISS 16–19 Sep. Dissection-day confounding is **not supported**. ERCC-mix confounding is confirmed via DIRECTIONS_REVIEW spot-check 3.
4. `data/results/clinical-axes/20261006T004825Z_g1g2/descriptive_signed_gene_effects.tsv` has per-mission gene effects for the axis genes (10 rows matching Podxl/Nid1/Slc12a3/Fn1). A5 is **partly seen**.
5. `docs/strategy/2026-10-07/DIRECTIONS_REVIEW.md` confirms the ERCC-mix pattern (ISS-T: FLT/VIV Mix 1, GC/BSL Mix 2; reversed in LAR) and the G2 values cited in §2.1.
