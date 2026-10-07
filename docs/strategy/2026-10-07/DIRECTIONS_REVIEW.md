# Adversarial review of the next-directions draft

Reviewer lenses: biostatistician/methodologist; senior space-biology/nephrology editor.
Inputs read: `DIRECTIONS_DRAFT.md`, `reader_priorart.md` (with addendum), `reader_leads.md`,
`reader_external_data.md`, and the Summary and §5 of the 2026-10-06 G1/G2 results doc. Five spot-checks
(listed at the end). No web searches.

## Verdict

| Direction | Verdict | One-line reason |
|---|---|---|
| D1 G3 threat package | **MODIFY** | (a) is the decisive test, but the gates aren't defined. Drop the "K0 sharpened" outcome branch: it rewards a post hoc, mediator-prone adjustment. |
| D2 OSD-457 | **MODIFY (demote, shrink)** | Already scored, so it is not a test. Keep it as a 2–3 day supplementary moderator row (WT primary). No main-text duration grid. |
| D3 paper | **MODIFY** | "Most reported responses do not recur" needs a denominator. Four in-house axes are not "the field's reported responses". |
| D4 same-animal cross-tissue | **KEEP, PROMOTE to Tier 1** | It is the only remaining *unseen* data that can falsify the handling/batch explanation, and it is cheap. |
| D5 proteomics | **KEEP, gated** | Run OSD-163 (21 MB) first, with a coverage gate. Run OSD-102 (3 GB) only if OSD-163 passes. |
| D6 second atlas | **DEFER / conditional drop** | The stop-list already removes podocyte from the title, so D6's own trigger never fires. G1 has 8/8 LOSO rebuilds. |
| D7 morphometry asks | **KEEP** | It is the only route to composition vs per-cell expression. Send in week 1 and gate nothing on it. |
| D8 GSE235042 | **DROP** | n=3/group, FPKM only, sex unknown, and requantifying from SRA has a real cost. No outcome would change a decision. |
| D9 watchlist | **KEEP** | Zero cost. |
| Stop-list | **KEEP, add two items** | Add: "podocyte rank 1 after sampling adjustment" as a claim, and any main-text days-since-return curve. |

## Major issues

**1. "Preregistration" after the numbers were seen is theatre unless the gate-setter is blind and
the claims are downgraded.** The readers have already computed (a), (b) and (c). Locking them now
fixes the code, not the inference. Fix:
- (i) Have the gates and decision table written by someone, or an agent, who has *not* read
  `reader_leads.md`. Hash the config and code via the rrrm2 CLI before execution, and commit the
  hash with a timestamp.
- (ii) Label G3 "registered sensitivity analysis (post hoc in origin)" everywhere, never
  "confirmatory".
- (iii) Put the real confirmatory weight on data nobody has looked at: D4 other organs and D5
  proteomics. Write *directional predictions* for those before download.
- (iv) Add a provenance ledger to the paper: targeted axes → 3/4 no-go → 49-set compartment scan →
  podocyte → G1/K0/G2 → reader post hoc checks → G3.

**2. D1(a) has no decision gate, and the likely outcome should be stated now.** The reader's
calibrated interval is already (−0.31, 1.68) with p ≈ 0.055. FLT−VIV is 0.28 (−0.45, 1.01). The
honest prior is that the "spaceflight-associated" label fails. Proposed frozen gate:
- **Primary estimand stays FLT−habitat GC.** That is the designed contrast. VIV differs in hardware,
  so FLT−VIV is flight plus habitat and cannot serve as the arbiter by itself.
- **Calibrated interval.** Inflate each mission's variance by a group-level τ²_g. Estimate τ²_g with a
  dependence-aware estimator: one contrast per independent control pair, or a multilevel model with
  mission as a random effect. Do not use the raw 35 contrasts, 24 of which come from OSD-253.
- **Label rule.** "Spaceflight-associated (robust to group-level variance)" requires the calibrated
  95% CI to exclude 0 **and** the FLT−VIV meta to have the same sign with ratio ≥ 0.5 of FLT−GC.
  Otherwise the label is "shift relative to habitat controls; not robust to control choice or
  group-level variance".
- **Apply to all 49 sets,** so the rank of podocyte under calibration is reported, not cherry-picked.

**3. OSD-771's primary contrast is batch-confounded, and the draft understates this.** Spot-check 3
below: in OSD-771 the ERCC mix is FLT/VIV = Mix 1 and GC/BSL = Mix 2 in ISS-T, reversed in LAR. FLT−GC
is therefore mix-confounded in both arms, and **FLT−VIV is the only mix-matched flight contrast**.
There, ISS-T podocyte is −0.17. OSD-771 also carries the largest structural effect (1.35). Fix: add a
frozen G3 item that substitutes FLT−VIV for OSD-771 only (and also drops OSD-771) and reports the
5-mission meta. D1(d) cannot "flag" this away: it is perfectly confounded within the mission, and only
contrast choice addresses it.

**4. D1(c) sampling covariates are a mediator/collider trap, and the "sharpened K0" branch must
go.**
- The cortex-PT, medulla and blood panels move with flight: PT −1.45 in OSD-513, blood −2.17 in
  MHU-3, and Siew reports tubule-density loss. Adjusting for them can remove real flight effects.
- Structural tracks PT within cells (r ≈ 0.33), so "structural collapses, podocyte survives" is the
  pattern you would expect from adjusting for a mediator correlated with structural.
- The panels were hand-picked after the results, and pooled within-cell slopes impose a common
  slope across missions. Per-mission slopes overfit and attenuate to 0.079.

Fix:
- Derive the panels by rule from an external atlas (top-N region-specific genes in Ransick or MKA,
  excluding target-set genes) before rerunning.
- Estimate the slopes in **control animals only**.
- Report the result as "robust / not robust" only. Rank stability is not inference. K0 remains the
  sole arbiter of separability.

**5. Thin margin and forking paths are not quantified anywhere in the plan.**
- Podocyte HS is the only one of 49 sets with lower CI > 0 (spot-check 1): p_mkh 0.042, empirical
  p 0.00095, within-family max-T 0.019.
- The family was chosen *after* the 4 targeted axes failed. Bonferroni over the 53 tests actually
  examined (4 axes + 49 sets) gives about 0.00095 × 53 ≈ 0.05, and that is before counting the gene
  map, tier and adjustment variants.
- 32 of 49 sets are positive, consistent with a broad upward offset.
- Leave-top-5 genes crosses zero.

Fix:
- State this arithmetic in the paper.
- Add a **matched-random-gene-set null**, unless K0 or G1 already contains one. Use about 10k sets
  matched on size, mean expression, GC and length, score each through the same meta, and report the
  podocyte percentile. This is simpler and more decisive for "specificity" than the covariate panels.

**6. D3's headline has no defined denominator, and the novelty is narrower than the draft says.**
Siew 2024 already claims structural remodeling, and its PODXL/NID1 vote count runs opposite to our
direction. With four in-house axes, "most reported renal responses" is an overreach that an editor
will catch. Fix: build a **prior-claim recurrence ledger** (see Missing directions).

**7. D2's days-since-return grid recreates the recovery story the stop-list retires.**
- The grid has four points from four labs, crossing strain, sex, preservation and euthanasia.
- MHU-3 mice underwent open-field, light/dark and Y-maze tests before isoflurane exsanguination
  (spot-check 5), and its IEG g of 1.50 is consistent with that. Half the cohort is Nrf2-KO.
- R+2 is from ibSLS; I could not confirm it in the OSDR metadata.

Fix: OSD-457 enters a supplementary moderator table (WT primary, KO separate, genotype-blocked
descriptive), with no fitted curve and no main-text figure.

**8. The sequencing has paper drafting outrun G3.** Weeks 1–3 draft Results while G3 can still
remove the attribution. Fix: in weeks 1–2 draft only the Introduction, Methods and the 4-axis
non-recurrence. Write no podocyte or structural Results prose until the G3 and ledger gates resolve.

## Missing directions

- **M1. Prior-claim recurrence ledger (Tier 1, about 1 week).**
  - Before scoring, enumerate the field's specific renal transcriptomic claims: Siew 2024's top
    pathways and most-frequent genes (PODXL, NID1, transporters, focal adhesion, tight junction,
    actin), Suzuki 2022 MHU-3 (bone mineralisation, blood pressure, lipid), and DCT remodeling.
  - Map each claim to a gene set by a frozen rule.
  - Score each with the existing G1 machinery and calibrated variance.
  - Verdict per claim: recurs / does not recur / not evaluable.
  - This gives the title its denominator, turns the Siew conflict into a tested estimand difference,
    and reuses existing code.
- **M2. Matched-random-set specificity null.** See issue 5.
- **M3. OSD-771 mix-matched substitution.** See issue 3.
- **M4. Group-level design as a field-level methods result.** In every mission, each condition is
  one dissection day and one processing batch, so the effective group-level n is 1 vs 1. The
  control-vs-control calibration is the most transferable contribution and should be a named
  Methods/Results section, not a sensitivity footnote.
- **M5. Staged, bounded execution.** This is a budget point, given this project's history of
  exhausting context. Run each direction as one bounded agent task with:
  - a written stop rule,
  - a maximum file-read list,
  - outputs written to a run directory, and
  - a 1-page result note as the only handoff.

  No agent re-reads the reader scratchpads.

## Recommended final ordering

| Order | Work | Gate or decision it feeds |
|---|---|---|
| 1 (wk 1) | Blind gate-writing and hash-lock for G3. Write D4/D5 directional predictions before download. Send the D7 asks. | Makes everything after it non-theatre |
| 2 (wk 1–2) | **G3 core:** (a) calibrated control choice across all 49 sets; M3 OSD-771 mix-matched substitution; M2 matched-random-set null. (b) GC/length and (c) externally derived, control-slope sampling covariates are compact sensitivities, "robust / not robust" only. | Decides "spaceflight-associated" and "podocyte-leaning" |
| 3 (wk 2–3) | **M1 prior-claim recurrence ledger** | Defines the paper's denominator and the Siew positioning |
| 4 (wk 2–3) | **D4 cross-tissue** (heart/brain, RRRM-2): structural, handling panels, control-group deviation axis | Kidney-specific vs systemic/batch; genuinely unseen |
| 5 (wk 1–6) | **D3 drafting**, gated: Introduction, Methods, 4-axis first; podocyte and structural text only after steps 2–4 | Paper |
| 6 (wk 3) | **D2 OSD-457**, 2–3 days, supplement only | Audit correction and descriptive context |
| 7 (wk 4) | **D5**: OSD-163 coverage gate, then OSD-102 if it passes | Second assay layer |
| — | D7 replies, D9 watch; D6 only if a reviewer demands it; D8 dropped | — |

**Single strongest figure for D3: the recurrence map.**
- Rows: the prior claims (M1) plus the 4 axes and the structural and podocyte sets.
- Columns: missions, with a calibrated meta column and I².
- Each row is marked by its verdict.

**Must NOT be in it:**
- podocyte as the headline or a separate hero panel;
- G2 "rebound" or recovery;
- the days-since-return grid;
- rank-after-sampling-adjustment;
- MHU-3 pooled as a replicate;
- glycocalyx or sialylation mechanism;
- any reader post hoc number presented without its G3 rerun.

**Venue.** npj Microgravity is realistic. Commun Biol needs M1 and M4 to land cleanly.

## Spot-checks performed

| # | Claim | How verified | Result |
|---|---|---|---|
| 1 | Podocyte HS g 0.689 (0.042, 1.336), rank 1 of 49, structural 0.549 | `compartment_context_meta_results.tsv`, sorted by t | **Confirmed.** It is the only set of 49 with CI_low > 0. p_mkh 0.042, empirical p 0.00095, max_t_fwer 0.019; structural max_t_fwer 0.129. 32 of 49 estimates are positive. |
| 2 | Structural I² is 0.67%, not 67.3% | Same file: i_squared alongside Q, and the endothelial rows for units | **Confirmed.** Structural Q = 4.03 on df 4 gives I² = 0.67. Endothelial is 72–75 in the same column, so the column is in percent. |
| 3 | ERCC mix confounded with condition, reversed between arms (OSD-771) | `a_OSD-771_…Illumina.txt`, crosstab of arm × group × mix | **Confirmed.** ISS-T: FLT/VIV Mix 1, GC/BSL Mix 2. LAR is reversed. FLT−VIV is the only mix-matched flight contrast (new implication). |
| 4 | G2 podocyte interaction −0.28 (−0.62, 0.06); structural p 0.009, FWER 0.024 | `recovery_persistence/persistence_contrasts.tsv` | **Confirmed.** Podocyte −0.281 (−0.622, 0.061), FL p 0.078, INCONCLUSIVE. Structural −0.349 (−0.626, −0.073), FL p 0.0088, max-T 0.024. |
| 5 | OSD-457: behaviour tests before dissection, R+2 | One OSDR API call (`/osdr/data/osd/meta/457`) | **Partly confirmed.** Open-field, light/dark and Y-maze tests, then isoflurane exsanguination, are stated; 31 d on the ISS. **R+2 is not stated** in the fields searched; it rests on ibSLS per the external-data reader. |
