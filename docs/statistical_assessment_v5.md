# Statistical Robustness & Biological-Claim Assessment — RRRM-2 Kidney Transcriptome (Manuscript v5)

**Reviewer pass:** May 2026 · **Production run reviewed:** `run_20260519_000547_2500g` · **Manuscript:** `latex_paper/manuscript_v5.tex` (and v4 for comparison)
**Code reviewed:** `src/networks/lar_reversal.py`, `scripts/run_lar_reversal_analysis.py`, `scripts/plot_lar_reversal_dashboard.py`, `tests/test_lar_reversal.py`, plus the state-space / mechanism-axis inputs that feed them.
**Verification:** all 15 unit tests re-run (pass); the LAR reversal vector statistics, component effects, and interaction tests were re-executed and re-derived independently. Numbers below were reproduced, not taken on trust.

---

## 1. Executive verdict

The project is **not** "we did DE and WGCNA and found this." It is a genuinely thoughtful cross-cohort contrast-vector study, and the manuscript is, in most places, honestly hedged. But the headline that v5 is built around — the **age-specific LAR reversal / three-program story** — is the weakest-supported part of the paper, and it is currently *over-promoted relative to what the statistics support*.

Three findings drive that verdict, all verified against the code and outputs:

1. **The single number the v5 Results section and the chat discourse treat as the sharpest result — young LAR matrix-minus-DCT effect = −0.782 — is not statistically distinguishable from zero** (Welch p = 0.36, BH q = 0.54). The "young LAR goes the other way" claim rests on a non-significant point estimate.
2. **The young state-space cosine of −0.993 has a bootstrap 95% CI of roughly [−0.99, +0.78]** — it spans almost the entire possible range and crosses zero. The manuscript reports the permutation p (0.0025) and q (0.0375) for this result but **omits the bootstrap CI**, which tells the opposite story. When the two resampling schemes disagree this hard, "robust" is not a defensible word.
3. **The statistically defensible part of the LAR analysis — the component arm×flight interactions — is exactly reproducible by a one-line ordinary-least-squares interaction model.** I verified this: `lm(score ~ flight*arm + age)` gives matrix interaction β = +1.738, p = 0.0006 and matrix−DCT β = +2.366, p = 0.011, matching the script's permutation output (q = 0.005 and 0.032) to within rounding. The permutation/bootstrap wrapper does not add discovery power; it re-prices a standard differential gene-set-score test.

None of this kills the project. The cross-cohort ISS-T/OSD-513 alignment (the v4 core) is real and modestly supported. But the claim hierarchy needs to be re-leveled: the LAR section should be demoted from "biological center of the paper" to "an arm-divergence observation with one suggestive, underpowered age interaction."

---

## 2. What the dataset can and cannot support

RRRM-2 (OSD-771) is **n = 5 per Age × Arm × Environment cell** (verified from `state_space_scores.tsv`: 16 cells × 5 = 80 RRRM-2 samples). The LAR reversal analysis uses only FLT and GC, so:

- Pooled arm effect: 10 FLT vs 10 GC per arm (age-balanced average of two 5-vs-5 contrasts).
- **Age-stratified arm effect: 5 FLT vs 5 GC.** Every "young LAR" / "old LAR" number in the manuscript and dashboard rests on a 5-vs-5 difference of means.

A 5-vs-5 mean difference can support a *large, consistent* main effect. It cannot support a **three-way Age × Arm × Flight interaction** as a confirmatory claim — that estimand divides the already-small signal across 8 cells. The manuscript's own Limitations section says exactly this ("age-stratified LAR results [are] descriptive even when the effect size is large"). The problem is that the **Results section and the abstract do not honor that caveat** — they present −0.782 and the −0.993 cosine as findings, and the title and claim hierarchy are reorganized around them.

So the honest framing of the data envelope:

| Estimand | Supportable? |
|---|---|
| ISS-T flight effect (matrix up, DCT down), pooled | Yes — large, consistent, externally echoed |
| LAR flight effect, pooled | Yes — and it is ≈ 0 (attenuation), see §4.5 |
| ISS-T vs OSD-513 directional alignment | Yes, modestly (bootstrap CI excludes 0; permutation p = 0.127) |
| Arm × Flight interaction (pooled), matrix axis | Yes — but it is ordinary DE on a gene-set score (§4.6) |
| **Age × Arm × Flight interaction (the young-LAR-reversal claim)** | **No — underpowered; the key number is n.s.** |

---

## 3. Statistical robustness — detailed findings

### 3.1 `matrix_minus_dct` is a redundant feature and it inflates the state-space cosine

`STATE_FEATURES` (`src/networks/lar_reversal.py:29`) contains five features, but I verified that one is a deterministic linear combination of two others:

```
matrix_minus_dct  ==  matrix_component − dct_transport_component
(max abs residual across all 219 samples = 8.9e-16)
```

So the "5-dimensional state space" is really **rank 4**, and the matrix and DCT directions are **double-counted** in every cosine and norm computed over `STATE_FEATURES`. Because the ISS-T flight effect loads almost entirely on the matrix↑/DCT↓ contrast, including `matrix_minus_dct` mechanically pushes `cos(v_LAR, v_ISS)` toward ±1. The young-cohort cosine of −0.993 is *partly an artifact of this redundancy*: when both arms' effect vectors are dominated by one shared axis that happens to point opposite ways, near-perfect anti-alignment is almost guaranteed in a rank-4 space.

**Fix:** drop `matrix_minus_dct` from any feature *set* used for cosine/projection (keep it only as a scalar readout), or compute the cosine on the genuinely independent components only.

### 3.2 The cosine is computed over far too few effective dimensions

`cos(v_LAR, v_ISS)` is reported for `state_space` (5 features, eff. rank 4), `mechanism_scores` (11 features), and `special_clock_s1p_rbm3` (12). These are not 5–12 *independent* measurements. `matrix_component`, `ecm_organization`, `fibrosis_tgfb_emt`, `integrin_cell_adhesion`, `mmp_adam_proteolysis` are all matrix/ECM gene-set scores and are strongly mutually correlated. The effective dimensionality of the mechanism vector is closer to **2–3** (a matrix/ECM/innate axis, a DCT/transport axis, maybe a clock axis). A cosine over ~2–3 effective dimensions is a high-variance statistic, and its sign is essentially set by the matrix axis. That is why the cosine "result" so closely tracks the gene-set-mean comparison the user already did by hand in chat.

### 3.3 The central problem: bootstrap and permutation disagree, and the manuscript reports only the favorable one

Re-derived from the actual code (`reversal_summary_for_features`, seed 20260520, 2000×2000):

| Feature set | Age | cos(LAR,ISS) | Bootstrap 95% CI | Permutation p | Manuscript reports |
|---|---|---|---|---|---|
| state_space | pooled | −0.493 | **[−0.962, +0.623]** | 0.258 | "weak … attenuation" ✓ |
| state_space | **YNG** | **−0.993** | **[−0.99, +0.78]** | **0.0025** | cosine, perm p, q — **CI omitted** |
| state_space | OLD | −0.133 | [−0.886, +0.672] | 0.439 | "weak, inconclusive" ✓ |
| mechanism_scores | pooled | −0.725 | [−0.879, −0.082] | 0.014 | cosine + CI + p ✓ |
| mechanism_scores | YNG | −0.766 | [−0.876, −0.114] | 0.012 | — |
| mechanism_scores | OLD | −0.495 | [−0.808, +0.282] | 0.106 | — |

The young state-space row is the engine of the v5 reframing, and it is **internally contradictory**: the permutation null says the FLT/GC labels are not exchangeable (p = 0.0025), but the bootstrap says the *direction* of the effect — the −0.993 cosine itself — is not a stable estimate (CI runs from −0.99 all the way to +0.78). Both came out of the same function call; the manuscript cites the permutation/q from that row and silently drops the CI.

Why they disagree: with n = 5 per cell, a bootstrap resample routinely triples one mouse and drops two others. The cosine is a ratio, so one or two influential mice flip it. The permutation null, by contrast, is narrow because shuffling labels within n = 5 cells produces small, structureless effect vectors — so *any* structured observed vector lands in the tail. The correct reading is **"there is a flight effect, but the −0.993 anti-alignment is driven by a couple of mice and is not a stable directional estimate."**

A defensible paper cannot report q = 0.0375 for a quantity whose own bootstrap CI crosses zero. Either both are shown side by side with the contradiction discussed, or the claim is downgraded.

### 3.4 The headline young-LAR effect size is not significant

From `lar_reversal_component_effects.tsv` (re-confirmed):

```
matrix_minus_dct  ISS-T  YNG   effect = +3.204   Welch p = 0.0046   q = 0.025   (significant)
matrix_minus_dct  LAR    YNG   effect = −0.782   Welch p = 0.361    q = 0.542   (NOT significant)
matrix_minus_dct  LAR    pooled effect = −0.112   Welch p = 0.865    q = 0.895   (≈ zero)
```

The young ISS-T shift is real. The young LAR −0.782 — the number the chat discourse called "GOING THE OTHER WAY" and the number v5's Results calls "moves in the opposite direction" — has p = 0.36. It is noise around zero. The *pooled* LAR effect is −0.11 (p = 0.87): textbook **Model A attenuation**, not reversal.

So the strongest honest statement is: *"LAR shows no detectable matrix−DCT flight response (pooled effect ≈ 0); the Arm×Flight interaction is significant because ISS-T responds and LAR does not."* That is attenuation. "Reversal," and especially "young-specific reversal," is not supported by the effect-level statistics.

### 3.5 The component interaction tests are ordinary least squares with extra steps (verified)

I ran a plain OLS model `score ~ flight + arm + age + flight:arm` on the RRRM-2 FLT/GC samples:

```
matrix_minus_dct        Arm×Flight interaction  β = +2.366   t = +2.68   p = 0.011
matrix_component        Arm×Flight interaction  β = +1.738   t = +3.76   p = 0.0006
dct_transport_component Arm×Flight interaction  β = −0.628   t = −1.02   p = 0.314
```

These match `lar_reversal_pathway_interactions.tsv` (matrix β = 1.738, q = 0.005; matrix−DCT β = 2.366, perm p = 0.019, q = 0.032; DCT q = 0.344) essentially exactly. **The permutation + bootstrap apparatus reproduces a standard linear-model interaction term** — the same term `limma`/`DESeq2` fit per gene. This is important for the "could DGEA do this?" question (§5): the *defensible* layer of the LAR analysis is differential gene-set testing. The custom machinery is a re-implementation, not a new capability.

A subtler point: because pooled LAR effects are ≈ 0 for nearly every feature, `interaction = ISS_effect − LAR_effect ≈ ISS_effect`. The interaction permutation test is therefore largely **detecting "ISS-T has a flight effect"** — pure DGEA — and `lar_opposes_iss` (set whenever the two signs merely differ, `lar_reversal.py:489`) flags small noise-level LAR values (−0.12, −0.15, −0.27) as "opposes," which over-reads attenuation as reversal.

### 3.6 The mechanism-switch decomposition is circular and carries no uncertainty

`mechanism_axis_definitions` (`lar_reversal.py:689`) hand-defines axes (`matrix_high_dct_low`, `circadian_dct_state`, `s1p_matrix`, `transport_rebound`, …) as weighted sums of **the same features** being decomposed. Projecting an effect vector onto hand-weighted combinations of its own components is descriptive bookkeeping, not a test. Three concrete problems:

- **No uncertainty.** `lar_mechanism_switch_scores.tsv` has point estimates only — no CI, no p, no q. The v5 "Panel D" numbers (circadian_dct ISS −1.94 / LAR +0.31; s1p_matrix +1.40 / −0.91; transport_rebound −0.35 / −0.91) are unquantified point estimates presented as a "three-program" finding.
- **Collinear, ill-conditioned basis.** The axes overlap heavily, so `np.linalg.lstsq` is ill-posed; `simple_projection_coefficient` and `multi_axis_regression_coefficient` disagree badly (transport_rebound ISS-T: −1.53 vs −0.35). Which one is "the" coefficient is arbitrary.
- **Poor fit.** For ISS-T, `residual_norm_all_axes = 4.88` vs `effect_norm = 5.92` → the six axes explain only ~32% of the effect-vector variance. The "decomposition into three clean programs" framing is not what the numbers show; two-thirds of the signal is residual.

### 3.7 `node2vec_displacement_magnitude_summary` cosine = +0.997 is a non-result

`node_summary_rows` computes `cos(D_ISS, D_LAR)` where `D` are **displacement *magnitudes*** — non-negative by construction. The cosine of two all-positive vectors is trivially near +1; the observed +0.9972 carries no information. It is labeled "exploratory" and the manuscript explains it away ("as expected for a non-directional movement metric"), which is fine — but it should not sit in `lar_reversal_vector_summary.tsv` as if it were a comparable measurement. Either drop it or compute a directional quantity.

### 3.8 Multiple-testing correction is applied across non-independent rows

`add_q_values` (`run_lar_reversal_analysis.py:156`) BH-corrects across all 17 rows of the vector summary. Those rows are not independent tests: `state_space` pooled/YNG/OLD are nested subsets of the same samples; `mechanism_scores`, `special_clock_s1p_rbm3`, and `all_gene_expression` are overlapping views of one expression matrix. BH assumes independence or positive dependence across *distinct* hypotheses; nested age scopes of one feature set violate the spirit of the family. The q-values are therefore optimistic, and they are also the q-values the manuscript quotes (e.g. young state-space q = 0.0375). Define the test family explicitly and correct within it, or report raw p with the family stated.

### 3.9 Cross-run input provenance

`run_lar_reversal_analysis.py` defaults mix runs: state scores from `run_20260519_000547_2500g`, WGCNA eigengenes from `run_20260505_remediated_2500g` (`:49`), cross-OSDR summary from `run_20260519_cosine_perm_null` (`:50`). The manifest records SHA-256 of each input, so it is auditable — but a single manuscript figure built from three different pipeline runs is a reproducibility smell and a reviewer will ask. Pin all inputs to one run, or document why each differs.

### 3.10 A shared RNG threads all resampling

One `np.random.default_rng(seed)` is created in `main()` and passed through every `reversal_summary_for_features` and `interaction_table` call in sequence (`run_lar_reversal_analysis.py:263` onward). Results are reproducible as a whole, but a single feature set's bootstrap is not independently reproducible, and reordering feature sets changes every downstream draw. Minor, but for a paper, give each analysis block its own seeded child generator (`rng.spawn`).

### 3.11 What is statistically solid

To be fair to the work — these parts hold up:

- **ISS-T pooled flight effect.** matrix +1.25 (p = 0.002), DCT −1.00 (p = 0.035), matrix−DCT +2.25 (p = 0.003). Large, consistent, age-robust at the pooled level.
- **ISS-T / OSD-513 alignment.** cosine +0.641, bootstrap CI [0.351, 0.814] excludes zero. The label-permutation p = 0.127 does not reach 0.05, and the manuscript says so plainly. "Robust under bootstrap, not significant under the external-label null" is a fair, correctly stated claim.
- **Pooled mechanism-score anti-alignment.** cosine −0.725, CI [−0.879, −0.082] excludes zero, perm p = 0.014. This is the *one* LAR-divergence number that survives its own bootstrap. It supports "LAR diverges from ISS-T at mechanism-score scale" — i.e. arm divergence — but **not** "reversal" and **not** "young-specific."
- **The unit tests pass** (15/15) and the manifest/SHA discipline is genuinely good practice.
- The manuscript's Methods, Limitations, and the OSD-253 "association not requirement" framing are careful and honest.

---

## 4. Are the biological claims sound?

Going down the v5 claim hierarchy (Results §"Claim hierarchy", lines 552–560):

**Claim 1 — "RRRM-2 cannot resolve the original age-rewiring contrast."** Sound and honestly stated. Good scientific hygiene to lead with a negative result.

**Claim 2 — "ISS-T aligns with OSD-513; LAR anti-aligns."** The ISS-T half is sound (bootstrap-supported, permutation-suggestive). The LAR "anti-aligns" half is weaker than stated: cosine −0.511 with CI [−0.726, −0.231] excludes zero, but permutation p = 0.178. "Anti-aligns" is acceptable as a directional description if paired with "not significant under the permutation null" — which the manuscript does do in the abstract. Acceptable, borderline.

**Claim 3 — "The aligned ISS-T direction is a matrix-high / DCT-low tubulointerstitial state."** Sound. This is the strongest, most defensible claim in the paper and should be the headline.

**Claim 4 — "Young LAR moves opposite to young ISS-T; old LAR does not."** **Not sound as stated.** The young LAR effect (−0.782) is not significant (p = 0.36); the young cosine (−0.993) has a bootstrap CI crossing zero. What the data support is: *ISS-T responds and LAR does not, and the failure-to-respond may be larger in young mice* — an attenuation-with-possible-age-modulation observation, not a reversal.

**Claim 5 — "Three separable programs (matrix/S1P/innate reversal, clock-DCT rebound, persistent transport disruption)."** **Not sound as a confirmatory claim.** It rests on the circular, uncertainty-free mechanism-switch decomposition (§3.6) that explains only ~32% of the effect variance. As a *hypothesis* it is interesting and worth stating — but it is currently written as a result.

**Claim 6 — "OSD-253: TLR4 association, not requirement."** Sound and carefully hedged. The strain-delta CIs cross zero and the manuscript says so.

**Claim 7 — "WGCNA / LIONESS / node2vec are annotation and prioritization."** Sound — correctly demoted, FDR-weak signals honestly reported.

**Net:** Claims 1, 3, 6, 7 are sound. Claim 2 is borderline-but-defensible. **Claims 4 and 5 — the entire reason v5 exists — are over-stated.** They should be reworded as hypotheses, or supported with the analyses in §6.

---

## 5. Are the claims novel — and could plain DGEA/DGEA recover them?

This is the question the user keeps circling, so directly:

**The matrix↑ / DCT↓ ISS-T signal: yes, DGEA recovers it.** The mechanism/state scores *are* gene-set mean expression scores; a FLT-vs-GC contrast on them is gene-set differential expression. The chat discourse already showed this by hand ("ISS-T matrix +0.201"). The component-effects table is Welch tests on gene-set scores. Verified equivalence in §3.5.

**The Arm × Flight interaction: yes, DGEA recovers it.** A one-line `lm(score ~ flight*arm + age)` reproduces the permutation interaction estimates and p-values exactly (§3.5). The custom permutation framework adds no discovery power over standard limma/DESeq2 with an interaction term.

**What is NOT trivially DGEA:**
- The **cross-cohort directional alignment** (ISS-T vs OSD-513 cosine, with bootstrap CI and a label-permutation null). This genuinely needs the vector framing — a per-gene DE table does not give you "the RRRM-2 ISS-T pathway vector points the same way as the OSD-513 pathway vector." This is the real methodological contribution and it is under-sold relative to the LAR section.
- The **explicit three-model competition** (attenuation vs reversal vs switch) as an *analysis design* is a good idea and is not standard DGEA. But as currently executed it does not actually adjudicate the models — the data land on Model A (attenuation) at the pooled level, and the script's own classifier returns three different model calls for three feature sets (`state_space`→A, `mechanism`→B, `special`→A; see `lar_reversal_model_calls.tsv`). An honest reading of the framework's own output is "attenuation, with one underpowered age interaction," not "reversal/switch."

**Novelty bottom line.** The novel, defensible contribution is the **cross-OSDR arm-specific remodeling-axis recurrence** — ISS-T captures a matrix/DCT remodeling direction that recurs in an independent flight cohort, while the live-return arm does not. That is more than "DE + WGCNA," and it is biologically interesting (it reframes spaceflight kidney response as arm/recovery-context-dependent). The "network rewiring discovery engine" and the "LAR three-program reversal" framings are *not* where the novelty is, and the v5 reframing unfortunately moved the spotlight onto the weakest material. The user's own instinct in the chat — "the novelty is the arm-specific remodeling geometry, especially the LAR decoupling phenotype" — is right; the execution just over-claimed the LAR half.

---

## 6. How to fix it — concrete, in priority order

**A. Re-level the claim hierarchy (no new computation).** Make Claim 3 (matrix-high/DCT-low recurrent state) the headline. Demote LAR to: *"the live-return arm shows no detectable matrix−DCT flight response (pooled effect ≈ 0); whether this is attenuation or active reversal, and whether it is age-modulated, is not resolved at n = 5."* Retitle away from "diverge along … programs" toward the v4-style recurrence framing or the "safer title" the chat itself proposed.

**B. Report the bootstrap CI everywhere the permutation p is reported.** Specifically, the young state-space row must show CI [−0.99, +0.78] next to p = 0.0025, with one sentence explaining the disagreement. Reviewers will find this; better to own it.

**C. Drop `matrix_minus_dct` from feature *sets*** used for cosine/projection (§3.1). Keep it only as a scalar. Re-run; expect the state-space cosines to move toward the (more honest, less extreme) mechanism-score values.

**D. Replace the mechanism-switch decomposition with a real test or label it exploratory.** Either attach bootstrap CIs to the projection coefficients and report the fraction of variance explained, or move the whole "three programs" panel to Supplementary and call it a hypothesis. Do not present a rank-deficient, 32%-fit projection as a result.

**E. Formalize the age interaction properly.** The chat already identified the right move: an explicit **age-stratified Arm × Flight interaction table** for `matrix_minus_dct`, `s1p_matrix`, `circadian_dct_state`, `transport_rebound`, with CIs and q-values — and then *report it as underpowered*. If the young-specific Age×Arm×Flight term is not significant (it almost certainly will not be at n = 5), say so. A clearly-stated underpowered interaction is publishable; a non-significant number presented as a finding is not.

**F. Run a genuine DGEA confirmation arm.** `limma`-voom or `DESeq2` with `~ age + arm + flight + arm:flight` on counts. Make the ISS-T vs LAR gene-level logFC scatter (the data already exist in `lar_reversal_gene_scatter.tsv` — note the quadrants are ~50/50: 4440 ISS-up/LAR-down vs 4345 ISS-down/LAR-up, which itself argues *against* a genome-wide reversal and *for* "LAR ≈ noise"). Showing DGEA agrees is a strength, not a weakness; it pre-empts the reviewer's first question.

**G. Fix `node2vec_displacement_magnitude_summary`** — drop it from the vector summary or replace the magnitude cosine with a signed/directional quantity (§3.7).

**H. Pin all inputs to one run** (§3.9) and give each resampling block its own child RNG (§3.10).

**I. Power/sensitivity statement.** Add an explicit note: with n = 5/cell, the minimum detectable Age×Arm×Flight interaction effect size is large; the young-LAR observation is hypothesis-generating and needs an independent cohort or a within-RRRM-2 replication design to confirm.

---

## 7. Verification log

- Environment: repo `venv` symlink is broken in this sandbox; used system Python 3 with `scipy 1.15.3`, `pandas`, `pytest` installed.
- `pytest tests/test_lar_reversal.py tests/test_tubulointerstitial_state.py tests/test_mechanism_axis.py` → **15 passed**.
- Re-ran `reversal_summary_for_features` for `state_space` and `mechanism_scores` (pooled/YNG/OLD, 2000 bootstrap × 2000 permutation, seed 20260520) → reproduced `lar_reversal_vector_summary.tsv` cosines, CIs, and permutation p-values to 3 decimals.
- Verified `matrix_minus_dct == matrix_component − dct_transport_component` (max abs residual 8.9e-16 over 219 rows).
- Re-derived component effects: young LAR `matrix_minus_dct` = −0.782, Welch p = 0.361, q = 0.542 (matches `lar_reversal_component_effects.tsv`).
- Independently fit OLS `score ~ flight*arm + age` on RRRM-2 FLT/GC: matrix interaction β = +1.738 (p = 0.0006), matrix−DCT β = +2.366 (p = 0.011), DCT β = −0.628 (p = 0.314) — matches the script's permutation interaction table.
- Confirmed RRRM-2 design is n = 5 per Age×Arm×condition cell (80 samples; FLT/GC contrasts are 10v10 pooled, 5v5 age-stratified).
- Sample-size and quadrant counts read directly from `state_space_scores.tsv` and `lar_reversal_gene_scatter.tsv`.

See `docs/annotated_issues_v5.md` for the file-and-line issue list.
