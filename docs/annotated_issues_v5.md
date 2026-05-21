# Annotated Issues — LAR Reversal Analysis & Manuscript v5

Companion to `docs/statistical_assessment_v5.md`. Severity: **[BLOCKER]** must fix before submission · **[MAJOR]** materially affects a claim · **[MINOR]** correctness/hygiene.

---

## Code issues

### [BLOCKER] Redundant feature inflates the state-space cosine
- **File:** `src/networks/lar_reversal.py:29-35` (`STATE_FEATURES`)
- `matrix_minus_dct` == `matrix_component − dct_transport_component` exactly (verified, residual 8.9e-16). The 5-feature "state space" is rank 4; matrix & DCT are double-weighted in every `cosine_similarity` / `vector_stats` / `projection_beta` call over `STATE_FEATURES`.
- **Effect:** pushes `cos(v_LAR, v_ISS)` artificially toward ±1, especially the young −0.993.
- **Fix:** remove `matrix_minus_dct` from `STATE_FEATURES`; keep it only as a scalar readout. Re-run the vector summary.

### [BLOCKER] Bootstrap CI and permutation p contradict each other for the young state-space result
- **File:** `src/networks/lar_reversal.py:341-354` (`reversal_summary_for_features`, bootstrap vs permutation blocks)
- For `state_space / YNG`: cosine −0.993, permutation p = 0.0025, **bootstrap 95% CI ≈ [−0.99, +0.78]** (crosses zero).
- The two resampling schemes disagree because n = 5/cell makes the cosine sample-driven. Manuscript reports the permutation/q, omits the CI.
- **Fix:** always report both; downgrade any claim where the bootstrap CI crosses zero.

### [MAJOR] `lar_opposes_iss` flags noise-level attenuation as reversal
- **File:** `src/networks/lar_reversal.py:489`
- `bool(np.sign(iss) != np.sign(lar) and abs(lar) > EPS and abs(iss) > EPS)` — set `True` whenever signs merely differ, so LAR values of −0.12 / −0.15 / −0.27 (not distinguishable from 0) are labeled "opposes ISS-T."
- **Fix:** require LAR to be significantly non-zero (CI excludes 0) before calling it opposition; otherwise label "attenuation."

### [MAJOR] Mechanism-switch decomposition is circular, ill-conditioned, and uncertainty-free
- **File:** `src/networks/lar_reversal.py:689-781` (`mechanism_axis_definitions`, `mechanism_switch_decomposition`)
- Axes are hand-weighted sums of the same features being decomposed; `np.linalg.lstsq` (`:766`) is ill-posed because axes are collinear (`simple_` vs `multi_axis_` coefficients disagree badly); output has no CI/p/q; `residual_norm_all_axes` (4.88) is ~82% of `effect_norm` (5.92) → axes explain only ~32% of variance.
- **Fix:** attach bootstrap CIs to projection coefficients and report variance explained, or move "three programs" to Supplementary as an explicitly exploratory hypothesis.

### [MAJOR] `node2vec_displacement_magnitude_summary` cosine is a meaningless +0.997
- **File:** `scripts/run_lar_reversal_analysis.py:345` + `src/networks/lar_reversal.py:114` (`node_summary_rows`)
- `D_ISS` / `D_LAR` are displacement *magnitudes* (non-negative); the cosine of two all-positive vectors is trivially ≈ +1.
- **Fix:** drop from `lar_reversal_vector_summary.tsv` or replace with a signed/directional metric.

### [MAJOR] BH correction applied across non-independent vector-summary rows
- **File:** `scripts/run_lar_reversal_analysis.py:156-163` (`add_q_values`)
- BH runs over all 17 rows; `state_space` pooled/YNG/OLD are nested subsets, and `mechanism`/`special`/`all_gene` are overlapping views of one matrix. q-values (incl. the cited young q = 0.0375) are optimistic.
- **Fix:** define the test family explicitly; correct within it (or report raw p with family stated).

### [MINOR] Cross-run input provenance
- **File:** `scripts/run_lar_reversal_analysis.py:47-50`
- State scores from `run_20260519_000547_2500g`, WGCNA from `run_20260505_remediated_2500g`, cross-OSDR from `run_20260519_cosine_perm_null`. Auditable via manifest SHAs, but one figure spans three runs.
- **Fix:** pin all inputs to a single run or document the divergence.

### [MINOR] Single shared RNG threads all resampling blocks
- **File:** `scripts/run_lar_reversal_analysis.py:263` (one `rng`, reused across every call at `:325-371`)
- Reproducible globally, but no block is independently reproducible; reordering feature sets changes all downstream draws.
- **Fix:** `rng.spawn()` a child generator per analysis block.

### [MINOR] Interaction ≈ ISS-T effect because pooled LAR ≈ 0
- **File:** `src/networks/lar_reversal.py:455-456` (`interaction_vec`)
- Since pooled LAR effects ≈ 0, `interaction = ISS − LAR ≈ ISS`; the permutation interaction test largely re-detects "ISS-T has a flight effect" (standard DGEA). Verified: OLS `lm(score ~ flight*arm + age)` reproduces the permutation interaction estimates/p-values.
- **Fix:** state in Methods that the component interaction is equivalent to a linear-model interaction term; do not present it as a network method.

### [MINOR] Unit tests cover math, not inference
- **File:** `tests/test_lar_reversal.py`
- Tests verify cosine/classification on toy inputs; none check false-positive rate of the permutation test or bootstrap-CI coverage.
- **Fix:** add a null-simulation test (random labels → permutation p uniform; effect ≈ 0).

---

## Manuscript issues (`latex_paper/manuscript_v5.tex`)

### [BLOCKER] Young-LAR headline rests on a non-significant number
- **Lines 102, 444, 477, 486, 556** — young LAR `matrix_minus_dct` = −0.782 is presented as "moves in the opposite direction." Actual Welch p = 0.361, q = 0.542. Pooled LAR effect = −0.112 (p = 0.865).
- **Fix:** reword to attenuation (pooled LAR ≈ 0); present the young-specific value only as an underpowered, non-significant observation.

### [BLOCKER] Young state-space cosine cited without its bootstrap CI
- **Line 477** — "cosine = −0.993, directional permutation p = 0.0025, BH q = 0.0375." Bootstrap CI [−0.99, +0.78] (same output row) is omitted.
- **Fix:** report the CI; discuss the bootstrap/permutation disagreement.

### [MAJOR] Title and claim hierarchy over-weight the LAR result
- **Lines 53-54** (title), **552-560** (claim hierarchy), abstract **99-103**.
- v5 reframes around "diverge along matrix/S1P, clock-DCT, and tubular-transport remodeling programs" — Claims 4 & 5 are the weakest material.
- **Fix:** make the matrix-high/DCT-low recurrent state (Claim 3) the headline; adopt a recurrence-framed or the chat's "safer" arm-divergence title.

### [MAJOR] "Three separable programs" presented as a result
- **Lines 77, 103, 488-492, 557, 567-577** — built on the circular, 32%-fit mechanism-switch decomposition.
- **Fix:** label as hypothesis; move Panel D logic to Supplementary or attach uncertainty.

### [MINOR] Age-stratified interaction asserted but not tabulated
- **Lines 486, 556** reference the age pattern; no age-stratified Arm×Flight interaction table with CIs/q exists (only descriptive effects + one pooled interaction).
- **Fix:** add the explicit age-stratified interaction table; report it as underpowered.

### [MINOR] No standalone DGEA confirmation
- A `limma`/`DESeq2` `~ age + arm + flight + arm:flight` model and the ISS-T-vs-LAR gene logFC scatter would pre-empt the "is this just DE?" question. (Scatter data already in `lar_reversal_gene_scatter.tsv`; quadrants are ~50/50, itself evidence against a genome-wide reversal.)
- **Fix:** add a DGEA confirmation paragraph + figure.
