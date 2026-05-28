# Manuscript v11 — Repository Audit

Scope: statistical robustness, biological claims, interpretability, results, plots,
figures. Cross-checked `latex_paper/manuscript_v11.tex` and `results_v11.tex`
against the v11 run root
`data/results/run_20260526_v11_dct1_phospho_mediation/`, the OSD-462 anchor run,
the regulator-activity run, the spatial-reference outputs, the public
perturbation-reference outputs, and the v11 source modules under `src/v11/`.

**Round 3 status (2026-05-28).** After the user's round-2 edits, a re-audit
verified every numerical claim against the May-28 TSV/JSON and identified 14
remaining items, of which 9 substantive ones have now been incorporated. The
manuscript and compendium both compile cleanly (27 and 17 pages). Section 0a
below lists the new fixes; the previous round-1 incorporation summary is
preserved as Section 0b.

## 0a. Round 3 fixes (this pass)

| Round-2 finding | Action taken |
|---|---|
| #1 DCT2-bottom-decile comparator under-reported | Added a side-by-side comparator table (manuscript Table 4 `tab:dct2compare`; compendium Table S\ref{tab:dct2comparesupp}) and one-sentence summaries in the Abstract, Results §3.4, and Conclusion. The supported claim is now "distal-nephron subtype-prior enrichment, strongest in the DCT1 top decile" rather than "DCT1-high parent-gene enrichment" without qualification. |
| #2 Compendium S5 vs S8.4/S8.5 duplication | Section S5 collapsed to a single forward-reference paragraph pointing to the perturbation-reference section. |
| #3 Compendium Table S1 missing perturbation modules | Provenance table now has a separate "Public perturbation-reference modules" row listing all four `src/v11/` perturbation modules. |
| #4 Naming inconsistency | Standardized to "public perturbation-reference analyses" across manuscript Methods §2.8, Results §3.7, Discussion §4.4, Author Contributions, Code Availability, compendium executive verdicts, perturbation section header, and figure caption. |
| #6 Author Contributions incomplete | Added "public perturbation-reference (low-K DCT, parent-protein-normalized phosphosite, IRI Visium DCT-transport, dDAVP/mpkDCT, KLHL3/CUL3)" to the I.S. contributions list. |
| #8 Matched-null "weakly" wording | Now stated as "supportive at FDR 0.10" with both the primary_p05 q=0.063 and the strict_q10 q=6.46e−4 explicit. |
| #9 Composition-aware ladder q values stale | Re-extracted from `h2_composition_adjusted_suppression_enrichment_single.tsv` and re-derived one-sided Fisher p-values from the contingency tables; M0 p=1.06e−7 matches one-sided exactly, M4 p=6.59e−4 matches one-sided exactly. Compendium Table S\ref{tab:h2ladder} is **not stale**: the ladder already uses directional Fisher. Item closed. |
| #10 One-site matched-null counter-result not in main text | Now reported in Results §3.4: "failed in the one-site-per-parent-gene reduction (observed mean DCT1 score among one-site suppressed parent genes was below the matched null; directional empirical p=0.996)". |
| #11 n=20 not in mediation description | Added "(n=20 animals with matched RNA and phosphoproteomics)" inline in Results §3.6. |
| #14 Duplicated opening sentence Results §3.7 vs Discussion §4.4 | Rewritten — Results §3.7 keeps the original, Discussion §4.4 now opens with "Among the five public perturbation references, the parent-protein-normalized re-test of the DCT1-high enrichment carried the strongest result". |

## 0a.1 Remaining round-2 items deferred (not fixed in this pass)

| Item | Reason | Recommendation |
|---|---|---|
| #5 Compendium parent-protein-normalized table missing composite-excluded row | The h2_occupancy enrichment TSV doesn't include a composite-excluded sensitivity for the parent-protein-normalized re-test — only `occupancy_p05`, `occupancy_q10`, `occupancy_exclude_anchor_genes`, `occupancy_exclude_ncc_sites`, `occupancy_single_site_p05`. Adding it requires re-running `h2_occupancy_normalized_phospho.py` with the new composite-exclusion filter. | Either re-run the module to produce the composite-excluded occupancy row, or add an explicit note to the table caption stating "composite-excluded sensitivity not run on the parent-protein-normalized version". |
| #7 Figure 5 panel D y-axis dominated by DCT-marker-high spike | Visual change; needs `publication_figures.py` modification to split panel D into two subpanels or annotate the day-14/week-6 DCT-adjacent points with their values. | Not a numerical bug; flag for the figure module owner. |
| #12 No `tests/v11/` directory | Reproducibility/regression-test gap. The May-28 re-run used seed=20260526 and the inputs are SHA-256-locked, but a future code change could silently break the contingency-table construction or the BH family. | Add fixture tests for the directional Fisher + BH (verifies #1), parent-gene Fisher and logistic, and a figure-snapshot test that catches the kind of blank-row bug that appeared and was later fixed in `v11_dct1_parent_gene_enrichment.pdf`. |
| #13 Methods §2.5 tie-breaker text not code-reviewed | The methods description of the one-site-per-parent-gene tie-breaker reads correctly, but the implementation in `src/v11/core_analysis.py` (or its helper) was not directly inspected for the exact tie-breaker order. | One-pass code review of the helper that emits `one_site_per_parent_gene` to confirm the implementation matches "lowest p, then more negative effect, then site identifier". |

---



Overall: the manuscript text is unusually well aligned with the underlying TSV/JSON
artifacts. The framing is appropriately hedged ("parent-gene enrichment, not
cell-of-origin"; "spatial contextualization, not validation"). The issues below
are tractable revisions, not headline-breakers.

## 0b. Round 2 incorporation: public perturbation-reference layer

| Branch | Run-root artifact | Headline number | Manuscript placement |
|---|---|---|---|
| Low-K DCT (GSE228367) | `perturbation/lowk_alignment_verdict.json` | Transport-target cosine OSD-462 −0.286 (CI −0.752 to 0.458; ρ=−0.420, q=0.218); RRRM-2 −0.311; OSD-513 −0.146 | Methods §2.8; Results §3.7 ("directional but not decisive"); Table 4 secondary; promotion rule not met |
| Occupancy-normalized DCT1 enrichment | `h2_occupancy/h2_occupancy_verdict.json` | Top-decile OR 1.52 (CI 1.36–1.70, q=3.73e−12); survives anchor, NCC, single-site, strict-q sensitivities | Abstract ("OR 1.52, q=3.73e−12"); Methods §2.8; Results §3.7 + Table S; Conclusion |
| IRI Visium DCT-transport spatial context | `spatial_reference/visium_dct_transport_verdict.json` | DCT-adjacent: day 14 −0.044 (p=0.005); week 6 −0.034 (p=0.018). DCT-marker-high: day 14 +0.199 (p=2.25e−15); week 6 +0.132 (p=2.70e−9) | Abstract; Methods §2.8; Results §3.7; Future Experiments (testable spatial prediction); Conclusion |
| PXD001729 dDAVP | `h2_pxd/` | 60 shared single sites; 0 shared transport-target sites; cosine −0.026 (CI −0.310 to 0.287) | Results §3.7 (coverage-bounded); Table 4 updated |
| KLHL3/CUL3 turnover | `h2_klhl3/` | 0 KLHL3 phosphosites; 1 CUL3 site; flat WNK/SPAK/NCC protein abundance | Results §3.7; Limitations; Future Experiments |
| Figure | `figures/v11/v11_perturbation_triangulation.pdf` | 4-panel: low-K cosines, target gene heatmap, occupancy forest, IRI DCT-transport by spatial context | Figure 5 in manuscript; Figure S6 in compendium |

The triangulation layer changes the manuscript's stance in three concrete ways:

1. **Abstract and Conclusion no longer claim DCT1 enrichment holds only after
   compartment-score adjustment.** They now also cite the parent-protein-normalized
   (occupancy) re-test at OR 1.52, q=3.73e−12 — a per-site normalization that is
   more direct than the M1 protein-covariate adjustment and gives a slightly
   stronger result.
2. **The upstream-mechanism question is now bounded rather than left open.** Intro
   question 5, Methods §2.8, Results §3.7, and Discussion §4.4 jointly explain
   that public references provide partial directional support (low-K) and
   coverage bounds (PXD, KLHL3) but do not identify the upstream WNK suppressor.
3. **The future spatial experiment now carries a specific prediction.** The IRI
   Visium DCT-transport-by-spot panel predicts DCT-adjacent transport drop with
   intact DCT-marker-high transport, which a single spaceflight spatial
   experiment can confirm or falsify.

---

## 1. Verified numbers (manuscript ↔ repo)

All headline quantities in Abstract, Results, and Tables 2–4 reconcile with the
generating TSV/JSON. Spot-checks against the canonical artifacts:

| Manuscript claim | Repo value | Source |
|---|---|---|
| OSD-462 RNA cosine 0.869 (CI 0.65–0.90) | 0.8690724689 (CI 0.6475–0.9028) | `osd462_anchor/results_summary.json::layer4_rna_gate` |
| OSD-513 ISS-T cosine 0.641 | 0.6408444531 | `run_20260518_201823_2500g/.../cross_osdr_alignment_summary.tsv` |
| ECM protein mean −0.107, p=0.012, concordance −0.40 | −0.1067, 0.0122, −0.3964 | `osd462_anchor/results_summary.json::layer1_protein` |
| NCC total protein +0.089 | +0.0888960 | same |
| KSEA SPAK/OSR1 z=−6.31 (p=2.8e−10), WNK z=−4.12 (p=3.7e−5) | exact match | `run_20260522_regulator_activity/osd462_kinase_activity_summary.tsv` |
| NCC Thr53 −0.851, Thr65 −0.790, Thr68 −0.794, Thr89 −0.930 | −0.8509, −0.7904, −0.7938, −0.9303 | `codex_review_osd462/.../phospho_axis_summary.tsv` |
| Composite Ser65;Thr68 = −1.563, p=9.3e−7 | −1.5633, 9.3e−7 | same |
| SPAK Ser366 −0.520; Ser383 −0.793 | −0.5202; −0.7930 | same |
| Endothelial vs NCC reg-phospho ρ=−0.762, p=0.0004 | matched in baseline lock | `v11_baseline_lock_summary.tsv` |
| DCT1 top-decile OR 1.51, q=1.13e−11 | OR 1.5053, q=1.126e−11 | `h2_enrichment/h2_dct1_sensitivity_summary.tsv` (primary_p05) |
| Anchor-excluded OR 1.49, q=3.68e−11 | OR 1.4893, q=3.68e−11 | same |
| Single-site-only OR 1.38, q=2.84e−6 | OR 1.3777, q=2.84e−6 | same |
| Full M4 OR 1.30, q=0.00158 | OR 1.2989, q=0.00158 | `h2_composition_adjusted_suppression_enrichment_single.tsv` |
| Mediation indirect (endothelial) −0.432, [−1.019,−0.036] | exact match | `h3_mediation_verdict.json` |
| Spatial niches (day14 inj-tub 0.827; DCT-adj 0.825; week6 DCT-adj 0.814; etc.) | 0.8271/0.8249/0.8137/0.8024/0.8019 | `spatial_reference/visium_niche_cosines.tsv` |
| DCT1 reference: 18,881/6,274 nuclei; Slc12a3 5.14/4.40; Pvalb 0.279/0.024; Trpv5 0.006/0.711; Calb1 1.01/2.97 | exact match | `external_qc/gse228367_marker_qc.tsv` |

The internal consistency between the manuscript, the companion compendium
(`results_v11.tex`), the v11 execution report (`docs/v11_execution_results.md`),
and the run-root TSV/JSON is high. SHA-256 checksums of all v11 inputs are
recorded in `manifests/v11_core_manifest.json` and the plan checksum
`d8b7465e09257d…` is locked in the compendium.

---

## 2. Statistical robustness — issues to address

### 2.1 OSD-513 bootstrap CI: 0.351 vs 0.373 mismatch (Results §3.1)
Manuscript: "cosine 0.641, bootstrap 95% CI 0.351 to 0.814."  
Repo: `0.6408 / CI 0.3733 to 0.8111` (`cross_osdr_alignment_summary.tsv`).
The point estimate matches exactly; the lower CI bound differs by ~0.02. Either
update the manuscript to 0.373–0.811, or, if 0.351 came from a different
bootstrap iteration/seed, point to the exact file. As written the CI cannot be
reproduced from the canonical artifact.

### 2.2 KSEA n=3 substrates per kinase (Results §3.2)
`osd462_kinase_activity_summary.tsv` shows SPAK_OSR1 z=−6.31 and WNK z=−4.12
were each computed against `n_substrates_quantified = 3`. The z-scores are tiny
because the background (n=21,083 sites) shrinks the standard error, but the
biological evidence is three sites per kinase (NCC T53/T65/T68 for SPAK; Stk39
S382/S383 + Oxsr1 S339 for WNK, with Stk39 S382 itself p=0.076 and only n=5
quantified channels). The manuscript already calls KSEA "a pathway-level
coherence check, not an unbiased kinome-wide discovery screen" (§2.7), but the
abstract still presents the z-scores as if they were primary discoveries. Add
the n_substrates value next to each z-score in Results (and Table S in the
compendium) so a reader is not misled by the apparent significance of the p
values.

### 2.3 Mediation model is an "approximate Bayesian OLS weak-prior fallback"
`h3_mediation_verdict.json::model` literally reads
`"approximate_bayesian_ols_weak_prior_fallback"`. This is a fallback because
the full Bayesian fit did not run on n=20. The Methods section calls this
"simple linear path models and bootstrap or approximate OLS uncertainty,"
which is accurate but understated; the manuscript should explicitly say that
the intended `brms`/SEM fit was not used and that intervals are approximate
posterior summaries from a normal-approximation fallback. The numerical
take-home — interval-crosses-zero for stromal and DCT, interval-excludes-zero
for endothelial and matrix-composite — is unaffected.

### 2.4 Phosphosite row dependence is acknowledged but not modeled
The Methods correctly note phosphosite rows are partially dependent through
shared parent genes, and the second-stage effect-level regression clusters SEs
by parent gene. The primary Fisher exact test for DCT1 enrichment, however,
treats sites as independent. The single-site-only sensitivity (one site per
parent gene) is the cleanest counterargument and already shows the top-decile
OR 1.38 survives. Recommend promoting the single-site-only number into the
abstract/Results paragraph (currently only top-quartile attenuation is
mentioned), and stating in the Discussion that the row-amplification
sensitivity primarily reduces effective n rather than effect size.

### 2.5 "Adjusted suppressed sites" are partially redefined per model
Each model M0–M5 defines its "suppressed" set by `β̂_F<0 & nominal p<0.05` in
that model's residual scale. This means the suppressed set size changes
across the ladder (M0: 2,430 → M4: 1,390 sites, see
`h2_composition_adjusted_suppression_enrichment_single.tsv`). The OR is a
fair comparison within a model, but the cross-model OR trend (1.40 → 1.55 →
1.36 → 1.24 → 1.30 → 1.38) partly reflects shifting denominators. The text
treats the ladder as a robustness check; explicitly say "set membership
changes by model" so reviewers do not read it as a single-set sensitivity.

### 2.6 Multiple testing scope is family-defined
BH is applied "within declared families". For the DCT1 enrichment family this
is well-defined and matches the compendium (single Fisher family with anchor,
NCC, and single-site sensitivities). The site-level BH for OSD-462
phosphosites (n=21,083) is separate. The manuscript should state once, near
the modelspec table, that there is no global FWER correction across families,
and list the families explicitly. The current Methods cover this implicitly
but it is buried.

### 2.7 KSEA q-value derivation
The KSEA table in the regulator-activity run shows `q_value = 5.5e-10` for SPAK
and `q = 3.7e-5` for WNK — i.e. q computed across two kinases only. With a BH
family of two, q ≈ p. This is honest in the file but not stated in the
manuscript; consider clarifying so reviewers understand q≈p here is not
fortuitous.

### 2.8 OSD-462 RNA gate used 11 pathways; cross-cohort vector used 9
`layer4_rna_gate.n_pathways = 11`, while §2.3 says "the primary recurrence
comparison used nine shared pathway features" (apoptosis, calcium handling,
DCT/NCC-WNK transport, ECM remodeling, fibrosis, inflammation, ion transport,
lipid metabolism, oxidative stress). The 11-pathway OSD-462 vector adds two
features (likely matrix/adhesion and endothelial/stromal context). Add one
sentence to §2.3 making this difference explicit.

---

## 3. Biological claims — points to tighten

### 3.1 NCC site residue letters in Table 2 are inconsistent
Table 2 lists "NCC Thr65" as a single site (effect −0.790) and
"NCC Ser65;Thr68" as a composite (effect −1.563). The residue letter at
position 65 cannot simultaneously be Thr (single) and Ser (composite). The
underlying `phospho_all_sites.tsv` does not carry a residue letter column —
only `site_position`. The residue letters in the manuscript Table appear to
have been added by hand. In mouse SLC12A3 (UniProt Q9DD25, long isoform), the
canonical SPAK/OSR1-activating cluster is T53/T58/S71 (rat T53/T58/T60); the
sites in OSD-462 are quantified at positions 53, 65, 68, 89 in whatever
isoform/protein-DB the proteomics core used. Fix in Table 2: either drop the
residue prefix (use "p53", "p65", "p68", "p89"), or add a footnote stating the
exact UniProt accession and isoform used and the residue at each position
under that accession. The composite label should be self-consistent with the
single-site labels.

### 3.2 Thr89 as "regulatory" is not sourced in the repo's curated substrate
file
`data/external/kinase_substrate/renal_kinase_substrate_core.tsv` lists the
canonical SPAK/OSR1 cluster as positions 53, 58, 65, 68 (citing
Pacheco-Alvarez 2006 and Richardson & Alessi 2008). Position 89 is included
as a "regulatory" site in Table 2 but does not appear in that curated list,
and is not standard in the WNK–SPAK–NCC literature. Either add a citation for
T89 as a regulatory site (some papers discuss adjacent N-terminal sites) or
relabel as "additional N-terminal phosphosite". The headline finding does not
depend on T89.

### 3.3 SPAK Ser382 has only n=5 (one-plex) and p=0.076
The compendium's WNK kinase output (z=−4.12) relies on Stk39 S382 (n_fl=5,
n_gc=5, p=0.076 — i.e. nominally non-significant), Stk39 S383 (n=10, p=2.8e−4),
and Oxsr1 S339 (n=10, p=0.87 in the spotted file). With one of three substrates
non-significant and another effectively flat in the per-site OLS, the WNK
z-score is dominated by Stk39 S383. The KSEA framing is fair, but Discussion
§4.2 should note this asymmetry rather than presenting WNK and SPAK output
suppression as symmetric findings.

### 3.4 Negative NCC anti-alignment vs PXD001729 is genuinely useful as a
boundary, but the abstract should not lean on it
The abstract currently quotes the result indirectly via "consistent with
reduced NCC activation." The compendium correctly states the dDAVP
anti-alignment cosine is −0.026 (CI −0.310 to 0.287) and Slc12a3/NCC has zero
shared transport target sites with OSD-462. That is properly framed as a
plausibility/coverage boundary, but the inability to test the vasopressin
anti-alignment hypothesis is a real limitation of the public-data layer and
deserves one sentence in the abstract or, better, the limitations section.

### 3.5 GSE228367 Welch test uses n=3 replicates per subtype
`gse228367_dct_prior_summary.json::statistical_caveat` says: "Welch tests use
three NK replicate-level mean-expression profiles per subtype; this is a
reference-prior analysis, not spaceflight evidence." The manuscript Methods
§2.5 says "DCT1 and DCT2 nuclei from normal-potassium controls were aggregated
by subtype and replicate" but does not state the replicate count. Add "n=3
biological replicates per subtype" so the strict-FDR DCT1-core count of 2
genes is interpretable as a power issue rather than a marker quality issue.

### 3.6 Tabula Muris Senis "did not map cleanly onto a simple accelerated
kidney-aging axis"
This claim is stated in Table 4 (Manuscript) but the supporting external
aging-axis projection is not included in the v11 run-root index. It exists in
prior runs (`tests/test_external_aging_axis.py` and earlier results). Either
cite the prior-run artifact in the compendium or include the v11 redo. As-is,
the claim is supported by older runs only.

---

## 4. Figures, plots, and tables — issues

### 4.1 v11_dct1_parent_gene_enrichment.pdf — "Single-site only" row is empty
Visual inspection of the forest plot shows the "Single-site only" panel is
rendered with no point estimate or CI for either DCT1 top-decile or top-quartile,
yet the data exist in
`h2_dct1_sensitivity_summary.tsv` (single_sites_only / fisher_dct1_top_decile
OR=1.38, q=2.84e−6; top_quartile OR=1.06, q=0.35). Fix in
`src/v11/publication_figures.py`: the panel layer is dropping the
single_sites_only family. This is the figure cited in Figure 4(a) of the
manuscript, so the reader sees a hole where the strongest within-row
sensitivity result should be plotted.

### 4.2 fig1_main_result_multipanel.pdf — Panel D "Full cosine = 0.55"
annotation does not match Results §3.1's "cosine 0.641"
The leave-one-out bar chart shows leave-one-pathway-out cosines around
0.40–0.70, but the legend dashed line is annotated "Full cosine = 0.55"
while the manuscript text and underlying TSV give the full
RRRM-2 ISS-T vs OSD-513 ISS-T cosine as 0.641. Either the figure is plotting
a different baseline (a different sub-vector or a centred version) or the
label is stale from an earlier run. Reconcile so the dashed line equals the
quoted 0.641 and label it accordingly.

### 4.3 fig_osd462_multiomics_dashboard.pdf is otherwise consistent
Panel A reports genome-wide Spearman 0.043 / Pearson 0.030 (matches JSON).
Panel B shows Slc12a3 +0.089, Stk39/SPAK small positive — matches. Panel C
shows the suppressed Slc12a3 53/65/68/89 sites and composites consistent with
Table 2. Panel D shows network-candidate translation is non-significant. No
changes needed.

### 4.4 Tables 2, 3, 4, S1–S6 reproduce the underlying numerics
All values cross-checked above. The Table-2 residue-letter issue is the only
substantive problem.

### 4.5 Figure 1 caption claims "live animal return is best treated as
recovery context"
This is supported by the data — OSD-513 LAR cosine is −0.511 (signed
anti-alignment) and the LAR arm did not pass the permutation null. The
caption is fine; only Panel D needs the relabel above.

---

## 5. Interpretability and reproducibility

### 5.1 Seeds and version pinning
- All v11 modules use `SEED = 20260526` (`core_analysis.py`,
  `h2_composition_aware_phospho.py`). Per-test sub-seeds are deterministic
  derivatives (`seed=20260526 + idx`, `seed=20260626 + idx`). Good.
- `manifests/v11_core_manifest.json` records SHA-256 of all inputs.
- `baseline/v11_baseline_input_manifest.json` also records SHA-256.
- The plan-checksum lock (`d8b7465e09257dd4…`) appears in both the compendium
  and the run-root, which lets reviewers verify the plan that was executed.
- `environment.yml`, `requirements.txt`, `install_r_packages.R` exist for
  reproducibility.

### 5.2 Pipeline entrypoint runs as advertised
`src/run_all_phases.py::phase_v11` references all the v11 modules. The
compendium states the dry-run succeeded for
`python -m src.run_all_phases --v11-only --dry-run --run-id run_20260526_v11_dryrun --v11-skip-spatial`.
The Xenium download checksum (`e6f5728ab7e47b499149f63b54e892c7`) matches the
figshare record per the compendium.

### 5.3 No v11-specific unit tests
`tests/` has 23 phase-level test files (`test_osd462_anchor.py`,
`test_regulator_activity.py`, etc.) that cover the upstream phases the v11
analysis depends on, but no `tests/v11/` directory and no direct tests of the
DCT1 enrichment Fisher logic, the composition-adjusted OLS, or the mediation
fallback. Recommend adding (a) a unit test that re-derives the top-decile
contingency table from a small fixture, (b) a snapshot test that the v11
figure module reproduces all panels (would have caught §4.1), and (c) a
deterministic mediation test against a synthetic fixture.

### 5.4 Manuscript code-availability statement is accurate
The "Data and Code Availability" block points to `src/` and `src/v11/` and the
correct entrypoint command. OSDR accessions (OSD-771, OSD-513, OSD-253,
OSD-462) and external accessions (GSE228367, PXD001729, GSE269622, GSE269719,
Tabula Muris Senis) are correctly listed.

### 5.5 Versioning of the manuscript itself
`latex_paper/` contains v2 through v11 plus several `manuscript_*_supplement*`
TeX files. The compile order is documented through the file names but not in
a README within `latex_paper/`. For reviewer reproducibility, a one-line note
in `latex_paper/README.md` indicating that `manuscript_v11.tex` +
`results_v11.tex` are the current pair would help.

---

## 6. Headline conclusions of the audit

1. The numerical claims in `manuscript_v11.tex` reproduce from the run-root
   artifacts to ~3 decimal places. The framing discipline is unusually
   strict — parent-gene enrichment, not cell localization; spatial
   contextualization, not validation; mediation as structure, not proof.
2. The most consequential correction is the **Table 2 residue-letter
   inconsistency** (Thr65 single vs Ser65 in the composite). This is a
   user-facing claim about the protein chemistry and should be fixed before
   submission.
3. The next most consequential is the **v11 DCT1 enrichment figure** dropping
   the single-site-only row. This is the strongest "not just NCC anchor"
   evidence and currently appears as a blank line in the rendered PDF.
4. The **OSD-513 lower-CI bound** (0.351 vs 0.373) and the **Figure 1D "Full
   cosine = 0.55"** label should both be reconciled to the canonical TSV.
5. The **KSEA z-scores** and the **mediation fallback** should each get one
   added sentence stating the n_substrates and the fallback nature; the
   numerical claims do not change.
6. **No v11 unit tests** is a real reproducibility gap; the catch in §4.1
   would have been automatic with a figure-snapshot test.
7. Reference accuracy: site-numbering and Welch-replicate caveats should be
   stated once in Methods to remove any chance of a reviewer reading the
   sub-1e−10 p-values as standalone evidence.

None of the above changes the headline result: a recurrent matrix/endothelial
RNA state with non-propagating protein abundance, sharply suppressed NCC
regulatory phosphorylation, and a DCT1-high parent-gene enrichment that
survives parent-protein and bulk-composition adjustment in the top decile.
