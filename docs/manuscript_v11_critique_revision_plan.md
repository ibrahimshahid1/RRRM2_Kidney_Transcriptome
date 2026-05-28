# Manuscript v11 Critique And Revision Plan

Reviewed artifact: `latex_paper/manuscript_v11.pdf`, using the synchronized source
`latex_paper/manuscript_v11.tex` plus v11 result artifacts under
`data/results/run_20260526_v11_dct1_phospho_mediation/`.

## Overall Verdict

The manuscript is scientifically much more disciplined than earlier versions
appear to have been. Its strongest feature is claim hygiene: it repeatedly
distinguishes DCT1-high parent-gene enrichment from DCT1 cell-of-origin
localization, and it treats spatial and mediation analyses as hypothesis
generating. That framing is essential and should be preserved.

The paper is not yet submission-ready. The central DCT1 top-decile result is
interesting, but the current manuscript still has several reviewer-visible
weaknesses:

1. One key sensitivity analysis is mislabeled or overinterpreted.
2. The KSEA and phosphosite-site labels need stricter biochemical accounting.
3. The statistical inference still leans heavily on phosphosite rows that are
   not independent biological units.
4. Spatial and mediation results are carefully hedged in text, but the abstract
   and narrative still risk sounding stronger than the data allow.
5. A figure-generation bug leaves the single-site sensitivity blank in the main
   DCT1 enrichment plot.
6. Reproducibility is good at the manifest level but weak at the unit-test level
   for the new v11-specific methods.

## Major Scientific Strengths

### 1. The paper has a clear and plausible biological spine

The strongest story is:

- public mouse kidney spaceflight RNA shows recurrent matrix/endothelial-high
  and DCT/NCC-WNK-low context;
- matched OSD-462 protein abundance does not mirror the RNA signature;
- the transporter signal is sharpest at regulatory phosphosite level;
- suppressed whole-kidney phosphosites are enriched among parent genes that are
  DCT1-high in an external native DCT1/DCT2 reference;
- public perturbation references narrow, but do not solve, the upstream WNK
  suppressor question.

That is a coherent observational multi-omic paper.

### 2. The most important claim is appropriately hedged

The abstract and Results repeatedly state that the DCT1 signal is parent-gene
enrichment in whole-kidney phosphoproteomics, not DCT1-localized
phosphoproteomics. This is exactly the right boundary. Keep that language.

### 3. The negative results improve credibility

The manuscript preserves null or bounded results: protein concordance is null,
continuous DCT1-gradient models are weak, low-K anti-alignment has wide CIs,
PXD001729 cannot test transport-target anti-alignment, and KLHL3/CUL3 turnover
is not resolvable. This makes the positive DCT1 top-decile enrichment more
credible because the paper is not trying to turn every analysis into support.

## Major Flaws And Risks

### 1. Critical reporting mismatch: "single-site" is not one-site-per-gene

Manuscript Methods say the sensitivity analyses "retained one representative
phosphosite per parent gene" (`manuscript_v11.tex`, line 174). The Results cite
the top-decile OR 1.38, q=2.84e-6 as the single-site sensitivity.

But the v11 core implementation defines `single_sites_only` as rows where
`site_position_str` is a single numeric position, not one representative site
per parent gene:

- `src/v11/core_analysis.py`, line 328: `is_single_site` is a regex for a
  single numeric site.
- `src/v11/core_analysis.py`, line 451: the sensitivity uses
  `phospho_prior[phospho_prior["is_single_site"]]`.

That filter removes composite phosphosite rows; it does not solve parent-gene
row dependence. The mapped "single-site" primary sensitivity still has 16,514
rows but only 5,109 unique parent genes. Therefore the current text overstates
how well row dependence has been controlled for the primary DCT1 enrichment.

Important nuance: the occupancy-normalized branch does implement a true
one-site-per-gene reduction in `src/v11/h2_occupancy_normalized_phospho.py`
lines 94-100 and uses it at lines 134-135. The primary unadjusted branch does
not.

Fix:

- Rename the current primary sensitivity to "single-position sites only" or
  "composite sites excluded."
- Add a true one-site-per-parent-gene sensitivity to the primary H2 enrichment,
  using the same selection rule as the occupancy branch or a pre-specified
  rule such as lowest phosphosite p value per parent gene.
- Recompute and report the true one-site-per-gene OR/q in the main DCT1 result.
- Update Methods, Results, Figure 4, and the compendium.

### 2. The DCT1 enrichment is row-level; parent-gene and site-count effects are still the main statistical vulnerability

The headline DCT1 top-decile enrichment is strong at row level:

- primary top decile: OR 1.51, q=1.13e-11;
- anchor excluded: OR 1.49, q=3.68e-11;
- full adjusted sensitivity: OR 1.30, q=0.00158;
- occupancy-normalized: OR 1.52, q=3.73e-12.

But the biological unit is ambiguous. A phosphosite row is not independent of
other rows on the same parent protein, and DCT1-high genes may have different
phosphosite observability or multi-site density. The manuscript acknowledges
this, but the main evidence still relies on a row-level Fisher test.

Fix:

- Add a parent-gene-level logistic or permutation analysis:
  - define whether each parent gene has at least one suppressed site;
  - model DCT1 top-decile status against suppression with covariates for protein
    abundance, peptide count, number of quantified phosphosites, and missingness;
  - use Firth/logistic regression or permutation if separation occurs.
- Add a cluster bootstrap over parent genes for the row-level OR.
- Add a site-count-stratified permutation preserving the number of sites per
  parent gene.
- Make the row-level Fisher result secondary to a parent-gene-aware confirmation.

### 3. Fisher exact test direction is inconsistent across branches

In the primary H2 code, `stats.fisher_exact(arr)` is called without an
`alternative`, so SciPy uses the two-sided default (`src/v11/core_analysis.py`,
line 411). In the occupancy and composition branches, Fisher tests use
`alternative="greater"`. The manuscript generally frames the test as
directional enrichment.

This probably does not change the headline top-decile conclusion, but it is a
reproducibility and interpretation blemish.

Fix:

- Decide whether the DCT1 enrichment family is directional or two-sided.
- Use the same alternative in all H2, occupancy, and composition tests.
- State the choice explicitly in Methods and table footnotes.

### 4. The KSEA result looks more decisive than its substrate count supports

The abstract reports KSEA SPAK/OSR1 z=-6.31 and WNK z=-4.12. The underlying
table shows only three quantified substrates per kinase against a 21,083-site
background. This is a coherence check, not an unbiased kinome-wide discovery.

The issue is strongest for WNK: Stk39 S383 is strong, Stk39 S382 is incomplete
and nominally weaker, and Oxsr1 S339 is essentially flat in the per-site table.
The WNK score therefore should not be presented as symmetric with the NCC/SPAK
site cluster.

Fix:

- Add `n_substrates_quantified=3` next to each KSEA z score in Results and
  compendium.
- In the abstract, say "targeted KSEA over three curated substrates per kinase"
  or move the z scores out of the abstract.
- Discuss WNK as weaker/asymmetric compared with the NCC phosphosite cluster.

### 5. NCC phosphosite labels need biochemical cleanup

Table 2 lists NCC "Thr65" and the composite "Ser65;Thr68." The residue at
position 65 cannot be both Thr and Ser in the same accession/isoform. The
underlying OSD-462 site table carries positions, not residue letters.

Also, the curated renal kinase-substrate file lists SPAK/OSR1 NCC positions
53, 58, 65, and 68, but not 89. The manuscript labels Thr89 as "regulatory,"
which needs a specific citation or a softer label.

Fix:

- Replace residue-letter labels with position-only labels such as NCC p53,
  p65, p68, p89 unless the exact accession and residue mapping are verified.
- Footnote the accession/isoform and site-numbering source.
- Relabel p89 as "additional N-terminal phosphosite" unless a regulatory
  citation is added.

### 6. OSD-513 recurrence numbers are internally inconsistent across artifacts

The manuscript reports OSD-513 cosine 0.641 with CI 0.351 to 0.814
(`manuscript_v11.tex`, line 257). The canonical 2,000-bootstrap summary in
`run_20260518_201823_2500g` reports 0.373 to 0.811. A separate cosine-permutation
run reports 0.351 to 0.814. Both appear to exist, but the manuscript does not
identify which run is canonical.

Fix:

- Pick one run as canonical for manuscript v11.
- Cite the exact artifact in the compendium.
- If the permutation-null run is canonical, say so because it explains the
  0.351 CI.

### 7. Figure 1D label is stale or computed from a different vector

The text reports full RRRM-2 ISS-T vs OSD-513 cosine 0.641, but Figure 1D labels
"Full cosine = 0.55." The figure may use a centered, filtered, or old vector,
but the caption does not explain this difference.

Fix:

- Regenerate Figure 1D so the dashed line matches the text, or relabel it as
  the exact statistic it represents.
- Add a small figure-note if panel D uses a different feature set from the main
  0.641 recurrence estimate.

### 8. The main DCT1 enrichment figure drops the single-site row

The rendered `v11_dct1_parent_gene_enrichment.png` leaves the "Single-site only"
row blank. The data exist under the label `single_sites_only`, but the plotting
code expects `single_site_only_p05`:

- data label: `single_sites_only`;
- plotting label expected at `src/v11/publication_figures.py`, line 178:
  `single_site_only_p05`.

Fix:

- Change the plotting code to the actual label or standardize the TSV label.
- After adding a true one-site-per-gene test, plot both "composites excluded"
  and "one site per gene" or just the latter.

### 9. Occupancy normalization is a useful robustness check, but not true occupancy

The occupancy branch computes:

`occupancy_effect = phosphosite flight effect - parent protein flight effect`.

The p-value threshold still comes from the raw phosphosite model
(`src/v11/h2_occupancy_normalized_phospho.py`, lines 116-124). That is fine as
a robustness check, but it is not a direct phospho-stoichiometry model with its
own uncertainty for the phospho-minus-protein contrast.

Fix:

- Rename as "parent-protein-normalized phosphosite effect" throughout; use
  "occupancy" only with a caveat.
- Compute uncertainty for the contrast if feasible by animal-level joint
  modeling or bootstrap.
- Do not imply biochemical site occupancy was measured directly.

### 10. Spatial analysis is visually persuasive but statistically fragile

The IRI Visium result is useful as a prediction generator. It is not validation.
The code compares spot-level distributions with Welch t tests
(`src/v11/spatial_dct_transport_check.py`, lines 115-129), but spots are not
independent biological replicates and the data are from IRI, not spaceflight.

The manuscript states this limitation, but the abstract gives the spatial
numbers enough prominence that readers may overread them.

Fix:

- Move spatial p values out of the abstract or explicitly say "spot-level,
  external IRI reference."
- Add animal/section-level aggregation if sample metadata allow it.
- Replace "predict" with "motivates the prediction" unless a true
  out-of-sample validation design is added.

### 11. Low-K comparison is mechanistically attractive but underpowered and gene-set-selected

The low-K result is directional only in a hand-focused 14-gene target subset,
with wide bootstrap intervals crossing zero. The genome-wide and DCT-prior
subsets are near zero. This is correctly not promoted to a primary mechanism,
but the paper still spends a lot of Results space on it.

Fix:

- Keep the low-K analysis in a compact "mechanism triage" paragraph.
- Add an explicit null sentence: genome-wide and DCT-prior subsets did not show
  stable anti-alignment.
- Avoid "potassium/chloride-like" in the title, abstract, or conclusion.

### 12. Mediation/path analysis should be demoted further unless the model details are explicit

The path analysis is cross-sectional, n=20, bulk-tissue, and uses an approximate
Bayesian OLS fallback rather than the intended full SEM/brms model. The text is
mostly careful, but Results currently says "consistent with negative indirect
paths" without telling the reader about the fallback.

Fix:

- State in Methods and Results that the full Bayesian/SEM fit was not used and
  the intervals are normal-approximation fallback summaries.
- Move the mediation forest to supplement unless it is essential to the main
  narrative.
- Frame as "covariance decomposition" rather than "mediation" in the main text.

### 13. Multiple-testing families are not explicit enough

The manuscript says BH was applied within defined families, but readers need to
know those families without searching the compendium. There is no global FWER
or global FDR across all exploratory analyses.

Fix:

- Add a short table listing each testing family, number of tests, correction
  scope, and whether it is primary, sensitivity, or exploratory.
- State that no global correction was applied across all families.

### 14. DCT1 reference prior needs more uncertainty and specificity analysis

The DCT1/DCT2 prior uses only three normal-potassium replicates per subtype.
Strict FDR marker sets are sparse. Percentile bins are reasonable, but the
paper should prove the DCT1 top-decile result is not an artifact of one scoring
definition.

Fix:

- Add sensitivity to:
  - DCT1-DCT2 difference score;
  - log2 ratio score;
  - rank-averaged score;
  - replicate leave-one-out DCT1 prior;
  - DCT1 top 5%, 10%, 15%, 20% thresholds.
- Add an explicit DCT2-high or CNT-associated comparator.

### 15. Reproducibility is strong for manifests but weak for tests

The v11 pipeline records manifests and checksums, which is good. But there are
no direct tests for the new v11 DCT1 enrichment logic, composition adjustment,
occupancy normalization, mediation fallback, or figure completeness.

Fix:

- Add tests for:
  - Fisher contingency-table construction on a tiny fixture;
  - true one-site-per-gene reduction;
  - consistency of all Fisher alternatives;
  - figure row completeness, which would have caught the blank single-site row;
  - occupancy effect calculation and caveat;
  - mediation fallback determinism.
- Add a `latex_paper/README.md` note saying `manuscript_v11.tex` and
  `results_v11.tex` are the current manuscript/compendium pair.

## Results-Specific Critique

### RNA recurrence

The recurrence result is plausible and useful, but it should be positioned as
pathway-level recurrence, not a universal or genome-wide spaceflight kidney
signature. OSD-462 uses 11 pathways in the anchor summary, while the cross-cohort
OSD-513 recurrence uses nine. The difference needs one sentence in Methods.

### RNA-protein mismatch

This is one of the paper's best results. The null protein-abundance layer
clarifies why the phosphoproteome matters. However, "matrix proteins moved
opposite their transcripts" should be interpreted cautiously because TMT protein
coverage, extracellular matrix extractability, and whole-kidney cellular
composition can all affect observed protein direction.

### NCC/SPAK/WNK phosphosite suppression

The NCC phosphosite cluster is the cleanest molecular anchor, especially because
total NCC protein is flat. This can stay central. The WNK KSEA interpretation
should be softened because the WNK substrate count is tiny and one substrate is
flat.

### DCT1-high parent-gene enrichment

This is the primary novel result and is worth keeping in the title. But after
fixing the single-site-per-gene issue, the title should remain precise:

Recommended title:

`DCT1-high parent-gene enrichment of flight-suppressed phosphosites in mouse spaceflight kidney`

Avoid:

`DCT1 phosphoproteomic suppression in spaceflight kidney`

### Public perturbation references

These analyses are useful boundaries, not mechanisms. The parent-protein-normalized
DCT1 enrichment belongs in the main text. Low-K, dDAVP, KLHL3/CUL3, and IRI
spatial context can be shorter in the main text and fuller in the compendium.

## Prioritized Revision Plan

### Phase 1: Fix correctness issues before any prose polishing

1. Implement a true primary one-site-per-parent-gene H2 sensitivity.
2. Rename the current `single_sites_only` output to "single-position sites only"
   or "composite sites excluded."
3. Standardize Fisher exact test alternatives across H2, composition, and
   occupancy branches.
4. Fix the Figure 4 blank row caused by the `single_site_only_p05` vs
   `single_sites_only` label mismatch.
5. Reconcile the OSD-513 CI and Figure 1D full-cosine label.
6. Correct NCC residue/site labels and the p89 regulatory wording.

### Phase 2: Strengthen the central DCT1 claim

1. Add parent-gene-level suppression models.
2. Add parent-gene cluster bootstrap or site-count-preserving permutations.
3. Add DCT1 prior score-definition sensitivity.
4. Add a DCT2/CNT comparator analysis.
5. Promote the parent-gene-aware result into the abstract only if it remains
   supportive.

### Phase 3: Tighten statistical reporting

1. Add a multiple-testing family table.
2. Add KSEA substrate counts and q-value family size.
3. Add explicit sample sizes for DCT reference pseudobulk replicates.
4. State that occupancy-normalized effects are not direct stoichiometry.
5. State that mediation uses an approximate OLS fallback.

### Phase 4: Rebalance the manuscript narrative

1. Keep the main text focused on four results:
   - RNA recurrence;
   - RNA-protein mismatch;
   - NCC/SPAK/WNK regulatory phosphosite suppression;
   - DCT1-high parent-gene enrichment with robust sensitivities.
2. Move most perturbation, spatial, live-return, aging, TLR4, and network
   material to the compendium unless it directly protects the main claim.
3. Shorten the abstract's spatial and low-K details.
4. Make the Discussion's first paragraph state the claim hierarchy plainly.

### Phase 5: Add reviewer-facing reproducibility support

1. Add v11-specific tests.
2. Update `latex_paper/README.md` to identify current manuscript and supplement.
3. Add a small artifact index table in the supplement that names the exact TSV
   or JSON source for every headline number.
4. Add a rendered-figure QA checklist.

## Suggested Claim Hierarchy After Revision

Primary claim:

> In matched OSD-462 whole-kidney phosphoproteomics, flight-suppressed
> phosphosites are enriched among parent genes that are DCT1-high in an external
> DCT1/DCT2 snRNA-seq reference, especially in the top decile, and this subset
> enrichment survives parent-protein and bulk-compartment sensitivity analyses.

Secondary claims:

> Across public mouse spaceflight kidney cohorts, terminal flight shows a
> pathway-level matrix/endothelial-high and DCT/NCC-WNK-low RNA context.

> In OSD-462, that RNA context does not propagate cleanly to protein abundance,
> while NCC/SPAK/WNK regulatory phosphosites are suppressed despite flat total
> NCC protein.

Exploratory claims:

> Bulk compartment scores, low-K references, IRI spatial references, dDAVP
> phosphoproteomics, and KLHL3/CUL3 checks generate mechanistic and spatial
> hypotheses but do not identify the upstream WNK suppressor or localize the
> lesion to DCT1 cells.

Forbidden claims:

- DCT1-specific phosphoproteomics.
- Cell-of-origin assignment from whole-kidney phosphoproteomics.
- Spatial validation of the spaceflight lesion.
- Causal mediation by endothelial remodeling.
- Vasopressin or KLHL3/CUL3 mechanism resolved by public data.
- Direct phosphosite occupancy/stochiometry measurement.

## Bottom Line

The manuscript's core idea is worth preserving: an established spaceflight
NCC/SPAK/WNK phosphoproteomic lesion sits inside a recurrent remodeling RNA
context and is enriched in a DCT1-high parent-gene subset. The most urgent
revision is not more analysis breadth; it is making the central enrichment
claim parent-gene-aware, fixing the mislabeled single-site sensitivity, and
cleaning up biochemical/statistical reporting so the paper can survive a
skeptical methods review.
