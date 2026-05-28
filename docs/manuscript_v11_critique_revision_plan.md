# Manuscript v11 Critique And Revision Plan

Reviewed artifact: `latex_paper/manuscript_v11.pdf`, using the synchronized
source `latex_paper/manuscript_v11.tex`, the companion compendium
`latex_paper/results_v11.tex`, and the v11 run artifacts under
`data/results/run_20260526_v11_dct1_phospho_mediation/`.

This is a fresh critique of the current v11 manuscript. Several flaws from the
older critique have been fixed: the primary H2 code now uses a directional
Fisher test, a one-row-per-parent-gene sensitivity exists, the DCT2 comparator
is reported, KSEA substrate counts are disclosed, NCC residue labels have been
softened to position-only labels, and there is now at least one v11-specific
test file. The remaining issues are more about claim hierarchy, biological
interpretation, statistical independence, and figure/readability polish.

## Overall Verdict

The paper has a coherent observational multi-omic story and much better claim
discipline than many public-data manuscripts. Its strongest contribution is not
discovering NCC dephosphorylation, which prior work already established, but
placing that endpoint into a recurrent cross-cohort RNA context and showing that
the OSD-462 transporter signal is sharper at the phosphorylation layer than at
the protein-abundance layer.

The central DCT1 result is interesting but still fragile as a headline. The
manuscript now correctly admits that the DCT2-bottom-decile comparator is also
enriched and is stronger than DCT1 in the two most row-dependence-robust
summaries. That admission is scientifically honest, but it creates a title and
abstract tension: the safest claim is distal-nephron subtype-prior enrichment,
not a DCT1-centered result. The current title still leads with DCT1.

My bottom-line recommendation is to revise the paper around this hierarchy:

1. Primary result: matched OSD-462 shows RNA-protein decoupling and NCC/SPAK/WNK
   regulatory-phosphosite suppression.
2. Primary enrichment result: suppressed whole-kidney phosphosite parent genes
   are enriched in distal-nephron subtype-prior bins, including DCT1-high and
   DCT2-leaning extremes.
3. DCT1-specific emphasis should be secondary and biologically motivated by the
   NCC/SPAK/WNK axis, not presented as exclusive or dominant after the DCT2
   comparator.
4. Spatial, path, low-K, dDAVP, and KLHL3/CUL3 analyses should remain
   hypothesis-generating boundary analyses.

## Major Strengths

### 1. The biological spine is plausible

The paper connects five observations into a coherent model:

- public mouse spaceflight kidney RNA repeatedly shows matrix/endothelial-high
  and DCT/NCC-WNK-low context;
- OSD-462 RNA recurs the RRRM-2 direction strongly;
- protein abundance does not mirror the RNA pattern cleanly;
- NCC/SPAK/WNK regulatory phosphosites are suppressed despite flat total NCC;
- parent genes of suppressed phosphosite rows are enriched in distal nephron
  subtype-prior bins from an external DCT reference.

That is a real story, and it is more useful than a generic differential
expression catalog.

### 2. Claim boundaries are mostly explicit

The manuscript repeatedly says that OSD-462 is whole-kidney phosphoproteomics,
that the DCT1 prior is transcriptomic, and that the result is parent-gene
prioritization rather than cell-of-origin localization. This language should be
kept.

### 3. Negative and bounded results improve credibility

The paper preserves non-supportive findings: continuous DCT1-gradient models
are weak, parent-gene logistic adjustment attenuates the top-decile result,
low-K anti-alignment is unstable, PXD001729 cannot test the key transport sites,
KLHL3/CUL3 turnover is not resolvable, and spatial data are external IRI
context. This discipline makes the positive results more believable.

## Major Flaws And Risks

### 1. The title overstates the DCT1 interpretation

The title is currently:

`DCT1-high parent-gene enrichment of flight-suppressed phosphosites in mouse spaceflight kidney`

That is technically true for one subset, but it no longer reflects the paper's
own safest interpretation. The abstract and Results state that the DCT2-bottom
decile is also enriched and is stronger than DCT1 in the one-row-per-parent-gene
and parent-gene Fisher analyses:

- one row per parent gene: DCT2-bottom OR 1.68 versus DCT1-top OR 1.49;
- parent-gene Fisher: DCT2-bottom OR 1.82 versus DCT1-top OR 1.52.

This means a skeptical reviewer will ask why the title names DCT1 rather than
distal-nephron subtype-prior enrichment or DCT-subtype extremes.

Fix:

- Retitle around "distal-nephron subtype-prior parent-gene enrichment" or
  "DCT-subtype-prior enrichment".
- Keep DCT1 in the subtitle or final clause because NCC/SPAK/WNK are
  DCT1-leaning in the reference.
- Do not let the first sentence of the abstract sound like DCT1 exclusivity.

Suggested title:

`Distal-nephron subtype-prior enrichment of flight-suppressed phosphosites in mouse spaceflight kidney`

Alternative:

`RNA-protein decoupling and distal-nephron phosphosite-prior enrichment in mouse spaceflight kidney`

### 2. The DCT1 result remains statistically vulnerable to observability and parent-gene structure

The row-level DCT1 top-decile Fisher result is strong, but phosphosite rows are
not independent biological units. The manuscript has added safeguards:

- composite/multi-position rows excluded;
- one representative row per parent gene;
- parent-gene Fisher;
- site-count-stratified permutation;
- parent-protein and compartment-score models.

Those help. The main unresolved weakness is that the covariate-adjusted
parent-gene logistic model for any nominally suppressed site is not supportive
for DCT1 top decile:

- DCT1 top decile logistic OR 1.14, q = 0.453.

The analogous DCT2-bottom-decile parent-gene logistic result is also only
borderline:

- DCT2 bottom decile logistic OR 1.27, q = 0.099.

This does not erase the enrichment, but it means the strongest independent-unit
model does not support the headline as strongly as the row and Fisher summaries
do.

Fix:

- Promote the logistic attenuation into the abstract or a boxed limitations
  sentence if the DCT1 title is kept.
- Add a parent-gene cluster bootstrap for the row-level odds ratio.
- Add a matched parent-gene permutation preserving quantified site count,
  peptide count, and parent abundance, not site count alone.
- Report whether the DCT1 and DCT2 extreme-bin effects survive in the same
  multivariable parent-gene model.

### 3. "One site per parent gene" is still slightly overphrased

The code selects one representative phosphosite row per gene using lowest
phosphosite p value, then more negative effect, then site identifier
(`src/v11/core_analysis.py`). This is a good row-dependence sensitivity. But it
is not literally "one single phosphosite site" if the selected row is a
composite/multi-position row. The manuscript's Methods are mostly accurate
("one representative row per gene"), but the abstract calls it a "true
one-site-per-parent-gene sensitivity."

Fix:

- Rename this everywhere to "one phosphosite row per parent gene".
- Add a stricter sensitivity: single-position rows only, then one
  representative row per parent gene.
- If that stricter sensitivity remains positive, it is the number that should
  appear in the abstract.

### 4. The DCT1/DCT2 reference prior is underpowered and threshold-dependent

GSE228367 is a useful native DCT reference, but the analysis rests on only three
normal-potassium replicates per subtype. The strict FDR marker inventory is
extremely sparse: two DCT1-core genes and zero DCT2-core genes in the compendium.
That supports the choice to use percentile bins, but it also means the
percentile threshold is doing a lot of work.

Fix:

- Add leave-one-replicate-out DCT1/DCT2 priors.
- Add score-definition sensitivity: mean difference, log2 ratio, rank-averaged
  score, and detection-aware score.
- Add threshold sensitivity across top/bottom 5%, 10%, 15%, 20%, and 25%.
- Show whether NCC/SPAK/WNK parent genes themselves sit in the claimed bins or
  are only DCT1-leaning outside the strict bin.

### 5. The primary suppressed-site definition uses nominal p < 0.05

The manuscript correctly says nominal p < 0.05 defines a directional enrichment
set, not individual-site significance. Still, the main enrichment depends on a
thresholded p-value set. This can select for high intensity, lower missingness,
more observable proteins, and multi-site parent proteins.

Fix:

- Add threshold-free analyses:
  - rank-based enrichment of DCT1/DCT2 prior versus signed phosphosite effect;
  - logistic or ordinal model using signed effect rank instead of nominal p;
  - competitive gene-set test at parent-gene level.
- Add effect-size-only sensitivity, for example negative effect below the 25th
  percentile, independent of p value.
- Report missingness and mean intensity distributions by DCT1/DCT2 bin.

### 6. Composition adjustment is useful, but not deconfounding

The composition-aware ladder is valuable, but it uses bulk RNA marker scores as
animal-level covariates in n=20 animals. These covariates may be highly
correlated with flight and may sit on the biological path. The suppressed-site
set also changes by model:

- M0 raw: 2,430 nominal suppressed sites;
- M4 full model: 1,390 nominal suppressed sites.

The OR ladder is therefore not tracking one fixed set through adjustment; it is
redefining the suppressed set on each residual scale.

Fix:

- State in Results, not just Methods, that set membership changes across M0-M5.
- Add a fixed-set analysis: take the M0 suppressed set and ask whether adjusted
  effects remain negative or enriched.
- Add covariate correlation/VIF diagnostics for flight, DCT identity,
  endothelial, stromal, and parent-protein abundance.
- Avoid phrases like "not explained by composition"; use "not eliminated by
  this bulk-composition sensitivity."

### 7. Parent-protein normalization is not stoichiometry

The parent-protein-normalized effect is:

`phosphosite flight effect - parent protein flight effect`

This is a useful robustness check, especially because NCC and SPAK parent
protein effects are slightly positive. But the uncertainty for the
phosphosite-minus-protein contrast is not propagated, and the protein and
phosphosite estimates come from different TMT tables with different coverage and
missingness. It is not direct biochemical occupancy.

Fix:

- Keep the manuscript's current "not direct stoichiometry" caveat.
- Rename all remaining "occupancy" language to "parent-protein-normalized".
- Add a bootstrap or animal-level paired model for the contrast where the same
  animals have both parent protein and phosphosite values.

### 8. KSEA is a coherence check, not strong kinase-output evidence

The manuscript now reports that targeted KSEA uses three curated substrates per
kinase. Good. But the z scores still look visually and rhetorically stronger
than the substrate evidence supports. WNK is especially asymmetric because the
score is driven by a tiny curated set, with individual substrate behavior that
is not uniformly strong.

Fix:

- Move KSEA z scores out of the abstract or keep only "targeted KSEA over three
  substrates per kinase was directionally supportive."
- Add a small substrate table beside the KSEA result with per-site effect, p,
  q, and number of finite channels.
- Treat NCC phosphosite suppression as the main phosphoproteomic fact; treat WNK
  KSEA as pathway coherence.

### 9. RNA recurrence is pathway-level and curated, not genome-wide validation

The RNA recurrence result is one of the manuscript's better pieces, but it is
based on curated gene-set/pathway scores and cosine similarity over small
feature vectors. Gene sets overlap biologically, and the OSD-513 recurrence uses
nine shared features while OSD-462 uses an 11-pathway anchor panel. The current
caption notes this, but the Methods should make the feature-set difference and
curation freeze more transparent.

Fix:

- Put gene-set member lists in a supplement table.
- Add a leave-one-gene-set-family-out analysis, not only leave-one-pathway-out.
- Bootstrap both the RRRM-2 reference and the external cohort, not only the
  external cohort, or clearly state that RRRM-2 uncertainty is treated as fixed.
- Add a concise genome-wide rank correlation or camera/fgsea sensitivity if
  feasible.

### 10. TMT protein and phosphoprotein models are simple relative-abundance models

The OSD-462 TMT analyses use log2 scaled signal-to-noise values with
median-centering and simple within-plex or OLS effects. This is reasonable for a
public reanalysis, but it has limitations:

- median centering can remove true global abundance shifts;
- missingness may be non-random;
- no empirical Bayes moderation is used for site-level phosphoproteomics;
- channel/plex effects and peptide/protein aggregation are simplified;
- phosphosite localization confidence and accession-specific residue mapping are
  not central in the current tables.

Fix:

- Add TMT QC panels: channel medians before/after normalization, missingness by
  group, plex balance, and sample PCA.
- Add a limma/MSstatsTMT-style sensitivity for protein and phosphosite effects
  if feasible.
- Report localization/ambiguity metadata for the NCC/SPAK/WNK sites.
- Keep "relative TMT abundance" wording throughout.

### 11. Spatial analysis is a prediction generator with pseudoreplication

The IRI Visium analysis is useful, but the spot-level t tests are not
animal-level replicated evidence. The code explicitly labels this as
"spot-level descriptive test; no animal-level spatial replication". The
manuscript states this caveat, but it still quotes p values prominently.

Fix:

- Move spot-level p values to the supplement.
- In the main text, report effect direction and magnitude without inferential
  language.
- If metadata permit, aggregate by section/animal before any testing.
- Change "found lower DCT transport" in the abstract to "observed a lower
  spot-level DCT-transport score in an external IRI reference."

### 12. Path analysis should be called covariance decomposition

The path model is cross-sectional, n=20, bulk-tissue, and uses an approximate
Bayesian OLS fallback rather than a full SEM/brms model. The Methods say this,
and the verdict JSON says it clearly. The Results still risk sounding like a
mediation result because of "indirect path" language.

Fix:

- Rename the section "Exploratory covariance decomposition".
- Keep the forest plot in the supplement unless the main text needs it.
- In the main text, say only that remodeling scores covary negatively with NCC
  regulatory phosphorylation after conditioning choices.
- Remove any causal verbs around endothelial or matrix remodeling.

### 13. Low-K and dDAVP branches are too large for their evidentiary weight

The low-K result is mechanistically attractive but supported only in a focused
14-gene subset with wide bootstrap intervals crossing zero. Genome-wide and
DCT-prior subsets are near zero or inconsistent. The dDAVP/mpkDCT comparison has
60 shared single sites but zero shared transport-target sites.

Fix:

- Compress low-K, dDAVP, and KLHL3/CUL3 into one "mechanism triage did not
  resolve upstream signal" paragraph in the main text.
- Keep detailed heatmaps/tables in the compendium.
- Avoid "potassium/chloride-like" as a named mechanism unless future targeted
  data support it.

### 14. The figure set needs cleanup before review

Current visual issues:

- `v11_dct1_parent_gene_enrichment.png`: the legend covers the strict-q red
  estimate/interval.
- `v11_parent_protein_composition_sensitivity.png`: legends cover M5 estimates
  in both panels.
- `v11_perturbation_triangulation.png`: the layout has excessive vertical white
  space, and panel D is dominated by the DCT-marker-high spike; the subtle
  DCT-adjacent decline is hard to see.
- `fig_osd462_multiomics_dashboard.png`: panel B annotation is cramped near the
  right edge.

Fix:

- Move legends outside plotting areas or into unused margins.
- Split perturbation Figure 5 into two figures or move the low-K/dDAVP material
  to supplement.
- For spatial DCT transport, use separate y axes or small multiples for
  DCT-marker-high and DCT-adjacent spots.
- Add a figure QA test that fails if key plotted rows have no visible artist or
  if a legend overlaps data bounds.

### 15. Reproducibility is better, but v11 test coverage is incomplete

There is now a v11-specific test file, `tests/test_v11_h2_enrichment.py`, which
checks Fisher directionality, the representative-row reducer, parent-gene table
collapse, and figure label handling. That is a meaningful improvement.

Remaining gaps:

- no direct tests for composition-adjusted OLS;
- no test for parent-protein-normalized contrast uncertainty or labels;
- no test for mediation fallback determinism;
- no visual snapshot or figure-completeness tests;
- no test that manuscript headline numbers match TSV/JSON sources.

Fix:

- Add v11 tests for composition adjustment, parent-normalized effects, spatial
  spot summary labels, and mediation fallback output.
- Add a headline-number reconciliation test that parses a small YAML/TSV table
  of expected manuscript numbers.
- Add `latex_paper/README.md` identifying `manuscript_v11.tex` and
  `results_v11.tex` as the current pair.

## Results-Specific Critique

### RNA recurrence

Useful and plausible, but it should be described as curated pathway-level
recurrence. It does not establish a universal spaceflight kidney transcriptomic
signature. Add stronger documentation of gene-set membership and uncertainty in
the RRRM-2 reference vector.

### RNA-protein mismatch

This is one of the paper's most persuasive results. Keep it central. Be careful
that "matrix proteins moved opposite their transcripts" may reflect TMT
coverage, ECM extractability, and relative abundance normalization, not only
biology.

### NCC/SPAK/WNK phosphosite suppression

The NCC regulatory-phosphosite cluster is the cleanest molecular anchor. Total
NCC flatness makes it more convincing. KSEA should remain secondary and targeted.
The p89 site is now appropriately softened as an additional N-terminal
phosphosite; keep it out of the curated regulatory-site core unless cited.

### DCT1/DCT2 parent-gene enrichment

This is the main novelty, but the current data support "distal-nephron
subtype-prior subset enrichment" better than "DCT1-centered enrichment." The
DCT1 top-decile result is real enough to keep, but the DCT2-bottom-decile
comparator must force a title and abstract reset.

### Composition and parent-protein robustness

These analyses are valuable sensitivity checks. They should be framed as "not
eliminated by these adjustments" rather than "independent of parent protein or
composition." The top-quartile attenuation and weak continuous-gradient models
are important boundaries.

### Public perturbation references

These analyses make the paper more honest by showing what public data cannot
resolve. The main text currently gives them more space than their evidentiary
weight deserves. Keep parent-protein-normalized enrichment in main; compress
low-K, dDAVP, KLHL3/CUL3, and spatial IRI into boundaries and predictions.

## Prioritized Revision Plan

### Phase 1: Reset the central claim

1. Change the title from DCT1-high enrichment to distal-nephron subtype-prior
   enrichment, or add a subtitle that immediately states DCT2-bottom-decile
   enrichment.
2. Rewrite the first abstract sentence so DCT1 is not the sole headline.
3. Make DCT2 comparator results part of the primary Results paragraph, not a
   boundary afterthought.
4. Change "true one-site-per-parent-gene" to "one phosphosite row per parent
   gene" unless a single-position-only representative analysis is added.

### Phase 2: Strengthen the core enrichment analysis

1. Add single-position-only plus one-row-per-parent-gene sensitivity.
2. Add parent-gene cluster bootstrap for DCT1 and DCT2 extreme-bin ORs.
3. Add multivariable parent-gene models that include DCT1-top and DCT2-bottom
   bins together.
4. Add matched parent-gene permutation preserving site count, peptide count,
   parent abundance, and missingness.
5. Add threshold-free signed-rank enrichment analyses.

### Phase 3: Harden the DCT reference prior

1. Add leave-one-replicate-out GSE228367 DCT1/DCT2 priors.
2. Add DCT1/DCT2 score-definition sensitivity.
3. Add percentile-threshold sensitivity across 5%, 10%, 15%, 20%, and 25%.
4. Report where NCC/SPAK/WNK anchor parent genes fall under each prior.

### Phase 4: Rebalance methods and statistics

1. Add a compact table of testing families with number of tests and correction
   scope.
2. State in Results that M0-M5 adjusted suppressed sets differ in size.
3. Add TMT QC and missingness summaries.
4. Keep KSEA as targeted coherence and move exact z/p values out of the
   abstract if space is tight.
5. Rename mediation/path results as covariance decomposition.

### Phase 5: Clean figures and tables

1. Move legends out of data regions in Figures 3-5.
2. Redesign perturbation Figure 5 or move most of it to supplement.
3. Add a DCT1-versus-DCT2 comparator forest plot as the key enrichment figure.
4. Add a small table showing row-level, one-row-per-gene, parent-gene Fisher,
   logistic, and permutation results side by side.
5. Add figure QA tests.

### Phase 6: Improve reproducibility

1. Add tests for composition adjustment, parent-normalized effects, spatial
   summaries, and mediation fallback.
2. Add a headline-number reconciliation test against run-root TSV/JSON files.
3. Add `latex_paper/README.md` documenting the current manuscript/compendium
   pair and build command.
4. Add a supplement artifact index mapping every headline number to its source
   file and column.

## Suggested Revised Claim Hierarchy

Primary claim:

> Matched OSD-462 whole-kidney multi-omics show that the spaceflight
> DCT/NCC-WNK transporter lesion resolves most clearly at the regulatory
> phosphosite layer rather than at protein abundance, inside a recurrent
> matrix/endothelial-high and DCT-transport-low RNA context.

Primary enrichment claim:

> Flight-suppressed whole-kidney phosphosite parent genes are enriched in
> distal-nephron subtype-prior extremes from an external DCT1/DCT2 snRNA-seq
> reference, including a DCT1-high top-decile subset and a DCT2-leaning
> bottom-decile comparator.

Secondary DCT1 claim:

> The DCT1-high top-decile result is biologically motivated by NCC/SPAK/WNK
> being DCT1-leaning in the reference, but it is not exclusive and does not
> localize phosphosite changes to DCT1 cells.

Exploratory claims:

> Bulk compartment scores, low-K DCT references, IRI spatial references,
> dDAVP/mpkDCT phosphoproteomics, and KLHL3/CUL3 coverage checks generate
> hypotheses and boundaries; they do not identify the upstream WNK suppressor,
> prove mediation, or spatially validate the spaceflight lesion.

Forbidden claims:

- DCT1-specific phosphoproteomics.
- DCT1 cell-of-origin assignment from whole-kidney data.
- DCT1 exclusivity.
- Spatial validation of spaceflight kidney tissue.
- Causal mediation by endothelial remodeling.
- Resolved potassium/chloride, vasopressin, or KLHL3/CUL3 mechanism.
- Direct phosphosite occupancy or stoichiometry measurement.

## Bottom Line

The manuscript is worth continuing. The RNA-protein mismatch and NCC regulatory
phosphosite suppression are strong enough to anchor the story. The DCT1-high
enrichment is interesting, but after the DCT2 comparator it should no longer be
the sole headline. The best next revision is a claim reset plus deeper
parent-gene-aware and DCT-prior sensitivity analysis, not a broader search for
more public datasets.
