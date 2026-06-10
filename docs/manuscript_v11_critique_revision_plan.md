# Manuscript v11 Critique And Revision Plan

Reviewed artifact: `latex_paper/manuscript_v11.pdf`, with the synchronized
source `latex_paper/manuscript_v11.tex`, the companion compendium
`latex_paper/results_v11.tex`, and v11 run artifacts under
`data/results/run_20260526_v11_dct1_phospho_mediation/`. This critique now
also incorporates the June 2 channel-centering QC artifacts:
`baseline/osd462_tmt_raw_channel_pattern_summary.tsv` and
`h2_enrichment/h2_dct_channel_centering_sensitivity.tsv`.

This critique is for the June 2 rewrite. It includes scientific, statistical,
methodological, figure, clarity, and prose-level issues.

## Bottom-Line Verdict

The rewrite is much more scientifically defensible than the earlier DCT1-heavy
version. The central claim is now closer to what the data can actually support:
flight-suppressed whole-kidney phosphosites are enriched among parent genes in
distal-nephron DCT1/DCT2 reference-prior extremes, with DCT1-high and
DCT2-leaning bins both involved. That is a publishable direction.

The main remaining problem is that the manuscript now reads partly like a
paper and partly like an internal analysis memo. It repeatedly tells the reader
what the claim is not, which checks were added, which branches were triaged,
which results were not promoted, and which appendix analyses are only context.
Those caveats are necessary, but they currently dominate the prose. The paper
needs to sound less like a defensive audit trail and more like a confident
scientific narrative with clearly separated primary, sensitivity, and appendix
evidence.

My recommended final hierarchy:

1. Primary biological contribution: OSD-462/RR-10 places the established
   NCC/SPAK/WNK dephosphorylation endpoint inside a recurrent remodeling RNA
   context, but the RNA context does not propagate cleanly to protein abundance.
2. Primary enrichment contribution: suppressed phosphosite parent genes are
   enriched in distal-nephron subtype-prior extremes, not exclusively DCT1.
3. DCT1 interpretation: biologically motivated by NCC/SPAK/WNK anchor genes,
   but weakened by matched parent-gene permutation and by stronger DCT2-leaning
   parent-gene-aware summaries.
4. Supporting extensions: kinome-wide KSEA recovery of WNK1 suppression and
   cross-cohort matrix/remodeling recurrence deserve real main-text space.
   Continuous DCT gradient, aldosterone, human urine, CMap, and IRI analyses
   are boundaries or hypotheses, not validation.

## Major Strengths

### 1. The central claim is now honest

The manuscript no longer tries to force a DCT1-only story. It directly reports
that DCT1 attenuates under matched parent-gene permutation (`p = 0.107`) while
the DCT2-leaning comparator remains supported (`p = 0.017`). That honesty makes
the paper more credible.

### 2. The matched multi-omic argument is strong

The clearest contribution is the RNA/protein/phosphoprotein mismatch. Flat NCC
protein with suppressed NCC/SPAK/WNK regulatory phosphosites is biologically
interesting and anchored in the right molecular layer.

### 3. The limitations are real, not decorative

The manuscript correctly states that whole-kidney phosphoproteomics cannot
assign cell of origin, that row-level phosphosite tests are not independent
biological units, that nominal `p < 0.05` defines an enrichment set rather than
site-level discovery, and that the public perturbation references do not solve
the upstream mechanism.

## Highest-Priority Problems

### 1. The abstract is too dense and too defensive

The abstract opens with a long block of odds ratios and q values. It contains
the right facts, but it reads like a statistical audit paragraph rather than a
scientific abstract. Phrases such as "the safest enrichment claim" and repeated
"rather than..." language are useful in a rebuttal, but they weaken the paper's
voice in the abstract.

Fix:

- Lead with the biological finding, not the sensitivity ladder.
- Keep only one DCT1 number, one DCT2 number, and the matched-permutation
  contrast.
- Move the full sensitivity ladder to Results/Table 3.
- Replace "safest claim" with "These results support..."

Suggested abstract spine:

> We reanalyzed public mouse kidney spaceflight multi-omics to ask how the
> established NCC/SPAK/WNK dephosphorylation phenotype relates to recurrent
> kidney remodeling and distal-nephron subtype programs. In matched RR-10/
> OSD-462 data, the transporter lesion was sharper at the phosphoproteomic
> layer than at protein abundance. Suppressed phosphosite parent genes were
> enriched in distal-nephron DCT1/DCT2 reference-prior extremes. DCT1-high
> genes were enriched in row and parent-gene analyses, but DCT1 support
> attenuated under matched parent-gene permutation, whereas the DCT2-leaning
> comparator remained supported. Thus, the data support distal-nephron
> subtype-prior enrichment anchored by NCC/SPAK/WNK biology, not DCT1
> cell-of-origin assignment.

### 2. The paper still rhetorically leans DCT1 while the robust signal is broader

The title is now appropriately broad, but the Results and Discussion still put
most of the biological weight on DCT1. The DCT2-leaning bin is not a side note:
it is stronger in several parent-gene-aware summaries:

- one row per parent gene: DCT2 OR 1.68 vs DCT1 OR 1.49;
- parent-gene Fisher: DCT2 OR 1.82 vs DCT1 OR 1.52;
- matched parent-gene permutation: DCT2 `p = 0.017` vs DCT1 `p = 0.107`;
- percentile sensitivity: DCT2 parent-gene enrichment remains strong through
  wider cutoffs, while DCT1 weakens at 20 percent.

Fix:

- Add a short biological paragraph explaining what the DCT2/CNT-leaning result
  could mean: DCT2/CNT calcium handling, aldosterone-sensitive/electrogenic
  distal transport, distal-nephron neighborhood effects, or a low-DCT1-score
  observability pattern. Keep magnesium/Trpm6 framed as DCT1-leaning.
- Do not treat DCT2 only as a "comparator"; call it a co-primary boundary of
  the distal-nephron subtype-prior result.
- Rename "DCT2-bottom-decile" in running prose. It is mathematically precise
  but confusing. Use "DCT2-leaning low-DCT1-score decile" once, then
  "DCT2-leaning decile."

### 3. The continuous and threshold-free evidence is not as clean as the bin evidence

The threshold-free table complicates the story. Row-level Spearman correlation
between DCT1 score and signed phosphosite effect is weak (`rho = -0.0085`,
`p = 0.226`). Parent-gene most-negative-site Spearman is positive
(`rho = 0.045`, `p = 0.001`), which is not supportive under the predicted
negative-gradient sign. The Mann-Whitney extreme-bin tests are supportive, but
the continuous gradient is not.

Fix:

- State clearly: "The evidence is strongest for extreme-bin enrichment, not
  for a continuous DCT1 gradient."
- Keep the now sign-aware `directional_read` label in
  `h2_dct_threshold_free_summary.tsv`, and explain in prose that the positive
  parent-gene Spearman coefficient is a boundary, not support.
- Avoid implying that higher DCT1 score broadly predicts stronger suppression.
  The data support top-decile subset enrichment.

### 4. Phosphosite TMT channel behavior needs explicit treatment

The raw-workbook QC artifacts show systematic lower phosphosite channel medians
for flight than ground control in both plexes, especially in the single-site and
composite-site phosphosite sheets. Protein channel medians are much more stable.
That is exactly the raw pattern a reviewer will worry about. However, the
primary phosphosite effect table feeding DCT-prior enrichment is channel
centered: `scripts/osd462/02_phospho_axis.py` explicitly calls
`compute_site_flight_effect_lm(..., channel_center=True)`, and the helper now
defaults to channel centering.

The empirical QC now answers both parts of the reviewer concern:

- Raw phosphosite scaled channel medians are lower in FL than GC: single-site
  sheets show FL-GC median shifts of about `-0.93` and `-0.59` across the two
  plexes; composite-site sheets show about `-0.97` and `-0.65`.
- Centering changes the global phosphosite effect baseline from a median
  uncentered effect of `-0.135` to a median centered effect of `+0.0068`, which
  is exactly what a TMT loading-like shift should do.
- Centered and uncentered site effects remain very highly rank-concordant
  (`Spearman rho = 0.998`), but the uncentered estimator calls far more sites
  nominally suppressed (`7893` vs `3134` across all compared sites), so centered
  estimation is the defensible primary estimator.
- The DCT-prior enrichment is not created by centering. The centered primary
  row-level test gives DCT1 top-decile OR `1.505`, `p = 4.74e-13`, and
  DCT2-leaning decile OR `1.293`, `p = 2.92e-6`. The uncentered sensitivity
  remains enriched, though attenuated at the row level: DCT1 OR `1.270`,
  `p = 3.94e-8`; DCT2 OR `1.169`, `p = 1.50e-4`.

Fix:

- Add a short QC paragraph and a supplement table using
  `baseline/osd462_tmt_raw_channel_pattern_summary.tsv`,
  `h2_enrichment/h2_dct_channel_centering_effect_comparison.tsv`, and
  `h2_enrichment/h2_dct_channel_centering_sensitivity.tsv`.
- Phrase the conclusion carefully: the raw FL-lower channel pattern exists,
  but the primary enrichment uses centered phosphosite effects and the
  enrichment survives the uncentered sensitivity. Therefore this is a QC issue
  to document, not a fatal artifact.

### 5. The recurrent RNA claim should emphasize matrix/remodeling as the positive result

The matrix/endothelial recurrence is stronger than the DCT/NCC-WNK transport
recurrence. The cross-cohort gene-set meta artifact reports matrix/ECM with
Stouffer `p = 0.00070`, while the DCT/NCC-WNK transport set is weak
(`p = 0.187`, median FDR about 0.948, median I2 about 73 percent). The main text
often says "matrix/endothelial-high and DCT-low" as if the two halves are
equally robust.

Fix:

- Give matrix/remodeling recurrence real main-text weight: it is one of the
  two strongest extension analyses and helps distinguish the paper from simply
  re-reporting Siew's NCC endpoint.
- Phrase the RNA result as "a recurrent matrix/endothelial remodeling context
  with directional DCT/NCC-WNK transport suppression" unless the specific
  pathway vector result is being discussed.
- Separate the recurrence evidence for matrix/remodeling from the recurrence
  evidence for DCT transport in one table.
- Do not let the RNA layer carry more DCT certainty than the phosphoproteomic
  layer.

### 6. Human urine concordance is too fragile for main-text "supportive" framing

The human section is disciplined in the artifacts, but it remains risky in the
paper. It uses one Twins flight subject, three physiological axes, four scored
directions, and one figure-level AQP2 direction mixed with machine-readable
Table S8 urine analytes. The sign test is non-significant (`p = 0.25`). Sodium,
urine volume, and magnesium are all physiologically confounded in spaceflight.

Fix:

- Move human concordance fully to appendix/supplement, or keep only one
  paragraph in the Discussion as "directional consistency with fluid-axis
  predictions."
- Do not call it validation.
- Do not include it in the abstract or conclusion except as a bounded appendix
  note.
- Keep AQP2 figure-level evidence separate from machine-readable Table S8
  analytes; do not make the reader feel a plotted protein direction and urine
  chemistry were analyzed at the same evidentiary level.

### 7. Aldosterone-axis evidence is directional context, not confirmation

The aldosterone-axis analysis is useful because it tests a biologically
plausible upstream explanation for distal-nephron transporter regulation. But
the result is directional rather than significant (`p = 0.06`). If the main
text treats it as a positive mechanistic result, reviewers can fairly object
that the paper is stacking near-misses around the phosphosite anchor.

Fix:

- Keep one sentence in Results or Discussion: aldosterone-axis evidence was
  directionally consistent but did not pass the prespecified significance
  threshold.
- Do not use aldosterone to explain the phosphoproteomic effect unless the
  wording is explicitly hypothesis-generating.
- Keep the detailed axis table in the supplement.

### 8. CMap is scope creep unless it is strictly appendix

The CMap screen is local, approximate, cell-line based, and currently uses an
uppercase mouse-to-human symbol assumption rather than a dedicated ortholog map.
It does not help the core kidney phosphoproteomics story enough to justify much
main-text space.

Fix:

- Remove CMap from the main conclusion.
- Keep a compact appendix table only after adding a real ortholog map.
- If no official CLUE query is run, label the result "local approximate
  connectivity screen" everywhere.
- Never mention treatments, countermeasures, or therapeutic implication.

### 9. Kinome-wide KSEA is a useful supporting positive, but not a full mechanism claim

The targeted KSEA is coherent but small: three quantified curated substrates
per kinase. The kinome-wide Johnson atlas screen is therefore valuable because
it recovers WNK1 suppression by an unbiased substrate-atlas route
(`z = -4.87`), orthogonal to the targeted NCC/SPAK/WNK anchor. That deserves
main-text mention. The boundary is that the atlas screen does not turn the
entire WNK-SPAK/OSR1-NCC pathway into an independently validated mechanism:
other members are weaker, missing, or not aligned enough to carry that claim.

Fix:

- Present WNK1 recovery in the main Results as an unbiased supporting positive.
- Keep targeted KSEA as pathway coherence, not independent kinase discovery.
- In Results, say the measured NCC/SPAK/WNK phosphosite rows are the primary
  evidence, with kinome-wide WNK1 suppression as orthogonal support.

### 10. Parent-protein-normalized enrichment is helpful but over-precise

The parent-protein-normalized analysis is the best robustness upgrade, but the
primary table still uses the raw phosphosite nominal p threshold. It subtracts
parent protein effect but does not fully propagate uncertainty in the
phosphosite-minus-protein contrast for all reported rows.

Fix:

- Describe it as "effect-normalized robustness," not occupancy.
- Put the paired animal-level phosphosite-minus-parent-protein model first if
  that is the more direct matched analysis.
- Avoid presenting `q = 1e-11` as if it resolves stoichiometry.

### 11. The compendium and figures are not synchronized with the manuscript rewrite

The manuscript was modified June 2, but `latex_paper/results_v11.tex` and the
main figures are May 29 artifacts. The compendium does not yet include the new
human concordance, CMap, aldosterone, or kinome-wide sections that the
manuscript now references, and it should now also include the channel-centering
QC artifacts. This creates a reproducibility and credibility gap.

Fix:

- Update `latex_paper/results_v11.tex` to include human concordance, CMap,
  aldosterone, kinome-wide KSEA, and channel-centering QC artifacts.
- Regenerate the figure set after the June 2 rewrite.
- Add a "current manuscript numbers" table mapping every headline number to
  the source TSV/JSON.
- Keep the repaired `analysis` column in
  `h2_dct_extreme_bin_cluster_bootstrap.tsv` so all-row and one-site-per-gene
  bootstrap rows remain distinguishable.

## Writing And Clarity Problems

### The prose sounds too internal

Terms like "Reach D", "Repair B", "branch", "promotion rule",
"claim discipline", "verdict", "RUN_ID", and "v11" are useful pipeline
language but should mostly stay out of the manuscript. They make the paper feel
like a project log.

Fix:

- In the paper: use "appendix screen", "sensitivity analysis", "predefined
  decision rule", "analysis module", or omit the label entirely.
- In Data and Code Availability: list the repo and analysis entrypoint, but do
  not dump many internal module names unless the journal requires it.
- Remove placeholder language such as "Authorship and collaboration structure
  will be revised with mentor and collaborator input" before circulation.

### The paper overuses defensive negation

The manuscript repeatedly says what analyses are "not": not validation, not
deconvolution, not discovery, not causal, not direct stoichiometry, not
treatment, not cell localization. These boundaries are correct, but the
repetition makes the paper sound insecure.

Fix:

- Put a single boundary sentence at the end of each relevant section.
- Use positive framing first: "This analysis prioritizes parent genes..." then
  one concise boundary: "...and does not localize phosphosites to DCT cells."
- Remove phrases like "real enough to motivate biology"; replace with
  "supports a DCT1-motivated but non-exclusive interpretation."

### The Results are too number-heavy

Many paragraphs are long strings of ORs, q values, CIs, and caveats. The reader
needs a stronger sentence telling them what the numbers mean.

Fix:

- Use one prose sentence for the finding, then table numbers.
- After each major table, add one interpretive sentence:
  "The pattern is therefore an extreme-bin enrichment, not a continuous DCT1
  gradient."
- Move sensitivity ladders out of paragraphs and into compact tables.

### The limitation section is a wall

The limitations paragraph is scientifically good but visually punishing. It is
too long for a reader to retain.

Fix:

Split it into four short paragraphs:

- public cohort/data limitations;
- phosphoproteomics/DCT-prior limitations;
- mechanism/spatial limitations;
- human/CMap appendix limitations.

### Some terminology is hard on readers

"DCT-subtype-prior parent-gene enrichment" is accurate but dense. It needs one
plain definition early:

> We use "DCT-subtype prior" to mean an external DCT1/DCT2 snRNA-seq ranking
> assigned to phosphosite parent genes; it is not a cell-of-origin call.

Then use shorter terms afterward: "DCT-prior enrichment", "DCT1-high genes",
and "DCT2-leaning genes."

## Novelty Framing

The paper should not try to become novel by inflating every extension into a
positive result. That would make the manuscript easier to attack. The novelty
is stronger if it is framed as a disciplined public multi-omics boundary map
around the known RR-10 kidney phenotype:

- Siew established the kidney transporter/phosphorylation phenotype. This
  manuscript adds the multi-layer reanalysis showing where that endpoint does
  and does not propagate across RNA, protein, phosphoprotein, public recurrence,
  kinase-substrate priors, and human-facing context.
- The two strongest extensions are the unbiased WNK1 KSEA recovery and the
  cross-cohort recurrence of the matrix/remodeling program.
- The bounded results are still useful: no broad continuous DCT gradient,
  aldosterone directionality without significance, human urine directionality
  without power, and CMap hypotheses only.

That is a more novel claim than "we re-found NCC suppression": it says which
parts of the public multi-omics ecosystem reproduce, which do not, and where
the next experiment should go.

## Figure Critique

### Main DCT-prior figure

The figure shows DCT1 top decile, DCT2-leaning decile, and DCT1 top quartile,
but it does not show the matched parent-gene permutation or logistic
attenuation that now governs the interpretation. The bottom note/legend layout
also collides visually.

Fix:

- Make the key figure a DCT1-vs-DCT2 evidence ladder:
  row Fisher, one-row-per-gene, parent-gene Fisher, logistic, matched
  permutation, threshold-free extreme-bin test.
- Drop DCT1 top quartile from the main panel or move it to supplement.
- Add direct labels for DCT1 and DCT2 odds ratios/p-values where possible.

### Parent-protein/composition figure

The figure is visually busy and still DCT1-centered. It does not help readers
understand the now broader distal-nephron result.

Fix:

- Separate parent-protein normalization from composition adjustment.
- Add DCT2-leaning comparator where relevant.
- Move the continuous-gradient panel to supplement unless it is explicitly used
  to show a negative result.

### Perturbation figure

This figure mixes low-K RNA anti-alignment, target-gene heatmaps,
parent-protein-normalized enrichment, and IRI spatial prediction. Those are
different evidence types. The figure reads like an analysis dashboard.

Fix:

- Split into two figures or demote to supplement.
- Main text can keep a small "public references did not identify a single
  upstream mechanism" table.
- Keep IRI spatial prediction as a focused figure only if it is central to the
  future-experiment pitch.

## Revision Plan

### Pass 1: Claim Hierarchy

1. Keep the broad title.
2. Rewrite the abstract around three findings: RNA/protein mismatch,
   phosphoproteomic NCC/SPAK/WNK lesion, distal-nephron subtype-prior
   enrichment.
3. Remove human/CMap from the abstract and downshift it in the conclusion.
4. Make DCT2-leaning enrichment a co-primary boundary of the DCT-prior result.

### Pass 2: Statistical Tightening

1. Add a phosphosite TMT QC paragraph reporting the new raw-channel and
   centered-vs-uncentered sensitivity artifacts.
2. Explain the threshold-free Spearman boundary: the row-level Spearman is
   weak and the parent-gene Spearman sign is not supportive, even though
   extreme-bin Mann-Whitney tests are supportive.
3. Add a clean parent-gene evidence table with DCT1 and DCT2 side by side:
   Fisher, logistic, cluster bootstrap, matched permutation, threshold-free,
   percentile sweep.
4. Rephrase RNA recurrence so matrix/endothelial recurrence is strongest and
   DCT transport recurrence is directional/contextual unless specifically
   supported by the pathway-vector contrast.
5. Clarify that parent-protein-normalized tests are robustness analyses, not
   stoichiometry or occupancy.

### Pass 3: Prose Cleanup

1. Remove internal labels from main text: "Reach", "Repair", "promotion rule",
   "v11", "verdict", "RUN_ID".
2. Replace repeated negative caveats with one boundary sentence per section.
3. Shorten Results paragraphs and push sensitivity ladders into tables.
4. Split the limitations section into readable paragraphs.
5. Replace "DCT2-bottom-decile" in prose with "DCT2-leaning decile" after a
   one-time definition.

### Pass 4: Figures And Supplement

1. Regenerate figures after the June 2 manuscript rewrite.
2. Redesign the main DCT-prior figure around the evidence ladder, including
   matched permutation and logistic attenuation.
3. Split or demote the perturbation figure.
4. Update `latex_paper/results_v11.tex` with human concordance, CMap,
   aldosterone, kinome-wide KSEA, channel-centering QC, and the current
   headline-number index.
5. Add a short supplement table showing all source artifact paths for headline
   numbers.

## Recommended Final Claim Language

Use:

> Flight-suppressed whole-kidney phosphosites are enriched among parent genes
> in distal-nephron DCT1/DCT2 reference-prior extremes. The DCT1-high subset is
> biologically anchored by NCC/SPAK/WNK genes, but the DCT2-leaning comparator
> is also enriched and is more stable in several parent-gene-aware summaries.
> The result therefore supports distal-nephron subtype-prior enrichment, not
> DCT1 exclusivity or DCT cell-of-origin assignment.

Avoid:

> DCT1 phosphosites are suppressed in spaceflight kidney.

> Human data validate the mouse DCT model.

> CMap identifies countermeasures.

> Composition adjustment rules out remodeling.

> Kinome-wide KSEA confirms the SPAK/OSR1 mechanism.

## Overall Priority

The paper does not need more novelty by overreaching. It needs a cleaner,
less defensive, better synchronized version of the story it now has. The
novelty is the disciplined integration: recurrent remodeling RNA context,
matched RNA/protein/phosphoprotein mismatch, and distal-nephron subtype-prior
phosphosite enrichment with explicit parent-gene and observability boundaries.
