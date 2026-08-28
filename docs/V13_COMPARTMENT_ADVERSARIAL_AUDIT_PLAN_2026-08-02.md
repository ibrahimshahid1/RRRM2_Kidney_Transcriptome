# V13 kidney-compartment adversarial audit

**Frozen:** 2026-08-02  
**Status:** post-hoc adversarial audit; primary podocyte and mesenchymal/stromal
results were already known. This is not a preregistration and cannot convert an
exploratory result into confirmatory evidence.

## Question

After applying one analysis to all major kidney compartments, does any
flight-labelled-channel phosphosite signal survive marker-specificity,
observability, normalization, parent-protein, contributor-stability, and
protein-class controls?

The inferential object is an externally annotated **parent protein** in
whole-kidney phosphoproteomics. It is not a cell of origin, a clinical kidney
injury endpoint, or proof of a spaceflight effect. Condition is perfectly
aliased with reporter-tag blocks in both plexes.

## What is already complete and will not be recomputed unnecessarily

The definitive v13 run already provides 63,504 balanced within-plex label
assignments, three parent-gene statistics, site-quality profiles, uncentered
and summed-signal profiles, parent-protein subtraction, and complete
leave-one-gene-out inference for the original capped atlas sets. Existing
reporter-position and phospho-versus-protein block diagnostics are also
retained.

NCC/SPAK cannot be used as a publication-grade positive control: Stage 0 found
zero isolated qualified canonical phosphoforms. A moderated fourth gene score
is not added post hoc because gene-specific exact-null standardization already
calibrates the declared scores, and another estimator cannot repair the tag
alias.

## Frozen marker tiers

The Mouse Kidney Atlas is kidney-only. It can support *within-kidney
compartment specificity*, not organ-level “kidney restriction.” For each of
ten compartments, genes must first have target mean CPM >= 1, target detection
in at least half of source studies, and target expression at least twofold
above every other mapped kidney compartment.

Five sets are then emitted per compartment:

1. `all_enriched`: every gene passing the twofold rule, with no top-200 cap;
2. `within_kidney_not_broad`: fewer than four non-target compartments meet
   CPM >= 1 and source-study detection >= 0.5;
3. `high_specificity`: target/max-other >= fourfold and target study detection
   >= 0.75;
4. `broad_enriched`: at least four non-target compartments meet the expression
   and detection thresholds;
5. `scaffold_excluded`: `all_enriched` after removing a frozen union of KEGG
   adhesion, ECM, focal-adhesion, actin, tight-junction, and vascular
   smooth-muscle genes.

The atlas `DCT_general` label is reported as `DCT1_context`; it is not true
subtype-resolved DCT1. The `DCT2_CNT_atlas` label is reported as
`DCT2_CNT_context`. The atlas fibroblast/stroma/pericyte grouping is reported
as `mesenchymal_stromal`.

## One exact family, in both directions

All 50 compartment-by-tier sets plus the broad structural/scaffold control are
one family. The primary profile is the qualified, singly modified,
non-composite, localization-score >= 19 universe. The primary statistic is the
competitive mean standardized parent-gene suppression score. All 63,504
balanced assignments are enumerated.

The run reports upper-tail, lower-tail, one-sided maxT, and two-sided
max-|T| familywise values. The max-|T| result across all 51 sets is the frozen
multiplicity decision for the primary median score. The signed-mean and
one-sided-maxmean sensitivities repeat the same 51-set family separately; they
are not a joint correction across all three score definitions. A negative
statistic is not silently converted into a new “increased phosphorylation”
story.

## Secondary observability null

For each of the ten `all_enriched` compartment sets and the scaffold control,
annotation labels are shuffled within propensity strata derived from:

- quantified site count;
- median phosphopeptide signal;
- missingness;
- phosphopeptide count;
- parent-protein availability;
- parent-protein abundance; and
- parent-protein peptide count.

There are 10,000 matched draws. The eleven p-values form one BH family. This
secondary gene-label null asks whether an annotation is unusual among similarly
observable genes; it does not replace the animal-label permutation.

## Contributor and sensitivity requirements

For each compartment the audit emits every leave-one-gene-out statistic,
positive-contribution concentration, top contributors, scaffold overlap,
uncentered direction, summed-signal direction, and parent-protein-adjusted
direction. Grey60 overlap is descriptive and cannot validate a phosphosite
result unless the overlapping genes are themselves observable and shifted.

A candidate survives the post-hoc adversarial gate only if all of the
following hold:

1. `all_enriched` is positive with max-|T| FWER <= 0.05;
2. median, signed mean, and one-sided maxmean statistics are all positive;
3. every leave-one-gene-out statistic remains positive;
4. the high-specificity tier has at least eight observable genes and remains
   positive;
5. the seven-covariate observability-matched BH q <= 0.05;
6. uncentered and parent-protein-adjusted directions remain positive;
7. no gene supplies >10% and the top ten supply <=50% of total positive Z.

Passing this gate would justify only a statement about concentration among
compartment-enriched parent proteins. It cannot override reporter-block alias,
whole-kidney dilution, lack of cell localization, or absence of renal
physiology, urine, histology, GFR, albuminuria, electrolyte, or injury-marker
measurements.

## Clinical-assay translation

The user-proposed clinical assay framework is used only to build a descriptive
observability map. Tissue RNA/protein/phosphopeptide coverage is not equated
with the corresponding serum or urine assay. No cancer analysis is performed
because the repository has neither renal-mass imaging nor tumor tissue.

## Stop rule

If no compartment survives the frozen adversarial gate, compartment mining in
OSD-462 stops. Remaining credible directions must then be design/methods work,
matched-library robustness, or new data acquisition with urine, histology,
physiology, and segment-resolved validation—not another post-hoc compartment
narrative.
