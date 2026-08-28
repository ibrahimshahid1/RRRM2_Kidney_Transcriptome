# Prospective DCT2/CNT–ASDN phosphoproteomic reanalysis

**Status:** Frozen replacement analysis specification  
**Frozen:** 2026-07-28  
**Starting repository commit:** `1d165c4`  
**Historical status:** Written after exploratory v11/v12 results were known  
**Controlling config:** `config/dct_asdn_phospho_reanalysis.yaml`

## 1. Primary question

After excluding prespecified canonical NCC/DCT1 machinery, are
flight-associated parent-gene phosphosite effects disproportionately shifted
toward lower phosphorylation among independently defined DCT2/CNT-transition or
aldosterone-sensitive distal-nephron (ASDN) programs?

The biological object is a set of **parent proteins associated with an external
cell-expression program**. Whole-kidney phosphoproteomics cannot assign a
phosphosite to a DCT2 or CNT cell, so neither the analysis nor the manuscript may
claim cell-of-origin localization.

## 2. Experimental estimand and unavoidable design limitation

The experimental contrast is ten space-flight versus ten ground-control WT
female kidneys from OSD-462/RR-10, split evenly across two TMTpro plexes.

The primary phosphosite effect is the plex-adjusted flight-minus-ground
difference in log2 official scaled reporter signal after channel median
centering. Uncentered official-scaled signal and channel-centered unscaled
signal-to-noise sums are prespecified sensitivities. The primary effect is a
**relative, site-specific phosphorylation estimand**. It is not an estimate of
absolute global phosphorylation.

Condition occupies the same reporter-channel ranges in both plexes. Reporter
position is therefore aliased with condition. No downstream model or label
permutation can remove that physical design limitation. The limitation must be
reported, and a strong biological headline requires concordance across the
predeclared normalization sensitivities and an orthogonal validation rationale.

## 3. Gate 0: sample and phosphopeptide provenance

No biological enrichment result may be promoted unless all of the following are
emitted and pass:

1. Exact sample, animal, condition, genotype, sex, plex, reporter, raw-file, and
   inclusion table.
2. Reconciliation of legacy `iTRAQ` metadata with the detailed TMTpro protocol.
3. NCC/SPAK peptide table containing accession/isoform, modified sequence,
   residue identity, site position, localization score, quantified-spectrum
   count, reporter completeness, and overlap with composite features.
4. Canonical regulatory-site annotation based on residue identity and literature,
   not position alone.
5. Separation of NCC Y65/Y68 features from the canonical
   SPAK/OSR1-regulated Ser/Thr set.
6. Deduplication of overlapping single/composite phosphoforms.

The primary global analysis uses uniquely localized, singly modified
phosphoforms with localization score at least 19. A localization-score threshold
of 13 and a deduplicated multi-modified analysis are sensitivities.

Filtering is label-blind: at least 16 of 20 finite reporter values and at least
6 of 10 finite values in each plex. Filtering is frozen before label
permutation.

## 4. Flight-blind subtype and compartment references

### 4.1 Discovery reference

Use the three normal-potassium GSE228367 biological replicates. Sum raw counts
within each biological replicate and subtype. Fit paired DCT2-versus-DCT1
pseudobulk differential expression with edgeR quasi-likelihood or an equivalent
count model.

Record log2 fold change, FDR, counts per million, detection fraction, and
direction in every biological-replicate pair.

### 4.2 Gene-set definitions

- **DCT1 core:** DCT1>DCT2 log2FC at least 1, FDR at most 0.05, detection at
  least 20%, same direction in all three pairs, and passing the frozen combined
  distal-specificity/breadth exclusion.
- **Strict DCT2-peaked:** DCT2>DCT1 log2FC at least 1, FDR at most 0.05,
  detection at least 20%, same direction in all three pairs, externally
  directionally validated, and passing the frozen combined
  distal-specificity/breadth exclusion.
- **DCT2/CNT transition:** DCT2>DCT1 log2FC at least 0.5, same direction in all
  three pairs, DCT2 detection at least 20%, retained or increased in CNT in the
  validation reference, and passing the frozen combined
  distal-specificity/breadth exclusion.
- **ASDN:** the literature-defined set frozen in the controlling config. It is
  not called DCT2-specific.

Strict-set FDR thresholds are not relaxed if a set is small or empty. The
transition set uses a predeclared effect/consistency/validation rule because
DCT2 is biologically transitional and the discovery reference has only three
replicates.

### 4.3 Independent validation and broad-expression flag

Validate subtype direction in GSE150338 or an equivalent independent
distal-nephron reference. Use a whole-kidney atlas to distinguish DCT2-peaked,
DCT2/CNT-shared, ASDN-shared, truly multi-compartment broad genes, and genes
that fail the stricter distal-specificity screen.

The validation reference defines cell-expression identity only. It is not
another spaceflight result.

Build comparator sets for proximal tubule, thick ascending limb, DCT1,
DCT2/CNT, principal cells, intercalated cells, podocytes, endothelium,
fibroblasts, and immune cells using the same flight-blind rules.

## 5. Site and parent-gene effects

For each eligible unique phosphoform, calculate the flight-minus-ground
contrast within each plex and average the two estimable plex contrasts with
equal weight. With complete balanced data this is equivalent to the flight
coefficient from `log2 reporter signal ~ flight + plex`; the explicit
within-plex form keeps the exchangeability structure visible when a site has
missing reporters.

Positive parent-gene scores mean lower phosphorylation in flight:

- **Primary:** median of `-beta_site` within each parent gene.
- **Sensitivity 1:** arithmetic mean of `-beta_site`.
- **Sensitivity 2:** one-sided maxmean, the mean of
  `max(-beta_site, 0)`.

No primary score uses a nominal-p membership gate or selects the lowest-p site.

For each parent gene, standardize the observed score against that gene's own
animal-label-permutation distribution. This calibrates site count, missingness,
and within-gene site structure.

## 6. Experimental null and set statistic

Permute five flight and five ground labels within each plex, using the same
permuted assignment for every phosphosite. There are
`choose(10,5)^2 = 63,504` balanced assignments, so the primary analysis uses
exact enumeration.

For every assignment:

1. Refit all site effects.
2. Recompute parent-gene scores.
3. Recompute gene-specific standardized scores.
4. Recompute every set statistic and exclusion tier.

The primary competitive statistic is:

`mean(Z_gene in set) - mean(Z_gene in eligible non-set background)`.

The self-contained set mean is secondary. Empirical p values use the exact
permutation distribution with the observed assignment included.

The two primary sets, DCT2/CNT-transition and ASDN, form one max-T
family-wise-error-controlled family. Kidney-compartment comparisons form a
separate prespecified family with BH correction.

## 7. Canonical-axis exclusions

The primary "beyond the canonical axis" analysis removes:

- `Slc12a3`
- `Stk39`
- `Oxsr1`
- `Wnk4`

A stricter sensitivity additionally removes `Wnk1` and all manually selected
NCC/SPAK anchor phosphoforms.

The broad-expression sensitivity excludes only target-expressed genes detected
across at least the configured number of unrelated kidney compartments. The
stricter combined distal-specificity/breadth flag remains a reference-signature
construction rule and is not substituted for this sensitivity. No post-result
manual scaffold blacklist is allowed.

## 8. Robustness and claim gates

A DCT2/CNT or ASDN claim passes only if:

1. At least eight phosphoproteome-observable genes remain after the primary
   canonical-axis exclusion.
2. The competitive effect is positive.
3. The exact max-T-adjusted permutation p is at most 0.05.
4. Direction is unchanged in the raw-signal and official-scaled-signal
   normalization sensitivities.
5. Direction remains positive after broad-expression exclusion.
6. Direction remains positive in every leave-one-gene-out analysis.
7. An independently derived subtype signature gives a concordant direction.
8. The result is not equally or more strongly explained by unrelated kidney
   compartment sets.

The strict DCT2/CNT title requires both the transition set and the independent
validation signature to pass. If only ASDN passes, the title and claim must say
ASDN rather than DCT2. If neither passes, the DCT2/CNT extension is rejected.

## 9. Secondary analyses

The following are supportive only:

- Nominal-p contingency tables.
- Matched annotation-label permutation.
- Motif-based KSEA.
- Parent-protein subtraction.
- Composition adjustment.
- RNA recurrence and RNA–protein concordance.

KSEA must be described as motif compatibility from the same phosphosite data,
not independent validation. Lower phosphorylation outside functionally verified
sites must not be equated automatically with lower parent-protein activity.

## 10. Required outputs

Every run writes:

- Input, config, code, and output hashes.
- Git commit and dirty-state summary.
- Complete sample/channel table.
- Peptide/residue provenance table.
- Reference-gene tables and validation results.
- Eligible-site and parent-gene score tables.
- Exact permutation results.
- Exclusion, normalization, broad-expression, and leave-one-gene-out results.
- Kidney-compartment forest-plot data.
- A machine-readable claim-decision file.

The manuscript is rewritten only after the claim-decision file is produced.
