# V13 reanalysis outcome and manuscript decision

**Date:** 2026-07-29  
**Frozen question:** After excluding the established NCC/DCT1 axis, is
flight-associated phosphosite suppression disproportionately concentrated in
independently defined DCT2/CNT-transition or aldosterone-sensitive
distal-nephron programs?

## Decision

The reanalysis does **not** support the v12 claim that phosphorylation
suppression extends beyond NCC/DCT1 into DCT2/CNT or the
aldosterone-sensitive distal nephron (ASDN).

- The DCT2/CNT-transition hypothesis is **non-evaluable**, not null: only 5 of
  27 frozen genes are observable after the primary exclusions, below the
  prespecified minimum of 8. The independently derived validation signature is
  similarly under-covered (5 of 17).
- The ASDN set has a positive channel-centered statistic, but it **fails the
  prespecified selectivity criterion**. A podocyte annotation has a larger
  statistic, and the secondary observability-matched gene-annotation test is
  not significant. No direct ASDN-versus-podocyte contrast was tested.
- The machine-readable statistical tier is therefore **`neither`**: the
  DCT2/CNT component is non-evaluable and the ASDN component fails.
- Biological publication promotion is separately blocked by perfect
  condition-to-reporter-tag-block aliasing and by the absence of an isolated,
  assay-qualified canonical NCC/SPAK phosphoform.

The appropriate product from this branch is a compact
provenance/statistical-boundary report. It is not a biology manuscript with a
“beyond NCC/DCT1,” DCT2/CNT, ASDN, canonical-site, or transporter-activity
headline.

## What was repaired

### 1. Assay and sample provenance

The exact OSD-462 MS design is now documented at the biological-sample and
reporter-channel level.

- The workbooks contain 30 wild-type female B6129SF2/J left-kidney samples in
  two 15-plexes: 5 baseline, 5 flight, and 5 ground-control samples per plex.
- The primary comparison is 10 flight versus 10 ground-control animals.
- The detailed protocol specifies TMTpro labeling, a +304.207-Da tag mass, and
  TMTpro isotope-impurity correction. A legacy repository field instead says
  iTRAQ. The assay is reported as **TMTpro isobaric-tag
  proteomics/phosphoproteomics, with a disclosed legacy metadata
  inconsistency**.
- Each plex assigns baseline to reporter positions 1--5, flight to 6--10, and
  ground control to 11--15. The pattern is repeated without a cross-plex label
  swap. Flight status therefore cannot be separated from a systematic
  reporter-position or tag-block effect.

### 2. Central phosphoform identity

The NCC/SPAK features were traced to accession, peptide, modified residues,
localization evidence, and reporter completeness.

- The NCC T53-indexed row comes from a peptide carrying both T53 and Y65
  modifications.
- The NCC features indexed at positions 65 and 68 map to Y65 and Y68; they are
  not the canonical NCC N-terminal regulatory sites.
- The SPAK S383-indexed row comes from an S382/S383 co-modified peptide.
- No isolated NCC T53/T58/S71 or SPAK T243/S383 phosphoform satisfies the
  frozen qualification rule.

Consequently, the MS assay does not establish isolated canonical-site
dephosphorylation, NCC inactivation, SPAK activity loss, phosphosite occupancy,
or downstream electrolyte physiology. Prior antibody work on RR-10 is relevant
literature context, but it uses the same biological experiment and is not an
independent replication of this reanalysis.

### 3. Flight-blind distal-nephron references

Subtype annotations were rebuilt without consulting OSD-462 flight effects or
phosphoproteomic observability.

- GSE228367 raw integer counts were aggregated by the three biological
  replicates and analyzed by paired pseudobulk DCT2-versus-DCT1 modeling.
- GSE150338 supplied an independent fine-subtype direction check and
  microdissected DCT-to-CNT retention check.
- The integrated Mouse Kidney Atlas supplied whole-kidney specificity,
  multi-compartment breadth, and prespecified kidney-compartment comparators.
- Frozen sets contain 27 DCT2/CNT-transition, 17 independently derived
  DCT2/CNT-validation, 29 ASDN, 12 strict-DCT2, and 6 DCT1-core genes.
- The minimum phosphoproteomic coverage of 8 genes was frozen before the new
  flight-label-permutation results were interpreted and was not relaxed when
  DCT2/CNT yielded only 5 observable genes.

### 4. Correct inferential unit and null

The primary analysis no longer selects phosphosites by nominal \(p<0.05\) and
does not treat correlated phosphosite rows as independent biological
observations.

- The label-blind primary universe contains 8,021 qualified phosphosite
  features mapping to 3,524 fixed-null-eligible parent genes.
- A site effect is the equal-weight mean of within-plex flight-minus-ground
  log2 reporter-signal differences.
- The primary parent-gene score is
  \(\operatorname{median}(-\widehat{\beta}_{site})\); positive values indicate
  lower relative phosphopeptide signal in the flight-labeled channels.
- Every gene is standardized against its own balanced-label null.
- A set statistic is the difference between mean gene-specific \(Z\) in the
  set and the eligible non-set background.
- The complete pipeline is recalculated for all
  \(\binom{10}{5}^2=63{,}504\) balanced within-plex label assignments.
- DCT2/CNT and ASDN form one maxT family. Major kidney compartments form a
  prespecified Benjamini--Hochberg family.
- The primary canonical exclusion removes `Slc12a3`, `Stk39`, `Oxsr1`, and
  `Wnk4`; the strict exclusion additionally removes `Wnk1`.

The label permutation is a conditional exchangeability calibration. It does
not remove or identify a systematic tag-position effect because the observed
biological condition is aliased with tag position.

## Exact primary outcome

After the primary canonical-axis exclusion:

| Program | Observable genes | Competitive statistic | Conditional exact \(p\) | Adjustment | Interpretation |
|---|---:|---:|---:|---:|---|
| DCT2/CNT transition | 5 | -- | -- | -- | Non-evaluable; frozen minimum is 8 |
| ASDN | 16 | 0.718 | 0.0291 | maxT 0.0291 | Positive conditional shift, but fails selectivity |
| Podocyte | 44 | 1.112 | 0.000236 | BH \(q=0.00213\) | Largest prespecified compartment statistic |
| Fibroblast | 34 | 0.707 | 0.00375 | BH \(q=0.0169\) | Positive comparator enrichment |
| Thick ascending limb | 16 | 0.588 | 0.0462 | BH \(q=0.134\) | Not significant after comparator-family correction |

Because DCT2/CNT did not meet its coverage threshold, ASDN was the only
evaluable member of the two-set primary family. Its maxT value therefore equals
its marginal exact permutation value; it should not be presented as adjustment
across two successfully evaluated hypotheses.

The comparator analysis does not prove a podocyte or fibroblast mechanism.
The larger podocyte statistic makes ASDN fail the prespecified operational
selectivity criterion; it is not a direct statistical ASDN-versus-podocyte
contrast. All annotations attach to parent proteins measured in whole-kidney
tissue and do not determine cell of origin.

## Robustness interpretation

The robustness ladder separates directional persistence from inferential
stability.

- ASDN remains directionally positive when individual observable genes are
  removed, so the sign is not carried by a single gene.
- Alternative aggregation is material: the signed-mean parent-gene statistic
  is 0.729 (\(p=0.0323\)), whereas the one-sided maxmean statistic is 0.613
  (\(p=0.0615\)).
- Adding `Wnk1` to the canonical exclusion attenuates the ASDN statistic to
  0.632 (\(p=0.0575\)).
- Using uncentered official signal attenuates the ASDN statistic to 0.202
  (\(p=0.216\)).
- The secondary observability-matched annotation test is not significant
  (\(p=0.247\)); it does not distinguish ASDN from genes with similar site
  count, intensity, and missingness.
- Parent-protein subtraction yields a larger standardized ASDN statistic in a
  reduced universe, but loses 930 parent genes and leaves the podocyte
  statistic slightly larger. It does not establish phosphosite occupancy or
  ASDN selectivity.
- Official-scaled and summed-signal estimates are algebraically related in the
  source workbook. Agreement between them is a provenance/normalization
  sensitivity, not independent validation.
- The deduplicated multi-modified/composite profile remains positive after the
  primary exclusion (statistic 0.649, \(p=0.0478\)). The true
  multi-compartment broad-expression exclusion remains positive but is no
  longer conventionally significant (12 genes; statistic 0.538,
  \(p=0.0696\)).
- Multi-modified/composite, localization, de-duplication, and broad-expression
  profiles test measurement and annotation choices. None can repair failed
  compartment selectivity, tag-block aliasing, or missing isolated canonical
  phosphoforms.

## Claims permitted and prohibited

### Permitted

> Within the channel-centered analysis, parent proteins in the predefined ASDN
> set had higher standardized phosphosite-suppression scores than the remaining
> observable proteins after canonical-axis exclusion. The shift was not
> ASDN-selective, DCT2/CNT coverage was insufficient for the frozen test, and
> the contrast is inseparable from reporter-tag block.

> Only five genes from the frozen DCT2/CNT-transition signature were
> phosphoproteomically observable. The DCT2/CNT hypothesis was therefore
> non-evaluable; this neither supports the proposed extension nor establishes
> absence of a DCT2/CNT response.

### Prohibited

- “Spaceflight suppresses NCC/SPAK regulatory phosphorylation.”
- “Suppression extends beyond NCC/DCT1.”
- “DCT2/CNT phosphoregulation was suppressed.”
- “ASDN pathway activity decreased.”
- “The result is robust across normalization.”
- “Parent-protein adjustment proves phosphorylation-specific regulation.”
- “Podocyte phosphorylation was suppressed.”
- “DCT2/CNT showed no effect.”
- Any individual nominal ASDN phosphosite presented as a discovery without
  residue-function evidence and appropriate multiplicity control.

## Disposition of the twelve-revision history

There is no stronger hidden DCT2/CNT result in an earlier draft. V11 introduced
the subtype question; v12 promoted the same thresholded result into a stronger
headline while its own score and adjusted-model sensitivities weakened the
claim.

The strongest genuinely distinct positive result elsewhere in the repository
is the 48-gene Grey60 ECM/cell-migration RNA module: it shows an
ISS-terminal-arm-associated eigengene shift in young and old RRRM-2 mice and is
strongly preserved in the ground-control reference network. It should be
developed, if desired, as a separate RNA paper with collection/recovery context
and limited external recurrence stated plainly. It should not be appended to
the late-distal report.

Retire from the main late-distal manuscript: PPAR sex differences,
LIONESS/node2vec, age rewiring, live-return reversal, TLR4 and S1P mechanism
branches, KSEA as “independent confirmation,” human urine, CMap, IRI spot-level
tests, dDAVP, low-potassium analogy, aging, and exhaustive composition/network
dashboards. These analyses do not repair the assay or subtype-identifiability
problems.

## Recommended next experiment

A biology-forward follow-up requires all three:

1. counterbalanced or cross-plex label-swapped isobaric-tag allocation;
2. targeted measurements that isolate the canonical phosphoforms of interest;
3. DCT1/DCT2/CNT-resolved sampling, paired with transporter function,
   electrolytes, renin/aldosterone, and blood-pressure physiology.

Without those data, further narrative expansion or additional external
reference mining cannot make the DCT2/CNT claim identifiable.

## Authoritative artifacts

- Stage 0 design/provenance:
  `data/results/run_20260728_osd462_stage0/reporting/`
- Frozen subtype reference:
  `data/processed/v13_dct_reference/`
- Definitive exact inference:
  `data/results/run_20260729_v13_continuous_phospho_exact_final/`
- Read-only figures and claim report:
  `data/results/run_20260729_v13_continuous_phospho_exact_final/reporting/`
