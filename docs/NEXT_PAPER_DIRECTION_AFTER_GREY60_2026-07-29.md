# Next paper direction after Grey60

**Status:** recommendation; not yet an executed analysis  
**Priority:** OSD-462 matched-library robustness first  
**Fallback:** cross-mission control-choice sensitivity second

Both the DCT2/CNT-ASDN extension and the Grey60 recurrent-module proposition
have now failed their locked claim gates. The next project should exploit a
design feature that makes a question identifiable rather than select another
nominal pathway.

## Priority 1 — OSD-462 matched-library robustness

### Question

> When the same kidney animals are profiled by UPX 3-prime-tag, polyA mRNA,
> and total-RNA sequencing, how much does library preparation change the
> inferred flight effect?

This is a measurement and inference paper, not a claim that one preparation
reveals the true biological answer.

### Why it is the strongest current direction

- The same 10 flight and 10 ground-control animals were measured three times.
- Biological differences are therefore held fixed across preparations.
- Grey60 already changes sign across preparations.
- Genome-wide flight-effect correlations were low in the preliminary audit.
- The repository contains raw counts, sample maps, VST matrices, official
  differential-expression tables, and assay metadata for all three
  preparations.

### Frozen inputs

- `data/external/osdr/OSD-462/GLDS-462_rna_seq_RSEM_Unnormalized_Counts_UPX_GLbulkRNAseq.csv`
- `data/external/osdr/OSD-462/GLDS-462_rna_seq_RSEM_Unnormalized_Counts_mRNA_GLbulkRNAseq.csv`
- `data/external/osdr/OSD-462/GLDS-462_rna_seq_RSEM_Unnormalized_Counts_totRNA_GLbulkRNAseq.csv`
- the corresponding three `SampleTable`, `VST_Counts`,
  `differential_expression`, and metadata assay files.

BL and VIV animals should be QC/sensitivity contexts, not part of the primary
FLT-versus-GC interaction.

### Primary estimand

Construct a common protein-coding gene universe before calculating contrasts.
After within-preparation normalization, fit a repeated-measures model for each
gene:

```text
expression ~ flight + preparation + flight:preparation + animal
```

The primary object is the preparation-dependent deviation in the flight
effect, not the number of FDR-significant genes.

The global endpoint should summarize the preparation-by-flight interaction
across the common transcriptome and calibrate it by permuting FLT/GC labels at
the **animal** level while rerunning normalization and modeling.

Report:

- pairwise Spearman and cosine concordance of flight-effect vectors;
- bootstrap confidence intervals;
- interaction variance relative to sampling variance;
- same-animal score correlations;
- official-OSDR-versus-harmonized-pipeline concordance.

### Prespecified secondary programs

Use a small frozen family:

- Grey60/ECM remodeling;
- DCT/NCC transport;
- endothelial;
- stromal/fibroblast;
- immune/injury;
- oxidative stress.

Test continuous program effects and preparation interactions. Correct only
across this locked family. DEG counts and overlap diagrams are descriptive.

### Required sensitivities

- protein-coding genes only;
- several common CPM thresholds;
- common dynamic-range matching;
- raw-count limma-voom versus VST/rank-based effects;
- leave-one-animal-out;
- robust regression;
- removal of low-count and preparation-specific biotypes;
- sample fingerprint and aliquot-identity verification.

### Hard go/no gates

All are required:

1. Sample and aliquot identity passes, with a complete 20-animal-by-3-preparation
   block.
2. The global preparation-by-flight endpoint has animal-label-permutation
   p<0.01.
3. At least one frozen biological program has an interaction interval
   excluding zero and genuinely opposite-signed preparation effects.
4. The result survives every leave-one-animal analysis and common-gene
   filtering.
5. Discordance is not restricted to low-expression, noncoding, or one
   normalization-specific gene class.

Retire the paper if harmonization raises pairwise effect concordance to at
least 0.70 and removes coherent interactions, or if the result is caused by
sample mismatch, one animal, low-expression features, or only the official
pipeline.

### Existing code to reuse

- `scripts/osd462/00_harmonize.py`
- `scripts/osd462/_common.py`
- `scripts/osd462/04_rna_recurrence.py`
- `scripts/osd462/08_stage0_provenance_audit.py`
- `src/multiomics/osd462_stage0.py`
- `src/v11/recurrence_meta.py`

A new repeated-preparation model and whole-profile animal-permutation engine
will be required.

### Paper only if it passes

A defensible working title would be:

> RNA library preparation alters inferred spaceflight kidney responses in
> matched biological samples

The paper would contain four results:

1. matched-sample and assay identity;
2. genome-wide preparation-by-flight heterogeneity;
3. biological-program consequences;
4. robustness and official-pipeline comparison.

## Priority 2 — cross-mission control-choice sensitivity

Run this only after resolving OSD-462 preparation instability.

### Question

> How much does the choice of hardware, vivarium, basal, or rerun control alter
> the estimated kidney response attributed to flight?

### Studies

- OSD-771: ISS-terminal and live-return arms, age retained;
- OSD-253: strain, duration, original blue-light GC, and rerun white-light GC
  kept separate;
- OSD-462: one biological study with nested preparation contexts;
- OSD-513: GC and VIV;
- OSD-163: GC and basal;
- OSD-102 excluded from the primary endpoint because the local design has no
  alternative reference group.

### Primary estimand

Within each study, estimate the change in a continuous flight effect produced
by changing the reference group:

```text
reference shift = (flight - control A) - (flight - control B)
```

Retain the covariance induced by reusing the same flight animals. Never pool
expression matrices or count OSD-462 preparations and OSD-253 cells as
independent studies.

Use the same six frozen programs as Priority 1. Report continuous effects,
reference shifts, confidence intervals, contrast-vector concordance, and
direction changes. FDR-status flips are secondary.

### Hard go/no gates

All are required:

1. At least three independent missions contain valid alternative controls.
2. A global control-choice test has empirical p<0.01.
3. Reference choice explains at least 25% of non-sampling pathway-effect
   variance.
4. At least two frozen pathways show a nonzero reference shift in at least two
   independent missions, or a replicated sign reversal with a direct
   interaction.
5. Results survive leave-one-study-out and exclusion of OSD-253 rerun/batch
   contrasts and nonsynchronous basal comparisons.

Retire it if median flight-contrast concordance is at least 0.80, observed
reference shifts match bootstrap noise, or OSD-253 alone drives the result.

### Existing code to reuse

- `scripts/run_contrast_vector_framework.py`
- `config/contrast_vector_framework.yaml`
- `src/networks/contrast_vectors.py`
- `src/validation/osd_external_validation.py`
- `src/v11/recurrence_meta.py`

## Biological fallback — design-aware heterogeneity map

If the matched-library paper fails, the best biological fallback is a
cross-mission kidney response **heterogeneity and boundary** analysis across
OSD-102, OSD-163, OSD-253, OSD-462, OSD-513, and OSD-771.

It must ask which programs vary with terminal versus live-return collection,
duration, strain, sex, preparation, and control design. It cannot be framed as
a universal signature. Its likely contribution is a map of what does *not*
generalize, which is scientifically useful but less novel than the
matched-library design.

## Directions to retire

- another post-hoc Grey60 subset or replacement module;
- PPAR sex differences;
- LIONESS/node2vec classification;
- TLR4 causality from OSD-253 strain contrasts;
- LAR “reversal”;
- DCT2/CNT or ASDN claims from the current OSD-462 phosphoproteome;
- a recurrent terminal Grey60 biology paper;
- RNA-to-protein propagation until reporter-tag aliasing and missing external
  proteomes are resolved.

## If a strong biology paper is the requirement

The current public-data repository is no longer enough. The clean experiment
would add segment-resolved kidney evidence: validated NCC/SPAK phosphosite
mapping, DCT1/DCT2/CNT-resolved protein or phosphoprotein measurements, and
renal electrolyte/renin-aldosterone physiology. Without new data, the honest
publishable space is measurement robustness, design sensitivity, and
cross-mission heterogeneity.

