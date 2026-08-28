# Strongest current-data paper: cross-mission podocyte-associated kidney RNA

Date: 2026-08-11  
Decision: **proceed as the strongest current-data biological manuscript**  
Scope: five terminal mouse spaceflight missions; bulk kidney RNA  
Not in scope: a methods paper, clinical biomarker prediction, AKI/CKD diagnosis,
or a mechanistic claim about podocyte injury

## Recommended title

> **A recurrent podocyte-associated kidney transcript program across five mouse spaceflight missions**

Acceptable more conservative alternative:

> **Cross-mission mouse kidney transcriptomics identifies a podocyte-associated expression shift after terminal spaceflight**

## One-sentence paper

Across five terminal mouse missions, flight kidneys showed higher abundance of
a high-specificity podocyte-associated transcript program, while prespecified
tubular-injury, fibrosis/maladaptive-repair, and distal-transport programs did
not recur; bulk tissue cannot determine whether the podocyte-associated shift
reflects cellular representation, kidney sampling, or coordinated
within-podocyte transcription.

## Why this is the strongest surviving paper

This result is stronger than the repository's other candidate stories because
it combines:

1. five distinct terminal-flight studies rather than one assay or one mission;
2. one biological animal as the analysis unit and one estimate per mission;
3. blocked whole-pipeline label permutations with family-wise correction;
4. comparison against a complete frozen kidney-compartment family rather than
   promotion of a hand-picked winning cell type;
5. a moderate pooled effect with a small-mission confidence interval excluding
   zero and no estimated between-mission heterogeneity;
6. recognizable podocyte biology in the gene-level leading edge; and
7. explicit negative results for the tubular-injury, fibrosis, distal-transport,
   Grey60, NCC-novelty, and DCT2-extension alternatives.

The manuscript's novelty is biological: a compartment-associated renal response
that recurs across mouse missions and differs from the tubule-centered story
that motivated the repository. The statistical machinery is support, not the
subject of the paper.

## Prior-art position

Do not claim that podocyte-associated molecules have never been mentioned in
spaceflight kidney work. Direct inspection of the 2024 Cosmic Kidney Fig. 5
shows that its PODXL/NID1 consensus includes nominal direction calls from mouse
kidney RNA as well as kidney proteomics; the accompanying prose separately
notes postflight human urine proteins. That study also reported glomerular
disease-ontology enrichment, while centering its direct spaceflight
interpretation on transport dephosphorylation and nephron/DCT remodeling. A
later RR-1/RR-3 transcriptomic comparison emphasized strain, lipid metabolism,
ECM, and TGF-beta signaling.

A forward-citation audit through 2026-08-11 found no published recurrent
podocyte-transcript program across independent mouse spaceflight missions. The
defensible novelty is therefore the mission-level, compartment-referenced
result: a distributed podocyte-associated RNA program is higher on average
across five terminal mouse missions even though generic fibrosis and
distal-transport transcript programs do not recur. This should be written as a
program-level extension and a different estimand, not as a correction of Siew
or a claim that prior work never examined glomerular biology.

## Exact evidence

### Primary clinically anchored family

All four scores were oriented so positive values represented the prespecified
adverse direction. The negative barrier-identity estimate therefore means
higher podocyte/barrier-marker abundance in flight.

| Program | Pooled Hedges g (mKH 95% CI) | maxT FWER | I2 | Result |
|---|---:|---:|---:|---|
| Glomerular barrier identity loss | -0.716 (-1.361, -0.072) | 0.00180 | 0.0% | Opposite-direction podocyte-associated lead |
| HAVCR1/LCN2 tubular injury | 0.570 (-0.095, 1.234) | 0.0263 | 8.0% | Suggestive but fails the small-mission interval criterion |
| Fibrosis/maladaptive repair | 0.311 (-0.732, 1.355) | 0.785 | 62.5% | No recurrent response |
| Distal transport identity loss | -0.153 (-1.319, 1.012) | 0.986 | 68.8% | No recurrent response |

### Complete kidney-compartment family

Among 49 evaluable compartment-by-specificity sets, only the
high-specificity podocyte program passed the 100,000-permutation max-|T|
family:

- pooled g = 0.689;
- mKH 95% CI = 0.042 to 1.336;
- maxT FWER = 0.0189;
- I2 = 0%;
- four of five mission estimates were positive.

Mission effects were -0.117 (OSD-102), 1.166 (OSD-163), 0.676 (OSD-253),
0.541 (OSD-462), and 1.174 (OSD-771). The prediction interval crossed zero,
so recurrence in a future mission is not guaranteed.

### Distributed biological coherence

The high-specificity atlas set contained 157 genes; 142 were observable in all
five missions. Of those 142:

- 105 had a positive pooled estimate;
- 69 were positive in at least four missions; and
- 23 were positive in all five missions.

Recognizable podocyte/glomerular machinery near the positive leading edge
included Nphs2, Col4a3, Myo1e, Wt1, Nphs1, Podxl, Crb2, Lmx1b, Mafb, Synpo,
Ptpro, and Plce1. Individual gene p-values are descriptive and must not be
presented as separately corrected discoveries.

### Essential claim-separation result

The six canonical barrier genes were strongly correlated with a disjoint
high-specificity podocyte score within missions (r = 0.77-0.90). The pooled
flight coefficient attenuated from 0.507 (0.131, 0.883) to 0.070 (-0.440,
0.581) after conditioning on that disjoint score. The result is therefore a
broader podocyte-associated program, not selective regulation of six barrier
genes.

The complete family-wise result remained after removing all canonical barrier
genes. With Nphs1, Nphs2, Synpo, Ptpro, Magi2, and Wt1 removed, pooled g was
0.675 (mKH 95% CI 0.027 to 1.322; max-|T| FWER = 0.0233; I2 = 0%). Removing
Podxl as well gave g = 0.676 (0.029 to 1.323; max-|T| FWER = 0.0229). Cd2ap
was requested for the expanded exclusion but was not a member of the source
high-specificity set. This gate therefore passes: the result is not generated
by the canonical barrier core.

Replacing the C57BL/6J OSD-253 arm with its C3H/HeJ arm also retained the
core-disjoint result (g = 0.893, 95% CI 0.005 to 1.782; family-wise p = 0.033),
although heterogeneity rose to 42.7%. This is a within-mission, post-hoc
cross-strain sensitivity, not a sixth independent mission.

## Manuscript architecture

### Introduction

Use three short moves:

1. Spaceflight kidney work has emphasized tubular transport, radiation, and
   broad remodeling, but whole-kidney studies span several missions and renal
   compartments.
2. Clinically grounded renal programs provide a finite external biological
   frame: glomerular barrier identity, tubular epithelial injury,
   fibrosis/maladaptive repair, and distal transport.
3. Ask which tissue programs recur across terminal mouse missions, then test
   the lead against a complete kidney-compartment reference.

Do not introduce NCC phosphorylation, DCT2 scoring, Grey60, KSEA, CMap, human
urine, or the twelve-draft project history.

### Results 1: Expected kidney-injury programs do not recur uniformly

Report all four frozen axes together. This is not a generic negative-results
section; it establishes that the signal is not simply "every kidney injury
program rises after flight." State that tubular injury is suggestive but does
not pass the complete evidence rule, while fibrosis and distal transport are
heterogeneous and null.

### Results 2: A high-specificity podocyte-associated program recurs

Report the 49-set family, the five mission effects, the pooled estimate,
family-wise p-value, heterogeneity, and prediction interval. Show all renal
compartments so the reader can see that podocytes were not promoted after
hiding a stronger endothelial, fibroblast, immune, proximal-tubule, TAL, DCT,
or collecting-duct result.

### Results 3: The result is distributed across podocyte machinery

Show gene-direction coherence and a compact leading-edge heat map. Emphasize
the set-level result. Include leave-one-mission, leave-one-gene, median-score,
RNA-preparation, OSD-163 mapping-rate, and matched-random-panel sensitivities.
Make removal of the six- and eight-gene canonical barrier cores the central
adversarial test. Add the PODXL/NID1 forced-gene result and the strict
same-atlas-tier matching audit. State that neither named gene is an individual
cross-mission discovery and that the common-support target-minus-matched
interval crosses zero despite its blocked-label support.

### Results 4: Bulk RNA leaves cellular representation unresolved

Show the canonical barrier score versus the disjoint podocyte proxy and the
attenuated adjusted estimate. Do not add an inferential comparison with exact
neighboring glomerular compartments using the current atlas. The exact-label
audit found podocytes in several source studies, but glomerular endothelial,
mesangial, and parietal epithelial labels each occurred in only one source
study; the high-specificity glomerular-endothelial set also contained fewer
than eight genes. A negative result for those controls would therefore use
inadequately replicated marker definitions. End with the three unresolved
explanations:

- a larger podocyte RNA fraction or altered podocyte representation;
- cortical/glomerular sampling differences; or
- coordinated cell-intrinsic transcription within podocytes.

Do not select one explanation without direct evidence.

### Discussion

Lead with the observation, not its caveats: a recurrent, podocyte-associated
bulk-kidney RNA shift was the only compartment result to pass the frozen family.
Then explain why this redirects attention beyond the tubule-centered literature.
State plainly that higher marker abundance is not synonymous with protection,
activation, injury, or increased podocyte number. Finish with one decisive
validation experiment: quantitative podocyte/glomerular histology plus urinary
albumin-to-creatinine or protein-to-creatinine in a future flight cohort.

## Main figures

1. **Five missions and four renal programs.** Cohort schematic plus forest plot
   of the four frozen axes.
2. **Kidney-compartment specificity.** Full compartment-family forest plot plus
   mission-level podocyte estimates.
3. **Distributed podocyte-program coherence.** Leading-edge heat map, gene sign
   counts, and canonical-core exclusion.
4. **Robustness and biological ambiguity.** Leave-one-out/technical sensitivity
   dashboard, disjoint-proxy adjustment, and a minimal interpretation schematic.

## Keep out of the main paper

- OSD-462 phosphoproteomics and NCC/SPAK;
- DCT2/CNT and ASDN phosphosite enrichment;
- Grey60/WGCNA/network topology;
- CMap, IRI spatial projection, dDAVP, aging, and low-potassium analyses;
- OSD-656 human urine proteins;
- live-return/recovery cohorts as if they were an independent time course; and
- the invalid gene-wise Stouffer fibrosis p = 0.0007.

OSD-253 C3H/HeJ and recovery-arm projections may appear as explicitly post-hoc
supplementary context. They must not be counted as additional independent
missions.

## Claim boundaries

### Supported

> Flight kidneys had higher bulk-tissue abundance of a high-specificity
> podocyte-associated transcript program across five terminal mouse missions.

> The average shift was compartment-selective within the executed frozen
> kidney-atlas family and distributed across many podocyte-associated genes.

### Not supported

- podocyte injury, protection, activation, hypertrophy, or increased number;
- albuminuria, filtration-barrier dysfunction, altered GFR, AKI, or CKD;
- prediction of a urinary biomarker;
- a selective change in canonical barrier genes independent of the broader
  podocyte-associated program; or
- a mechanism involving microgravity, radiation, hemodynamics, or hormonal
  regulation.

## Submission decision

This is the strongest manuscript available from the current repository and is
worth writing as a concise biological cross-mission paper. The decisive
canonical-core exclusion gate passed. The neighboring-glomerular-cell control
could not be tested credibly with the current reference and therefore remains
a limitation rather than a submission gate:

1. **Passed:** the high-specificity podocyte result retained family-wise support
   after removal of canonical barrier genes.
2. **Passed:** forcing `Podxl` and `Nid1` into the program left the estimate
   essentially unchanged; neither gene was individually significant.
3. **Passed for the observed missions:** strict same-atlas-tier matching and
   q95 common-support trimming retained blocked-label support; the
   target-minus-matched mKH interval crossed zero, limiting future-mission
   generalization.
4. **Not evaluable with the current atlas:** exact glomerular endothelial,
   mesangial, and parietal epithelial identities lack cross-study marker
   replication, so cortical/glomerular sampling remains unresolved.

The title and text must therefore retain "podocyte-associated" and must not
claim a cell-intrinsic or podocyte-specific change.

## Reproducible source outputs

- `data/results/run_20260811_clinical_renal_axes_cross_mission/primary_meta_results.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/compartment_context_meta_results.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/podocyte_gene_coherence_summary.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/podocyte_gene_coherence_meta.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/barrier_proxy_adjustment_meta.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/disjoint_podocyte_family/compartment_context_meta_results.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/disjoint_podocyte_family_c3h/compartment_context_meta_results.tsv`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/podxl_nid1_forced_sensitivity/`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/strict_podocyte_matching/`
- `docs/PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md`
- `docs/PODOCYTE_STRICT_MATCHING_AUDIT_2026-08-11.md`
- `figures/clinical_renal_axes_cross_mission/`
