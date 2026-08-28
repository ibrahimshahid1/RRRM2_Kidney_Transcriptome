# Cross-mission clinical renal-axis analysis: final decision report

Date: 2026-08-11  
Status: completed targeted retrospective analysis  
Primary lock: `config/clinical_renal_axes_cross_mission.yaml`

## Executive decision

The proposed paper—“which clinically anchored renal injury axes recur after
spaceflight?”—is a **no-go as an injury-biomarker paper**. Neither the
HAVCR1/LCN2 tissue-injury panel, fibrosis/maladaptive-repair panel, nor distal
transport panel passed the complete cross-mission evidence standard.

The analysis did reveal one genuinely distinct lead: **bulk-kidney abundance
of highly podocyte-specific transcripts was higher after terminal spaceflight
in four of five mouse missions**. This result survived the frozen
kidney-compartment family, observability-matched random panels, alternative RNA
preparations, mission removal, gene removal, and the flagged OSD-163 mapping-rate
sensitivity.

That lead is not evidence of podocyte injury, improved glomerular filtration,
albuminuria, or barrier failure. Adjustment for a disjoint podocyte-marker
program absorbed the six-gene barrier-core effect. The observable result is
therefore a **broader podocyte-associated tissue-expression shift**, which bulk
RNA cannot separate into altered cell representation, sampling composition,
or coordinated within-podocyte transcription.

The practical decision is:

- Do not convert the current results into a KIM-1/NGAL, fibrosis, DCT, or
  glomerular-disease biomarker manuscript.
- Proceed with a concise cross-mission biological paper centered on the
  podocyte-associated bulk-kidney transcript program.
- Treat direct glomerular histology, cell-resolved expression, ultrastructure,
  and urinary albumin/protein normalized to creatinine as the decisive next
  validation, not as measurements that the present repository contains.

## Question and analysis unit

The frozen question was:

> Which clinically anchored renal tissue programs show directional engagement
> during terminal spaceflight, and which recur across independent mouse
> missions?

The primary terminal-flight synthesis used one mission-level estimate from
each of OSD-102, OSD-163, OSD-253, OSD-462, and OSD-771. Duration strata in
OSD-253 and age strata in OSD-771 were combined within mission and were not
counted as independent studies. OSD-513 and the OSD-771 LAR arm were retained
only as recovery moderators.

Animal-level signed gene scores were calculated before group comparison.
Mission-level Hedges g values were synthesized with REML and modified
Hartung–Knapp uncertainty. The primary family used 100,000 label permutations
within each mission's age/duration exchangeability blocks and maximum-|T|
family-wise error control over the four frozen axes.

The lock is retrospective, not a preregistration. It was created after earlier
repository work had shown ECM-like and distal-transport patterns. The
glomerular and HAVCR1/LCN2 panel results had not been inspected before the lock.
During a label-blind observability audit, the CPM eligibility threshold was
amended from 1.0 to 0.1 before any new group effect was calculated; the reason
is recorded in the configuration.

## Primary results

All scores were oriented so positive values represented the prespecified
adverse direction. Thus, the negative glomerular estimate means **higher**, not
lower, barrier/podocyte-marker expression in flight.

| Frozen program | Pooled Hedges g (mHK 95% CI) | maxT FWER | I² | Decision |
|---|---:|---:|---:|---|
| Glomerular barrier identity loss | -0.716 (-1.361, -0.072) | 0.00180 | 0.0% | Reject “identity loss”; investigate the opposite-direction podocyte-associated shift |
| Tubular epithelial injury (Havcr1/Lcn2) | 0.570 (-0.095, 1.234) | 0.0263 | 8.0% | No-go: permutation signal, but the small-mission mHK interval crosses zero |
| Fibrosis/maladaptive remodelling | 0.311 (-0.732, 1.355) | 0.785 | 62.5% | No-go: null and heterogeneous |
| Distal transport identity loss | -0.153 (-1.319, 1.012) | 0.986 | 68.8% | No-go: null, heterogeneous, and directionally inconsistent |

The prediction interval also crossed zero for the glomerular panel (-1.455,
0.023). This matters: the average effect was coherent in the observed mission
set, but the analysis does not guarantee the same direction in a new mission.

### Mission-level pattern

- Glomerular adverse-orientation effects were approximately 0.001, -0.981,
  -1.288, -0.913, and -0.413 in OSD-102, OSD-163, OSD-253, OSD-462, and
  OSD-771, respectively. Four missions therefore showed higher podocyte/barrier
  marker abundance; OSD-102 was essentially null.
- Tubular injury was positive in four missions but negative in OSD-163. Neither
  Havcr1 nor Lcn2 had a pooled gene-level interval excluding zero, and Lcn2
  contributed 60.9% of the absolute two-gene signal.
- Fibrosis ranged from -0.715 in OSD-462 to +1.500 in OSD-771. The repair and
  ECM subdomains were each null.
- Distal transport ranged from -1.498 in OSD-462 to +0.719 in OSD-102. Its DCT
  and ASDN subdomains were each null.

These results retire the older gene-wise Stouffer matrix/ECM p=0.0007 as
primary evidence. That analysis treated correlated genes as inferential units;
the animal- and mission-level analysis does not reproduce a general recurrent
fibrosis program.

## Adversarial audit of the podocyte lead

### Evidence that the signal is real at the annotation level

1. The signed-median result was similar: g=-0.739, mHK CI -1.386 to -0.092,
   maxT FWER=0.00130.
2. The common-observable-gene analysis was identical to the primary result.
3. Adding Podxl and Cd2ap retained the direction and interval exclusion:
   g=-0.704, CI -1.345 to -0.063.
4. Every leave-one-mission estimate retained the direction (-0.579 to -0.855).
5. Every leave-one-core-gene estimate retained the direction and at least 73%
   of the full magnitude. Nphs2 supplied 34.2% of the absolute signed gene
   effect, but removing Nphs2 still gave g=-0.524.
6. Replacing the OSD-462 total-RNA preparation with mRNA or UPX retained the
   result: g=-0.673 and -0.686, with both mHK intervals excluding zero.
7. Using the OSD-253 white-light rerun control retained the direction
   (g=-0.785), although its interval crossed zero (-1.591, 0.021).
8. Residualizing the flagged OSD-163 uniquely mapped-read percentage retained
   the meta-effect: g=-0.684, CI -1.327 to -0.042.
9. Against 10,000 six-gene panels matched on RNA observability, expression,
   variability, and kidney-cell-type breadth, the barrier-core expression
   increase was unusual (two-sided matched-panel p=0.0062). Results were stable
   to matching pools of 20, 50, 100, and 200 candidates per target.
10. In a separately frozen 49-set compartment family, only the
    high-specificity podocyte program passed: g=0.689, mHK CI 0.042 to 1.336,
    maxT FWER=0.0189, I²=0%. Four of five missions were positive.
11. The result survived removal of Nphs1, Nphs2, Synpo, Ptpro, Magi2, and Wt1
    while rerunning a 51-set, 100,000-permutation family: g=0.675, CI 0.027 to
    1.322, maxT FWER=0.0233, I²=0%. Removing Podxl as well produced essentially
    the same result (g=0.676, CI 0.029 to 1.323, FWER=0.0229).
12. Replacing the C57BL/6J OSD-253 arm with its C3H/HeJ arm retained the
    core-disjoint result (g=0.893, CI 0.005 to 1.782, FWER=0.0330), although
    heterogeneity increased to 42.7%. This is a within-mission sensitivity,
    not an additional independent mission.

### Evidence that limits the biological claim

The six core genes were highly correlated with a disjoint high-specificity
podocyte score within missions (Pearson r=0.77–0.90). In a blocked regression:

| Model | Flight coefficient (mHK 95% CI) | I² |
|---|---:|---:|
| Unadjusted, HC3 | 0.507 (0.131, 0.883) | 0.0% |
| Adjusted for disjoint podocyte proxy, HC3 | 0.070 (-0.440, 0.581) | 76.0% |

This does not prove a compositional artifact: the proxy could represent
podocyte abundance, coordinated podocyte state, or both. It does show that the
six barrier genes do not change specifically beyond the broader podocyte
program. Accordingly:

- allowed: “higher bulk-kidney abundance of a podocyte-associated transcript
  program”;
- not allowed: “podocyte activation,” “podocyte protection,” “barrier
  improvement,” “glomerular injury,” “albuminuria,” or “filtration change.”

An exact-label atlas audit could not supply a strong neighboring-glomerular
negative control. Podocytes occurred in several atlas source studies, whereas
exact glomerular-endothelial, mesangial, and parietal-epithelial labels each
occurred in only one; the high-specificity glomerular-endothelial definition
also fell below the eight-gene evaluability rule. A negative flight comparison
would therefore be weak evidence. General cortical/glomerular sampling remains
an unresolved explanation and is stated as such rather than assigned a p-value.

The result is also secondary in the larger project history. The broad
compartment scan followed the targeted barrier result, and the all-enriched,
broad-enriched, scaffold-excluded, and within-kidney-not-broad podocyte tiers
did not pass the full 49-set family. The high-specificity tier did. That makes
this a strong hypothesis, not a clean confirmatory endpoint.

The matched-panel null is supportive rather than exact. The barrier genes sit
near the extreme of the atlas abundance/specificity distribution, so even their
nearest candidate pools generally had lower atlas expression. Pool-size
stability reduces concern about an arbitrary nearest-neighbour choice, but the
matched p-value should not be described as eliminating all marker-selection
effects.

## Human urine context

OSD-656 cannot validate the mouse result. Its inflammatory urine panel measured
only 4 of the 34 frozen tissue-axis genes: LCN2, CCL2, EGF, and TIMP1. It did not
measure HAVCR1, any of the six barrier-core markers, any ECM gene, or any distal
transport gene.

Relative to each subject's averaged preflight baseline, at R+1:

- LCN2 was higher in all three subjects with an available sample (mean paired
  NPQ delta +1.124).
- CCL2 was higher in all three (mean +1.011).
- EGF and TIMP1 were mixed.

R+45 was mixed. By R+82, CCL2, EGF, and TIMP1 were lower in all four subjects;
LCN2 was lower in three and essentially flat in the fourth.

No p-values were calculated. This is, at most, transient postflight LCN2/CCL2
urine context. It has four crew members, no nonflight control, no inflight
sample, relative NPQ rather than a clinical concentration, and no urine
creatinine correction. Urinary protein and kidney-tissue RNA are not equivalent
measurements.

## Relationship to OSD-462 phosphoproteomics

The earlier podocyte-annotated phosphosite result is not a second validation:

- OSD-462 contributes animals to both the RNA and phosphoproteomic results.
- Whole-kidney parent-protein annotations do not localize a phosphosite to a
  podocyte.
- Condition is aliased with reporter-tag position in the OSD-462 MS design.
- The frozen all-enriched podocyte phosphoproteomic test did not pass its full
  family; selected high-specificity and observability-matched variants did.

It may be reported as exploratory assay-layer context only. P-values should not
be combined, and “independent” or “orthogonal confirmation” is not defensible.

## Direct-measurement availability audit

No available spaceflight cohort currently supplies a direct matched test of the
podocyte/barrier interpretation.

- RR-10/OSD-462 has the same bulk omics used here. Its published morphometry is
  NCC-positive DCT morphometry, not quantitative glomerular or podocyte
  morphometry. Qualitative cleared-kidney imaging reported no obvious gross
  lesion, but that cannot resolve podocyte foot processes or the glomerular
  basement membrane.
- OSD-457 provides bulk kidney RNA and selected systemic measurements, but no
  public albumin/protein urine endpoint, podocyte stain, permeability test, or
  quantitative glomerular morphometry was identified.
- OSD-771 used a right-kidney coronal section for RNA. The public study records
  blood collection but no urine renal chemistry, kidney histology, or imaging
  endpoint.
- OSD-656 lacks albumin, total protein, urine creatinine, cystatin C, nephrin,
  podocin, synaptopodin, and WT1.
- OSD-575 contains human serum chemistry from four Inspiration4 crew and can
  provide filtration context, but eGFR is neither a podocyte measurement nor a
  permeability endpoint.

The Cosmic Kidney study's increased urine Protein:Cr and glomerular capillary
microthrombi came from the separate NSRL22A ground-GCR experiment approximately
six months after radiation exposure. That is useful radiation-model context,
not a spaceflight replication, and it contains no podocyte-specific stain or
ultrastructure.

## What the repository establishes and does not establish

The current evidence supports:

> Across five terminal mouse missions, spaceflight was associated with higher
> bulk-kidney expression of a high-specificity podocyte-associated program; the
> canonical barrier-gene signal was not separable from that broader program,
> and neither podocyte injury nor filtration-barrier dysfunction was measured.

The repository does not establish:

- a recurrent fibrosis/maladaptive-repair response;
- recurrent distal-transport transcript loss;
- a definitive tissue HAVCR1/LCN2 injury response;
- urinary KIM-1, NGAL, albumin, or fibrosis-biomarker changes;
- altered GFR, creatinine, cystatin C, or urine output;
- podocyte number, foot-process architecture, or glomerular permeability;
- AKI, CKD, renal cancer, or a kidney morphological lesion.

## Go/no-go by manuscript concept

| Manuscript concept | Decision | Reason |
|---|---|---|
| Cross-mission clinical kidney-injury axes | **No-go** | No adverse biological axis passed the complete evidence standard |
| Recurrent renal fibrosis/remodelling | **No-go** | Null average, I²=62.5%, and opposite mission effects |
| Tissue NGAL/KIM-1 response | **No-go** | Suggestive direction only; mHK interval crosses zero; two-gene endpoint is Lcn2-heavy |
| Distal nephron transcript adaptation | **No-go** | Null and heterogeneous; does not rescue the NCC paper |
| Podocyte/barrier disease biomarker | **No-go** | No barrier function or injury measurement; core effect is absorbed by broader podocyte identity |
| Recurrent podocyte-associated tissue program | **Go, with a narrow bulk-tissue claim** | Strongest distinct cross-mission result; bulk RNA cannot resolve its cellular meaning |

## Recommended biological paper

The strongest current-data manuscript is:

> **A recurrent podocyte-associated kidney transcript program across five mouse
> spaceflight missions**

Its central claim is higher bulk-kidney abundance of a podocyte-associated
program, not podocyte injury or barrier dysfunction. The full manuscript
architecture and claim boundaries are in
`docs/PODOCYTE_CROSS_MISSION_PAPER_BLUEPRINT_2026-08-11.md`.

The next prospective biological validation should be one direct test of the
podocyte interpretation:

1. **Archived-tissue validation:** quantify WT1-positive podocyte nuclei per
   glomerulus and per tuft area; stain NPHS1/nephrin, NPHS2/podocin, and
   synaptopodin; measure glomerular tuft area and cortical sampling fraction.
2. **Functional barrier assay:** urinary albumin-to-creatinine or
   protein-to-creatinine, ideally during flight and at terminal collection.
3. **Ultrastructure:** electron-microscopy foot-process width/effacement and
   glomerular-basement-membrane measurements.
4. **Cell resolution:** a correctly labelled kidney single-nucleus or spatial
   dataset that distinguishes more podocytes from higher expression per
   podocyte. The currently corrected RRRM-1 single-cell sample labels are
   unknown and cannot be used for this purpose.

The current data are sufficient for the narrow biological paper above. A direct
assay would determine whether the program reflects podocyte number, sampling,
or cell-intrinsic state and would permit a stronger mechanistic or disease claim.

## Reproducible outputs

- Primary results: `data/results/run_20260811_clinical_renal_axes_cross_mission/primary_meta_results.tsv`
- Mission effects: `data/results/run_20260811_clinical_renal_axes_cross_mission/primary_mission_effects.tsv`
- Age/duration strata: `data/results/run_20260811_clinical_renal_axes_cross_mission/primary_stratum_effects.tsv`
- Complete compartment family: `data/results/run_20260811_clinical_renal_axes_cross_mission/compartment_context_meta_results.tsv`
- Matched-panel audit: `data/results/run_20260811_clinical_renal_axes_cross_mission/barrier_matched_panel_summary.tsv`
- Podocyte-proxy adjustment: `data/results/run_20260811_clinical_renal_axes_cross_mission/barrier_proxy_adjustment_meta.tsv`
- Canonical-core exclusion: `data/results/run_20260811_clinical_renal_axes_cross_mission/disjoint_podocyte_family/compartment_context_meta_results.tsv`
- C3H/HeJ sensitivity: `data/results/run_20260811_clinical_renal_axes_cross_mission/disjoint_podocyte_family_c3h/compartment_context_meta_results.tsv`
- Gene-level coherence: `data/results/run_20260811_clinical_renal_axes_cross_mission/podocyte_gene_coherence_summary.tsv`
- OSD-163 QC adjustment: `data/results/run_20260811_clinical_renal_axes_cross_mission/technical_qc_covariate_sensitivity.tsv`
- OSD-656 context: `data/results/run_20260811_clinical_renal_axes_cross_mission/osd656_urine_context_summary.tsv`
- Figures: `figures/clinical_renal_axes_cross_mission/`
- Provenance manifest: `data/results/run_20260811_clinical_renal_axes_cross_mission/manifest.json`

Main commands:

```bash
venv/bin/python scripts/clinical_axes/run_cross_mission.py
venv/bin/python scripts/clinical_axes/run_sensitivities.py
venv/bin/python scripts/clinical_axes/run_compartment_context.py --permutations 100000
venv/bin/python scripts/clinical_axes/run_compartment_context.py --permutations 100000 --add-podocyte-disjoint-variants
venv/bin/python scripts/clinical_axes/run_podocyte_gene_coherence.py
venv/bin/python scripts/clinical_axes/run_matched_panel_null.py
venv/bin/python scripts/clinical_axes/run_barrier_specificity_adjustment.py
venv/bin/python scripts/clinical_axes/run_osd656_urine_context.py
venv/bin/python scripts/clinical_axes/plot_results.py
```
