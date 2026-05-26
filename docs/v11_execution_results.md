# V11 Execution Results

Date: May 26, 2026

Run root: `data/results/run_20260526_v11_dct1_phospho_mediation/`

Plan executed from: `docs/v11_execution_research_plan.md`

Plan checksum recorded in run: `d8b7465e09257dd4444fc2abb969a02066f2d360e4fbc04056791a60c9755e3f`

## Executive Verdict

The v11 analysis makes the manuscript meaningfully more novel than v10, but it does not turn the paper into a causal mechanistic proof.

The strongest upgraded claim is:

> Spaceflight-suppressed OSD-462 phosphosites are enriched in a DCT1-high reference-prior subset, while the NCC/SPAK/WNK regulatory lesion remains a whole-kidney phosphoproteomic observation and the remodeling-to-transporter link is consistent with, not proof of, partial mediation.

After the composition-aware follow-up, this should be sharpened to:

> The DCT1-high top-decile phosphosite subset remains enriched for flight-suppressed sites after conservative parent-protein and bulk-composition adjustment, but the broader continuous DCT1-prior gradient and one-shot flight-by-DCT1 interaction are weak.

The safest manuscript frame is therefore:

> A DCT1-prioritized, phosphoproteome-informed, mediation-specified hypothesis of spaceflight kidney remodeling.

I would avoid "DCT1-anchored mechanism" unless the title or abstract immediately says the anchoring is reference-prior based.

## What Ran

Core scripts:

- `src/v11/build_gse228367_dct_prior.R`
- `src/v11/core_analysis.py`
- `src/v11/h2_composition_aware_phospho.py`
- `src/v11/spatial_reference_projection.py`
- `src/v11/publication_figures.py`

Pipeline entrypoint:

- `python -m src.run_all_phases --v11-only --run-id <run_id>`

Core outputs:

- `baseline/v11_baseline_lock_summary.tsv`
- `external_qc/gse228367_marker_qc.tsv`
- `external_qc/pxd001729_phosphosite_effects.tsv`
- `dct_prior/gse228367_dct1_vs_dct2_de.tsv`
- `dct_prior/osd462_phosphosite_dct1_prior.tsv`
- `h2_enrichment/h2_dct1_sensitivity_summary.tsv`
- `h2_pxd/h2_pxd001729_ddavp_antialignment_summary.tsv`
- `h2_klhl3/h2_klhl3_cul3_effects.tsv`
- `h2_composition_adjusted/h2_composition_adjusted_suppression_enrichment_single.tsv`
- `h2_composition_adjusted/h2_composition_effect_level_dct1_ladder_single.tsv`
- `h2_composition_adjusted/h2_composition_site_fixed_interaction_ladder_single.tsv`
- `h3_mediation/h3_mediation_model_summary.tsv`
- `h3_mediation/h3_mediation_power_simulation.tsv`
- `spatial_reference/spatial_reference_projection_verdict.json`
- `figures/v11/v11_h2_primary_dct1_enrichment.png`
- `figures/v11/v11_h2_composition_robustness.png`
- `figures/v11/v11_h3_mediation_forest.png`
- `figures/v11/v11_spatial_reference_projection.png`
- `figures/v11/v11_xenium_annotation_inventory.png`

## H1: Recurrent Remodeling State

H1 remains supported and should stay as the v10 backbone, not be inflated.

Key locked quantities:

| Component | Result | Interpretation |
| --- | ---: | --- |
| OSD-462 RNA recurrence | pathway cosine 0.869 | OSD-462 recurs the RRRM-2 matrix-high / DCT-low direction |
| OSD-462 protein abundance | null | RNA remodeling does not propagate cleanly as protein abundance |
| NCC regulatory phospho score | -0.832 | NCC regulatory phosphosites are suppressed while total NCC protein is flat |
| KSEA SPAK/OSR1 | -6.31 | inferred kinase-output down |
| KSEA WNK | -4.12 | inferred WNK-output down |
| Endothelial vs NCC phospho | rho = -0.762 | strongest compartment anti-correlation with NCC regulatory phosphorylation |

This is still an observational cross-omic decoupling result, but it is a strong one. It aligns with the Cosmic Kidney paper's report that spaceflight induces renal transporter dephosphorylation and nephron remodeling, but v11 adds a cleaner hypothesis structure around the bulk RNA context and phospho lesion.

## GSE228367 DCT Reference

The GSE228367 DCT reference loaded successfully.

Inventory:

- DCT1 object: 18,881 nuclei.
- DCT2 object: 6,274 nuclei.
- Replicates available: three normal-potassium controls per subtype.
- Analysis used replicate-level pseudobulk summaries where possible.

Marker sanity check:

| Gene | DCT1 mean | DCT2 mean | Read |
| --- | ---: | ---: | --- |
| Slc12a3 | 5.14 | 4.40 | present in both, higher in DCT1 |
| Pvalb | 0.279 | 0.024 | DCT1-high |
| Trpm6 | 2.40 | 1.90 | DCT1-high |
| Klhl3 | 2.84 | 2.65 | mild DCT1-high |
| Trpv5 | 0.006 | 0.711 | DCT2-high |
| Calb1 | 1.01 | 2.97 | DCT2-high |
| Stk39 | 1.98 | 1.54 | DCT1-leaning |
| Wnk4 | 1.72 | 1.40 | DCT1-leaning |

Important statistical boundary:

The strict FDR-defined DCT1 core is too sparse for binary gene-set testing: only 2 genes passed the strict rule, and 0 genes passed DCT2 core under the same rule. This is not a failure of the reference dataset; it is mostly an n=3 pseudobulk power issue plus the biological fact that DCT1 and DCT2 are related subsegments with gradients rather than cleanly separated compartments.

Manuscript consequence:

Use the continuous DCT1 enrichment score and DCT1-high percentile prior as the main reference-prior object. Do not claim that a large FDR-defined DCT1 marker set drives the result.

## H2: DCT1-Prioritized Phosphosite Suppression

Primary H2 result: partially positive.

The cleanest evidence is enrichment of suppressed phosphosites in DCT1-high percentile bins.

| Test | Suppressed set | Result | Read |
| --- | --- | ---: | --- |
| Matched-null mean DCT1 score | p < 0.05 suppressed sites | q = 0.083 | weak continuous support at FDR 0.10 |
| Fisher DCT1 top quartile | p < 0.05 suppressed sites | OR = 1.14, q = 0.016 | supported |
| Fisher DCT1 top decile | p < 0.05 suppressed sites | OR = 1.51, q = 1.13e-11 | supported |
| DCT1 top quartile, anchor genes excluded | p < 0.05 suppressed sites | OR = 1.13, q = 0.0196 | survives anchor-gene exclusion |
| DCT1 top decile, anchor genes excluded | p < 0.05 suppressed sites | OR = 1.49, q = 3.68e-11 | survives anchor-gene exclusion |
| DCT1 top quartile, NCC sites excluded | p < 0.05 suppressed sites | OR = 1.13, q = 0.0196 | survives NCC-site exclusion |
| Single-site-only DCT1 top decile | p < 0.05 suppressed sites | OR = 1.38, q = 2.84e-6 | survives, but top-quartile signal weakens |

What this supports:

Suppressed phosphosites in OSD-462 are not just the known NCC/SPAK/WNK anchor sites. They are enriched among parent genes that are DCT1-high in the external DCT reference.

What this does not support:

It does not prove the spaceflight phosphosites came from DCT1 cells. OSD-462 is whole-kidney phosphoproteomics. The result is a reference-prior enrichment, not DCT-isolated spaceflight phosphoproteomics.

Best claim label:

> DCT1-prioritized phospho suppression is supported at reference-prior resolution.

Avoid:

> DCT1-specific phosphoproteomic remodeling.

## H2 Composition-Aware Robustness

This follow-up asks the reviewer-obvious question: is the DCT1-prior phosphosite result just whole-kidney DCT dilution, parent-protein abundance, or endothelial/stromal expansion?

Short answer: no for the DCT1-high top-decile subset, but yes/unclear for a broad continuous DCT1-gradient claim.

Single-site phosphosite model ladder:

| Model | Adjustment | DCT1 top-decile suppression OR | q value | Read |
| --- | --- | ---: | ---: | --- |
| M0 | raw flight effect | 1.40 | 6.35e-7 | raw DCT1-high subset enrichment |
| M1 | parent protein abundance | 1.55 | 2.14e-10 | not explained by parent protein abundance |
| M2 | DCT score | 1.36 | 2.17e-5 | survives estimated DCT abundance adjustment |
| M3 | endothelial + stromal scores | 1.24 | 0.00545 | attenuated, still positive |
| M4 | parent protein + DCT + endothelial + stromal | 1.30 | 0.00158 | full conservative model passes |
| M5 | parent protein + composition PC1 | 1.38 | 6.53e-5 | PC-based composition adjustment passes |

The top-quartile result is less stable. In the full M4 model, DCT1 top quartile gives OR = 1.05, q = 0.282, while DCT1 top decile remains positive. This means the robust signal is concentrated in the most DCT1-prioritized parent-gene subset.

The more aggressive continuous models do not pass:

| Test | Key term | Coefficient | q value | Read |
| --- | --- | ---: | ---: | --- |
| Two-stage M0 DCT1-only | DCT1 prior | -0.00110 | 0.877 | weak negative, not significant |
| Two-stage M4 full | DCT1 prior | +0.00274 | 0.877 | no continuous suppression gradient |
| One-shot LM0 site-fixed | flight x DCT1 prior | -0.00145 | 0.665 | direction negative, weak |
| One-shot LM4 full site-fixed | flight x DCT1 prior | -0.00200 | 0.665 | direction negative, weak |

Best manuscript claim:

> The DCT1-high top-decile enrichment of flight-suppressed phosphosites persisted after conservative adjustment for parent protein abundance and estimated DCT/endothelial/stromal composition.

Avoid:

> We deconvolved DCT-specific phosphoproteomics from whole-kidney data.

Interpretation boundary:

The animal-level DCT, endothelial, and stromal scores are estimated from bulk RNA. They are useful sensitivity adjustments, but not true deconvolution of phosphoproteomic cell of origin. Also, composition adjustment may remove part of the biological pathway if endothelial/stromal expansion or DCT loss is itself on the causal path from flight to phosphosite suppression.

## PXD001729 DCT-Lineage Phosphoproteomics

PXD001729 loaded and mapped, but the anti-alignment hypothesis did not pass in a useful way.

Key result:

| Comparison | Shared sites | Shared genes | Cosine | 95% bootstrap CI | Read |
| --- | ---: | ---: | ---: | --- | --- |
| all shared single sites | 60 | 49 | -0.026 | -0.310 to 0.287 | essentially null / unstable |
| transport target shared sites | 0 | 0 | NA | NA | not testable |

Target coverage boundary:

- PXD001729 contains phosphosites for Wnk1, Wnk4, Stk39, Oxsr1, and Nedd4l.
- It does not provide shared changed target sites for the key NCC/SPAK/WNK transport subset used in OSD-462.
- Slc12a3/NCC did not map as a shared transport phosphosite for this anti-alignment test.

Manuscript consequence:

PXD001729 should be used as a DCT-lineage phosphoproteome plausibility reference, not as evidence that spaceflight opposes dDAVP/NCC activation. The anti-alignment claim should be removed or moved to a negative/coverage-boundary paragraph.

Best wording:

> The mpkDCT dDAVP reference confirms that DCT-lineage cells carry a rich phosphoproteomic signaling space, but public overlap with OSD-462 is insufficient to test a targeted vasopressin anti-alignment model for NCC/SPAK/WNK sites.

## KLHL3/CUL3 Mechanism Check

This remains not testable from the public data.

Key coverage:

| Gene | OSD-462 phosphosites | PXD phosphosites | Protein effect | Read |
| --- | ---: | ---: | ---: | --- |
| Klhl3 | 0 | 0 | +0.013 | no phosphosite evidence |
| Cul3 | 1 | 0 | -0.001 | no useful turnover signal |
| Wnk1 | 12 | 12 | -0.024 | WNK1 phospho observed, protein mostly flat |
| Wnk4 | 8 | 1 | +0.087 | WNK4 phospho observed, protein not depleted |
| Stk39/SPAK | 5 | 3 | +0.081 | phospho suppression despite protein not depleted |
| Slc12a3/NCC | 16 | 0 | +0.089 | phospho suppression with protein flat/up |

Interpretation:

The public data do not distinguish KLHL3/CUL3-mediated WNK degradation from ionic/osmotic WNK-SPAK suppression. The lack of KLHL3 phosphosite coverage and lack of ubiquitinomics are decisive limitations. If anything, flat WNK/SPAK/NCC protein abundance makes a pure turnover story less directly supported by these data.

Best manuscript use:

Name KLHL3/CUL3-WNK turnover as a candidate future mechanism, not as a result.

## H3: Mediation-Specified Causal Skeleton

The mediation model ran as an approximate Bayesian linear-model fallback, not a full `brms`/SEM fit. The result is useful as causal-structure specification, not causal proof.

Outcome: animal-level NCC regulatory phosphorylation score.

Predictor: flight status.

Candidate mediators: endothelial, stromal/fibroblast, DCT identity, and matrix/endothelial composite scores.

| Mediator | Indirect median | 95% interval | P(indirect < 0) | Read |
| --- | ---: | --- | ---: | --- |
| endothelial | -0.432 | -1.019 to -0.036 | 0.984 | consistent with negative indirect path |
| stromal/fibroblast | -0.295 | -0.843 to 0.025 | 0.964 | suggestive, interval crosses 0 |
| DCT identity | -0.187 | -0.658 to 0.075 | 0.922 | composition-linked, interval crosses 0 |
| matrix/endothelial composite | -0.408 | -0.990 to -0.033 | 0.985 | consistent with negative indirect path |

The direct flight effect remains negative after mediator adjustment in all models. That matters: the data look more like partial mediation or shared composition-linked remodeling than a simple "remodeling fully explains phospho suppression" model.

Future-n planning:

| Mediator | n for approx 80% interval-exclusion power |
| --- | ---: |
| endothelial | about 40 animals |
| matrix/endothelial composite | about 40 animals |
| stromal/fibroblast | about 60 animals |
| DCT identity | about 80 animals |

These are approximate Sobel-style simulations using the observed effect scale. They should be used for planning, not as formal design guarantees.

Best manuscript claim:

> Animal-matched RR-10 data are consistent with a negative indirect path through endothelial or matrix/endothelial remodeling, while the persistent direct flight effect and bulk-tissue composition confounding prevent causal interpretation.

Avoid:

> Endothelial remodeling mediates NCC dephosphorylation.

## External Spatial Reference Contextualization

This analysis is now implemented as an external kidney injury/repair spatial reference, not as spaceflight spatial validation.

Downloaded sources:

- GSE269622 Visium Space Ranger archives: sham, 4 h, 12 h, day 2, day 14, and week 6 mouse IRI kidneys.
- GSE269719 processed Xenium AnnData object: 1,374,915 cells, 299 genes, checksum verified against the figshare MD5.

Primary Visium result:

| IRI timepoint | OSD-462 pathway cosine vs sham | RRRM-2 pathway cosine vs sham | OSD-462 gene-level cosine vs sham | RRRM-2 gene-level cosine vs sham |
| --- | ---: | ---: | ---: | ---: |
| 4 h | 0.627 | 0.709 | 0.041 | 0.079 |
| 12 h | 0.780 | 0.699 | 0.063 | 0.054 |
| day 2 | 0.792 | 0.667 | 0.061 | 0.090 |
| day 14 | 0.797 | 0.720 | 0.083 | 0.070 |
| week 6 | 0.801 | 0.687 | 0.028 | 0.275 |

Read:

The spaceflight remodeling vector aligns strongly with IRI repair-stage pathway programs, especially day 14 to week 6, but genome-wide gene-vector cosines are much smaller. This supports a pathway-level injury/repair contextualization, not a claim that spaceflight kidney is globally transcriptionally identical to IRI.

Top OSD-462 spatial niche matches by pathway cosine:

| Timepoint | Niche | Cosine |
| --- | --- | ---: |
| day 14 | injured tubule-enriched spots | 0.827 |
| day 14 | DCT-adjacent spots | 0.825 |
| week 6 | DCT-adjacent spots | 0.814 |
| day 14 | fibro-inflammatory repair-enriched spots | 0.807 |
| week 6 | injured tubule-enriched spots | 0.802 |
| week 6 | endothelial-enriched spots | 0.802 |

Secondary Xenium inventory:

The Xenium object is appropriate for annotation/neighborhood context only because it is a 299-gene targeted panel, not a whole-transcriptome projection source. The object includes `celltype_plot`, `time`, and `CN` annotations. Relevant counts include DCT cells = 48,774, endothelial cells = 129,261, fibroblasts = 160,328, injured proximal tubule cells = 146,097, distal tubule niche cells = 77,103, and fibro-inflammatory niche cells = 105,536.

Best manuscript wording:

> Because no RR-10 spatial transcriptomic dataset is publicly available, we used GSE269622/GSE269719 only as an external kidney injury/repair spatial reference atlas. We projected the spaceflight bulk RNA remodeling vector onto Visium-derived spatial gene-expression programs and used Xenium annotations to contextualize whether aligned programs corresponded to tubule-injury, interstitial, endothelial, fibro-inflammatory, or DCT-adjacent neighborhoods. This analysis does not localize the spaceflight lesion directly; it asks which known kidney injury/repair spatial niches the bulk spaceflight remodeling state most resembles.

Avoid:

> Spatial validation of the spaceflight DCT lesion.

## Manuscript Action

The paper can now be reframed around three core hypotheses plus a spatial-reference supplement, with tighter claim labels:

1. H1: recurrent matrix/endothelial-high and DCT-low RNA remodeling state.
2. H2: DCT1-prioritized phosphosite suppression at reference-prior resolution, with top-decile enrichment surviving parent-protein and composition-aware sensitivity models.
3. H3: remodeling-linked versus parallel-response causal structure specified by animal-matched mediation models.
4. Supplement: external IRI spatial reference contextualization, explicitly not spaceflight spatial validation.

This is more novel than v10 because it adds a falsifiable DCT1 prior and an explicit causal skeleton. It is still not a true mechanism paper because the decisive experiments are absent: DCT-enriched spaceflight phosphoproteomics, ubiquitinomics, perturbation/rescue, and actual RR-10 spatial transcriptomics.

Recommended title:

> DCT1-prioritized phospho-regulatory suppression in spaceflight kidney remodeling

Recommended abstract sentence:

> By integrating a DCT-enriched snRNA-seq reference, a DCT-lineage phosphoproteomic reference, and animal-matched RR-10 multi-omics, we refine the established NCC/SPAK lesion as a DCT1-prioritized whole-kidney phosphoproteomic phenotype and test whether matrix/endothelial remodeling is statistically linked to, or parallel with, transporter deactivation.

## Evidence Ladder

| Layer | Dataset | Result | Claim strength |
| --- | --- | --- | --- |
| RNA remodeling | RRRM-2, OSD-513, OSD-462 | recurrent matrix/endothelial-high and DCT-low direction | strong |
| Protein abundance | OSD-462 | targeted protein concordance null | strong negative |
| Phospho lesion | OSD-462 | NCC/SPAK/WNK regulatory phospho down | strong, established anchor |
| DCT subtype prior | GSE228367 -> OSD-462 | suppressed phosphosites enriched in DCT1-high parent genes | new, reference-prior positive |
| Composition-aware H2 | OSD-462 matched animals + RNA compartment scores | DCT1 top-decile enrichment survives parent-protein and composition adjustment; continuous interaction does not pass | robust subset positive |
| DCT phospho reference | PXD001729 -> OSD-462 | anti-alignment not stable/testable | negative/coverage boundary |
| KLHL3/CUL3 | OSD-462 + PXD001729 | not distinguishable without ubiquitinomics | negative/coverage boundary |
| Mediation | OSD-462 matched animals | endothelial/composite indirect path consistent with negative mediation | suggestive, composition-confounded |
| Spatial | GSE269622/GSE269719 | pathway-level alignment to IRI repair-stage and DCT-adjacent/injured-tubule niches | contextual supplement, not validation |

## Sources Grounding The Interpretation

Local sources in `transcriptomic_texts/` used:

- `s41467-024-49212-1.pdf`: Cosmic Kidney disease paper.
- `s41526-025-00465-0 (3).pdf`: spaceflight strain-dependent kidney transcriptomics.
- `s41598-023-50195-0.pdf`: kidney fibrosis single-cell injury reference.

Online/source metadata checked:

- Cosmic Kidney disease, Nature Communications 2024: https://www.nature.com/articles/s41467-024-49212-1
- GSE228367 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228367
- Su et al. DCT-enriched snRNA-seq, JASN 2024: https://pubmed.ncbi.nlm.nih.gov/38238903/
- Cheng et al. mpkDCT dDAVP phosphoproteomics, Scientific Reports 2015: https://www.nature.com/articles/srep12829
- WNK-SPAK-NCC and CUL3/KLHL3 review context: https://www.nature.com/articles/s41440-020-0437-x
- KSEA method reference: https://academic.oup.com/bioinformatics/article/33/21/3489/3892392
- GSE269622 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269622
- GSE269719 GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269719
- IRI Xenium/Visium spatial reference, Nature Communications 2025: https://www.nature.com/articles/s41467-025-62599-9
- IRI spatial analysis code: https://github.com/TheHumphreysLab/IRI_Xenium_Visium
