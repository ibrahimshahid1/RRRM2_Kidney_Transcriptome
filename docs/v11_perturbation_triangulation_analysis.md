# V11 Perturbation Triangulation Analysis

## Purpose

This analysis was added to ask whether public perturbation datasets can move the
manuscript closer to the unresolved mechanistic question:

> Why is WNK-SPAK-NCC regulatory phosphorylation suppressed in mouse spaceflight kidney?

The existing v11 manuscript already showed that flight-suppressed OSD-462
whole-kidney phosphosites are enriched among parent genes that are DCT1-high in
a native DCT reference. That result refines the lesion toward DCT1-prioritized
parent genes, but it does not identify the upstream signal that suppresses the
WNK-SPAK-NCC pathway.

The perturbation-triangulation pass therefore asked whether the spaceflight
signature resembles any public mechanistic reference:

- potassium/chloride-WNK regulation
- vasopressin/cAMP signaling
- injury-repair remodeling
- KLHL3/CUL3-mediated WNK turnover

The goal was not to prove mechanism. The goal was to decide whether any public
perturbation gives a strong enough match to become a main mechanistic result, or
whether the manuscript should report the public-data ceiling clearly and use the
results to motivate additional experimental data.

## Decision Rule

The highest-priority test was GSE228367 low-potassium DCT anti-alignment.

Low potassium normally activates the DCT/NCC-WNK adaptive program. Spaceflight
shows suppression of NCC/SPAK/WNK regulatory phosphorylation. Therefore, a strong
negative alignment between the low-K DCT response and the spaceflight DCT
transport signature would support a potassium/chloride-like DCT deactivation or
failed low-K adaptation hypothesis.

The planned decision rule was:

| Low-K result | Action |
| --- | --- |
| Strong anti-alignment, cosine less than -0.3, CI excluding zero, and target genes moving oppositely | Promote low-K anti-alignment to a primary mechanistic result alongside DCT1 enrichment. |
| Weak, unstable, or null anti-alignment | Treat low-K as mechanism-triage evidence, not a primary result. Move emphasis to parent-protein-normalized phosphosite enrichment and spatial prediction. |

## Analysis Stack

| Priority | Analysis | Main question | Output role |
| --- | --- | --- | --- |
| 1 | GSE228367 low-K DCT anti-alignment | Does spaceflight oppose the native DCT low-potassium adaptive program? | Mechanism triage; possible primary result if strong. |
| 2 | OSD-462 occupancy-normalized phosphosite effects | Does DCT1-high parent-gene enrichment survive parent-protein normalization? | Main robustness upgrade to H2. |
| 3 | GSE269622 IRI DCT transport spatial check | Does injury-repair remodeling co-occur with DCT transport suppression in DCT-adjacent spatial niches? | Spatial prediction/context. |
| 4 | PXD001729 dDAVP/mpkDCT phosphoproteomics | Does the flight phosphosite lesion resemble vasopressin/cAMP remodeling? | Secondary reference; expected to be weak. |
| 5 | KLHL3/CUL3/WNK turnover check | Is there public evidence for KLHL3/CUL3-mediated WNK turnover? | Future-data rationale. |

## Branch 1: GSE228367 Low-K DCT Anti-Alignment

### Why this was done

GSE228367 is the most relevant public transcriptomic perturbation for the
upstream WNK-SPAK-NCC question because it profiles native DCT-enriched nuclei
under normal potassium and potassium-restricted conditions. If spaceflight
suppression of NCC/SPAK/WNK reflects altered potassium/chloride sensing, then
the spaceflight DCT transport program might oppose the canonical low-K DCT
activation response.

### Inputs

- `data/external/dct_reference/GSE228367/GSE228367_RAW.tar`
  - GEO raw archive containing filtered 10x matrices:
    - `NK1`, `NK2`, `NK3`
    - `KD1`, `KD2`, `KD3`
- `data/results/run_20260526_v11_dct1_phospho_mediation/dct_prior/gse228367_dct1_vs_dct2_de.tsv`
  - Native DCT1/DCT2 prior used to define DCT1-prior and DCT2-prior gene subsets.
- `data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_flight_effects.tsv`
  - OSD-462 RNA flight effect vector.
- `data/results/run_20260519_000547_2500g/contrast_vectors/mechanism_axis/tubulointerstitial_state/lar_reversal/lar_reversal_gene_scatter.tsv`
  - RRRM-2 ISS-T RNA flight effect vector.
- `data/results/run_20260522_regulator_activity/rna_effects/OSD513_gene_stat.tsv`
  - OSD-513 RNA flight statistic vector.

### How it was done

The six GSE228367 filtered 10x matrices were loaded, library-size normalized,
log-transformed, and summarized as sample-level pseudobulk mean expression.

For each gene:

```text
lowK_effect_g = mean(log-normalized expression in KD1-KD3)
              - mean(log-normalized expression in NK1-NK3)
```

The low-K vector was compared with each spaceflight RNA vector using:

- cosine similarity
- Spearman correlation
- bootstrap confidence intervals by resampling genes

The comparisons were run across:

- all genes overlapping the DCT prior
- DCT1-prior top-decile genes
- DCT2-prior top-decile genes
- DCT1-prior top-quartile genes
- DCT2-prior top-quartile genes
- focused DCT/WNK transport target genes

The focused transport target set was:

```text
Slc12a3, Stk39, Oxsr1, Wnk1, Wnk4, Kcnj10, Kcnj16,
Clcnkb, Bsnd, Klhl3, Cul3, Ppp1ca, Ppp1r1a, Calb1
```

### Main results

The strongest biologically relevant signal was in the 14 transport target genes,
but it did not meet the pre-specified primary-result threshold.

| Spaceflight vector | Gene subset | Cosine | 95% CI | Spearman rho | q | Interpretation |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| OSD-462 RNA | transport targets, n=14 | -0.286 | -0.752 to 0.458 | -0.420 | 0.218 | Directional anti-alignment, not decisive |
| RRRM-2 ISS-T RNA | transport targets, n=14 | -0.311 | -0.780 to 0.325 | -0.341 | 0.306 | Directional anti-alignment, not decisive |
| OSD-513 RNA | transport targets, n=14 | -0.146 | -0.768 to 0.506 | -0.169 | 0.695 | Weak |

For OSD-462 RNA, 10 of the 14 focused DCT/WNK genes moved opposite the low-K
response, including `Slc12a3`, `Wnk1`, `Klhl3`, and `Cul3`. However, the small
target-gene set had wide bootstrap intervals, and most individual low-K effects
were not FDR-significant in the three-versus-three pseudobulk design.

### Interpretation

The low-K comparison gives partial directional support at the DCT/WNK target
gene level, but it should not be promoted to a primary mechanism claim. The
clean manuscript wording is:

> A native DCT low-potassium reference showed directional anti-alignment with
> the spaceflight DCT/WNK target-gene signature, but the effect was not stable
> enough to identify potassium/chloride-WNK regulation as the upstream mechanism.

### Output files

- `data/results/run_20260526_v11_dct1_phospho_mediation/perturbation/lowk_dct_pseudobulk_effects.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/perturbation/lowk_dct_alignment_summary.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/perturbation/lowk_target_gene_table.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/perturbation/lowk_dct1_dct2_specificity.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/perturbation/lowk_alignment_verdict.json`

## Branch 2: OSD-462 Parent-Protein-Normalized Phosphosite Effects

### Why this was done

The main v11 H2 claim is that flight-suppressed whole-kidney phosphosites are
enriched among DCT1-high parent genes, especially in the DCT1 top decile.

A reviewer could ask whether this is just parent protein abundance. To address
that directly, each phosphosite flight effect was normalized against the flight
effect of its parent protein.

### Inputs

- `data/results/run_20260526_v11_dct1_phospho_mediation/dct_prior/osd462_phosphosite_dct1_prior.tsv`
  - OSD-462 phosphosite effects with DCT1/DCT2 parent-gene priors and parent
    protein flight effects.

### How it was done

For each phosphosite:

```text
occupancy_effect_s = phosphosite_flight_effect_s
                   - parent_protein_flight_effect_g
```

This is a parent-protein-normalized robustness metric. It is not a direct
biochemical phosphorylation stoichiometry measurement, because it uses
whole-kidney TMT protein and phosphosite effect estimates rather than direct
site occupancy.

DCT1 enrichment was then retested for:

- DCT1 top-decile parent genes
- DCT1 top-quartile parent genes
- nominal p < 0.05 suppressed sites
- strict q < 0.10 suppressed sites
- anchor-gene-excluded sensitivity
- NCC-site-excluded sensitivity
- single-site-per-gene sensitivity

### Main results

This was the strongest new result from the triangulation pass.

| Analysis | DCT1 subset | OR | 95% CI | q |
| --- | --- | ---: | ---: | ---: |
| occupancy p < 0.05 | top decile | 1.52 | 1.36 to 1.70 | 3.73e-12 |
| occupancy p < 0.05 | top quartile | 1.17 | 1.06 to 1.28 | 0.00113 |
| strict q < 0.10 | top decile | 1.61 | 1.42 to 1.82 | 3.73e-12 |
| strict q < 0.10 | top quartile | 1.21 | 1.08 to 1.34 | 0.000683 |
| exclude anchor genes | top decile | 1.51 | 1.34 to 1.69 | 9.75e-12 |
| exclude NCC sites | top decile | 1.51 | 1.35 to 1.69 | 9.75e-12 |
| single-site per gene | top decile | 1.50 | 1.19 to 1.90 | 8.42e-4 |

Target-site examples:

| Site | Raw phosphosite effect | Parent protein effect | Occupancy-normalized effect |
| --- | ---: | ---: | ---: |
| `Slc12a3/NCC Thr53` | -0.851 | +0.089 | -0.940 |
| `Slc12a3/NCC Thr65` | -0.790 | +0.089 | -0.879 |
| `Slc12a3/NCC Thr68` | -0.794 | +0.089 | -0.883 |
| `Slc12a3/NCC Thr89` | -0.930 | +0.089 | -1.019 |
| `Stk39/SPAK Ser366` | -0.520 | +0.081 | -0.601 |
| `Stk39/SPAK Ser383` | -0.793 | +0.081 | -0.874 |

### Interpretation

This supports the statement:

> The DCT1 top-decile enrichment of flight-suppressed whole-kidney phosphosites
> persists after parent-protein normalization.

This is stronger than the original enrichment result because it directly
addresses the possibility that the DCT1 signal is trivially explained by parent
protein abundance.

The wording should still avoid saying "DCT1 phosphoproteomics" or "cell of
origin." The result is still parent-gene-prioritized whole-kidney
phosphoproteomics.

### Output files

- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_occupancy/h2_occupancy_site_effects.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_occupancy/h2_occupancy_dct1_enrichment.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_occupancy/h2_occupancy_target_sites.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_occupancy/h2_occupancy_verdict.json`

## Branch 3: GSE269622 IRI Spatial DCT Transport Check

### Why this was done

The existing v11 spatial-reference analysis showed that the spaceflight bulk RNA
state resembles late injury/repair spatial programs in an external IRI Visium
atlas. That helped with context, but it did not ask whether DCT transport
suppression occurs in those spatial neighborhoods.

This branch asked whether external IRI repair niches show DCT transport
suppression in DCT-adjacent regions.

### Inputs

- `data/external/spatial_reference/GSE269622_Visium/`
  - Visium Space Ranger archives for:
    - sham
    - 4 h IRI
    - 12 h IRI
    - day 2 IRI
    - day 14 IRI
    - week 6 IRI
- Existing marker definitions in `src/v11/spatial_reference_projection.py`

### How it was done

For each Visium spot, a DCT transport score was computed as the mean standardized
expression of:

```text
Slc12a3, Stk39, Wnk1, Wnk4, Oxsr1,
Kcnj10, Kcnj16, Clcnkb, Bsnd, Calb1
```

Spots were grouped into:

- all spots
- DCT-marker-high spots
- DCT-adjacent spots
- non-DCT-adjacent spots
- marker-enriched injured tubule, fibroblast/interstitial, endothelial, immune,
  and fibro-inflammatory repair niches

DCT-adjacent spots were defined as spots located near DCT-marker-high spots, but
not themselves in the DCT-marker-high top quartile.

### Main results

Late IRI showed lower DCT transport score in DCT-adjacent spots compared with
sham DCT-adjacent spots:

| Spatial group | Timepoint | Difference versus sham | p | Spots |
| --- | --- | ---: | ---: | ---: |
| DCT-adjacent spots | day 14 | -0.0436 | 0.0050 | 932 vs 1134 |
| DCT-adjacent spots | week 6 | -0.0344 | 0.0182 | 1574 vs 1134 |

In contrast, DCT-marker-high spots showed higher DCT transport score at day 14
and week 6:

| Spatial group | Timepoint | Difference versus sham | p |
| --- | --- | ---: | ---: |
| DCT-marker-high spots | day 14 | +0.199 | 2.25e-15 |
| DCT-marker-high spots | week 6 | +0.132 | 2.70e-9 |

### Interpretation

The spatial-reference result should be framed carefully:

> In an external IRI Visium atlas, late repair-stage DCT-adjacent neighborhoods
> showed lower DCT transport score, while DCT-marker-high spots retained or
> increased DCT transport score. This suggests a testable spatial prediction:
> future spaceflight kidney sections should determine whether DCT/NCC-WNK loss
> is intrinsic to DCT cells, concentrated in DCT-adjacent repair neighborhoods,
> or both.

This is not spaceflight spatial validation. It is a spatial reference that makes
the future experiment more concrete.

### Output files

- `data/results/run_20260526_v11_dct1_phospho_mediation/spatial_reference/visium_dct_transport_by_niche.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/spatial_reference/visium_dct_transport_timepoint_summary.tsv`
- `data/results/run_20260526_v11_dct1_phospho_mediation/spatial_reference/visium_dct_transport_spot_scores.tsv.gz`
- `data/results/run_20260526_v11_dct1_phospho_mediation/spatial_reference/visium_dct_transport_verdict.json`

## Branch 4: PXD001729 dDAVP/mpkDCT Phosphoproteomics

### Why this was included

PXD001729 is a cultured DCT-lineage phosphoproteomic dataset after acute dDAVP
stimulation. It is relevant because vasopressin/cAMP can regulate distal nephron
transport signaling, but it is not native kidney tissue and not spaceflight.

### Result

The shared phosphosite map was too limited to support a strong vasopressin/cAMP
mechanism:

- shared single phosphosites: 60
- shared transport-target phosphosites: 0
- missing transport targets included:
  - `Slc12a3`
  - `Stk39`
  - `Oxsr1`
  - `Wnk1`
  - `Wnk4`
  - `Klhl3`
  - `Cul3`

### Interpretation

The dDAVP/mpkDCT reference does not currently explain the OSD-462 spaceflight
transporter phosphosite lesion.

### Output files

- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_pxd/h2_pxd001729_ddavp_antialignment_verdict.json`
- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_pxd/pxd001729_osd462_shared_phosphosites.tsv`

## Branch 5: KLHL3/CUL3/WNK Turnover

### Why this was included

WNK kinases can be regulated by KLHL3/CUL3-mediated degradation. If public
OSD-462 data contained useful KLHL3/CUL3/WNK abundance, phosphorylation, or
ubiquitin-remnant evidence, it might distinguish WNK turnover from ionic/osmotic
suppression.

### Result

- KLHL3 Ser433 was not detected in the OSD-462 targeted phosphosite table.
- Public OSD-462 does not include ubiquitinomics.
- Current public data cannot distinguish KLHL3/CUL3-mediated WNK turnover from
  ionic/osmotic WNK-SPAK suppression.

### Interpretation

This remains a future-data question. A real test would need DCT-enriched protein
or phosphoprotein data, WNK protein turnover readouts, ubiquitin-remnant data, or
targeted perturbation.

### Output files

- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_klhl3/h2_klhl3_cul3_interpretation.md`
- `data/results/run_20260526_v11_dct1_phospho_mediation/h2_klhl3/h2_klhl3_cul3_effects.tsv`

## Figure Output

A compact perturbation-triangulation figure was generated:

- `data/results/run_20260526_v11_dct1_phospho_mediation/figures/v11/v11_perturbation_triangulation.png`
- `data/results/run_20260526_v11_dct1_phospho_mediation/figures/v11/v11_perturbation_triangulation.pdf`

Panels:

| Panel | Content |
| --- | --- |
| A | Low-K DCT response versus OSD-462, OSD-513, and RRRM-2 RNA vectors. |
| B | Focused DCT/WNK target-gene direction heatmap. |
| C | Parent-protein-normalized DCT1 enrichment odds ratios. |
| D | IRI Visium DCT transport score by timepoint and spatial context. |

## Runner Integration

The following modules were added to the v11 phase in `src/run_all_phases.py`:

```text
src.v11.perturbation_gse228367_lowk
src.v11.h2_occupancy_normalized_phospho
src.v11.spatial_dct_transport_check
src.v11.perturbation_triage_summary
```

The full v11 dry-run command successfully lists these steps:

```bash
./venv/bin/python -m src.run_all_phases \
  --v11-only \
  --skip-r \
  --dry-run \
  --run-id run_20260526_v11_dct1_phospho_mediation
```

## Final Interpretation

The latest triangulation pass changes the manuscript in three ways.

First, it prevents overclaiming. The obvious potassium/chloride-WNK public
reference does not cleanly explain the spaceflight phenotype. Low-K DCT
anti-alignment is directional for transport targets but not stable enough to
serve as a primary mechanistic result.

Second, it strengthens H2. Parent-protein-normalized phosphosite effects remain
enriched among DCT1 top-decile parent genes. This directly addresses the concern
that DCT1 enrichment is merely a parent protein abundance artifact.

Third, it sharpens the spatial prediction. External IRI Visium data suggest that
late repair-stage DCT-adjacent neighborhoods can show lower DCT transport score.
That makes the future spaceflight spatial experiment clearer: test whether
DCT/NCC-WNK suppression is DCT-intrinsic, DCT-adjacent remodeling-associated, or
both.

The manuscript should therefore emphasize:

> Spaceflight-associated whole-kidney phosphosite suppression remains enriched
> among DCT1-high parent genes after parent-protein normalization. Public
> perturbation references do not identify a single upstream WNK-suppressing
> mechanism, but low-K DCT and IRI spatial references generate testable
> mechanistic and spatial predictions for future animal-matched spaceflight
> kidney experiments.

