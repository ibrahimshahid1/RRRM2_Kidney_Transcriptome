# V11 Execution Research Plan

Status: execution plan / pre-registration draft

Date: May 2026

Source: `/Users/ibrahimshahid/Downloads/spaceflight_kidney_v11_research_plan.pdf`

Working title: DCT1-anchored phospho-regulatory suppression in spaceflight kidney

Core framing: v11 should move the project from a cross-omic context map toward a
DCT1-anchored, phosphoproteome-informed, mediation-specified mechanistic
hypothesis. It must not claim causal proof, DCT-isolated spaceflight
phosphoproteomics, or new discovery of NCC dephosphorylation.

## 1. Purpose

v10 is a rigorous context-mapping study around an established phenotype. It shows
that a recurrent matrix-high / DCT-low RNA state surrounds the established
spaceflight NCC/SPAK phosphorylation lesion, that the RNA state does not carry
through as protein abundance, and that the transporter lesion resolves at the
phospho-activity layer.

The v11 goal is to answer the question v10 leaves open:

Can the established transporter lesion be refined into a DCT1/NCC-centered
phospho-regulatory phenotype, and can the relationship between the kidney-wide
remodeling state and NCC/SPAK phospho-suppression be formally specified?

Every v11 analysis must target one of two gaps:

1. Causal structure: is remodeling statistically upstream of, downstream of, or
   parallel to transporter deactivation?
2. Functional consequence / biological specificity: is the transporter lesion
   DCT1/NCC-centered rather than generic whole-kidney phospho noise?

## 2. Three-hypothesis architecture

### H1: Spaceflight induces a recurrent kidney-wide remodeling state

Status: already established in v10. Preserve, do not inflate.

Existing evidence:

- RRRM-2 ISS-T, OSD-513, and OSD-462 RNA recur along a matrix-high / DCT-low
  direction.
- OSD-462 RNA vs RRRM-2 ISS-T pathway-vector cosine is 0.87 with 95% CI
  approximately 0.65 to 0.90.
- OSD-462 protein abundance does not match the RNA remodeling direction.
- The measurable transporter lesion appears at the NCC/SPAK/WNK
  phospho-activity layer.

V11 addition:

- Situate the endothelial/stromal component as a composition axis of the DCT-low
  signal.
- Optionally project the RNA remodeling state onto an external IRI spatial atlas,
  but only as a spatial reference, not as spaceflight spatial evidence.

Claim label:

- Supported / established / reproduced.
- Spatially referenced only if optional Step 7 is completed.

### H2: The transporter lesion is DCT1/NCC-centered

Status: primary v11 upgrade.

Core logic:

1. GSE228367 provides a native DCT RNA reference with control DCT1 and DCT2
   nuclei.
2. PXD001729 provides a cultured DCT-lineage phosphoproteomic reference under
   dDAVP/vasopressin stimulation.
3. OSD-462 provides the actual whole-kidney spaceflight phosphoproteomic anchor.

Main test:

For each OSD-462 phosphosite, assign a DCT1 enrichment score from the GSE228367
DCT1 vs DCT2 reference. Test whether OSD-462 flight-suppressed phosphosites are
enriched in DCT1-high parent genes.

Important wording:

- Say "DCT1-prioritized", "DCT1-informed", or "DCT1-centered at the reference
  prior level".
- Do not say "DCT1-isolated" or "DCT1-proven".

Claim label if positive:

- DCT1-localized at RNA reference level.
- DCT-lineage plausible at phospho reference level.
- Spaceflight-suppressed at whole-kidney phospho level.
- DCT specificity is inferred, not directly measured.

Claim label if null:

- DCT1 specificity not supported.
- Report null honestly and do not tune thresholds to rescue H2.

### H3: Remodeling state and transporter lesion are formally related

Status: hypothesis specification, not mechanistic proof.

Question:

Does the endothelial/stromal remodeling state statistically mediate the relation
between flight and NCC regulatory phosphorylation, or are remodeling and NCC
phospho-suppression better treated as parallel flight responses?

What the model can deliver:

- Posterior estimate of an indirect path through endothelial/stromal state.
- Wide credible intervals expected at n=20.
- Directional probability for the indirect path.
- Power / sample-size estimate for a future experiment.

What it cannot deliver:

- Causal proof.
- Time ordering.
- Separation of true biological mediation from bulk tissue composition dilution.
- DCT cell-autonomous phospho suppression.

Claim label:

- "Mediation-specified" or "causal structure bounded".
- Avoid "mediation-tested" in title or abstract unless wording clearly says the
  test is underpowered and observational.

## 3. Dataset inventory and local paths

### Existing v10 datasets

| Dataset | Role | Local status |
| --- | --- | --- |
| OSD-771 / RRRM-2 | discovery and factorial bridge | existing pipeline outputs |
| OSD-513 | independent RNA recurrence cohort | existing local data |
| OSD-253 | strain / TLR4 context cohort | existing local data |
| OSD-462 / RR-10 | matched RNA, proteome, phosphoproteome anchor | existing local data |
| Tabula Muris Senis kidney | external aging reference | existing local data |

### New v11 datasets

| Dataset | Role | Local path | Status |
| --- | --- | --- | --- |
| GSE228367 | DCT1/DCT2 RNA reference | `data/external/dct_reference/GSE228367/` | DCT1 and DCT2 control RDS files downloaded |
| PXD001729 | mpkDCT dDAVP phosphoproteomics reference | `data/external/phosphoproteomics/PXD001729/` | supplementary PDF and Tables S2-S6 downloaded |
| GSE269719 / GSE269622 | optional spatial IRI reference | `data/external/spatial_reference/GSE269719_GSE269622/` | placeholder only; raw archives not downloaded |

Critical framing rule:

GSE228367 and PXD001729 are reference scaffolds, not spaceflight evidence. Every
sentence using these datasets should say "in the DCT reference dataset" or "in
DCT-lineage cells" rather than implying direct spaceflight DCT measurement.

## 4. Pre-registration and run discipline

Before running new v11 analyses:

1. Treat this file as the pre-registration draft.
2. Freeze analysis definitions for H1, H2, and H3.
3. Compute and record a SHA-256 checksum of this plan in the first v11 run
   manifest.
4. Use one v11 run root:

   `data/results/run_YYYYMMDD_v11_dct1_phospho_mediation/`

5. Every output table must include a manifest with:

   - input file paths
   - input SHA-256 hashes where feasible
   - package versions
   - filtering thresholds
   - random seeds
   - analysis date

Do not inspect H2 results and then change thresholds. Threshold changes after
inspection must be labeled exploratory sensitivity analyses.

## 5. Step 0: baseline lock from v10

Purpose:

Fix the v10 quantities that v11 uses as inputs. Do not re-derive the v10
biological story while implementing H2/H3.

Inputs:

- OSD-462 phosphosite flight effects from the existing OSD-462 anchor run.
- OSD-462 protein effects from the existing OSD-462 anchor run.
- OSD-462 per-animal NCC regulatory phosphorylation score from phenotype anchor.
- OSD-462 per-animal compartment scores from cell-type marker panels.
- Existing RRRM-2 / OSD-513 / OSD-462 pathway-vector recurrence outputs.

Checks:

- OSD-462 RNA recurrence gate remains PASS.
- Protein concordance remains null for matrix and DCT gene sets.
- NCC regulatory phosphosites remain specifically suppressed while total NCC is
  flat.
- KSEA WNK and SPAK/OSR1 positive control remains inferred activity down.

Outputs:

- `v11_baseline_lock_summary.tsv`
- `v11_baseline_input_manifest.json`

## 6. Step 1: verify downloaded external datasets

### Step 1a: GSE228367 verification

Local inputs:

- `data/external/dct_reference/GSE228367/GSE228367_CONTROL_DCT1.rds.gz`
- `data/external/dct_reference/GSE228367/GSE228367_CONTROL_DCT2.rds.gz`

Already completed:

- Downloaded DCT1 and DCT2 control RDS files.
- Verified gzip integrity.
- Recorded SHA-256 hashes in `SHA256SUMS`.

Required before analysis:

1. Load both objects in R.
2. Identify object class: Seurat, SingleCellExperiment, or another R object.
3. Extract expression matrix and metadata.
4. Confirm baseline/control condition is normal potassium or equivalent control.
5. Confirm DCT1 markers are high in DCT1 object:

   - Slc12a3
   - Pvalb
   - Trpm6
   - Klhl3

6. Confirm DCT2 / connecting-transition markers are relatively higher in DCT2
   where expected:

   - Scnn1a
   - Trpv5
   - Calb1, noting that Calb1 may be mixed and should not be forced.

7. Record the number of cells/nuclei, animals/replicates if available, and
   metadata fields.

Important statistical rule:

If biological replicate or animal IDs exist, DCT1 vs DCT2 differential analysis
should be pseudobulked by replicate. Do not treat 25,000 cells as 25,000
independent biological samples. If replicate metadata are missing or unusable,
cell-level marker analysis may be used only as a reference-prior construction
with an explicit pseudoreplication caveat.

Outputs:

- `gse228367_object_inventory.tsv`
- `gse228367_marker_qc.tsv`
- `gse228367_load_manifest.json`

### Step 1b: PXD001729 verification

Local inputs:

- `data/external/phosphoproteomics/PXD001729/41598_2015_BFsrep12829_MOESM1_ESM.pdf`
- `data/external/phosphoproteomics/PXD001729/41598_2015_BFsrep12829_MOESM2_ESM.xls`
- `data/external/phosphoproteomics/PXD001729/41598_2015_BFsrep12829_MOESM3_ESM.xls`
- `data/external/phosphoproteomics/PXD001729/41598_2015_BFsrep12829_MOESM4_ESM.xls`
- `data/external/phosphoproteomics/PXD001729/41598_2015_BFsrep12829_MOESM5_ESM.xls`
- `data/external/phosphoproteomics/PXD001729/41598_2015_BFsrep12829_MOESM6_ESM.xls`

Required before analysis:

1. Identify which sheet/table contains quantified phosphosites and dDAVP effect
   sizes.
2. Extract:

   - gene symbol
   - phosphosite residue
   - phosphopeptide / localization information
   - dDAVP vs vehicle log2 fold change or equivalent effect
   - p value or significance label
   - site localization confidence where available

3. Confirm whether these target genes/sites appear:

   - Slc12a3 / NCC
   - Stk39 / SPAK
   - Oxsr1 / OSR1
   - Wnk1
   - Wnk4
   - Klhl3
   - Cul3
   - Nedd4l
   - Sgk1

4. Record the count of total phosphosites and the count of dDAVP-altered sites.

Outputs:

- `pxd001729_table_inventory.tsv`
- `pxd001729_phosphosite_effects.tsv`
- `pxd001729_target_site_coverage.tsv`
- `pxd001729_load_manifest.json`

## 7. Step 2: build DCT1/DCT2 reference prior from GSE228367

### Step 2a: DCT1 vs DCT2 differential program

Inputs:

- Verified GSE228367 control DCT1 expression object.
- Verified GSE228367 control DCT2 expression object.

Preferred analysis:

1. Normalize expression consistently across DCT1 and DCT2 objects.
2. If replicate metadata exist, create pseudobulk profiles per replicate and
   subtype.
3. Fit DCT1 vs DCT2 differential expression.
4. Export continuous effect sizes and uncertainty:

   - log2 mean DCT1 / mean DCT2
   - standard error if model supports it
   - test statistic
   - p value
   - FDR
   - mean expression in DCT1
   - mean expression in DCT2
   - detection fraction in DCT1
   - detection fraction in DCT2

Fallback analysis:

If pseudobulk is impossible, compute a cell-level marker table and label it as a
reference-prior marker analysis rather than a biological replicate-level DE test.

Core gene-set definitions:

- DCT1 core: log2FC > 1, FDR < 0.05, detected in a meaningful fraction of DCT1
  nuclei.
- DCT2 core: log2FC < -1, FDR < 0.05, detected in a meaningful fraction of DCT2
  nuclei.
- DCT shared: expressed in both DCT1 and DCT2 without strong subtype bias.
- Unresolved: low expression, unmapped, or inconsistent.

Marker sanity checks:

- Slc12a3, Pvalb, Trpm6, and Klhl3 should be DCT1-high or at least compatible
  with DCT1 identity.
- Scnn1a and Trpv5 should lean DCT2/late-DCT where present.
- Calb1 should be treated empirically; do not force it into DCT1 or DCT2.

Outputs:

- `gse228367_dct1_vs_dct2_de.tsv`
- `gse228367_dct1_core_genes.tsv`
- `gse228367_dct2_core_genes.tsv`
- `gse228367_dct_shared_genes.tsv`
- `gse228367_dct_prior_summary.json`

### Step 2b: build DCT1 enrichment prior for OSD-462 proteins and phosphosites

Inputs:

- `gse228367_dct1_vs_dct2_de.tsv`
- OSD-462 protein table from existing anchor run.
- OSD-462 phosphosite table from existing anchor run.
- Project gene symbol / Ensembl mapping resources.

For each OSD-462 protein, assign:

- DCT1 enrichment score
- DCT1 core membership
- DCT2 core membership
- DCT shared membership
- DCT unresolved / unmapped
- mean DCT expression where available
- mapping confidence

For each OSD-462 phosphosite, assign:

- parent gene
- phosphosite identifier
- DCT1 enrichment score of parent gene
- DCT1 / DCT2 / shared class
- OSD-462 flight effect
- OSD-462 p value and FDR
- site role where known: NCC regulatory, NCC non-regulatory, SPAK/WNK pathway,
  other

Coverage checks:

- Fraction of OSD-462 proteins with DCT prior annotation.
- Fraction of OSD-462 phosphosites with DCT prior annotation.
- Coverage of Slc12a3, Stk39, Oxsr1, Wnk1, Wnk4, Klhl3, Cul3, Nedd4l, and Sgk1.

Outputs:

- `osd462_protein_dct1_prior.tsv`
- `osd462_phosphosite_dct1_prior.tsv`
- `osd462_dct1_prior_coverage.tsv`
- `dct1_enrichment_prior_v1.tsv`

## 8. Step 3: H2 main DCT1 enrichment test

Question:

Are OSD-462 spaceflight-suppressed phosphosites enriched in DCT1-high parent
genes?

Primary definition:

- Suppressed phosphosite = OSD-462 flight effect < 0 and phosphosite p < 0.05.
- Background = all quantified OSD-462 phosphosites with mapped parent genes and
  available DCT1 prior.

Primary tests:

1. Fisher exact test:

   - rows: suppressed vs not suppressed
   - columns: DCT1 core vs not DCT1 core

2. Mann-Whitney test:

   - compare continuous DCT1 enrichment scores for suppressed vs not suppressed
     phosphosites.

3. Matched-null sensitivity:

   - random phosphosite sets matched on parent protein abundance, phosphosite
     quantification coverage, and parent-gene expression where possible.

Specificity controls:

1. Repeat for DCT2 core membership.
2. Repeat for DCT shared membership.
3. Repeat excluding anchor genes:

   - Slc12a3
   - Stk39
   - Oxsr1
   - Wnk1
   - Wnk4

4. Repeat excluding all NCC phosphosites.
5. Repeat with stricter suppressed definition:

   - flight effect < 0 and FDR < 0.1

Falsification rule:

If DCT1 enrichment of suppressed phosphosites is not above null at FDR <= 0.1,
the DCT1-specificity claim is not supported. Do not tune thresholds to force a
positive result.

Claim labels:

| Result | Label |
| --- | --- |
| DCT1 enriched, DCT2 not enriched, sensitivity survives anchor-gene exclusion | DCT1-prioritized phospho suppression supported |
| DCT1 enriched only before anchor-gene exclusion | DCT1-prioritized anchor genes, broader DCT1 enrichment not established |
| DCT1 and DCT2 both enriched | DCT/DCT-shared phospho suppression, not subtype-specific |
| No enrichment | DCT1 specificity not supported |

Outputs:

- `h2_dct1_phosphosite_enrichment_summary.tsv`
- `h2_dct1_phosphosite_enrichment_background.tsv`
- `h2_dct1_sensitivity_summary.tsv`
- `h2_dct1_enrichment_verdict.json`

## 9. Step 4: PXD001729 dDAVP anti-alignment test

Question:

Does the OSD-462 spaceflight phosphosite direction oppose a DCT-lineage
vasopressin/dDAVP activation direction?

### Step 4a: build PXD001729 dDAVP direction vector

Inputs:

- `pxd001729_phosphosite_effects.tsv`

Procedure:

1. Standardize phosphosite identifiers to gene + residue.
2. Build dDAVP vs vehicle effect vector.
3. Flag dDAVP-altered sites where significance is available.
4. Restrict a targeted DCT transport/signaling subset:

   - Slc12a3
   - Stk39
   - Oxsr1
   - Wnk1
   - Wnk4
   - Klhl3
   - Cul3
   - Nedd4l
   - Sgk1
   - Kcnj10
   - Kcnj16
   - Trpm6
   - Pvalb
   - Calb1

Outputs:

- `pxd001729_ddavp_direction.tsv`
- `pxd001729_dct_target_direction.tsv`

### Step 4b: map PXD001729 to OSD-462

Inputs:

- `pxd001729_ddavp_direction.tsv`
- OSD-462 phosphosite effects.

Procedure:

1. Map shared phosphosites by gene + residue.
2. Record all unmapped target sites.
3. Check whether NCC/SPAK/WNK sites are shared.
4. Do not force site matches across ambiguous isoforms without recording
   uncertainty.

Outputs:

- `pxd001729_osd462_shared_phosphosites.tsv`
- `pxd001729_osd462_mapping_failures.tsv`

### Step 4c: anti-alignment statistic

Primary statistic:

- Cosine similarity between dDAVP effect vector and OSD-462 flight effect vector
  over shared phosphosites.

Uncertainty:

- Bootstrap shared phosphosites with replacement, 1000 or more draws.

Interpretation:

| Result | Label |
| --- | --- |
| cosine < 0, CI excludes 0 | anti-aligned with DCT activation direction |
| cosine < 0, CI crosses 0 | suggestive anti-alignment, underpowered |
| cosine approximately 0 | no directional relationship |
| cosine > 0 | not consistent with suppression of dDAVP activation direction |

Important caveat:

PXD001729 is an mpkDCT cell-line dDAVP stimulus. It is a biological plausibility
reference, not validation of spaceflight DCT behavior.

Outputs:

- `h2_pxd001729_ddavp_antialignment_summary.tsv`
- `h2_pxd001729_ddavp_antialignment_bootstrap.tsv`
- `h2_pxd001729_ddavp_antialignment_verdict.json`

## 10. Step 5: KLHL3 / CUL3 exploratory mechanism check

Question:

Can existing data distinguish KLHL3/CUL3-driven WNK turnover from ionic/osmotic
suppression of WNK activity?

Answer expected:

Probably not definitively. Public data lack ubiquitinomics. This is exploratory.

Checks:

1. In OSD-462, inspect:

   - Klhl3 RNA effect
   - Cul3 RNA effect
   - KLHL3 protein effect if quantified
   - CUL3 protein effect if quantified
   - KLHL3 phosphosites if quantified
   - WNK1 / WNK4 total protein effects

2. In PXD001729, inspect:

   - KLHL3 phosphosites
   - CUL3 phosphosites
   - WNK1 / WNK4 phosphosites

3. Specifically check whether KLHL3 Ser433 or any known WNK-regulatory KLHL3
   site is present.

Interpretation:

- If KLHL3 regulatory phosphorylation increases and WNK protein decreases, this
  is consistent with WNK turnover but not proof.
- If WNK protein is flat and WNK/SPAK/NCC phospho is down, ionic/osmotic or
  kinase-activity suppression remains plausible.
- If sites are absent, report not testable.

Outputs:

- `h2_klhl3_cul3_site_coverage.tsv`
- `h2_klhl3_cul3_effects.tsv`
- `h2_klhl3_cul3_interpretation.md`

Claim label:

- Exploratory candidate observation only.
- No mechanism claim without ubiquitinomics or perturbation.

## 11. Step 6: Bayesian mediation / causal-structure specification

Do not run this step until the H2 DCT1 prior is complete.

Question:

Are the OSD-462 data more consistent with partial mediation through an
endothelial/stromal remodeling state, or with parallel flight responses?

Data:

- n=20 matched OSD-462 flight and ground animals.
- X = flight status.
- Y = per-animal NCC regulatory phosphorylation activity score.
- Candidate mediator M = one per model.

Candidate mediators:

1. Endothelial compartment score.
2. Stromal/fibroblast compartment score.
3. DCT identity score as a composition-control mediator.
4. Matrix/endothelial composite score.
5. DCT1-prior-weighted remodeling score if Step 2 supports construction.

Model form:

Mediator model:

`M_i = alpha_m + a * X_i + error_m`

Outcome model:

`Y_i = alpha_y + c_prime * X_i + b * M_i + error_y`

Indirect effect:

`a * b`

Direct effect:

`c_prime`

Total effect:

`c_prime + a * b`

Implementation:

- Preferred: Bayesian SEM in `brms` or `blavaan`.
- Acceptable fallback: Bayesian linear models with posterior draws of a and b,
  then draw-wise multiplication for the indirect effect.
- Standardize X, M, and Y where appropriate so path coefficients are comparable.
- Priors:

  - Normal(0, 1) for standardized path coefficients.
  - Half-Normal(0, 0.5) for residual SDs.

- 4 chains.
- 4000 iterations.
- 2000 warmup.
- Rhat < 1.01.

Report:

- posterior median of a, b, c_prime, and a*b
- 95% credible interval
- P(indirect effect < 0) or the pre-registered directional probability
- posterior predictive checks
- sensitivity to mediator choice

Decision table:

| Posterior result | Label | Do not say |
| --- | --- | --- |
| indirect effect < 0 and CI excludes 0 | consistent with mediation | demonstrates mediation |
| indirect effect < 0 and CI includes 0 | suggestive of mediation, underpowered | supports mediation |
| indirect effect approximately 0 | parallel responses cannot be ruled out | mediation absent |
| indirect effect > 0 | inconsistent with mediation hypothesis | hide or soften |

Power / future-n deliverable:

Use posterior simulations to estimate the sample size required to detect the
observed indirect effect with 80% power under a future matched design.

Outputs:

- `h3_mediation_input_scores.tsv`
- `h3_mediation_model_summary.tsv`
- `h3_mediation_posterior_draws.tsv.gz`
- `h3_mediation_power_simulation.tsv`
- `h3_mediation_verdict.json`

## 12. Step 7: optional spatial reference projection

Do this only after H2 and H3 core analyses are complete.

Sources:

- GSE269719: Xenium mouse kidney IRI/repair spatial reference.
- GSE269622: Visium mouse kidney IRI/repair spatial reference.

Local status:

- Not downloaded in the core setup because raw archives are large and the
  analysis is optional.

Purpose:

Use an IRI spatial atlas as a reference to ask which injury/repair spatial niche
the spaceflight bulk RNA vector most resembles.

Procedure:

1. Download processed annotations only if available. Avoid raw FASTQ unless
   absolutely required.
2. Extract cell types / niches:

   - DCT
   - endothelial
   - fibroblast / stromal
   - macrophage / immune
   - injured tubular
   - peri-tubular interstitial

3. Compute pathway or gene-set signatures for the existing 11-pathway panel.
4. Compare RRRM-2 ISS-T and OSD-462 RNA flight vectors against spatial niche
   signatures by cosine similarity.
5. Identify best matching timepoint and niche.

Claim label:

- Spatial reference projection.
- Not spaceflight spatial transcriptomics.
- Not evidence that a spatial niche exists in RR-10 unless actual sections are
  profiled.

Outputs:

- `spatial_reference_dataset_inventory.tsv`
- `spatial_reference_pathway_signatures.tsv`
- `spatial_reference_projection_summary.tsv`
- `spatial_reference_projection_verdict.json`

## 13. Integrated v11 claim ladder

| Layer | Dataset | Claim | Confidence label |
| --- | --- | --- | --- |
| RNA context recurrence | RRRM-2, OSD-513, OSD-462 | matrix-high / DCT-low state recurs | established |
| Protein abundance | OSD-462 | matrix and DCT RNA signals do not propagate as protein abundance | established negative |
| Phospho-regulatory suppression | OSD-462 | NCC regulatory cluster and SPAK/WNK activity are down | established / reproduced from Siew et al. |
| DCT1 subtype prior | GSE228367 -> OSD-462 | suppressed phosphosites are DCT1-prioritized if H2 passes | new falsifiable reference-prior result |
| DCT phospho reference | PXD001729 -> OSD-462 | spaceflight direction anti-aligns with dDAVP activation if Step 4 passes | cell-line plausibility reference |
| Composition / compartment | marker panels + DCT1 prior | DCT-low includes interstitial expansion component | bounded claim |
| Mediation structure | OSD-462 matched n=20 | consistent with partial mediation or parallel responses, with wide CI | hypothesis specified |
| Spatial context | GSE269719/GSE269622 optional | resemblance to IRI spatial niche | optional reference only |
| Physiology | absent | animal-matched functional outcome not available | future work |

## 14. Manuscript restructure

### Proposed three-hypothesis Results layout

1. H1: recurrent spaceflight kidney remodeling context.
2. H1 boundary: cross-omic decoupling and protein-abundance null.
3. H2: DCT1/DCT2 reference prior from GSE228367.
4. H2: DCT1 enrichment of OSD-462 flight-suppressed phosphosites.
5. H2: PXD001729 dDAVP anti-alignment as DCT-lineage phospho plausibility.
6. H2 boundary: DCT specificity inferred, not directly measured.
7. H3: matched OSD-462 mediation / parallel-response model.
8. H3 boundary: underpowered, cross-sectional, composition-confounded.
9. Integrated evidence ladder and named future experiments.

### Title guidance

Better:

- Spaceflight kidney remodeling is DCT1-prioritized and cross-omically decoupled:
  a phosphoproteome-informed causal hypothesis

Safer:

- DCT1-anchored phospho-regulatory suppression in spaceflight kidney remodeling

Avoid:

- "Mediation-tested mechanism"
- "Causal mechanism"
- "DCT-specific spaceflight phosphoproteomics"

### One-sentence contribution

Spaceflight kidney remodeling is not simply a bulk RNA injury signature: a
recurrent matrix/endothelial-high and DCT-low RNA context surrounds an
established NCC/SPAK phosphorylation lesion that does not propagate as protein
abundance, maps through external references to a DCT1/NCC-centered
phospho-regulatory hypothesis, and yields an explicit, underpowered mediation
model distinguishing remodeling-linked from parallel transporter deactivation.

## 15. Stop rules

Primary stop rule:

If H2 Step 3 fails to show DCT1 enrichment above null at FDR <= 0.1, do not
claim DCT1-prioritized phospho suppression. Report the null and either:

1. submit the paper as v10 plus a DCT-reference null result, or
2. keep H3 as future work rather than running a mediation model built on a failed
   DCT1 prior.

Secondary stop rules:

- If PXD001729 shares too few phosphosites with OSD-462 for a stable cosine,
  replace the anti-alignment claim with a target-site coverage table.
- If mediation credible intervals are wide, use "suggestive / underpowered" and
  lead with the power calculation.
- If optional spatial data require raw reprocessing beyond scope, skip spatial.
- Do not add new omics/network layers unless they directly address H2 or H3.

## 16. Future experiments to name explicitly

### Priority 1: DCT-enriched spaceflight phosphoproteomics

Design:

- NCC-Cre-INTACT or equivalent DCT-enriched isolation from spaceflight or a
  carefully chosen analog.
- FANS-sort DCT nuclei/cells.
- Split material for snRNA-seq and phosphoproteomics where feasible.

Readout:

- DCT-specific phosphoproteomic map of flight vs ground.

Why it matters:

- Directly resolves whether NCC/SPAK phospho-suppression is DCT cell-autonomous
  or a whole-kidney bulk effect shaped by interstitial expansion.

### Priority 2: WNK-SPAK-NCC functional perturbation / rescue

Design A:

- mpkDCT cells or kidney organoids under simulated microgravity.
- Ubiquitin remnant or K48 ubiquitin enrichment.
- Test whether WNK1/WNK4 ubiquitination increases.

Design B:

- Simulated microgravity with chloride or ionic-state clamping.
- Test whether NCC phosphorylation is rescued.

Why it matters:

- Distinguishes KLHL3/WNK turnover from ionic/osmotic kinase suppression.

### Priority 3: spaceflight kidney spatial transcriptomics

Design:

- Visium or Xenium on RR-10 or equivalent flight kidney sections.
- Pair with pNCC/pSPAK immunofluorescence on adjacent sections.

Readout:

- Geographic relation between matrix/endothelial/stromal expansion and DCT/NCC
  segments.

Why it matters:

- Resolves DCT suppression vs segment loss vs interstitial dilution more directly
  than bulk RNA.

## 17. Proposed timeline

| Week | Step | Deliverable |
| --- | --- | --- |
| 1 | pre-registration and baseline lock | plan checksum, baseline manifest |
| 1-2 | load and verify GSE228367 | object inventory, marker QC |
| 1-2 | load and verify PXD001729 | phosphosite table inventory |
| 2-3 | DCT1 vs DCT2 reference prior | DCT1/DCT2 DE and core sets |
| 3-4 | map DCT prior to OSD-462 | protein/phosphosite prior tables |
| 4-5 | H2 enrichment tests | DCT1 enrichment verdict |
| 5-6 | PXD001729 anti-alignment | dDAVP anti-alignment verdict |
| 6 | KLHL3/CUL3 check | exploratory mechanism table |
| 7-8 | H3 mediation model | posterior summary and power estimate |
| 9 | apply claim labels | H2/H3 integrated verdict |
| 9-10 | optional spatial projection | spatial reference supplement if feasible |
| 10-14 | manuscript rewrite | v11 draft |

## 18. Final honesty envelope

This study cannot become a full mechanistic paper from public data alone.

It can become a sharper, more novel paper that:

- preserves the v10 cross-omic decoupling result,
- adds DCT1/DCT2 subtype specificity at reference-prior resolution,
- uses DCT-lineage phosphoproteomics to test biological plausibility,
- formally specifies mediation vs parallel-response structure,
- reports nulls and wide intervals openly,
- and names the exact experiment required for causal proof.

That is meaningfully more than v10 while staying inside the evidence envelope.
