# Flight-blind DCT subtype reference rebuild

## Purpose

This workflow defines the gene sets consumed by
`config/dct_asdn_phospho_reanalysis.yaml` without reading OSD-462 flight or
phosphoproteomic results. The production handoff is:

`data/processed/v13_dct_reference/frozen_gene_sets.tsv`

The frozen marker rules are in
`config/dct_subtype_reference_freeze_v1.yaml`. The ASDN genes in that file are
required to match the production analysis config exactly.

## Reference roles

1. **GSE228367 discovery:** the author-separated normal-potassium DCT1 and
   DCT2 objects supply cell barcodes, subtype membership, and `Rep` metadata
   only. Their `RNA@counts` slots are corrected and approximately 99.5%
   fractional, so they are rejected for edgeR. Every subtype barcode is
   matched to the official NK1–NK3 filtered 10x H5 matrices in
   `GSE228367_RAW.tar`; integer UMI counts are then summed by biological
   replicate. A paired edgeR quasi-likelihood model
   (`~ replicate + subtype`) estimates DCT2 minus DCT1. Individual nuclei are
   never treated as independent replicates.
2. **GSE150338 segment validation:** official microdissected DCT, CNT, and CCD
   replicate TPM profiles are used only to test whether a discovery candidate
   is retained from DCT into CNT. Technical aliquots are averaged within
   mouse. DCT is not relabeled as DCT1 or DCT2.
3. **Fine-subtype validation:** the official author-integrated GSE150338
   Seurat object supplies raw RNA counts, three `Batch_10X` replicates, and
   integrated clusters. The DCT1/DCT2/CNT cluster map is the pre-existing
   marker-checked mapping in `scripts/build_hybrid_reference.R`; it is frozen
   in the config and re-audited with canonical markers. It is not represented
   as an author-provided cell-type column.
   A separate `DCT2_CNT_external_validation` signature is derived entirely
   from these GSE150338 fine-subtype and segment data, followed by the frozen
   combined distal-specificity/breadth exclusion. It never consults GSE228367
   membership and is marked as an external-validation set rather than part of
   the primary hypothesis family.
4. **Whole-kidney specificity:** the integrated 141,401-cell Mouse Kidney
   Atlas supplies a true multi-compartment breadth flag, a separate combined
   distal-specificity/breadth flag, and all major-compartment comparator sets.
   Full `raw.X` (16,711 genes) is used instead of the
   3,000-variable-gene `counts` layer. `Celltype_finest` is primary and
   `Celltype` is the fallback label. Its `Origin` field identifies eight
   source studies—not donors or animals—so source-study-stratified detection
   is a reproducibility check rather than animal-level inference.

## Frozen rules

- `DCT1_core`: DCT1-enriched in the paired GSE228367 discovery, directionally
  replicated in an independent fine-subtype reference, and passing the combined
  distal-specificity/breadth exclusion.
- `strict_DCT2_peaked`: DCT2-enriched over DCT1 in discovery, directionally
  confirmed in the independent fine-subtype object, replicate-consistent,
  adequately detected, and passing the combined distal-specificity/breadth
  exclusion. DCT2-versus-CNT effects are retained as characterization rather
  than added as a post-freeze gate.
- `DCT2_CNT_transition`: DCT2-enriched over DCT1 in GSE228367, retained in CNT
  under the predeclared GSE150338 non-inferiority margin, and passing the
  combined distal-specificity/breadth exclusion.
- `DCT2_CNT_external_validation`: independently defined from GSE150338 genes
  that are DCT2-enriched over DCT1, retained in CNT in both the fine-subtype
  and microdissected-segment references, and passing the combined
  distal-specificity/breadth exclusion. This is a validation-only signature and
  is not in the primary multiplicity family.
- `ASDN`: literature-frozen categories copied exactly from the production
  analysis config.
- Major comparator sets are frozen for proximal tubule, TAL, DCT1, a curated
  late-DCT/CNT context set, principal and intercalated cells, podocytes,
  endothelium, fibroblasts, and immune cells.

The broad-expression table distinguishes target expression, true
multi-compartment breadth, and lack of distal specificity. `target_expressed`
requires maximum distal-target CPM at least 1. `broadly_expressed` is true only
when a target-expressed gene is detected above the configured floor in at least
four unrelated compartments. This is the flag used by the phosphoproteomic
broad-expression sensitivity.

The separate `non_distal_specific_or_broad` compatibility flag preserves the
earlier combined rule: a target-expressed gene is flagged when it is truly broad
or when its maximum unrelated-compartment expression is at least half its
maximum distal-target expression. Data-driven DCT signatures continue to use
this stricter specificity exclusion, so separating the meanings does not
silently broaden the frozen signatures. A gene absent from the distal target
but specific to one unrelated compartment is false for both flags.

The strict DCT1/DCT2 discovery thresholds are absolute log2 fold change at
least 1, BH FDR at most 0.05, target-cell detection at least 20%, and the same
direction in all three biological-replicate pairs. The transition discovery
rule follows the production plan: DCT2-minus-DCT1 log2 fold change at least
0.5, DCT2 detection at least 20%, and the same direction in all three pairs;
it does not add an unplanned discovery-FDR gate. External DCT-to-CNT retention
requires both an observed log2 fold change of at least -0.5 and a one-sided
lower confidence bound above -1 with BH FDR at most 0.05.

## Commands

```bash
venv/bin/python scripts/subtype_reference/01_inventory_inputs.py \
  --config config/dct_subtype_reference_freeze_v1.yaml \
  --output-dir data/processed/v13_dct_reference/inventory

Rscript scripts/subtype_reference/02_gse228367_pseudobulk.R \
  --config config/dct_subtype_reference_freeze_v1.yaml \
  --output-dir data/processed/v13_dct_reference/gse228367

Rscript scripts/subtype_reference/02b_gse150338_segment_validation.R \
  --config config/dct_subtype_reference_freeze_v1.yaml \
  --output-dir data/processed/v13_dct_reference/gse150338

Rscript scripts/subtype_reference/02c_gse150338_fine_subtype_pseudobulk.R \
  --config config/dct_subtype_reference_freeze_v1.yaml \
  --output-dir data/processed/v13_dct_reference/gse150338_fine

venv/bin/python scripts/subtype_reference/03_atlas_pseudobulk.py \
  --config config/dct_subtype_reference_freeze_v1.yaml \
  --input data/external/single_cell_atlases/mouse_kidney_atlas/mka.h5ad \
  --output-dir data/processed/v13_dct_reference/mouse_kidney_atlas \
  --reference-label mouse_kidney_atlas \
  --source-study-column Origin \
  --cell-type-column Celltype_finest \
  --fallback-cell-type-column Celltype

venv/bin/python scripts/subtype_reference/04_build_frozen_signatures.py \
  --config config/dct_subtype_reference_freeze_v1.yaml \
  --output-dir data/processed/v13_dct_reference \
  --discovery data/processed/v13_dct_reference/gse228367/gse228367_dct2_vs_dct1_paired_pseudobulk.tsv \
  --validation data/processed/v13_dct_reference/gse150338_fine/gse150338_fine_subtype_validation.tsv \
  --segment-validation data/processed/v13_dct_reference/gse150338/gse150338_segment_transition_validation.tsv \
  --atlas data/processed/v13_dct_reference/mouse_kidney_atlas/atlas_compartment_expression.tsv.gz
```

## Current input status and limitations

- Both author-separated GSE228367 control objects are present and their
  SHA-256 hashes match the frozen manifest. Their RDS count slots are not used:
  25,155 of 25,155 subtype-member barcodes map exactly to the official filtered
  10x matrices, and all six reconstructed pseudobulks contain integer counts.
- Both official GSE150338 references are present and checksum-verified: the
  microdissected segment TPM matrix and the author-integrated Seurat RData.
  Fine cell types are a frozen marker-based mapping of integrated clusters,
  not an embedded author cell-type field.
- The Mouse Kidney Atlas is present and checksum-verified (MD5
  `6a2e53086f100acffe99c054d3e8b7e9`). Full `raw.X` contains 16,711 genes;
  132,275 of 141,401 cells (93.55%) map to a frozen compartment, yielding 65
  retained source-study-by-compartment pseudobulks. Every required compartment
  is represented in 4–8 source studies. The atlas supports cross-compartment
  specificity, not animal-level p-values.

The workflow itself performs no downloads. Any future missing input generates
an explicit inventory and signature-gate record.

## Frozen handoff

The completed flight-blind build contains 6 `DCT1_core`, 12
`strict_DCT2_peaked`, 27 `DCT2_CNT_transition`, 17
`DCT2_CNT_external_validation`, and 29 literature-frozen `ASDN` genes. The
independent external set has `analysis_role=external_validation`; it is
available for a directional-concordance gate but is not part of the two-set
primary multiplicity family. These counts were frozen without consulting
phosphoproteome observability or flight effects.
