# Deferred: second independent atlas validation of the podocyte marker set (G1 item 9)

**Status:** DEFERRED by owner decision, 2026-10-06. Not started; nothing below has been computed.
**Recorded in:** `config/clinical_axes_podocyte_sensitivities.yaml` (`items.9_second_atlas.status: DEFERRED`).
**Plan source:** work package WP-A9 of the 2026-10-05 implementation plan, with facts verified on
2026-10-06 during the Mouse Kidney Atlas rebuild.

## 1. What this item would establish

The 157-gene `podocyte__high_specificity` set was defined from one reference, the integrated Mouse
Kidney Atlas (MKA). "Independently defined" currently means flight-label-blind, not reference-robust.
Item 8 (leave-one-atlas-source-study-out) probes robustness *within* MKA. Item 9 asks whether a
podocyte set built the same way from a **different atlas, made by a different lab with different
animals**, reproduces both the gene membership and the cross-mission flight shift.

## 2. Why it was deferred

- K0 (`docs/` revision conclusion, `config/revisions/clinical-axes.yaml`) already demoted the podocyte
  program from headline: the podocyte flight coefficient attenuates 0.198 → 0.112 (score units; HC3 CI
  crosses zero) after adjusting for structural drift. The story is now "shared structural drift with a
  podocyte-leaning tail", so reference robustness of the podocyte *label* is lower-stakes.
- Effort: about 2–3 days of coding plus label curation; the candidate atlas's podocyte labels need
  remapping (section 4).

**Un-defer if any of these becomes true:** "podocyte" stays in the paper's title or main claims; a
reviewer questions single-reference marker definition; item 8 shows the set is fragile to dropping a
source study.

## 3. Reference choice: Tabula Muris Senis (TMS) kidney

Independence from MKA was checked against MKA's eight `Origin` source studies (github.com/nrclaudio/MKA):
Wu19 GSE119531, Miao21 GSE157079, Park18 GSE107585, Kirita20 GSE139107, Dumas20 E-MTAB-8145,
Conway20 GSE140023, Hinze21 GSE145690, Janosevic21 GSE151658. Neither Tabula Muris nor TMS is among them.
TMS is a different lab, 10x v2 plus Smart-seq2, C57BL/6JN, both sexes, 1–30 months, with mouse-level
donors, so replication can use the animal as the unit (MKA only supports study-level units).

### Files already downloaded (gitignored, `data/external/single_cell_atlases/tms_kidney_cellxgene/`)

Source: CELLxGENE collection `0b9d8a04-bb9d-44da-aa27-705bb65b54eb` (TMS, Nature 2020,
doi:10.1038/s41586-020-2496-1).

| File | Size | SHA-256 | Content |
|---|---|---|---|
| `cfd31d35-ceb2-4990-b0b5-3434ecbb7967.h5ad` | 343,537,513 B | `02d37b0cc51c2aaaf2cc27c3a2ddadd4d97bafebe247c37d64a363226ef5c6b4` | 10x 3′ v2; 21,647 cells × 17,943 genes (Ensembl IDs, `feature_name` in var); `raw.X` integer counts; 16 donors (11 M / 5 F) aged 1/3/18/21/24/30 m |
| `7769f88b-96e2-46e8-a2b2-2b7282804da7.h5ad` | 33,234,927 B | `01317ea79bc33dfb927bca35930f4aec8d4138ce4aa0a60f19609a1c2158cc3e` | Smart-seq2; 1,833 cells × 21,025 genes; 14 donors; 10 cell types; **no podocytes** |

obs columns (10x): `cell_type` (23 ontology labels), `free_annotation` (19), `donor_id`, `sex`, `age`,
`subtissue`, `method`.

Alternative sources (equivalent content, different packaging) on public S3:
`https://czb-tabula-muris-senis.s3.amazonaws.com/Data-objects/tabula-muris-senis-droplet-processed-official-annotations-Kidney.h5ad`
(776,671,418 B); FACS kidney (63,459,198 B); all-tissue raw droplet object (4,056,529,818 B).

### Known pitfall: podocyte labelling

In the CELLxGENE 10x object the ontology `cell_type` labels **447 of 486** podocytes as
`plasmatocyte`. The usable label is `free_annotation == "Epcam    podocyte"` (486 cells; 124 from
female donors). The compartment mapping must therefore use `free_annotation` for podocytes, never
`cell_type`, and the alias audit must show this explicitly.

The repo's config already anticipates a female-only TMS fallback,
`tms_kidney_female_ALLDATASETS_counts_innerGenes.h5ad`
(`config/dct_subtype_reference_freeze_v1.yaml:61-64`, `required: false`). That file is not a published
object: it is the female cells of both TMS kidney datasets, inner-joined on genes. Restricting to
females leaves only 124 podocytes across at most 5 donors, so item 9 should use **both sexes** with sex
recorded, not the female-only object.

## 4. Implementation plan (preserved)

**New files:**
- `scripts/subtype_reference/05_tms_kidney_pseudobulk.py`
- `scripts/clinical_axes/run_second_atlas_validation.py`
- `tests/test_second_atlas_validation.py`
- registry stage `second-atlas-validation` (`optional: true`, stage-level `requires` on the TMS file and the frozen tiers)

**Steps:**
1. **Counts.** Use `raw.X` and assert the counts are integer-like. If they are not, fall back to the
   all-tissue raw droplet object subset to `tissue == Kidney`.
2. **Compartments.** Map labels to compartments with the frozen alias regexes
   (`dct_subtype_reference_freeze_v1.yaml`, `whole_kidney_compartment_aliases`), using
   `free_annotation` for podocytes (section 3). Write an audit table of every source label →
   compartment and every unmapped label.
3. **Pseudobulk.** The replication unit is `donor_id`; a pseudobulk needs at least 25 cells. Reuse the
   CPM, detection and mean/median logic of `03_atlas_pseudobulk.py`; detection fraction is computed over
   donors instead of source studies.
4. **Evaluability gate (preregistered; fixed before any scoring).** The podocyte compartment needs
   at least 3 donors with at least 25 podocytes each. If not, pool by age group and require at least 3
   units of at least 25 cells. If that also fails, item 9 is NOT_EVALUABLE.
   Compartments absent from TMS (expected: DCT2/CNT) are dropped from the max-other denominator, and
   this is disclosed.
5. **Tiers.** Call `build_marker_tiers` with a copied v13 config whose `set_test.primary_family` is
   restricted to the TMS-available sets (the function raises unless the config lists exactly its set
   names). Thresholds are unchanged: mean CPM ≥ 1, detection ≥ 0.75, ≥ 4-fold above every other
   modeled compartment.

**Validation metrics:**
- **(a) Definitional overlap.** Jaccard and hypergeometric overlap between the MKA and TMS podocyte HS
  sets over the gene universe shared by both atlases. Replication rate: the fraction of the 157 MKA
  genes that are at least 2-fold podocyte-enriched in TMS.
- **(b) Flight scoring.** Score the five terminal missions with three sets: TMS podocyte HS, the
  consensus set (MKA ∩ TMS), and the TMS-only set. Report the pooled g (REML/mHK) and K0 P1 for each.
- **(c) Family rank.** Rank of the podocyte set within the TMS-defined family; descriptive, seed
  offset +900.

**Gates (already in the preregistration YAML):**
- The TMS podocyte HS estimate has the sign of the reference with retention ≥ 0.50.
- The consensus set estimate is > 0.

Run under the baseline and runner-up gene maps (reconstruction sensitivity), as for every G1 item.

**Fallback references** if the TMS podocyte compartment is NOT_EVALUABLE: Ransick 2019 (GSE129798) or
Denisenko 2020 (GSE141115) from GEO. The MKA-extended atlas (303,791 cells) is not suitable: it adds
injury and PKD populations, and its labels were transferred from the baseline atlas.

**Estimated cost:** 2–3 days of coding; about 15 minutes of compute; downloads already done (377 MB).

## 5. Risks to carry forward

- Podocyte capture in droplet data is low (486 cells), so the evaluability gate may fail.
- TMS has no DCT2/CNT compartment, which makes the max-other denominators differ from MKA's.
- TMS animals are older on average (up to 30 months); age-related podocyte loss could shift detection.
- Platform differences (10x v2 versus the MKA studies' chemistries) affect detection rates; thresholds
  are deliberately not re-tuned.
