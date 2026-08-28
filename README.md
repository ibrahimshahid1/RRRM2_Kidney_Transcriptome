# RRRM2 Kidney Transcriptome and OSD-462 Phosphoproteomic Reanalysis

This repository contains the code, frozen configurations, tests, provenance
audits, exact permutation outputs, figures, and manuscript for a reanalysis of
mouse-kidney spaceflight omics.

The current analysis is v13. It supersedes the v11/v12 claim that
flight-associated phosphorylation suppression extends beyond the canonical
NCC/DCT1 axis into DCT2/CNT or the aldosterone-sensitive distal nephron
(ASDN).

## Current scientific conclusion

The frozen v13 question was:

> After excluding the established NCC/DCT1 axis, is flight-associated
> phosphosite suppression disproportionately concentrated in independently
> defined DCT2/CNT-transition or ASDN programs?

The definitive exact result is:

- **DCT2/CNT transition: non-evaluable.** Five frozen genes were observable
  after the primary exclusions; the prespecified minimum was eight.
- **ASDN: failed.** Sixteen genes produced a positive channel-centered
  competitive statistic (`0.718`; conditional exact `p = maxT = 0.0291`), but
  the result failed the prespecified selectivity criterion. A podocyte
  annotation had a larger statistic (`1.112`; `p = 2.36e-4`; BH
  `q = 0.00213`), and the observability-matched ASDN test was not significant
  (`p = 0.247`). No direct ASDN-versus-podocyte contrast was tested.
- **Overall statistical tier: `neither`.** This preserves DCT2/CNT as
  under-covered rather than calling it null, while recording that ASDN failed
  a frozen claim gate.
- **Publication promotion: blocked.** Condition is perfectly aligned with
  reporter-tag blocks in both plexes, and the MS workbook contains no
  isolated, assay-qualified canonical NCC/SPAK phosphoform.

The defensible conclusion is therefore a technical and null boundary:

> OSD-462 does not support a selective extension of phosphosite suppression
> into DCT2/CNT or ASDN. Its reporter allocation and phosphoform provenance
> also prevent attribution to flight or isolated canonical NCC/SPAK
> regulation.

This is not evidence that DCT2/CNT biology is unchanged. The frozen DCT2/CNT
test lacked sufficient phosphoproteomic coverage.

## Why the previous headline was retired

Stage 0 resolved three blocking provenance issues:

1. The detailed protocol specifies **TMTpro** labeling, a `+304.207`-Da tag
   mass, and TMTpro isotope-impurity correction. A legacy repository field
   says iTRAQ. The manuscript now reports “TMTpro isobaric-tag
   proteomics/phosphoproteomics” and discloses the inconsistency.
2. The two 15-plexes each place baseline samples in reporter positions 1--5,
   flight in 6--10, and ground control in 11--15. There is no cross-plex label
   swap. A flight coefficient is therefore inseparable from a systematic
   reporter-position effect.
3. The NCC T53-indexed feature is a T53/Y65 co-modified peptide, and the SPAK
   S383-indexed feature is an S382/S383 co-modified peptide. NCC positions 65
   and 68 map to tyrosines. Zero isolated NCC T53/T58/S71 or SPAK T243/S383
   phosphoforms passed the frozen qualification rule.

Prior antibody work on RR-10 remains relevant literature context, but it uses
the same biological experiment and is not an independent validation of this
reanalysis.

## What changed statistically

The old nominal-`p < 0.05` phosphosite contingency endpoint was replaced.

- The primary unit is one **parent gene**, not a phosphosite row.
- The primary site universe is label-blind: localization score at least 19,
  singly modified analytical phosphoforms, no composite rows, and fixed
  completeness requirements.
- A site effect is the equal-weight mean of the two within-plex
  flight-minus-ground log2 contrasts.
- The primary gene score is
  `median(-site effect)`, with positive values representing lower relative
  phosphopeptide signal in flight-labeled channels.
- Every gene is standardized against its own balanced-label null.
- A set is tested competitively as mean gene-specific `Z` in the set minus
  mean `Z` in the eligible background.
- The complete site, gene, and set calculation is rerun for all
  `choose(10, 5)^2 = 63,504` balanced assignments.
- DCT2/CNT transition and ASDN form one maxT family. Major kidney-compartment
  comparators form one Benjamini--Hochberg family.
- The primary canonical exclusion removes `Slc12a3`, `Stk39`, `Oxsr1`, and
  `Wnk4`; the strict exclusion additionally removes `Wnk1`.

This conditional permutation calibration preserves site dependence, site
count, missingness, and measured signal structure. It does not remove the
observed condition--reporter-tag alias.

## Flight-blind expression references

The subtype and compartment sets were built without reading OSD-462 flight
effects or phosphoproteomic observability.

- **GSE228367:** raw integer counts reconstructed from official 10x matrices;
  DCT1 and DCT2 counts summed by three biological replicates; paired edgeR
  pseudobulk model.
- **GSE150338:** independent fine-subtype direction check and microdissected
  DCT-to-CNT retention check.
- **Mouse Kidney Atlas:** whole-kidney specificity, true multi-compartment
  breadth, and prespecified comparator sets.

The frozen build contains:

| Set | Frozen genes | Observable after primary exclusions |
|---|---:|---:|
| DCT2/CNT transition | 27 | 5 |
| DCT2/CNT external validation | 17 | 5 |
| ASDN | 29 | 16 |
| strict DCT2 | 12 | not a primary set |
| DCT1 core | 6 | not a primary set |

The frozen membership file has SHA-256:
`ed64f0e20361d21ebc4d2483fce383f5f82e94bb7ec43b66156826c8eb86b945`.

## Key exact results

Primary canonical-axis exclusion, median parent-gene score:

| Program | Observable genes | Statistic | Exact p | Adjusted p |
|---|---:|---:|---:|---:|
| DCT2/CNT transition | 5 | not evaluated | not evaluated | not evaluated |
| ASDN | 16 | 0.718 | 0.0291 | maxT 0.0291 |
| Podocyte | 44 | 1.112 | 0.000236 | BH q 0.00213 |
| Fibroblast | 34 | 0.707 | 0.00375 | BH q 0.0169 |
| Thick ascending limb | 16 | 0.588 | 0.0462 | BH q 0.134 |
| DCT1 | 9 | 0.187 | 0.301 | BH q 0.452 |

Important ASDN sensitivities:

| Profile | Genes | Statistic | Exact p |
|---|---:|---:|---:|
| Primary centered scaled signal | 16 | 0.718 | 0.0291 |
| Strict exclusion adding `Wnk1` | 15 | 0.632 | 0.0575 |
| Deduplicated multi-modified/composite | 16 | 0.649 | 0.0478 |
| Exclude true multi-compartment-broad genes | 12 | 0.538 | 0.0696 |
| Official scaled signal, uncentered | 16 | 0.202 | 0.216 |
| Summed S/N, centered | 16 | 0.789 | 0.0247 |
| Parent-protein subtraction | 11 | 1.337 | 0.000772 |

The parent-protein profile loses 930 genes and still has a slightly larger
podocyte statistic (`1.342`). Scaled and summed S/N are algebraically related
in the source workbook, so their agreement is not independent validation.

## Repository map

### Authoritative specifications and decisions

- [`config/dct_asdn_phospho_reanalysis.yaml`](config/dct_asdn_phospho_reanalysis.yaml)
- [`config/dct_subtype_reference_freeze_v1.yaml`](config/dct_subtype_reference_freeze_v1.yaml)
- [`docs/DCT_ASDN_PHOSPHO_REANALYSIS_PLAN_2026-07-28.md`](docs/DCT_ASDN_PHOSPHO_REANALYSIS_PLAN_2026-07-28.md)
- [`docs/V13_FINAL_REANALYSIS_AND_MANUSCRIPT_DECISION_2026-07-29.md`](docs/V13_FINAL_REANALYSIS_AND_MANUSCRIPT_DECISION_2026-07-29.md)
- [`docs/MANUSCRIPT_V1_V12_CROSS_VERSION_AUDIT_2026-07-29.md`](docs/MANUSCRIPT_V1_V12_CROSS_VERSION_AUDIT_2026-07-29.md)

### Implementation

- `src/multiomics/osd462_stage0.py`
- `src/subtype_reference/`
- `src/v13/`
- `scripts/osd462/08_stage0_provenance_audit.py`
- `scripts/osd462/09_stage0_manuscript_reporting.py`
- `scripts/subtype_reference/`
- `scripts/v13/`

### Tests

- `tests/test_osd462_stage0.py`
- `tests/test_osd462_stage0_reporting.py`
- `tests/test_subtype_reference_builder.py`
- `tests/test_continuous_phospho_inference.py`
- `tests/test_v13_corrected_engine_contracts.py`
- `tests/test_v13_reporting.py`
- `tests/test_v13_provenance.py`

### Manuscript

- Source:
  [`latex_paper/manuscript_v13.tex`](latex_paper/manuscript_v13.tex)
- Rendered PDF:
  [`latex_paper/manuscript_v13.pdf`](latex_paper/manuscript_v13.pdf)

The v13 manuscript was rebuilt from scratch as an eight-page
provenance/statistical-boundary report. It has no table of contents, no
executive summary, two main figures, compact captions, and an abstract that
fits on the first page.

## Reproduction

The complete analysis requires the external data bundle under `data/external/`.
Large external inputs and generated matrices may be excluded from Git.

### Stage 0 provenance

```bash
venv/bin/python scripts/osd462/08_stage0_provenance_audit.py
venv/bin/python scripts/osd462/09_stage0_manuscript_reporting.py
```

### Flight-blind reference build

See
[`docs/DCT_SUBTYPE_REFERENCE_REBUILD_2026-07-28.md`](docs/DCT_SUBTYPE_REFERENCE_REBUILD_2026-07-28.md)
for the complete Python and R command sequence.

### Exact continuous inference

```bash
PYTHONPATH=. venv/bin/python scripts/v13/run_continuous_phospho_inference.py \
  --config config/dct_asdn_phospho_reanalysis.yaml \
  --output-dir data/results/run_20260729_v13_continuous_phospho_exact_final
```

### Read-only report and provenance manifest

```bash
MPLCONFIGDIR=/tmp/rrrm2-mpl-cache \
XDG_CACHE_HOME=/tmp/rrrm2-xdg-cache \
venv/bin/python scripts/v13/build_inference_report.py \
  --input-dir data/results/run_20260729_v13_continuous_phospho_exact_final

venv/bin/python scripts/v13/finalize_inference_provenance.py \
  --run-dir data/results/run_20260729_v13_continuous_phospho_exact_final
```

### Manuscript

```bash
cd latex_paper
latexmk -pdf -interaction=nonstopmode -halt-on-error manuscript_v13.tex
```

## Definitive outputs

The exact run is:

`data/results/run_20260729_v13_continuous_phospho_exact_final/`

Core artifacts:

- `claim_tier.tsv`
- `claim_gates.tsv`
- `set_level_permutation_inference.tsv`
- `parent_gene_null_calibration.tsv`
- `leave_one_gene_out.tsv`
- `matched_observability_sensitivity.tsv`
- `manual_anchor_exclusion_audit.tsv`
- `fixed_universe_attrition.tsv`
- `normalization_*audit.tsv`
- `manifest.json`
- `provenance_manifest.json`
- `reporting/claim_decision_summary.md`
- `reporting/primary_compartment_enrichment.pdf`
- `reporting/leading_parent_gene_matrix.pdf`

The exact sample/channel and phosphoform tables are under:

`data/results/run_20260728_osd462_stage0/reporting/`

## Testing

Run:

```bash
venv/bin/python -m pytest -q
```

Current full-suite status:

```text
199 passed, 0 failed
```

## Claim boundaries

The current repository does not support:

- “Spaceflight suppresses NCC/SPAK regulatory phosphorylation” as an OSD-462
  MS result.
- “Suppression extends beyond NCC/DCT1.”
- “DCT2/CNT phosphoregulation was suppressed.”
- “ASDN pathway activity decreased.”
- “Podocyte phosphorylation was suppressed.”
- transporter flux, electrolyte, aldosterone, renin, or blood-pressure
  conclusions.
- cell-of-origin localization from whole-kidney phosphoproteomics.

The next decisive experiment requires counterbalanced or label-swapped
isobaric tags, isolated-site targeted measurements, segment-resolved sampling,
and matched renal physiology.
