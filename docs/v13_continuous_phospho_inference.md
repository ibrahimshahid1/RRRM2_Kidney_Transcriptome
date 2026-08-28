# V13 continuous parent-gene phosphoproteomic inference

## Purpose

This analysis answers one prospectively frozen question:

> After excluding the established NCC/DCT1 axis, is flight-associated
> phosphosite suppression disproportionately concentrated in independently
> defined DCT2/CNT-transition or aldosterone-sensitive distal-nephron programs?

The implementation is
[`src/v13/continuous_phospho_inference.py`](../src/v13/continuous_phospho_inference.py).
The executable wrapper is
[`scripts/v13/run_continuous_phospho_inference.py`](../scripts/v13/run_continuous_phospho_inference.py).
All production choices are read directly from
[`config/dct_asdn_phospho_reanalysis.yaml`](../config/dct_asdn_phospho_reanalysis.yaml).
The inference code does not build or alter the external subtype signatures.

## Primary estimand

No nominal phosphosite p-value threshold is used to define the primary
endpoint.

1. Load the official scaled signal from `siteQuant_360`.
2. Log2-transform it and median-center each channel within plex.
3. Apply the label-blind site filter:
   localization score at least 19, one modified residue per analytical
   phosphoform, no composite feature, at least 16/20 finite values, and at
   least 6/10 finite values in each plex.
4. De-duplicate by accession plus modified-peptide/phosphoform identity.
5. For every balanced label assignment, calculate FL-GC separately in each
   plex. Each plex must have at least 3 finite reassigned FL and 3 finite
   reassigned GC values. The site effect is the equal-weight mean of the two
   valid plex contrasts.
6. Collapse sites to a parent-gene score:

   `median(-site FL-GC effect)`.

   Positive values therefore mean lower phosphorylation in flight. Signed mean
   and one-sided maxmean are sensitivity scores.
7. Calibrate every parent gene against its own label-permutation null:

   `Z_g = (G_observed - mean(G_null)) / sd(G_null)`.

8. Test a frozen set competitively:

   `mean(Z in set) - mean(Z in eligible non-set background)`.

   The self-contained set mean is secondary.

The primary randomization p-value is obtained by rerunning the complete site,
gene, and set statistic under every balanced within-plex assignment. With two
10-channel plexes and 5 FL/5 GC labels in each, there are
`choose(10, 5)^2 = 63,504` assignments.

## Fixed inference universe under missingness

The 6/10 label-blind completeness rule does not mathematically guarantee that
every site will satisfy 3/3 under every reassignment. The code therefore:

- records the number of valid site and parent-gene statistics for every
  assignment;
- records each site's fraction of valid assignments; and
- fixes each score's inferential gene universe to genes with a finite observed
  Z and a valid parent-gene statistic in every assignment used for null
  calibration.

Set membership consequently does not vary across permutations. In sampled
smoke mode, that guarantee applies only to the sampled assignments; smoke
outputs are never confirmatory results.

## Frozen multiplicity and robustness

- DCT2/CNT-transition and ASDN are one primary family. Their competitive
  statistics use maxT family-wise adjustment.
- Major kidney compartment comparators use Benjamini-Hochberg adjustment in
  their declared family.
- Canonical exclusions are nested: the primary exclusion removes `Slc12a3`,
  `Stk39`, `Oxsr1`, and `Wnk4`; the strict exclusion additionally removes
  `Wnk1`.
- Leave-one-gene-out statistics rerun the permutation comparison for every
  member.
- The secondary annotation null matches set size within bins of site count,
  intensity, and missingness. It is restricted to the primary profile, parent
  score, canonical exclusion, and frozen primary sets; it does not replace the
  animal-label null.
- True multi-compartment broad-expression exclusion, the stricter
  non-distal-or-broad diagnostic, localization score 13, position-level
  de-duplication, a deduplicated multi-modified/composite phosphoform universe,
  uncentered official-scaled signal, summed S/N signal, and parent-protein
  subtraction are explicit sensitivity profiles.

Composite localization strings are parsed component by component and the
minimum component score must pass. The Stage-0 table is also checked for
isolated qualified NCC/SPAK anchor phosphoforms. There are none in the current
workbook, so the configured strict feature removal is an explicit, recorded
no-op rather than an unimplemented sensitivity.

Parent-protein subtraction is recalculated within every label assignment from
the protein workbook. A fixed observed protein coefficient is not subtracted
from permuted phosphosite effects.

## Normalization equivalence

The code audits the ratio between official scaled signal and summed S/N for
every row and plex. In the current workbook this ratio is effectively constant
within a row/plex, so uncentered log2 FL-GC contrasts are algebraically
equivalent. Scaled-versus-summed output is therefore a provenance audit, not
independent validation. Channel-centered versus uncentered is the meaningful
normalization sensitivity.

No normalization resolves the study's condition-to-reporter-tag aliasing. The
manifest carries that limitation explicitly.

## Claim gates

`claim_gates.tsv` reports pass, fail, or non-evaluable for each frozen gate:

- at least eight observable genes;
- positive competitive effect;
- primary-family maxT below alpha;
- every leave-one-gene-out competitive effect positive;
- centered-versus-uncentered and official-scaled-versus-summed-S/N direction
  concordance;
- one complete leave-one-gene-out result for every observable set member;
- positive true multi-compartment broad-expression-exclusion result;
- positive strict canonical-exclusion result;
- an independently derived `DCT2_CNT_external_validation` direction for a
  DCT2/CNT claim; and
- no unrelated kidney compartment with an equal or larger competitive effect.

If `DCT2_CNT_external_validation` is absent, the DCT2/CNT gate is
non-evaluable—not failed. `claim_tier.tsv` reports exactly one of:
`DCT2_CNT`, `ASDN_only`, `neither`, or `non_evaluable`.

The statistical tier is separate from publication promotion. The machine
output cannot permit a promoted biological title while condition remains
perfectly aliased with reporter-tag blocks. A title claiming extension from a
canonical NCC/SPAK anchor is additionally blocked because Stage 0 found no
isolated qualified anchor phosphoform. When DCT2/CNT is under-covered and ASDN
fails, the tier is `neither`: no beyond-axis claim is supported, while the
DCT2/CNT component is still described as non-evaluable rather than null.

## Run commands

Deterministic smoke run:

```bash
venv/bin/python scripts/v13/run_continuous_phospho_inference.py \
  --config config/dct_asdn_phospho_reanalysis.yaml \
  --output-dir data/results/run_YYYYMMDD_v13_continuous_phospho_smoke \
  --smoke
```

The smoke mode uses the observed assignment plus 31 deterministic balanced
null assignments and 100 matched-observability draws. It exercises every
selected profile and output contract but must not be cited.

Exact production run:

```bash
venv/bin/python scripts/v13/run_continuous_phospho_inference.py \
  --config config/dct_asdn_phospho_reanalysis.yaml \
  --output-dir data/results/run_YYYYMMDD_v13_continuous_phospho_exact
```

Run selected profiles while debugging:

```bash
venv/bin/python scripts/v13/run_continuous_phospho_inference.py \
  --config config/dct_asdn_phospho_reanalysis.yaml \
  --output-dir /tmp/v13_primary_smoke \
  --profile primary \
  --smoke
```

Do not launch or interpret the exact run until the flight-blind
`frozen_gene_sets.tsv` and its independent validation set pass their own
reference-builder gates.

## Output contract

| Artifact | Meaning |
|---|---|
| `set_level_permutation_inference.tsv` | Competitive and self-contained set statistics, empirical p-values, maxT, and comparator BH |
| `parent_gene_null_calibration.tsv` | Observed gene scores, gene-specific null mean/SD, Z, observability, and fixed-universe eligibility |
| `leave_one_gene_out.tsv` | Per-member deletion statistics and empirical p-values |
| `matched_observability_sensitivity.tsv` | Secondary gene-annotation null |
| `site_universe_<profile>.tsv` | Exact sites/phosphoforms used in each profile and permutation-validity frequency |
| `site_filter_audit.tsv` | Label-blind exclusion counts |
| `assignment_validity.tsv` | Valid site/gene counts for every assignment |
| `fixed_universe_attrition.tsv` | Site/gene attrition imposed by all-assignment estimability |
| `normalization_*audit.tsv` | Scaled-versus-summed algebraic-equivalence diagnostics |
| `claim_gates.tsv` | Machine-readable gate status |
| `claim_tier.tsv` | Statistical tier plus separate publication-promotion decision and blockers |
| `manual_anchor_exclusion_audit.tsv` | Stage-0-qualified anchor count and strict-exclusion equivalence check |
| `manifest.json` | Input hashes, commit, frozen estimands, limitations, and run mode |
| `provenance_manifest.json` | Post-run code, input, output, commit, and complete dirty-state hashes; generated after reporting |

After generating the read-only report, finalize provenance so every delivered
artifact is hashed without recursively hashing the provenance file itself:

```bash
venv/bin/python scripts/v13/finalize_inference_provenance.py \
  --run-dir data/results/run_YYYYMMDD_v13_continuous_phospho_exact
```

## Definitive production outcome

The corrected exact production run is
`data/results/run_20260729_v13_continuous_phospho_exact_final/`.
It completed all nine declared profiles and all 63,504 balanced assignments.
The machine-readable decision is:

- overall statistical tier: `neither`;
- DCT2/CNT: `non_evaluable` (5 observable genes; minimum 8);
- ASDN: `fail` despite a positive primary statistic, because the podocyte
  comparator is stronger;
- publication promotion: blocked by reporter-tag aliasing and zero isolated
  Stage-0-qualified canonical NCC/SPAK phosphoforms.

The read-only figures, robustness tables, and claim summary are under the
run's `reporting/` directory. `provenance_manifest.json` hashes the delivered
run outputs and the code and inputs used to generate them.

## Interpretation boundary

The unit attached to a subtype set is a parent protein, not a cell. Even a
passing DCT2/CNT set result means that parent proteins associated with an
external DCT2/CNT program carry disproportionate phosphorylation suppression
in whole-kidney tissue. It does not localize those sites to DCT2/CNT cells.
