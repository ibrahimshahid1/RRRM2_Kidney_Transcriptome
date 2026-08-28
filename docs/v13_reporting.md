# V13 inference reporting

`src/v13/reporting.py` is a read-only reporting layer for a completed v13
continuous phosphoproteomic inference directory. It does not refit sites,
recalculate parent-gene scores, rerun permutations, alter frozen memberships,
change multiplicity correction, or reevaluate claim gates.

Run it after either a smoke or exact inference:

```bash
MPLCONFIGDIR=/tmp/rrrm2-mpl-cache \
venv/bin/python scripts/v13/build_inference_report.py \
  --input-dir data/results/<v13-inference-run>
```

The default destination is `<v13-inference-run>/reporting`. A separate
destination can be supplied with `--output-dir`.

## Outputs

- `claim_decision_summary.tsv` and `.md`: one row per frozen primary claim,
  including the claim tier, exact gate status, primary statistic, empirical
  p-value, family-adjusted p-value, failed or non-evaluable gates, statistical
  eligibility, the separate publication-promotion status and blockers, design
  warning, and permitted-claim boundary. This keeps DCT2/CNT coverage
  non-evaluability distinct from an ASDN claim-gate failure.
- `primary_compartment_enrichment.tsv`, `.png`, `.svg`, and `.pdf`: the declared
  primary competitive test across the nephron and kidney-compartment
  comparators. The plotted horizontal bars are permutation-null intervals,
  not effect-size confidence intervals. Sets below the frozen minimum
  observable-gene count are shown explicitly as non-evaluable.
- `leading_parent_gene_matrix.tsv`, `.png`, `.svg`, and `.pdf`: observable ASDN and
  DCT2/CNT members ranked by primary parent-gene Z. The table retains raw
  score, site count, parent-protein-subtracted sensitivity where available,
  frozen subtype/ASDN annotations, and the true multi-compartment
  broad-expression flag. The flag is derived first from the run-internal
  `exclude_multicompartment_broad_expression` profile so that a subsequently
  rebuilt external reference cannot change an old report. The stricter
  `exclude_non_distal_or_broad_expression` profile remains a separate
  diagnostic and is not mislabeled as a broad-scaffold exclusion.
- `robustness_profiles_exclusions.tsv` and `.md`: the frozen primary score and
  competitive test across every executed profile and canonical-axis
  exclusion. Missing set rows are represented as non-evaluable rather than
  silently omitted.
- `report_manifest.json`: source run, smoke/exact status, output paths, and an
  explicit statement that inference was not recomputed.

Smoke reports carry a visible warning and must not be cited as inferential
results. Exact label permutation still cannot remove the OSD-462
condition-to-reporter-tag alignment recorded in the inference manifest.
Reporting also surfaces the summed-S/N normalization-direction gate and
leave-one-gene-out completeness gate directly from `claim_gates.tsv`; it never
recomputes either.
