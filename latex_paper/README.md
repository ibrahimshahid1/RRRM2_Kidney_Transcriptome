# LaTeX Manuscript Artifacts

## Current v11 Pair

- `manuscript_v11.tex` / `manuscript_v11.pdf`: current main manuscript draft, titled around distal-nephron subtype-prior enrichment.
- `results_v11.tex` / `results_v11.pdf`: companion results compendium for the same v11 artifact set.
- `v11_headline_number_index.tsv`: source-artifact, filter, and column map for reported headline numbers.
- RNA recurrence curation artifacts live under the v11 run root in `baseline/rna_recurrence_gene_set_members.tsv`, `baseline/osd462_rna_recurrence_leave_one_family.tsv`, and `baseline/osd462_rna_recurrence_paired_pathway_bootstrap.tsv`.
- Raw OSD-462 TMT QC summaries live under the same run root in `baseline/osd462_tmt_channel_qc.tsv` and `baseline/osd462_tmt_missingness_by_condition.tsv`.

The v11 manuscript pulls publication figures first from:

`../data/results/run_20260526_v11_dct1_phospho_mediation/figures/v11/`

and then falls back to older local figure folders for legacy panels.

## Build

From this directory:

```bash
latexmk -pdf manuscript_v11.tex
```

```bash
latexmk -pdf results_v11.tex
```

The end-to-end analysis entry point is:

```bash
python -m src.run_all_phases --v11-only --run-id RUN_ID
```
