# LaTeX Manuscript Artifacts

## Current v11 Pair

- `manuscript_v11.tex` / `manuscript_v11.pdf`: current main manuscript draft.
- `results_v11.tex` / `results_v11.pdf`: companion results compendium when regenerated.

The v11 manuscript pulls publication figures first from:

`../data/results/run_20260526_v11_dct1_phospho_mediation/figures/v11/`

and then falls back to older local figure folders for legacy panels.

## Build

From this directory:

```bash
latexmk -pdf manuscript_v11.tex
```

The end-to-end analysis entry point is:

```bash
python -m src.run_all_phases --v11-only --run-id RUN_ID
```
