# WGCNA Pivot Manuscript: Figure And Rerun Plan

This plan supports `latex_paper/manuscript_wgcna_pivot.tex`.

## Main Rerun

Run the WGCNA-mode pipeline from the repository root:

```bash
venv/bin/python -m src.run_all_phases \
  --network-method wgcna \
  --run-id run_YYYYMMDD_wgcna_5000g \
  --wgcna-genes 5000 \
  --wgcna-pres-perms 1000 \
  --external-validation \
  --external-studies OSD-102,OSD-513,OSD-163,OSD-253
```

Current caveat: `src/run_all_phases.py` does not yet call the GC-reference preservation script, so run this manually after the WGCNA phase:

```bash
Rscript src/networks/wgcna_gc_reference_preservation.R \
  data/results/<RUN_ID>/phase1_residuals/Rtech.tsv.gz \
  data/results/<RUN_ID>/phase1_residuals/meta_phase1.tsv.gz \
  data/results/<RUN_ID>/wgcna/gc_reference \
  1000
```

Then rerun the WGCNA follow-up and figures:

```bash
RRRM_RESULTS_DIR=data/results/<RUN_ID> venv/bin/python -m src.networks.wgcna_followup

venv/bin/python -m src.visualization.wgcna_publication_figures \
  --results_dir=data/results/<RUN_ID>
```

## Secondary Benchmark Reruns

These are not primary discovery analyses; they support the negative topology-rewiring benchmark.

Leakage-safe LIONESS-vs-expression benchmark:

```bash
venv/bin/python -m src.validation.enhanced_cv \
  --phase2_dir=data/results/run_20260505_remediated_2500g/networks \
  --rtech=data/results/run_20260505_remediated_2500g/phase1_residuals/Rtech.tsv.gz \
  --meta=data/results/run_20260505_remediated_2500g/phase1_residuals/meta_phase1.tsv.gz \
  --outdir=data/results/run_20260505_remediated_2500g/phase8_validation \
  --max_genes=2500 \
  --topk=80 \
  --n_perms=1000
```

Direct differential-correlation diagnostic:

```bash
venv/bin/python -m src.statistics.direct_coexpression_test \
  --rtech=data/results/run_20260505_remediated_2500g/phase1_residuals/Rtech.tsv.gz \
  --meta=data/results/run_20260505_remediated_2500g/phase1_residuals/meta_phase1.tsv.gz \
  --phase2_dir=data/results/run_20260505_remediated_2500g/networks \
  --rewiring_dir=data/results/run_20260505_remediated_2500g/phase3_rewiring \
  --outdir=data/results/run_20260505_remediated_2500g/phase10_direct_coexpr \
  --n_perms=1000
```

## Main Figures

Use these as the manuscript core:

1. `fig1_module_trait_heatmap.png`: WGCNA module-trait model.
2. `fig2_eigengene_grey60.png`: lead grey60 module eigengene by group.
3. `fig6_contrast_forest.png`: FLT-vs-GC simple-effect contrasts.
4. `fig5_grey60_hub_genes.png`: grey60 hub genes by kME.
5. A new GC-reference preservation figure should be added. Current manuscript uses a table because no dedicated figure exists yet.
6. `fig7_external_validation.png`: external WGCNA module projection.

Supplementary figures:

- `supplementary/figS1_eigengene_salmon.png`
- `supplementary/figS2_eigengene_green.png`
- `supplementary/figS3_eigengene_pink.png`
- `supplementary/figS4_eigengene_royalblue.png`
- Phase 8 network-vs-expression benchmark table/plot.
- Direct differential-correlation diagnostic summary.
- Old remediated LIONESS inference tables.

## Abstract Backbone

Use this as the abstract thesis:

> Spaceflight is associated with renal remodeling, but whether kidney transcriptomic responses reflect altered co-expression topology or shifts in preserved transcriptional programs remains unclear. RRRM-2 WGCNA identifies an ISS-T-specific module activity response led by grey60/ME17, a 48-gene ECM/cell-migration module. The module increases under flight in ISS-T young and old mice but not in LAR, and GC-reference preservation shows its co-expression structure remains strongly preserved. Negative LIONESS and direct differential-correlation benchmarks do not support broad topology rewiring. The detectable RRRM-2 kidney response is therefore best described as ISS-T-specific activation of preserved remodeling modules.

## Claims To Keep

- RRRM-2 shows ISS-T-specific module activity shifts.
- Grey60/ME17 is the lead ECM/cell-migration module.
- GC-reference preservation supports activity shift, not topology disruption.
- LIONESS/node2vec and direct differential-correlation benchmarks are negative.
- External WGCNA projection detects ECM modules in OSD-513, but this is context support rather than universal replication.

## Claims To Avoid

- Do not claim broad kidney co-expression topology rewiring.
- Do not claim grey60 topology is disrupted.
- Do not lead with PPAR/lipid as the main RRRM-2 module result.
- Do not claim OSD-102 validates the grey60 response.
- Do not claim OSD-253 cleanly maps onto RRRM-2 ISS-T/LAR; it is a duration/context cohort with caveats.
