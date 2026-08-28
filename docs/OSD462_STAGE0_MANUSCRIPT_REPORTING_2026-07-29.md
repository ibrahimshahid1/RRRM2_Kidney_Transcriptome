# OSD-462 Stage 0 manuscript reporting

[`scripts/osd462/09_stage0_manuscript_reporting.py`](../scripts/osd462/09_stage0_manuscript_reporting.py)
is a read-only reporting layer over the frozen Stage 0 artifacts. It does not
read raw workbooks, recompute phosphosite effects, change assay qualification,
or modify manuscript files.

Run:

```bash
venv/bin/python scripts/osd462/09_stage0_manuscript_reporting.py
```

Inputs:

- `data/results/run_20260728_osd462_stage0/osd462_sample_plex_channel_design.tsv`
- `data/results/run_20260728_osd462_stage0/osd462_ncc_spak_phosphosite_provenance.tsv`

Outputs:

- `data/results/run_20260728_osd462_stage0/reporting/osd462_exact_sample_channel_inclusion.tsv`
  contains one row per biological sample/reporter channel. It explicitly marks
  the 20 FL/GC primary-contrast samples, the 10 BL context-only samples,
  presence in both MS modalities, ISA agreement, assay name, raw archives,
  condition/tag-block aliasing, and label-swap status.
- `data/results/run_20260728_osd462_stage0/reporting/osd462_slc12a3_stk39_phosphoform_provenance_compact.tsv`
  retains every Slc12a3 and Stk39 feature from the Stage 0 provenance audit,
  with source row, residue-resolved feature, observed peptide phosphoform,
  localization evidence, reporter completeness, canonical components, and an
  explicit isolated-canonical qualification and reason.
- `figures/v13/osd462_reporter_tag_map.{png,pdf,svg}` shows the exact TMTpro
  tag order, sample assignment, contiguous BL/FL/GC blocks, identical
  assignments across both plexes, and absence of a label swap. PNG is rendered
  at 600 dpi; PDF and SVG are vector outputs.

The script fails rather than silently drawing a figure if the two modalities
do not contain identical 30-sample maps, a plex does not contain the expected
5/5/5 BL/FL/GC blocks, an ISA mapping fails, a label swap appears, provenance
rows are dropped, or the frozen isolated-canonical count differs from zero.

## Interpretation boundary

The reporter-tag map is a design disclosure, not a correction for confounding.
Because the same condition-specific tag blocks are reused in both plexes, the
FL–GC coefficient cannot distinguish condition from a systematic reporter-tag
block effect. Likewise, the compact phosphoform table reports the T53- and
S383-indexed co-modified features as context; it does not convert them into
isolated canonical-site or phospho-occupancy measurements.
