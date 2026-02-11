# RRRM2 Kidney Transcriptome — Plot Outputs

Generated: 2026-02-11T03:36:38

## Summary

- **Total plots generated**: 23
- **Total tables generated**: 22

## Phase Breakdown

### phase0_deconvolution
- Plots: 6
- Tables: 3

### phase3_rewiring
- Plots: 4
- Tables: 5

### phase5_silent_shifters
- Plots: 3
- Tables: 5

### phase6_uncertainty
- Plots: 3
- Tables: 3

### phase6_regression
- Plots: 7
- Tables: 6

## File Formats

- **PNG**: High-resolution (300 DPI) figures for publication
- **TSV**: Tab-separated tables for downstream analysis

## Usage

To regenerate these plots:
```bash
python scripts/plot_output.py --out_root plots --tag <your_tag>
```