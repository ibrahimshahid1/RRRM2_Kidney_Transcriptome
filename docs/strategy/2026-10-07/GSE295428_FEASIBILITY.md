# GSE295428 (RRRM-1 single-cell) kidney feasibility: label-blind inventory

**Date:** 2026-10-07
**Question:** Does GSE295428 contain enough confidently identified podocytes per animal to support an
animal-level pseudobulk test of the frozen podocyte HS and structural programs?
**Blinding:** no flight-vs-control effect was computed. Only cell counts, tissue identity and podocyte
counts per sample were examined.

## Data

- GEO GSE295428 per-sample 10x matrices (`barcodes` / `features` / `matrix.mtx`), downloaded
  2026-10-07.
- Kidney samples GSM8947832–39, titled "Kidney, HC, Mouse 1–4" and "Kidney, LAR, Mouse 1–4".
- Spleen-LAR samples GSM8947977–80, used for comparison.
- BALB/cAnNTac females. On 2026-04-22 GEO reset these samples' `tissue`, `age` and `treatment` fields
  to "unknown".

## Method

- **Tissue identity:** fraction of cells detecting any marker of a family.
  - kidney proximal tubule: Lrp2, Slc34a1, Slc22a6, Kap
  - TAL: Umod, Slc12a1
  - DCT/CD: Slc12a3, Aqp2
  - immune: Ptprc, Cd79a, Cd3e
  - B cell: Cd79a, Ms4a1, Cd19
  - T cell: Cd3e, Cd3g
  - thymocyte: Rag1, Dntt
  - plus endothelial and other families
- **Podocytes:**
  - loose rule: Nphs1 or Nphs2 detected, and at least 2 of {Nphs1, Nphs2, Podxl, Ptpro, Synpo, Wt1,
    Mafb, Clic5};
  - strict rule: the same, with at least 4 of those markers.
- **Duplicate check:** barcode overlap between the "Kidney, LAR" and "Spleen, LAR" files of the same
  mouse.

## Results

| Sample | Cells | Median UMI | Podocytes, loose / strict | Proximal-tubule fraction | Immune fraction | B / T fraction |
|---|---|---|---|---|---|---|
| Kidney-HC mouse 1 | 3,290 | 5,049 | 56 / 14 | 0.92 | 0.06 | 0.00 / 0.01 |
| Kidney-HC mouse 2 | 3,931 | 4,391 | 53 / 14 | 0.71 | 0.06 | 0.00 / 0.01 |
| Kidney-HC mouse 3 | 3,917 | 4,535 | 46 / 23 | 0.93 | 0.10 | 0.00 / 0.02 |
| Kidney-HC mouse 4 | 4,016 | 4,286 | 37 / 14 | 0.75 | 0.06 | 0.00 / 0.01 |
| Kidney-LAR mouse 1 | 3,156 | 1,401 | 0 / 0 | 0.00 | 0.96 | 0.58 / 0.31 |
| Kidney-LAR mouse 2 | 1,350 | 5,521 | 2 / 0 | 0.00 | 0.93 | 0.60 / 0.30 |
| Kidney-LAR mouse 3 | 5,689 | 2,743 | 6 / 0 | 0.00 | 0.98 | 0.59 / 0.36 |
| Kidney-LAR mouse 4 | 3,139 | 2,768 | 1 / 0 | 0.00 | 0.98 | 0.61 / 0.34 |

**Barcode overlap between "Kidney, LAR" and "Spleen, LAR", per mouse:** 99.9%, 99.6%, 99.9% and 99.8%
of the smaller file.

## Conclusion

1. **The four "Kidney, LAR" samples are spleen, not kidney.**
   - They contain about 60% B cells, about 30–36% T cells, red-pulp macrophages, and no
     tubular cells.
   - They are near-duplicates of the "Spleen, LAR" libraries.
   - The actual RRRM-1 live-return kidney single-cell data are **not in GSE295428**. This is presumably
     why GEO reset the tissue labels.
2. **The habitat-control kidney samples are genuine kidney**, with 14–23 strict (37–56 loose)
   podocytes per animal.
3. **GSE295428 cannot test a flight or LAR effect in kidney podocytes.** Even if the true LAR kidney
   files were deposited, 14–23 strict podocytes per animal from about 4,000-cell whole-kidney
   libraries, with n = 4 vs 4 and a pool confound (LAR mice 1–2 are in Pool 4), would support at best a
   very low-precision descriptive pseudobulk.
4. **Action:** report the duplicate upload to GEO and the RRRM-1 authors. The project would benefit if
   the correct LAR kidney libraries were released. Until then, this path is closed; the habitat-control
   kidney cells could at most serve as a BALB/c reference for marker detectability.

The scripts used are in the session scratchpad (`scripts_sc/podo_inventory.py`,
`scripts_sc/lar_identity.py`). They are diagnostic only, and the method above is sufficient to
reproduce them.
