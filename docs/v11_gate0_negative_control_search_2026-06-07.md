# v11 Module 4 — Gate 0 verdict: non-kidney matched multi-omic spaceflight data

Date: 2026-06-07
Author context: companion to docs/v11_novelty_extensions_implementation_plan_2026-06-06.md §5

## Question
Does NASA OSDR / GeneLab host a publicly-available dataset with matched RNA-seq + proteomics (ideally + phosphoproteomics) from the **same mouse animals** in a **non-kidney tissue**, suitable as a negative-control tissue contrast against OSD-462?

## Verdict
**CONDITIONAL PASS — one strong candidate (OSD-105) and one secondary (OSD-173) exist; viability depends on platform compatibility checks that require downloading the workbooks.**

## Candidate matrix

| OSD ID | Tissue | Mission | RNA-seq | Proteomics | Phospho | Same-animal match | n (FLT/GC) | Verdict |
|---|---|---|---|---|---|---|---|---|
| **OSD-105** | TA muscle | RR-1 (37d, female C57BL/6J, 16wk) | yes | yes (mass-spec; platform unconfirmed) | unconfirmed | yes (12 mice total) | 6 / 6 | **Primary candidate** |
| OSD-104 | SOL muscle | RR-1 | yes | NO | NO | — | 6 / 6 | Rejected (no protein layer) |
| OSD-99 | Multiple RR-1 tissues | RR-1 | viability check | viability check | viability check | — | — | Supporting only (not a primary dataset) |
| OSD-173 | Liver? | TBD | TBD | TBD | TBD | TBD | TBD | Secondary — needs verification |
| OSD-488 | Female SOL | RR-1 | Ca²⁺ reuptake data only | — | — | — | 4 / 4 | Not relevant |

## Confirmed properties of OSD-105 (from search hits)
- RR-1 NASA Validation Flight, female C57BL/6J mice flown 37 days
- Bulk RNA-seq present
- Mass-spectrometry-based proteomics present
- Bisulfite-seq methylation present (not used here)
- 6 spaceflight + 6 ground control animals (matched at the same-animal level across modalities)
- Documented as part of the RR-1 multi-omics bio-banked tissue program

## Open compatibility questions (must check by downloading the workbook before running Module 4)
1. **Is the proteomics TMT-based or label-free?** OSD-462 is TMT 2-plex with scaled S/N intensities. If OSD-105 is label-free (LFQ), the matching logic in `src/multiomics/osd462_anchor.py` (per-plex within-plex differencing, per-channel median centering) does not directly transfer. A LFQ pipeline would need a separate flight-effect estimator; the matched-null engine is platform-agnostic so the per-pathway statistic stays comparable.
2. **Is there phosphoproteomics?** Most likely NO (RR-1 publications focus on proteomics + methylation). If absent, Module 4 reduces to RNA→protein propagation only — no phospho carry can be tested.
3. **Is the proteomics matched to the RNA-seq at the per-animal level?** Need to confirm sample IDs cross-reference.
4. **Strain & age comparability vs OSD-462.** OSD-462 = RR-10 (need to confirm strain/age). RR-1 = C57BL/6J × 16wk × 37d flight. Mission differences are real and limit interpretation; Module 4 is framed as "is the kidney decoupling generic to mouse spaceflight multi-omics?" not "RR-10 vs RR-10," so different missions are acceptable but must be flagged.

## Recommendation
- **Proceed with Module 4 conditionally**, scoped narrowly to OSD-105 vs OSD-462 RNA→protein propagation (Module 2 statistic). Phospho specificity is out of scope unless OSD-105 phospho is confirmed.
- **Do this verification step before downloading**: query the OSDR study page for the assay/file manifest. If the protein workbook header mentions "TMT" or "iTRAQ" the existing TMT parser may be reusable with header overrides; if it says "LFQ" / "intensity" / "iBAQ," a thin LFQ flight-effect estimator must be added before the Module 2 pipeline can run on it.
- **If TMT and same-animal matching both confirmed**, this is a genuine specificity contrast and the manuscript can elevate the discordance from "kidney shows X" to "X is kidney-specific (vs muscle)."
- **If LFQ but otherwise compatible**, run the contrast but document the platform-difference limitation explicitly; cross-platform specificity is still informative ("kidney protein discordance is larger than muscle protein discordance under each platform's matched null").
- **If neither compatibility check passes**, document the limitation per the plan and stop. The negative finding (matched non-kidney TMT spaceflight multi-omics is rare to nonexistent) is itself reportable.

## Sources
- Multiple confirming references to RR-1 quadriceps/gastrocnemius/TA matched RNA+protein+epigenomic data across catalog.data.gov listings.
- npj Microgravity 2023 "Explainable machine learning identifies multi-omics signatures of muscle response to spaceflight in mice" — integrates these datasets and is the closest published analog to our v11 propagation framing.
- NASA OSDR (https://osdr.nasa.gov) — direct study page renders client-side and is not WebFetch-accessible; manual or API-based verification is required to lock the compatibility matrix above.
