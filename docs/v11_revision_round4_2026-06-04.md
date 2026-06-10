# v11 revision — round 4 (figure fixes + final tone/consistency), 2026-06-04

Build: `latexmk` exit 0, 34 pages, abstract 326 words, 22/22 citations resolve, no undefined refs. Three figures re-rendered from the pipeline and verified visually.

## The "8,64" artifact — diagnosed and fixed
It was **not a citation**. `pdftotext` showed the literal text was "Canonical recurrence estimates" followed by "0.64" — the first forest-plot point's **cosine data-label colliding with the panel title**. The script title was already clean, so the committed PDF was a stale render. Fixed by moving the cosine labels to the right of each error bar (`fig_cross_cohort_recurrence`) and re-rendering `fig1_main_result_multipanel`. Verified: title and labels are now cleanly separated (cosines 0.64 / −0.51 / 0.87 to the right of their points).

## Figure panels removed (re-rendered)
- **Composition figure** (`v11_parent_protein_composition_sensitivity`): dropped the continuous-gradient **Panel B** (a null result already summarized in one sentence); now a single full-width adjustment-ladder panel (M0–M5, DCT1 decile + quartile).
- **Perturbation figure** (`v11_perturbation_triangulation`): dropped the redundant parent-protein-normalized **Panel C** (covered by the composition figure / evidence ladder); Panel D (IRI Visium) now spans the bottom row, relabeled C. Caption updated to A/B/C.

## TMT-QC figure + spec tables → formal Supplementary Information
The appendix is now "Supplementary Information" with **S-numbering**: the two specification tables are **Tables S1–S2** and the TMT channel-centering QC figure is **Figure S1** (Section S1). In-text pointers updated; all cross-references resolve.

## Final tone / consistency fixes
- **Abstract** now leads the DCT result with **DCT2 surviving the permutation** ("led by the DCT2-leaning bin … survived …; the DCT1-high bin was strongly enriched at the row level but did not survive that permutation, p=0.107"). DCT1 is no longer front-loaded.
- **§4.1**: removed the undercut ("functional consequence remains an experimental prediction rather than a measured endpoint"); the confident chain stands, and the missing-data point is folded into "what remains to be measured is the magnitude … not the existence of the chain."
- **"concordant" disambiguated**: the term is now reserved for the **RNA–protein** relationship. Cross-cohort uses → "agreement/agree in direction/shared-direction" (Methods §2.3, §3.1); the §3.3 animal-level use → "group-level agreement"; the §3.8 human-vs-mouse use → "directionally aligned." The named "urine concordance" analysis keeps its term (separate translational domain).

## Still genuinely outstanding (your call)
- **New schematics #2/#3** you'd asked about earlier: a data-layer summary (TikZ) and a DCT1/DCT2 ranked-score panel from `gse228367_gene_stats_overall.tsv` — not yet added; say the word.
- **Figure 2B legibility** (pathway-label overlap) — can fix with ggrepel-style nudging on that figure next if you want it.
- **DESeq2/limma/edgeR per-cohort attribution** and the **age-by-flight stability metric** — still not cleanly documented in the repo; left honest/general.
