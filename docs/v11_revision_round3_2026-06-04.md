# v11 revision — round 3 (second critique list + §3.8 follow-ups), 2026-06-04

Build: `latexmk` exit 0, **34 pages**, abstract **320 words** (was 398 → 334 → 320), 22/22 citations resolve, all new labels (`fig:axis`, `app:specs`, `tab:modelspec`, `fig:tmtqc`) resolve. Every statistic preserved; added numbers were all sourced from repo artifacts, none invented.

## Tone — lead with the finding, qualify after
- **Abstract** now leads the DCT-prior result with the **DCT2-leaning bin (survives the matched permutation)** as the headline and states DCT1 honestly (OR 1.51 at the row level, **does not survive**, p=0.107).
- §3.4 / §4.3 / Table 6 caption / Fig 5 caption now use **one consistent DCT1 framing**: supported by most sensitivity tests, not by the matched permutation.
- §3.5 — removed the repeated "bulk-composition sensitivities, not cell-type-resolved" refrain (kept once, in Limitations).
- §3.7 — all four perturbation branches rewritten to **lead with the finding**; the vasopressin/KLHL3 paragraph now frames the flat WNK/SPAK/NCC abundance as a **meaningful mechanistic negative** (argues against degradation, points to activity-level regulation).
- §4.1 — strengthened the functional chain ("possible" removed): reduced NCC phosphorylation → reduced activation → coupled K/Ca/Mg handling is a known chain; only its **magnitude** is unmeasured here.

## §3.8 — kept and reframed (per your two follow-ups)
Retitled "Network, cross-species, and translational constraints." Opens on the two constraints instead of "retained as context." The LIONESS/node2vec null is now stated confidently as a **methodological integrity check**; the Twins urine concordance is framed as a **preliminary translational signal** (three concordant physiological axes) with an explicit statement of what would validate it (machine-readable inflight human urine/distal-nephron data across subjects). CMap stays brief as hypothesis-generation.

## Biology now explained
- **New Figure 1 (TikZ schematic)** of the WNK → SPAK/OSR1 → NCC axis, plus an Introduction paragraph explaining the N-terminal activation gate and the **DCT1 vs DCT2** distinction (DCT1 = WNK–SPAK–NCC electroneutral NaCl axis; DCT2/CNT = aldosterone-sensitive + Ca/Mg).
- §3.3 explains the **bulk-dilution** mechanism (expanding interstitial/endothelial cells capture more reads, lowering apparent DCT levels).
- §3.7 states the **low-K → NCC activation** logic that makes the anti-alignment meaningful.

## Statistical accuracy (sourced from repo)
- **KSEA atlas substrate counts** added: WNK1 (101), WNK3 (208), SPAK/STLK3 (302), OSR1 (221) — lets the reader weigh reliability and shows SPAK/OSR1 have *more* substrates yet aren't suppressed.
- **VIF concern addressed**: the DCT1 top-decile OR stays >1 and significant under each covariate added individually (M1 1.55, M2 1.36, M3 1.24), so the result isn't a collinear-full-model artifact (from `h2_composition_adjusted_suppression_enrichment_single.tsv`).
- **cosine 0.87 vs meta p=0.19 reconciled** in §3.1: the high anchor cosine is whole-panel concordance dominated by remodeling, and OSD-462 is the matched anchor (not an independent recurrence cohort); the transport set is isolated in the meta-analysis, where it is the weaker half.

## Prose / structure
- §3.2 now **interprets** the cross-layer table (the matrix/ECM proteins-opposite-transcripts surprise) instead of transcribing it; the endothelial↔NCC-phospho relationship is flagged as the strongest animal-level finding and foreshadowed.
- Methods **Tables 1 & 2 moved to an Appendix**; **Figure 7 (TMT QC) moved to the Appendix**; §2.11 global-FDR justification cut to one sentence.
- "concordant" disambiguated (defined at first use; §3.3 reworded to "moved together").

## Citations
- **Johnson et al. 2023** kinase atlas added (ref) and cited in §2.7 and §3.2.
- **Siew 2024 cited in Results §3.2** with the honest framing that the NCC phosphosite suppression is a same-animals (RR-10/OSD-462) **re-derivation** of the antibody-validated endpoint, not independent validation.
- Twins **figure-level AQP2** limitation made explicit at first cite (§2.10).
- Title is now **finding-forward**.

---

## Still outstanding — needs the figure-rendering pipeline or your call
These require re-running `src/v11/publication_figures.py` (and the R DCT-prior script) on the data; I did the LaTeX-level moves and caption fixes but did not regenerate image contents:

1. **Figure 2A "8,64" superscript** — the script title is clean (`publication_figures.py:294` sets just "Canonical recurrence estimates", no superscript), so this is a **stale render artifact**; re-rendering Figure 2 should clear it. Worth confirming there isn't a second annotation adding it.
2. **Figure 2B** label legibility (ggrepel/nudge) — re-render.
3. **Figure 8 split** — Panel B (14-gene heatmap) and Panel D (IRI Visium) into standalone figures; Panels A/C to supplement — re-render.
4. **Figure 6 Panel B** (continuous-gradient) removal → full-width Panel A — re-render.
5. **Two more new schematics** you asked for — a data-layer summary (I can draw this in TikZ next) and a DCT1/DCT2 ranked-score panel (needs the `D_g` values from `gse228367_gene_stats_overall.tsv`); say the word and I'll add them.
6. **DESeq2/limma/edgeR per-cohort attribution** — the repo uses limma for the LIONESS edge regression and OSDR-provided DE tables for OSD-462; a clean per-cohort tool map isn't documented, so I left the standard-practice citation rather than assert specifics.
7. **Age-by-flight stability metric** (from round 2) — still no clean quantitative value in the repo; the text reads "did not hold a consistent sign under bootstrap resampling."
