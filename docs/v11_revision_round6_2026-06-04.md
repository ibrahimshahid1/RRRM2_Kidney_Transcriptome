# v11 revision — round 6 (final punch-list), 2026-06-04

Build: `latexmk` exit 0, 35 pages, 27/27 citations resolve. All items from the latest review addressed.

- **§4.2 DCT1 caveat de-tripled.** The third in-section restatement is gone; the permutation status is now one sentence ("the DCT1-high bin is strong at the row level but does not survive the matched permutation while the DCT2-leaning bin does, so the supported claim is … led by the DCT2-leaning extreme, with DCT1 biologically anchored but more conditional").
- **Parent-protein OR discrepancy resolved.** §4.2 now names both analyses and both ORs together — per-site normalization (OR 1.52) and the paired animal-level model (OR 1.56) — and states explicitly they are two distinct analyses whose ORs differ slightly but agree in direction.
- **§3.8 made factual.** The interpretive clauses ("not a stable post-flight signature," "did not support an accelerated-aging reading") are removed from Results; the interpretation lives in the Discussion Context-analyses block.
- **Limitations cut ~55%** (≈500 → 223 words), three tight paragraphs, every load-bearing caveat retained (whole-kidney/parent-gene prior, DCT1-attenuates-DCT2-stable, composition overadjustment, IRI-not-spaceflight, coverage gaps, exploratory human/CMap).
- **Figure 5C (dashboard Panel C) now spatially separates regulatory vs non-regulatory phosphosites** — NCC N-terminal + SPAK/OSR1 regulatory sites grouped above a dashed divider (all suppressed, left), non-regulatory NCC sites below (near zero), with right-side group labels. Re-rendered from `05_plot_dashboard.py` and verified; the "regulatory sites all fall left" pattern is now legible without cross-referencing the table.

## Remaining (your call, flagged as lower priority)
- **Methods §2.1–2.10 collapse to ~3 sections.** Still the one structural item outstanding; it's mostly header consolidation + paragraph joining, but several in-text "§2.x" references and the per-family testing table would need their cross-references checked, so it's best done deliberately for the submission version. I can take it on whenever you want.
