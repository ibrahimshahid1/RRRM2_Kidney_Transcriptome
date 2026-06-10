# v11 revision — round 5 (hero result, two new figures, biology + structure), 2026-06-04

Build: `latexmk` exit 0, 35 pages, abstract **300 words**, 27/27 citations resolve. Two new figures built from real repo data and verified visually.

## Two new figures (the ones you said were missing)
- **Cross-layer directional heatmap** (new, in §3.2): an 11-pathway × 3-layer (RNA/protein/phospho) grid, each cell colored by within-layer flight effect (z) and marked ↑/↓/≈, ordered by RNA effect, with the RNA→protein discordance now visually unmistakable (ECM red at RNA, blue at protein; DCT/NCC-WNK blue at RNA, red at protein). Built from `cross_layer/osd462_cross_layer_pathway_matrix.tsv`. Table 1 is retained for exact counts; the figure carries the impact.
- **Per-animal endothelial vs NCC-regulatory-phosphorylation scatter** (new, in §3.3): the n=20 "money shot," flight (red) vs ground control (blue), regression line, ρ = −0.762, p = 0.0004 — and the script recomputed ρ to exactly −0.762 from `h3_mediation/h3_mediation_input_scores.tsv`, confirming the reported value.

## RNA–protein mismatch promoted to the hero result (§3.2)
Added the explicit framing: a reader with RNA alone would conclude matrix-remodeling proteins accumulate, but the protein layer moves the opposite way, so the RNA signature is a tissue-context readout, "actively misleading about protein remodeling" — only protein/phospho recover the lesion. The compartment-dilution mechanism now carries a deconvolution citation \citep{newman2015cibersort}, and the matrix protein test is named (matched gene-set null).

## Figure 1 schematic now shows activation state
Redrawn to distinguish **NCC–P (active → Na⁺/Cl⁻ ON)** from **NCC (inactive → OFF)**, with the red "spaceflight: −P" arrow shifting active→inactive at constant total protein. Verified rendered.

## Discussion restructured + biology
- **Merged §4.1/4.2** (main interpretation + prior work) and **§4.3/4.4** (subtype-prior + upstream mechanism) into continuous blocks.
- Removed "interpretive as much as confirmatory"; the contribution is now stated plainly: *first matched RNA+protein+phospho decomposition of the spaceflight kidney response, with a non-obvious result (RNA and protein disagree; phospho aligns)*.
- **New DCT2-biology paragraph**: why a DCT2-leaning bin can be the permutation-robust one — flight may suppress a broader aldosterone-sensitive distal program (SGK1–NEDD4L→ENaC, TRPV5/CALB1 calcium handling, WNK–SPAK→ROMK/KCCs), so DCT2/CNT-expressed substrates lose phosphorylation without NCC being DCT2-specific. Frames it as a testable hypothesis.
- **Conclusion rewritten** as a field-level "so what": bulk transcriptomics alone will misread the remodeling RNA signature and miss the phospho event, so the field should prioritize phospho / DCT-resolved / spatial phosphoproteomics paired with electrolyte physiology — not more bulk RNA.

## Introduction biology + citations
- Foundational **WNK–SPAK/OSR1–NCC** citations added \citep{vitari2005wnk,moriguchi2005wnk,richardson2008wnk}, plus the **Gitelman-syndrome physiological hook** (salt wasting, hypotension tendency, hypokalemia/hypomagnesemia) \citep{gamba2005cct} — the spaceflight-medicine stakes.
- One sentence on **what Cosmic Kidney Disease left unresolved** (RNA→protein propagation; which layer the signal concentrates in).
- Hypothesis phrasing fixed ("would not be recapitulated," not "not necessarily").

## Abstract restructured (326 → 300 words)
Findings-forward arc; dropped "in matched public data"; the layer-decomposition summary moved to sentence 3; the I² stat removed; the RNA–protein mismatch ("a reader with RNA alone would wrongly infer protein-level remodeling") and the DCT2-survives-permutation result each get a dedicated early sentence; contribution stated ("first matched cross-layer decomposition").

## Methods QC/citations
- DCT prior nuclei counts added: **18,881 DCT1 and 6,274 DCT2 nuclei** (three replicates).
- Gene-set Z-score acknowledged as a standard per-sample standardized mean, not novel.
- KSEA atlas usage clarified (full Johnson et al. atlas substrate set, ≥3 quantified substrates, e.g. 101 for WNK1 / 209 for WNK3; custom implementation).
- **Figure 2 moved from Methods into Results §3.1** (beside the recurrence figure).

## Remaining (the larger lifts, your call)
1. **Full Methods §2.1–2.10 collapse to ~3 sections** — I merged the Discussion pairs; the Methods 10→3 merge is the next big restructure (mostly header removal + paragraph joining, but worth doing carefully).
2. **Discussion §4.5/4.6 merge and a ~40% Limitations cut** — partially addressed; the remaining merge + trim is straightforward but sizable.
3. **In-image figure re-renders**: separate regulatory vs non-regulatory phosphosites in the OSD-462 phospho panel; visually de-emphasize the intermediate rows of the DCT evidence ladder so the matched-permutation row dominates. Both are `publication_figures.py` edits + re-render, like the panels I already fixed.
4. **§3.3 → §3.2 reorder** (compartment dilution before the mismatch table) — I added the dilution mechanism into §3.2's hero paragraph, which softens the inversion, but the full reorder is still open.
