# v11 revision — round 7 (accuracy, tone, figure honesty), 2026-06-04

Build: `latexmk` exit 0, 35 pages, 28/28 citations resolve.

## Critical fix
- **Abstract factual bug corrected.** The sentence implying matrix RNA fell is replaced: "Matrix/ECM RNA rose while its protein moved in the opposite direction, whereas DCT/NCC-WNK RNA fell while total NCC and related protein abundance stayed flat-to-slightly-positive." Now matches Results exactly.

## Tone / positioning
- **Title** softened from "does not propagate" to **"Spaceflight kidney remodeling RNA is decoupled from protein abundance, while NCC/WNK transporter suppression resolves at phosphorylation"** (no longer nitpickable by the one near-zero concordant pathway).
- **Conclusion** de-aggressived: "transcriptomic surveys … will misread/will miss" → "bulk RNA alone is insufficient to infer protein remodeling or transporter activation state." Keeps the point without sounding anti-RNA-seq.
- **Defensive caveats trimmed** (e.g., CMap "context-level class descriptions, not treatment recommendations" → "we interpret CMap outputs only as perturbagen-class context").

## Statistical / scientific honesty
- **§3.5 retitled** "DCT1 enrichment persists under abundance and composition adjustment but remains observability-sensitive" (no longer implies DCT1 is robust overall).
- **DCT2-leaning guardrail** added: defined by lowest DCT1 score, captures Trpv5/Calb1 at its expected end but may include broadly expressed genes — labeled "DCT2-leaning," not DCT2-specific.
- **Targeted-KSEA p-values de-emphasized**: the 3-substrate $z$ scores are kept as directional coherence; the extreme $p$ values ($2.8\times10^{-10}$ etc.) are removed and explicitly "not interpreted inferentially."
- **Channel-centering tradeoff** stated: median-centering may also remove a real global dephosphorylation, so the interpretable signal is the subtype-prior *enrichment*, not absolute global dephosphorylation magnitude.
- **Fig 4 caption** softened to "no pathway shows strong, nontrivial RNA-to-protein propagation."

## Figures
- **Fig 5A relabeled** "Cross-cohort RNA → matched-anchor protein" (honest about RRRM-2 RNA vs OSD-462 protein), with a caption note; re-rendered.
- **Fig 6B relabeled** "matrix recurs, transport heterogeneous" (was "transport does not"); re-rendered and verified.
- **Fig 7 (endothelial↔NCC scatter)**: caption now reports the within-condition correlation (pooled ρ=−0.67; GC −0.73, flight −0.20), showing ρ=−0.762 is not just flight/control separation.
- **Fig 5C** (last round) regulatory vs non-regulatory phosphosites are spatially separated; caption updated.

## Citations
- Added **Subramanya & Ellison 2014** (distal convoluted tubule review) to the DCT2/CNT mechanistic paragraph (with Richardson & Alessi for WNK–SPAK).

## Remaining (larger structural work, flagged across rounds)
1. **Methods §2.1–2.10 → ~3 sections** (needs in-text §-reference checks).
2. **Dataset table** (assay/mission/sex-strain/n/role/accession) and a **gene-set provenance supplementary table** (every member + per-pathway citation).
3. **Concept/terminology box** early in the paper.
4. **Table 5** reorder so the parent-gene evidence ladder is primary and raw phosphosite-row counts move to supplement.
5. **Move to supplement**: Fig 5 Panel D (network control) and the perturbation-triangulation figure; trim §3.7 detail; trim repeated numbers in the Discussion.
6. **Fig 1** add DCT1/DCT2/CNT context; **Fig 2** plot CI bars in Panel A; **Fig 3** state score units/SEM-vs-CI; **Fig 9** add DCT2 or demote.
7. Full **Results reorder** into the five-beat narrative you outlined.

These are the remaining items; I can take the Methods collapse + the two tables (dataset, gene-set provenance) next, since those are the highest-impact for an editor.
