# v11 revision — round 2 (implementing the section-by-section critique), 2026-06-04

Worked through the full section-by-section critique. Prose/structure only; every statistic, citation, figure, and table number is unchanged except where a value was *added* from the repo (all sourced, none invented) or a duplicated value was consolidated. Build: `latexmk` exit 0, 31 pages, 21/21 citations resolve, abstract 398 → **330 words**.

## Abstract
- Re-inverted the arc: opens on the established endpoint + the **gap** ("Whether the recurring transcriptomic response tracks that endpoint across protein and phosphorylation layers has not been tested systematically in matched public data") instead of a generic remodeling sentence.
- Added the interpretive note to the transport `p=0.19` (limited power **and** cross-cohort heterogeneity, median `I²=73%` — value already in Results).
- Added the mechanistic thread to the RNA–protein mismatch ("bulk dilution by shifting interstitial/endothelial compartments and post-transcriptional regulation").
- Removed the defensive "a DCT-subtype prior is … not a cell-of-origin assignment" clarification (belongs in Methods/Results, not the abstract).
- New ending lands the punchline: phosphorylation is the readout that **survives multi-layer decomposition** and yields a specific spatial prediction.
- Spelled out NCC (sodium-chloride cotransporter); trimmed back to 330 words.

## Introduction
- Smoothed the stone→tubulopathy transition ("shifted the emphasis from stone physiology to active tubular dysfunction") and folded the obvious "rodent studies are useful" sentence into it.
- **NCC spelled out at first use** with gene symbol (`Slc12a3`).
- De-restated the methods sentence (dropped "addressed at three levels").
- Added the **consequences of the localization fork** (DCT-intrinsic → direct transporter target; remodeling-associated → secondary to interstitial/endothelial expansion).

## Methods (data-grounded additions in **bold**)
- §2.1: added the dataset-hierarchy framing (anchor / recurrence / context / reference) and made **why RRRM-2 is primary** explicit (richest design; supports within-mission contrasts + age strata) and **why OSD-513 is held separate** (recurrence tested on data not used to define the pattern).
- §2.2: gene-set **curation provenance** — literature-curated panels of canonical members (not imported ontology terms), symbols resolved against an Ensembl-backed id map. (Sourced from `config/gene_sets.yaml`.)
- §2.3: acknowledged the **2,000-bootstrap vs 5,000-permutation asymmetry** (finer tail resolution for small p); added a cosine-interpretation note (read against nulls / overlapping CIs, not a fixed cutoff).
- §2.4: **TMT missingness quantified** — "negligible (median and maximum per-channel missing fraction 0.0)" from the emitted QC file `osd462_tmt_missingness_by_condition.tsv`; replaced the defensive "not limma" with the actual rationale (2 plexes, few channels → little shrinkage benefit).
- §2.5: reframed marker checks as a **post-hoc** sanity check that didn't feed back into the analysis; **justified the one-sided Fisher** (directional hypothesis; applied symmetrically to DCT1 and DCT2 bins).
- §2.6: named the **biological hypothesis behind each M0–M5 rung**.
- §2.9: tagged each perturbation branch with its **pre-specified role** (mechanism triage / robustness / spatial prediction / secondary + future-data).
- §2.11: added the explicit reason **no global FDR** is applied (non-independent re-tests; a global correction would penalize robustness checks).
- Table 1: composition family scope now specified (M0–M5 × decile/quartile). Table 2: enrichment unit clarified (phosphosite rows or parent genes).

## Results
- §3.1 now **leads with the biology**; dropped the awkward "canonical 2,000-bootstrap run"; "distinguished" not "separated the two halves."
- §3.2: added the key interpretive sentence where flat NCC is shown (**NCC is activated by phosphorylation, so flat protein ≠ flat function**); NCC-protein control made explicit; NRF2 logic cleaned up.
- §3.3: folded the repeated bulk-RNA sentence into the result; **"underpowered" → "not significant, consistent with limited power at n=20."**
- §3.4: added a plain-language handle for "phosphosite-row level"; sharpened the DCT1→DCT2 transition (DCT1 attenuates under matching, so the robust DCT2 result is informative).
- **Moved the TMT channel-centering subsection from §3.4 into §3.5** (composition sensitivities), as you flagged.
- §3.6: added what would actually distinguish the remodeling-linked vs parallel forks (spatial / time-course).
- §3.7 (low-K): now states plainly that **none of the three cohorts met the pre-specified promotion threshold** (cosine < −0.3, CI excluding zero) → hypothesis, not primary claim.
- §3.8: the LIONESS/node2vec null is now interpreted as a **useful negative control**.

## Discussion / Conclusion
- §4.2 reframed from "adds context to [16]" to **interpretive** ("constrains how the endpoint should be read"), citing the mismatch + the unbiased kinome-wide WNK corroboration.
- Consolidated the duplicated "direct test would require…" sentence (kept in §4.3, removed from §3.5).
- §4.5: added the weak external prior from the IRI repair reference (DCT-transport drop concentrated in DCT-adjacent spots → leans remodeling-linked, hedged).
- §4.6: added a synthesis sentence (the context analyses mainly **exclude** alternatives, narrowing the model).
- §4.7: added why the CMap screen is retained at all.
- §4.8: **prioritized** the future experiments (spatial/DCT-enriched phosphoproteomics is decisive; KLHL3/ubiquitinomics is lower because it probes one mechanism).
- Conclusion **rewritten to lead with the punchline**, and the WNK1/WNK3-vs-SPAK/OSR1 asymmetry now gets its own interpreted sentence.

## Deliberately left (with reason)
- **Limitations (§4.7) length / overlap with Results sensitivities** — kept thorough; this is the canonical limitations section and trimming honest caveats works against the rigor you want. Optional consolidation if a journal pushes on length.
- **Terminology sweep** ("DCT-subtype-prior" ≈ "distal-nephron subtype-prior", 9 vs 10 uses) — not blanket-replaced because the variants are partly meaning-bearing (the DCT1/DCT2 axis vs the broader distal-nephron frame). Recommend standardizing by hand on "distal-nephron subtype-prior" for the headline and "DCT1/DCT2 prior" only for the reference itself.
- **Compound-noun handles** — left as-is to avoid churn; a defined short handle on first use (e.g. "the matched permutation") could be introduced if you want.

## Needs your input (I would not invent these)
1. **Age-by-flight stability metric.** No clean quantitative value in the repo, so the text now reads "did not hold a consistent sign under bootstrap resampling." If you have a bootstrap direction-agreement fraction or Jaccard, drop it in.
2. **"Figure 3D" network-candidate panel.** I added the network-null interpretation to §3.8 text, but couldn't verify which figure panel you meant — confirm the panel reference if you want it cited inline.
3. **Compendium dependency** (12+ "companion results compendium" pointers). Fine for a preprint; for journal submission these likely need to become numbered supplementary items or be described in-text. Strategy call for you.
