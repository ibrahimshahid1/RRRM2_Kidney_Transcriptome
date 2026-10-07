# Strategy reader: Prior art and novelty

_Research reader's notes (2026-10-07). Written under a ~15 search/fetch budget. pmc.ncbi.nlm.nih.gov, ebi.ac.uk (Europe PMC), scidok, and sainsburywellcome.org were all blocked by the egress proxy, so NO full text was read. Everything below rests on search-result abstracts and snippets. Items marked UNVERIFIED must be checked before being cited or relied on. Not a decision document._

## Summary

- The only paper confirmed to bear directly on the project is Siew et al. 2024 (Nat Commun), which pooled 11 spaceflight mouse missions. Its confirmed headline findings are renal transporter dephosphorylation, DCT size expansion with loss of tubule density, and GCR-induced damage. I could not confirm what it says about glomeruli, podocytes, PODXL, kidney weight, or per-mission reproducibility. The project's podocyte claims therefore cannot yet be cleared against it.
- Grandke et al. 2026 (Nat Commun) is a multi-organ miRNA paper whose title says "ECM and developmental pathways" and "age-independent tissue adaptation". It is the closest prior art for a shared structural/ECM theme across organs (at miRNA level, not mRNA). Kidney-specific content unverified.
- The RRRM-2 design (ISS-T at 55-58 d, LAR 32 d flight plus 24 d recovery) is confirmed from the GeneLab dataset record. No RRRM-2 kidney paper was found in search.
- DCT remodeling is PRIOR (Siew 2024). Everything else is either unverified or not found; absence from a ~15-query search is weak evidence of novelty.
- No found paper frames "most reported renal responses do not recur across missions" or tests exposure-duration dependence of podocyte scores.

## Key papers

1. **Siew K, et al. 2024. Cosmic kidney disease: an integrated pan-omic, physiological and morphological study into spaceflight-induced renal dysfunction. Nat Commun 15:4923 (DOI 10.1038/s41467-024-49212-1; PMC11167060).**
   - Data: 11 spaceflight mouse and 5 human missions, 1 simulated-microgravity rat mission, 4 simulated-GCR mouse missions; transcriptomic, proteomic, epigenomic, epiproteomic, metabolomic, metagenomic, clinical chemistry, histology, 3D imaging, miRNA-ISH, tissue weights.
   - Confirmed findings: transporter dephosphorylation; DCT expansion with reduced overall tubule density; renal damage/dysfunction under Mars-roundtrip-equivalent simulated GCR; coined "cosmic kidney disease".
   - UNVERIFIED: glomerular/podocyte (PODXL) results, kidney weight, GCR-related glomerular findings, which OSD datasets, any recovery data. Must read full text (PMC blocked here; try publisher site or Crick/Salford repository PDF).
2. **Grandke F, et al. 2026. MiRNAs shape mouse age-independent tissue adaptation to spaceflight via ECM and developmental pathways. Nat Commun (DOI 10.1038/s41467-026-68737-1; PMID 41644562; published 5 Feb 2026).**
   - Data: 686 small-RNA samples, female mice, 13 solid organs, 3 and 8 months old, at least 3 weeks on ISS vs ground controls.
   - Findings (abstract-level): spaceflight effects in systemic tissue-remodeling pathways along the fat-liver-pancreas axis and in heart, brain, spleen, thymus; MIR-17/92 and MIR-1/133 families. Whether kidney (and OSD-913) is in the 13 organs, and whether it reports ECM/cytoskeleton at mRNA level, is UNVERIFIED. Search did not return "OSD-913" or "RRRM" strings.
3. **Suzuki T, et al. 2022. Gene expression changes related to bone mineralization, blood pressure and lipid metabolism in mouse kidneys after space travel. Kidney Int.** (MHU-3, Nrf2 WT/KO, ~31 d ISS; DOI not retrieved.) Only the title and the Nrf2-homeostasis framing were seen. Findings on podocyte/glomerular genes, DCT, or recovery are UNVERIFIED. Related MHU-3 papers: Nrf2 and metabolic response during and after spaceflight (Commun Biol 2021, PMC8660801); Nrf2 alleviates spaceflight-induced immunosuppression and thrombotic microangiopathy (Commun Biol 2023, PMC10457343). The "during and after" wording suggests some post-flight data exist in MHU-3; check.
4. **Hindlimb unloading kidney literature.** Search surfaced: (a) Alteration of vasopressin-aquaporin system in hindlimb unloading mice, Front Physiol 2025 (DOI 10.3389/fphys.2025.1535053): copeptin and AQP2 up at 1 week, AQP2 mRNA/protein and copeptin down at 4 weeks (a duration-dependent reversal, but for the collecting duct, not podocytes); (b) a snippet that HU mouse kidneys show Bowman's-capsule widening and increased glomerular area with ER-stress involvement (source paper not identified; possibly the PMID 41946794 / GSE235042 paper, but this was NOT confirmed; PMID lookup returned nothing relevant). Whether that paper reports podocyte transcript changes is UNVERIFIED.
5. **RRRM-2 dataset record (GeneLab/data.gov).** 40 female C57BL/6NTac mice, 12 wk and 29 wk at launch; ISS-T sacrifice at 55-58 d; LAR returned live at 32 d and sacrificed after 24 d recovery; ground, vivarium, basal controls; kidney bulk RNA-seq (ribodepleted), 5 young + 5 old per group. No associated peer-reviewed RRRM-2 kidney paper was found.
6. **Cross-mission meta-analyses.** Found: (a) NASA GeneLab-derived microarray comparison of mouse and human altered-gravity datasets (PMC11053165), emphasizing highly reproducible genes and some musculoskeletal overlap; (b) a systems-biology analysis of seven rodent datasets proposing TGF-beta1 as a common master regulator and 13 shared miRNAs (likely Beheshti et al. 2018, Sci Rep, s41598-018-22613-1; UNVERIFIED which organs/kidney). Neither addresses kidney or a reproducibility-failure framing as far as snippets show.
7. **Early RR-1 cytoskeleton paper** (heart and lung, PMC5955502): cytoskeletal protein content unchanged but mRNA of cytoskeletal genes changed. Marginal precedent for cytoskeletal transcript drift without protein change.

## Verdict table

Labels are provisional given that no full text was read.

| # | Claim | Verdict | Papers / basis |
|---|---|---|---|
| i | Podocyte/glomerular transcript abundance up across missions | NOVEL (provisional) | Not seen in any snippet. Siew 2024 pooled 11 missions and may contain glomerular transcript results; must read before claiming. HU Bowman's-capsule/glomerular-area widening is morphological support, not transcript prior art |
| ii | Shared low-heterogeneity structural/ECM/cytoskeletal drift | PARTLY PRIOR | Grandke 2026 (ECM/tissue-remodeling across organs, miRNA level); RR-1 cytoskeletal mRNA change (heart/lung); Beheshti-type shared-regulator work. No kidney mRNA, cross-mission, low-heterogeneity framing found |
| iii | Glomerular barrier genes higher in flight (opposite of identity loss) | NOVEL (provisional) / possible CONFLICT | No paper seen claiming identity loss, but Siew 2024 "damage and dysfunction" narrative and HU glomerular widening could be read as opposed; check Siew full text for podocyte gene direction |
| iv | Recovery after live return | PARTLY PRIOR | MHU-3 "during and after spaceflight" papers (Commun Biol 2021/2023, Suzuki 2022) hint at post-flight data; RRRM-2 design is public but no kidney recovery paper found. Unverified |
| v | Exposure-duration dependence | PARTLY PRIOR | HU time course (1 vs 4 wk, AQP2/copeptin) shows duration dependence for a different compartment; no spaceflight podocyte duration analysis found |
| vi | "Most renal responses do not recur across missions" | NOVEL (provisional) | Not found. Closest are GeneLab cross-study reproducibility comparisons (musculoskeletal) and Siew's pooled-mission approach, which emphasises consistency rather than failure to recur |
| vii | DCT density/remodeling | PRIOR | Siew 2024 (DCT expansion, reduced tubule density). Do not claim |

## Gaps/forward-search still owed

1. Read Siew 2024 full text and supplementary: glomerular/podocyte findings, PODXL, kidney weight, GCR glomerular data, dataset list (OSD IDs), and any per-mission DE concordance. Highest priority; decides i, iii, vi.
2. Identify and read the PMID 41946794 / GSE235042 hindlimb-unloading paper.
3. Read Grandke 2026: is kidney among the 13 organs, is it OSD-913, and what ECM/cytoskeletal statements are made.
4. Find the Suzuki 2022 Kidney Int DOI and findings; check MHU-3 post-flight kidney data.
5. Search for any RRRM-1 multi-organ single-cell paper (OSD-918 series) and any RRRM-2 kidney paper or bioRxiv preprint; forward-search 2026 citations of Siew 2024.
6. Search other recovery studies (e.g. STS-era post-flight rodent kidney, RR live-return) and GeneLab AWG cross-mission kidney analyses (Beheshti, Cekanaviciute, Overbey, da Silveira 2020 Cell).
7. Retry blocked hosts via alternative mirrors (Europe PMC, publisher site) once egress allows.

## Sources

- [Siew et al. 2024, Nat Commun (PMC11167060)](https://pmc.ncbi.nlm.nih.gov/articles/PMC11167060/) (search snippet only)
- [Siew 2024 record, Salford repository](https://salford-repository.worktribe.com/output/3335990/cosmic-kidney-disease-an-integrated-pan-omic-physiological-and-morphological-study-into-spaceflight-induced-renal-dysfunction)
- [Siew 2024, Crick publications](https://www.crick.ac.uk/research/publications/cosmic-kidney-disease-an-integrated-pan-omic-physiological-and-morphological-study-into-spaceflight-induced-renal-dysfunction)
- [Grandke et al. 2026 (PMC12876965)](https://pmc.ncbi.nlm.nih.gov/articles/PMC12876965)
- [Grandke et al. 2026 PDF, Saarland repository](https://scidok.sulb.uni-saarland.de/bitstream/20.500.11880/42057/1/s41467-026-68737-1.pdf)
- [Hindlimb unloading vasopressin-aquaporin, Front Physiol 2025](https://www.frontiersin.org/journals/physiology/articles/10.3389/fphys.2025.1535053/full)
- [MHU-3 Nrf2 metabolic response, Commun Biol 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8660801)
- [MHU-3 Nrf2 immunosuppression/TMA, Commun Biol 2023](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10457343/)
- [RRRM-2 kidney transcriptional profiling dataset record](https://catalog.data.gov/dataset/transcriptional-profiling-of-kidney-tissue-from-mice-flown-on-the-rodent-research-referenc)
- [GeneLab-derived microarray meta-comparison (PMC11053165)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11053165/)
- [Seven-dataset rodent systems-biology analysis (Sci Rep 2018)](https://link.springer.com/10.1038/s41598-018-22613-1)
- [RR-1 heart/lung cytoskeleton (PMC5955502)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5955502/)
- Repo context read: docs/PRIOR_ART_VERDICT_2026-07-29.md; docs/strategy/2026-10-07/reader_history.md (top 40 lines)

## Coordinator addendum (2026-10-07): Siew et al. 2024, verified from the full text

Source: Europe PMC full-text XML, PMC11167060 (doi:10.1038/s41467-024-49212-1).

**What Siew 2024 did**
- Pooled 11 spaceflight-exposed mouse missions, 5 human missions, 1 simulated-microgravity rat study
  and 4 simulated-GCR mouse studies.
- Integrated omics by pathway enrichment and vote counting across 24 datasets. This is not a
  per-mission effect-size meta-analysis.

**What it reports, against each of our claims**

- **(ii) Structural drift: PRIOR in qualitative form.**
  - Siew states that "structural remodelling of the nephron occurs within a month of spaceflight
    exposure".
  - The enriched pathways include focal adhesion, tight junction, gap junction, actin and cell cycle.
  - Those overlap our structural scaffold control, which is defined from the KEGG adherens junction,
    ECM-receptor, focal adhesion, actin, tight junction and VSMC terms.
  - What remains novel: the calibrated cross-mission effect size, heterogeneity (I² ≈ 0),
    family-wise inference, and the finding that the shift is not separable from a podocyte-associated
    program.
- **(i)/(iii) Podocyte HS up and barrier genes up: CONFLICT in direction.**
  - Siew lists NID1 and PODXL among "the most frequently downregulated gene products across kidney
    datasets".
  - Our per-mission estimates go up: the Podxl effect is positive in the C57BL/6 and C3H cohorts.
  - The project's K1 audit already frames this as a different estimand (a vote count over datasets
    and layers, versus an animal-level standardized set score).
  - It must be presented that way, not as a correction.
- **(vii) DCT remodeling: PRIOR.** In RR-10, DCT density is lower and the retained tubules are larger.
- **Kidney weight: PRIOR.** Kidney weight relative to bodyweight rises with simulated microgravity and
  with microgravity plus GCR. Complete weight data existed only for BNL-1/2/3 and RR-23.
- **Glomerular injury: PRIOR, but in GCR only.** Microthrombi and thrombotic microangiopathy appear in
  27% of GCR-exposed females, with proteinuria in GCR vs sham. This is radiation, not standard ISS
  spaceflight.
- **(vi) Reproducibility framing: NOVEL.** Siew's design (enrichment and vote counting, mostly female
  animals) does not test per-mission recurrence. A preregistered "what recurs" analysis directly
  complements it.

**Implication for the paper.** Position it as a quantitative, preregistered cross-mission estimate
that tests which of the field's reported renal responses recur. Siew's structural-remodeling claim
largely does recur, as structural drift with a podocyte-leaning top. Transporter and identity-loss
axes do not. Our podocyte-associated *abundance* runs opposite to Siew's PODXL vote count; the
difference is in the estimand.
