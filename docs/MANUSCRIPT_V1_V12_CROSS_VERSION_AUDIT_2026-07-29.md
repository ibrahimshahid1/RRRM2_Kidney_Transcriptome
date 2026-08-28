# Manuscript v1-v12 cross-version audit

**Date:** 2026-07-29  
**Scope:** `manuscript_draft.tex`, `manuscript_wgcna_pivot.tex`, and `manuscript_v2.tex` through `manuscript_v12.tex`, checked against the Stage 0 OSD-462 provenance audit and the locked exact v13 continuous-phosphoproteomic run.

## Executive conclusion

The repository does not contain a stronger, hidden DCT2/CNT result in an earlier manuscript. Before v11, “DCT” was a curated bulk-RNA transport score; it was not a subtype-resolved phosphorylation endpoint. V11 introduced the genuinely new question, but v12 promoted essentially the same thresholded, low-DCT1 enrichment result into a stronger headline while simultaneously documenting that the result was score-dependent and null in the adjusted parent-gene models.

The revision history is therefore a sequence of scientifically useful falsifications and rhetorical pivots, not twelve independent discoveries:

1. a cohort-confounded PPAR/sex claim;
2. a credible ISS-T-associated, preserved ECM/cell-migration RNA module;
3. an unstable age-rewiring hypothesis replaced by a matrix-high/DCT-low cross-cohort RNA direction;
4. an overcalled live-return reversal corrected to attenuation;
5. a multi-omic architecture built around an already published RR-10 NCC result;
6. a thresholded subtype-prior analysis promoted beyond what its inferential unit and score sensitivity supported.

The current repository boundary is stricter still. Stage 0 found no isolated, assay-qualified canonical NCC/SPAK phosphoform, and flight status is perfectly aliased with reporter-tag blocks. In the exact v13 run, the DCT2/CNT-transition set is non-evaluable because only 5 genes are phosphoproteomically observable (frozen minimum 8). ASDN has 16 observable genes and a positive competitive shift after canonical exclusions (empirical/maxT \(p=0.0291\)), but it fails the prespecified selectivity gate because a podocyte annotation has a larger statistic; no direct ASDN-versus-podocyte contrast was tested. The corrected machine-readable tier is therefore `neither`: the DCT2/CNT component remains explicitly non-evaluable, while the ASDN component fails. Neither a DCT2/CNT title nor an ASDN “beyond-axis” claim is permitted.

## Corpus identity and audit method

`latex_paper/main.tex` is an incomplete journal scaffold with placeholder text, not the first scientific manuscript. The substantive lineage begins with `latex_paper/manuscript_draft.tex`. For every revision, the title, abstract, declared primary endpoint, Methods estimand, Results, and Discussion claim were checked directly. The final disposition also uses:

- `docs/OSD462_STAGE0_PROVENANCE_AUDIT_2026-07-28.md`;
- `docs/V13_MANUSCRIPT_DECISION_MEMO_STAGE0_2026-07-29.md`;
- `data/results/run_20260729_v13_continuous_phospho_exact_final/claim_tier.tsv`;
- `data/results/run_20260729_v13_continuous_phospho_exact_final/claim_gates.tsv`;
- `data/results/run_20260729_v13_continuous_phospho_exact_final/set_level_permutation_inference.tsv`.

## Verified manuscript lineage

| Revision | Headline and visible primary endpoint | Inferential unit | What actually changed, and what is not defensible |
|---|---|---|---|
| **Original draft** | “Sex-dependent PPAR signaling…”; male OSD-513 PPAR enrichment, female OSD-102 null, and a comparison of LIONESS/network classification with expression. | Per-gene expression/edge statistics, gene-set enrichment, and mouse-held-out classification. | The sex claim is not identified: sex is confounded with mission, strain, duration, processing, and cohort, and no sex-by-flight interaction is estimated. Calling male \(q=0.079\) replication “strong” is inflated. The useful result is the negative method comparison: sample-specific network features did not outperform expression. |
| **WGCNA pivot** | “Mission-arm-dependent activation of a preserved ECM/cell-migration module…”; grey60/ME17 becomes the main result. | Module eigengene per mouse under the factorial model; GC-reference network preservation; per-animal external module scores. | This is the strongest distinct positive result in the early manuscripts: a 48-gene ECM/cell-migration module shifts in ISS-T young (\(q=6.1\times10^{-5}\)) and old (\(q=0.0108\)), not LAR; ECM enrichment OR 29.8, \(q=0.0138\); GC-reference preservation \(Z_{\rm summary}=20.6\). “Activation” should be “eigengene shift,” and ISS-T is entangled with recovery/dissection context. External corrected support occurs only in OSD-513, not universally. |
| **v2** | “ISS-T-specific activity shifts…”; expands WGCNA, LIONESS, node2vec, direct differential correlation, hubs, and external projections. | Same mouse-level eigengene analysis plus gene/edge/network summaries. | Correctly changes the interpretation from topology rewiring to an activity shift in a preserved module. It then overextends hub membership into TLR4/TGF-\(\beta\) mechanism and overreads a single \(n=10\) LAR-young classifier comparison. No new endpoint improves on grey60. |
| **v3** | “ISS-T-aligned extracellular-matrix and tubular-transport remodeling…”; abandons age rewiring for cross-cohort pathway-vector recurrence. | Nine curated pathway-effect dimensions; cosine between cohort contrast vectors; bootstrap resampling of the external cohort while holding RRRM-2 fixed. | The rejection of the unstable age-rewiring vector is scientifically important. The ISS-T/OSD-513 cosine of 0.641 is driven mainly by matrix/fibrosis, with a smaller DCT component. The bootstrap omits RRRM-2 reference uncertainty, the dimensions are correlated and curated, and “recurrence” is stronger than the design warrants. |
| **v4** | “ISS-T, but not LAR, aligns…”; adds a matrix-high/DCT-low sample score and mechanism axes. | Contrast-vector cosine plus a per-mouse derived state score. | The manuscript explicitly reports an OSD-513 fixed-reference label-permutation \(p=0.127\) yet calls the alignment claimable. Worse, the state axis is selected from pathways concordant in RRRM-2 and OSD-513 and then tested back in those animals; its later q-values are not independent validation. Correlations with constituent mechanism scores are partly part-whole correlations. |
| **v5** | “Terminal-flight and live-return…diverge…”; promotes young-LAR reversal and three mechanism programs. | Low-dimensional state/mechanism/module vectors, cosines, and per-mouse component scores. | The young-LAR cosine near \(-0.993\) is unstable because the LAR effect magnitude is near zero, and many hand-defined axes expand analytic flexibility. This is the clearest example of a directional statistic being mistaken for evidence of a resolved biological response. |
| **v6** | “Cross-cohort contrast vectors identify…”; revises reversal to attenuation/divergence and removes a redundant derived component from vector geometry. | Same vector and mouse-score units, with a cleaner estimand. | This is a healthy correction: pooled LAR is near zero; young-LAR \(p=0.361,\ q=0.542\), and its vector interval crosses zero. Mechanism divergence remains exploratory. The OSD-513 “recurrence” and outcome-trained state score remain statistically weak. |
| **v7** | “Cross-cohort transcriptomics and OSD-462 multi-omics refine…”; recurrent RNA, discordant/flat protein, and decreased NCC/SPAK phosphosite features. | Cohort-level RNA pathway vector; gene/set-level protein summaries; phosphosite-row OLS with plex covariate. | This introduces the multi-omic architecture, not a new biological endpoint: the NCC phenotype was already reported from RR-10. The version treats integer-position features as regulatory sites and calls OSD-462 independent even though it is the same biological experiment. Stage 0 later invalidates the isolated canonical-site interpretation. |
| **v8** | “A matrix-high/DCT-low direction…aligns with the published Cosmic Kidney Disease NCC phenotype”; explicitly presents a reproducibility rather than discovery paper. | Same RNA, protein, and phosphosite units, organized around four internally declared hypotheses. | The four decisions are a valuable audit record: ECM RNA-to-protein concordance falsified; DCT RNA-to-protein concordance falsified/total NCC flat; phospho-axis supported under the then-used feature mapping; network-candidate translation falsified. “Orthogonal” or “independent reproduction” is nevertheless wrong because the same OSD-462 measurements and same RR-10 animals underlie the analysis. The declaration is manuscript prose, not independently verified preregistration. |
| **v9** | “A recurrent transcriptomic context for distal tubule transporter suppression…”; treats NCC as a positive-control anchor. | Same multi-omic units, with emphasis on context mapping and null boundaries. | This is the best-calibrated pre-subtype framing: the contribution is RNA context, not rediscovery of NCC. It is still undermined retrospectively by Stage 0 residue/provenance findings and cannot serve as a fallback biological paper. |
| **v10** | “Cross-omic decoupling…”; adds matched-animal correlations, composition scores, KSEA, and regulator activity. | Animal-level correlations (\(n=20\)); phosphosite/substrate rows; bulk-RNA compartment scores. | Near-zero RNA/protein agreement is evidence of little detected concordance, not biological “decoupling.” The striking endothelial-phospho correlation is mainly flight/control separation; within-group correlations are weak. The intended non-regulatory phosphosite control correlates at least as strongly in some analyses. KSEA reuses the same phosphosite effects and is neither independent nor orthogonal; kinase-specific causal language is unsupported. |
| **v11** | “Layer-resolved…implicate distal-nephron phosphoregulatory suppression”; first asks whether suppressed sites map to DCT subtype extremes. | Primary Fisher test treats nominally suppressed phosphosite rows as observations; parent-gene Fisher/logistic/bootstrap/matched permutation are sensitivities. | This is the first genuinely new question, but the primary unit is wrong for the headline. “DCT2-leaning” means low DCT1-minus-DCT2 mean expression, not positive DCT2 specificity or cell of origin. The DCT2 matched permutation is already present (\(p\approx0.017\) under the abundance score), while the joint adjusted logistic result is not significant (\(q\approx0.104\)). |
| **v12** | “Spaceflight suppresses distal-nephron phosphoregulation beyond the canonical NCC/DCT1 axis”; the v11 subtype result becomes the title, abstract, model, and conclusion. | Same nominal-\(p<0.05\) site-row Fisher endpoint remains visually primary; gene-level and threshold-free methods remain secondary. | This is mainly a claim-level revision, not a new analysis stage. The added score ladder weakens the claim: DCT2-leaning OR 1.77 by abundance difference, 1.85 detection-aware, 1.33 specificity, and 1.00 under rank averaging; specificity-matched permutation \(p=0.062\); adjusted parent-gene logistic models are null for every score. The title asserts robustness that the Results refute. |

## What is scientifically useful but missing or demoted in v12?

There is no hidden direct DCT2/CNT evidence. There are, however, useful results and audit boundaries:

1. **Grey60 is the strongest earlier positive result.** It is a compact, preserved ECM/cell-migration RNA module with an ISS-T-associated eigengene shift and coherent hubs. It is absent from the v11-v12 main story. It could support a separate, narrowly framed transcriptomic paper, provided “ISS-T-associated” replaces “flight-activated” and the collection/recovery confounding and limited external recurrence are central.
2. **The network-negative result is unusually clean.** LIONESS/node2vec and direct differential-correlation features did not consistently outperform expression or translate into protein/phosphoprotein hits. This is useful methodological evidence, but not a kidney mechanism.
3. **The rejected age and LAR claims are valuable falsifications.** Age-rewiring vectors were unstable, and the live-return reversal dissolved when effect magnitude and uncertainty were examined. These should remain in provenance or a supplement, not be revived.
4. **The matrix/ECM RNA direction is more reproducible than the DCT RNA direction.** V12 itself ultimately reports five-cohort matrix recurrence \(p=7.0\times10^{-4}\), whereas DCT/NCC-WNK is heterogeneous and non-significant (\(p=0.19,\ I^2\approx73\%\)). It can motivate the tissue context, but it cannot validate phosphorylation or DCT2/CNT attribution.
5. **The v8 four-hypothesis record is worth preserving.** Its three falsifications prevent the story from silently drifting back toward protein abundance, network candidates, or universal cross-omic propagation. The historical “phospho-axis supported” decision must now be annotated as superseded by Stage 0 rather than presented as validated.

## Current v13 boundary

Stage 0 resolves several ambiguities that v7-v12 did not:

- the detailed protocol supports **TMTpro** labeling, while a legacy repository field says iTRAQ; the correct wording is “TMTpro isobaric-tag proteomics/phosphoproteomics, with a disclosed legacy metadata inconsistency”;
- the workbooks contain 30 wild-type samples across two plexes (baseline, flight, and ground), while the primary comparison uses 10 flight versus 10 ground samples; the p21-null arm is not part of this workbook contrast;
- reporter-tag blocks and flight status are perfectly aliased in both plexes, with no label swap, so label permutation cannot separate a biological flight effect from a tag-block effect;
- the NCC T53-indexed feature is a T53/Y65 co-modified peptide, and the SPAK S383-indexed feature is an S382/S383 co-modified peptide;
- no isolated canonical NCC T53/T58/S71 or SPAK T243/S383 feature passed the frozen assay-qualification rule.

The exact v13 run replaces nominal-site selection with 8,021 label-blind phosphosite features, 3,524 fixed-null-eligible parent genes, one continuous score per gene, and all 63,504 within-plex label assignments. Its decision is:

| Frozen primary set | Exact outcome | Interpretation |
|---|---|---|
| DCT2/CNT transition | 5 observable genes; minimum is 8; gate `non_evaluable`. | The dataset cannot adjudicate the declared DCT2/CNT hypothesis. This is not evidence for a null and cannot support a subtype title. |
| ASDN | 16 observable genes; positive competitive statistic 0.718; maxT \(p=0.0291\); leave-one-gene-out and declared direction sensitivities pass. | There is an ASDN-associated shift, but it fails the prespecified selectivity criterion. The podocyte comparator has a larger statistic (1.112; empirical \(p=2.36\times10^{-4}\), comparator BH \(q=0.00213\)); fibroblast is also enriched (statistic 0.707; \(p=0.00375,\ q=0.0169\)). No direct ASDN-versus-podocyte contrast was tested. The ASDN claim gate therefore fails. These are parent-protein annotations in whole kidney, not cell localization. |
| Overall tier | `neither`; DCT2 component `non_evaluable`; ASDN component `fail`; DCT2 title false; ASDN-beyond-axis claim false. | A biological late-distal subtype paper is a no-go under the frozen decision rule. The valid output is a provenance/methods/null-boundary report, or no manuscript from this endpoint. |

## Branch disposition for the v13 outcome

Here, **Keep** means main text of a provenance/methods or null-boundary report—not permission to write a DCT2/CNT biological paper.

| Branch | Disposition | Reason |
|---|---|---|
| Exact sample/plex/tag map, TMTpro/legacy-iTRAQ reconciliation, residue and phosphoform audit | **Keep** | This is the most publication-grade new work and defines what the assay can and cannot identify. |
| Continuous one-gene-one-score, exact within-plex permutation framework and frozen claim gates | **Keep** | It fixes the primary statistical-unit and nominal-threshold problems and reports the non-evaluable outcome transparently. |
| Whole-kidney compartment comparison | **Keep** | It makes ASDN fail the prespecified late-distal selectivity criterion; no direct ASDN-versus-comparator contrast was tested. Describe annotations, not cell-of-origin effects. |
| ASDN positive but non-selective shift | **Context** | Numerically real under the declared set test, but it fails the comparator gate and cannot carry title, abstract, or conclusion. |
| Flat total NCC abundance | **Context** | Permits only “no detected whole-kidney abundance decrease”; it does not establish activation state or exclude altered turnover/degradation. |
| Matrix-high/ECM RNA recurrence | **Context** | A short motivation paragraph or separate RNA paper; not validation of phosphosite biology. |
| Grey60 preserved ECM/cell-migration module | **Supplement / separate-paper candidate** | Coherent and distinct, but unrelated to the frozen late-distal endpoint and limited by arm/collection context. |
| Network methods fail to beat expression; network candidates fail multi-omic translation | **Supplement** | Useful negative method result and audit trail, not mechanism. |
| V8 four-hypothesis decision record | **Supplement** | Preserves falsifications; the old H3 decision must be marked superseded by Stage 0. |
| KSEA/WNK motif compatibility | **Supplement at most** | Same underlying phosphosite data, overlapping motif assignments, no kinase-specific in-vivo causality, no independence. |
| Composition M0-M5 ladders and matched-animal correlations | **Supplement at most** | Small \(n\), collinearity, post-treatment covariates, reselected sets, and pooled group separation prevent causal or localization claims. |
| PPAR sex dimorphism; TLR4 dependence; LAR reversal; age mimicry; clock-DCT/S1P subdivisions | **Retire** | Confounded, unstable, post-hoc, or explicitly falsified by later revisions. |
| Canonical NCC/SPAK regulatory dephosphorylation as an OSD-462 MS result | **Retire** | Stage 0 found only co-modified features and no isolated assay-qualified canonical phosphoform. Prior antibody work on the same RR-10 material is context, not independent replication. |
| V11-v12 nominal-\(p<0.05\) low-DCT1 Fisher analysis and “beyond NCC/DCT1” title | **Retire** | Superseded by the continuous exact analysis; wrong primary unit, indirect subtype definition, score dependence, and no passing v13 claim tier. |
| Human urine, low-K, IRI spatial spots, dDAVP, CMap, aging, live return, exhaustive pathway/network branches | **Retire from main; archive as compendium if desired** | None repairs subtype attribution, assay aliasing, or phosphoform identity. Several add pseudoreplication or remote biological analogies. |

## Prose and document-form evolution

The AI-like quality is not merely vocabulary. From v4 onward, the manuscripts increasingly embed an audit memo inside the paper: each result is followed by a qualification, a defense of why the qualification does not erase it, a narrowed restatement, and a future experiment. V7-v12 add executive summaries, claim hierarchies, roadmap figures, extensive “supported claim” language, and repeated epistemic verdicts in the abstract, Results, Discussion, captions, and conclusion. The same few contrasts—RNA as context, phosphorylation as the informative layer, DCT2 as modest but robust—are restated until the caveat-management becomes more visible than the evidence.

V12 is the sharpest example: the abstract attempts to announce the extension, litigate four scoring systems and several null models, and preserve the extension anyway. A journal rewrite should not edit that architecture sentence by sentence. It should be rebuilt from the permitted endpoint outward: one result per Results sentence; interpretation in Discussion; no table of contents or executive summary; no “claim adjudication” slogans; short captions; and only analyses that change the conclusion.

## Final judgment

The strongest genuinely distinct, earlier repository result is grey60, not DCT2/CNT. The latest exact phosphoproteomic result is broader compartment-associated suppression with a suggestive but non-selective ASDN component, not a distal-nephron-specific extension. The established NCC biology remains important literature context, but this MS reanalysis cannot make it the new result or even use it as a secure isolated-site fallback.

Accordingly, the revision history supports one of two honest products:

1. a compact provenance/methods/null-boundary report centered on assay identifiability, reporter-tag aliasing, exact continuous inference, and the failed/non-evaluable subtype gates; or
2. a separate RNA manuscript centered on the ISS-T-associated preserved ECM module and matrix remodeling, with the DCT/NCC RNA direction explicitly treated as heterogeneous context.

It does not currently support a v13 biological paper claiming suppression beyond NCC/DCT1 into DCT2/CNT or the ASDN.
