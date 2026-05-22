# RRRM-2 Kidney Transcriptome Project Closeout Packet

## Purpose

This document gives a concise, honest way to explain where the project landed after the original OSD-771 proposal, the later network analyses, and the cross-OSDR rescue attempts. It is written for mentor meetings, poster/presentation framing, and deciding whether to close or archive the project.

The point is not to pretend the project became the original proposal. It did not. The point is to preserve what was actually learned and stop the project from becoming an undefined failure.

---

## One-Sentence Finding

Public bulk mouse kidney spaceflight transcriptomes can recover a coarse terminal-flight renal remodeling signature, but they do not support a confident DCT-specific age-dependent network-rewiring mechanism with the available sample sizes, metadata, and validation layers.

---

## What The Original Proposal Claimed

The original proposal asked whether OSD-771 could reveal how spaceflight and aging reshape distal convoluted tubule biology, especially the NCC/WNK ion-transport hub, using condition-specific co-expression graphs, LIONESS/sample-specific networks, node2vec embeddings, and hub-shift scores.

The hoped-for result was a clean mechanistic claim:

> Spaceflight and aging rewire DCT/NCC/WNK-centered kidney networks, and hub-shift genes reveal molecular drivers of nephron remodeling relevant to calcium homeostasis and kidney stone risk.

That claim is not supported strongly enough by the final analyses.

---

## What Was Actually Done

The project did not stop at one failed analysis. It stress-tested the original idea across several analytical resolutions:

- OSD-771/RRRM-2 preprocessing, residualized expression, pathway scoring, and age/arm/flight contrasts.
- DCT/NCC/WNK-focused marker and pathway scoring.
- WGCNA module analysis for co-expression structure.
- LIONESS-derived sample-specific network features.
- node2vec/Procrustes embedding displacement and candidate hub-shift ranking.
- Silent-shifter and local-network candidate analyses.
- Leakage-safe comparisons of network features against expression-only baselines.
- Cross-OSDR recurrence/context analyses using OSD-513, OSD-253, and earlier explored OSDR cohorts.
- External aging-axis projection using Tabula Muris Senis kidney.
- Reframed state-space analysis of matrix, DCT transport, immune/context, preservation/stress, and mechanism scores.

The analysis also generated formal documentation:

- `latex_paper/manuscript_v6.tex` and `latex_paper/manuscript_v6.pdf`
- `latex_paper/results.tex` and `latex_paper/results.pdf`
- `latex_paper/methods.tex` and `latex_paper/methods.pdf`
- `docs/statistical_assessment_v5.md`
- `docs/biology_first_framing.md`

---

## What Held Up

The strongest recoverable biological signal is a terminal-flight kidney remodeling state.

In RRRM-2 ISS-T and directionally in OSD-513, flight is associated with:

- higher extracellular-matrix / matrix-remodeling biology;
- higher adhesion/proteolysis and related injury-repair context;
- association with innate/S1P/macrophage-context signatures;
- lower distal convoluted tubule transport signatures, especially DCT/NCC/WNK-linked biology.

The most defensible biological language is:

> Terminal-flight mouse kidney transcriptomes show a matrix-high and DCT/NCC-WNK-low tubulointerstitial remodeling state. Live animal return does not cleanly preserve that state, suggesting attenuation or redirection of matrix-associated injury-repair biology after return/reloading, while distal transport recovery remains unresolved.

This is a state-level transcriptomic result, not a direct demonstration of fibrosis, DCT dysfunction, kidney stone risk, or causal TLR4/S1P mechanisms.

---

## What Did Not Hold Up

The original DCT-specific network-rewiring mechanism did not become statistically resolvable.

Unsupported or under-supported claims:

- A confident age-dependent DCT/NCC/WNK network-rewiring mechanism.
- A robust hub-shift gene discovery claim from node2vec/LIONESS.
- A confirmed young-specific live-return reversal.
- A clean three-program LAR switch model.
- Network-derived features outperforming expression-only baselines in small-N prediction.
- Strict TLR4 requirement for the remodeling state.
- Simple accelerated kidney aging.

The negative result is not just that one method failed. The broader lesson is that the available public bulk spaceflight transcriptomic datasets are too small, heterogeneous, and metadata-limited for the fine-grained network-mechanism claim that the proposal wanted.

---

## Best Presentable Version

### Title Option 1

Bulk OSDR Kidney Transcriptomes Support Terminal-Flight Renal Remodeling Signatures But Not DCT-Specific Network-Rewiring Claims

### Title Option 2

Resolution Limits of Network-Mechanism Discovery in Public Mouse Kidney Spaceflight Transcriptomes

### Title Option 3

From DCT Hub Shifts to Renal Remodeling: A Cross-OSDR Stress Test of Mouse Kidney Spaceflight Transcriptomes

### Core Abstract

NASA OSDR mouse kidney transcriptomes provide rare access to spaceflight-associated renal molecular responses, but their ability to support fine-grained network-mechanism claims remains uncertain. We analyzed RRRM-2/OSD-771, a bulk RNA-seq cohort crossing age, mission arm, and flight/control status, using residualized expression, curated pathway scores, WGCNA modules, LIONESS sample-specific network summaries, node2vec embedding displacement, cross-cohort contrast vectors, and external OSDR context cohorts. The original hypothesis was that spaceflight and aging would induce DCT/NCC/WNK-centered network rewiring detectable through hub-shift genes. This claim was not statistically stable at primary module/pathway resolution and was not resolved by sample-specific network or embedding methods. Instead, the strongest recoverable signal was a terminal-flight renal remodeling state: extracellular-matrix, adhesion/proteolysis, and innate/S1P-associated biology increased while DCT/NCC/WNK transport signatures decreased. This matrix-high/DCT-low state was strongest in RRRM-2 ISS terminal-dissection samples and recurred directionally in OSD-513. Live animal return did not show the same pooled state shift, suggesting attenuation or redirection of the terminal-flight remodeling program rather than confirmed reversal. OSD-253 supported TLR4/innate association within a broader remodeling context but did not establish strict TLR4 dependence. These results argue that public bulk kidney spaceflight transcriptomes are useful for state-level hypothesis generation, but insufficient on their own for DCT-specific network-mechanism discovery without segment-resolved profiling, matched histology, physiology, or better-powered designs.

---

## Two-Minute Spoken Version

The project started with a stronger hypothesis than the data could support. I originally proposed that OSD-771 would let me detect age- and spaceflight-dependent DCT/NCC/WNK network rewiring using co-expression graphs, LIONESS, node2vec, and hub-shift genes.

After running the analyses, that specific mechanism did not hold up. The within-cohort aging directions were unstable, the network-derived features did not produce a robust gene-level mechanism, and the available external OSDR cohorts were too heterogeneous to validate a DCT-specific hub-shift claim.

What did hold up was a narrower transcriptomic state. RRRM-2 ISS terminal-dissection samples showed a matrix-high and DCT/NCC-WNK-low remodeling pattern, and OSD-513 pointed in a similar direction. That suggests terminal-flight kidney tissue may enter an injury-repair or tubulointerstitial remodeling state, with extracellular-matrix and adhesion/proteolysis programs rising while distal transport signatures fall. Live animal return did not preserve that pooled state shift, which suggests attenuation or redirection after return/reloading, but not a confirmed reversal.

So the final project is best framed as a boundary study: public bulk OSDR kidney transcriptomes can support coarse renal remodeling signatures, but not the fine-grained DCT network-rewiring mechanism I originally proposed.

---

## Email Draft To Dr. Casaletto

Subject: RRRM-2 kidney project update and honest reframing

Dear Dr. Casaletto,

I wanted to give you a candid update on where the RRRM-2 kidney transcriptome project landed. The original proposal focused on detecting DCT/NCC/WNK-centered age- and spaceflight-dependent network rewiring in OSD-771 using co-expression graphs, LIONESS, node2vec, and hub-shift genes.

After running the full analysis and several external OSDR checks, I do not think that original mechanistic claim is statistically supportable from the available data. The main limitation is not one failed method, but a mismatch between the fine-grained DCT network-rewiring claim and the bulk, small-N, heterogeneous public datasets available for validation.

The strongest defensible result is narrower: RRRM-2 ISS terminal-dissection samples show a matrix-high / DCT/NCC-WNK-low renal remodeling state, and OSD-513 shows directional recurrence of the same state. Live animal return does not show the same pooled remodeling shift, suggesting attenuation or redirection after return/reloading rather than confirmed reversal. OSD-253 supports TLR4/innate association within a broader remodeling context, but not strict TLR4 dependence.

I have documented the analysis as a reframed/negative study rather than forcing the original claim. My current view is that the project can be presented as a stress test of what public bulk kidney spaceflight transcriptomes can and cannot resolve: they can support state-level remodeling hypotheses, but not confident DCT-specific network-mechanism discovery without segment-resolved profiling, histology/physiology, or a better-powered design.

I would appreciate your advice on whether this is worth polishing into a cautious poster/preprint-style artifact, or whether I should close it as a documented negative result and pivot.

Best,
Ibrahim

---

## What Not To Say

Avoid:

- "We discovered the mechanism of DCT remodeling in spaceflight."
- "node2vec identified hub-shift genes driving kidney injury."
- "LAR reverses the ISS-T state."
- "This proves fibrosis."
- "This proves TLR4 dependence."
- "This proves accelerated kidney aging."
- "All methods independently validated the same story."

Use instead:

- "The original DCT network-rewiring claim was not statistically resolvable."
- "The strongest supported signal is a terminal-flight matrix-high/DCT-low remodeling state."
- "Network methods provided annotation and candidate prioritization, not primary evidence."
- "The project defines the resolution limit of the available public datasets."

---

## Poster Structure

1. **Question**
   Can public mouse kidney spaceflight transcriptomes resolve DCT/NCC/WNK network remodeling under spaceflight and aging?

2. **Approach**
   OSD-771/RRRM-2 plus external OSDR cohorts; pathway scoring, WGCNA, LIONESS, node2vec, cross-cohort contrast vectors, and mechanism/state scoring.

3. **Original Hypothesis Test**
   Age-dependent DCT network-rewiring was not stable enough to support.

4. **Recoverable Signal**
   ISS-T shows matrix/ECM/adhesion/proteolysis up and DCT/NCC/WNK transport down; OSD-513 directionally recurs.

5. **Live Return**
   LAR does not show the same pooled state shift; best read as attenuation/redirection, not confirmed reversal.

6. **Boundary Finding**
   Public bulk OSDR kidney datasets support state-level remodeling hypotheses but not fine-grained DCT hub-shift mechanisms.

7. **Future Work**
   Segment-resolved profiling, matched histology, distal-nephron NCC/WNK protein/phosphorylation assays, urine/serum chemistry, and better-powered recovery designs.

---

## If The Project Is Closed

Closing the project should mean archiving it with a documented outcome, not deleting it.

Minimum closeout artifacts:

- `latex_paper/manuscript_v6.pdf`
- `latex_paper/results.pdf`
- `latex_paper/methods.pdf`
- `docs/project_closeout_packet.md`

The closeout sentence is:

> I tested whether public OSDR kidney transcriptomes could support a DCT-centered spaceflight/aging network-rewiring mechanism. They could not. The strongest recoverable signal was a terminal-flight renal remodeling/transport-suppression state, which I documented as a reframed negative result and boundary analysis.

That is enough to close the loop without pretending the year did not happen.

