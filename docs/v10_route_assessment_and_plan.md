# v10 Route — Honest Assessment and Implementation Plan

**Status:** decision document · May 2026
**Companion docs:** `docs/statistical_assessment_v5.md`, `docs/osd462_multiomics_analysis_plan.md`

---

## 1. Where the project actually stands

The analysis is rigorous, verified, and honestly reported. That is real and it is not in question. What is in question is *novelty*, and the situation is fixed by facts that no further computation changes:

- The phenotype — spaceflight NCC/DCT transporter dephosphorylation and DCT remodeling — is established by the Cosmic Kidney Disease study (Siew et al. 2024), using the **same RR-10 / OSD-462 data** this project anchors to.
- v9 is a "state map around an established phenotype." Your own second-opinion AI rated it ~4.5/10 novelty. That is roughly right.
- Your v10 plan adds regulator-activity layers. The same AI rated the *maximal computational-only* version at ~6–6.5/10, and was, if anything, generous.

So the honest baseline: this is a sound project whose novelty ceiling, by computation alone, is moderate. That is a property of the datasets and of the fact that Siew got there first — not a property of your effort or skill.

## 2. The pattern worth naming

This project has gone v4 → v9, and each round ends the same way: a new layer, brief satisfaction, then the same dissatisfaction returns. v10 is the next iteration of that loop. It is worth being honest about why the loop does not close: the thing being chased is *discovery-level novelty*, and discovery-level novelty is not present in these datasets to be found. Every layer so far has been competent, and every layer has left the core unchanged. v10 will do the same — it will produce a better paper and it will not produce the feeling you are looking for.

That is not an argument against finishing. It is an argument for deciding, deliberately, what "finished" means — before writing more code.

## 3. What v10 (as you drafted it) will and will not do

Your v10 plan is genuinely well-disciplined. The wording control ("candidate upstream organizers," "prioritization not discovery," keeping null axes, evidence grades) is good scientific hygiene, and the instinct that **activity-first beats network-topology** is correct — that is real growth in the analysis thinking. KSEA on the OSD-462 phosphoproteome is the right call: NCC biology *is* phosphorylation biology.

But be clear-eyed about two things:

1. **v10 does not escape the reproducibility framing — it decorates it.** Your own v10 central claim begins "The *established* spaceflight DCT/NCC phenotype is embedded in...". "Established" = not yours. You said you do not want a paper that is "about reproducibility of another paper's findings"; v10 is still that paper, with a regulator-activity context map added around it.
2. **The strongest v10 layer largely re-derives Siew.** KSEA on OSD-462 phospho will return "WNK/SPAK activity down" as its top hit — that is the NCC dephosphorylation finding restated as a kinase score. That is a good *positive control*, but it is not new. Whatever else KSEA surfaces is exploratory, hypothesis-level nomination.

## 4. Recommended route: decide the destination, do one tight expansion, then stop

### Step 0 — Decide the destination (before any code)

There are two honest destinations. Pick one now.

- **Destination A — a finite, submittable, modest-novelty paper.** Space-biology / methods venues exist for exactly this: *npj Microgravity*, *Life Sciences in Space Research*, *Frontiers in Physiology*, *PLOS ONE*, *GigaScience*. A ~5–6/10-novelty cross-cohort + regulator-activity + honest-negatives paper is publishable there. This is reachable computationally, from here, within a bounded effort.
- **Destination B — high novelty.** This is only reachable with data you do not have: spatial transcriptomics or single-cell kidney, or histology/wet-lab perturbation. That is not a script. It is a mentor/collaborator conversation. No v10, v11, or v12 computational layer substitutes for it.

This plan implements **Destination A**. If your honest answer is B, stop here and have the collaborator conversation instead — that is the higher-value action and it is not something to delay behind more analysis.

### Step 1 — One final analytical expansion (a trimmed v10)

Two layers, not six. Each is bounded, low-risk, and standard.

**Layer A — Kinase-activity inference on OSD-462 phosphoproteomics.** KSEA (or decoupleR with a kinase-substrate prior; PTM-SEA as a cross-check). Positive control: WNK/SPAK/OSR1 must come back suppressed. Report other coherent kinase axes (SGK1, MAPK/ERK, p38/JNK, FAK/SRC, IKK/NF-κB-related) as *exploratory nominations*. One output table: `osd462_kinase_activity_summary.tsv`. This is the highest-value new layer.

**Layer B — TF / pathway activity inference on RNA.** decoupleR with PROGENy (pathway activity) and DoRothEA (TF activity), mouse versions, run *within each cohort* (RRRM-2 ISS-T, RRRM-2 LAR, OSD-513, OSD-253, OSD-462 RNA) — never by raw cross-study pooling. Report cross-cohort recurrence of TGF-β/SMAD, NF-κB, JAK-STAT, MAPK, hypoxia, NRF2, TEAD/YAP. One output table: `rna_regulator_activity_summary.tsv`.

**Optional Layer C — MOFA+ on OSD-462 only.** OSD-462 is the one cohort with RNA + protein + phospho on the same mice. A single MOFA+ model asking whether one latent factor co-loads DCT-RNA-down, NCC/SPAK-phospho-down, and ECM-RNA. Contained, genuinely interesting, but still within the RR-10 cohort. Do it only if Layers A/B finish cleanly and you have appetite.

### Step 2 — Integration

One mechanism-context table and one figure. The table ranks candidate axes (WNK/SPAK/NCC; S1P/adhesion/integrin; TGF-β/ECM; TLR/innate/macrophage; stress/preservation; hypoxia/NRF2) with columns for RNA recurrence, OSD-462 phospho support, OSD-253 context, and an explicit evidence grade. **Keep the null boundaries in the figure** — network-candidate translation failed, protein-abundance concordance failed. The negative results are part of the product and are the most genuinely useful thing the project produced.

### Step 3 — v10 manuscript, honest framing

Title and claims must not say "discovery" or "mechanism." The honest claim is the one you already drafted: *a regulator-activity context map around an established phenotype, with explicit null boundaries.* Cite Siew throughout as the source of the phenotype. Submit to a Destination-A venue.

### Step 4 — The stop rule (write this down and mean it)

**There is no v11.** "Done" is defined now: after Layers A and B (and optionally C), the integration figure, and the v10 manuscript, the project is submitted. If the post-v10 dissatisfaction returns — and based on the v4→v9 history it will — that feeling is not a signal to add another layer. It is the signal that the project is complete and the remaining gap is novelty that this data cannot provide.

## 5. What I would cut from your v10, and why

- **Perturbation-signature matching (LINCS/CMap) — cut, or demote to a clearly-labeled exploratory supplement.** LINCS is human cell lines; you have a mouse whole-kidney DCT signature. Cross-species + cross-tissue + cross-context matching yields generic hits ("TGF-β stimulation," "stress," HDAC inhibitors) that look like mechanism but are not. For a paper whose strength is rigor, this layer adds a soft, fishable section that *undercuts* the rigor. High effort, low defensible yield.
- **Histology / source-data "audit" — cut as an analysis layer.** You have no tissue and no raw images. Auditing Siew's published figures is literature review, not analysis. Keep it as one honest Discussion paragraph citing Siew's NCC/pNCC histology — not a layer with a table.

Trimming these is not timidity; it is novelty-per-hour discipline. Six layers at ~6/10 is the same paper as two layers at ~6/10, with four extra weeks of code to write, verify, and defend.

## 6. The only route to genuine novelty

For completeness, because it is the honest answer to the question you keep really asking: the project becomes genuinely novel only with **data you do not currently have** — spatial transcriptomics or single-cell kidney (to localize the ECM/S1P/innate signal relative to remodeled DCTs), or histology / wet-lab perturbation (to test causality). Each is a collaboration, not a computation. Your manuscript affiliation already says "mentor feedback pending" — the highest-value next action, if novelty is non-negotiable, is that conversation, framed as: *"I have a verified cross-cohort computational analysis around the Cosmic Kidney phenotype; what tissue or single-cell resource could turn it into a discovery?"*

## 7. Concrete first actions

1. Answer Step 0 in one sentence: Destination A or B.
2. If A: implement Layer A (KSEA on OSD-462 phospho) — it is bounded and it is the one layer with real value. Then Layer B.
3. If B: stop coding; book the mentor/collaborator conversation.
4. Either way: write down the Step 4 stop rule where you will see it.

The work you have done is solid and verified. The most useful thing now is not a tenth layer — it is a decision.
