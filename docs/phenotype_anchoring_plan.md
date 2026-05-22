# Phenotype-Anchoring Reframe — Plan

**Status:** plan · May 2026 · **Companion:** `docs/regulator_activity_prioritization_plan.md`, `docs/v10_route_assessment_and_plan.md`

## Idea

Reframe the paper around a **measured phenotype endpoint** rather than a moved
pathway score — the way the retinal paper predicted measured 4-HNE / TUNEL.
The endpoint here is **loss of NCC/SPAK regulatory phosphorylation** (a measured
transporter-activity readout), with DCT remodeling as a second, established
anchor. The question becomes: *does the DCT/NCC-low RNA state align with a
measured NCC-activity phenotype, animal by animal?*

## Honest tiering of "phenotype" — state this explicitly in the paper

| Tier | Endpoint | Status | What it is |
|---|---|---|---|
| 1 | NCC/SPAK regulatory phosphorylation | **available, animal-matched** | a *molecular activity* phenotype (mass spec) — strong, but still an omic layer of the same tissue homogenate, not imaging/histology |
| 2 | DCT remodeling (segment size, density, pNCC) | **established by Siew 2024** | a *tissue morphology* phenotype — cite as literature anchor unless numeric source data are obtained |
| 3 | Urine/serum electrolytes, stone-risk | **not available** | a *physiology* phenotype — name explicitly as the missing rung / future work |

Be candid: this is **molecular-activity-phenotype anchoring**, one rung below the
retinal paper's imaging anchor. That honesty is a strength, not a weakness.

## Feasibility (verified)

RR-10 RNA-seq and phosphoproteomics were run on the **same animals**: all 30
phospho samples (`RR10_KDN_WT_*`) match RNA samples by animal ID. FL + GC gives
~20 animals with matched RNA + phospho; +BL gives ~30. Per-sample phospho
intensities (`siteQuant_360` channels) and per-sample RNA (VST matrix) are both
available, so a sample-level comparison is constructible.

## Score definitions

**NCC activity phenotype score (per animal).** Pre-specified regulatory sites
only: Slc12a3 pThr53, pThr58, pThr65, pThr68 and SPAK/Stk39 regulatory Ser
sites; Slc12a3 pThr89 included as a sensitivity variant. Per-animal score =
mean of z-scored per-channel phospho intensities across those sites. Where total
Slc12a3 protein is quantified per animal, also compute a protein-normalized
variant (phospho minus total-protein z) and report both — total NCC protein was
shown flat, so this should change little, which is itself the point.
Non-regulatory NCC sites (Ser96/120/122/124) are kept as a **negative-control
score**, not part of the phenotype.

**RNA DCT/NCC-low state score (per animal).** Per-sample mean z of the
DCT/NCC-WNK transport gene set from the OSD-462 VST matrix (and the
matrix-minus-DCT state score as a secondary). **Circularity control:** compute a
variant with Slc12a3 transcript removed from the gene set — the RNA score must
not be carried by the one gene whose phospho is the endpoint.

## The comparison — done rigorously

Join the two scores by animal. Then, in order:

1. **Group level (primary, robust).** Both scores, FL vs GC. The defensible
   claim lives here: flight animals are lower on the RNA DCT/NCC-low state *and*
   lower on NCC regulatory phosphorylation.
2. **Sample level, all animals.** Spearman correlation of RNA-state vs
   NCC-activity score. Report it — but flag that it is inflated by the BL/FL/GC
   group means.
3. **Sample level, condition-adjusted (the real test).** Correlation of the two
   scores *after removing condition means* (within-condition / partial
   correlation). This asks whether the RNA–phospho link holds beyond the shared
   flight response. With ~7–10 animals per condition it is underpowered — report
   the estimate and say so plainly. Do **not** claim per-animal "prediction"
   unless this survives.
4. **Negative control.** Repeat with the non-regulatory NCC phosphosite score;
   it should *not* track the RNA state as strongly. And repeat the RNA score
   with Slc12a3 transcript removed; the relationship should persist.

Honest claim ceiling: "animals and groups lower in the DCT/NCC-low RNA state are
lower in measured NCC regulatory phosphorylation; the relationship is clearest
at the group level, the condition-adjusted test is underpowered." That is a
phenotype-anchored statement and it is defensible.

## Phenotype-anchor figure

One multi-panel figure: (A) RNA state recurrence across RRRM-2 / OSD-513 /
OSD-462; (B) NCC activity phenotype — regulatory pNCC/SPAK down, total NCC
protein flat, non-regulatory sites flat; (C) animal-matched RNA-state vs
NCC-activity scatter, colored by condition, with condition-adjusted fit; (D) DCT
morphology direction from Siew 2024 (literature anchor, clearly labelled); (E)
evidence ladder: RNA state → phospho activity → tissue morphology → physiology,
with rungs 1–2 filled and rung 3 marked as future work.

## Reframed central claim

> Using the established spaceflight DCT/NCC tubulopathy as a phenotype anchor,
> we identify RNA programs that co-vary with measured NCC/SPAK regulatory
> phosphorylation suppression — animal-matched within RR-10 — and with known DCT
> remodeling across public mouse kidney spaceflight datasets.

"Established," "co-vary," "measured," "known DCT remodeling [cite Siew]" — all
honest. It is a genuine step up from "we define a transcriptomic context."

## What this is and is not

- **Is:** a measured-endpoint reframe; an animal-matched RNA↔activity linkage
  that no prior analysis here had; an honest evidence ladder.
- **Is not:** new biology (the NCC phenotype is Siew's), and not a
  tissue/physiology anchor (rungs 2–3 are literature / absent). The phospho
  endpoint and the RNA state are two omic layers of the same homogenate — say so.
- **No new data required**; this reorganizes existing OSD-462 data around a
  phenotype spine. If Siew numeric source data become available, Tier 2 upgrades
  from literature anchor to quantitative anchor.

## Implementation sketch

Add to `src/multiomics/`: per-animal score construction from `siteQuant_360`
channels and the OSD-462 VST matrix; the matched join; group / all-sample /
condition-adjusted comparisons; negative-control variants. One script under
`scripts/regulator_activity/`, unit tests for the scoring and the
condition-adjustment, a manifest, and the phenotype-anchor figure. Estimated:
one focused module + script, comparable in size to the OSD-462 anchor layer.
