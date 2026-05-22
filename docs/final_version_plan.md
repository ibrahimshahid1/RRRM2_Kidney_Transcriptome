# Final Version — Consolidated Plan and Status

**Status:** consolidated plan + implementation status · May 2026
**Supersedes (kept for history):** `v10_route_assessment_and_plan.md`,
`regulator_activity_prioritization_plan.md`, `phenotype_anchoring_plan.md`.
**Destination:** Destination A — a finite, honest, modest-novelty paper for a
space-biology venue (npj Microgravity / Life Sciences in Space Research), with
Dr. Casaletto as mentor; usable as an STS/Regeneron project.

## 1. The frame

Phenotype anchoring is the organizing spine. The paper hangs off a **measured
endpoint** — loss of NCC/SPAK regulatory phosphorylation — and an evidence
ladder, rather than off a moved pathway score.

**Central claim.** *Using the established spaceflight DCT/NCC tubulopathy as a
phenotype anchor, we identify RNA programs that co-vary with measured NCC/SPAK
regulatory-phosphorylation suppression — animal-matched within RR-10 — and with
known DCT remodeling across public mouse kidney spaceflight datasets.*

**Evidence ladder.** RNA state recurrence → NCC phospho-activity (measured,
Tier 1) → DCT morphology (Siew 2024, Tier 2 literature anchor) → physiology
(Tier 3, not available; named as future work).

## 2. Components and status

| Component | Role | Status |
|---|---|---|
| v9 cross-cohort RNA state recurrence | RNA-state spine | done (verified earlier) |
| OSD-462 multi-omic anchor (v9 layers 1–4) | protein null + phospho activity | done, verified |
| **Layer A — KSEA kinase activity** | quantifies the Tier-1 phenotype | **done, verified** — WNK/SPAK positive control passes (z = −6.31 / −4.12) |
| **Layer B — TF/pathway activity** | supporting regulator nomination | **code done, unit-tested; one network step runs on your machine** |
| **Phenotype-anchor matched comparison** | animal-matched RNA ↔ activity | **done, verified — see §3** |
| Integration figure + table | the A–E phenotype-anchor figure | to build |
| v-final manuscript | honest write-up | to write |

Code: `src/multiomics/regulator_activity.py`, `src/multiomics/phenotype_anchor.py`,
`scripts/regulator_activity/`, `tests/test_regulator_activity.py` +
`tests/test_phenotype_anchor.py` (15 tests pass).

## 3. Phenotype-anchoring results (run `run_20260522_phenotype_anchor`)

Animal-matched within RR-10 (20 FL+GC animals; RNA and phospho on the same
kidneys). RNA score = mean-z of the DCT/NCC-WNK transport gene set; NCC activity
score = mean-z of regulatory phosphosites (Slc12a3 T53/T58/T65/T68, SPAK
S382/S383). Both signed so higher = more DCT/NCC function.

- **Group level — concordant.** Flight lowers the DCT transport RNA program
  (−0.32) and NCC regulatory phosphorylation (−1.30). The regulatory cluster is
  suppressed ~4.3× more than non-regulatory NCC sites (−1.30 vs −0.30) —
  magnitude specificity holds here.
- **Per-animal — concordant direction, underpowered.** RNA and NCC activity
  covary positively (Spearman ≈ +0.33–0.37) but not significantly
  (p ≈ 0.11–0.16), and the condition-adjusted test is likewise n.s. — exactly
  the underpowered-at-n=20 outcome anticipated.
- **Slc12a3-removed control — passes.** The link persists when Slc12a3
  transcript is dropped from the RNA score (Spearman +0.30); it is not one gene
  covarying with its own phospho.
- **Non-regulatory-site control — honest wrinkle.** Non-regulatory NCC
  phosphosites also correlate per-animal with the RNA score (+0.55, p = 0.011) —
  the per-animal covariation is **not** exclusive to regulatory sites (plausibly
  all NCC sites partly track per-animal DCT-cell content / NCC abundance). The
  specificity evidence is therefore the **group-level magnitude** (regulatory
  suppressed ~4× more), not the per-animal correlation. Report this openly.

**Honest claim ceiling:** group-level concordance between the DCT transport RNA
program and measured NCC regulatory phosphorylation is supported; an
animal-level predictive link is suggestive but underpowered. This is a
phenotype-anchored, suggestive, defensible result — and labelled as such.

## 4. What is left

1. Run Layer B on a network-connected machine (warm the decoupler cache, then
   `run_regulator_activity.py --rna-effects ...`). Suggestive TF/pathway
   candidate lists are an acceptable deliverable.
2. Build the A–E phenotype-anchor figure.
3. Write the final manuscript on the phenotype-anchor frame; cite Siew 2024 as
   the source of the phenotype; keep all four hypothesis outcomes and both
   controls (including the non-regulatory-site wrinkle) visible.
4. **Stop rule:** this is the last version. After the figure and manuscript, the
   project is submitted.

## 5. Honest scope (carry into the manuscript)

- The NCC dephosphorylation phenotype is established (Siew 2024, same RR-10
  data); this work reproduces and contextualizes it, and does not discover it.
- The Tier-1 phenotype is a *molecular activity* endpoint (mass spec), one rung
  below the imaging anchor of the retinal paper; Tier 2 is literature; Tier 3 is
  absent. State the ladder explicitly.
- Suggestive results and ranked candidate lists are intended deliverables;
  the requirement is honest labelling of confidence, not certainty.
