# Stage 0 go/no-go — protocol eligibility and effect-size anchor

**Date:** 2026-07-29
**Run:** `data/results/run_20260729_stage0_protocol_inventory/`
**Scripts:** `scripts/stage0/protocol_inventory.py`,
`scripts/stage0/coverage_confounding_gate.py`,
`scripts/stage0/axis_effect_size_anchor.py`

**Verdict: conditional go, but not for the four-axis ranking.** The data support
a narrower two-axis contrast. Details below.

---

## 1. Transcript integrity is degraded in every cohort, and RIN does not show it

GeneLab's published QC files contain RSeQC gene-body coverage, so the metric
Lai Polo et al. (2020) identified as diagnostic is downloadable — no realignment
needed. Ratio of the 5–20 % bin to the 80–95 % bin; 1.0 is uniform.

| Cohort | n | 5′/3′ ratio | SD | Range | Mean RIN | Library |
|---|---:|---:|---:|---|---:|---|
| OSD-102 | 12 | 0.837 | 0.013 | 0.815–0.864 | 9.43 | ribo-depletion |
| OSD-253 | 112 | 0.686 | 0.179 | 0.369–1.146 | 9.33 | ribo-depletion |
| OSD-163 | 18 | 0.639 | 0.035 | 0.569–0.693 | 8.33 | ribo-depletion |
| OSD-513 | 30 | 0.533 | 0.036 | 0.459–0.590 | 8.77 | ribo-depletion |

Every cohort shows 5′ loss, spanning a 1.6× range between cohorts, while RIN sits
at 8.3–9.4 throughout. That is Lai Polo's core observation reproduced in kidney:
**RNA quality metrics look fine and coverage does not.** OSD-253 is also
internally heterogeneous (SD 0.179, two read lengths, a "Ground Control Rerun"
factor level), so it is multi-batch.

OSD-462 and OSD-771 have no local QC file and are **unevaluated**. OSD-771 is the
primary RRRM-2 dataset, so this gap must be closed before anything is locked.

## 2. One cohort is confounded at the measurement layer

Within-cohort flight-versus-ground contrast on the same coverage ratio:

| Cohort | n FLT/GC | Coverage g | p | Gate |
|---|---|---:|---:|---|
| OSD-102 | 6/6 | +0.474 | 0.374 | pass |
| OSD-163 | 6/12 | +0.363 | 0.439 | pass |
| **OSD-253** | **20/92** | **−0.556** | **0.0023** | **CONFOUNDED** |
| OSD-513 | 10/20 | +0.290 | 0.404 | pass |

OSD-253 flight samples have systematically worse 5′ coverage than their controls.
Any flight estimate from that cohort is confounded by measurement.

Note the three passes are underpowered — g of 0.29–0.47 is not small, merely
non-significant at these sample sizes. "Pass" means *no difference detected*.

## 3. Observed effect sizes — the power anchor

Per-sample mean z across axis genes, Hedges g for flight vs ground control,
DerSimonian–Laird random effects. Gene sets are the repository's existing
`config/gene_sets.yaml` panels, **not** the frozen externally-sign-learned
panels, so these are design-calibration numbers rather than the study.

| Axis | Pooled g | 95 % CI | I² | p |
|---|---:|---|---:|---:|
| **Fibrosis** | **0.798** | 0.306 – 1.291 | 9.3 % | 0.0015 |
| ECM remodeling | 0.639 | 0.104 – 1.174 | 21.9 % | 0.019 |
| DCT/NCC-WNK | −0.307 | −1.364 – 0.751 | **77.8 %** | 0.57 |

Per cohort, the distal axis runs −0.29, +0.69, +0.26, **−2.01** (OSD-513).

## 4. Leave-one-cohort-out

| Axis | Full | −OSD-102 | −OSD-163 | −OSD-253 | −OSD-513 |
|---|---:|---:|---:|---:|---:|
| Fibrosis | 0.798 | 0.905 | 0.888 | 0.803 *(CI crosses 0)* | 0.618, I²=0 |
| ECM | 0.639 | 0.752 | 0.744 | 0.656 *(crosses 0)* | 0.441 *(crosses 0)* |
| DCT/NCC-WNK | −0.307 | −0.325 | −0.622 | −0.540 | **+0.229, CI −0.270–0.728, I²=0** |

Three things follow.

**Fibrosis is robust.** It survives removal of any single cohort, and removing
OSD-513 leaves g = 0.618 with I² = 0 — positive, significant, homogeneous.

**ECM is not.** It loses significance when either OSD-253 or OSD-513 is dropped.

**The distal-transport signal is one cohort.** Remove OSD-513 and the axis
becomes g = +0.229, CI −0.270 to 0.728, I² = 0 %. Three cohorts agree on
approximately nothing; a single live-return cohort at g = −2.01 generates the
entire pooled effect and all 78 % of the heterogeneity.

## 5. The two problem cohorts are the two driving the story

OSD-513 has the worst coverage (0.533) and the most extreme effect on all three
axes. OSD-253 fails the confounding gate outright. Exactly the pattern Lai Polo
predicts, and it means the apparent distal-transport effect and the ECM effect
both rest disproportionately on the least trustworthy measurements.

## 6. Verdict

**The four-axis ranking is not viable.** Fibrosis at 0.798 sits right at the
detectable threshold for this design; ECM at 0.639 is below it and not
LOO-robust. Pairwise separation among four axes needs more power than detecting
one, so the ranking would be uninterpretable.

**A two-axis contrast is viable**, and it is sharper than the ranking anyway:

> A fibrosis/remodeling program is reproducibly engaged across mouse kidney
> spaceflight cohorts (g = 0.80, I² = 9 %, robust to leave-one-cohort-out),
> whereas the distal-nephron transport program — the dominant published claim —
> is not (g = +0.23, CI −0.27 to 0.73, I² = 0 % after removing the single
> anomalous cohort that generates the entire apparent effect).

The negative is load-bearing here because after removing OSD-513 the distal
interval is *tight and homogeneous*, which is equivalence-quality evidence
against effects larger than about |g| = 0.73 — not merely a non-significant
result.

## 7. Required before locking

1. Pull QC metrics for **OSD-462 and OSD-771** and run both gates. OSD-771 is the
   primary dataset and is currently unevaluated.
2. Decide OSD-253's disposition. It fails the coverage gate; excluding it makes
   fibrosis borderline (CI −0.030 to 1.635). That tension must be resolved by a
   rule, not by preference.
3. Replace `config/gene_sets.yaml` panels with the externally sign-learned frozen
   panels. These numbers move once panels change.
4. Recompute p-values exactly; the current ones use a normal approximation.
5. **Declare prior knowledge in the preregistration.** The fibrosis and distal
   axes were examined before this analysis and the results were known. State it.

## 8. Boundary

Design calibration, not biology. The axes are annotation panels in bulk tissue;
none of this localises anything, and every cohort carries the collection and
preservation confounds documented in §1–2.
