# Preregistration: podocyte-set sensitivities (G1) and recovery persistence (G2)

**Locked:** 2026-10-06, by the commit that adds this file together with the two machine-readable
preregistrations below. Every new G1/G2 script refuses to run on real data unless its YAML is committed
and unchanged, and records the YAML's sha256 and the git commit in its manifest.

| File | Governs |
|---|---|
| `config/clinical_axes_podocyte_sensitivities.yaml` | G1 items 1–8 (item 9 deferred) |
| `config/clinical_axes_recovery_persistence.yaml` | G2 recovery/persistence and its margin sweep |
| `config/gene_map_reconstruction.yaml` | the reconstructed gene-ID baseline both depend on |

This document restates the rules for readers. If it and the YAML disagree, the YAML governs.

## Why these analyses, and their status

The 157-gene podocyte high-specificity (HS) result (pooled g = 0.689, mHK 95% CI 0.042–1.336, max-T
FWER 0.0189, 5 terminal missions, 83 animals) was the only one of 49 compartment sets to clear the
family. K0 then showed it is not separable from broad structural drift: the podocyte coefficient
attenuates 0.198 → 0.112 score units (HC3 CI −0.031 to 0.255) after adjusting for the disjoint structural
score. The project's working description is therefore "shared structural drift with a podocyte-leaning
tail", reported as a bounded observation.

- **G1** asks whether that number is robust. The sensitivities the paper blueprint listed as complete
  had actually been run on the 4 clinical axes and the 6-gene barrier core, not on the 157-gene set.
- **G2** asks whether the terminal shift survives live-animal return and Earth recovery.

Both are post hoc relative to the original 2026-08-11 lock and are preregistered here, before any of
their new quantities was computed.

## Gene-ID map baseline (decided before any new computation)

The ~14k id_map behind the 2026-08-11 numbers was not tracked and is lost. A replacement was chosen by a
rule fixed in advance: smallest maximum absolute deviation D from 27 already-published target values,
requiring the published 49-set family.

| Map | D | Evaluable sets | Role |
|---|---|---|---|
| Ensembl 116 minus 9 non-coding symbols (`id_map_reconstructed.tsv`) | 0.00044 | 49 | **baseline** |
| raw Ensembl 116 (`id_map.tsv`) | 0.05238 | 49 | runner-up: reconstruction sensitivity |
| protein-coding only | 0.00243 | 48 | disqualified (family changes) |
| OSDR annotation only | 0.01329 | 48 | disqualified (family changes) |

The nine exclusions were reverse-engineered from two published rows, so agreement on those rows is
partly by construction. Values not used in the fit still reproduce: K0, structural, the disjoint variants,
barrier adjustment and the gene-coherence counts. **Every verdict is computed under both the baseline and
the runner-up map. A verdict that differs between them is reported as "not robust to gene-ID
reconstruction", and the more conservative verdict governs.**

## G1 — sensitivities of the 157-gene podocyte HS set

Default gate: the estimate has the reference's sign AND estimate / reference ≥ 0.50. Reference values are
read from the Phase-0 reproduction run, never typed in.

| Item | Analysis | Gate |
|---|---|---|
| 1 | signed-median scoring across the 49-set family | median pooled estimate > 0 |
| 2 | genes eligible in all five missions | pooled estimate > 0 (equals the K4 strict-matching target, already known: 0.6926) |
| 3 | leave one mission out | all 5 estimates > 0 |
| 4 | gene influence | every leave-one-gene estimate keeps direction with ≥ 50% retention; no gene > 25% of the signed effect |
| 5 | OSD-462 mRNA and UPX preparations | default gate |
| 6 | OSD-253 white-light rerun control | default gate |
| 7 | OSD-163 mapping-rate adjustment (label-blind residualization) | default gate |
| 8 | leave one atlas source study out (8 rebuilds) | default gate for every QC-passing rebuild with ≥ 8 podocyte HS genes; < 8 is NOT_EVALUABLE |
| 9 | second independent atlas | **deferred**; see `docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md` |

The overall label is ROBUST_DIRECTION if items 1–7 all pass, FRAGILE (naming the failing items) if any
fails, and INCOMPLETE if none fails but one is not evaluable. Item 8 decides DEFINITION_ROBUST. The
structural sets and the K0 models are reported as descriptive comparators.

**K0 rule.** Nothing in G1 can upgrade K0. Its verdict rests on the HC3 interval for P1, which crosses 0.
The corrected Freedman–Lane joint permutation is the K0 primary scheme from now on, and the legacy 0.061 is
kept for audit. A corrected family-wise p below 0.05 may be reported only as "P1 is unusual within {P1, S1,
D} under a corrected null, but its interval crosses zero; still a bounded observation."

## G2 — persistence through live return and recovery

**Design.** The primary estimand is within OSD-771 (RRRM-2):

- 40 animals: the ISS-terminal and live-animal-return (LAR) arms, each with its own ground controls, 5
  per arm × age × condition cell.
- Model: `score ~ 0 + arm×age blocks + flight_ISS-T + flight_LAR`, OLS with HC3 errors.
- Inference: Freedman–Lane residual permutation within blocks, 100,000 permutations, the same
  permutation for every set.
- Primary target: podocyte HS. Key secondary: the disjoint structural set and the podocyte − structural
  difference.

The two recovery cohorts are not pooled (k = 2 gives df = 1). OSD-513 is descriptive only.

**Classification, on the margin M = θ_LAR − 0.5·θ_ISS-T:**

| Label | Rule |
|---|---|
| NOT_EVALUABLE_TECHNICAL | the LAR coverage metric is flagged by the frozen QC rule |
| NOT_EVALUABLE_NO_TERMINAL_EFFECT | the HC3 95% CI of θ_ISS-T includes 0, or its sign differs from the pooled terminal estimate |
| PERSISTS | Freedman–Lane one-sided p(M > 0) ≤ 0.05 AND the HC3 90% CI lower bound of M > 0 |
| DECAYS | Freedman–Lane one-sided p(M < 0) ≤ 0.05 AND the HC3 90% CI upper bound of M < 0 (sublabel REVERSED if the θ_LAR 95% CI upper bound < 0) |
| INCONCLUSIVE | otherwise, including any disagreement between Freedman–Lane and HC3 |

**Margin sweep (owner request).** κ ∈ {0, 0.25, 0.5, 0.75, 1} is computed and reported for every run,
whatever the primary label. Only κ = 0.5 sets the classification. κ = 0 asks whether any shift is
retained (θ_LAR > 0). κ = 1 asks whether nothing decayed, which is the same as the interaction test. The
Fieller 90% lower bound on the ratio θ_LAR/θ_ISS-T is reported as the largest margin compatible with
persistence at one-sided 5%.

**Confounds stated in every output.** Recovery is aliased with:

- shorter flight (32 versus 53–56 days);
- euthanasia on Earth versus on the ISS, and the handling that goes with it;
- re-entry and landing.

There is also winner's curse: the podocyte set was selected using OSD-771 among other missions, which
inflates θ_ISS-T and biases the ratio downward. OSD-513 animals were sampled about 1 day after landing,
with a 3-day offset for the ground controls.

**Power** (simulated with the full decision rule, 4,000 simulations × 999 permutations; SD units):

| True θ_ISS-T | True ratio ρ | P(PERSISTS) | P(DECAYS) | P(INCONCLUSIVE) | P(NOT_EVALUABLE) |
|---|---|---|---|---|---|
| 0.689 | 1 | 0.009 | 0.006 | 0.229 | 0.756 |
| 0.689 | 0 | 0.001 | 0.069 | 0.180 | 0.750 |
| 1.174 | 1 | 0.117 | 0.001 | 0.512 | 0.369 |
| 1.174 | 0.5 | 0.009 | 0.032 | 0.586 | 0.372 |
| 1.174 | 0 | 0.001 | 0.213 | 0.421 | 0.366 |
| 1.174 | −0.5 | 0.000 | 0.495 | 0.135 | 0.370 |

The minimum detectable interaction at 80% power is about 1.8 SD. **The expected outcome is
INCONCLUSIVE or NOT_EVALUABLE. Either label reflects an underpowered design and is not evidence for or
against persistence.** A non-significant interaction is never reported as "persisted".

## Known before lock (provenance disclosure)

- **Already published:** the 4-axis results and their sensitivities; the 49-set primary result and
  its mission effects; the disjoint, C3H and Podxl/Nid1 variants; K0 under the legacy permutation; the
  K4 target 0.6926 (= item 2).
- **Under the runner-up map only:** the corrected K0 joint and Freedman–Lane FWERs (P1 0.0536 / 0.0563).
- **Already known recovery effects:** the LAR and OSD-513 4-axis recovery effects, and the ISS-T
  podocyte mission effect (1.174).
- **Flight-blind only:** the atlas leave-one-study-out tier rebuilds (gene membership only; all 8 pass
  QC, Jaccard 0.79–0.97 versus the full set).
- **Not computed by anyone before lock:** any G1 item 1 or 3–8 flight effect on the 157-gene set, and any
  LAR or OSD-513 effect for a compartment set, the podocyte set, the structural sets or the K0 contrast.
