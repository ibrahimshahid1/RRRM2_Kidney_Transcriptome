# Podocyte-set sensitivities (G1) and recovery persistence (G2): results

**Date:** 2026-10-06
**Run:** `data/results/clinical-axes/20261006T004825Z_g1g2` (gitignored; regenerate with the CLI commands in §7)
**Preregistration:** `docs/CLINICAL_AXES_G1_G2_PREREGISTRATION_2026-10-06.md`, locked by commit `e574feb`

## Summary

1. **The podocyte high-specificity (HS) result is directionally robust (G1: ROBUST_DIRECTION).**
   - All seven preregistered sensitivities pass: median scoring, common genes, leave-one-mission,
     gene influence, OSD-462 preparation, OSD-253 control choice and OSD-163 mapping rate.
   - All eight leave-one-atlas-study-out marker rebuilds pass (item 8, DEFINITION_ROBUST).
   - No verdict changes between the reconstructed and raw Ensembl-116 gene maps.
   - Robustness is of *direction and magnitude*, not of family-wise significance. Several variants have
     intervals that cross zero. Removing the 5 largest positive contributors (of 145 contributing genes)
     moves the lower CI bound below zero.
2. **K0 stands: a shared structural drift with a podocyte-leaning tail.**
   - The corrected permutation lowers P1's family-wise p from 0.061 to 0.042, but only on the baseline
     map; it is 0.056 on the runner-up.
   - P1's HC3 interval still crosses zero, in the primary analysis and in every sensitivity variant.
   - By preregistered rule the verdict is not upgraded.
3. **Recovery persistence is INCONCLUSIVE for the podocyte set (G2, primary).**
   - Within RRRM-2 (OSD-771), the podocyte shift is clear in ISS-terminal mice (arm g = 1.17) and small
     in live-return mice after 24 days on Earth (arm g = 0.14).
   - The design cannot separate persistence from decay. The point estimate retains about 10% of the
     terminal shift, with a Fieller 90% interval for the ratio of −0.66 to 1.08.
   - This was the preregistered expected outcome, and no margin in the sweep gives a decision for the
     podocyte set.
4. **The structural scaffold shift is smaller after live return (key secondary).**
   - For the disjoint structural set the flight × arm interaction is −0.35 score units (95% CI −0.63 to
     −0.07), Freedman–Lane p = 0.009, preregistered max-T over {podocyte, structural, difference} = 0.024.
   - Its 0.5-margin label stays INCONCLUSIVE because the Holm-adjusted permutation p is 0.053.
   - Recovery is aliased with shorter flight (32 vs 53–56 d) and on-Earth euthanasia, so this is *not*
     evidence of recovery as such.

## 1. Provenance and reproduction

- **Code.** Lock commit `e574feb` (2026-10-06 00:48:10 UTC); the stages ran at `f04f84d`, which changed
  only manifest date serialization after the first attempt failed on writing its manifest.
- **Integrity.** Every new stage recorded the committed preregistration sha256s
  (G1 `71eac97c…`, G2 `93f5bdb9…`), `prereg_last_commit = e574feb`, and a clean working tree.
- **Push timing.** A GitHub access outage delayed the push of the lock commit to about 00:55 UTC. By
  then two G1 stages had run (00:50–00:53): the median and common-gene compartment families. The local
  commit order and the manifests establish that the lock came first; the public timestamp is a few
  minutes later for those two stages.
- **Inputs.**
  - OSDR data from the NASA public S3 mirror (94 files, sizes verified against the OSDR file API).
  - Mouse Kidney Atlas `mka.h5ad` from Zenodo 10.5281/zenodo.17395591, MD5 `6a2e5308…` (matches the
    project fingerprint).
  - Rebuilt tiers `frozen_compartment_tiers.tsv` sha256 `f080a0d9…`: 157 podocyte HS genes, 49
    evaluable sets.
  - Gene map: reconstructed baseline `id_map_reconstructed.tsv` (`d36462df…`); runner-up raw Ensembl
    116 (`98e622bc…`). See `config/gene_map_reconstruction.yaml`.

**Reproduction in this run (baseline map):**

| Quantity | Published | This run |
|---|---|---|
| podocyte HS g (mHK 95% CI), max-T FWER | 0.689 (0.042, 1.336), 0.0189 | 0.689 (0.042, 1.336), 0.0192 |
| podocyte mission effects | −0.117 / 1.166 / 0.676 / 0.541 / 1.174 | identical |
| structural control g, FWER | 0.549 (−0.095, 1.193), 0.128 | 0.549 (−0.095, 1.193), 0.129 |
| K0 P0 / P1 / S1 (score units, HC3) | 0.198 / 0.112 (−0.031, 0.255) / 0.031 | 0.198 / 0.112 (−0.031, 0.255) / 0.031 |
| K0 P1 max-T (legacy scheme) | 0.061 | 0.0605 |
| 4 axes | −0.716 / 0.570 / 0.311 / −0.154 | identical |
| 49-set family vs Phase-0 reference run | — | max abs difference 0 |

**Correction to earlier wording.** Under the reconstructed map the structural control ranks **2nd** by
|t|, not 4th. Ranks 2–5 lie within 0.02 of each other (structural 0.549, podocyte broad_enriched 0.545,
podocyte all_enriched 0.545, principal broad_enriched 0.530). The published order depended on two
non-target podocyte tiers (0.567 / 0.554) that no reconstructed map reproduces exactly. The substantive
point, that structural genes track the podocyte signal closely, is unchanged.

## 2. G1: sensitivities of the 157-gene podocyte HS set

Reference: g = 0.689 (0.042, 1.336). Default gate: same sign AND retention ≥ 0.50.

| Item | Variant | g (mHK 95% CI) | Missions + | Retention | Gate |
|---|---|---|---|---|---|
| 1 | signed-median scoring | 0.631 (−0.014, 1.276) | 4/5 | 0.92 | PASS |
| 2 | genes eligible in all 5 missions | 0.693 (0.044, 1.341) | 4/5 | 1.01 | PASS |
| 3 | omit OSD-102 | 0.846 (0.035, 1.657) | 4/4 | 1.23 | PASS (all 5) |
| 3 | omit OSD-163 | 0.614 (−0.184, 1.412) | 3/4 | 0.89 | |
| 3 | omit OSD-253 | 0.692 (−0.232, 1.616) | 3/4 | 1.00 | |
| 3 | omit OSD-462 | 0.739 (−0.189, 1.668) | 3/4 | 1.07 | |
| 3 | omit OSD-771 | 0.550 (−0.292, 1.392) | 3/4 | 0.80 | |
| 4 | leave-one-gene (min over 145 contributing genes) | 0.666 (0.020, 1.312) | — | 0.97 | PASS |
| 5 | OSD-462 mRNA preparation | 0.751 (0.100, 1.402) | 4/5 | 1.09 | PASS |
| 5 | OSD-462 UPX preparation | 0.703 (0.055, 1.351) | 4/5 | 1.02 | PASS |
| 6 | OSD-253 white-light rerun control | 0.802 (0.095, 1.509) | 4/5 | 1.16 | PASS |
| 7 | OSD-163 mapping-rate residualization | 0.669 (0.023, 1.314) | 4/5 | 0.97 | PASS |
| 7 | OSD-163 ANCOVA (descriptive) | 0.732 (0.013, 1.450) | 4/5 | 1.06 | — |

**Overall: ROBUST_DIRECTION.** It is the same under the runner-up gene map (`idmap_sensitivity.tsv`:
no flips).

**Family-level variants** (49 sets, 100,000 blocked permutations each):

| Variant | podocyte HS | Rank | max-T FWER |
|---|---|---|---|
| mean (primary) | 0.689 | 1 | 0.019 |
| signed median | 0.631 | 1 | 0.045 |
| common genes | 0.693 | 1 | 0.018 |

**Gene influence (item 4, plus descriptive diagnostics).**

- No gene carries more than 3.6% of the signed effect; the largest is Mgat5. The top 10 genes together
  carry 24%, and the effective number of positive contributors is 67.
- Canonical podocyte genes are among the top 37 contributors (Nphs2 and Col4a3 rank 6th and 7th; Myo1e,
  Wt1, Nphs1, Crb2 and Podxl follow), alongside Mgat5, Arhgef16, Aplp1, Dcdc2a and Mapt.
- **The leave-top-k curve shows the marginal margin:**

| Top-k positive contributors removed | g | Lower CI bound |
|---|---|---|
| 2 | 0.646 | 0.002 |
| 5 | 0.598 | −0.044 |
| 15 | 0.451 | −0.193 |
| 37 | 0.196 | −0.428 |

  Removing ranked contributors always lowers the estimate, so this is an influence diagnostic, not a
  test. It shows a signal spread across many genes whose interval clears zero only narrowly.

**Comparators (descriptive, not gated).**

- The structural sets behave like the podocyte set under every variant: full structural set 0.549, and
  0.680–0.341 across leave-one-mission.
- Omitting OSD-771 hurts the structural set most (0.33–0.34) and the podocyte set least (0.55).
- The structural set is slightly *stronger* under the OSD-462 mRNA/UPX preparations (0.72–0.74, intervals
  just above zero).

**K0 under the variants.** P1 stays positive in every variant: 0.083 for OSD-462 UPX up to 0.149 for the
OSD-253 rerun control. Its HC3 interval crosses zero in every variant. The direct difference D is
positive in every variant, and its interval also crosses zero in every variant.

### Item 8: leave-one-atlas-source-study-out marker rebuild

All 8 rebuilds pass QC and are evaluable. Thresholds are unchanged.

| Omitted study | Contributes podocytes | HS genes | Jaccard vs full | g (95% CI) | max-T | Rank | Retention |
|---|---|---|---|---|---|---|---|
| Conway20 | no | 151 | 0.93 | 0.718 (0.068, 1.368) | 0.012 | 1 | 1.04 |
| Dumas20 | no | 153 | 0.97 | 0.700 (0.053, 1.348) | 0.017 | 1 | 1.02 |
| Janosevic21 | no | 151 | 0.89 | 0.752 (0.101, 1.404) | 0.007 | 1 | 1.09 |
| Park18 | no | 140 | 0.89 | 0.728 (0.080, 1.377) | 0.010 | 1 | 1.06 |
| Hinze20 | yes | 155 | 0.79 | 0.669 (0.026, 1.313) | 0.025 | 2 | 0.97 |
| Kirita20 | yes | 162 | 0.83 | 0.697 (0.049, 1.346) | 0.017 | 1 | 1.01 |
| Miao21 | yes | 159 | 0.82 | 0.588 (−0.052, 1.229) | 0.077 | 1 | 0.85 |
| Wu19 | yes | 162 | 0.85 | 0.695 (0.028, 1.362) | 0.024 | 1 | 1.01 |

**DEFINITION_ROBUST: 8/8 pass.** The result is identical under the runner-up map (retention 0.85–1.09).

- Dropping a podocyte-contributing study makes the "≥75% of studies" rule effectively 3-of-3 and changes
  membership most (Jaccard 0.79–0.85), but the flight shift is retained.
- K0 P1 on the rebuilt sets ranges from 0.088 to 0.134. Every interval crosses zero.
- The atlas labels the Hinze study "Hinze20"; its README cites it as Hinze 2021.

Item 9 (a second independent atlas) is deferred; see `docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md`.

## 3. K0 permutation correction

The original max-|T| null drew independent label permutations per model, so the family maximum was taken
over mismatched shuffles. It also held the structural covariate fixed while permuting flight. The script
now reports three schemes; Freedman–Lane joint is the primary going forward.

| P1 family-wise p | Baseline map | Runner-up map |
|---|---|---|
| legacy independent streams | 0.0605 | 0.0827 |
| joint label permutation | 0.0388 | 0.0536 |
| **Freedman–Lane joint** | **0.0415** | **0.0563** |

The verdict follows the preregistered rule.

- P1's HC3 interval (−0.031, 0.255) crosses zero.
- The corrected p is below 0.05 only on the baseline map, so it is not robust to gene-ID reconstruction.
- **Permitted statement:** "P1 is unusual within {P1, S1, D} under a corrected null, but its interval
  crosses zero; still a bounded observation."
- Leave-one-mission estimates of the direct difference D range from 0.058 to 0.100, and every interval
  crosses zero.
- **Side observation from the K0 calibration simulations:** the parametric REML/mHK p is strongly
  conservative at k = 5 (0.8% rejection at nominal 5%). The permutation p-values are the calibrated ones.

## 4. G2: recovery persistence

**Design (preregistered).**

- Within OSD-771: 40 animals, ISS-terminal (ISS-T) and live-animal-return (LAR) arms, each arm with its
  own ground controls, blocked by arm × age.
- LAR mice flew 32 days and then spent 24 days on Earth. ISS-T mice flew 53–56 days and were euthanized
  on orbit.
- Model: OLS with HC3 errors; Freedman–Lane permutation within blocks, 100,000 permutations; per-arm
  exact enumeration (63,504 allocations).
- **Technical gate passed:** the LAR primary coverage metric is balanced. LAR flight animals do have
  higher RIN (g = 0.93, p = 0.045), a secondary metric that does not gate but should be kept in mind.

**Primary and key-secondary results (baseline map; identical labels under the runner-up):**

| Target | ISS-T arm g | LAR arm g | Interaction θ_L − θ_T, score units (95% CI) | FL p | Margin θ_L − 0.5θ_T (90% CI) | ρ̂ (Fieller 90%) | Label |
|---|---|---|---|---|---|---|---|
| **podocyte HS (primary)** | 1.17 (exact p 0.009) | 0.14 (p 0.73) | −0.281 (−0.622, 0.061) | 0.078 | −0.125 (−0.336, 0.087) | 0.10 (−0.66, 1.08) | **INCONCLUSIVE** |
| structural disjoint | 1.36 (p 0.004) | −0.43 (p 0.30) | −0.350 (−0.626, −0.073) | 0.009; max-T 0.024 | −0.213 (−0.407, −0.020) | −0.28 (−1.24, 0.42) | INCONCLUSIVE* |
| podocyte − structural | 0.16 (p 0.69) | 0.60 (p 0.16) | 0.069 (−0.265, 0.403) | 0.65 | — | — | NOT_EVALUABLE_NO_TERMINAL_EFFECT |
| P1-type adjusted podocyte (secondary) | θ_T 0.219 (−0.170, 0.607) | θ_L 0.083 | −0.135 | 0.50 (HC3) | — | — | NOT_EVALUABLE_NO_TERMINAL_EFFECT |

\* Structural disjoint would be DECAYS without the preregistered Holm adjustment across the two
key-secondary margins: unadjusted FL p = 0.026, Holm 0.053; HC3 90% upper bound −0.020. The flat-prior
summary gives P(ρ ≥ 0.5) = 0.03.

**Margin sweep (secondary, descriptive; requested by the owner):**

| κ | podocyte HS margin (90% CI) | Label | structural disjoint margin (90% CI) | Label |
|---|---|---|---|---|
| 0 (any retention) | 0.031 (−0.153, 0.215) | INCONCLUSIVE | −0.077 (−0.258, 0.104) | INCONCLUSIVE |
| 0.25 | −0.047 (−0.238, 0.144) | INCONCLUSIVE | −0.145 (−0.330, 0.039) | INCONCLUSIVE |
| **0.5 (primary)** | −0.125 (−0.336, 0.087) | INCONCLUSIVE | −0.213 (−0.407, −0.020) | INCONCLUSIVE* |
| 0.75 | −0.203 (−0.446, 0.041) | INCONCLUSIVE | −0.281 (−0.491, −0.072) | DECAYS |
| 1 (no decay) | −0.281 (−0.565, 0.003) | INCONCLUSIVE | −0.350 (−0.580, −0.119) | DECAYS |

Your hypothesis held: the 0.5 margin was INCONCLUSIVE, as predicted. The sweep shows that **no** margin
decides the podocyte set. The data are compatible both with substantial retention and with full decay;
at κ = 1 the 90% upper bound misses zero by 0.003. For the structural set, the sweep rules out retention
of ≥ 75% of the terminal shift.

**Context and exploratory analyses (descriptive; no label is a finding):**

- **49-set family.**
  - 37 of 48 exploratory sets have a negative interaction (smaller shift in LAR).
  - 13 sets pass BH q < 0.05 and 6 pass max-T < 0.05: four of the five endothelial tiers, TAL
    high-specificity and mesenchymal/stromal broad.
  - 31 sets have no established ISS-T effect.
- **Pattern concordance.** Per-arm effects across the 49 sets are *anti-correlated* between ISS-T and LAR
  (Pearson r = −0.65). The preregistered one-sided test for positive concordance gives p = 0.92. The
  negative correlation suggests a rebound pattern rather than simple fading. It is an unplanned
  direction and is descriptive only.
- **4 clinical axes.** The fibrosis axis has a negative interaction (FL p = 0.003). The other three have no
  ISS-T effect within OSD-771.
- **OSD-513** (live return: 37 days in flight, euthanized about 1 day after landing, ground controls 3
  days later; separate study; descriptive only):
  - podocyte HS g = 1.56 (exact p 0.003, max-T 0.057);
  - structural disjoint g = −0.18;
  - podocyte − structural g = 1.53 (max-T 0.009);
  - concordance with the pooled terminal pattern r = 0.37 (p 0.18).

  OSD-513 therefore shows a podocyte-specific shift that the terminal missions do not (K0). With k = 1
  and a different study design, this cannot be pooled or interpreted as recovery.
- **Exposure length (OSD-253, descriptive).**
  - Podocyte HS: day-25 g = 0.04, day-75 g = 1.73.
  - Structural: day-25 g = −0.19, day-75 g = 1.62.

  This dose-like pattern makes the LAR interpretation harder: LAR's shorter flight alone could explain a
  smaller shift.

**Power context.** Before the lock, the simulated probability of a PERSISTS label under full persistence
was 0.12 at the observed ISS-T effect size. INCONCLUSIVE was the expected outcome.

## 5. What may and may not be claimed

**May be claimed:**

- "Across five terminal missions, a podocyte-associated program shifts upward in bulk kidney (g = 0.69),
  and the direction and magnitude survive median scoring, gene influence, leave-one-mission, RNA
  preparation, control choice, mapping-rate adjustment and atlas leave-one-study-out marker rebuilds."
- "The shift is not separable from broad structural drift (K0): a shared structural drift with a
  podocyte-leaning tail."
- "Within RRRM-2, the terminal podocyte shift was not clearly present after 32 d flight + 24 d Earth
  recovery (LAR g = 0.14 vs ISS-T 1.17), but the design could not distinguish persistence from decay
  (interaction −0.28 [−0.62, 0.06]; MDE ≈ 1.8 SD)."
- "The structural scaffold shift was smaller in the LAR arm than the ISS-T arm (interaction p = 0.009;
  family-wise 0.024). Recovery, shorter exposure and on-orbit versus on-Earth euthanasia are aliased."

**May not be claimed:**

- that the podocyte program is separable, specific, or "the headline";
- that the signal "persisted" or "recovered";
- that the interval-crossing variants are "significant";
- anything about OSD-513 as recovery evidence;
- that the 0.5-margin label for the structural set is DECAYS.

## 6. Open items

- Item 9 (second atlas) is deferred; it should be un-deferred if "podocyte" stays in the title.
- A design that separates recovery from exposure length would need a terminal arm matched to LAR's
  32-day flight.
- The LAR flight-group RIN imbalance could be checked with a RIN-adjusted sensitivity. This was not
  preregistered.

## 7. Reproduce

```bash
python scripts/clinical_axes/build_reconstructed_id_map.py          # baseline gene map
RUN=data/results/clinical-axes/<UTC>_g1g2
for s in cross-mission compartment-context compartment-context-median compartment-context-common \
         podocyte-disjoint podocyte-scaffold-specificity podocyte-sensitivities atlas-loso-markers \
         atlas-loso-markers-idmap recovery-persistence recovery-persistence-idmap recovery-persistence-power; do
  python -m rrrm2 run clinical-axes --run-dir $RUN --stage $s --only
done
```

Key output hashes (sha256 prefix):

| File | sha256 prefix |
|---|---|
| `compartment_context_meta_results.tsv` | `a5a4e2d5` |
| `podocyte_scaffold_permutation.tsv` | `d0b93a83` |
| `podocyte_sensitivities/podocyte_sensitivity_summary.tsv` | `dfdc23f1` |
| `podocyte_sensitivity_gates.json` | `538e88be` |
| `atlas_loso/atlas_loso_summary.tsv` | `354a9ce0` |
| `recovery_persistence/persistence_contrasts.tsv` | `b80aa621` |
| `persistence_margin_sweep.tsv` | `5bdc59c1` |
