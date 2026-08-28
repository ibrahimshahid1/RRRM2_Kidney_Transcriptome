# Strict atlas-selection matching audit of the podocyte RNA result

**Date:** 2026-08-11  
**Status:** post-hoc adversarial sensitivity  
**Question:** does the five-mission high-specificity podocyte result disappear
when the comparison genes are required to occupy the same atlas
abundance/specificity stratum?

## Verdict

**K4 passes. The support does not collapse under the stricter reference.**

The stricter audit is stronger than the earlier six-gene all-gene matched
panel null in two ways:

1. every comparator gene was selected by the same frozen
   `high_specificity` atlas rule as the podocyte markers; and
2. matching directly included the atlas variables used to create that tier,
   as well as flight-label-blind bulk-expression abundance and variance.

The result remains a bulk-kidney podocyte-*associated* expression result. This
audit does not distinguish cell representation, tissue sampling, or
within-podocyte regulation.

## Frozen target and candidate pool

- Target: the 142 `podocyte__high_specificity` genes observable in all five
  primary terminal missions.
- Candidates: 543 genes from the other frozen high-specificity kidney
  compartments that were also observable in all five missions.
- No flight label or flight-effect estimate entered target filtering,
  candidate filtering, matching, or common-support trimming.
- Matching variables:
  - median VST abundance across missions;
  - median log VST variance across missions;
  - log podocyte/target-compartment atlas CPM;
  - log maximum non-target-compartment CPM;
  - log2 target-to-maximum-other specificity;
  - target-compartment source-study detection fraction; and
  - number of non-target compartments in which the gene was detected.

## Analysis 1: balanced same-tier random panels

For every podocyte target gene, candidate pools were restricted to its 10 or
25 closest other-compartment high-specificity markers. Every random panel had
142 unique genes. A panel was retained only when its mean differed from the
target by no more than 0.25 robust standard deviations on **every** matching
variable.

| Candidate pool | Accepted panels | Acceptance | Target g | Target mKH 95% CI | One-sided p | Two-sided p | t-based p |
|---|---:|---:|---:|---:|---:|---:|---:|
| 10 nearest | 10,000 | 52.8% | 0.6926 | 0.0443 to 1.3409 | 0.000100 | 0.000100 | 0.000100 |
| 25 nearest | 10,000 | 9.79% | 0.6926 | 0.0443 to 1.3409 | 0.000100 | 0.000100 | 0.000100 |

No sampled same-tier panel reached the observed target estimate or absolute
t statistic. The reported Monte Carlo value is therefore
\((0+1)/(10{,}000+1)\), not zero.

For the 10-neighbour audit, mean target-to-panel imbalances ranged from 0.011
to 0.213 robust SD in absolute value, and the enforced maximum was 0.250. The
largest residual mean differences were lower specificity (-0.213), lower
atlas target-compartment abundance (-0.181), higher bulk variance (0.168), and
higher bulk abundance (0.161) among matched panels.

## Analysis 2: one-to-one optimal same-tier contrast

A separate analysis selected a single minimum-total-distance comparator for
each target, without replacement. It then tested the within-animal difference
between the podocyte panel score and its matched comparator by 100,000
flight-label permutations within the declared mission strata.

| Target definition | Genes | Target-minus-matched g | mKH 95% CI | mKH p | Blocked empirical p | maxT across two schemes |
|---|---:|---:|---:|---:|---:|---:|
| Full common-observable tier | 142 | 0.5420 | -0.1654 to 1.2495 | 0.1005 | 0.01480 | 0.01665 |
| q95 common-support trim | 133 | 0.5003 | -0.2225 to 1.2232 | 0.1270 | 0.02641 | 0.02956 |

For the stricter common-support analysis, a target was retained only if its
nearest non-podocyte candidate distance was no greater than the 95th
percentile of cross-compartment nearest-neighbour distances among candidate
markers. This removed nine unusually difficult-to-match targets:
`Asic2`, `Ceacam2`, `Magi2`, `Nphs2`, `Opcml`, `Srgap1`, `Tdrd5`, `Thsd7a`,
and `Wt1`.

After trimming, absolute panel-mean differences were at most 0.142 robust SD
across all seven variables. The largest assigned pair distance was 1.625. The
q95 rule is a target-overlap rule rather than a hard pairwise caliper; the
one-to-one constraint can assign an individual pair slightly beyond the
1.567 target-overlap threshold.

The q95 target-minus-comparator contrast was positive in OSD-102, OSD-163,
OSD-253, and OSD-462 and negative in OSD-771. Thus the observed-mission label
test supports selectivity, whereas the mKH interval does not establish that
the target-minus-comparator contrast generalizes to a wider population of
future missions.

## Remaining limitations

- This was designed after the podocyte result was known and is a sensitivity
  analysis, not prospective confirmation.
- The available same-tier pool is itself unbalanced: endothelial markers are
  the largest category and supplied 81 of 133 q95 optimal matches. That is a
  consequence of covariate proximity and frozen tier sizes, but it means the
  comparator is not a uniformly balanced mixture of kidney compartments.
- Nine podocyte targets, including `Nphs2`, `Magi2`, and `Wt1`, lacked q95
  common support. The positive result after removing them is reassuring, but
  no matching procedure can invent overlap where the reference contains
  none.
- Random-panel p values are annotation-null sensitivities. The blocked
  animal-label test of the target-minus-matched score is the stronger
  inferential result.

## Reproduction

```bash
PYTHONPATH=. venv/bin/python scripts/clinical_axes/run_strict_podocyte_matching_audit.py
PYTHONPATH=. venv/bin/pytest -q tests/test_clinical_axes_matching.py tests/test_clinical_axes_statistics.py tests/test_clinical_axes_integration.py
```

Machine-readable outputs are under
`data/results/run_20260811_clinical_renal_axes_cross_mission/strict_podocyte_matching/`.
