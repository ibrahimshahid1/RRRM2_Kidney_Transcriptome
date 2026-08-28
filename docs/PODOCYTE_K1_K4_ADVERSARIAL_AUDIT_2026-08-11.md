# Podocyte paper K1–K4 adversarial audit

**Date:** 2026-08-11  
**Status:** completed; post-hoc biological and statistical audit  
**Decision:** K1 changes the literature framing, while K2–K4 support retaining
the cross-mission podocyte-associated program as the repository's strongest
current biological paper.

## Executive verdict

| Check | Verdict | Consequence |
|---|---|---|
| K1: what Siew's PODXL/NID1 statement measures | **Mixed, and not urine-only.** The sentence juxtaposes Roscosmos urine proteomics with a separate cross-mission, cross-omic kidney consensus that includes mouse kidney RNA. | Do not claim that Siew saw only urine proteins. Present the new result as a program-level extension and a different mission-level estimand, with gene-level directional discordance noted. |
| K2: later competing claim | **No recurrent spaceflight podocyte-transcript program found** among forward citations screened through 2026-08-11. | The exact cross-mission program claim remains novel, although glomerular observations and individual markers are not novel. |
| K3: force `Podxl` and `Nid1` into the program | **Pass.** The set estimate is essentially unchanged; neither gene is individually negative overall. | The positive program is not concealing a two-gene contradiction. Do not present either gene as an individual discovery. |
| K4: same-tier strict matching | **Pass for the five observed missions.** The result survives common-support trimming and exact blocked-label testing. | The atlas abundance/specificity extreme does not explain the result. Generalization of the target-minus-matched contrast to future missions remains uncertain. |

## K1. Direct read of Siew Fig. 5 and Supplementary Data 3

### What the prose actually does

The PODXL/NID1 sentence appears inside a paragraph about proteins enriched in
human Roscosmos urine after flight. It then makes a second statement: NID1 and
PODXL were among the most frequently downregulated gene products **across
kidney datasets**, citing Fig. 5 and Supplementary Data 3. These are two
different measurements placed next to one another:

1. postflight human urine-protein abundance; and
2. a pan-mission kidney-tissue gene-product consensus.

The Fig. 5 caption confirms that the ranking score counts nominal Wald
`p < 0.05` direction calls and uses only kidney-specific transcriptome and
proteome datasets in the score. Epiproteomic, epigenomic, plasma, and exosome
datasets are shown only for visual comparison. Therefore the downregulation
statement is neither a Roscosmos-urine result nor a cross-species visualization
artifact.

### Direct categorical calls

Reading the Fig. 5 heatmap directly:

| Gene | Kidney-tissue calls contributing to the reported pattern | Supplementary Data 3 score |
|---|---|---:|
| `NID1` | Down in the BNL-3 kidney proteome and in mouse kidney RNA from RR-10, RR-23, RR-7 C3H at 75 days, and RR-7 C57 at 75 days | -25: five down, zero up |
| `PODXL` | Down in mouse kidney RNA from RR-7 C3H at 25 days, RR-10, MHU-3, RR-7 C3H at 75 days, and RR-7 C57 at 75 days; up in two mouse kidney proteomic columns | -21: five down, two up |

Supplementary Data 3 does not contain the per-dataset matrix; it contains the
gene-product consensus and score. Its Notes sheet defines the score as the
direction sum multiplied by the number of observed sets. `NID1` is row 3 of
the Downregulated sheet and `PODXL` is row 7.

### Correct interpretation against the present analysis

There is a real gene-level directional discrepancy with some of the same
missions, particularly RR-7/OSD-253. It is not, however, a contradiction
between like-for-like endpoints:

- Siew's display is a vote-count of nominally significant categorical calls
  across proteomic and transcriptomic datasets;
- the present primary endpoint is a continuous, 158-gene, animal-level
  expression score, followed by one effect per mission and blocked-label
  permutation; and
- the present gene-specific confidence intervals for `Podxl` and `Nid1` both
  include zero.

The manuscript should therefore say that prior work reported categorical
downregulation of individual PODXL/NID1 gene products, whereas the current
analysis identifies a broader recurrent podocyte-associated RNA program. It
should not advertise a correction of Siew, and it should not describe Siew's
observation as urine-only.

Primary source: [Siew et al. 2024](https://www.nature.com/articles/s41467-024-49212-1).

## K2. Forward-citation search

Forward citations to Siew 2024 and Finch/da Silveira 2025 were enumerated with
Crossref, OpenAlex, and Semantic Scholar on 2026-08-11. After version-pair
deduplication, 51 distinct Siew-citing publications were screened by title and
abstract; available full text was searched for `podocyt*`, `glomerul*`,
`PODXL`, `NID1`, `NPHS1`, `NPHS2`, and `SYNPO`. The three forward citations to
Finch/da Silveira were also screened.

No publication was found that reports recurrent elevation of a
high-specificity podocyte-associated bulk-kidney transcript program across
independent mouse spaceflight missions. Relevant neighboring findings remain:

- Siew's individual PODXL/NID1 and GCR glomerular observations;
- Finch/da Silveira's flight-downregulated `Wnt11` and discussion of
  glomerular development, without a podocyte-cell program claim;
- GCR-associated endothelial injury/thrombotic microangiopathy; and
- hindlimb-unloading glomerular morphology, which is simulated microgravity
  rather than flight and does not report a podocyte transcript program.

This is a targeted forward-citation audit, not proof that no uncited or
unindexed work exists. The precise novelty claim remains defensible as of the
search date.

Key sources: [Siew et al. 2024](https://www.nature.com/articles/s41467-024-49212-1),
[Finch/da Silveira et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11937539/),
and [Hoorn and Gameiro 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11630797/).

## K3. PODXL/NID1 forced-gene sensitivity

`Podxl` was already a member of the 157-gene high-specificity podocyte set;
`Nid1` was added. Both genes passed the frozen `0.1 CPM` eligibility rule in
every flight and control sample set across all five missions.

| Set | Pooled Hedges g | mKH 95% CI | mKH p | Blocked empirical p | Two-set maxT p | I2 |
|---|---:|---:|---:|---:|---:|---:|
| Original 157 genes | 0.6890 | 0.0416 to 1.3364 | 0.0418 | 0.00088 | 0.00091 | 0% |
| Forced 158 genes | 0.6859 | 0.0383 to 1.3335 | 0.0424 | 0.00097 | 0.00100 | 0% |

The two-set maxT value is a sensitivity-family calibration; it does not replace
the original 49-set compartment-family FWER of 0.0189.

### Individual genes

| Gene | OSD-102 | OSD-163 | OSD-253 | OSD-462 | OSD-771 | Pooled g (mKH 95% CI) | Direction count |
|---|---:|---:|---:|---:|---:|---:|---:|
| `Podxl` | -0.466 | 0.872 | 0.739 | 0.098 | 1.302 | 0.511 (-0.335 to 1.357) | 4 positive / 1 negative |
| `Nid1` | -0.378 | 0.086 | 0.857 | -0.529 | 1.975 | 0.389 (-0.887 to 1.666) | 3 positive / 2 negative |

Neither gene is an individually supported discovery; `Nid1` is notably
heterogeneous (`I2 = 72.7%`). The hypothesized failure mode—both genes negative
while the program is positive—did not occur.

Artifacts:

- `scripts/clinical_axes/run_podocyte_podxl_nid1_sensitivity.py`
- `data/results/run_20260811_clinical_renal_axes_cross_mission/podxl_nid1_forced_sensitivity/`

## K4. Strict same-atlas-tier matched audit

The stricter audit used 142 podocyte markers observable in every mission and
543 non-podocyte genes selected by the identical frozen high-specificity atlas
rule. Matching used seven flight-label-blind variables spanning bulk
abundance/variance, atlas abundance, specificity, source-study detection, and
kidney-compartment breadth.

### Balanced random panels

For both 10-nearest and 25-nearest candidate pools, 10,000 unique 142-gene
panels were retained only when every aggregate covariate imbalance was at most
0.25 robust SD. The target estimate was `g = 0.6926` (mKH 95% CI 0.0443 to
1.3409). No random panel was as extreme by estimate or absolute t statistic;
the finite-simulation p was `(0 + 1) / (10,000 + 1) = 0.000100`.

### Optimal one-to-one contrast and common support

| Target definition | Genes | Target-minus-matched g | mKH 95% CI | mKH p | Blocked empirical p | Two-scheme maxT |
|---|---:|---:|---:|---:|---:|---:|
| Full common-observable | 142 | 0.5420 | -0.1654 to 1.2495 | 0.1005 | 0.01480 | 0.01665 |
| q95 common-support trim | 133 | 0.5003 | -0.2225 to 1.2232 | 0.1270 | 0.02641 | 0.02956 |

The q95 trim removed `Asic2`, `Ceacam2`, `Magi2`, `Nphs2`, `Opcml`, `Srgap1`,
`Tdrd5`, `Thsd7a`, and `Wt1`. Mean imbalance after trimming was at most 0.142
robust SD across all seven variables.

Thus K4 passes for the observed missions. The mKH interval on the
target-minus-matched contrast still crosses zero, so the contrast is not shown
to generalize across a wider population of future mission designs. The control
pool is also endothelium-heavy: endothelial genes supplied 81 of 133 trimmed
matches. This arose from the frozen tier sizes and covariate proximity rather
than flight outcomes, but it remains a limitation.

Detailed report: `docs/PODOCYTE_STRICT_MATCHING_AUDIT_2026-08-11.md`.

## Manuscript consequence

The strongest defensible paper remains:

> **A recurrent podocyte-associated kidney transcript program across five
> mouse spaceflight missions**

The lead claim should be limited to higher bulk-kidney abundance of an
atlas-defined podocyte-associated RNA program across the five analyzed terminal
missions. The paper must not infer podocyte injury, protection, number,
filtration-barrier failure, albuminuria, or a urinary biomarker.

The Siew relationship should be written as follows:

> Prior pan-omic work placed PODXL and NID1 among recurrently downregulated
> kidney gene products and also detected postflight urinary changes in their
> proteins. Here, a continuous animal-level synthesis across five mouse
> missions identified a broader podocyte-associated tissue-RNA program whose
> average abundance was higher after flight; neither PODXL nor NID1 alone was
> a significant cross-mission effect.

That sentence acknowledges the direct overlap without manufacturing either an
agreement or a clean contradiction.

## Verification

```bash
venv/bin/python -m pytest -q \
  tests/test_clinical_axes_statistics.py \
  tests/test_clinical_axes_integration.py \
  tests/test_clinical_axes_matching.py
```

Result: **15 passed**.
