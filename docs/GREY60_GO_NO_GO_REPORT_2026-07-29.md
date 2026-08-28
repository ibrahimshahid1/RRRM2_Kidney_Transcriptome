# Grey60 adversarial go/no report

**Decision date:** 2026-07-29  
**Decision:** **NO-GO — retire Grey60 as a standalone recurrent-module paper**  
**Locked specification:** `docs/GREY60_ADVERSARIAL_REANALYSIS_PLAN_2026-07-29.md`  
**Production results:** `data/results/run_20260729_grey60_adversarial/`

## Executive conclusion

Grey60 is a real and unusually coherent **within-OSD-771 expression score**.
It is not, however, supported as a stable compact coexpression module or as a
recurrent terminal-flight kidney program.

Two failures are decisive:

1. **Flight-blind network reconstruction failed in all 27 specifications.**
   The primary BSL+VIV network placed only 12 frozen Grey60 genes in a broad
   511-gene module (Jaccard 0.022), and that module's projected ISS-terminal
   effect was null (blocked-permutation p=0.831).
2. **Compatible external recurrence failed.** Only two of four terminal-study
   point estimates were positive. The modified Hartung-Knapp terminal
   meta-estimate was Hedges g=0.163 (95% CI -0.658 to 0.984, p=0.573).

The locked rule required Gates A-E all to pass. Gates B, D, and E failed.
Gate D failed narrowly and is cautionary rather than independently fatal;
Gates B and E retire the branch by themselves.

## Gate decisions

| Gate | Question | Result | Decision |
|---|---|---:|---|
| A | Does the frozen score retain selection-adjusted internal evidence? | maxT p=0.00235 for the strongest expected ISS stratum; pooled bootstrap CI 0.956 to 1.760 | PASS |
| B | Is Grey60 recovered as a compact flight-blind module? | 12/48 overlap in the primary 511-gene module; Jaccard 0.022; projected p=0.831; 0/27 grid recoveries | **FAIL** |
| C | Is the result robust to mice, genes, hubs, and marker removal? | 48/48 gene contrasts positive; every leave-one-mouse/gene effect positive; maximum contribution 3.05% | PASS |
| D | Does it clear technical, composition, and baseline falsification? | 89.0% of effect retained after adjustment, but baseline absolute g=0.563 exceeded the locked 0.50 ceiling | **FAIL** |
| E | Does it recur in compatible terminal kidney cohorts? | 2/4 positive; meta g=0.163, 95% CI -0.658 to 0.984 | **FAIL** |
| F | Is it more than a matched generic state within OSD-771? | 82.3% retained after independent endothelial/fibroblast-score adjustment; matched-set percentile 1.00 | PASS, wording only |

## What survived

### A strong, distributed OSD-771 score

The primary 48-gene mean-z score increased in the ISS-terminal arm:

- pooled young/old flight-minus-control effect: 1.369 score units;
- young: 1.714;
- old: 1.023;
- pooled animal-bootstrap 95% interval: 0.956 to 1.760.

The selection-adjusted result is concentrated in the young stratum. Young
survived the 220-test maxT family (p=0.00235), whereas the concordant old
stratum did not independently survive that family (maxT p=0.334).

The result was not caused by one animal, hub, or gene:

- all 48 gene-wise ISS-terminal contrasts were positive;
- the exact sign-test p was 3.55e-15;
- no gene contributed more than 3.05% of the total absolute contrast;
- all leave-one-gene and leave-one-mouse effects remained positive;
- removing the top ten historical kME genes retained an effect of 1.275;
- removing all 17 Grey60 genes overlapping independent atlas-defined
  endothelial/fibroblast marker sets retained an effect of 1.417.

Technical and bulk-composition adjustment also did not erase the OSD-771
contrast. The raw-VST estimate fell from 1.470 to 1.309 after technical and
composition covariates, retaining 89.0% of its magnitude. The blocked
Freedman-Lane p was approximately 1e-4 and the maximum VIF was 3.86.

These results support:

> OSD-771 contains a large, distributed, ISS-terminal-associated
> vascular/stromal/ECM expression shift that is not explained by a few genes
> or by the measured bulk-composition covariates.

They do not establish a recurrent module, kidney specificity, cellular
localization, or causation by microgravity.

## Why the module claim failed

### Grey60 does not reconstruct as a compact flight-blind network object

The 48 genes were frozen, then WGCNA was rerun using only BSL+VIV samples so
that FLT/GC outcomes could not determine the reference network. The primary
5,000-gene, minimum-size-30, merge-height-0.25 network produced:

- only 35 of the 48 frozen genes in its variance-selected universe;
- a best-matching 511-gene blue module;
- 12 overlapping Grey60 genes;
- Jaccard overlap 0.0219;
- overlap q=0.000788;
- projected ISS-terminal effect 0.0354;
- blocked-permutation p=0.831.

The small overlap q value does not rescue compactness. A large 511-gene module
can be statistically enriched for 12 genes while remaining almost entirely
different from the proposed 48-gene object.

None of 27 prespecified WGCNA settings recovered Grey60. All best-overlap
modules had a positive projected point estimate, but zero met the recovery
criteria and the primary projection was null.

The older preservation result also referred to a different object. The
reported Zsummary of approximately 20.57 belongs to a broad 439-gene GC-only
green module containing 40 Grey60 genes. Compact combined-network Grey60 had
Zsummary approximately 7.80. WGCNA color labels are run-specific and cannot be
used as module identity.

The defensible interpretation is therefore:

> Grey60 is an unusually flight-responsive subset of a broader
> vascular/stromal program, not an independently recovered compact
> coexpression module.

### The baseline gate gave a caution, not a mechanistic explanation

Gate D failed because the BSL-versus-VIV negative-control Hedges g was 0.563,
narrowly above the locked absolute-g ceiling of 0.50. Its corresponding
factorial effects were imprecise and none survived maxT correction
(baseline maxT p=1.0).

This is not strong evidence that baseline context explains the OSD-771 flight
contrast. It is evidence that the score is somewhat responsive to housing,
handling, collection, or related context and therefore should not be described
as flight-specific. Gate B and Gate E remain the decisive failures.

## External recurrence

The frozen 48-gene score was projected without refitting membership. OSD-513
was designated in advance as a live-return moderator and excluded from the
terminal synthesis.

| Study | Context | Hedges g | 95% CI | Direction |
|---|---|---:|---:|---|
| OSD-102 | terminal | -0.342 | -1.483 to 0.800 | negative |
| OSD-163 | terminal | 0.625 | -0.539 to 1.789 | positive |
| OSD-253 | terminal, C57, original GC | 0.526 | -0.423 to 1.476 | positive |
| OSD-462 | terminal, three-preparation within-animal mean | -0.112 | -0.989 to 0.765 | negative |
| Terminal random-effects meta | terminal | 0.163 | -0.658 to 0.984 | unresolved |
| OSD-513 | live return | 2.456 | 1.200 to 3.712 | positive, excluded |

The terminal meta-analysis used REML with modified Hartung-Knapp uncertainty:
p=0.573, I-squared=0%, and maximum study weight 33.2%. I-squared=0% is not
evidence of biological equivalence here; four small, imprecise studies have
little power to detect heterogeneity. Unmodified Hartung-Knapp and
DerSimonian-Laird sensitivities also had intervals crossing zero.

OSD-253 was not a consistent duration replication. With a fixed score
reference, its original-control effects were -0.007 at approximately 25 days
and 1.283 at approximately 75 days. The pooled study estimate was driven by
the longer-duration cell.

OSD-462 is especially informative because the biological animals are matched
across RNA preparations:

| Preparation | Hedges g | 95% CI | permutation p |
|---|---:|---:|---:|
| UPX 3-prime tag | -0.091 | -0.968 to 0.786 | 0.838 |
| polyA mRNA | 0.975 | 0.042 to 1.908 | 0.0239 |
| total RNA | -1.162 | -2.117 to -0.207 | 0.0149 |

The same animals therefore yield near-null, positive, and negative Grey60
effects depending on RNA preparation. Averaging the three preparations in the
locked recurrence test prevents cherry-picking; it does not make the
underlying discordance biologically interpretable.

The strong OSD-513 result cannot rescue recurrence. OSD-513 is male and
post-return, whereas the primary recurrence question concerned compatible
terminal cohorts. It may reflect sex, recovery, timing, or other design
differences.

## Final scientific boundary

The repository supports this statement:

> In OSD-771, ISS-terminal flight is associated with a coherent,
> composition-robust vascular/stromal/ECM expression score, particularly in
> young mice.

It does not support:

- Grey60 is a preserved compact coexpression module;
- Grey60 recurs as a general terminal-flight kidney response;
- the response is kidney-specific rather than generic injury/remodeling;
- flight rewires the Grey60 network;
- the score localizes to one kidney cell type;
- microgravity, rather than the bundled mission and collection context, causes
  the score.

Grey60 may remain in a supplement, a design-heterogeneity analysis, or a
methods case study. It should not become the next standalone biology paper.

## Audit amendments

The first external execution exposed a technical-replicate suffix parsing bug
before a final decision was issued; it was corrected and the three OSD-462
preparations were verified to align to 20 unique animals.

An independent code audit also identified two conservative improvements:

1. OSD-253 primary and rerun-control sensitivities now use one fixed
   cross-scenario score reference rather than silently changing gene weights.
2. The primary meta-analysis uses modified Hartung-Knapp uncertainty, which
   cannot become narrower than the conventional random-effects interval when
   the small-study residual scale is below one.

Both corrections weakened the pooled result slightly and did not alter the
NO-GO decision. Complete 48-gene coverage and expected FLT/GC sample counts are
now asserted for every external matrix.

## Recommended next action

Proceed only with the locked **OSD-462 matched-library robustness** audit
described in `docs/NEXT_PAPER_DIRECTION_AFTER_GREY60_2026-07-29.md`.

That project asks an identifiable question using the same animals:

> How much does RNA library preparation alter the inferred kidney response to
> spaceflight?

It should run before a cross-mission control-choice paper because unresolved
OSD-462 preparation instability changes how OSD-462 can enter any larger RNA
synthesis.

