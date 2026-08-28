# State of play — 2026-07-29

One page. What's established, what's closed, what's still open. No
recommendations, no next steps. This exists so the record isn't only in
someone's head.

---

## Established (defensible with numbers)

**OSD-771 stromal/ECM score.** Pooled ISS-terminal effect 1.369 (young 1.714,
old 1.023), animal-bootstrap 95% CI 0.956–1.760. All 48 gene-wise contrasts
positive, exact sign-test p = 3.55e-15, no gene above 3.05% of the contrast, all
leave-one-gene and leave-one-mouse effects positive. 89.0% retained after
technical and composition adjustment; 82.3% retained after independent
endothelial/fibroblast adjustment. Scope: OSD-771, ISS-terminal arm, entangled
with collection context (baseline-vs-vivarium g = 0.563).

**Cross-cohort asymmetry.** Matrix/ECM signed Stouffer p = 7.03e-4 (median
I² = 63.8%); DCT/NCC-WNK p = 0.187 (median I² = 73.3%, unstable leave-one-out
sign). The matrix program recurs across cohorts; the transport program does not.

**OSD-462 design.** Condition perfectly aliased with reporter-tag block in both
plexes, no cross-plex swap (baseline 126–128C, flight 129N–131N, ground
131C–133C). Verified directly from ISA metadata.

**OSD-462 phosphoform identity.** Positions 65/68 are tyrosines; the T53- and
S383-indexed features come from co-modified peptides; zero isolated,
qualification-passing canonical NCC/SPAK phosphoforms.

**OSD-462 intensity gradient.** Gene-level suppression Z correlates with median
phosphopeptide signal at Spearman −0.633 (official-scaled) versus −0.143
(summed S/N). Lowest signal decile mean Z = +2.36, highest = −1.09.

**OSD-462 preparation dependence.** Same 20 animals: g = −0.09 (UPX), +0.98
(polyA), −1.16 (total RNA).

**OSD-462 reporter-position effect is small.** Within each block the five
channels hold five exchangeable animals, so a within-block trend is free of
biological confounding. Across all six block × plex estimates (~14–15k sites
each) no slope is significant: p = 0.080, 0.085, 0.615, 0.754, 0.654, 0.835.
Pooled slope +0.0072 log2 per channel step. Extrapolated across the five steps
between the flight and ground block centres this predicts −0.036 log2, against
an observed flight-minus-ground difference of −0.179 (plex 1) and −0.157
(plex 2) — roughly 20–23%, and not distinguishable from zero. The observed
contrast is 4–5× larger than a linear positional trend can produce.

Two caveats. A *step-shaped* tag effect at block boundaries would be invisible
to a within-block slope, so the aliasing is bounded, not resolved. And within a
block, position is confounded with animal identity (channels were filled in
animal-ID order), so the estimate is "position plus any animal-ID-ordered
nuisance." It remains the best available bound.
Run: `data/results/run_20260729_osd462_reporter_position/`.

**The flight-block dip is phospho-specific.** The protein workbook measures the
same animals in the same channel layout without Fe-NTA enrichment, so a
handling, loading or reporter-chemistry effect upstream of the phospho/protein
split should appear in both layers.

| Layer | Plex | n | Flight − ground (log2) | Fraction of features negative |
|---|---:|---:|---:|---:|
| Phospho | 1 | 14,959 | −0.179 | 85.6% |
| Phospho | 2 | 14,133 | −0.157 | 81.0% |
| Protein | 1 | 7,427 | −0.021 | 57.3% |
| Protein | 2 | 7,321 | −0.002 | 53.2% |

Paired within parent protein across 3,869 proteins quantified in both layers:
mean phospho-minus-protein difference −0.145 log2, bootstrap 95% CI
[−0.151, −0.140]; 3,356 of 3,869 proteins more negative at the phospho layer
(86.7%), exact sign test p < 1e-300.

The block profile is also non-monotone in both plexes — baseline high, flight
low, ground high — which is not the shape a positional drift produces.

What this excludes: shared upstream handling or loading, and reporter-position
chemistry (isotopic impurity does not care whether a peptide is
phosphorylated). What it does **not** exclude: anything specific to the phospho
arm downstream of the split — Fe-NTA enrichment efficiency, the separate
phospho labelling reaction (tc882-883 versus tc884-885), or phosphopeptide
loading. Condition remains aliased with tag block in both layers, so this is
still a flight-associated contrast, not a causal one.

Also note the effect is highly uniform (>80% of features in one direction),
which is equally consistent with a global technical shift in the phospho arm
and with a genuinely global biological one. The published Cosmic Kidney
observation of several thousand altered phosphosites and a general phosphatase
increase is the biological version of the same pattern.
Run: `data/results/run_20260729_osd462_layer_block_shift/`.

**ccRCC coordination boundary.** Glycolysis coordination transition robust at
patient level (7/7, exact p = 0.008); PT-identity coordination not (1/7,
p = 0.19). PT-identity loss prognostic as a bulk expression module (HR 0.64),
surviving adjustment for FBP1 and HIF/glycolysis. Per-patient network-
coordination score adds no prognostic value over mean expression.

## Closed

- **DCT2/CNT extension.** Non-evaluable at the frozen coverage gate (5 of 27
  observable) and negative in direction in every profile (−1.04 to −1.17).
- **ASDN extension.** Statistic 0.718, exact p = 0.0291 unadjusted, but
  p = 0.263 against an intensity-matched null. Fails selectivity.
- **Canonical NCC/SPAK dephosphorylation as an OSD-462 MS result.** No isolated
  qualified phosphoform. Prior antibody work on the same RR-10 animals is
  context, not replication.
- **Grey60 as a recurrent module.** 0/27 flight-blind reconstructions, Jaccard
  0.022, projected p = 0.831; terminal meta g = 0.163 (−0.658 to 0.984). The
  Zsummary ≈ 20.6 belongs to a 439-gene GC-only module, not compact Grey60
  (Z ≈ 7.8).
- **Spaceflight kidney ECM/TGF-β remodeling as a novel claim.** Published
  March 2025, da Silveira et al., npj Microgravity, doi:10.1038/s41526-025-00465-0.
- **Network features derived from X.** Capped by the DPI argument; empirically
  failed to beat expression under leave-one-cohort-out evaluation.
- **v12 and v13 as biology manuscripts.**

## Prior art that closed off proposed methods work

- Equivalence + differential testing for expanded analyte classification:
  QuEStVar (2024).
- Disattenuation of cross-layer correlation: Franks 2017, Upadhya & Ryan 2022.
- Power versus gene-set purity: pathway power analysis (2012).
- Curated versus data-derived signatures: BMC Bioinformatics (2020).
- Fit-and-test-on-same-data invalidity: the selective-inference and data-thinning
  literature.
- TMT reporter interference and design remedies: Brenes 2019, optTMT 2026.

## Open

- The paired small-n / large-n demonstration across RRRM-2 and TCGA KIRC. Both
  halves exist; nothing joins them.
- The silent-shifter construct — rewiring plus positive expression equivalence —
  which is the one network statistic not covered by Proposition 4.
- The within-block reporter-position diagnostic. Days of work, not done.
- OSDR design/detectability audit. Not started.
- `email_to_casaletto_draft.md`, dated 2026-06-19, unsent and now out of date.

## Provenance gap

`latex_paper/*` is gitignored except manuscript_v11, and the July run
directories are untracked. Freezing claims cannot currently be verified from
repository history.
