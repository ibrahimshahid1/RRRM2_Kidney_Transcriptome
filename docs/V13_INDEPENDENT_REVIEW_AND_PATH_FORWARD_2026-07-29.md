# Independent review of the v13 reanalysis, and what to do with this project

**Date:** 2026-07-29
**Scope:** review of `V13_FINAL_REANALYSIS_AND_MANUSCRIPT_DECISION_2026-07-29.md`,
`V13_MANUSCRIPT_DECISION_MEMO_STAGE0_2026-07-29.md`,
`MANUSCRIPT_V1_V12_CROSS_VERSION_AUDIT_2026-07-29.md`, the Stage 0 provenance
audit, and `run_20260729_v13_continuous_phospho_exact_final/`.
**New analysis produced by this review:**
`scripts/v13/intensity_confound_audit.py` →
`data/results/run_20260729_v13_intensity_confound_audit/`.

---

## 1. Bottom line

The v13 work is correct, and its conclusion holds. I re-derived the load-bearing
numbers from the emitted artefacts and re-read the raw ISA metadata by hand. The
late-distal branch is closed.

One correction, in the direction of *less* support rather than more: the ASDN
result that v13 reports as "positive but non-selective" does not survive a
control that was not run. Gene-level suppression scores are strongly
intensity-dependent (Spearman −0.63 against median phosphopeptide signal). When
gene sets are tested against an intensity-matched competitive null, ASDN's
statistic loses most of its magnitude and its significance
(**p = 0.263** in the primary profile, versus 0.029 unmatched). ASDN should be
reported as *not distinguishable from intensity-matched background*, not as a
positive-but-unselective shift.

Two further points the v13 documents leave implicit:

- The DCT2/CNT set is not merely under-covered; its five observable genes point
  **against** the hypothesis in every quantification profile (statistic −1.04 to
  −1.17). "Non-evaluable" is the right formal verdict against the frozen n ≥ 8
  gate, but the observed direction should be stated so no one revisits it hoping
  for a hidden positive.
- The only annotation that survives both intensity matching and the choice of
  quantification column is **podocyte**. This is not a new finding to chase —
  see §4 — but it is the honest description of what is left in the phospho
  layer.

The project's problem was never statistical technique. It is that the OSD-462
phosphoproteomics cannot yield an identifiable biological claim, for two reasons
that no reanalysis can repair, and the descriptive claim it *can* yield has
already been published on the same animals. §5 covers this. §6 gives three
products that are real.

---

## 2. What I verified

| Claim in the v13 documents | Status | How |
|---|---|---|
| Condition perfectly aliased with reporter-tag block; no cross-plex swap | **Verified independently** | Read `a_OSD-462_phosphoprotein-profiling_mass-spectrometry_Orbitrap Eclipse.txt` directly. All 30 rows: baseline → 126, 127N, 127C, 128N, 128C; flight → 129N, 129C, 130N, 130C, 131N; ground → 131C, 132N, 132C, 133N, 133C. Identical map in plex 1 and plex 2. |
| 30 WT female samples, two 15-plexes, 5/5/5 per plex, primary contrast 10 flight vs 10 ground, no p21-null arm in the MS workbooks | **Verified** | Same file; sample names `RR10_KDN_WT_{BSL,FLT,GC}_*`. |
| NCC positions 65/68 are tyrosines, not canonical regulatory residues; T53 and S383 rows come from co-modified peptides; zero isolated qualified canonical phosphoforms | **Consistent and internally documented** | Stage 0 emits peptide sequence, accession, AScore, redundancy and reporter completeness per row. Mouse NCC T53/T58/S71 and SPAK T243/S383 are the correct canonical references. |
| DCT2/CNT 5 of 27 observable; ASDN 16 observable, statistic 0.718, exact p = 0.0291; podocyte 1.112, BH q = 0.0021; fibroblast 0.707, q = 0.0169 | **Verified** | `reporting/primary_compartment_enrichment.tsv`, `set_level_permutation_inference.tsv`. |
| All 63,504 balanced within-plex assignments enumerated | **Verified** | `n_null_valid = 63504` throughout. |
| No hidden stronger DCT2/CNT result in v1–v12 | **Agreed** | The cross-version audit is accurate and unusually candid; v11 introduced the question and v12 only promoted the same thresholded result. |

The provenance discipline here is good — frozen gene sets are copied into the
run directory, gates are machine-readable, the coverage threshold was not
relaxed when DCT2/CNT came in under it. That last point matters: it is the
single strongest piece of evidence that this run was not steered toward a
desired answer.

---

## 3. The intensity confound (new)

### What the v13 pipeline controls, and what it does not

Each parent gene is standardised against its own balanced-label null, which
correctly removes gene-specific *variance* differences arising from site count
and missingness. It does not remove a systematic *shift* in the observed
effects that tracks phosphopeptide signal, because the permutation null is
centred at zero for every gene by construction.

That shift is large:

| Signal decile (primary profile) | Mean median log2 signal | Mean gene Z |
|---:|---:|---:|
| 1 (lowest) | 2.504 | **+2.36** |
| 5 | 2.705 | +0.13 |
| 10 (highest) | 2.772 | **−1.09** |

Spearman(gene Z, median signal) = **−0.633** over 3,524 eligible genes. The
gradient spans ~3.4 Z units — roughly three times the largest compartment
statistic in the paper. Any gene set weighted toward low-signal phosphoproteins
will show apparent "suppression" for purely technical reasons.

### Re-testing every frozen set against an intensity-matched null

Competitive statistic recomputed exactly as in v13, but with the gene-label
null drawn within 20 strata of median phosphopeptide signal (20,000
permutations; primary canonical-axis exclusion applied):

| Gene set | n | Statistic | Intensity-matched p — primary | S/N-sum centred | Scaled uncentred |
|---|---:|---:|---:|---:|---:|
| Podocyte | 44 | +1.11 | **0.0054** | **0.00005** | **0.0148** |
| Fibroblast | 34 | +0.71 | **0.0041** | **0.0032** | 0.277 |
| Immune | 20 | +0.71 | **0.0163** | **0.0193** | 0.591 |
| **ASDN** | 16 | +0.72 | **0.263** | 0.035 | 0.576 |
| Thick ascending limb | 16 | +0.59 | 0.621 | 0.087 | 0.402 |
| DCT1 | 9 | +0.19 | 0.126 | 0.464 | 0.185 |
| Endothelial | 35 | +0.35 | 0.441 | 0.101 | 0.943 |
| Principal cell | 34 | +0.07 | 0.926 | 0.321 | 0.967 |
| Proximal tubule | 60 | −0.32 | 0.995 | 0.718 | 0.995 |
| Intercalated cell | 32 | −0.43 | 0.997 | 0.856 | 0.994 |
| **DCT2/CNT transition** | 5 | **−1.12** | 0.981 | 0.900 | 0.995 |

Three consequences.

1. **ASDN fails a control it should have passed.** 71 % of its statistic in the
   primary profile is attributable to intensity composition, and it is
   significant in only one of three quantification profiles. The v13
   observability-matched test (p = 0.247) was pointing at exactly this but was
   filed as a secondary; it should be primary. The permitted-claim sentence in
   the final memo ("parent proteins in the predefined ASDN set had higher
   standardised phosphosite-suppression scores than the remaining observable
   proteins") is not supportable once intensity is matched, and should be
   struck or explicitly conditioned.
2. **The normalisation profiles do not agree as much as the memo implies.** The
   intensity–Z correlation is −0.633 under official-scaled values but −0.143
   under summed signal-to-noise; the uncentred profile eliminates fibroblast,
   immune and ASDN entirely. Calling this a "provenance/normalization
   sensitivity" understates it: the site-level effect estimates are
   substantially determined by the quantification column chosen.
3. **Podocyte is the only survivor** across all three profiles under intensity
   matching. Note what the podocyte marker set contains — Synpo, Magi2, Iqgap2,
   Myo1e, Parva, Nck2, Pak1, Ablim3, Arhgef18, Tmod3, Mtss1, Plekhh2 — an
   almost purely cytoskeletal/adhesion, phosphorylation-dense set, alongside a
   fibroblast set that is similarly contractile (Acta2, Myl9, Tagln, Cnn1,
   Myh11, Tpm2). The proximal-tubule set, at the opposite end, is metabolic
   enzymes and solute carriers. The compartment axis is plausibly a
   protein-class axis.

### One methodological note on set construction

The frozen DCT2/CNT-transition set is Anxa11, Cdh16, Cdo1, Gata3, Hif1a, Hoxd4,
Hoxd9, Pax2, Peli2, Phactr1, Tbx2, Tbck, Trpv4, Slc2a9 and several
`Gm`/`Rik` transcripts — developmental transcription factors and low-expression
genes. The biologically canonical late-distal markers (Trpv5, Calb1, S100g,
Scnn1b/g) were partitioned into ASDN. So the DCT2/CNT non-evaluability is partly
a consequence of non-overlapping set design plus n = 3 pseudobulk, not purely of
phosphoproteomic coverage. This does not change the verdict — a better set would
still be sparsely observed — but it should be recorded, because a reviewer or
collaborator will ask whether the set was a fair test of the hypothesis.

---

## 4. Do not chase podocyte

The temptation is obvious: podocyte survives everything, spaceflight glomerular
change is topical, there is a figure already. Resist it, for three reasons.

1. It is still inside the tag-aliased contrast. Nothing about a podocyte
   annotation makes the flight coefficient causally interpretable.
2. It is annotation, not localisation — the same objection that killed DCT2/CNT
   applies verbatim.
3. It is very likely pre-empted and probably not even a compartment effect. The
   Cosmic Kidney Disease study already reported, from this material, "several
   thousand significantly altered phosphosites, suggesting a general increase in
   phosphatase and decrease in kinase activities." A broad downward phospho
   shift concentrated on cytoskeletal/adhesion scaffolds is a refinement of a
   published observation, not a discovery.

Chasing it would be the thirteenth instance of the pattern in §7.

---

## 5. Why there is no v14

Two identifiability failures, neither fixable by analysis:

- **Design.** Flight occupies reporter channels 129N–131N and ground occupies
  131C–133C in both plexes, with no label swap. Isobaric-tag isotope impurity
  and interference are systematically structured by channel position, so the
  flight coefficient and a tag-block effect are the same parameter. Plex
  covariates, median centring, and label permutation all leave it intact —
  permutation preserves the dependence structure but cannot reconstruct a
  randomisation that never happened.
- **Attribution.** Whole-kidney phosphoproteomics identifies a phosphosite's
  parent protein. No annotation, atlas, or specificity score converts that into
  a cell of origin. Every subtype claim from this assay is a claim about gene
  annotations, permanently.

And a third, external:

- **Pre-emption.** NCC dephosphorylation in RR-10 was published with antibody
  validation, and the global phosphosite shift was published too — both from the
  same animals ([Nat Commun 2024](https://www.nature.com/articles/s41467-024-49212-1);
  [Clin Kidney J 2024](https://academic.oup.com/ckj/article/17/12/sfae329/7849774)).
  The fallback position is occupied.

Your "F = ma" objection is correct, and it cuts the way you didn't want it to.
A laborious re-derivation of a known result is a weak paper — which is an
argument against publishing the NCC finding as novelty, not an argument for
inflating the DCT2 finding to fill the gap. The escape from that trap is a
different question, not a bolder claim on the same data.

---

## 6. Three products that are real

### A. OSD-462 data-quality / reanalysis note — do this first

Highest value per unit effort, lowest risk, ~85 % already written in this
repository. It requires no new data and no new positive result.

Content:

1. The condition-by-reporter-tag map, with the observation that no cross-plex
   swap exists and therefore what is and is not estimable from this dataset.
2. The TMTpro-versus-legacy-iTRAQ metadata inconsistency.
3. The residue-indexing problem: positions 65/68 are tyrosines; the T53 and
   S383 rows are co-modified peptides; zero isolated qualified canonical
   phosphoforms. Include the peptide-level appendix Stage 0 already emits.
4. The intensity-dependent gradient in the processed workbook values and its
   dependence on quantification column (§3) — a direct, quantitative caution
   for anyone doing set enrichment on these tables.
5. A short "what this dataset can and cannot answer" section.

Why it's worth doing: OSD-462 is public and will be reused. These are exactly
the errors a downstream user makes silently. It converts three months of audit
work into something citable, and the honest framing is a strength rather than a
liability. Venue: npj Microgravity, *Life Sciences in Space Research*, or a
preprint paired with a GeneLab curation ticket on the assay-label field. Also
worth emailing GeneLab about the metadata field regardless of publication —
that's a service contribution with no downside.

Effort: 2–4 weeks, mostly writing. The v13 manuscript is already close; it needs
its scope narrowed from "null-boundary report on a biological hypothesis" to
"reanalysis note on a public dataset," which is a shorter and more natural
document.

### B. The RNA / stromal paper — the one branch with a live biological question

This is the only branch that is free of the isobaric-tag confound (RNA-seq is
not isobarically labelled), has independent cohorts, and has a positive result
that has survived several rounds of adversarial review.

Two components, and the second is the interesting one:

- The preserved 48-gene Grey60 ECM/cell-migration module with an
  ISS-terminal-arm-associated eigengene shift (young q = 6.1 × 10⁻⁵, old
  q = 0.0108; ECM enrichment OR 29.8, q = 0.0138; GC-reference preservation
  Z_summary = 20.6).
- The cross-cohort asymmetry: matrix/ECM RNA recurrence is reproducible across
  five cohorts (p = 7.0 × 10⁻⁴) while the distal-transport RNA response is
  heterogeneous and non-significant (p = 0.19, I² ≈ 73 %).

The second point is a corrective to the field rather than a me-too result: the
spaceflight-kidney story that gets the attention (distal transport) is the less
reproducible one across cohorts, and the unglamorous one (matrix remodelling) is
the reproducible one. That is a defensible, non-obvious claim with a clear
message, and it does not require any phosphoproteomics.

It must pass the adversarial go/no-go before writing:

1. Reconstruct the module and contrast exactly from committed code.
2. Show it is not driven by collection arm, dissection/recovery context, age
   imbalance, or cell composition — ISS-T is entangled with collection context
   and this is the main threat.
3. Random-effects meta-analysis over cohort-level standardised effects, with a
   confidence interval per cohort, replacing the weighted Stouffer p-value.
4. Bootstrap both cohorts, not just the external one — the current interval is
   too narrow because the RRRM-2 reference vector is held fixed.
5. Leave-one-gene and leave-one-cohort sensitivities; check the module is not
   generic injury/ECM activation that would appear in any stressed kidney.
6. Say "eigengene shift," not "activation."

If it fails, retire it and stop. Do not rescue it.

### C. LayerScore methods paper — already spun off, now has its worked example

The v13 machinery is a genuinely reusable contribution: continuous parent-gene
scores, whole-pipeline label permutation with exact enumeration, per-gene null
standardisation, frozen coverage gates, machine-readable claim tiers, and now
intensity-matched competitive nulls. It has something most methods papers lack —
a fully documented case where the framework correctly refused to certify a
hypothesis its authors wanted. Fold §3 in as a worked demonstration of why
per-gene variance standardisation is insufficient without intensity matching.

---

## 7. What to retire, and the pattern worth naming

Retire from all main texts: PPAR sex dimorphism, LIONESS/node2vec, age rewiring,
LAR reversal, TLR4 and S1P branches, KSEA as independent confirmation, human
urine, CMap, IRI spot-level tests, dDAVP, low-potassium analogy, aging
projection, composition M0–M5 ladders, the cohort-level MR Spearman, and
canonical NCC/SPAK dephosphorylation as an OSD-462 MS result. The v13 documents
already say this; it is right.

The pattern underneath twelve revisions: the dataset stayed fixed and the
hypothesis moved. PPAR → networks → age → live return → multi-omic anchor →
DCT1 → DCT2/CNT → ASDN. Each pivot was locally reasonable and several were
scientifically admirable corrections. But the informative content of
RRRM-2 + OSD-462 was largely extracted by about v8, and everything after has
been re-slicing the same variance — which is why each successive headline is
weaker than the one before, and why the prose had to become more defensive to
carry it.

That is also the answer to your prose complaint. The Codex estimate of ~8.5/10
for conspicuous AI drafting is fair, but the litigating-every-paragraph style is
downstream of the claim, not the writing. When a claim needs five qualifications
to survive, every paragraph has to relitigate it. Products A and B both have
claims that survive a single plain sentence, and their prose will get short
almost automatically. Do not do a line-edit pass on v12 or v13 prose — rewrite
from a blank document once the scope is fixed.

---

## 8. Immediate next actions

1. Add §3 to the v13 record: the intensity gradient, the intensity-matched
   re-test, the struck ASDN permitted-claim sentence, and the observed negative
   direction of the five DCT2/CNT genes. Artefacts are already in
   `data/results/run_20260729_v13_intensity_confound_audit/`.
2. Re-scope v13 from "null-boundary biology report" to "OSD-462 reanalysis and
   data-quality note" (product A). Rewrite from a blank file.
3. Email GeneLab about the assay-label field and the residue indexing.
4. Run the Grey60 go/no-go as a locked, prespecified audit before any writing.
   One decision, one document, no narrative until it passes.
5. Commit the v13 run directories and generators. The `.gitignore` exclusion of
   the July outputs is the one remaining provenance hole, and it undercuts every
   "frozen before testing" statement in the record.
6. Do not open a v14 of the late-distal manuscript.

A defensible DCT2/CNT paper needs: counterbalanced or label-swapped isobaric
allocation, targeted assays that isolate the phosphoforms of interest,
segment-resolved sampling, and matched electrolyte/aldosterone/blood-pressure
physiology. That is a new experiment, not a new analysis. Worth proposing;
not worth simulating.
