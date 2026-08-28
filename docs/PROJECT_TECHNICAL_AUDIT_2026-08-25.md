# RRRM2_Kidney_Transcriptome — Deep Technical Project Audit

**Audit date:** 2026-08-25
**Auditor basis:** direct inspection of the repository at `main` (HEAD `1d165c4`), including source, configs, tests (executed), run outputs, manuscripts, planning documents, and commit history. Not based on README claims.
**Repository state at audit:** 67 commits (2025-12-15 → 2026-07-17); 19 modified and 63 untracked paths; last file-level activity 2026-08-11; test suite executed on 2026-08-25 → **225 passed, 0 failed, 26.4 s**.

---

## 1. What this project is

RRRM2_Kidney_Transcriptome is a **computational research project in statistical genomics**, not a software product. Its object of study is the renal response of mice to spaceflight, and its output is a sequence of inferential claims about that response — most of which the project itself has tested and retired.

### Origin: this did not begin as a reanalysis

The current README calls the project a "reanalysis." That word describes the last quarter of its life and should not be read backward onto the whole. The founding artifact is an ASGSR-style poster/proposal dated **9 December 2025** (`node2vec_asgsr_osd771_kidney_transcriptomics.pdf`), six days before the first commit. The first README states the motivating gap in the project's own terms:

> The Gap: Transcriptomic mechanisms underlying these changes remain unclear, specifically regarding genes that rewire networks without changing expression levels.

That is an original hypothesis about an unexamined class of signal, not a check on a published result. It was pursued on a newly released dataset chosen for a specific structural reason — per the project's own retrospective, RRRM-2/OSD-771 is "the only kidney dataset that crosses age, mission arm, and flight in a single design." The category the work set out to find, **silent shifters**, was named by the project, not borrowed. The word "reanalysis" first appears in `docs/` on 2026-07-28, with the v13 rebuild.

The accurate three-part shape is therefore:

| Period | What the project was |
|---|---|
| Dec 2025 – May 2026 | **Original method-driven investigation.** New dataset, self-identified gap, novel application of graph embeddings to a co-expression rewiring question, author-named target category, conference proposal |
| May – Jul 2026 | **Self-audit and repeated pivot.** The original method failed its own prespecified validation; the dataset stayed fixed while the question moved |
| Jul 2026 – present | **Reanalysis *and* original synthesis.** The v13 phospho work is genuinely a reanalysis of OSD-462; the cross-mission layer is not a reanalysis of any published claim — by the project's own forward-citation audit, no published equivalent exists |

Calling the whole project a reanalysis understates it in a specific way: it implies the question came from someone else's paper. It did not.

It spans four activity types that are genuinely entangled:

1. **Data analysis / bioinformatics.** Bulk RNA-seq and isobaric-tag phosphoproteomics from NASA Open Science Data Repository (OSDR) missions.
2. **Statistical methodology.** Purpose-built exact and blocked permutation inference, meta-analysis, and null-calibration machinery that does not exist off the shelf in the form used here.
3. **Forensic data auditing.** Reconstruction of assay provenance from raw ISA metadata and MS workbooks, which overturned the project's own prior headline.
4. **Scientific writing.** Thirteen numbered manuscript versions plus supplements, audits, posters, and a retrospective.

### The problem it addresses

**As posed at the start:** ordinary differential expression measures how much each gene changes, not how genes change what they are co-regulated *with*. The project's founding claim was that spaceflight-associated renal biology includes regulatory reorganization invisible to DE, and that graph embeddings of co-expression networks could measure it — quantifying, per gene, how far its network context moved between flight and control, and whether flight amplifies, dampens, or redirects the ordinary aging trajectory.

**Field context, not the trigger:** the flagship result (Siew et al. 2024, *Nat Commun*, doi:10.1038/s41467-024-49212-1) reports distal-nephron remodeling and NCC dephosphorylation after spaceflight. That work supplied the biological anchor genes and the tubule-centered expectation the project tested against; it did not supply the question.

**As the problem now stands:** after the rewiring premise failed its own validation and the OSD-462 phospho anchor failed provenance audit, the operative question became **which reported renal spaceflight effects actually recur when the analysis unit is one animal and one mission rather than one gene or one dataset, and how large an effect can be excluded for those that do not.**

### Inputs

| Input | Detail |
|---|---|
| Primary cohort | OSD-771 / GLDS-674 — RRRM-2 mouse kidney bulk RNA-seq, **n = 80** female C57BL/6NTac, 57,278 genes. Full factorial: Age (Young 16 wk / Old 34 wk) × Arm (ISS-Terminal / Live-Animal-Return) × Environment (FLT / HGC / VIV / BSL), 5 mice per cell (`config/metadata_design.yaml`) |
| External missions | OSD-102, OSD-163, OSD-253, OSD-462, OSD-513 bulk RNA-seq from OSDR (`data/external/osdr/`) |
| Phosphoproteomics/proteomics | OSD-462 (RR-10) two 15-plex isobaric workbooks, 20 wild-type female mice |
| Single-cell references | Mouse Kidney Atlas; Tabula Muris Senis kidney; Chen atlas; GSE228367; GSE150338 |
| Spatial references | GSE269622 Visium, GSE269719 Xenium |
| Prior-knowledge resources | Johnson 2023 Ser/Thr kinome atlas, Yaron-Barir 2024 Tyr atlas, LINCS/CMap (GSE92742), decoupleR priors, KEGG/GO/Reactome sets |
| Human context | OSD-656 Inspiration4 urine; NASA Twins Study supplementary tables |

### Outputs

Frozen YAML specifications; timestamped run directories (81 under `data/results/`) containing TSV result tables, permutation nulls, PDF/PNG/SVG figures, and JSON manifests carrying SHA-256 hashes of inputs, code, outputs, and git state; decision documents in `docs/`; and LaTeX manuscripts.

### What distinguishes it from a simpler implementation of the same idea

Three things, concretely:

- **The inference is built, not called.** The v13 endpoint is a full enumeration of all 63,504 balanced within-plex label assignments, with each parent gene standardized against *its own* null distribution — computed by re-running site contrasts, gene aggregation, and set statistics for every assignment. There is no library that does this; `src/v13/continuous_phospho_inference.py` (3,048 lines) does it in two chunked passes.
- **Claim gates are code.** `config/*.yaml` files encode minimum observable-gene counts, selectivity criteria, and go/no-go rules *before* the run; the engine emits `claim_tier.tsv` and `claim_gates.tsv` as computed artifacts. The current tier is the machine-produced string `neither`.
- **Negative controls are frozen into the design.** `broad_structural_scaffold_control__all` exists in the compartment family for the express purpose of killing the compartment framing if it wins. That control now ranks 4th of 49 — and the project's own documentation says so rather than hiding it.

---

## 2. The technical workflow as it actually exists

There are **four parallel workflows**, built in sequence. Only the fourth is the current development frontier; the first is retired; the second and third are complete and locked.

### Workflow A — RNA network rewiring (Dec 2025 – May 2026). Retired.

```
raw counts (n=80) → MuSiC deconvolution (R) → DESeq2 VST + SVA residualization (R)
  → gene panel selection (2,500 protein-coding, biotype-filtered, anchor force-include)
  → cell-standardized shared skeleton (Ledoit-Wolf shrinkage → partial correlation → top-k=80)
  → LIONESS sample-specific edge weights → limma edge-wise regression → predicted networks
  → PecanPy node2vec embeddings (multi-seed, signed pos/neg channels)
  → orthogonal Procrustes alignment on configured anchors → cosine-distance rewiring
  → "silent shifters" (high rewiring ∧ |log2FC| < 0.3 ∧ DE FDR > 0.2)
  → permutation + bootstrap inference → GO/KEGG/Reactome grounding
```

Implemented in `src/preprocessing/`, `src/networks/` (16 modules, 5,976 LOC), `src/statistics/`, `src/enrichment/`, orchestrated by `src/run_all_phases.py` (1,413 lines, phases {0,1,2,3,5,6,7}). Manual intervention: R steps run separately; anchors curated in `config/anchor_genes.yaml`.

Where it stops: the leakage-safe classifier validation (`src/validation/enhanced_cv.py`, all transforms inside folds) **returned negative** — network features did not beat expression baselines. `METHODOLOGY.md` states this plainly: "The audited run is negative." The retrospective (`latex_paper/network_journey_retrospective.tex`, §"why the network never beat differential expression") is the postmortem.

### Workflow B — OSD-462 Stage 0 provenance audit (July 2026). Complete.

```
OSD-462 ISA metadata + MS workbooks
  → sample/plex/reporter-tag design table (src/multiomics/osd462_stage0.py, 1,110 lines)
  → assay-label reconciliation (protocol TMTpro +304.207 Da vs legacy field "iTRAQ")
  → reporter-position map (baseline 126–128C, flight 129N–131N, ground 131C–133C, both plexes)
  → peptide-level phosphoform qualification (localization score, singly-modified, composite detection)
  → isolated_canonical_assay_features() → count of qualified NCC/SPAK phosphoforms
Output: run_20260728_osd462_stage0/, figures/v13/osd462_reporter_tag_map.*
```

Three findings, each verified in the run outputs:
1. Assay chemistry is internally inconsistent between protocol and metadata field.
2. **Condition is perfectly aliased with reporter-tag block in both plexes, with no cross-plex label swap.** A flight coefficient is therefore not separable from a systematic tag-position effect.
3. **Zero** isolated, qualification-passing canonical NCC/SPAK regulatory phosphoforms exist in the workbook. The feature indexed as NCC T53 is a T53/Y65 co-modified peptide; the SPAK S383 feature is S382/S383 co-modified; positions 65 and 68 are tyrosines, not the canonical regulatory Thr/Ser.

This audit **retracted the project's own v11/v12 headline.**

### Workflow C — v13 exact continuous phospho inference (July 2026). Complete and locked.

```
frozen config (config/dct_asdn_phospho_reanalysis.yaml)
  + flight-blind gene sets (config/dct_subtype_reference_freeze_v1.yaml, SHA-256 ed64f0e0…)
  → label-blind site universe (loc score ≥ 19, singly modified, no composites, completeness rules)
     → 8,021 phosphosite features, 3,524 parent genes
  → enumerate_balanced_labels(): C(10,5)² = 63,504 assignments, observed forced to row 0
  → PASS 1 (_first_pass_calibration): per-gene null mean/SD over all assignments, chunked (1,024)
  → PASS 2 (_second_pass_inference): competitive set statistic = mean gene Z in set − mean Z in background
  → maxT family {DCT2/CNT, ASDN}; BH family {compartment comparators}
  → claim gates → claim_tier.tsv
  → src/v13/reporting.py → figures + claim_decision_summary.md
  → finalize_inference_provenance.py → provenance_manifest.json (input/code/output SHA-256 + git state)
```

Nine analysis profiles run in parallel (primary, loc-13, site-position collapse, dedup multi-modified, two broad-expression exclusions, uncentered, summed S/N, parent-protein subtracted). Run directory: `data/results/run_20260729_v13_continuous_phospho_exact_final/` (74 MB, 41 files).

Result: **`claim_tier = neither`.** DCT2/CNT non-evaluable (5 observable genes vs frozen minimum 8); ASDN nominally positive (0.718, exact p = maxT = 0.0291) but failed the selectivity gate because podocyte scored higher (1.112, p = 2.36e-4) and the observability-matched ASDN test was null (p = 0.247). Publication promotion blocked by the design/provenance gate.

Manuscript v13 (`latex_paper/manuscript_v13.tex`, 8 pages) reports exactly this.

### Workflow D — cross-mission clinical renal axes (August 2026). Complete run; paper not written.

```
5 terminal missions (OSD-102, -163, -253, -462, -771), 83 animals (12/12/19/20/20)
  → per-mission VST matrices, symbol mapping, CPM ≥ 0.1 eligibility (src/clinical_axes/data.py)
  → animal-level signed gene z-scores → equal-weight subdomain → axis score (statistics.py)
  → Hedges g per exchangeability stratum → inverse-variance fixed effect within mission
  → REML random effects across missions + modified Hartung–Knapp CI + prediction interval
  → 100,000 blocked whole-pipeline label permutations, max-|T| FWER over the 4 frozen axes
  → adversarial layer: leave-one-mission, leave-one-gene, signed median, common-observable,
     3 RNA preparations, OSD-253 rerun control, OSD-163 mapping-rate residualization,
     10,000 observability-matched panels, 49-set frozen compartment family,
     strict same-tier optimal matching (Hungarian assignment)
Output: run_20260811_clinical_renal_axes_cross_mission/ (94 files),
        figures/clinical_renal_axes_cross_mission/ (3 figures)
```

Primary result, verified directly from `primary_meta_results.tsv` and the decision report:

| Frozen axis | Pooled g (mHK 95% CI) | maxT FWER | I² | Decision |
|---|---:|---:|---:|---|
| Glomerular barrier identity loss | −0.716 (−1.361, −0.072) | 0.0018 | 0.0% | Reject "identity loss"; opposite-direction lead |
| Tubular injury (Havcr1/Lcn2) | 0.570 (−0.095, 1.234) | 0.0263 | 8.0% | No-go — CI crosses zero |
| Fibrosis / maladaptive repair | 0.311 (−0.732, 1.355) | 0.785 | 62.5% | No-go — null, heterogeneous |
| Distal transport identity loss | −0.153 (−1.319, 1.012) | 0.986 | 68.8% | No-go — null, heterogeneous |

Compartment family (verified by reading `compartment_context_meta_results.tsv` directly, 49 sets):

| Rank | Set | g | mHK CI | maxT FWER | I² |
|---:|---|---:|---:|---:|---:|
| 1 | podocyte high_specificity | 0.689 | 0.042 – 1.336 | **0.0189** | 0% |
| 2 | podocyte all_enriched | 0.567 | −0.078 – 1.213 | 0.105 | 0% |
| 3 | podocyte broad_enriched | 0.554 | −0.090 – 1.198 | 0.120 | 0% |
| **4** | **broad_structural_scaffold_control__all** | **0.549** | −0.095 – 1.193 | 0.128 | 67.3% |
| 5 | podocyte scaffold_excluded | 0.537 | −0.108 – 1.181 | 0.149 | 0% |
| 6 | principal_cell broad_enriched | 0.530 | −0.102 – 1.162 | 0.143 | 0% |

Only row 1 clears the family. But the declared structural negative control is row 4 — above the podocyte tier built by *removing* structural genes. This is the open question the project is currently sitting on.

---

## 3. Capability audit

### A. Implemented and demonstrated

- **Exact conditional permutation inference over the full balanced label space** (63,504 assignments), with per-gene null standardization, chunked two-pass computation, nine parallel analysis profiles, and machine-emitted claim gates. Run artifacts and hashes exist.
- **Cross-mission random-effects meta-analysis under 100,000 blocked whole-pipeline permutations** with REML τ², modified Hartung–Knapp intervals, prediction intervals, and max-|T| FWER. Implemented in a batched/vectorized form (`_batch_tau2_reml`, `_batch_reml_mkh` in `src/clinical_axes/statistics.py`) so that the full REML fit is recomputed for every permutation rather than approximated.
- **MS assay-provenance forensics.** Reporter-tag design reconstruction from ISA metadata; peptide-level phosphoform qualification; detection of co-modified peptides and residue mis-indexing.
- **Flight-blind reference construction.** GSE228367 raw-count reconstruction from 10x matrices, replicate-summed DCT1/DCT2, paired edgeR pseudobulk; GSE150338 fine-subtype and microdissection checks; Mouse Kidney Atlas compartment tiers. Built without reading flight effects (`src/subtype_reference/reference_builder.py`, 788 lines).
- **Observability-matched null panels.** Genes matched on abundance, variance, atlas specificity, source-study detection fraction, and non-target compartment breadth, with balanced unique panels and optimal (Hungarian) assignment (`src/clinical_axes/matching.py`).
- **Provenance manifesting.** SHA-256 over inputs, code files, every output, plus git commit and porcelain-status hash.
- **Reproducibility harness.** 225 tests, executed clean at audit time.
- **Manuscript production.** 13 numbered versions plus supplements; v13 compiles to an 8-page PDF.

### B. Implemented but preliminary or approximate

- **The podocyte cross-mission result.** Real, adversarially tested eleven ways, and it survives all of them — but the effect is moderate, the mHK CI barely excludes zero (0.042–1.336), the prediction interval crosses zero, and the declared structural control sits 0.14 below it against CI widths of ~1.3. The project's own §3b amendment states the honest reading: a band of structural/compartment sets drifting positive at g ≈ 0.45–0.57 with the podocyte tier at the top.
- **The compartment lock is retrospective.** `plan_changelog.md` says so explicitly for both the 2026-08-02 and 2026-08-11 locks. These are post-hoc specifications written before *computing new results*, not preregistrations predating discovery — and the documents refuse the stronger word.
- **OSD-462 reporter-position bounding.** Within-block slopes are non-significant across all six block × plex estimates (pooled +0.0072 log2/step, predicting −0.036 vs observed −0.179/−0.157). This bounds a *linear* positional trend at ~20–23% of the effect. A step-shaped block effect remains invisible to it. Position is also confounded with animal-ID fill order. The project labels this "bounded, not resolved."

### C. Partially implemented

- **Cell-type deconvolution.** MuSiC pipeline, hybrid reference build, DCT_early/DCT_late splitting, and three generations of sanity checks (`data/results/deconvolution_sanity_check*.tsv`) exist. Composition is used as a residualization covariate and in composition-aware phospho robustness tests, but deconvolved proportions are not the basis of any surviving claim.
- **`scripts/run_full_pipeline.py`** — 157 lines, described in `FILE_SUMMARY.md` as "Legacy placeholder (not currently used)." Accurate.
- **`src/v11/cmap_screen.py`** (LINCS/CMap) and **`kinome_atlas_ksea.py`** — implemented and tested, but explicitly labelled appendix-tier / hypothesis-generating; they feed no current claim.
- **OSD-656 human urine context** (`src/clinical_axes/osd656_context.py`) — implemented, run, and deliberately non-inferential by design; the module docstring forbids using it as validation.

### D. Currently under active implementation

Nothing is mid-edit. The last file modification is 2026-08-11 16:27. The frontier is defined but paused:

- **K0** — regress the high-specificity podocyte score on flight adjusting for `broad_structural_scaffold_control__all`, blocked by mission, HC3, same permutation scheme, plus the reciprocal. **Not implemented.** The nearest machinery, `scripts/clinical_axes/run_barrier_specificity_adjustment.py`, hardcodes barrier-core as outcome and the disjoint podocyte tier as proxy (lines 68–83); K0 requires generalizing outcome/adjuster selection. This is a small change to an existing, working script.
- **The podocyte manuscript.** Fully blueprinted (`docs/PODOCYTE_CROSS_MISSION_PAPER_BLUEPRINT_2026-08-11.md` — title, one-sentence paper, evidence tables, prior-art position). No `.tex` exists.

### E. Planned or documented but not implemented

- Glomerular morphometry validation (WT1⁺ nuclei per glomerulus and per tuft area; nephrin/podocin/synaptopodin staining) — specified precisely, routed as a collaboration ask through the NASA PI, contingent on external tissue and an external group.
- Repository hygiene for reviewer verification: committing run directories, dropping the `latex_paper/*` gitignore exclusion, tagging frozen configs. Scheduled at weeks 7–9 in the plan; not done.
- Forward-citation searches from Melica 2024; the Siew/da Silveira ones are done.
- Standalone publication of the retrospective, sequenced *after* the primary paper.

Explicitly retired in `docs/POST_V13_DIRECTION_DECISION_2026-08-11.md` §9: no v14, no further OSD-462 compartment mining (2026-08-02 stop rule), no Grey60 rescue, no LayerScore method paper, no OSD-462 data note.

---

## 4. What was actually built

**The exact permutation inference engine** (`src/v13/continuous_phospho_inference.py`). Designing the estimand is the substantive part: replacing a nominal-p phosphosite contingency endpoint with one continuous score per parent gene (`median(−site effect)`), standardized against a gene-specific null that preserves that gene's own site count, missingness pattern, and within-protein site dependence. Then making it computable: two passes, chunked at 1,024 assignments, with the observed assignment pinned to row 0 and verified by assertion.

**The cross-mission meta-analytic estimator** (`src/clinical_axes/`, 2,274 LOC across 6 modules). Choosing one biological animal as the unit and one estimate per mission — rather than treating correlated genes as inferential units — is the design decision that retired the project's own earlier `p = 7e-4` matrix/ECM result. The batched REML implementation exists because a whole-pipeline permutation requires refitting τ² 100,000 times.

**The Stage-0 forensic audit.** This is the single highest-value piece of work in the repository and it is destructive of the author's own prior claim. Finding that the T53-indexed feature is co-modified, that positions 65/68 are tyrosines, and that condition is perfectly aliased with reporter-tag block required reading the ISA metadata and MS workbooks at peptide level rather than trusting summary tables.

**Flight-blind reference construction.** Rebuilding DCT1/DCT2 references from GSE228367 raw 10x matrices with replicate-summed pseudobulk and paired edgeR, plus GSE150338 microdissection validation, specifically so that the marker sets could not have been influenced by the flight effects they would later be tested against.

**Observability-matched and same-tier-matched nulls.** Two generations: 10,000 six-gene panels matched on RNA observability/expression/variability/breadth; then a stricter audit where every comparator was drawn from the same frozen `high_specificity` atlas rule and matched on seven variables including atlas specificity ratio and detection breadth, with common-support trimming and optimal unique assignment.

**The claim-gate architecture.** Encoding go/no-go criteria in frozen YAML and emitting the tier as a computed artifact is a design decision that removes the author's discretion at the point where discretion is most dangerous.

**The network pipeline** (Workflow A) — skeleton construction, LIONESS, edge regression, node2vec/Procrustes rewiring, silent shifters, leakage-safe CV. Built, audited, found negative, documented as negative.

**Division of labor.** The repository is single-author. `agents_instruction.md` (34 KB), `docs/plan_changelog.md`, and `docs/owner_decisions.md` show the working pattern: Ibrahim writes the specification, estimand, guardrails, thresholds, and acceptance criteria; AI coding agents implement against that specification; decisions are recorded as owner decisions. One external interaction is visible in git — PR #1 (`copilot/fix-marker-validation-issue`, 2026-02-17), a marker-validation and two-stage deconvolution change, merged after code review. A NASA PI is referenced as available for tissue access and publishing help; no PI-authored artifact appears in the repository.

---

## 5. Technical depth

**Enumerating the exact null instead of sampling it.** With 5 flight and 5 ground animals in each of two plexes, the balanced label space is C(10,5)² = 63,504 — small enough to enumerate, large enough that naively recomputing site contrasts → gene aggregation → set statistics 63,504 times over 8,021 sites × 3,524 genes is expensive. The solution is a two-pass chunked design where pass 1 accumulates per-gene sums and sums-of-squares to get the null moments and pass 2 reuses them. *Status: complete and correct as implemented; validated by contract tests in `tests/test_v13_corrected_engine_contracts.py`.*

**Gene-specific null calibration under structured missingness.** Phosphosites are unevenly observed across genes, and sites on the same protein are dependent. Standardizing each gene against its own permutation null — rather than a pooled null — preserves site count, missingness pattern, and dependence structure. *Status: complete; this is the methodological core of v13.*

**Refitting REML inside a permutation loop.** A "whole-pipeline" permutation means the mission-level Hedges g, the inverse-variance within-mission pooling, the REML τ², and the Hartung–Knapp adjustment must all be recomputed for each of 100,000 label draws. Done naively this is intractable in Python; `_batch_reml_score` / `_batch_tau2_reml` / `_batch_reml_mkh` vectorize the τ² root-find across permutations. *Status: complete.*

**Non-identifiability that no estimator can fix.** OSD-462's condition is perfectly aliased with reporter-tag block. The project's response — bounding the linear component, showing the phospho-vs-protein layer asymmetry, and then declaring the claim blocked regardless of the p-value — is the correct one. *Status: correctly diagnosed, correctly declared unresolvable with these data. Genuinely unresolved and unresolvable without a label-swapped experiment.*

**Separating a compartment signal from a protein-class signal.** If broadly-expressed structural genes drive the effect, "podocyte program" is the wrong name for it. The frozen scaffold control was built for this, and it now ranks 4th of 49. The literal kill criterion ("beats every compartment set") is not triggered; the spirit of it is in question. *Status: **unresolved**. This is the K0 gap, and the project has correctly identified it as the load-bearing open question rather than declaring victory on the FWER.*

**Composition vs. per-cell transcription in bulk RNA.** Higher podocyte-marker abundance can mean more podocytes, different tissue sampling, or more transcription per podocyte. *Status: acknowledged as unresolvable in bulk RNA; correctly deferred to morphometry.*

**Debugging a pipeline that kept producing attractive results.** The commit history and audit documents record the pattern: `docs/statistical_assessment_v5.md` (25 KB) and `annotated_issues_v5.md` produced the May 2026 remediation, which removed post-hoc focused permutation, replaced unsigned Stouffer aggregation as primary inference, moved all transforms inside CV folds, tightened the silent-shifter definition to require contrast-matched DE, and forbade selecting Procrustes anchors from observed low-rewiring genes (which would have been circular). Each of these is a self-inflicted correction that reduced the strength of an existing claim.

---

## 6. Technical maturity

The right question for this project: **could its analyses support a defensible scientific conclusion, and how far can its outputs be carried?**

**Yes, for negative and boundary claims — at manuscript strength.** The v13 manuscript exists, compiles, and states a conclusion that the run artifacts support exactly: OSD-462 does not support a selective extension of phosphosite suppression into DCT2/CNT or ASDN, and its design and phosphoform provenance prevent attribution to flight or isolated canonical NCC/SPAK regulation. Every number in the abstract traces to a hashed table in the locked run directory.

**Provisionally, for the one positive claim.** The cross-mission podocyte result has more adversarial testing behind it than most published spaceflight-omics results: five missions, animal-level units, blocked whole-pipeline permutation with FWER control, a 49-set frozen comparator family, two generations of matched nulls, and leave-one-out at both mission and gene level. What it lacks is (a) K0, and (b) any measurement that separates cell number from per-cell expression. It is a strong hypothesis at manuscript-draft strength, and the blueprint is honest that the compartment scan followed the targeted result.

**No, for mechanism or physiology.** The README's "Claim boundaries" section lists seven things the repository does not support, including transporter flux, electrolyte, aldosterone, renin, and blood-pressure conclusions, and cell-of-origin localization from whole-kidney phosphoproteomics. That list is accurate.

**Reproducibility, with one real caveat.** Code, configs, tests, and provenance manifests are in good shape — 225 tests pass from a clean checkout, and every locked run carries input/code/output hashes. But `.gitignore` excludes `data/results/`, `data/external/`, `data/processed/`, and all of `latex_paper/*` except v11; the last commit is 2026-07-17 while the current analyses date to 2026-08-11. **The claim that configurations were frozen before testing is currently unverifiable from git history.** The project's own plan flags this as "an afternoon of work" and schedules it at weeks 7–9.

---

## 7. Intended system, direction, and what "done" means

### A. What exists now, as part of the intended endpoint

The cross-mission estimator (`src/clinical_axes/`), the frozen 4-axis family, the 49-set compartment family, the adversarial suite, the three figures, and the decision report are not preparation for the intended paper — they *are* its results section. The blueprint says the decision report "is most of a draft." The Stage-0 material is also load-bearing in the new context: it justifies why OSD-462 enters the RNA synthesis on total-RNA with mRNA/UPX as within-animal sensitivities, and why its phosphoproteomics is excluded from the evidence base.

### B. What is being built

One test (K0 + its reciprocal) and one manuscript. K0 decides the title, not whether the paper exists — the blueprint carries two headline branches that share introduction, methods, eligibility gates, four-axis null results, and boundary statements:

- **If K0 survives:** *A recurrent podocyte-associated kidney transcript program across five mouse spaceflight missions.*
- **If K0 collapses:** *Most reported renal responses to spaceflight do not recur across independent mouse missions* — with the structural drift reported as the one positive and explicitly labelled non-attributable.

### C. What remains

1. K0 implementation and run (small — generalize an existing script).
2. Manuscript drafting (~4 weeks planned).
3. Repository verifiability work: commit run directories, un-ignore `latex_paper/`, tag frozen configs.
4. Remaining forward-citation search (Melica 2024).
5. PI contact and the morphometry collaboration ask.
6. Internal review → bioRxiv → journal submission.

### Definition of "done"

A researcher opens the preprint and finds: a mission-level effect estimate with a confidence interval and a prediction interval for each of four clinically anchored renal programs across five independent mouse spaceflight missions, computed from one score per animal and one estimate per mission; the eligibility gates (5′/3′ coverage, OSD-253's confounded arm, OSD-462's reporter alias and preparation instability) that justify which cohorts entered; one positive compartment-level result reported alongside the frozen structural control that nearly matches it; the specific boundary statement that bulk RNA cannot separate cell representation from per-cell transcription; and a precisely specified morphometry experiment (WT1⁺ nuclei per glomerulus and per tuft area, plus nephrin/podocin/synaptopodin) that would resolve it. They can clone the repository, run `pytest`, re-execute the locked run from the frozen config, and reproduce every number by hash.

Concretely, "done" means someone can answer a question they cannot answer today: **which reported renal spaceflight effects survive when each mission contributes one estimate, and how large an effect can be excluded for the ones that do not?** No existing work provides this — Siew 2024 pools across missions without mission-level effect-size synthesis; da Silveira 2025 compares two studies by strain.

### Is the current architecture consistent with that endpoint?

Yes, and this is a change from the project's earlier state. Workflow A (networks) was a dead end and is architecturally irrelevant to the endpoint — it survives as a retired module tree and a retrospective. Workflows B–D are all load-bearing. `src/clinical_axes/statistics.py` is deliberately free of repository paths, cohort names, gene panels, and biological interpretation, which is what makes it reusable across axis families; the compartment family, the strict-matching audit, and (with a small change) K0 all plug into the same primitives. Reaching the endpoint requires *adding* one test and writing prose, not redesigning anything.

One structural tension is worth naming: the project's own diagnosis is that across v1–v13 the **dataset was fixed and the question floated** — PPAR → rewiring → age interaction → LAR reversal → multi-omic anchor → DCT1 → DCT2/CNT → ASDN → Grey60 → four clinical axes. Ten hypotheses, one corpus. Workflow D inverts this (fixed question, five cohorts), which is why it is the right frontier. But the podocyte headline was itself found by scanning a 49-set family after the targeted axes failed, so the *newest* result was again produced by the old pattern. The blueprint says this in those words, which is the correct handling.

---

## 8. Where the chain stops

```
scientific question   ██████████ implemented (frozen, four axes + compartment family)
cohort assembly       ██████████ implemented (5 missions, 83 animals, eligibility gates)
preprocessing         ██████████ implemented (VST, symbol mapping, CPM eligibility, QC audit)
estimator             ██████████ implemented (animal → mission → REML/mHK → blocked permutation)
adversarial audit     ██████████ implemented (11 checks + strict matching + K1–K4)
specificity resolution ███░░░░░░░ K0 NOT RUN  ← the chain stops here
manuscript            ░░░░░░░░░░ blueprinted, no .tex
external validation   ░░░░░░░░░░ specified; requires tissue + external group
preprint / publication ░░░░░░░░░░ target Dec 2026 – Jan 2027
```

The parallel phospho chain (Workflows B–C) runs to completion — question → audit → exact inference → claim gate → manuscript v13 — and terminates at a documented boundary rather than a gap. The network chain (Workflow A) terminates at a negative validation result and a retrospective.

---

## 9. Gaps between current state and endpoint

**K0: podocyte-vs-structural-scaffold specificity.**
*Missing:* a regression of the high-specificity podocyte score on flight adjusting for `broad_structural_scaffold_control__all`, blocked by mission, HC3, same permutation scheme, plus the reciprocal.
*Why it matters:* it decides whether the paper's title names podocytes. The declared negative control ranks 4th of 49, and the FWER result answers "is this row unusual within the family," not "is this a podocyte signal rather than a structural one."
*Status:* clearly planned, authorization-blocked, not implemented.
*Where it belongs:* inside the project — generalize `run_barrier_specificity_adjustment.py` to parameterize outcome and adjuster.
*Closure artifact:* a `podocyte_scaffold_adjustment_meta.tsv` giving adjusted and reciprocal-adjusted estimates with mHK intervals, in the same format as the existing `barrier_proxy_adjustment_meta.tsv` (which showed 0.507 → 0.070 for the barrier case).

**Cell representation vs. per-cell transcription.**
*Missing:* any measurement that distinguishes more podocytes from more expression per podocyte.
*Why it matters:* the paper's central biological claim is unresolvable without it; bulk RNA cannot do it in principle.
*Status:* specified in detail, requires external tissue and an external group; not scheduled as a precondition for the preprint.
*Where it belongs:* experimental measurement, outside the project — Sarder's morphometry pipeline (already used in the flagship paper), on sections that already exist.
*Closure artifact:* WT1⁺ nuclei per glomerulus and per tuft area, plus nephrin/podocin/synaptopodin staining, on flight vs. ground sections.

**Reviewer-verifiable freeze history.**
*Missing:* commits for the run directories and configs; `latex_paper/*` is gitignored except v11; last commit predates the current analyses by ~4 weeks.
*Why it matters:* every "frozen before testing" statement in the manuscript is currently unverifiable from repository history, and a reviewer will check.
*Status:* known, scheduled, trivial.
*Closure artifact:* a tagged commit containing `config/clinical_renal_axes_cross_mission.yaml` and the run directory, dated before the results commit.

**The podocyte manuscript.**
*Missing:* the `.tex`.
*Status:* blueprinted in full; both headline branches specified; figures exist.
*Closure artifact:* `latex_paper/manuscript_podocyte_cross_mission.tex` compiling to a preprint PDF.

**Prediction beyond the observed missions.**
*Missing:* nothing implementable — the prediction interval crosses zero and a sixth mission would be required.
*Status:* correctly reported as a limitation rather than treated as a gap.

---

## 10. Validation and evidence

| Mechanism | What it validates | What it does not | Strength |
|---|---|---|---|
| 225 unit/contract tests (executed clean, 2026-08-25) | Statistical primitives, config consistency, engine contracts, fold safety, provenance schema | That the biological conclusion is correct | Strong for code; silent on science |
| Exact enumeration of all 63,504 balanced assignments | The conditional p-value is exact, not Monte Carlo | Exchangeability itself — which the reporter alias breaks | Strong, but conditional on an assumption the project says is violated |
| 100,000 blocked whole-pipeline permutations + maxT | Family-wise error over the frozen axes | Whether the frozen family was the right family | Strong within the family |
| Observability-matched panels (10,000 draws; pools of 20/50/100/200) | The effect is unusual against genes matched on expression/variance/breadth (p = 0.0062) | Whether matching variables capture the relevant nuisance | Good |
| Strict same-tier matching (Hungarian, 7 covariates, 142 targets vs 543 candidates) | The atlas abundance/specificity extreme does not explain the result | Generalization of the target-minus-matched contrast to new missions | Good; the strongest single check |
| Leave-one-mission / leave-one-gene | No single mission or gene carries the effect (all retained direction; Nphs2 = 34.2%, removal still g = −0.524) | Correlated-gene structure across the panel | Good |
| Three RNA preparations of the same 20 animals | Result holds under total-RNA, mRNA, UPX (g = −0.673, −0.686) | Library-prep discordance in general | Good |
| Frozen structural negative control | That a protein-class explanation was *tested* | Which of podocyte vs. structural drift carries the signal — **it ranks 4th of 49** | Currently the weakest link |
| Provenance manifests (SHA-256 inputs/code/outputs/git) | That a given table came from a given code state and input | That the config predates the result — git history does not currently show this | Strong in file, weak in repository |
| Flight-blind reference construction | Marker sets were not selected using the flight effects they test | Whether the atlas itself is right for mouse kidney | Good |
| Prior-art audit (51 Siew-citing papers screened, Crossref/OpenAlex/Semantic Scholar) | No published recurrent podocyte-transcript program across missions | Systematic-review completeness — the project says it is not one | Adequate and honestly scoped |
| Leakage-safe cross-validation (all transforms in-fold) | That network features do **not** beat expression baselines | — | Strong; a validated negative |

**Outputs with little or no independent validation:** everything downstream of bulk RNA compartment attribution. There is no histology, no cell-resolved expression, no protein-level podocyte measurement, no urinary albumin/creatinine, and no independent sixth mission.

**Validation that only becomes relevant later:** morphometry and ultrastructure are meaningless until the K0 question fixes what the claim is about.

---

## 11. Important technical decisions and project evolution

**9 Dec 2025 — the founding proposal.** An ASGSR-style poster/proposal (`node2vec_asgsr_osd771_kidney_transcriptomics.pdf`) predates the git history by six days. It specifies the method concretely — node2vec biased random walks (p = 0.25, q = 4.0, walk length 80, 200 walks/node, window 10, ten seeds), Procrustes alignment on housekeeping/ribosomal anchors, rewiring read off as cosine distance — and names four aims, the fourth of which is validating the embeddings through cross-condition predictive modeling. That fourth aim is the one that later kills the premise; the project promised the test in advance and then ran it.

**Dec 2025 – Feb 2026 — the network premise.** node2vec embeddings on sample-specific networks to find "silent shifters": genes that rewire without changing expression. Deconvolution via MuSiC with a hybrid Tabula Muris Senis + Chen reference; DCT split into early/late (commit `250187f`); DCT residualized for anti-confounding (`b096e3c`). PR #1 added a two-stage deconvolution and marker guards after code review.

**Mar 2026 — first writeup, first cracks** (`latex_paper/main.tex`). The retrospective's §"Anatomy of a silent shifter (and why the definition was ill-posed)" records the core problem: the definition mixed a rewiring statistic with a null DE result without requiring the DE contrast to match the rewiring contrast.

**Apr–May 2026 — the audits land.** `docs/statistical_assessment_v5.md` (25 KB) and `annotated_issues_v5.md` drove commit `35ab40c` ("Remediate methodology and add validation modules"). Effects: post-hoc focused permutation removed entirely; unsigned Stouffer aggregation demoted from primary inference to a signed incident-edge t-statistic calibrated by within-stratum permutation; silent shifters redefined to require contrast-matched DE; Procrustes anchors moved from observed low-rewiring genes (circular) to a config file with a ≥20 mapped-anchor hard failure; ribosomal/translation genes excluded from anchors because translation biology is a result of interest, not a neutral reference.

**May 10 — the cross-validation that would not cooperate.** Fold-safe rerun of the network-vs-expression classifier came back negative. This is the moment the network premise stops being viable, and it is recorded as such rather than buried.

**May 12–20 — three pivots in eight days.** WGCNA module activity (`manuscript_wgcna_pivot.tex`); then the contrast-vector framework (`agents_instruction.md`, locked May 2026) reformulating the question as a projection of the flight aging vector onto the control aging vector with β/cosθ/ρ interpretation categories and five explicit guardrails; then the LAR-reversal detour and the v5 reckoning.

**May 22 – Jul 2026 — the phosphoproteomic pivot.** OSD-462 multi-omics enters (`125e2a2`). v11 builds DCT1 phospho-mediation, kinome KSEA, CMap screening, human concordance, and a proteome observability-bias audit. v12 promotes the DCT2 subtype result to the headline.

**Jul 28–29 — Stage 0 invalidates the anchor.** The provenance audit finds the reporter alias, the assay-label inconsistency, and zero isolated canonical phosphoforms. The v11/v12 headline is retracted by its own author. v13 is rebuilt from scratch as an 8-page provenance/boundary report with a new estimand, and returns `claim_tier = neither`.

**Jul 29 — Grey60 no-go.** The last standing positive RNA result (48-gene ECM/cell-migration module) is put through six gates. Gates B (flight-blind module reconstruction: 12/48 overlap, Jaccard 0.022, projected p = 0.831, 0/27 grid recoveries), D (baseline falsification g = 0.563 vs 0.50 ceiling), and E (external recurrence: meta g = 0.163, CI −0.658 to 0.984) fail. Retired.

**Aug 2 — compartment-wide stop rule.** A final post-hoc audit across ten kidney compartments. Four tier/score rows pass their score-specific family; none passes the complete adversarial gate. Further compartment mining is retired by rule.

**Aug 11 — the reframe.** The cross-mission clinical-axes analysis returns no-go on all four axes plus one podocyte lead. The direction document diagnoses the structural error — dataset fixed, question floating, ten hypotheses on one corpus — and inverts it: fix the question, float the evidence base. Same day, the §3b amendment self-corrects the podocyte enthusiasm after reading `compartment_context_meta_results.tsv` directly and finding the structural control at rank 4.

**Net effect:** the project has retired more claims than it has kept. That is visible in the artifacts, not just asserted.

---

## 12. Failures, limitations, unresolved problems

### Major — affect what the current work can claim

1. **Podocyte vs. structural drift is unresolved** (K0). Rank 4 of 49 for the declared negative control; point estimates differ by 0.14 against CI widths of ~1.3.
2. **OSD-462 condition is perfectly aliased with reporter-tag block.** Non-identifiable from these data. Blocks the entire phospho layer from causal attribution regardless of statistics.
3. **Bulk RNA cannot separate cell representation from per-cell transcription.** Structural, not fixable in software.
4. **The prediction interval crosses zero** (−1.455, 0.023 for the glomerular panel). Coherence across five observed missions is not a guarantee for a sixth.
5. **Adjustment for a disjoint podocyte proxy absorbs the barrier-core effect** — verified directly: 0.507 (CI 0.162–0.854) → 0.070 (CI −0.379–0.533), with I² rising 0% → 72%. The six barrier-core genes do not move beyond the broader program.
6. **Locks are retrospective.** Stated in `plan_changelog.md` for the 2026-08-02 and 2026-08-11 specifications; the compartment scan followed the targeted result.
7. **Freeze history is not verifiable from git.** Last commit 2026-07-17; the August analyses, configs, and run directories are uncommitted or gitignored.
8. **Small n throughout.** 5 mice per design cell in OSD-771; 12–20 animals per mission; five missions in the synthesis; k = 5 for every REML fit.

### Secondary

- Documentation drift: `README.md` reports "199 passed" against an actual 225; `FILE_SUMMARY.md` references `src/utils.py`, which does not exist (`src/common.py` does); `METHODOLOGY.md` carries a scope notice because most of its content describes the retired pipeline.
- `scripts/run_full_pipeline.py` is a 157-line unused placeholder alongside the real orchestrator `src/run_all_phases.py`.
- Repository hygiene: `.DS_Store` files, a stray `.swp`, `__pycache__` directories tracked in `scripts/`, `permutedStats-actualModules.RData` (564 KB) at root, `venv/` (1.2 GB) and `tmp/` (43 MB) in the working tree, and zipped run archives under `data/results/`.
- `LICENSE` is a zero-byte file.
- Environment drift: `requirements.txt` pins Python 3.11 while `environment.yml` specifies Python 3.10; several dependencies are commented out with notes about conflicts.
- Deconvolution has three generations of sanity-check outputs and no single authoritative one.

### Expected incompleteness — not failures

- No podocyte manuscript: drafting has not started, by plan.
- No histology or morphometry: requires external tissue and an external group, and is deliberately not a precondition for the preprint.
- No sixth mission: none exists.
- K0 not run: awaiting owner authorization per §12 of the direction document.

---

## 13. Actual technical ambition

The ambition has changed once, and the change is the substance of the project's second year.

**Originally:** demonstrate that co-expression network rewiring carries renal spaceflight signal that differential expression misses, and deliver a ranked set of "silent shifter" genes as candidate drivers of distal-tubule remodeling. This was an original bet on an unexamined signal class, pursued on a dataset chosen for its unique three-way design. It was tested by the project's own prespecified classifier validation and did not survive.

**Now:** the project is trying to build **a mission-level reproducibility layer for spaceflight renal biology** — an estimator that takes independent rodent missions, produces one effect estimate per mission with a confidence interval and a prediction interval, controls family-wise error across a prespecified set of biological programs, and reports both what recurs and what can be excluded.

That is a different object from what the field currently has. Siew 2024 pools ~11 mouse missions pan-omically but does not synthesize mission-level effect sizes with FWER control; da Silveira 2025 compares two studies by strain. The forward-citation audit found no cross-mission animal-level meta-analysis of spaceflight kidney.

The existing components are the layer: `src/clinical_axes/` is the estimator, the frozen configs are the program definitions, the adversarial suite is the credibility layer, the Stage-0 material is the eligibility layer, and the retired workflows are the specificity control that shows the pipeline returns nulls. When the chain completes, someone can ask "does effect X recur across missions, and if not, how large an effect can be excluded?" and get a numeric answer with an interval.

**How much of that is in the repository today:** the machinery, entirely. The application to four clinical axes and 49 compartments, entirely. The one open specificity question, not yet. The write-up, not yet. The experimental validation, not at all and not intended to be.

**Is development moving toward it?** Yes, and the trajectory corrected in August. The first year fixed the dataset and floated the question. Workflow D fixes the question and floats the evidence base, which is the design the ambition requires. The remaining risk is that the newest positive result was itself found by scanning a 49-set family — the pattern the reframe was meant to end — and the project's own documents say so.

---

## 14. Project scale

| Measure | Value |
|---|---|
| Git-tracked files | 253 (184 `.py`, 23 `.R`, 22 `.md`, 10 `.yaml`) |
| Python in `src/` | 39,443 lines across 16 subpackages, 95 modules |
| Largest modules | `v13/continuous_phospho_inference.py` 3,048; `v11/core_analysis.py` 2,004; `v13/reporting.py` 1,529; `run_all_phases.py` 1,413; `v13/compartment_adversarial_audit.py` 1,314; `multiomics/osd462_stage0.py` 1,110 |
| Analysis scripts | ~80 under `scripts/`, organized into 8 workflow subdirectories |
| R code | 27 scripts, 8,088 lines (deconvolution, DESeq2/limma, WGCNA, pseudobulk references) |
| Tests | 44 files, 225 tests, 5,586 lines — all passing at audit |
| Frozen configurations | 16 YAML files (`config/`) |
| Decision & planning documents | 65 files in `docs/` (~1.2 MB) |
| Manuscript versions | 13 numbered (v1–v13) plus supplements, methods, results compendia, critique reports, posters, a retrospective, and a preregistration — 300 files in `latex_paper/` (~51 MB) |
| Timestamped run directories | 81 under `data/results/` |
| Data volume | 14 GB |
| Commits | 67 (2025-12-15 → 2026-07-17) |
| Primary cohort | 80 mice × 57,278 genes, 2 × 2 × 4 factorial |
| Cross-mission synthesis | 5 missions, 83 animals (12/12/19/20/20) |
| Phospho analysis universe | 8,021 sites → 3,524 parent genes |
| Exact permutations | 63,504 (complete enumeration) |
| Monte Carlo permutations | 100,000 blocked, per axis family |
| Frozen gene-set families | 4 clinical axes; 49-set compartment family; 51-set coherence family; 10 compartments × breadth tiers |
| Matched-null draws | 10,000 (observability); balanced same-tier panels over 142 targets vs 543 candidates |

---

## 15. External dependencies vs. project work

| Dependency | What it does | What the project supplies |
|---|---|---|
| NumPy / pandas / SciPy | Array math, tabular I/O, distributions | The permutation design, chunking strategy, per-gene null calibration, and every statistic computed on top |
| statsmodels | Regression fitting, HC3 covariance | The estimand (animal → mission → meta), the blocking structure, the exchangeability strata |
| — (no meta-analysis library) | — | **REML τ² estimation, modified Hartung–Knapp intervals, prediction intervals, and their batched forms are implemented in-repo** (`src/clinical_axes/statistics.py`) |
| PecanPy | node2vec random walks and embedding training | Graph construction from the skeleton, signed pos/neg channel design, multi-seed protocol, Procrustes alignment and rewiring metric |
| limma / edgeR (via R) | Empirical-Bayes moderation, pseudobulk DE fitting | Design matrices, contrast definitions, covariate policy (`validate_covariate_policy`), replicate-summed pseudobulk construction |
| MuSiC | Cell-type deconvolution solve | Reference assembly (TMS + Chen hybrid), DCT subtype splitting, marker curation, sanity-check design |
| DESeq2 / SVA | VST normalization, surrogate variables | Which covariates are removed, and the policy forbidding double-adjustment downstream |
| WGCNA | Module detection | The 27-specification flight-blind grid, the preservation criteria, the Jaccard/projection gates |
| scikit-learn | Classifiers, PCA, scaling | The fold-safe protocol that puts *every* transform inside folds — the design decision that turned the classifier result negative |
| gseapy | Preranked GSEA | Gene-set curation, direction registry, the `context_detected` vs `replicated` distinction |
| SciPy `linear_sum_assignment` | Hungarian assignment | The matching covariates, common-support trimming, balanced-panel construction |
| Mouse Kidney Atlas / GSE228367 / GSE150338 | Expression reference data | Compartment tier rules, specificity thresholds, the flight-blind construction protocol, the frozen membership file (SHA-256 `ed64f0e0…`) |
| AI coding assistance | Implementation against written specifications | The estimands, guardrails, thresholds, claim gates, negative controls, kill criteria, and every scientific decision — recorded in `agents_instruction.md`, `docs/plan_changelog.md`, and `docs/owner_decisions.md` before implementation |

The pattern throughout: standard tools do the numerical work; the project defines the question, the unit of analysis, the null, the family, the gates, and the boundary.

---

## 16. Factual evidence for later writers

### Technical independence

- The scientific question changed ten times across v1–v13; each change is a documented decision with a dated specification, a config, and a run directory. The current question was written by the author, not assigned.
- Every claim gate, kill criterion, and negative control in the repository was specified by the author before the corresponding run.
- The August reframe (`POST_V13_DIRECTION_DECISION_2026-08-11.md`) diagnoses the project's own structural error — dataset fixed, question floating — and redesigns around it.

### Technical depth

- Exact enumeration of C(10,5)² = 63,504 balanced within-plex label assignments with per-gene null standardization, implemented as a chunked two-pass computation.
- Batched REML τ² estimation with modified Hartung–Knapp intervals recomputed inside 100,000 blocked whole-pipeline permutations.
- Peptide-level phosphoform forensics that identified co-modified peptides and tyrosine mis-indexing behind features labelled as canonical NCC T53 and SPAK S383.
- Detection of a perfect condition-to-reporter-tag-block alias from ISA metadata, followed by a quantitative bound on the linear positional component (pooled slope +0.0072 log2/channel step; predicted −0.036 vs observed −0.179/−0.157) and an explicit statement that a step-shaped effect remains unbounded.
- Seven-covariate strict same-tier matching with optimal unique assignment and common-support trimming.

### Persistence and debugging

- The May 2026 remediation removed post-hoc focused permutation, replaced unsigned Stouffer aggregation as primary inference, moved all CV transforms inside folds, tightened the silent-shifter definition, and eliminated circular Procrustes anchor selection — each reducing the strength of an existing claim.
- The fold-safe classifier rerun (May 10) returned negative and ended the network premise; the result is documented as a negative rather than discarded.
- Grey60 was put through six gates; three failed; it was retired.
- The §3b same-day amendment: after writing a positive-headline direction document, the author read `compartment_context_meta_results.tsv` directly, found the frozen structural control at rank 4 of 49, and amended the document to say the podocyte framing might not survive.

### Research judgment

- Retracted the project's own v11/v12 headline after Stage 0, and rebuilt v13 from scratch as a boundary report.
- Retired the gene-wise Stouffer matrix/ECM `p = 7e-4` result on the grounds that it treated correlated genes as inferential units.
- Refused to claim barrier failure despite a significant glomerular result, because astronaut literature reports *reduced* urinary albumin excretion in flight.
- Refused to frame the direction disagreement with a ~100-author consortium paper as a correction, after reading their Fig. 5 and Supplementary Data 3 directly and finding the claim was neither urine-only nor cross-species (K1).
- Declared the OSD-462 phospho claim blocked by design regardless of its p-value.
- Wrote a stop rule (2026-08-02) retiring further compartment mining, and a "stays retired" list (§9 of the direction document) enumerating seven abandoned directions.

### Initiative

- Assembled an external evidence base spanning six OSDR missions, three single-cell atlases, two spatial references, two kinome atlases, LINCS/CMap, and human urine cohorts — none of which was part of the original single-dataset scope.
- Built a cross-mission meta-analytic estimator that, by the project's own forward-citation audit, does not exist in the published spaceflight kidney literature.
- Wrote `network_journey_retrospective.tex` — a dated trail of claims made, tested, and retracted by the same author — and sequenced it deliberately *after* the primary paper.

### Collaboration and mentorship

- Single-author repository. One merged external PR (`copilot/fix-marker-validation-issue`, 2026-02-17) addressing marker validation and two-stage deconvolution, with code-review responses.
- A NASA PI is referenced as available for tissue access, suggestions, and publishing help; no PI-authored artifact is present.
- The identified validation path routes through the existing Siew consortium (Sarder's morphometry pipeline, Roufosse's pathology, Walker-Samuel's imaging) on sections already cut, rather than a new tissue request.

### Tangible outputs

13 numbered manuscripts (v13 compiles to an 8-page PDF); a preregistration document; a poster; 81 timestamped run directories with hashed provenance manifests; 3 publication figures for the cross-mission analysis plus the v13 reporter-tag map and two v13 result figures; 16 frozen configurations; 225 passing tests; a frozen gene-set membership file with a published SHA-256; ~39 k lines of Python and ~8 k lines of R; 65 decision and audit documents.

---

## 17. Factual characterization for later use

1. Began the project from a self-identified gap rather than an assigned or derivative question: designed and proposed (ASGSR, 9 Dec 2025) an original study measuring co-expression network *rewiring* — regulatory reorganization invisible to differential expression — in the only mouse kidney dataset crossing age, mission arm, and flight in one factorial design, and named the target category ("silent shifters").
2. Built a cross-mission meta-analytic estimator for mouse spaceflight kidney RNA that computes one signed score per animal, one Hedges g per mission, REML random-effects pooling with modified Hartung–Knapp intervals, and family-wise error control over 100,000 blocked whole-pipeline label permutations — applied to five missions and 83 animals across four prespecified clinical renal programs and a 49-set kidney-compartment family.
3. Built an exact conditional permutation engine for isobaric-tag phosphoproteomics that enumerates all 63,504 balanced within-plex label assignments and standardizes each of 3,524 parent genes against its own null, preserving site count, missingness, and within-protein dependence.
4. Reconstructed OSD-462's assay provenance from ISA metadata and MS workbooks, establishing that condition is perfectly aliased with reporter-tag block in both plexes and that zero isolated, qualification-passing canonical NCC/SPAK phosphoforms exist — and retracted the project's own prior headline on that basis.
5. Built the RRRM-2 network-rewiring pipeline (deconvolution → residualization → shrinkage skeleton → LIONESS → edge regression → node2vec → Procrustes rewiring → permutation inference), then designed the leakage-safe cross-validation that showed network features do not outperform expression baselines, and retired the approach.
6. Constructed flight-blind kidney subtype and compartment references from GSE228367 raw 10x matrices, GSE150338 microdissection data, and the Mouse Kidney Atlas, frozen by SHA-256 before any flight-effect testing.
7. Designed a claim-gate architecture in which go/no-go criteria are frozen in YAML and the decision tier is emitted as a computed artifact — currently `neither` for the v13 phospho analysis and no-go for all four cross-mission clinical axes.
8. Identified a recurrent podocyte-associated bulk-kidney transcript program elevated in four of five terminal spaceflight missions (g = 0.689, mHK CI 0.042–1.336, maxT FWER = 0.0189, I² = 0%), and independently identified that a frozen broadly-expressed structural control ranks 4th of 49 in the same family — the specificity test that would resolve it is specified and not yet run.
9. Wrote 13 manuscript versions, of which v13 is a completed 8-page provenance-and-boundary report whose every reported number traces to a hash-manifested run directory.
10. Currently extending the cross-mission layer into a preprint targeted for December 2026 – January 2027, with two prespecified headline branches selected by one outstanding specificity test.

---

## 18. Canonical project fact sheet

**Project** — RRRM2_Kidney_Transcriptome

**Project type** — Computational research: statistical genomics, with substantial custom statistical-methods implementation, data-provenance forensics, and scientific writing. Not a software product; not an experimental or hardware project.

**Origin** — Not a reanalysis. Began 9 Dec 2025 from an ASGSR proposal predating the repository, testing an original hypothesis — that spaceflight-associated renal regulation includes co-expression network rewiring invisible to differential expression — on RRRM-2/OSD-771, chosen as the only mouse kidney dataset crossing age, mission arm, and flight in one factorial design. The target category, "silent shifters," was named by the project. Reanalysis framing enters only in July 2026, and covers the OSD-462 phospho work specifically.

**Original objective** — Quantify per-gene co-expression rewiring between flight and control, determine whether flight amplifies, dampens, reverses, or redirects the ordinary aging trajectory, and prioritize silent shifters as candidate drivers of distal-tubule remodeling. Retired after the project's own prespecified classifier validation returned negative.

**Current objective** — Determine which reported renal responses to spaceflight actually recur across independent mouse missions when the unit of analysis is one animal and one mission, and quantify what can be excluded for those that do not.

**Current state** — Three analysis workflows complete and locked (Stage-0 provenance audit, v13 exact phospho inference, cross-mission clinical axes); one retired with a documented negative result (network rewiring). One manuscript published-ready (v13, 8 pages, compiled). One specificity test (K0) outstanding, blocking a title decision on the next manuscript, which is blueprinted but undrafted. Last file activity 2026-08-11; last commit 2026-07-17.

**What currently exists** — `src/clinical_axes/` (cross-mission estimator, 2,274 LOC); `src/v13/` (exact permutation engine + reporting + compartment audit, 5,914 LOC); `src/multiomics/osd462_stage0.py` (assay forensics, 1,110 LOC); `src/subtype_reference/` (flight-blind references); `src/v11/` (9,595 LOC of superseded phospho analyses, tests still green); `src/networks/`, `src/statistics/`, `src/validation/`, `src/enrichment/`, `src/preprocessing/` (retired RNA pipeline, 39 k LOC total in `src/`); 16 frozen YAML configs; 225 passing tests; 81 run directories with SHA-256 provenance manifests; 13 manuscripts; 65 decision documents; 3 cross-mission figures.

**What is currently under implementation** — Nothing mid-edit. Frontier defined and paused: K0 (podocyte score adjusted for `broad_structural_scaffold_control__all`, plus reciprocal) and the podocyte cross-mission manuscript.

**Intended completed system** — A published mission-level reproducibility layer for spaceflight renal biology: frozen program definitions, eligibility gates, an animal→mission→meta estimator with FWER control, an adversarial suite, and a precisely specified validation experiment — reproducible from the repository by hash.

**Definition of "done"** — A preprint reporting, for four clinically anchored renal programs across five independent mouse missions, one mission-level effect estimate each with confidence and prediction intervals; the eligibility gates that determined cohort inclusion; one positive compartment-level result reported alongside the frozen structural control that nearly matches it; the explicit boundary that bulk RNA cannot separate cell representation from per-cell transcription; and a morphometry experiment specified precisely enough to execute. Reproducible: clone → `pytest` → re-run from frozen config → match hashes.

**Inputs** — OSD-771/GLDS-674 (n = 80, 57,278 genes, 2×2×4 factorial); OSD-102/-163/-253/-462/-513 bulk RNA-seq; OSD-462 two 15-plex TMTpro phospho/proteomics (20 mice); Mouse Kidney Atlas, Tabula Muris Senis, Chen atlas, GSE228367, GSE150338; GSE269622/GSE269719 spatial; Johnson 2023 + Yaron-Barir 2024 kinome atlases; LINCS/CMap; OSD-656 human urine.

**Core process** — Frozen YAML specification → label-blind eligibility and scoring → animal-level signed scores → mission-level Hedges g → REML/mHK meta-analysis (RNA) *or* exact enumeration of all balanced within-plex assignments with per-gene null standardization (phospho) → blocked/exact permutation with maxT or BH → adversarial sensitivity suite → machine-emitted claim gates → figures, reports, and a hashed provenance manifest.

**Outputs** — `claim_tier.tsv` / `claim_gates.tsv`; per-axis and per-set meta results with intervals, τ², I², prediction intervals, and FWER; permutation nulls; matched-panel nulls; leave-one-out tables; coverage and attrition audits; publication figures (PDF/PNG/SVG); `provenance_manifest.json`; decision documents; LaTeX manuscripts.

**Major components** — (1) cross-mission clinical-axis estimator; (2) exact phospho permutation engine; (3) Stage-0 assay-provenance auditor; (4) flight-blind reference builder; (5) adversarial/matching audit suite; (6) reporting and provenance layer; (7) retired RNA network-rewiring pipeline; (8) manuscript corpus.

**Key technologies and methods** — Hedges g; REML random effects with modified Hartung–Knapp; prediction intervals; blocked whole-pipeline label permutation with max-|T| FWER and BH; exact conditional permutation by complete enumeration; per-gene null standardization; competitive gene-set statistics; observability- and specificity-matched null panels with Hungarian assignment; Ledoit-Wolf shrinkage partial correlation; LIONESS; limma/edgeR; node2vec (PecanPy) with orthogonal Procrustes; WGCNA; MuSiC deconvolution; DESeq2 VST + SVA; preranked GSEA; KSEA; SHA-256 provenance manifesting.

**My major technical contributions** — Estimand design (unit of analysis, null construction, family definitions); the exact enumeration engine and its two-pass chunked implementation; the batched REML/mHK permutation machinery; the assay-provenance forensics that retracted the project's own headline; flight-blind reference construction; the claim-gate architecture and its frozen negative controls; the matched-null designs; the leakage-safe CV that produced a validated negative; the full manuscript corpus and decision record.

**Strongest demonstrated capabilities** — Exact conditional inference over a complete balanced label space; whole-pipeline permutation meta-analysis at 100 k draws with REML refit; detecting and quantifying a design confound (condition ↔ reporter-tag block) that invalidates an entire measurement layer; running a positive result through eleven independent adversarial checks and reporting where it weakens.

**Major current limitations** — Podocyte-vs-structural specificity unresolved (K0); OSD-462 condition non-identifiable from reporter-tag block; bulk RNA cannot separate composition from per-cell transcription; prediction interval crosses zero; barrier-core effect absorbed by the disjoint podocyte proxy (0.507 → 0.070); locks are retrospective; freeze history not verifiable from git; small n throughout (k = 5 missions, 5 mice per design cell).

**Major unfinished components** — K0 test; podocyte manuscript; repository verifiability commit; morphometry validation (external); one remaining forward-citation search.

**Validation performed** — 225 tests passing (executed 2026-08-25); exact enumeration; 100 k blocked permutations with FWER; 10 k observability-matched panels (p = 0.0062); strict same-tier matching (142 targets, 543 candidates, 7 covariates); leave-one-mission and leave-one-gene; three RNA preparations of the same animals; OSD-253 rerun control; OSD-163 mapping-rate residualization; 49-set frozen compartment family; SHA-256 provenance manifests; prior-art audit over 51 forward citations.

**Current development frontier** — K0 and the two-branch manuscript. Paused since 2026-08-11 pending owner authorization.

**Intended endpoint** — bioRxiv preprint December 2026, journal submission January 2027 (npj Microgravity or a nephrology venue); the retrospective published separately afterward; morphometry pursued as collaboration, not as a submission precondition.

**Most important next technical steps** — (1) Generalize `run_barrier_specificity_adjustment.py` to parameterize outcome and adjuster; run K0 and its reciprocal. (2) Draft the sections common to both branches (introduction, eligibility, methods, four-axis nulls, boundary) — they do not depend on K0. (3) Commit run directories and configs, un-ignore `latex_paper/`, tag the freeze. (4) Contact the PI regarding morphometry on existing sections.

**Tangible artifacts and results** — `latex_paper/manuscript_v13.pdf` (8 pp., compiled); `figures/clinical_renal_axes_cross_mission/` (3 figures); `figures/v13/osd462_reporter_tag_map.*`; `data/results/run_20260729_v13_continuous_phospho_exact_final/` (41 files, 74 MB, hash-manifested); `data/results/run_20260811_clinical_renal_axes_cross_mission/` (94 files); `data/results/run_20260728_osd462_stage0/`; `data/results/run_20260729_grey60_adversarial/`; frozen gene-set membership SHA-256 `ed64f0e20361d21ebc4d2483fce383f5f82e94bb7ec43b66156826c8eb86b945`; `latex_paper/network_journey_retrospective.pdf`; `docs/layerscore_design_preregistration_2026-07-19.pdf`.

**Useful quantitative context** — 39,443 LOC Python in `src/` + ~8,000 LOC R + 5,586 LOC tests (225 tests); 253 tracked files; 67 commits over 8 months; 81 run directories; 14 GB data; 16 frozen configs; 65 decision documents; 13 manuscript versions; 5 missions / 83 animals in the primary synthesis; 63,504 exact permutations; 100,000 blocked permutations; 8,021 phosphosites → 3,524 parent genes; 49-set compartment family.

---

## 19. Concise descriptions

### One sentence

RRRM2_Kidney_Transcriptome began as an original investigation into whether spaceflight reorganizes kidney gene co-regulation in ways differential expression cannot see, and — after that premise failed its own validation and an assay-provenance audit invalidated the project's phosphoproteomic anchor — is now a mission-level reproducibility layer under development, built on an exact permutation engine for isobaric-tag phosphoproteomics and an animal-level cross-mission meta-analytic estimator, that determines which reported renal spaceflight effects recur across independent mouse missions and how large an effect can be excluded for those that do not.

### Three sentences

The project set out in December 2025 to test an original hypothesis — that spaceflight-associated renal biology includes co-expression network rewiring invisible to differential expression — on the only mouse kidney dataset crossing age, mission arm, and flight in one design, and is now intended to become a published mission-level reproducibility layer for spaceflight renal biology: frozen program definitions, cohort eligibility gates, an estimator that yields one effect size per mission with confidence and prediction intervals, family-wise error control, and a specified validation experiment. What exists is the machinery and its application — an exact conditional permutation engine over all 63,504 balanced label assignments for OSD-462 phosphoproteomics, a cross-mission estimator applied to five missions and 83 animals across four clinical renal programs and a 49-set compartment family, a Stage-0 assay audit that established the OSD-462 flight contrast is non-identifiable from reporter-tag position, flight-blind reference construction, 225 passing tests, hashed provenance manifests, and thirteen manuscript versions of which v13 is a completed 8-page boundary report. What remains is one specificity test that decides whether the surviving positive result is a podocyte signal or the upper tail of a broad structural drift, the manuscript itself, a repository commit that makes the config freeze verifiable to reviewers, and glomerular morphometry that only an external group with tissue can perform.

### Technical paragraph

RRRM2_Kidney_Transcriptome began as an original study of whether spaceflight reorganizes kidney gene co-regulation — building a pipeline from cell-type deconvolution through sample-specific co-expression networks, node2vec embeddings, and Procrustes-aligned rewiring to find genes whose network context shifts without changing expression. That premise was retired when the project's own leakage-safe cross-validation showed network features did not outperform expression baselines. The work now determines which reported renal responses to spaceflight recur when each mission contributes a single effect estimate. Three analysis layers are complete. A Stage-0 forensic audit reconstructed OSD-462's sample, reporter-channel, and phosphoform provenance from ISA metadata and MS workbooks, establishing that condition is perfectly aliased with reporter-tag block in both 15-plexes and that no isolated, qualification-passing canonical NCC or SPAK phosphoform exists — which retracted the project's own v11/v12 headline. A v13 inference engine then enumerated all 63,504 balanced within-plex label assignments over 8,021 phosphosites and 3,524 parent genes, standardizing each gene against its own null and emitting a machine-computed claim tier of `neither`. A cross-mission estimator computes animal-level signed scores, mission-level Hedges g, REML pooling with modified Hartung–Knapp intervals, and 100,000 blocked whole-pipeline permutations with max-|T| control; applied to five missions and 83 animals, it returned no-go on all four frozen clinical axes and one surviving positive — a high-specificity podocyte transcript program elevated in four of five missions (g = 0.689, FWER = 0.0189, I² = 0%). An earlier network-rewiring pipeline was built, audited, found not to outperform expression baselines, and retired. Active development is one specificity test — whether the podocyte result survives adjustment for a frozen broadly-expressed structural control that ranks 4th of 49 in the same family — followed by a preprint targeted for December 2026.

---

*Audit produced 2026-08-25 by direct repository inspection. Test suite executed during the audit (225 passed). All quantitative claims above were read from run outputs, configs, source, or git history rather than from README or planning prose; where a document and the artifacts disagree, the artifacts are reported.*
