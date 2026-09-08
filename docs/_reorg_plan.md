# File reorganization plan

**Written:** 2026-08-28. **Baseline:** `26aad1a`, working tree clean at session
start, `venv/bin/python -m pytest -q` → **225 passed, 25 warnings in 24.08s**
(re-verified during this analysis, not quoted from README — README's "199
passed" is stale).

This document is a plan only. Nothing in the repository was moved, deleted, or
edited to produce it. Every claim below is followed by the command that
produced it.

---

## 0. Live-concurrency warning — read this before executing anything

The working tree was clean at 19:20. By 19:29 `git status --porcelain` showed:

```
?? config/revisions/
?? docs/_commentary_scripts.md
?? docs/_commentary_src.md
?? pyproject.toml
?? rrrm2/
?? tests/test_rrrm2_cli.py
?? tests/test_rrrm2_registry.py
```

Another wave is building an experiment registry and CLI **right now**. Three
consequences bind this plan:

1. **`config/revisions/*.yaml` is a machine-checked path registry.**
   `tests/test_rrrm2_registry.py:38-40` is
   `assert stage.entry_path.exists(), f"{rev_id}/{stage.id}: missing {stage.entry}"`
   over every `entry:` in every revision file, and `:45-50` asserts every
   `spec_config:` and every `docs:` entry exists. Moving a registered script or
   a registered doc turns the suite **red**.
2. **`docs/_commentary_scripts.md` and `docs/_commentary_src.md` are path +
   exact-line-range catalogues** ("Line ranges are inclusive and exact… A later
   wave applies the replacements"). Every `git mv` in this plan invalidates the
   keys of that wave's worklist.
3. **`rrrm2/paths.py` discovers the repo root by**
   `_MARKERS = ("config/revisions", "src", "scripts")`. `src/` and `scripts/`
   must continue to exist as root directories. This plan keeps both.

**Sequencing recommendation:** let the commentary wave land first, then execute
this plan, then let the commentary wave re-key. If that is not possible,
execute Batch 0 only (deletes; touches no `.py`/`.R` file the commentary
catalogues).

---

## 1. Summary counts

| Category | Count |
|---|---:|
| Files proposed to **move** | **107** |
| — root-level, git-tracked | 4 |
| — root-level, untracked | 3 |
| — `src/` whole packages → `legacy/` (Batch L1) | 26 |
| — `src/` partial-package files → `legacy/` (Batch L2, optional) | 11 |
| — `scripts/` loose top-level → subdirs | 63 |
| Files/paths proposed to **delete** | **57** |
| **gitignore-only** changes | **9 rule edits** (0 files moved) |
| **High-risk** moves | **9** |
| Moves blocked outright (do not attempt) | **22 `src/` files** (§3.3) |
| Documented reproduction commands broken | **0 in `README.md`**; 8 elsewhere (§6) |

**Biggest breakage trap** — see §7.

---

## 2. Repo-root inventory and classification

Evidence commands:
`ls -la`; `git ls-files | grep -v /`; `git check-ignore -v <f>`;
`git ls-files --error-unmatch <f>`; `grep -rnF "<f>" README.md docs/ config/ scripts/ src/ tests/ Dockerfile .dockerignore`.

### 2.1 Every root-level entry, classified

| Path | tracked | ignored | referenced by | verdict |
|---|---|---|---|---|
| `.gitignore`, `LICENSE`, `README.md`, `requirements.txt`, `environment.yml`, `Dockerfile`, `.dockerignore` | yes | no | build/CI | **keep at root** |
| `agents_instruction.md` | **yes** | no | `config/hyperparameters.yaml:120`, `config/contrast_vector_framework.yaml:3`, `config/revisions/contrast-vectors.yaml` `docs:[0]`, `docs/plan_changelog.md:79`, `docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:235,393` | **keep at root** — 5 refs, one of them registry-asserted |
| `install_limma.R`, `install_r_packages.R` | yes | no | `docs/v11_novelty_extensions_implementation_plan_2026-06-06.md:57` | **keep at root** (env bootstrap) |
| `email_to_casaletto_draft.md` | **yes** | no | `docs/STATE_OF_PLAY_2026-07-29.md:134` | **move → `docs/`** |
| `manuscript_v11_prose_audit.md` | **yes** | no | none found | **move → `docs/`** |
| `manuscript_v12_prose_audit.md` | **yes** | no | none found | **move → `docs/`** |
| `latex/repository_statistical_audit_2026_04_24.tex` | **yes** | no | none found (`grep -rn 'latex/' README.md docs/ config/` → 0 non-`latex_paper` hits) | **move → `docs/`**, then `rmdir latex/` |
| `latex_paper/` (2 tracked of 296 on disk) | partial | rule `latex_paper/*` + `!manuscript_v11.*` | `README.md:183,185` point at **manuscript_v13** | **keep**; fix gitignore (§4) |
| `figures/` (9 tracked) | yes | no | v13 + clinical-axes figures | **keep at root** |
| `config/`, `data/`, `docs/`, `scripts/`, `src/`, `tests/` | yes | no | everything | **keep at root** |
| `CLOUD_DEPLOYMENT.md` | no | `/CLOUD_DEPLOYMENT.md` | documents `src/run_all_phases.py` (3 commands) | **gitignore only** (update if L1 runs) |
| `FILE_SUMMARY.md` | no | `/FILE_SUMMARY.md` | **`config/revisions/network-rewiring.yaml` `docs:[1]` → asserted by `tests/test_rrrm2_registry.py:50`** | **keep at root — do NOT move or delete** |
| `METHODOLOGY.md` | no | `/METHODOLOGY.md` | **`config/revisions/network-rewiring.yaml` `docs:[0]` → registry-asserted** | **keep at root — do NOT move or delete** |
| `RESULTS.md` | no | `/RESULTS.md` | `scripts/audit_fix_status.py:57` (`REPO / "RESULTS.md"`) | **keep at root** |
| `COMMENT_ARCHIVE.md` (122 KB) | no | `COMMENT_ARCHIVE.md` | historical command archive; ~20 stale commands | **gitignore only** |
| `walkthrough.md` | no | `/walkthrough.md` | none | **DELETE** — content is `# OLS Hardware Driver Integration — Walkthrough`, about OpenLiveStacker C++ drivers and an "Achernar Rust daemon". Belongs to a **different project**. |
| `results.txt` (18 KB) | no | `/results.txt` | none | **DELETE** — generated `PIPELINE: SYNTHESIZED RESULTS … run_20260408 vs run_20260209`, superseded Workflow-A output |
| `chrome_paper_doi_links.txt` (9.4 KB) | no | `/chrome_paper_doi_links.txt` | none | **move → `tmp/`** (already-ignored dir). Content is a 2026-05-14 DOI collection — research value, keep it, just not at root |
| `.chrome_paper_doi_links.txt.swp` | no | `*.swp` | none | **DELETE** (vim swapfile) |
| `permutedStats-actualModules.RData` (564 KB) | no | `*.RData` **and** `/permutedStats-...` | none | **move → `data/external/legacy_wgcna/`**. Not delete: it is a WGCNA permutation intermediate that was never in git history, so deletion is irreversible. |
| `node2vec_asgsr_osd771_kidney_transcriptomics.pdf` (231 KB) | no | `/node2vec_...pdf` | **`docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:15` and `:383`** — "the founding artifact", "the founding proposal" | **move → `docs/` and track**; update the 2 citations. Never delete: it is cited twice as the project's origin document. |
| `.DS_Store` (41 copies repo-wide, `find . -name .DS_Store -not -path './venv/*' -not -path './.git/*' \| wc -l` → 41) | no | `.DS_Store` | none | **DELETE all** |
| `.idea/` (7 files) | no | `.idea/` | none | **DELETE** (IDE state) |
| `.pytest_cache/` | no | `.pytest_cache/` | pytest | **gitignore only** (regenerated) |
| `.claude/` | no | `.claude/` | agent state | **gitignore only** |
| `deconv_music_out/` | no | `/deconv_music_out/` | `docs/v11_novelty_extensions_implementation_plan_2026-06-06.md:57` describes it as *"the empty `deconv_music_out/` target"* | **DELETE** — 0 files on disk; the one reference is prose noting it is empty |
| `tms_kidney_female_counts/` | no | ignored | none | **DELETE** — 0 files |
| `tms_kidney_ytrfemale_counts/` | no | ignored | none | **DELETE** — 0 files |
| `plots/` | no | `/plots/` | output root in `src/visualization/network_diagnostics.py:337`, `scripts/plot_output.py` | **DELETE** — contains only `.DS_Store`; recreated on demand |
| `notebooks/` (3 files, 0 tracked) | no | 2 explicit rules | none | **gitignore only** — replace 2 rules with `/notebooks/` |
| `tmp/` (43 MB), `transcriptomic_texts/` (182 MB), `venv/` (1.2 GB) | no | ignored | — | **gitignore only** |

### 2.2 Non-root junk found while inventorying

| Path | evidence | verdict |
|---|---|---|
| `src/networks/.edge_regression.py.swp`, `.embeddings.py.swp`, `.lioness.py.swp` | `ls -la src/networks/.*.swp` → 3 × 16 KB, May 17 | **DELETE** |
| `src/plots/` (2 dirs, 7 files) | Created by the bug at `src/visualization/network_diagnostics.py:34` — `REPO_ROOT = Path(__file__).resolve().parents[1]` resolves to `src/`, not the repo root, and `:337` writes `REPO_ROOT / "plots" / …`. `.gitignore` line `/src/plots/` papers over it. Untracked. | **DELETE** + fix the index (§5) |
| `src/data/results/` (6 dirs, 2 files) | Same failure mode; `.gitignore` has `/src/data/results/`. Untracked. | **DELETE** |
| `scripts/data/processed/networks/phase2/regression/` | `find scripts/data -type f \| wc -l` → **0**. Empty dir tree from a run with the wrong CWD. Untracked. | **DELETE** |
| `scripts/run_permutations.py` | `wc -c` → **0 bytes**; `git ls-files -s` → `e69de29…` (the empty blob). Tracked, referenced nowhere. | **`git rm`** |

---

## 3. Legacy archival — Workflow A

### 3.1 How the boundary was determined (not guessed)

An AST import graph was built over all 244 `.py` files in `src/`, `scripts/`,
`tests/` (`ast.Import` / `ast.ImportFrom`, relative imports resolved by level),
then a **transitive reachability closure** was computed from these seeds:

- everything in `tests/`
- `scripts/{osd462,v13,v11,clinical_axes,subtype_reference,stage0,celltype,regulator_activity,grey60,integration}/`
- `src/{v11,v13,clinical_axes,subtype_reference,multiomics,grey60}/`

A `src/` file reachable from that closure **cannot move** without editing its
consumers.

### 3.2 The chain that decides it — `src/networks/` is load-bearing for v13

```
src/v13/continuous_phospho_inference.py:34
    from src.multiomics.osd462_anchor import PLEX1, PLEX2, TmtTable, parse_tmt_sheet
src/multiomics/osd462_anchor.py:534   (inside aligned_pathway_cosine, function-local)
    from src.networks.contrast_vectors import cosine
```

The locked v13 engine depends on a Workflow-A module, via an import that is
**not at module top level**. `tests/test_osd462_anchor.py:110` does call
`aligned_pathway_cosine`, so the suite would catch it — but a top-level-imports-only
audit would not. See §7.

### 3.3 Per-file verdict, Workflow-A candidate packages

`BLOCKED` = reachable from tests or current workflows. `MOVABLE` = not
reachable, transitively.

| File | Verdict | Importers (class) |
|---|---|---|
| `src/networks/contrast_vectors.py` | **BLOCKED** | `tests/test_contrast_decomposition.py` [TEST], `scripts/osd462/04_rna_recurrence.py` [CUR], `src/multiomics/osd462_anchor.py` [CUR, lazy], + 4 legacy |
| `src/networks/cross_osdr_projection.py` | **BLOCKED** | `tests/test_cross_osdr_alignment.py`, `src/v11/aldosterone_axis.py` [CUR], `src/v11/recurrence_meta.py` [CUR] |
| `src/networks/edge_regression.py` | **BLOCKED** | `tests/test_lioness_methodology.py` |
| `src/networks/external_aging_axis.py` | **BLOCKED** | `tests/test_external_aging_axis.py` |
| `src/networks/lar_reversal.py` | **BLOCKED** | `tests/test_lar_reversal.py` |
| `src/networks/lioness.py` | **BLOCKED** | `tests/test_lioness_methodology.py` |
| `src/networks/manuscript_decision.py` | **BLOCKED** | `tests/test_decision_tree.py` |
| `src/networks/mechanism_axis.py` | **BLOCKED** | `tests/test_mechanism_axis.py` |
| `src/networks/procrustes.py` | **BLOCKED** | `tests/test_procrustes_anchors.py` |
| `src/networks/shared_topology.py` | **BLOCKED** | `tests/test_config_consistency.py`, `tests/test_lioness_methodology.py` |
| `src/networks/stability_test.py` | **BLOCKED** | `tests/test_stability_test.py` |
| `src/networks/tubulointerstitial_state.py` | **BLOCKED** | `tests/test_tubulointerstitial_state.py` |
| `src/networks/alternative_methods.py` | **BLOCKED (transitive)** | `src/validation/enhanced_cv.py` → `tests/test_enhanced_cv_fold_safety.py` |
| `src/statistics/full_pipeline_permutation.py` | **BLOCKED** | `tests/test_full_pipeline_permutation.py` |
| `src/statistics/full_regression.py` | **BLOCKED** | `tests/test_full_regression_aggregation.py` |
| `src/statistics/permutation_bootstrap.py` | **BLOCKED** | `tests/test_config_consistency.py`, `tests/test_permutation_modes.py` |
| `src/statistics/silent_shifters.py` | **BLOCKED** | `tests/test_silent_shifters.py` |
| `src/validation/cross_validation.py` | **BLOCKED** | `tests/test_config_consistency.py`, `tests/test_lioness_methodology.py` |
| `src/validation/enhanced_cv.py` | **BLOCKED** | `tests/test_enhanced_cv_fold_safety.py` |
| `src/validation/external_replication.py` | **BLOCKED** | `tests/test_external_protocol.py` |
| `src/validation/osd_external_validation.py` | **BLOCKED** | `tests/test_osd_external_validation.py` |
| `src/validation/sample_features.py` | **BLOCKED (transitive)** | `src/validation/cross_validation.py`, `enhanced_cv.py` |
| `src/preprocessing/` — all 11 | MOVABLE | none current |
| `src/enrichment/` — all 4 | MOVABLE | internal only |
| `src/visualization/` — all 5 | MOVABLE | internal only (`publication_figures.py:936` lazily imports `src.v11.publication_figures`, which **stays** — absolute import survives the move) |
| `src/markers/` — all 3 | MOVABLE | `scripts/discover_dct_markers.py` (also moving) |
| `src/data/` — 2 `.py` | MOVABLE | `scripts/build_id_map.py` (also moving) |
| `src/networks/embeddings.py`, `wgcna_followup.py`, `wgcna_analysis.R`, `wgcna_gc_reference_preservation.R` | MOVABLE | `scripts/phase3_node2vec_embedding.py` (also moving) |
| `src/statistics/bootstrap_decomposition.py`, `direct_coexpression_test.py`, `interaction_metrics.py`, `differential_expression.R` | MOVABLE | none |
| `src/validation/continuous_target.py`, `multi_study_pool.py`, `wgcna_external_validation.py` | MOVABLE | none |
| `src/run_all_phases.py` | MOVABLE **but** referenced by `Dockerfile:73`, `config/revisions/network-rewiring.yaml:42,51` | move with edits |

**Verdict: `src/networks/`, `src/statistics/`, `src/validation/` cannot be
archived as packages.** 22 of their files are test-reachable. Archiving them
would require editing 18 test files (`test_lioness_methodology`,
`test_procrustes_anchors`, `test_silent_shifters`, `test_config_consistency`,
`test_enhanced_cv_fold_safety`, `test_permutation_modes`,
`test_full_regression_aggregation`, `test_full_pipeline_permutation`,
`test_external_protocol`, `test_osd_external_validation`,
`test_contrast_decomposition`, `test_cross_osdr_alignment`, `test_decision_tree`,
`test_external_aging_axis`, `test_lar_reversal`, `test_mechanism_axis`,
`test_stability_test`, `test_tubulointerstitial_state`) **plus** the v13 chain
in §3.2. Not recommended.

### 3.4 Proposed legacy layout

```
legacy/
  workflow_a_network_rewiring/
    __init__.py                 (NEW — required for `python -m legacy…`)
    run_all_phases.py
    preprocessing/              (whole package, 11 files)
    enrichment/                 (whole package, 4 files)
    visualization/              (whole package, 5 files)
    markers/                    (whole package, 3 files)
    idmap/                      (from src/data/, 2 files — renamed to avoid a
                                 second directory called `data`)
    networks/   __init__.py NEW (4 WGCNA/embedding files)      [Batch L2]
    statistics/ __init__.py NEW (4 files)                      [Batch L2]
    validation/ __init__.py NEW (3 files)                      [Batch L2]
    scripts/
      contrast_vectors/
  scratch/
```

**Recommendation:** execute **Batch L1** (whole packages + `run_all_phases.py`)
and **skip Batch L2**. Splitting `src/networks`, `src/statistics`,
`src/validation` across two trees for 11 files makes the layout *harder* to
read than leaving them, and each surviving package still needs its `__init__.py`
in place. If L2 is skipped, add a one-line header comment to those 11 files
instead. L2 is specified below for completeness, marked optional.

Absolute `from src.common import …` and `from src.networks.contrast_vectors import …`
inside moved files keep working unchanged — repo root is on `sys.path` for every
entry point. Only imports of modules that *themselves moved* need rewriting.

---

## 4. THE MOVE TABLE

Risk key: **none** = no consumer of any kind; **low** = consumers exist and the
test suite or a registry test catches a mistake; **high** = a consumer exists
that nothing automatically checks, or a frozen artifact must be edited.

### 4.1 Root-level

| current path | proposed path | rationale | risk | must be updated alongside |
|---|---|---|---|---|
| `email_to_casaletto_draft.md` | `docs/email_to_casaletto_draft_2026-06-19.md` | correspondence draft, not a root artifact | low | `docs/STATE_OF_PLAY_2026-07-29.md:134` |
| `manuscript_v11_prose_audit.md` | `docs/manuscript_v11_prose_audit.md` | audit prose belongs with the other audits | none | — |
| `manuscript_v12_prose_audit.md` | `docs/manuscript_v12_prose_audit.md` | same | none | — |
| `latex/repository_statistical_audit_2026_04_24.tex` | `docs/repository_statistical_audit_2026_04_24.tex` | sole occupant of a root dir; unreferenced | none | `rmdir latex/` after |
| `node2vec_asgsr_osd771_kidney_transcriptomics.pdf` | `docs/founding_proposal_2025-12-09_node2vec_asgsr_osd771.pdf` | cited twice as the founding artifact but currently untracked and gitignored | **high** | `docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:15` and `:383`; drop `.gitignore:/node2vec_…pdf`; `git add` the new path |
| `permutedStats-actualModules.RData` | `data/external/legacy_wgcna/permutedStats-actualModules.RData` | 564 KB binary at root; never tracked, so not recoverable if deleted | none | — (target dir is inside ignored `data/external/`) |
| `chrome_paper_doi_links.txt` | `tmp/chrome_paper_doi_links.txt` | literature collection, deliberately gitignored; keep content, unclutter root | none | drop `.gitignore:/chrome_paper_doi_links.txt` |

### 4.2 `src/` → `legacy/` — Batch L1 (whole packages)

| current path | proposed path | rationale | risk | must be updated alongside |
|---|---|---|---|---|
| `src/preprocessing/` (`__init__.py`, `common_gene_universe.py`, `data_alignment.py`, `deconvolution_sanity.py`, `export_counts.py`, `qc_variance.py`, `deconvolution.R`, `deconvolution_sensitivity.R`, `export_phase1.R`, `multi_study_harmonization.R`, `residualization.R`) | `legacy/workflow_a_network_rewiring/preprocessing/` | Workflow A phases 0–1; zero current/test consumers | low | `run_all_phases.py` module string `src.preprocessing.qc_variance` + 4 `run_rscript` paths; depth fixes (§5) |
| `src/enrichment/` (`__init__.py`, `biological_grounding.py`, `gene_set_loader.py`, `regression_enrichment.py`) | `legacy/workflow_a_network_rewiring/enrichment/` | Workflow A phase 7 | low | `run_all_phases.py:` `src.enrichment.biological_grounding`; internal `src.enrichment.gene_set_loader` imports in `biological_grounding.py` and `regression_enrichment.py` (+ the docstring at `gene_set_loader.py:13`) |
| `src/visualization/` (`__init__.py`, `network_diagnostics.py`, `publication_figures.py`, `publication_plots.py`, `wgcna_publication_figures.py`) | `legacy/workflow_a_network_rewiring/visualization/` | Workflow A plotting; the v13 figures live in `src/v13/reporting.py` | low | 3 `run_all_phases.py` module strings; depth fixes (§5) |
| `src/markers/` (`__init__.py`, `discover_dct.py`, `discover_markers.py`) | `legacy/workflow_a_network_rewiring/markers/` | Workflow A marker discovery | low | 2 `run_all_phases.py` module strings; `scripts/discover_dct_markers.py` wrapper (also moving) |
| `src/data/__init__.py`, `src/data/build_id_map.py` | `legacy/workflow_a_network_rewiring/idmap/` | Ensembl→symbol map builder, Workflow A only | low | `run_all_phases.py` × 2 (`src.data.build_id_map` at `:206,:211`); `scripts/build_id_map.py` wrapper (also moving); `src/networks/shared_topology.py:526,603` and `src/enrichment/biological_grounding.py:257` mention it in **error-message text only** — update for accuracy, not correctness |
| `src/run_all_phases.py` | `legacy/workflow_a_network_rewiring/run_all_phases.py` | the Workflow A orchestrator, 1,413 lines, status `retired` | **high** | `Dockerfile:73` `ENTRYPOINT ["python3", "src/run_all_phases.py"]`; `config/revisions/network-rewiring.yaml:42` and `:51` (**registry-asserted**); `CLOUD_DEPLOYMENT.md:161,178,255`; `FILE_SUMMARY.md:128,131,134`; its own 13 module strings + 6 rscript paths (§4.2a); `parent.parent` → `parents[2]` |

**4.2a — the 19 strings inside `run_all_phases.py` that must be rewritten**
(the other 13 module strings point at modules that **stay** and must be left alone):

*rewrite* — `src.data.build_id_map` (×2), `src.markers.discover_dct`,
`src.markers.discover_markers`, `src.preprocessing.qc_variance`,
`src.networks.embeddings`†, `src.networks.wgcna_followup`†,
`src.statistics.direct_coexpression_test`†, `src.statistics.interaction_metrics`†,
`src.enrichment.biological_grounding`, `src.validation.wgcna_external_validation`†,
`src.visualization.publication_figures`, `src.visualization.publication_plots`,
`src.visualization.wgcna_publication_figures`; and the six
`run_rscript(...)` paths `src/preprocessing/{deconvolution,deconvolution_sensitivity,residualization,export_phase1}.R`,
`src/networks/wgcna_analysis.R`†, `src/statistics/differential_expression.R`†.
(† = only if Batch L2 also runs.)

*leave alone* — `src.networks.{shared_topology,lioness,edge_regression,procrustes}`,
`src.statistics.{silent_shifters,permutation_bootstrap,full_pipeline_permutation,full_regression}`,
`src.validation.{cross_validation,enhanced_cv,osd_external_validation}`,
`src.v11.{spatial_dct_transport_check,spatial_reference_projection}`.

### 4.3 `src/` → `legacy/` — Batch L2 (partial packages, OPTIONAL, not recommended)

| current path | proposed path | risk | must be updated alongside |
|---|---|---|---|
| `src/networks/embeddings.py` | `legacy/workflow_a_network_rewiring/networks/embeddings.py` | low | `scripts/phase3_node2vec_embedding.py` (moving); new `networks/__init__.py`; `src/networks/__init__.py:24` `try: from . import embeddings` — swallowed by its own `except ImportError: pass`, verified by reading the file |
| `src/networks/wgcna_followup.py` | `…/networks/wgcna_followup.py` | low | `run_all_phases.py`; `.parent.parent.parent` → `parents[3]` |
| `src/networks/wgcna_analysis.R`, `wgcna_gc_reference_preservation.R` | `…/networks/` | low | `run_all_phases.py:601` rscript path |
| `src/statistics/bootstrap_decomposition.py` | `…/statistics/` | none | — (its `src.networks.contrast_vectors` import stays valid) |
| `src/statistics/direct_coexpression_test.py` | `…/statistics/` | low | `run_all_phases.py:646` |
| `src/statistics/interaction_metrics.py` | `…/statistics/` | low | `run_all_phases.py:471`; `src/statistics/__init__.py` try/except |
| `src/statistics/differential_expression.R` | `…/statistics/` | low | `run_all_phases.py:569` |
| `src/validation/continuous_target.py`, `multi_study_pool.py`, `wgcna_external_validation.py` | `…/validation/` | low | `run_all_phases.py:659`; `wgcna_external_validation.py`'s `src.validation.osd_external_validation` import **stays valid** (that module does not move) |

### 4.4 `scripts/` regrouping — all 65 loose tracked top-level entries

`git ls-files scripts | awk -F/ 'NF==2'` → **65**. Every one is placed below.

#### 4.4a Workflow A phase pipeline → `legacy/workflow_a_network_rewiring/scripts/`

| current path | rationale | risk | must be updated alongside |
|---|---|---|---|
| `scripts/phase2_edge_regression.py` | 8-line wrapper for `src.networks.edge_regression` | low | — (target module stays) |
| `scripts/phase3_node2vec_embedding.py` | wrapper for `src.networks.embeddings` | low | `run_phase3_pipeline.py:106` hardcoded `"scripts/phase3_node2vec_embedding.py"`; import path if L2 runs |
| `scripts/phase3_node_rewiring_from_deltas.py` | phase 3.1 | low | `run_phase3_pipeline.py:94`; `parents[1]`→`[3]` |
| `scripts/phase3_procrustes_rewiring.py` | phase 3.3 | low | `run_phase3_pipeline.py:125`; `parents[1]`→`[3]` |
| `scripts/phase4_anchor_qc_report.py` | phase 4 prereg doc | low | `run_phase4_7_pipeline.py:69`; `parents[1]`→`[3]` |
| `scripts/phase5_build_silent_shifters_strict.py` | wrapper | low | `run_phase4_7_pipeline.py:116` |
| `scripts/phase5_derive_interaction_persistence.py` | phase 5 derived | low | `run_phase4_7_pipeline.py:66`; `parents[1]`→`[3]` |
| `scripts/phase6_edge_regression_full.py` | wrapper | none | — |
| `scripts/phase6_perm_bootstrap_node_rewiring.py` | wrapper | low | `run_phase4_7_pipeline.py:73` |
| `scripts/phase7_grounding_fast.py` | phase 7 | low | `run_phase4_7_pipeline.py:93`; `parents[1]`→`[3]` |
| `scripts/build_phase2_skeleton.py` | wrapper for `src.networks.shared_topology` | low | `run_phase2_pipeline.py:65` (comment only) |
| `scripts/compute_phase2_lioness_on_skeleton.py` | wrapper for `src.networks.lioness` | none | — |
| `scripts/run_phase1_networks.py` | **already dead**: `runpy` → `ImportError: cannot import name 'lioness_correlation_edges' from 'src.networks.lioness'`; also imports the deleted `src.statistics.rewiring_metrics` | none | `.parent.parent`→`parents[3]` |
| `scripts/run_phase2_pipeline.py` | phase 2 orchestrator | low | `.parent.parent` at `:35`→`parents[3]` |
| `scripts/run_phase3_pipeline.py` | phase 3 orchestrator | low | 3 hardcoded `"scripts/…"` strings; `parents[1]`→`[3]` |
| `scripts/run_phase4_7_pipeline.py` | phases 4–7 orchestrator | low | `parents[1]`→`[3]` (`SCRIPTS = REPO_ROOT / "scripts"` at `:26` must become the legacy scripts dir) |
| `scripts/run_full_pipeline.py` | `FILE_SUMMARY.md`: "Legacy placeholder (not currently used)"; `docs/provenance_audit_2026-06-06.md:110`: "every phase is a `logger.info("[Placeholder] …")`" | low | `.parent.parent / 'src'`→`parents[3] / 'src'`; `docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:197,425`; `docs/provenance_audit_2026-06-06.md:110` |
| `scripts/pick_anchors.py` | **already dead**: `ModuleNotFoundError: No module named 'src.statistics.rewiring_metrics'` (deleted in `beab419`) | none | `.parent.parent`→`parents[3]` |
| `scripts/make_consensus_anchors.py` | Procrustes anchor selection | none | — |
| `scripts/plot_output.py` | all-phase publication plots | low | `parents[1]`→`[3]`; `src/visualization/publication_plots.py:2,766` mention it in docstring/help text |
| `scripts/plot_skeleton_diagnostics.py` | phase 2 diagnostics | low | `parents[1]`→`[3]`; `src/visualization/network_diagnostics.py:1,13` docstring |
| `scripts/batch_embeddings.sh` | node2vec batch driver | none | — |
| `scripts/align_metadata_to_counts.py` | GLDS-674 metadata alignment (phase 0) | low | `.parent.parent`→`parents[3]` |
| `scripts/export_raw_counts.py` | raw-count export for MuSiC | low | `scripts/run_deconvolution.R:144` and `src/preprocessing/deconvolution.R:149` both `stop("… Run scripts/export_raw_counts.py")` |
| `scripts/build_id_map.py` | 6-line wrapper for `src.data.build_id_map` | low | `os.path.dirname(__file__), ".."` → `"..", ".."`; import path (target moves in L1) |
| `scripts/discover_dct_markers.py` | 6-line wrapper for `src.markers.discover_dct` | low | same two edits |
| `scripts/audit_fix_status.py` | self-audit of the 14 remediation fixes | low | `parents[1]`→`[3]`; keeps reading `REPO / "RESULTS.md"` — verify after the index change |
| `scripts/DESeq2.R` | VST source | none | — |
| `scripts/run_gene_level_DE.R` | limma DE | none | `src/statistics/differential_expression.R:1` header comment |
| `scripts/export_phase1_to_python.R` | R→Python export | none | `scripts/phase1_to_python.R:1`, `src/preprocessing/export_phase1.R:1` header comments |
| `scripts/phase1_to_python.R` | duplicate of the above | none | — |
| `scripts/run_deconvolution.R` | MuSiC deconvolution driver | none | — |

#### 4.4b Contrast-vector era → `legacy/workflow_a_network_rewiring/scripts/contrast_vectors/`

| current path | rationale | risk | must be updated alongside |
|---|---|---|---|
| `scripts/run_contrast_vector_framework.py` | the May-2026 contrast-vector orchestrator; `config/revisions/contrast-vectors.yaml` marks it workflow **A / retired** | **high** | `config/revisions/contrast-vectors.yaml:37` (**registry-asserted**); `src/run_all_phases.py:933`; `scripts/audit_stability_gates.py:103` (message text); `config/contrast_vector_framework.yaml:6` comment; `docs/NEXT_PAPER_DIRECTION_AFTER_GREY60_2026-07-29.md:196`; `parents[1]`→`[4]`; its own 3 sub-script paths at `:1000,:1037,:1061` |
| `scripts/audit_stability_gates.py` | stability-gate auditor for the above | **high** | `config/revisions/contrast-vectors.yaml:46` (**registry-asserted**); `parents[1]`→`[4]` |
| `scripts/run_lar_reversal_analysis.py` | LAR reversal; `docs/MANUSCRIPT_V1_V12_CROSS_VERSION_AUDIT:94` lists "LAR reversal" under **Retire** | **high** | `run_contrast_vector_framework.py:1061`; `docs/annotated_issues_v5.md:32,37,42,47`; `docs/statistical_assessment_v5.md:4`; `parents[1]`→`[4]` |
| `scripts/run_mechanism_axis_prioritization.py` | mechanism-axis prioritization | low | `run_contrast_vector_framework.py:1000`; `parents[1]`→`[4]` (note `:448` `outdir.parents[1]` is an **output**-path index — leave it alone) |
| `scripts/run_tubulointerstitial_state_analysis.py` | tubulointerstitial state space | low | `run_contrast_vector_framework.py:1037`; `parents[1]`→`[4]` |
| `scripts/plot_lar_reversal_dashboard.py` | LAR dashboard | low | `docs/statistical_assessment_v5.md:4`; `parents[1]`→`[4]` |
| `scripts/plot_manuscript_v4_figures.py` | manuscript v4 figures (superseded by v13) | none | `parents[1]`→`[4]` |

#### 4.4c One-off scratch → `legacy/scratch/`

All twenty have **zero** references in `docs/`, `config/`, `src/`, `scripts/`,
`tests/` (verified per file). Risk **none** unless noted.

`analyze_dct_variance.py`, `audit_reference_labels.py`, `check_dct_universe.py`,
`check_reference_counts.py`, `compare_two_runs.py`,
`debug_deconvolution_phase0.py`, `explore_chen_atlas.R`,
`extract_cell_types_from_csv.py`, `extract_cell_types_to_file.py`,
`extract_ecm_edges.py`, `inspect_10x_reference.py`, `inspect_sce_structure.R`,
`map.R`, `map_top_hits.py`, `pathway_presence.py`, `run_pseudobulk_test.R`,
`sanity_marker_panel.R`, `targeted_fdr.py`,
`module_convergence_drilldown.py` (`parents[1]`→`[2]`),
`module_convergence_fisher.py` (`parents[1]`→`[2]`).

#### 4.4d New experiment subdirs

| current path | proposed path | rationale | risk | must be updated alongside |
|---|---|---|---|---|
| `scripts/spaceflight_kidney_classifier.py` | `scripts/sf_classifier/spaceflight_kidney_classifier.py` | 2026-07-03 cross-cohort classifier; `docs/ANALYSIS_TRIAGE_AND_NETWORK_REFRAME_2026-07-29.md:148` calls its machinery live and reusable ("the one positive result in this repository has exactly this shape"). **Not scratch.** | low | `docs/_commentary_scripts.md:848`; script uses `ROOT` computed by its own idiom — verify `data/results/run_20260701_sf_classifier` still resolves |
| `scripts/sf_classifier_component_decomp.py` | `scripts/sf_classifier/` | same class | low | `docs/_commentary_scripts.md:799` |
| `scripts/sf_classifier_negctrl_figure.py` | `scripts/sf_classifier/` | same class | low | `docs/_commentary_scripts.md:823` |
| `scripts/osd253_strain_module_projection.py` | `scripts/grey60/osd253_strain_module_projection.py` | its **only** consumer is `scripts/grey60/04_run_external.py:15` (`from scripts.osd253_strain_module_projection import …`) — this is Grey60 machinery, not a shared utility | **high** | `scripts/grey60/04_run_external.py:15` import path → `scripts.grey60.osd253_strain_module_projection`. **No test covers `scripts/grey60/04_run_external.py`** (`tests/test_grey60_adversarial.py` imports `src.grey60.adversarial` only), so this breakage is silent — verify by hand |

#### 4.4e Stays where it is

| path | why |
|---|---|
| `scripts/build_hybrid_reference.R` | **Recommend NOT moving.** It is pinned inside a frozen spec: `config/dct_subtype_reference_freeze_v1.yaml:91` `cluster_map_source: scripts/build_hybrid_reference.R`, and cited at `docs/DCT_SUBTYPE_REFERENCE_REBUILD_2026-07-28.md:33` ("it is frozen in the config"). Moving it means editing a frozen provenance record for cosmetic gain. *If* it is moved, the natural home is `scripts/subtype_reference/` and both lines must change — risk **high**. Note: editing that YAML does **not** invalidate the published SHA-256 `ed64f0e2…`; `shasum -a 256 config/dct_subtype_reference_freeze_v1.yaml` → `a79a7be8…`, so `ed64f0e2…` is the generated membership file under `data/results/`, not the config. |

#### 4.4f Untracked entries under `scripts/`

| path | ignored by | verdict |
|---|---|---|
| `scripts/sample.py` | `.gitignore:108 /scripts/sample.py` | **DELETE** (scratch, never tracked) |
| `scripts/phase7_2` | `.gitignore:109 /scripts/phase7_2` | **DELETE** (scratch, never tracked) |
| `scripts/.DS_Store` | `.DS_Store` | **DELETE** |

---

## 5. `parents[N]` depth-sensitivity table

This is the silent-failure surface: `Path(__file__).resolve().parents[N]` does
not raise when `N` is too small, it just returns a directory that is *not* the
repo root, and downstream paths then resolve under `src/` or `scripts/`.

Depth rule: a file `a/b/c.py` needs `parents[2]`; `a/b.py` needs `parents[1]`.

Found with `grep -rn --include='*.py' "parents\[" src scripts tests` **and**
`grep -rn -E "\.parent\.parent|dirname\(__file__\)"` — the second grep matters,
because 10 files use the `.parent.parent` spelling and are invisible to the first.

### 5.1 Files whose index changes under this plan

| file (current) | proposed path | current expr | required expr |
|---|---|---|---|
| `scripts/audit_fix_status.py:8` | `legacy/workflow_a_network_rewiring/scripts/` | `parents[1]` | `parents[3]` |
| `scripts/phase3_node_rewiring_from_deltas.py:21` | same | `parents[1]` | `parents[3]` |
| `scripts/phase3_procrustes_rewiring.py:21` | same | `parents[1]` | `parents[3]` |
| `scripts/phase4_anchor_qc_report.py:27` | same | `parents[1]` | `parents[3]` |
| `scripts/phase5_derive_interaction_persistence.py:22` | same | `parents[1]` | `parents[3]` |
| `scripts/phase7_grounding_fast.py:28` | same | `parents[1]` | `parents[3]` |
| `scripts/plot_output.py:31` | same | `parents[1]` | `parents[3]` |
| `scripts/plot_skeleton_diagnostics.py:34` | same | `parents[1]` | `parents[3]` |
| `scripts/run_phase3_pipeline.py:22` | same | `parents[1]` | `parents[3]` |
| `scripts/run_phase4_7_pipeline.py:25` | same | `parents[1]` | `parents[3]` |
| `scripts/align_metadata_to_counts.py:14` | same | `.parent.parent` (≡`parents[1]`) | `parents[3]` |
| `scripts/pick_anchors.py:19` | same | `.parent.parent` | `parents[3]` |
| `scripts/run_phase1_networks.py:19` | same | `.parent.parent` | `parents[3]` |
| `scripts/run_phase2_pipeline.py:35` | same | `.parent.parent` | `parents[3]` |
| `scripts/run_full_pipeline.py:12` | same | `.parent.parent / 'src'` | `parents[3] / 'src'` |
| `scripts/build_id_map.py:4` | same | `dirname(__file__), ".."` | `dirname(__file__), "..", ".."` |
| `scripts/discover_dct_markers.py:4` | same | `dirname(__file__), ".."` | `dirname(__file__), "..", ".."` |
| `scripts/audit_stability_gates.py:13` | `…/scripts/contrast_vectors/` | `parents[1]` | `parents[4]` |
| `scripts/run_contrast_vector_framework.py:19` | same | `parents[1]` | `parents[4]` |
| `scripts/run_lar_reversal_analysis.py:16` | same | `parents[1]` | `parents[4]` |
| `scripts/run_mechanism_axis_prioritization.py:18` | same | `parents[1]` | `parents[4]` |
| `scripts/run_tubulointerstitial_state_analysis.py:15` | same | `parents[1]` | `parents[4]` |
| `scripts/plot_lar_reversal_dashboard.py:17` | same | `parents[1]` | `parents[4]` |
| `scripts/plot_manuscript_v4_figures.py:16` | same | `parents[1]` | `parents[4]` |
| `scripts/module_convergence_drilldown.py:6` | `legacy/scratch/` | `parents[1]` | `parents[2]` |
| `scripts/module_convergence_fisher.py:11` | `legacy/scratch/` | `parents[1]` | `parents[2]` |
| `src/run_all_phases.py:14` | `legacy/workflow_a_network_rewiring/` | `.parent.parent` (≡`parents[1]`) | `parents[2]` |
| `src/preprocessing/deconvolution_sanity.py:11` | `…/preprocessing/` | `parents[2]` | `parents[3]` |
| `src/enrichment/regression_enrichment.py:22` | `…/enrichment/` | `parents[2]` | `parents[3]` |
| `src/visualization/publication_plots.py:31` | `…/visualization/` | `parents[2]` | `parents[3]` |
| `src/visualization/publication_figures.py:963` | same | `parents[2]` | `parents[3]` |
| `src/networks/wgcna_followup.py:23` *(L2 only)* | `…/networks/` | `.parent.parent.parent` (≡`parents[2]`) | `parents[3]` |

### 5.2 Files that are ALREADY WRONG — fix while moving, and say so

| file | current | resolves to | evidence it fires | required |
|---|---|---|---|---|
| `src/visualization/network_diagnostics.py:34` | `parents[1]` | `src/` | `:337` writes `REPO_ROOT / "plots" / …`; **`src/plots/20260407_182630_skeleton_diagnostics/` exists on disk with 7 files**, and `.gitignore` carries a `/src/plots/` rule to hide it | `parents[3]` at the new depth |
| `src/preprocessing/data_alignment.py:14` | `.parent.parent` | `src/` | `:16-21` build `data/raw/...`, `data/processed/...` under `src/`; `.gitignore` carries `/src/data/results/`, and `src/data/results/` exists with `phase5_silent_shifters_strict/` and `phase7_grounding/` outputs | `parents[3]` at the new depth |

**Fixing these changes runtime behaviour** (outputs move from `src/data`,
`src/plots` to `data/`, `plots/`). That is the correct behaviour, but it is a
behaviour change, not a no-op. Do it deliberately, in its own commit, and drop
the two `.gitignore` cover-up rules at the same time.

### 5.3 Files with `parents[N]` that this plan does NOT move — leave untouched

`src/common.py:21` (`parents[1]`, depth 1 — correct);
all of `scripts/{osd462,v11,v13,clinical_axes,grey60,regulator_activity,stage0,subtype_reference,celltype,integration}/*` (`parents[2]`, depth 2 — correct);
`src/{v11,v13,multiomics}/*` (`parents[2]` — correct);
`tests/*` (`parents[1]` — correct).
Moving `scripts/osd253_strain_module_projection.py` into `scripts/grey60/`
changes its depth 1→2, but that file contains **no** `parents[N]` expression
(`grep -n "parents\[" scripts/osd253_strain_module_projection.py` → no match);
it relies on CWD = repo root, which is unchanged.

---

## 6. Documented-command impact

Extraction: `grep -rn -E "(python|Rscript|python3)[^ ]* +[^ ]*(scripts|src)/[A-Za-z0-9_./-]+\.(py|R)" README.md docs/*.md *.md`.

### 6.1 Commands that KEEP WORKING — zero changes needed

| doc:line | command | why safe |
|---|---|---|
| `README.md:200` | `venv/bin/python scripts/osd462/08_stage0_provenance_audit.py` | `scripts/osd462/` untouched |
| `README.md:201` | `venv/bin/python scripts/osd462/09_stage0_manuscript_reporting.py` | untouched |
| `README.md:213` | `PYTHONPATH=. venv/bin/python scripts/v13/run_continuous_phospho_inference.py --config config/dct_asdn_phospho_reanalysis.yaml …` | `scripts/v13/` and `config/` untouched |
| `README.md:223` | `venv/bin/python scripts/v13/build_inference_report.py …` | untouched |
| `README.md:226` | `venv/bin/python scripts/v13/finalize_inference_provenance.py …` | untouched |
| `README.md:234` | `latexmk -pdf … manuscript_v13.tex` | `latex_paper/` untouched |
| `README.md:269` | `venv/bin/python -m pytest -q` | rootdir/CWD unchanged; `src/` still exists |
| `docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md:322–330` (9 commands) | `scripts/clinical_axes/*` | untouched |
| `docs/DCT_SUBTYPE_REFERENCE_REBUILD_2026-07-28.md:101,105,109,113,117,126` (6) | `scripts/subtype_reference/*` | untouched |
| `docs/OSD462_STAGE0_PROVENANCE_AUDIT_2026-07-28.md:11`, `docs/OSD462_STAGE0_MANUSCRIPT_REPORTING_2026-07-29.md:11` | stage-0 scripts | untouched |
| `docs/PODOCYTE_STRICT_MATCHING_AUDIT_2026-08-11.md:115` | `scripts/clinical_axes/run_strict_podocyte_matching_audit.py` | untouched |
| `docs/v13_continuous_phospho_inference.md:149,162,170,204`, `docs/v13_reporting.md:12` | `scripts/v13/*` | untouched |
| `docs/v11_layer_specificity_analysis_log_2026-06-07.md:361,364` | `scripts/v11/*` | untouched |

**Every command in `README.md`'s Reproduction section survives this plan.**

### 6.2 Commands and file references that BREAK

| doc:line | reference | breaks because | fix |
|---|---|---|---|
| `CLOUD_DEPLOYMENT.md:161` | `!python src/run_all_phases.py --phases 2 3 5 6 7 --run-id terra_run_001` | L1 moves `run_all_phases.py` | rewrite path |
| `CLOUD_DEPLOYMENT.md:178` | `python src/run_all_phases.py --dry-run` | same | rewrite |
| `CLOUD_DEPLOYMENT.md:255` | `python src/run_all_phases.py \ …` | same | rewrite |
| `FILE_SUMMARY.md:128` | `python src/run_all_phases.py` | same | rewrite |
| `FILE_SUMMARY.md:131` | `python src/run_all_phases.py --phases 2 3 5` | same | rewrite |
| `FILE_SUMMARY.md:134` | `python src/run_all_phases.py --skip-r` | same | rewrite |
| `FILE_SUMMARY.md` §"Reorganized Module Structure" | whole `src/preprocessing`, `src/networks`, `src/statistics`, `src/enrichment`, `src/validation`, `src/visualization` inventory + the `python -m src.networks.shared_topology` usage block | L1/L2 | rewrite the section; it is already stale (it documents `src/utils.py`, which does not exist) |
| `docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:197,425` | `scripts/run_full_pipeline.py` | §4.4a | rewrite path (2 refs) |
| `docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:15,383` | `node2vec_asgsr_osd771_kidney_transcriptomics.pdf` | §4.1 | rewrite path (2 refs) |
| `docs/provenance_audit_2026-06-06.md:110` | `scripts/run_full_pipeline.py` | §4.4a | rewrite |
| `docs/NEXT_PAPER_DIRECTION_AFTER_GREY60_2026-07-29.md:196` | `scripts/run_contrast_vector_framework.py` | §4.4b | rewrite |
| `docs/statistical_assessment_v5.md:4` | `scripts/run_lar_reversal_analysis.py`, `scripts/plot_lar_reversal_dashboard.py` | §4.4b | rewrite (2 refs) |
| `docs/annotated_issues_v5.md:32,37,42,47` | `scripts/run_lar_reversal_analysis.py:<line>` | §4.4b | rewrite (4 refs); line numbers stay valid |
| `docs/STATE_OF_PLAY_2026-07-29.md:134` | `email_to_casaletto_draft.md` | §4.1 | rewrite |
| `docs/DCT_SUBTYPE_REFERENCE_REBUILD_2026-07-28.md:33` | `scripts/build_hybrid_reference.R` | only if §4.4e is overridden | leave the file in place instead |
| `COMMENT_ARCHIVE.md` (~20 commands: `:22,42,431,456,476,507,530,550,569,585,626,694,900,1073,1100,1147,1670,1829,2090,2119–2122,2992`) | mixed | untracked, gitignored, explicitly historical | **do not update** — record it as knowingly stale |
| `docs/_commentary_scripts.md` (65 `## scripts/…` headings + line ranges) | every moved script | concurrent wave's worklist | coordinate (§0) |
| `docs/_commentary_src.md` (all `## src/…` headings for moved packages) | every moved `src/` module | same | coordinate (§0) |

### 6.3 Non-doc references that break

| location | reference | fix |
|---|---|---|
| `Dockerfile:73` | `ENTRYPOINT ["python3", "src/run_all_phases.py"]` | rewrite to the legacy path, or repoint at a current workflow |
| `config/revisions/network-rewiring.yaml:42,51` | `entry: src/run_all_phases.py` | rewrite — **`tests/test_rrrm2_registry.py::test_stage_entries_exist_on_disk` fails otherwise** |
| `config/revisions/contrast-vectors.yaml:37,46` | `entry: scripts/run_contrast_vector_framework.py`, `scripts/audit_stability_gates.py` | rewrite — same test |
| `config/contrast_vector_framework.yaml:6` | comment naming the orchestrator | rewrite (cosmetic) |
| `config/dct_subtype_reference_freeze_v1.yaml:91` | `cluster_map_source: scripts/build_hybrid_reference.R` | only if §4.4e is overridden |
| `scripts/grey60/04_run_external.py:15` | `from scripts.osd253_strain_module_projection import …` | rewrite to `scripts.grey60.…` — **no test covers this** |
| `src/run_all_phases.py` × 19 strings | see §4.2a | rewrite |
| `scripts/run_contrast_vector_framework.py:1000,1037,1061` | `REPO_ROOT / "scripts" / "<name>.py"` × 3 | rewrite to the contrast_vectors subdir |
| `scripts/run_phase3_pipeline.py:94,106,125` | `"scripts/phase3_*.py"` CWD-relative × 3 | rewrite |
| `scripts/run_phase4_7_pipeline.py:26` | `SCRIPTS = REPO_ROOT / "scripts"` | repoint at the legacy scripts dir |
| `src/preprocessing/deconvolution.R:149`, `scripts/run_deconvolution.R:144` | `stop("… Run scripts/export_raw_counts.py")` | rewrite (message text) |

### 6.4 Test imports that constrain the plan

`grep -rn 'import\|from' tests/` resolved through the AST graph. Every module
path referenced by `tests/`:

- `src.clinical_axes.{analysis,data,matching,osd656_context,statistics}` — **stay**
- `src.multiomics.{celltype_panels,osd462_anchor,osd462_stage0,phenotype_anchor,regulator_activity}` — **stay**
- `src.v11.{aldosterone_axis,cmap_screen,core_analysis,dct_continuous_gradient,h2_composition_aware_phospho,human_concordance,kinome_atlas_ksea,matched_null,observability_audit,publication_figures,recurrence_meta,rna_protein_propagation}` — **stay**
- `src.v13.{compartment_adversarial_audit,continuous_phospho_inference,reporting}` — **stay**
- `src.subtype_reference.reference_builder`, `src.grey60.adversarial` — **stay**
- `src.networks.{contrast_vectors,cross_osdr_projection,edge_regression,external_aging_axis,lar_reversal,lioness,manuscript_decision,mechanism_axis,procrustes,shared_topology,stability_test,tubulointerstitial_state}` — **stay** (12 files)
- `src.statistics.{full_pipeline_permutation,full_regression,permutation_bootstrap,silent_shifters}` — **stay**
- `src.validation.{cross_validation,enhanced_cv,external_replication,osd_external_validation}` — **stay**
- `scripts.clinical_axes.run_compartment_context` — **stay** (`tests/test_clinical_axes_integration.py:8`)
- config paths: `config/{clinical_renal_axes_cross_mission,hyperparameters,dct_asdn_phospho_reanalysis,dct_subtype_reference_freeze_v1}.yaml` — **all stay in `config/`**

Note there is **no `conftest.py`, no `tests/__init__.py`, and (until today) no
`pytest.ini`**. `import src.*` works only because `python -m pytest` puts the
CWD on `sys.path[0]`. The plan does not change that, but do not "simplify" the
documented invocation to bare `pytest`.

---

## 7. The single biggest breakage trap

**`config/revisions/*.yaml` is a machine-asserted path registry that did not
exist when this reorganization was scoped, and it is enforced by a test.**

```
tests/test_rrrm2_registry.py:38  def test_stage_entries_exist_on_disk(rev_id):
tests/test_rrrm2_registry.py:40      assert stage.entry_path.exists(), f"{rev_id}/{stage.id}: missing {stage.entry}"
tests/test_rrrm2_registry.py:48      assert (REPO_ROOT / rev.spec_config).exists(), rev.spec_config
tests/test_rrrm2_registry.py:50      assert (REPO_ROOT / doc).exists(), f"{rev_id}: missing doc {doc}"
```

`rrrm2/registry.py:36` resolves `entry_path` as `REPO_ROOT / self.entry`, so
every one of the **41 registered `entry:` paths** and every `docs:` entry across
nine revision files is a hard constraint. Three specific consequences that a
path-grep of `README.md` and `docs/` would miss entirely:

1. `network-rewiring.yaml` lists `METHODOLOGY.md` and `FILE_SUMMARY.md` in its
   `docs:` array. Both are **gitignored root files**. They look like deletable
   root junk; deleting either turns the suite red.
2. `contrast-vectors.yaml` lists `agents_instruction.md` as `docs:[0]`. Same
   trap for the largest root markdown file.
3. `network-rewiring.yaml` registers `src/run_all_phases.py` twice — the single
   file this plan most wants to archive.

**Runner-up trap**, and the one specific to code rather than metadata: the
function-local import at `src/multiomics/osd462_anchor.py:534`,
`from src.networks.contrast_vectors import cosine`. The locked v13 inference
engine reaches into the retired Workflow-A network package through an import
that never appears at module top level, so any audit built on top-level imports
concludes `src/networks/` is fully retired. It is not.

**Third**, and the one this brief specifically predicted: eleven
depth-sensitive files spell their repo-root computation as `.parent.parent`
rather than `parents[N]` (§5) — including `src/run_all_phases.py:14` itself —
so a `grep "parents\["` sweep silently misses them, and two files
(`src/visualization/network_diagnostics.py:34`,
`src/preprocessing/data_alignment.py:14`) already have the wrong index today,
with `src/plots/` and `src/data/results/` on disk as the proof.

---

## 8. Delete / gitignore table

### 8.1 Deletes

| path | tracked | ignored | referenced | recoverable from git? | command |
|---|---|---|---|---|---|
| `walkthrough.md` | no | `/walkthrough.md` | none | no | `rm walkthrough.md` — content is an OpenLiveStacker/Rust-daemon walkthrough from an unrelated project |
| `results.txt` | no | `/results.txt` | none | no | `rm results.txt` — regenerable Workflow-A run comparison |
| `.chrome_paper_doi_links.txt.swp` | no | `*.swp` | none | no | `rm` |
| `src/networks/.edge_regression.py.swp`, `.embeddings.py.swp`, `.lioness.py.swp` | no | `*.swp` | none | no | `rm src/networks/.*.swp` |
| 41 × `.DS_Store` | no | `.DS_Store` | none | no | `find . -name .DS_Store -not -path './venv/*' -not -path './.git/*' -delete` |
| `.idea/` (7 files) | no | `.idea/` | none | no | `rm -rf .idea` |
| `deconv_music_out/` | no | yes | prose only (`docs/v11_novelty…:57` says it is empty) | n/a (empty) | `rmdir` |
| `tms_kidney_female_counts/` | no | yes | none | n/a (empty) | `rmdir` |
| `tms_kidney_ytrfemale_counts/` | no | yes | none | n/a (empty) | `rmdir` |
| `plots/` | no | `/plots/` | as an output root | n/a (only `.DS_Store`) | `rm -rf plots` — recreated on demand |
| `scripts/data/` tree | no | no rule, but 0 files | none | n/a | `rm -rf scripts/data` — stray empty dirs from a wrong-CWD run |
| `src/plots/` | no | `/src/plots/` | none | no | `rm -rf src/plots` — artifact of the §5.2 bug |
| `src/data/results/` | no | `/src/data/results/` | none | no | `rm -rf src/data/results` — same |
| `scripts/sample.py` | no | `/scripts/sample.py` | none | no | `rm` |
| `scripts/phase7_2` | no | `/scripts/phase7_2` | none | no | `rm` |
| `notebooks/sample.py`, `notebooks/notebook00.ipynb`, `notebooks/__pycache__` | no | 2 explicit rules | none | no | `rm -rf notebooks` (optional) |
| **`scripts/run_permutations.py`** | **yes** | no | none | **yes** | `git rm scripts/run_permutations.py` — 0 bytes, blob `e69de29…`, tracked since 2025-12-15, never populated |

**57 paths total.** Only one is a tracked deletion.

### 8.2 gitignore-only changes (no file moves)

| rule | action | reason |
|---|---|---|
| `/results.txt` | remove | file deleted |
| `/walkthrough.md` | remove | file deleted |
| `/permutedStats-actualModules.RData` | remove | file moved into already-ignored `data/external/`; the generic `*.RData` rule still covers it |
| `/chrome_paper_doi_links.txt` | remove | file moved into already-ignored `tmp/` |
| `/node2vec_asgsr_osd771_kidney_transcriptomics.pdf` | remove | file moved into `docs/` and tracked (it is cited twice as the founding artifact) |
| `/notebooks/notebook00.ipynb`, `/notebooks/sample.py` | replace with `/notebooks/` | simpler, and covers `notebooks/__pycache__` |
| `/src/plots/`, `/src/data/results/` | remove **after** the §5.2 index fixes land | these rules exist only to hide the output of two path bugs |
| `latex_paper/*` + `!latex_paper/manuscript_v11.tex` + `!latex_paper/manuscript_v11.pdf` | add `!latex_paper/manuscript_v13.tex` and `!latex_paper/manuscript_v13.pdf` | `README.md:183,185` link to v13, which is currently **excluded from git** — a reviewer cloning the repo gets no manuscript. `docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md` §3E lists "dropping the `latex_paper/*` gitignore exclusion" as outstanding repo hygiene |
| `/scripts/sample.py`, `/scripts/phase7_2` | keep | harmless once the files are gone; cheap insurance |

---

## 9. Ordered execution script

Run `venv/bin/python -m pytest -q` after every batch. Expected: **225 passed**
(or 225 + whatever the concurrent registry wave adds). Commit per batch.

### Batch 0 — deletes and gitignore (touches no `.py`/`.R`; suite stays green)

```bash
cd /Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome

find . -name .DS_Store -not -path './venv/*' -not -path './.git/*' -delete
rm -f .chrome_paper_doi_links.txt.swp src/networks/.edge_regression.py.swp \
      src/networks/.embeddings.py.swp src/networks/.lioness.py.swp
rm -f walkthrough.md results.txt
rm -rf .idea plots src/plots src/data/results scripts/data
rmdir deconv_music_out tms_kidney_female_counts tms_kidney_ytrfemale_counts
rm -f scripts/sample.py scripts/phase7_2
git rm scripts/run_permutations.py

mkdir -p data/external/legacy_wgcna
mv permutedStats-actualModules.RData data/external/legacy_wgcna/
mv chrome_paper_doi_links.txt tmp/
```

**Follow-up edits in this batch:** apply the `.gitignore` changes of §8.2
*except* the `/src/plots/` and `/src/data/results/` removals (those wait for
Batch 4).

### Batch 1 — root-level document moves

```bash
git mv email_to_casaletto_draft.md docs/email_to_casaletto_draft_2026-06-19.md
git mv manuscript_v11_prose_audit.md docs/
git mv manuscript_v12_prose_audit.md docs/
git mv latex/repository_statistical_audit_2026_04_24.tex docs/
rmdir latex
mv node2vec_asgsr_osd771_kidney_transcriptomics.pdf \
   docs/founding_proposal_2025-12-09_node2vec_asgsr_osd771.pdf
git add docs/founding_proposal_2025-12-09_node2vec_asgsr_osd771.pdf
```

**Follow-up edits:** `docs/STATE_OF_PLAY_2026-07-29.md:134`;
`docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md:15` and `:383`. No test is affected
(no revision file lists any of these in `docs:` — verified by reading all nine).

### Batch 2 — `scripts/` regrouping that does not touch Workflow A

```bash
mkdir -p scripts/sf_classifier
git mv scripts/spaceflight_kidney_classifier.py   scripts/sf_classifier/
git mv scripts/sf_classifier_component_decomp.py  scripts/sf_classifier/
git mv scripts/sf_classifier_negctrl_figure.py    scripts/sf_classifier/

git mv scripts/osd253_strain_module_projection.py scripts/grey60/
```

**Follow-up edits, mandatory:**
- `scripts/grey60/04_run_external.py:15` →
  `from scripts.grey60.osd253_strain_module_projection import (…)`.
  **No test exercises this file** — verify manually with
  `PYTHONPATH=. venv/bin/python -c "import importlib.util,sys; …"` or a `--help` run.
- The three `sf_classifier` scripts compute their own `ROOT`; confirm
  `data/results/run_20260701_sf_classifier` still resolves after the depth
  change (they use `os.path.join(ROOT, "data", "results", …)` at
  `spaceflight_kidney_classifier.py:35`, `sf_classifier_component_decomp.py:16`,
  `sf_classifier_negctrl_figure.py:19`).
- `docs/_commentary_scripts.md:799,823,848` headings.

### Batch 3 — Workflow A scripts → `legacy/`

```bash
mkdir -p legacy/workflow_a_network_rewiring/scripts/contrast_vectors legacy/scratch
L=legacy/workflow_a_network_rewiring/scripts

git mv scripts/phase2_edge_regression.py scripts/phase3_node2vec_embedding.py \
       scripts/phase3_node_rewiring_from_deltas.py scripts/phase3_procrustes_rewiring.py \
       scripts/phase4_anchor_qc_report.py scripts/phase5_build_silent_shifters_strict.py \
       scripts/phase5_derive_interaction_persistence.py scripts/phase6_edge_regression_full.py \
       scripts/phase6_perm_bootstrap_node_rewiring.py scripts/phase7_grounding_fast.py \
       scripts/build_phase2_skeleton.py scripts/compute_phase2_lioness_on_skeleton.py \
       scripts/run_phase1_networks.py scripts/run_phase2_pipeline.py \
       scripts/run_phase3_pipeline.py scripts/run_phase4_7_pipeline.py \
       scripts/run_full_pipeline.py scripts/pick_anchors.py \
       scripts/make_consensus_anchors.py scripts/plot_output.py \
       scripts/plot_skeleton_diagnostics.py scripts/batch_embeddings.sh \
       scripts/align_metadata_to_counts.py scripts/export_raw_counts.py \
       scripts/build_id_map.py scripts/discover_dct_markers.py \
       scripts/audit_fix_status.py scripts/DESeq2.R scripts/run_gene_level_DE.R \
       scripts/export_phase1_to_python.R scripts/phase1_to_python.R \
       scripts/run_deconvolution.R  $L/

git mv scripts/run_contrast_vector_framework.py scripts/audit_stability_gates.py \
       scripts/run_lar_reversal_analysis.py scripts/run_mechanism_axis_prioritization.py \
       scripts/run_tubulointerstitial_state_analysis.py \
       scripts/plot_lar_reversal_dashboard.py scripts/plot_manuscript_v4_figures.py \
       $L/contrast_vectors/

git mv scripts/analyze_dct_variance.py scripts/audit_reference_labels.py \
       scripts/check_dct_universe.py scripts/check_reference_counts.py \
       scripts/compare_two_runs.py scripts/debug_deconvolution_phase0.py \
       scripts/explore_chen_atlas.R scripts/extract_cell_types_from_csv.py \
       scripts/extract_cell_types_to_file.py scripts/extract_ecm_edges.py \
       scripts/inspect_10x_reference.py scripts/inspect_sce_structure.R \
       scripts/map.R scripts/map_top_hits.py scripts/pathway_presence.py \
       scripts/run_pseudobulk_test.R scripts/sanity_marker_panel.R \
       scripts/targeted_fdr.py scripts/module_convergence_drilldown.py \
       scripts/module_convergence_fisher.py  legacy/scratch/
```

**Follow-up edits, mandatory — the suite will not go green without them:**
1. All 26 depth fixes in §5.1 for these files.
2. `config/revisions/contrast-vectors.yaml:37` and `:46` → new paths.
   *(`tests/test_rrrm2_registry.py::test_stage_entries_exist_on_disk` fails otherwise.)*
3. `src/run_all_phases.py:933` → new path for `run_contrast_vector_framework.py`.
4. `scripts/run_contrast_vector_framework.py:1000,1037,1061` → `contrast_vectors/` sub-paths.
5. `run_phase3_pipeline.py:94,106,125`; `run_phase4_7_pipeline.py:26`.
6. Message-text refs: `scripts/run_deconvolution.R:144`,
   `src/preprocessing/deconvolution.R:149`, `audit_stability_gates.py:103`,
   `config/contrast_vector_framework.yaml:6`.
7. Doc refs: §6.2 rows for `run_full_pipeline`, `run_contrast_vector_framework`,
   `run_lar_reversal_analysis`, `plot_lar_reversal_dashboard`.

`scripts/build_hybrid_reference.R` is deliberately left in `scripts/` (§4.4e).
After this batch `scripts/` holds exactly that one loose file plus the eleven
experiment subdirectories.

### Batch 4 — Workflow A `src/` packages → `legacy/` (Batch L1)

```bash
mkdir -p legacy/workflow_a_network_rewiring
W=legacy/workflow_a_network_rewiring

git mv src/preprocessing  $W/preprocessing
git mv src/enrichment     $W/enrichment
git mv src/visualization  $W/visualization
git mv src/markers        $W/markers
git mv src/data           $W/idmap          # `src/data/results/` already deleted in Batch 0
git mv src/run_all_phases.py $W/run_all_phases.py
```

**Follow-up edits, mandatory:**
1. Create `legacy/__init__.py` and `legacy/workflow_a_network_rewiring/__init__.py`
   (empty) so `python -m legacy.workflow_a_network_rewiring.<pkg>.<mod>` works.
2. `run_all_phases.py`: `.parent.parent` → `parents[2]`; rewrite the 13 module
   strings and 6 rscript paths listed in §4.2a.
3. Depth fixes §5.1 for `deconvolution_sanity.py`, `regression_enrichment.py`,
   `publication_plots.py`, `publication_figures.py`.
4. **Behaviour-changing fixes §5.2**: `network_diagnostics.py:34` and
   `data_alignment.py:14`. Then remove the `/src/plots/` and
   `/src/data/results/` rules from `.gitignore`.
5. `src.enrichment.gene_set_loader` imports inside `biological_grounding.py`,
   `regression_enrichment.py`, and the docstring at `gene_set_loader.py:13`.
6. `scripts/build_id_map.py` and `scripts/discover_dct_markers.py` (now in
   `legacy/.../scripts/`): both the `dirname(__file__), ".."` fix and the
   `src.data.build_id_map` / `src.markers.discover_dct` import paths.
7. `Dockerfile:73` ENTRYPOINT.
8. `config/revisions/network-rewiring.yaml:42` and `:51`
   *(registry-asserted — the suite fails otherwise).*
9. `CLOUD_DEPLOYMENT.md:161,178,255`; `FILE_SUMMARY.md:128,131,134` and its
   module-structure section.
10. Error-message text at `src/networks/shared_topology.py:526,603` and
    `src/enrichment/biological_grounding.py:257` (accuracy only).

**Do not proceed to a Batch 5 that moves `src/networks/`, `src/statistics/`, or
`src/validation/`.** §3.3 lists the 22 blocked files and the 18 test files that
would have to be rewritten, plus the v13 dependency in §3.2.

### Batch 5 (OPTIONAL, not recommended) — partial-package files, Batch L2

Only if the owner accepts splitting three packages across two trees. Moves and
follow-ups are enumerated in §4.3 and the `†` entries of §4.2a. Requires three
new `__init__.py` files under `legacy/workflow_a_network_rewiring/`.

---

## 10. Open questions and things this plan deliberately does not decide

- **`scripts/build_hybrid_reference.R`** — evidence for moving: it is a
  reference-construction script and every other reference script lives in
  `scripts/subtype_reference/`. Evidence against: `config/dct_subtype_reference_freeze_v1.yaml:91`
  pins its path inside a *frozen* spec, and `docs/DCT_SUBTYPE_REFERENCE_REBUILD_2026-07-28.md:33`
  says the mapping "is frozen in the config". I recommend leaving it, but this
  is an owner call, not a technical one.
- **Batch L2** — genuinely ambiguous. Moving 11 files splits three packages;
  leaving them keeps three retired modules in `src/`. Evidence both ways is in
  §3.3 and §4.3.
- **`node2vec_…pdf` tracking** — moving it into `docs/` makes a 231 KB binary
  tracked for the first time. Justified by the two audit citations, but it is a
  git-content change, not just a reorganization.
- **Whether `.idea/` and `notebooks/` should be deleted at all** — both are
  purely local; deleting them is a convenience, not a correctness improvement.
- **`docs/_commentary_scripts.md` / `docs/_commentary_src.md` re-keying** — owned
  by the concurrent wave (§0), not by this plan.
