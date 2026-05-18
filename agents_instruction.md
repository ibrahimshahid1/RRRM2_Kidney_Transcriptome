# Cross-OSDR Network-Contrast Framework — Working Plan

**Status:** Working reference for Claude Code agents implementing the new analysis pivot.
**Owner:** Ibrahim Shahid.
**Date locked:** May 2026.
**Working scope:** Minimum-viable (RRRM-2 + OSD-513 + OSD-253). Full-scope extension only after MV is shipped and signal is confirmed.

This document is the canonical reference for every agent working on this analysis. Read it before opening a code editor. When in doubt about a methodological choice, this document overrides earlier drafts (`manuscript_v2.tex`, `rrrm2_network_wgcna_consolidated_reference.pdf`, `remediation_guide.tex`). Earlier documents are kept for historical context; this one supersedes them.

---

## 1. The scientific question (do not rewrite this)

**Plain English:** After accounting for kidney cell composition, does spaceflight *amplify*, *dampen*, *reverse*, or *redirect* age-related kidney transcriptional network change in mouse? Separately for ISS-T (terminal on-orbit) and LAR (live animal return).

**Formal estimand (contrast of contrasts):**

For each arm $a \in \{\text{ISS-T}, \text{LAR}\}$ and each scale of network summary $\mathbf{N}(\cdot)$:

$$
\mathbf{A}^a_{GC} \;=\; \mathbf{N}(\text{Old}, a, GC) - \mathbf{N}(\text{Young}, a, GC) \quad \text{(control aging vector)}
$$

$$
\mathbf{A}^a_{FLT} \;=\; \mathbf{N}(\text{Old}, a, FLT) - \mathbf{N}(\text{Young}, a, FLT) \quad \text{(flight aging vector)}
$$

$$
\boldsymbol{\Delta}^a_{\text{AgeRewire}} \;=\; \mathbf{A}^a_{FLT} - \mathbf{A}^a_{GC} \quad \text{(within-arm Flight} \times \text{Age interaction)}
$$

$$
\boldsymbol{\Delta}_{\text{Arm}} \;=\; \boldsymbol{\Delta}^{\text{ISS-T}}_{\text{AgeRewire}} - \boldsymbol{\Delta}^{\text{LAR}}_{\text{AgeRewire}} \quad \text{(three-way interaction)}
$$

**Decomposition:** project $\mathbf{A}^a_{FLT}$ onto $\mathbf{A}^a_{GC}$:

$$
\beta_a = \frac{\mathbf{A}^a_{FLT} \cdot \mathbf{A}^a_{GC}}{\mathbf{A}^a_{GC} \cdot \mathbf{A}^a_{GC}}
\quad,\quad
\cos\theta_a = \frac{\mathbf{A}^a_{FLT} \cdot \mathbf{A}^a_{GC}}{\|\mathbf{A}^a_{FLT}\| \, \|\mathbf{A}^a_{GC}\|}
\quad,\quad
\mathbf{R}_a = \mathbf{A}^a_{FLT} - \beta_a \mathbf{A}^a_{GC}, \quad \rho_a = \frac{\|\mathbf{R}_a\|}{\|\mathbf{A}^a_{FLT}\|}
$$

**Interpretation taxonomy:**

| Statistic | Value | Verbal interpretation |
|---|---|---|
| $\beta_a$ | $>1$ | flight **amplifies** control aging direction |
| $\beta_a$ | $0 < \beta_a < 1$ | flight **dampens** control aging direction |
| $\beta_a$ | $<0$ | flight **reverses** control aging direction |
| $\beta_a$ | $\approx 0$ | flight **redirects** (not aligned with aging) |
| $\rho_a$ | small | flight aging mostly along control aging axis |
| $\rho_a$ | large | flight aging mostly redirected into a new direction |

This is the entire intellectual point of the paper. Every method choice below serves it.

---

## 2. Statistical framework: the five guardrails

Naive computation of $\beta$ and $\cos\theta$ at $n=5+5$ per cell will fail. The guardrails below are not optional; they are the framework.

### 2.1 Guardrail A — Stability test gate (run first, before anything else)

Before *any* projection or decomposition, test whether $\mathbf{A}^a_{GC}$ is itself a stable estimated direction.

**Procedure:**

1. For each arm $a$, compute $\mathbf{A}^a_{GC}$ on the full RRRM-2 GC samples.
2. Bootstrap-resample mice within each Age × Arm × GC cell, $B = 2000$ times.
3. For each bootstrap replicate $b$, recompute $\mathbf{A}^{a,(b)}_{GC}$.
4. Compute the angular stability:

$$
\text{stability}_a = \text{median}_b \cos(\mathbf{A}^a_{GC}, \mathbf{A}^{a,(b)}_{GC})
$$

5. Compute the 95% angular envelope: the 2.5th percentile of the bootstrap cosines with the full-sample direction.

**Pass criterion (pre-registered):**

- **Edge-level** $\mathbf{N}$: stability median $\geq 0.30$ AND 2.5th percentile $\geq 0.10$ — *if edge-level fails, edge-level projection is DEMOTED to exploratory.*
- **Gene-level node-aggregate** $\mathbf{N}$: stability median $\geq 0.40$ AND 2.5th percentile $\geq 0.20$ — *if gene-level fails, gene-level projection is DEMOTED to exploratory.*
- **Module-level** $\mathbf{N}$: stability median $\geq 0.60$ AND 2.5th percentile $\geq 0.40$ — *if module-level fails, abandon the within-RRRM-2 projection layer entirely and use only external-axis projection (Guardrail E).*

**Output artifact:** `data/results/<run>/contrast_vectors/agc_stability_report.tsv` — must be checked before any downstream projection runs.

### 2.2 Guardrail B — Resolution hierarchy (primary → secondary → exploratory)

Project at multiple resolutions and lead with the resolution that *passes* the stability gate. The resolution hierarchy, ordered from most-robust to least-robust at $n=5/$cell:

| Resolution | Vector | Default role |
|---|---|---|
| **Module activity** (eigengene or mean-of-module) | 20-30 dim | **Primary** |
| **Pathway activity** (curated 8-panel + GO/Reactome 50-set) | 50-100 dim | **Primary** |
| **Gene-level node summaries** (differential connectivity, $\Delta$kME, centrality change) | 2500 dim | **Secondary** |
| **LIONESS-derived per-sample node score aggregate** | 2500 dim | **Secondary** |
| **Edge-level** (per-edge $\Delta z$) | $\sim 180{,}000$ dim | **Exploratory only** |

Do not promote a lower-resolution result to primary unless it passes the stability gate.

### 2.3 Guardrail C — Shrinkage / precision weighting

Raw dot products weight all features equally. Noisy features should contribute less.

For each feature $k$ (edge, gene, or module) at each level, estimate the within-condition variance $\hat\sigma^2_k$ from bootstrap resampling of GC samples. Define weight:

$$
w_k = \frac{1}{\hat\sigma^2_k + \epsilon}
$$

where $\epsilon$ is a small floor (e.g. the 10th percentile of $\hat\sigma^2$). Then compute weighted projection:

$$
\beta_a^{(w)} = \frac{\mathbf{A}^a_{FLT}{}^\top W \, \mathbf{A}^a_{GC}}{\mathbf{A}^a_{GC}{}^\top W \, \mathbf{A}^a_{GC}}
\quad,\quad
\cos\theta_a^{(w)} = \frac{\mathbf{A}^a_{FLT}{}^\top W \, \mathbf{A}^a_{GC}}{\sqrt{\mathbf{A}^a_{FLT}{}^\top W \mathbf{A}^a_{FLT}} \, \sqrt{\mathbf{A}^a_{GC}{}^\top W \mathbf{A}^a_{GC}}}
$$

where $W = \text{diag}(w_1, \ldots, w_K)$.

Report both unweighted and weighted versions; lead with weighted in the manuscript.

### 2.4 Guardrail D — Bootstrap the whole pipeline (no point estimates)

Every output — $\beta_a, \cos\theta_a, \rho_a, \|\boldsymbol{\Delta}^a_{\text{AgeRewire}}\|$ — must be reported with a bootstrap confidence interval, not a single number.

**Procedure:**

1. Stratified bootstrap resample mice within each Age × Arm × EnvGroup cell, $B = 2000$ times.
2. For each bootstrap replicate, recompute the full chain: residualization (or held-out residuals if computed once), backbone projection (held fixed), $\mathbf{A}^a_{GC}$, $\mathbf{A}^a_{FLT}$, $\beta_a$, $\cos\theta_a$, $\rho_a$, $\boldsymbol{\Delta}^a_{\text{AgeRewire}}$.
3. Report median and 95% percentile CI for each statistic at each resolution.
4. **Categorize the interpretation by the bootstrap distribution, not by the point estimate.** If 70% of bootstrap replicates fall in the "amplify" range ($\beta > 1$) and 25% in "dampen" ($0 < \beta < 1$), the honest claim is "predominantly amplification with non-trivial dampening uncertainty" — NOT "amplification."

**Output artifact:** `data/results/<run>/contrast_vectors/<arm>/<resolution>/bootstrap_decomposition.tsv` with columns `statistic, point_estimate, median, ci_low, ci_high, frac_amplify, frac_dampen, frac_reverse, frac_redirect`.

### 2.5 Guardrail E — External aging axis (the cleanest fix when available)

The cleanest version of the analysis defines the aging axis from an **independent** kidney aging reference, not from RRRM-2 controls.

**Implementation:**

1. Use Tabula Muris Senis kidney female bulk-or-pseudobulk aging signature as $\mathbf{A}^{\text{ext}}_{\text{age}}$ at gene/module/pathway level.
2. Project RRRM-2 flight effects $\mathbf{F}^a = \mathbf{N}(FLT, a) - \mathbf{N}(GC, a)$ (pooled over age, separately per arm) onto $\mathbf{A}^{\text{ext}}_{\text{age}}$:

$$
\beta^{\text{ext}}_a = \frac{\mathbf{F}^a \cdot \mathbf{A}^{\text{ext}}_{\text{age}}}{\mathbf{A}^{\text{ext}}_{\text{age}} \cdot \mathbf{A}^{\text{ext}}_{\text{age}}}
$$

3. Same decomposition: $\beta^{\text{ext}}_a$ classifies whether flight is aging-like (positive), anti-aging-like (negative), or orthogonal in each arm.

**Pass criterion:** if the within-RRRM-2 aging-axis stability test (Guardrail A) passes at module/pathway level, report both within-RRRM-2 and external-axis projections side by side. If within-RRRM-2 fails, use external-axis as the *only* aging-axis projection and explicitly state this in the methods.

This is the single most powerful clarity move available. Lean on it.

---

## 3. Dataset roster and roles

Locked roster for **minimum viable scope** (MV). All other datasets are optional extension only.

| Dataset | Accession | Design | MV role |
|---|---|---|---|
| RRRM-2 | OSD-771 | Female C57BL/6NTac kidney; full Age × Arm × EnvGroup factorial; n=80 (5/cell) | **Primary factorial bridge** |
| RR-23 | OSD-513 | Male C57BL/6J kidney; FLT vs HGC vs VGC; 16–17 wk; 38 day flight; n=9/group | **Powered FLT-vs-GC flight vector** for projection recurrence test |
| RR-7 | OSD-253 | Female C57BL/6J + C3H/HeJ kidney; 25-day and 75-day timepoints; original GC and 2021 white-LED rerun GC | **Strain/duration/light sensitivity** |
| Tabula Muris Senis kidney | external | Single-cell + bulk kidney aging atlas across ages | **External aging axis** (Guardrail E) |

Extension scope (defer until MV ships):

| Dataset | Accession | Role |
|---|---|---|
| RR-1 | OSD-102 | Female C57BL/6J flight vector (context replication) |
| RR-3 | OSD-163 | BALB/c flight vector (cross-strain context, Finch et al. 2025) |
| RR-10 | OSD-462 | Mechanistic context (p21-null subset, WT-only flight vector if usable) |
| MHU-3 | OSD-457 | Multi-tissue Nrf2-KO and WT (stress-response genotype context) |
| RRRM-1 | OSD-913 | Kidney miRNA-seq (regulatory layer; do not use as mRNA co-expression evidence) |

**Pooling rule (absolute):** do NOT pool raw expression across studies into a single correlation network. Pool contrast-level statistics (signed Stouffer Z, random-effects meta-analysis) only after within-study contrast vectors are estimated.

---

## 4. Pre-registered decision thresholds

These thresholds are locked before any contrast estimation runs. Do not retroactively adjust them.

### 4.1 Stability gates (Guardrail A)

See §2.1. Passing thresholds:
- Edge level: median cos ≥ 0.30, 2.5th pct ≥ 0.10
- Gene level: median cos ≥ 0.40, 2.5th pct ≥ 0.20
- Module/pathway: median cos ≥ 0.60, 2.5th pct ≥ 0.40

### 4.2 Significance thresholds (decomposition)

- $\beta_a$ is "significant" if its 95% bootstrap CI does NOT include 0 AND does NOT include 1 (so the interpretation category is unambiguous).
- $\cos\theta_a$ is "significant" if its 95% bootstrap CI does NOT include 0 AND has $|\cos\theta_a| > 0.2$ at the lower bound.
- $\rho_a$ "redirectedness" is "substantial" if median $\rho_a \geq 0.5$.

### 4.3 Cross-OSDR recurrence (external projection)

For each external study $s$ with FLT-vs-GC flight vector $\mathbf{F}_s$, compute $\cos(\mathbf{F}_s, \mathbf{F}_{\text{RRRM-2}, a})$ where $\mathbf{F}_{\text{RRRM-2}, a}$ is the RRRM-2 FLT-vs-GC flight vector (age-pooled, arm-specific).

**Recurrence claim** requires:
- Point estimate $\cos > 0.2$
- 95% bootstrap CI (over study-specific resampling) lower bound $> 0$
- Same sign as the within-RRRM-2 vector direction

### 4.4 Interpretation category assignment

For $\beta_a$ category assignment (amplify/dampen/reverse/redirect), report the *fraction* of the bootstrap distribution falling in each category:

- "amplify": $\beta > 1$
- "dampen": $0 < \beta < 1$
- "reverse": $\beta < 0$
- "redirect": $|\beta| < 0.1$ AND $\rho > 0.5$

The headline category is assigned only if it captures > 70% of the bootstrap distribution. Otherwise report as "mixed (X% amplify, Y% dampen, Z% reverse, W% redirect)."

---

## 5. Phase-by-phase implementation plan

### Phase 0 — Repository state audit and harmonization scaffolding

**Deliverable:** clean source tree, harmonized metadata schema, deconvolution + residualization re-run consistent across all MV datasets.

**Tasks:**

1. Audit existing source files for compatibility with new framework:
   - `src/networks/lioness.py` — already updated to `raw_ranknorm` (✓)
   - `src/networks/edge_regression.py` — already handles `lioness_edges.npy` (✓)
   - `src/networks/embeddings.py` — supports `signed_split` (✓; used as exploratory only in new framework)
   - `src/networks/wgcna_analysis.R` — keep, role changes to module annotation
   - `src/statistics/full_regression.py` — keep, label-permutation aggregator (✓)
   - `src/statistics/permutation_bootstrap.py` — keep, hierarchical FDR (✓)
   - `src/validation/external_replication.py` — keep, ranked-GSEA external testing (✓)

2. Build harmonized metadata table at `data/processed/multi_study_metadata.tsv` with columns: `study, sample_id, strain, sex, age_weeks, age_group, arm, envgroup, mission_duration_days, recovery_days, dissection_location, light_protocol, platform, read_length, library_batch, RIN`.

3. Pull OSD-513 and OSD-253 VST / count matrices to `data/external/osdr/<accession>/`. Verify ISA metadata is parsable.

4. Add Tabula Muris Senis kidney aging signature at `data/external/aging_reference/tms_kidney_female_aging_axis.tsv` — per-gene aging effect estimates (Old vs Young) on the same Ensembl ID space as the RRRM-2 panel.

5. New module: `src/networks/contrast_vectors.py` (skeleton; full implementation in Phase 4).

6. New module: `src/networks/stability_test.py` (skeleton; implementation in Phase 3).

7. Update `config/hyperparameters.yaml` with new section `contrast_vector_framework:` containing the pre-registered thresholds from §4.

8. Run `scripts/audit_fix_status.py` to confirm prior remediation (14 fixes) still passes.

**Acceptance:** all four MV datasets have residualized expression matrices, harmonized metadata, and a frozen Ensembl ID space at `data/processed/multi_study_gene_universe.tsv`.

### Phase 1 — Multi-study preprocessing and residualization

**Deliverable:** residualized expression matrices for each MV dataset on the common gene universe, with provenance.

**Tasks:**

1. Residualize each dataset independently. **Do not pool raw expression.** Per-study residualization model:

$$
X_{g,i,s} = \mu_{g,s} + C_{i,s} \alpha_{g,s} + T_{i,s} \gamma_{g,s} + \varepsilon_{g,i,s}
$$

where $C$ is estimated kidney cell-type composition (from MuSiC against TMS kidney reference) and $T$ is study-specific technical covariates (batch, library kit, RIN, read length, sequencing run, sex if applicable, dissection location if applicable).

2. Save residuals at `data/processed/<study>/residuals.tsv.gz` with metadata stamp.

3. Sanity check: per-study residualized PC1 should not separate by batch but should still capture treatment/age structure.

4. Common gene universe: intersect Ensembl IDs across all four MV studies, keeping only protein-coding genes with $\geq 1$ CPM in $\geq 20\%$ of samples in EVERY study.

**Acceptance:** four study-specific residual matrices on a frozen common gene universe (~10,000–15,000 genes expected after intersection).

### Phase 2 — Frozen backbone construction

**Deliverable:** a frozen partial-correlation skeleton on the common gene universe, built from controls-only pooled across studies (after batch correction at the residual level).

**Tasks:**

1. Build a control-only Ledoit–Wolf-shrinkage partial-correlation backbone on pooled control samples (HGC/VGC/Basal from each study). Pool only after residualization. Apply ComBat-seq with study as batch on the residual matrix before partial-correlation estimation.

2. Top-$k$ neighbor selection at $k = 80$.

3. Save backbone at `data/processed/networks/control_pooled_backbone/skeleton_edges.tsv` and `phase2_genes.txt`.

4. Sensitivity: also build per-study control backbones and report top-$k$ overlap. Backbone is "stable" if Jaccard ≥ 0.50 across study-specific backbones.

5. **Critical:** the backbone is locked before contrast estimation. Save commit hash + backbone hash to `data/processed/networks/control_pooled_backbone/manifest.json`.

**Acceptance:** one frozen control-only backbone with stability report. If per-study backbones disagree (Jaccard < 0.5), document and use the most conservative subset (recurring edges across all studies).

### Phase 3 — Control aging vector stability test (GATE)

**Deliverable:** stability report for $\mathbf{A}^a_{GC}$ at each resolution, gating subsequent phases.

**Tasks:**

1. Compute $\mathbf{A}^a_{GC} = \mathbf{N}(\text{Old}, a, GC) - \mathbf{N}(\text{Young}, a, GC)$ for $a \in \{\text{ISS-T}, \text{LAR}\}$ at each resolution:
   - Edge: per-edge $\Delta z = \mathrm{atanh}(r_{\text{Old}}) - \mathrm{atanh}(r_{\text{Young}})$ on the frozen backbone.
   - Gene-level: aggregate $|\Delta z|$ over edges incident on each gene.
   - Module: residualized eigengene difference per module (from RRRM-2-internal WGCNA on residuals).
   - Pathway: mean residualized expression difference per gene set in `config/gene_sets.yaml`.

2. Bootstrap $B = 2000$ replicates of mice within each Age × Arm × GC cell.

3. Compute angular stability per resolution per arm: median and 2.5th-percentile cosine of bootstrap replicate against full-sample direction.

4. Apply stability gates from §2.1. Write pass/fail per resolution to `data/results/<run>/contrast_vectors/agc_stability_report.tsv`.

5. **Branching decision encoded in code:** if module-level fails the gate, halt Phase 4 within-RRRM-2 projection and proceed to Phase 4 external-axis-only path.

**Acceptance:** stability report exists, branching decision is logged to `agc_stability_decision.json` with explicit pass/fail per resolution and arm.

### Phase 4 — Within-RRRM-2 contrast vector decomposition

**Deliverable:** per-resolution decomposition of $\mathbf{A}^a_{FLT}$ along $\mathbf{A}^a_{GC}$ with bootstrap CIs, classified by interpretation category.

**Tasks:**

1. For each resolution that passed the stability gate (Phase 3):
   - Compute $\mathbf{A}^a_{FLT}$.
   - Compute $\beta_a$ (unweighted), $\beta_a^{(w)}$ (variance-weighted per Guardrail C), $\cos\theta_a$, $\cos\theta_a^{(w)}$, $\rho_a$.
   - Bootstrap-recompute $B = 2000$ times. Save full bootstrap distribution to `data/results/<run>/contrast_vectors/<arm>/<resolution>/bootstrap.tsv`.
   - Apply interpretation thresholds (§4.4). Write category assignment per arm per resolution to `interpretation_categories.tsv`.

2. Permutation calibration: within each Age × Arm stratum, permute Old/Young labels $K = 5000$ times. Recompute $\beta_a, \cos\theta_a, \rho_a$ on permuted labels. Empirical $p$-value: fraction of permutations with $|\beta_a|$ at least as extreme as observed.

3. Save permutation distributions for the same statistics at the same paths.

**Acceptance:** per-arm per-resolution decomposition + bootstrap + permutation calibration. Each statistic has a point estimate, 95% CI, and empirical $p$-value.

### Phase 5 — Module/pathway annotation of redirected component

**Deliverable:** biological interpretation of which modules/pathways carry the redirected component $\mathbf{R}_a$ at each resolution.

**Tasks:**

1. For each arm, compute the redirected component $\mathbf{R}_a = \mathbf{A}^a_{FLT} - \beta_a \mathbf{A}^a_{GC}$ at gene level.

2. Test concentration of $\mathbf{R}_a$ in each WGCNA module: rank genes by $|R_{a,g}|$, test top decile vs each module by Fisher's exact, BH-correct across modules. Save to `data/results/<run>/contrast_vectors/<arm>/redirected_module_enrichment.tsv`.

3. Same for the eight curated pathway panel and for GO/Reactome pathways from the existing enrichment infrastructure.

4. Hub-gene interpretation: the top-loading genes of $\mathbf{R}_a$ are the biological candidates for what flight is "redirecting" toward. Report top 50 with biological annotation.

**Acceptance:** redirected-component enrichment table per arm, biological-annotation report.

### Phase 6 — External OSDR projection (cross-OSDR recurrence)

**Deliverable:** alignment of RRRM-2-derived flight network vectors with OSD-513 and OSD-253 flight vectors.

**Tasks:**

1. For each external study $s$, compute FLT-vs-GC flight vector at each resolution:

$$
\mathbf{F}_s = \mathbf{N}_s(FLT) - \mathbf{N}_s(GC/HGC)
$$

For OSD-253, compute four variants: C57BL/6J × {25d, 75d} × {original GC, 2021 white-LED rerun GC}; C3H/HeJ × {25d, 75d} × {original GC, 2021 white-LED rerun GC}. The light/strain/duration sensitivity is reported as a 2×2×2 = 8-cell matrix.

2. Compute three projections per external study per arm of RRRM-2:
   - $\cos(\mathbf{F}_s, \mathbf{F}_{\text{RRRM-2}, a})$ — flight direction recurrence
   - $\cos(\mathbf{F}_s, \boldsymbol{\Delta}^a_{\text{AgeRewire}})$ — alignment with the within-arm Flight × Age interaction direction
   - $\cos(\mathbf{F}_s, \mathbf{R}^a)$ — alignment with the redirected component

3. Bootstrap uncertainty over external-study sample resampling. Apply recurrence thresholds from §4.3.

4. Meta-analysis: combine study-level signed projections via signed Stouffer Z weighted by effective sample size. Leave-one-study-out robustness checks.

**Acceptance:** alignment matrices saved to `data/results/<run>/contrast_vectors/cross_osdr_alignment/`. Recurrence verdict (pass/fail per external study per resolution per RRRM-2 contrast type).

### Phase 7 — External aging axis projection (Guardrail E)

**Deliverable:** RRRM-2 flight effect projection onto Tabula Muris Senis kidney aging axis.

**Tasks:**

1. Define $\mathbf{A}^{\text{ext}}_{\text{age}}$ from TMS kidney female: per-gene effect of age on bulk-or-pseudobulk expression, restricted to common gene universe. Save to `data/external/aging_reference/tms_kidney_aging_axis.tsv`.

2. Compute RRRM-2 arm-specific flight vectors (age-pooled): $\mathbf{F}^a = \mathbf{N}_{\text{RRRM-2}}(FLT, a) - \mathbf{N}_{\text{RRRM-2}}(GC, a)$, $a \in \{\text{ISS-T}, \text{LAR}\}$.

3. Project: $\beta^{\text{ext}}_a$, $\cos\theta^{\text{ext}}_a$, $\rho^{\text{ext}}_a$. Bootstrap-resample TMS samples and RRRM-2 samples jointly.

4. Interpretation: same taxonomy as §1, applied to flight-vs-external-aging-axis alignment. Report whether flight is "aging-like" (positive alignment) or "anti-aging-like" (negative) or "orthogonal" in each arm.

5. This is the cleaner-than-within-RRRM-2 version of the headline question. **If Guardrail A (stability test) fails for module/pathway level, this becomes the only aging-axis projection in the paper.**

**Acceptance:** external-axis projection per arm with bootstrap CI and interpretation category.

### Phase 8 — Manuscript decision (Phase 6 decision tree)

**Deliverable:** decision artifact at `data/results/<run>/manuscript_decision.json` that maps observed results to manuscript implication per the pre-registered tree.

| Observed result | Manuscript implication |
|---|---|
| RRRM-2 within-cohort aging-axis projection passes stability AND $\beta_a$ or $\cos\theta_a$ is significant in $\geq 1$ arm, AND recurs in $\geq 1$ external cohort | **Strong network-aging paper.** Spaceflight modifies kidney aging network architecture. Headline framing: "Flight $X$ control aging direction in $Y$ arm" where $X$ is amplify/dampen/reverse/redirect. |
| RRRM-2 within-cohort is weak alone but external alignment (OSD-513 or TMS external axis) is significant | **Cross-OSDR flight-network biology paper.** RRRM-2 provides factorial bridge; external cohorts provide recurrence. Headline: "Flight aging-axis alignment is detectable cross-cohort even when single-cohort discovery is underpowered." |
| Only WGCNA module-activity layer is positive; no projection signal | **Modest module-activity paper** — write only if intellectually content with the smaller scope. Consider rolling forward to a later paper with expanded cohorts instead. |
| Nothing survives at any layer | **Do not publish as standalone biology.** Write up as methods + negative constraint paper at a methods venue, or fold into a future collaboration. |

The decision is *data-derived*, not pre-decided. The framework runs to completion; the decision artifact records which branch the data put us on.

---

## 6. File and module layout (what to create)

### 6.1 New source modules

```
src/networks/
    contrast_vectors.py          # Core: A^a_GC, A^a_FLT, beta, cos, rho, R; bootstrap; permutation
    stability_test.py            # Guardrail A: bootstrap stability of A^a_GC at each resolution
    cross_osdr_projection.py     # Phase 6: external study F_s vectors, cosine alignment, meta-analysis
    external_aging_axis.py       # Phase 7: TMS kidney aging axis projection

src/statistics/
    bootstrap_decomposition.py   # Bootstrap pipeline runner used by Phase 3-7

src/preprocessing/
    multi_study_harmonization.R  # Phase 1: per-study residualization + harmonized output
    common_gene_universe.py      # Phase 0/1: Ensembl intersection across studies

scripts/
    run_contrast_vector_framework.py   # Orchestrator for Phases 0-8
    audit_stability_gates.py            # Convenience: check that Phase 3 gates have been run
```

### 6.2 New configs

```
config/contrast_vector_framework.yaml   # Pre-registered thresholds from §4
                                        # Bootstrap iteration counts (B=2000, K_perm=5000)
                                        # Resolution hierarchy + stability gates
                                        # External cohort manifests (paths, sample IDs)

config/aging_reference.yaml             # TMS kidney aging axis: path, columns, expected sample counts
```

### 6.3 Output directory layout (per run)

```
data/results/<run_id>/contrast_vectors/
    agc_stability_report.tsv             # Phase 3 gate outputs
    agc_stability_decision.json           # Pass/fail per resolution × arm
    ISS-T/
        edge/                             # If gate passed
            bootstrap.tsv
            permutation.tsv
            point_estimates.json
        gene/
            ...
        module/
            ...
        pathway/
            ...
        redirected_module_enrichment.tsv
        redirected_pathway_enrichment.tsv
        redirected_top_genes.tsv
    LAR/
        (same structure as ISS-T)
    cross_osdr_alignment/
        OSD-513/
            F_vs_RRRM2_ISS-T.tsv          # cosine, bootstrap CI
            F_vs_RRRM2_LAR.tsv
            F_vs_DeltaAgeRewire.tsv
            F_vs_R_redirected.tsv
        OSD-253/
            (8 scenarios: strain × duration × control_scenario)
            scenario_matrix.tsv
        meta_analysis_summary.tsv
    external_aging_axis/
        tms_projection_ISS-T.tsv
        tms_projection_LAR.tsv
    manuscript_decision.json
```

### 6.4 New tests

```
tests/test_stability_test.py            # Synthetic data: stable vs unstable vectors trigger gate correctly
tests/test_contrast_decomposition.py    # Known geometry test: hand-built A_GC and A_FLT give expected β
tests/test_external_aging_axis.py       # TMS axis loads cleanly; projection respects gene-universe intersection
tests/test_cross_osdr_alignment.py      # External study cosine computation; bootstrap CI behaves
tests/test_decision_tree.py             # Manuscript decision artifact maps outcomes to branches correctly
```

---

## 7. What this framework explicitly is NOT

- Not a WGCNA-only paper. WGCNA modules are used for biological interpretation of $\mathbf{R}_a$ (Phase 5), not as the primary inferential object.
- Not a DE-only paper. DE is used for the silent-shifter recurrence (high network change + low DE) but not as the headline.
- Not a generic "we did everything and reported" paper. The framework is structured to commit to a specific scientific question and answer it at the resolution the data supports.
- Not a forced-discovery paper. The decision tree (§5 Phase 8) explicitly includes "no positive result" branches.
- Not a node2vec/Procrustes-centric paper. Node2vec embeddings are exploratory only; the contrast-vector decomposition is the primary geometry.
- Not a paper that pools raw expression across missions. All cross-OSDR work happens at contrast-statistic level after within-study estimation.
- Not silent-shifter-as-positive-category. Silent shifters are now defined as *recurring rank patterns across studies*, tested by rank aggregation, not as a single-study FDR-passing gene list.

---

## 8. Acceptance criteria for "MV is shipped"

The minimum viable analysis is shipped when all of the following exist:

1. `data/processed/multi_study_metadata.tsv` with all four MV studies.
2. `data/processed/<study>/residuals.tsv.gz` for OSD-771, OSD-513, OSD-253 (with C57 and C3H subsets).
3. `data/processed/networks/control_pooled_backbone/skeleton_edges.tsv` with frozen manifest.
4. `data/results/<run>/contrast_vectors/agc_stability_report.tsv` AND `agc_stability_decision.json`.
5. Within-RRRM-2 decomposition at every resolution that passed the gate, with bootstrap CIs and permutation p-values.
6. Redirected-component module/pathway enrichment per arm.
7. Cross-OSDR alignment with OSD-513 and the 8 OSD-253 scenarios.
8. External-axis (TMS) projection per arm.
9. `manuscript_decision.json` with the branch assignment.
10. Unit tests (the five listed in §6.4) pass.
11. `scripts/audit_fix_status.py` reports `present=14`.
12. A short results memo (1-2 pages) summarizing per-resolution decomposition outcomes for the manuscript author.

---

## 9. Engineering style and conventions for agents

- **Test before claim.** Every new function gets a unit test. Synthetic data tests for the geometric functions (e.g., `compute_beta_cos_rho(A_FLT, A_GC)` with known inputs) are non-negotiable.
- **Bootstrap iteration counts are config-driven.** Never hard-code $B$ or $K$ inside a function. Read from `config/contrast_vector_framework.yaml`.
- **Manifest everything.** Every run output writes a manifest with: git commit, config hash, input dataset hashes, output dataset hashes, parameter values, completion timestamp.
- **No silent fallbacks.** If `A_GC` fails the stability gate at gene level, do NOT silently fall back to module level. Halt, log the failure, and require explicit user re-run with `--bypass-stability` (which should be reserved for sensitivity analyses, not headline results).
- **Backbone is frozen.** Once the backbone is locked, do not modify it. If a new dataset is added, build a fresh backbone in a new directory.
- **Two-stage decision logging.** Phase 3 writes a decision JSON. Phase 4 reads it before running. Phase 8 reads all upstream decision JSONs before writing the manuscript decision.
- **Reproducibility is not optional.** Every run gets a unique run ID with timestamp. Every result is traceable to a commit + a config + a manifest.

---

## 10. Open questions for the human owner (Ibrahim) before MV starts

These are decisions only the human can make. Document the answers in `docs/owner_decisions.md` before any agent runs Phase 0.

1. Do you want the within-RRRM-2 control aging vector estimate to use ISS-T HGC + VIV pooled, or only HGC? (Pooling HGC + VIV gives larger N but introduces a potential housing-type confound.)
2. Do you want the OSD-253 sensitivity to lead with the 2021 white-LED rerun controls (clean light, dirty read-length) or the original GC (clean read-length, dirty light)? Both should be run; which is "primary" for the manuscript?
3. Is `Basal` included in the GC pool for the backbone, or only `HGC`+`VGC`? (Basal mice were sacrificed 2 days before launch; including them broadens the control reference but introduces a temporal offset.)
4. Tabula Muris Senis kidney aging axis: which age contrast specifically? (Young = 1–3 month vs Old = 18–24 month is standard; consider also Young = 3 month vs Mid = 18 month if that better matches RRRM-2's 16- and 34-week mice.)
5. Is the manuscript target still npj Microgravity? Or does the new framework push toward Communications Biology / NAR Genomics & Bioinformatics? (Affects writing tone but not analysis.)
6. Hard scope-creep budget: if MV ships clean, does extension to OSD-102 / OSD-163 happen in the same paper, or a follow-up? (Lock now; revisit only with explicit owner authorization.)
7. Wet-lab follow-up (CRISPRi grey60 hubs, TLR4 organoid work, etc.) — is this on the table for the current paper, or strictly future work? (Affects the discussion section's tone.)

---

## 11. Where this document lives and how to update it

This document is `docs/CROSS_OSDR_NETWORK_CONTRAST_PLAN.md`. It is the canonical reference for the analysis. Any agent or collaborator working on Phases 0–8 must read it first.

If the plan changes mid-stream:
1. Open a PR that modifies this document.
2. Note the change in `docs/plan_changelog.md` with date and rationale.
3. Do not modify in-flight phase outputs; spin a new run ID.

If the plan does not change but the implementation reveals an issue (e.g., a stability gate is too tight or too loose), the right move is to add a sensitivity-analysis branch to the affected phase, NOT to silently adjust the gate. Document the sensitivity branch in `docs/plan_changelog.md`.

---

## 12. Bottom line

The original methodology is not being discarded. It is being reframed from *"RRRM-2 alone discovers age-dependent network rewiring genes"* to *"RRRM-2 defines the only exact age-by-arm flight network contrast, and a cross-OSDR framework tests whether that network direction recurs across independent kidney spaceflight cohorts."*

That is closer to the original ambition than any of the intermediate pivots, and it avoids pretending that $n=5$ per factorial cell can support strong individual edge discovery by itself.

The five guardrails (stability gate, resolution hierarchy, shrinkage weighting, full-pipeline bootstrap, external aging axis) are what convert this from "another small-cohort fishing expedition" into a defensible quantitative test of a specific scientific question.

Minimum viable scope ships in 3–4 weeks of focused work assuming the existing repo machinery handles most of the inner computation. If MV produces a positive result, the paper writes itself. If MV produces a negative result, the framework itself is the contribution.

Do not deviate from this plan without updating it first.
