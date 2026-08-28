# Methods paper: analysis plan

**Date:** 2026-07-29
**Status:** draft plan, not yet frozen. Nothing here should be run before §7.

---

## 0. Novelty check — do this before anything else

Two quick searches show the space is not empty, and the plan is positioned
accordingly.

**Already published, do not claim as new:**

- Gene-set method power comparisons (GSEA vs global test vs Wald-type), including
  correlation-aware simulations.
- RNA-seq power/sample-size analysis with mean–variance dependence.
- Sample size effects on gene-set analysis reproducibility and specificity.
- WGCNA module stability via resampling — `rWGCNA`, SABRE, and module-validation
  benchmarks all exist.

**Not covered, and where this paper lives:**

1. **Power as a function of prior accuracy.** Existing comparisons treat the gene
   set as given. Nobody parameterises the overlap between your external \(w\) and
   the truth and asks where the crossover is. That overlap is the practically
   decisive quantity and it is never reported.
2. **The derivation-source axis.** Same-\(X\)-derived \(w\) versus external \(w\)
   versus transfer-derived \(w\), compared at fixed everything else. The gap
   between the first and third is the empirical cost of Theorem 3, and it has not
   been measured.
3. **Identifiability/abstention as a unifying frame** across layer carriage,
   network derivation, and annotation attribution.
4. **Three worked negative cases on a heavily reused public cohort.**

Step 0 is a real literature review, not two searches. Budget a week. If (1) and
(2) turn out to be covered, the paper reduces to the abstention framework plus
the empirical cases — still viable, but decide that consciously rather than
discovering it in review.

---

## 1. Thesis

> In small-\(n\), bundled-exposure, heavily reused omics cohorts, the binding
> constraint is not power but whether the estimand is identifiable from the
> available structure. A reanalysis framework should screen claim classes for
> identifiability before testing, and abstain when the screen fails.

Each demonstration tests a different failure mode: **carriage** (does an effect
at one layer transfer to the next), **derivation** (can a statistic built from
\(X\) carry more than \(X\)), **attribution** (can a whole-tissue measurement be
assigned to a compartment).

If that sentence cannot be written cleanly at drafting time, do not merge
LayerScore and the network note — ship them separately.

---

## 2. S1 — The simulation (the main new work)

The only part requiring substantial new computation, and the paper's positive
half. Without it the paper is purely negative.

### 2.1 Generative model — calibrated, not toy

- Take the RRRM-2 ground-control VST matrix; keep \(m = 5{,}000\) most variable
  genes.
- Fit \(X = \mu + \Lambda f + \varepsilon\) by PCA on GC samples only,
  \(K \in \{10,20\}\) factors, gene-specific residual variance taken from the
  real residuals. This inherits realistic correlation structure and mean–variance
  trend.
- Define "true modules" as gene blocks with high \(|\Lambda_{gk}|\) on factor
  \(k\). These are the objects an external prior is trying to match.

### 2.2 Alternative classes

Injected as a mean shift in the flight arm. **Total signal energy
\(\lVert\delta\rVert_2\) is held constant across A1–A3** so the comparison is
about the distribution of signal, not its amount. This is the design decision
that makes the result interpretable.

| Class | Structure | Purpose |
|---|---|---|
| A0 | no shift | type I error calibration |
| A1 | 20 genes, large \(\delta\) | sparse–strong; DE should win |
| A2 | 200 genes of one true module, small \(\delta\) | dense–weak–aligned; a good \(w\) should win |
| A3 | 200 random genes, small \(\delta\) | dense–weak–**diffuse**; the control that stops "module scores always win" |
| A4 | A1 + A2 | mixed |
| A5 | within-module correlation 0.3 → 0.6, **no mean shift** | second-order only; Proposition 1 made empirical |

### 2.3 Methods compared

| | \(w\) source | Tests |
|---|---|---|
| M1 | — (gene-wise limma + BH) | baseline |
| M2 | oracle module membership | upper bound |
| M3 | external \(w\) at overlap \(\rho \in \{1, .75, .5, .25, 0\}\) with truth | **prior-accuracy axis** |
| M4 | WGCNA eigengene fit on the *same* samples | the capped case |
| M5 | WGCNA eigengene fit on an *independent* draw | transfer case |
| M6 | differential correlation test | for A5 |
| M7 | node2vec → Procrustes → cosine rewiring | optional; the pipeline the theory targets |

\(n\) per group \(\in \{3, 5, 10, 20, 50\}\); 1,000 replicates per cell.

M7 is expensive and optional. Include it if feasible — it is the only way to
demonstrate the §5 circularity and Lemma 1 costs empirically rather than by
argument. Drop it first if time is short.

### 2.4 Outcomes

- Type I error at nominal 5% under A0.
- Power at 5% FWER (set-level) and FDR (gene-level) under A1–A5.
- **The phase boundary:** for each \((n, \text{class})\), the minimum prior
  overlap \(\rho^\star\) at which M3 beats M1.

### 2.5 Predictions — write these down before running

| | Prediction |
|---|---|
| P1 | M1 wins under A1 at every \(n\). |
| P2 | M2 and M3(\(\rho{=}1\)) win under A2, with the margin growing as \(n\) falls. |
| P3 | \(\rho^\star\) is high (0.5–0.7): imprecise priors do not help. |
| P4 | M4 is anticonservative under A0 and/or shows apparent power that does not transfer — Theorem 3's practical cost. |
| P5 | M5 is valid but attenuated relative to M2. |
| P6 | Under A5, **everything** has near-nominal power at \(n=5\). |
| P7 | Under A3, module scores do **not** beat DE. |

P4, P6 and P7 are the quotable results. P7 in particular is what keeps the paper
honest — a simulation in which the favoured method always wins is not evidence.

### 2.6 Deliverable

A phase diagram over \((n, \text{sparsity}, \rho)\) with the decision rule
stated in one line, e.g. *use a module score when the alternative is dense and
weak and your prior overlaps truth by more than \(\rho^\star\); otherwise use
gene-wise testing; never use a \(w\) fitted on the samples you test.*

---

## 3. S2 — Generalise the Grey60 non-recovery (cheap, high value)

**Question:** is 0/27 specific to Grey60, or generic to WGCNA modules at
\(n \approx 20\)?

**Method:** for *every* module in the original RRRM-2 run — not just grey60 —
freeze membership, re-run flight-blind on BSL+VIV across the same 27
specifications, compute per-module Jaccard recovery.

**Outcome:** the distribution of recovery across modules. If the median Jaccard
is low, the claim becomes *"WGCNA modules at this sample size do not reconstruct
under a held-out sample split"* — general and citable, rather than an anecdote
about one module.

**Second output:** how often the colour-matched module differs from the
overlap-matched module across specs. This turns the
\(Z_{\text{summary}} \approx 20.6\) mislabelling into a quantified general hazard
rather than a one-off error.

Machinery already exists. Effort: low.

---

## 4. S3 — RNA-preparation reliability (do this first)

The best new empirical idea in the plan, and the cheapest way to de-risk
LayerScore's premise.

UPX 3′ tag, polyA mRNA and total RNA on the **same 20 animals** is a natural
reliability experiment. Nobody planned it that way; it exists.

**Steps:**

1. Estimate reliability of the flight effect at the RNA layer: cross-preparation
   correlation of per-gene effects, and ICC across the three preparations.
2. Do the same for the protein layer wherever a split is available (split-half
   across animals, or across plexes).
3. Compute the attenuation implied: an observed cross-layer correlation is
   bounded by \(\sqrt{r_R \, r_P}\) times the true one.
4. Ask the decisive question: **is the observed RNA–protein \(\rho = -0.034\)
   consistent with a true correlation of 0.3 given the measured reliabilities?**

If yes, "the RNA program does not propagate to protein" is demonstrated to be
unidentifiable — with data, not assertion. That is the strongest single result
available from this material, and it converts the preparation sign-flips from an
embarrassment into a measurement.

**De-risking value:** if reliabilities come back *high*, attenuation cannot
explain the cross-layer null, and LayerScore's central premise weakens. Better to
learn that in week two than month five. **Run S3 before committing to the merge.**

---

## 5. S4 and S5 — supporting analyses

**S4 — within-block reporter-position diagnostic.** As specified in
`ANALYSIS_TRIAGE_AND_NETWORK_REFRAME_2026-07-29.md` §1. Five exchangeable animals
per block means a within-block positional trend is a pure tag effect. Six
estimates (3 blocks × 2 plexes) per modality. Primary home is the data note; the
methods paper cites it as the attribution-failure exhibit.

**S5 — transfer template.** Mostly writing. Formalise structure-from-A evaluated
in held-out-B against panels matched on size, expression and detection.
`run_20260701_sf_classifier/` supplies the leave-one-cohort-out machinery; Grey60
(terminal meta \(g = 0.163\), CI −0.658 to 0.984) is the worked negative
instance.

---

## 6. Sequencing and effort

| Order | Analysis | Effort | Why here |
|---|---|---|---|
| 1 | Step 0 literature review | 1 week | Decides whether §2 is novel |
| 2 | S3 reliability | 1–2 weeks | De-risks the merge; highest payoff per unit work |
| 3 | S2 module recovery | 1 week | Cheap; generalises a result already in hand |
| 4 | S4 reporter position | 3–5 days | Also unblocks the data note |
| 5 | **Freeze S1 design (§7)** | 1 day | — |
| 6 | S1 simulation | 6–10 weeks | The main work |
| 7 | S5 + drafting | 4–6 weeks | — |

Realistic solo timeline: data note out in 4–6 weeks running in parallel; methods
preprint 4–6 months. Do not compress S1 by cutting cells from the grid — cut M7
instead.

---

## 7. Freezing discipline

The credibility of this paper depends entirely on the analyses being frozen
before they run. This is the one thing that distinguishes it from the previous
twelve revisions, and a simulation is the easiest thing in science to tune until
it looks good.

Before S1 runs, commit and tag: generative model and calibration source,
alternative classes and their energy normalisation, the method list and \(\rho\)
grid, the \(n\) grid, replicate count, outcome definitions, and the P1–P7
predictions verbatim. Record the commit hash in the manuscript. Any post-hoc
addition goes in a clearly labelled exploratory section.

Also: commit the v13 run directories and remove the `latex_paper/*` gitignore
exclusion. Every "frozen before testing" statement in this project is currently
unverifiable from repository history, and that is fixable in an afternoon.

---

## 8. Kill criteria — decide now, not later

| Finding | Consequence |
|---|---|
| Step 0 shows (1) and (2) are already published | Drop S1's novelty claim; paper reduces to framework + three cases. Still viable, but reframe deliberately. |
| S3 shows high reliability at both layers | LayerScore's premise weakens. Do not paper over it — report it, and the merged paper's carriage section becomes a negative result about the *method*, not the data. |
| S1 shows the useful regime is empty | Acceptable, and arguably a stronger result: *the regime where these methods help is empirically rare.* Agree to this outcome in advance. |
| S1 shows module scores win broadly, including A3 | Suspect a bug. A method that wins under diffuse signal is detecting something other than structure. Debug before believing. |
| S2 shows most modules *do* recover flight-blind | Grey60's failure becomes specific rather than generic. Weaker, still reportable, and an honest correction to §3. |

The third row is the important one. Agreeing in advance that a negative phase
diagram is a publishable outcome is what prevents S1 from being tuned.
