# Analysis triage, and what network methodology can legitimately ask

**Date:** 2026-07-29
**Companion to:** `V13_INDEPENDENT_REVIEW_AND_PATH_FORWARD_2026-07-29.md`

Two questions: (1) is the proposed Tier 1–4 compartment battery worth running,
and (2) given the information bound, how does one form real questions with
network methodology?

---

## Part 1 — Triage of the proposed battery

### Tier 1: do not run

Every Tier 1 item — podocyte, fibroblast, endothelial, TAL, principal,
intercalated, and the full frozen compartment comparison — terminates at the
same sentence: *annotation enrichment inside a tag-aliased contrast*. The
battery is well designed. That is the problem: you could execute all of it
cleanly, get an unambiguous winner, and still not be able to write a biological
claim, because the controls address selection and observability while the
blocking failure is design (flight ≡ channels 129N–131N in both plexes) and
attribution (parent protein ≠ cell of origin).

The complete frozen compartment comparison has in any case already been run
(`primary_compartment_enrichment.tsv`) and re-run under intensity matching
(`run_20260729_v13_intensity_confound_audit/`). Repeating it with more
sensitivities does not change what it is.

### One item is degenerate — check before spending a day on it

**Direction-reversal control.** For the primary gene statistic
\(G_g=\operatorname{median}_s(-\hat\beta_s)\), the flight-increased score is
\(\operatorname{median}_s(\hat\beta_s)=-G_g\) exactly. The balanced within-plex
permutation set is closed under global label flip, so the null is symmetric
about zero and \(Z'_g=-Z_g\) exactly. The competitive set statistic is linear in
\(Z\), so the reversed statistic is exactly the negative of the forward one:
podocyte would be −1.11, fibroblast −0.71. The control cannot fail and tests
nothing.

It *is* informative for `one_sided_maxmean`, which treats the two tails
asymmetrically. If you run a direction-reversal control, run it there only.

### Four items worth running — but as evidence for the data-quality note

These do not license a biological claim. They characterise the artefact, which
is what product A needs.

**1. Reporter-position diagnostics — highest-value analysis remaining.**

There is a clean, unconfounded version of this that the plan does not state.
Within each treatment block, the five reporter channels hold five *biologically
exchangeable* animals of the same group. So a systematic trend across channel
position **within a block** is a pure tag-position effect with no biological
confound at all.

Concretely:

- For each of the three blocks (baseline 126–128C, flight 129N–131N, ground
  131C–133C), regress per-site log signal on within-block channel index.
- Do this separately for both plexes and both modalities (protein and
  phosphosite), so you have six independent estimates of the same nuisance
  slope.
- Test whether the within-block slopes agree in sign and magnitude.
- If they do, extrapolate the fitted positional trend across the block boundary
  and compute how much of the observed flight-minus-ground difference the
  positional effect alone would produce.

That last number is what makes the note citable rather than a complaint: it
converts "condition is aliased with tag position" from a design observation into
a measured effect with a magnitude and an uncertainty. Either outcome is worth
publishing — a detectable slope bounds the flight effect; a null slope is an
honest statement that the aliasing has no measurable magnitude at this depth,
which is genuinely useful to anyone reusing OSD-462.

Pair it with the intensity gradient already measured
(Spearman −0.63 official-scaled vs −0.14 summed S/N): both are properties of the
processed workbook that downstream users will hit silently.

**2. Protein-class control.** Build a cytoskeleton/adhesion/contractile set and
run it through the identical pipeline alongside the compartment sets. The
podocyte set (Synpo, Magi2, Iqgap2, Myo1e, Parva, Nck2, Pak1, Ablim3, Tmod3,
Mtss1) and the fibroblast set (Acta2, Myl9, Tagln, Cnn1, Myh11, Tpm2) are
cytoskeleton-dominated; the proximal-tubule set, at the opposite end of the
ranking, is metabolic enzymes and solute carriers. If a protein-class set beats
every compartment set, the compartment axis is a protein-class axis and every
compartment framing dies in one line. Cheap, and decisive either way.

**3. Restricted versus broad marker stratification.** One pass over the existing
sets. Feeds the note directly.

**4. Leading-edge coherence.** Already largely visible in
`leading_parent_gene_matrix.tsv`. Write it down rather than re-derive it.

### Skip

- NCC regulatory versus non-regulatory features — moot; the residue identities
  were wrong and no isolated canonical phosphoform qualified.
- Phosphoform-quality ladder, alternative parent-gene statistics,
  leave-one-gene-out, centred versus uncentred — all already in the exact run.
- Baseline as a biological validation group — already correctly rejected. Its
  only use is as a third channel block in diagnostic (1).

---

## Part 2 — Network methodology after the information bound

### The bound is right; its own §10 is undersold

The note is correct and I would not change the theory. Theorem 3 is exactly
right, Proposition 2 is a sharp and useful observation, and Proposition 3 is a
practical calculation that most people never do.

But it files the important thing under "caveat." Theorem 3 is a statement about
*information about \(Y\)*. Almost no analysis is limited by information — it is
limited by **power at fixed \(n\)**, which is a bias–variance question, not an
information question. §10 is where all the legitimate action is, not a footnote
to §9.

### Proposition 4 is the reframe

Proposition 4 says that a contrast on a linear network summary *is* a
differential-expression test on the composite feature \(Xw\). Read positively
rather than negatively, this means:

> The entire content of "network methodology" is the choice of \(w\).

So the well-posed question is never "does the network know something \(X\)
doesn't." It is:

> Is there a weight vector \(w\) that I could not have written down without the
> network, and does it concentrate power on a class of alternatives I actually
> care about at \(n=5\)–\(10\)?

That is answerable, falsifiable, and it is a methods question rather than a
biology question.

### Three regimes for \(w\)

| Regime | Capped by Thm 3? | Examples | Verdict |
|---|---|---|---|
| \(w=f(X)\) — network estimated from the same matrix | **Yes**, and worse in practice | LIONESS, node2vec on the co-expression graph, direct differential correlation, rewiring cosines | Dead. Prop 1, Lemma 1 and the §5 circularity explain precisely why v1–v2 failed. Do not revisit. |
| \(w=f(Z)\), \(Z\) external | **No** | Pathway/regulon graphs, kinase–substrate priors, cell-type atlases, interactome propagation | Legitimate. This is what v11–v13 have actually been doing without calling it network methodology. |
| \(w=f(X_{\text{other}})\) — structure from a *different* matrix | **No**, w.r.t. the \(X\) being tested | Grey60 defined in ground-control animals and tested in flight; module defined in cohort A, tested in cohort B | Legitimate, and the one positive result in this repository has exactly this shape. |

The third row is the sleeper. A module estimated in one cohort is not a function
of the cohort you evaluate it in, so the bound does not apply to the transfer
test. `run_20260701_sf_classifier/` already contains leave-one-cohort-out AUC by
panel, random-panel nulls, and a negative control — the machinery for this exists
and was built for a different purpose.

### Four question templates that survive the theorem

1. **Transfer.** Is structure estimated in condition/cohort A predictive in
   held-out B, beating random panels matched on size, expression and detection?
   Non-circular by construction; the estimand is stated cleanly; a negative
   answer is a result. Grey60 preservation (\(Z_{\text{summary}}=20.6\)) is a
   weak version of this — the strong version asks whether a GC-defined module
   *predicts the direction* of the flight response in a cohort it was never
   fitted on.

2. **Power / estimand (simulation).** For which classes of alternative does a
   fixed external \(w\) beat gene-wise testing at \(n=5\)–\(10\)? Sparse
   large-effect alternatives favour DE; dense weak low-rank alternatives favour
   a good \(w\); a wrong \(w\) is worse than DE everywhere. Proposition 1 is the
   negative half of this analysis. The positive half — where the crossover
   actually is — is unwritten and is a real contribution.

3. **Organisation / constraint.** Given an external graph \(G\), is the observed
   DE vector *smoother on \(G\)* than a degree- and expression-matched null?
   This does not ask the network to discover anything; it asks whether the
   response is localised in the interactome or diffuse. It uses external \(Z\),
   so the bound does not bite, and a diffuse answer is publishable and
   interesting.

4. **Design.** Which nodes or segments would most reduce uncertainty about a
   module state if measured? Honest because it makes no claim about existing
   data, and it is the correct foundation for the "next experiment" section of
   any of these manuscripts.

### The concrete recommendation

`network_information_bound.pdf` is about 90 % of a publishable methods note, and
this repository already contains its empirical companion: a real dataset in
which the network features were run properly, evaluated leakage-safely, and
failed to beat an expression baseline. Theorem 3 predicts that failure; the
classifier run measures it.

A note consisting of (i) the bound and its chain of propositions, (ii) the
measured failure with leave-one-cohort-out evaluation and matched random-panel
nulls, and (iii) a small simulation locating the crossover point from template 2,
is a genuine contribution. Negative methods results accompanied by a theorem
explaining them are rare, and this one closes the loop back to the question the
repository was started to answer.

It is arguably a stronger standalone than the LayerScore note, and the two
merge naturally: both are about what a derived statistic can and cannot carry.

### Have these already been tried?

Partly, and it is worth being precise about which part. The \(w=f(X)\) versions
were tried in the original draft and v2 and failed. What has *not* been done is
the reframing: treating that failure as the result, rather than as a setback on
the way to a biology paper. The theory that explains the failure was written in
June 2026, after the fact, and has never been connected back to the empirical
run that demonstrates it.
