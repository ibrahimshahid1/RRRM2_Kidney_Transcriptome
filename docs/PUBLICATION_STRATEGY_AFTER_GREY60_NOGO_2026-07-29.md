# Publication strategy after the Grey60 NO-GO

**Date:** 2026-07-29
**Trigger:** Grey60 locked go/no-go returned NO-GO. Gates B (flight-blind
recovery, 0/27) and E (terminal recurrence, g = 0.163, CI −0.658 to 0.984) are
independently fatal; Gate D failed narrowly and is cautionary.

---

## 1. The Grey60 audit was the right call, executed well

Worth stating plainly, because it is easy to read a NO-GO as another failure.
The audit froze 48 genes, re-ran WGCNA on BSL+VIV only so that FLT/GC could not
inform the reference network, swept 27 specifications, used modified
Hartung–Knapp for the terminal meta, caught a technical-replicate suffix parsing
bug before issuing a decision, accepted an independent code audit that made the
pooled result *weaker*, and kept it.

It also disposed of a claim the earlier drafts leaned on: the
\(Z_{\text{summary}} \approx 20.6\) preservation statistic belongs to a broad
439-gene GC-only module, not to compact Grey60 (\(Z \approx 7.8\)). WGCNA colour
labels are run-specific and were being treated as module identity. That
correction alone justifies the exercise.

This is a higher standard of evidence than most published spaceflight-omics
reanalysis meets. That fact is load-bearing for §4.

---

## 2. Why OSD datasets behave like this

Not because GeneLab is careless — publishing ISA-formatted, reusable data is
right. The problem is that reuse practice treats these like ordinary omics
cohorts, and structurally they are not.

**a. \(n\) is small and irreducible.** Five per group, set by payload mass, cage
capacity and crew time. No analysis fixes this.

**b. The exposure is a bundle, not a variable.** "Spaceflight" is microgravity +
radiation + elevated CO₂ + housing + light cycle + launch vibration + diet +
handling + carcass timing. Ground control matches some of it. Different missions
bundle it differently, so cross-cohort "replication" is replication of a
*mission*, not of a treatment. Gate D found exactly this: baseline-versus-
vivarium \(g = 0.563\) — the score responds to housing/handling context.

**c. Batch structure *is* the design.** With five per group you cannot randomise
technical layout without splitting groups or wasting capacity, so people do not.
OSD-462's reporter-channel aliasing is this disease in proteomic form; condition
becomes inseparable from batch permanently, and no reanalysis undoes it.

**d. Preservation and collection confounds are unrandomisable.** Terminal versus
live-return, on-orbit versus post-return euthanasia, freezing versus fresh
dissection, hours-to-days between death and preservation. These correlate with
condition rather than varying at random.

**e. Pipeline variance is comparable to the effect.** Panel D of the Grey60
figure is the cleanest demonstration in this repository: the *same 20 animals*
give \(g = -0.09\) (UPX 3′ tag), \(+0.98\) (polyA mRNA, p = 0.024), and
\(-1.16\) (total RNA, p = 0.015).

One nuance worth preserving rather than flattening: those three preparations
measure genuinely different RNA populations — total RNA carries non-polyadenylated
and intronic/nascent signal, 3′ tag has different length bias. So the discordance
is not purely technical noise. It means something sharper and more useful:
**"the flight RNA effect" is not a well-defined quantity until the preparation is
specified.** That is an estimand problem, not a data-quality complaint, and it is
the same lesson as everything else in this project.

**f. Effective multiplicity across the field is uncontrolled.** The same handful
of matrices are reanalysed by many groups looking for the same shape of result.
Apparent consistency in the literature is partly survivorship.

**The honest summary:** OSD data are well suited to description and hypothesis
generation, and structurally poorly suited to effect estimation, replication and
mechanistic attribution. Twelve revisions asked them for the second kind of
thing. That was a category error made in good faith, and it is the reason the
methods direction is not a consolation prize — it is the correct read of what
this material can support.

---

## 3. Do not run LayerScore and the network note as two papers

They look like different projects. They are the same paper.

| | LayerScore | Network bound |
|---|---|---|
| Question | Can an RNA-only effect be called non-propagating? | Can a derived network statistic out-inform \(X\)? |
| Limiting quantity | Reliability / attenuation of cross-layer estimates | Data-processing bound; \(O(1)\) per-edge variance at \(n=5\) |
| Correct output | Abstention (`indeterminate`) | Abstention (no gain unless \(Z \not\to X\)) |
| Worked example | OSD-462 | OSD-771 + 4 cohorts |

Both say: *the binding constraint is not power, it is whether the estimand is
identifiable from the available structure — and the correct output when it is
not is a principled refusal, not a weaker claim.* Two thin notes making the same
argument from one repository compete with each other and neither reaches
critical mass.

### The merged paper

**Thesis.** In small-\(n\), bundled-exposure, heavily-reused public omics
cohorts, claim classes should be screened for identifiability before testing,
and a framework should be able to abstain.

**Structure:**

1. **Framework.** LayerScore's decision rule and abstention class provide the
   actionable core — this is what a reader can adopt.
2. **Theory.** The information-bound propositions provide rigour: DPI chain,
   \(O(1)\) differential-correlation variance, linear-functional collapse
   (Prop 4), BH permutation floor, rotation non-identifiability, Procrustes
   circularity. Most methods papers in this space have no theorems.
3. **Demonstration I — derivation limits.** Grey60: 0/27 flight-blind
   recoveries, Jaccard 0.022, projected \(p = 0.831\); the \(Z_{\text{summary}}\)
   mislabelling; terminal meta \(g = 0.163\). Theorem 3 predicts it, the audit
   measures it.
4. **Demonstration II — layer limits.** OSD-462 RNA/protein/phospho, plus the
   `run_20260701_sf_classifier` leave-one-cohort-out result where network
   features failed to beat expression under leakage-safe evaluation.
5. **Demonstration III — attribution limits.** The v13 subtype test: frozen
   coverage gate, exact 63,504-assignment permutation, intensity-matched null,
   machine-readable tier `neither`. A framework refusing to certify its authors'
   preferred hypothesis is the paper's strongest single exhibit.
6. **Positive half.** The simulation locating where a fixed external \(w\) beats
   gene-wise testing at \(n = 5\text{–}10\). Without this the paper is only
   negative; with it, it tells people when to use these methods.
7. **Transfer template.** The one design that escapes the bound — structure from
   cohort A evaluated in held-out B against panels matched on size, expression
   and detection. Grey60 is the worked negative instance.

Section 6 is the only part requiring substantial new work, and it is
simulation, not new data.

---

## 4. Sequencing

**First — the OSD-462 data note (weeks, mostly writing).** Different genre and
different audience from the methods paper, so they do not compete, and it gives
the methods paper something to cite for provenance. Contents: TMTpro-versus-
legacy-iTRAQ metadata, the condition-by-reporter-tag map, the intensity gradient
(Spearman −0.63 official-scaled versus −0.14 summed S/N), the Y65/Y68 residue
indexing and co-modified peptide problem — **and now the same-animal RNA
preparation sign flips**, which belong here as a property of the dataset that
reusers will hit silently. Four independent ways OSD-462 misleads a downstream
user is a strong, useful note.

Add the within-block reporter-position diagnostic first; it is the one
measurement that turns the tag section from a design critique into a quantity.

**Second — the merged methods paper (months).** Everything above.

**Not now — anything biological from this material.** No v14, no podocyte, no
Grey60 rescue, no OSD-253 duration follow-up (−0.007 at 25 d versus +1.283 at
75 d is interesting and underpowered; note it, do not chase it).

---

## 5. What LayerScore loses and gains

Loses: standalone identity, and the framing where OSD-462 is a demo appended to
a general method.

Gains: theorems, three worked demonstrations instead of one, and a thesis with a
point of view rather than a proposed metric. A methods paper whose central claim
is "here is when to abstain, here is why, and here are three cases where we did"
is considerably harder to reject than one proposing another score.

The γ_P regression-dilution correction and the
`rna_only_supported` / `indeterminate` distinction survive intact and become the
framework's operational core.
