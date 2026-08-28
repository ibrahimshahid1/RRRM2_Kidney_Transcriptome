# Direction after v13 and the clinical-axes no-go

**Date:** 2026-08-11
**Status:** historical pre-audit recommendation; K1–K4 are now completed
**Trigger:** `CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md` returned no-go on
all four frozen axes, leaving one unexplained podocyte-associated lead.
**Constraint:** preprint by November 2026 – January 2027; ~3 months of work;
NASA PI available for tissue access, suggestions, and publishing help.

> **Post-audit note.** The K1–K4 predictions in this document are superseded by
> `PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md`. In brief: Siew's
> PODXL/NID1 consensus is not urine-only and includes mouse kidney RNA; no
> forward-citing recurrent podocyte-transcript claim was found; forcing
> `Podxl`/`Nid1` did not weaken the program; and strict same-tier matching
> retained support across the five observed missions. Read the completed audit
> before using the directional-conflict language below.

---

## 1. Decision

Write **one paper**, and make it a positive directional claim rather than a
negative audit:

> Podocyte-associated transcripts are consistently *higher*, not lower, in
> mouse kidney after terminal spaceflight: a mission-level reanalysis of five
> independent cohorts.

The four-axis null becomes the specificity control inside that paper, not the
headline. The validation experiment becomes the closing section and the basis
for a collaboration request, not a precondition for submission.

This is achievable by January on work that is already done.

---

## 2. The reframe: the no-go was a category error

The 2026-08-11 report concluded no-go **as an injury-biomarker paper**. That
conclusion is correct and should stand. It was then read as "no paper," which
does not follow.

Every product decision in this repository since May has been scored against one
rubric: *is this a biology discovery paper?* The answer has been no thirteen
times. But the rubric was never the only one available. The repository has for
a year been producing a different object — a mission-level estimate with an
interval attached — and that object is publishable on its own terms.

The retrospective already states the operative principle and then does not
apply it:

> A null is publishable when it is narrow. "We failed to detect" is not a
> finding; "we can exclude effects above this size" is.

The structural error is older and simpler. Across v1–v13 the **fixed parameter
was the dataset** and the **floating parameter was the question**: PPAR →
rewiring → age interaction → LAR reversal → multi-omic anchor → DCT1 →
DCT2/CNT → ASDN → Grey60 → four clinical axes. Ten hypotheses, one corpus. That
is the parameter strategy inverted. Fixing the question and letting the
evidence base float is the correction, and it is what §3 does.

---

## 3. What changed today: the podocyte lead is not orphaned

The lead looked uninterpretable because nothing external anchored it. That was
wrong. Siew et al. 2024 (*Nat Commun*, doi:10.1038/s41467-024-49212-1) — the
field's flagship, ~100 authors including much of NASA GeneLab — contains three
glomerular observations that the repository has not connected:

1. **A glomerular lesion in their GCR arm.** NSRL22A GCRsim mice showed
   significantly greater urinary protein excretion, described as "consistent
   with a glomerular lesion," plus microthrombi in glomerular capillary tufts
   and overt thrombotic microangiopathy in 27% (3/11) of females versus 0% of
   shams.
2. **Glomerular disease ontologies.** Nephrotic syndrome, membranous
   glomerulonephritis, and amyloidosis are among the enriched ontologies shared
   across their spaceflight and simulated-GCR experiments.
3. **A directional conflict on a named gene.** They report Podocalyxin (PODXL,
   the podocyte glycocalyx constituent) and Nidogen-1 (NID1) as "two of the
   most frequently downregulated gene products across kidney datasets."

The repository's mission-level result points the other way. The high-specificity
podocyte program is positive in four of five missions (g = 0.689, mHK CI
0.042–1.336, maxT FWER = 0.0189, I² = 0%), and the barrier-core direction is
retained when Podxl and Cd2ap are explicitly added (g = −0.704 in adverse
orientation, i.e. higher expression).

So the glomerulus is already present in the flagship's data — but assigned to
the *radiation* arm, while the *spaceflight* headline stayed tubular
(DCT expansion, NCC dephosphorylation). The repository's contribution is to
show a recurrent podocyte-associated signal in the spaceflight arm across
missions, and to disagree with the reported direction.

**This is a much better paper than "no axis passed."** It is positive,
specific, gene-named, falsifiable, and it engages the flagship on its own
material rather than merely failing to reproduce it.

### Why the direction disagreement is defensible

The PODXL claim rests on frequency of downregulation *across datasets* — a
vote-count over correlated corpora. The 2026-08-11 report already identifies
that inference mode as the reason the older matrix/ECM `p = 7e-4` was retired:

> That analysis treated correlated genes as inferential units; the animal- and
> mission-level analysis does not reproduce a general recurrent fibrosis
> program.

The same objection applies to a cross-dataset frequency count. The repository's
estimator — animal-level signed scores, mission-level Hedges g, REML with
modified Hartung–Knapp, 100,000 blocked permutations under maxT — is the
stronger design for a directional question. That argument is available and
should be made explicitly, not implied.

---

## 3b. AMENDMENT (same day): the frozen structural control nearly matches the lead

Added after inspecting `compartment_context_meta_results.tsv` directly. This
qualifies §3 substantially and takes priority over it.

Ranked by estimate across the 51-set family:

| Rank | Set | g | mHK CI | maxT FWER |
|---|---|---:|---:|---:|
| 1 | podocyte high_specificity | 0.689 | 0.042 – 1.336 | **0.0189** |
| 2 | podocyte all_enriched | 0.567 | −0.078 – 1.213 | 0.105 |
| 3 | podocyte broad_enriched | 0.554 | −0.090 – 1.198 | 0.120 |
| **4** | **broad_structural_scaffold_control__all** | **0.549** | −0.095 – 1.193 | 0.128 |
| 5 | podocyte scaffold_excluded | 0.537 | −0.108 – 1.181 | 0.149 |
| 6 | principal_cell broad_enriched | 0.530 | −0.102 – 1.162 | 0.143 |

`broad_structural_scaffold_control__all` is the frozen union of KEGG adhesion,
ECM, focal-adhesion, actin, tight-junction, and vascular smooth-muscle genes.
It exists for exactly one purpose, stated in
`ANALYSIS_TRIAGE_AND_NETWORK_REFRAME_2026-07-29.md`:

> If a protein-class set beats every compartment set, the compartment axis is a
> protein-class axis and every compartment framing dies in one line.

It does not beat every compartment set — podocyte high_specificity exceeds it,
and only that row clears the family. So the criterion as literally written is
not triggered. But the control ranks **4th of 51**, above every other podocyte
tier, and critically **above `podocyte__scaffold_excluded`** — the tier built by
removing those very structural genes. Principal cell sits alongside it. All
these rows have tau² = 0 and I² = 0.

The honest description of the data is therefore not "a podocyte program is
elevated." It is: **a cluster of broadly-expressed structural and compartment
sets drifts positive at g ≈ 0.45–0.57, and the high-specificity podocyte tier
sits at the top of that cluster at 0.689.** Whether 0.689 is a podocyte signal
or the upper tail of a general structural drift is not resolved by these
numbers, and the multiplicity result alone does not settle it — passing a
family-wise threshold that a negative control narrowly misses is a weak
discrimination when the point estimates differ by 0.14 against CI widths of
~1.3.

Note also that this pattern — several unrelated sets moving together, zero
heterogeneity, effects clustered in a narrow band — is the signature of a shared
global, compositional, or normalization-level shift. That is the same class of
explanation the report already entertains for the podocyte proxy.

### K0 — the test that decides whether §3 exists at all

Run **before** K1–K4, and before writing anything:

> Regress the high-specificity podocyte score on flight while adjusting for the
> `broad_structural_scaffold_control__all` score, blocked by mission, HC3, with
> the same permutation scheme.

The machinery exists: `scripts/clinical_axes/run_barrier_specificity_adjustment.py`
already does exactly this shape of adjustment (it produced the 0.507 → 0.070
collapse when barrier-core was adjusted for the disjoint podocyte proxy).

| K0 outcome | Consequence |
|---|---|
| Podocyte coefficient survives adjustment with an interval excluding zero | §3 stands. Proceed to K1–K4. The paper is real and the structural control becomes a strength — a declared negative control that was run and reported. |
| Coefficient collapses toward zero, as barrier-core did | The podocyte framing dies, but the paper does not. The headline becomes the cross-mission reproducibility result (§5 alt), with the structural drift reported as the one positive and explicitly labelled non-attributable. |
| Ambiguous — attenuates but retains direction | Report as a bounded observation, not a headline. Do not promote it. |

Also run the reciprocal: regress the structural control score on flight
adjusting for the podocyte score. If both attenuate symmetrically, neither set
is the carrier and the shared-drift explanation is the parsimonious one.

**Do not skip K0 because the FWER passed.** The FWER answers "is this row
unusual within the family." It does not answer "is this row a podocyte
signal rather than a structural one," and that is the question the paper's
title makes.

---

## 4. Kill criteria — resolve these before writing a word

The last two novelty pitches in this repository died on contact with search.
Do not skip this. **K0 in §3b comes first and can retire all of these.**

| # | Check | Kills the paper if |
|---|---|---|
| K1 | Is Siew's PODXL/NID1 statement about **mouse kidney tissue RNA**, or about human/Roscosmos **urine proteomics** and cross-layer gene products (their Fig. 5 / Supplementary Data 3)? Read the figure and supplement directly. | It is urine-proteomic or cross-species only. Then there is no contradiction, only a different measurement — reframe as extension, not disagreement. **Most likely outcome; plan for it.** |
| K2 | Forward-citation search from Siew 2024 and da Silveira 2025 for any published spaceflight glomerular or podocyte transcript claim. | Someone has already reported recurrent podocyte-program elevation. |
| K3 | Re-run the podocyte program with Podxl and Nid1 **forced in**, and report the per-gene mission-level effects for both. | Podxl and Nid1 are individually negative while the program is positive — that is a real internal inconsistency and must be reported prominently, not buried. |
| K4 | Confirm the high-specificity tier result is not an artifact of the atlas abundance/specificity extreme flagged in the matched-panel caveat. | Matched-panel support collapses under a stricter pool. |

K1 is the load-bearing one. If K1 says "no contradiction," the paper still works
— the title becomes *A recurrent podocyte-associated kidney program across five
spaceflight missions*, and Siew's GCR glomerular lesion becomes converging
support rather than a foil. Either way there is a paper. Decide the framing from
K1's answer rather than choosing it now.

---

## 5. What the paper contains

1. **Design and eligibility.** Stage-0 gates, 5′/3′ coverage, OSD-253's
   confounded arm, OSD-462's reporter-tag alias and preparation instability.
   This is why the estimate is trustworthy and it is already written.
2. **The four-axis frozen family, reported honestly.** Fibrosis null and
   heterogeneous (g = 0.311, I² = 62.5%). Distal transport null
   (g = −0.153, I² = 68.8%). Tubular injury crosses zero. This is the
   specificity control: it shows the pipeline is capable of returning nulls.
3. **The podocyte result** with all ten adversarial checks — signed median,
   common-observable genes, expanded panel, leave-one-mission, leave-one-gene,
   three RNA preparations, OSD-253 rerun control, OSD-163 mapping-rate
   residualization, observability-matched panels (p = 0.0062), and the frozen
   49-set compartment family.
4. **The boundary, stated in the abstract.** Adjustment for a disjoint podocyte
   proxy absorbs the barrier-core effect (0.507 → 0.070). Bulk RNA cannot
   separate altered cell representation from coordinated per-cell transcription.
   Prediction interval crosses zero. No barrier function was measured.
5. **The resolving experiment**, specified precisely enough that someone with
   tissue can run it (§7).

Figures already exist in `figures/clinical_renal_axes_cross_mission/`. The
decision report is most of a draft.

---

## 6. Sequencing to a January preprint

| Weeks | Work |
|---|---|
| 1–2 | K1–K4. Read Siew Fig. 5 + Supplementary Data 3 directly. Forward-citation search. Re-run with Podxl/Nid1 forced in. |
| 2–3 | Send the PI a one-page summary: the mission-level podocyte result, the direction question, the proposed morphometry. Ask two things — is this already known internally, and who holds sections. This runs in parallel and gates nothing. |
| 3–7 | Draft. Reuse the decision report as the results skeleton. |
| 7–9 | Commit the run directories, drop the `latex_paper/*` gitignore exclusion, tag the frozen config. Every "frozen before testing" claim is currently unverifiable from repository history — a reviewer will check, and it is an afternoon of work. |
| 9–12 | Internal review by the PI, revise, post to bioRxiv. |

Target: bioRxiv by December, journal submission January. npj Microgravity is
the obvious venue given the prior publication, though the direction
disagreement may argue for a nephrology venue instead — ask the PI.

---

## 7. Validation: ask the consortium, do not file an NBISC request

The 2026-08-11 report's four validation options assumed acquiring tissue and
finding someone to assay it. That is the slow path and it is unnecessary as a
first move.

The capability already exists inside the Siew consortium:

- **Automated glomerular morphometry** — Sarder's group (U. Florida, Medicine-
  Nephrology / Intelligent Critical Care) performed the AI-based tubule and
  nuclear-density quantification in the flagship paper. Glomerular and podocyte
  morphometry is the same pipeline pointed at a different structure.
- **Renal pathology** — Roufosse (Imperial) read the histology.
- **3D imaging** — Walker-Samuel (UCL Centre for Advanced Biomedical Imaging).
- **The tissue and sections already exist** and have been through H&E, Masson,
  MSB, PAS, Picrosirius red, and Von Kossa.

The ask is therefore small and specific: **WT1-positive nuclei per glomerulus
and per tuft area, plus nephrin/podocin/synaptopodin staining, on sections that
are already cut.** That is a request an established group can absorb, and it is
exactly what distinguishes more podocytes from more expression per podocyte.

Route it through the PI. An NBISC request remains the fallback if the
collaboration does not materialize, but filing one now would consume the
timeline for a result that arrives after the preprint anyway.

---

## 5 alt. If K0 collapses: the same paper, different headline

The two branches share most of their structure, which is the point. §5 items 1,
2 and 4 are identical either way. Only the headline and one results section
change:

> Most reported renal responses to spaceflight do not recur across independent
> mouse missions: a mission-level reanalysis of five terminal cohorts.

Contents: the eligibility gates, the four-axis frozen family returning null
(fibrosis g = 0.311 with I² = 62.5%; distal transport g = −0.153 with
I² = 68.8%; tubular injury crossing zero), and one positive — a coordinated
structural transcript drift at g ≈ 0.45–0.57 spanning podocyte, principal-cell
and the frozen scaffold control alike, which no compartment annotation
resolves and which is reported as such.

This directly tests the field's two standing claims — the distal-transport/NCC
axis (Siew 2024) and the ECM/TGF-β axis (da Silveira 2025) — with a stronger
design than either, and reports what survives. The estimator sensitivity is
itself a result: the earlier Stage-0 synthesis gave fibrosis g = 0.798 robust
to leave-one-cohort-out, and the animal-level mission-level estimator gives
0.311 and null. Explaining that gap is a contribution, not an embarrassment.

**Consequence for scheduling: start drafting before K0 returns.** The
introduction, methods, eligibility, and null-axis results are common to both
branches. K0 decides the title and one section, not the manuscript.

---

## 8. The retrospective

`network_journey_retrospective.tex` is genuinely unusual — a dated, versioned
trail of claims made, tested, and retracted by the same author. Post it as a
standalone preprint or long-form piece *after* the primary paper, not before.
Publishing the postmortem first invites reviewers to read the primary paper as
the fourteenth attempt rather than as the first properly specified one.

---

## 9. Stays retired

No v14. No further OSD-462 compartment mining (the 2026-08-02 stop rule holds).
No Grey60 rescue. No S1 simulation. No LayerScore as a new method. No
matched-library robustness paper — library-preparation discordance is well
documented (PMC3747248, PMC4620295).

**No OSD-462 data note.** Owner decision, 2026-08-11. It is removed from the
portfolio entirely and is not a fallback.

One consequence to handle rather than ignore: OSD-462 is one of the five
missions in the primary synthesis, so a reviewer will ask why a cohort with
condition perfectly aliased to reporter-tag block and a −0.09 / +0.98 / −1.16
preparation sign flip is eligible at all. The Stage-0 material therefore does
not disappear — it moves into the eligibility and sensitivity subsections of the
main paper, where it is load-bearing:

- OSD-462 enters the RNA synthesis on total-RNA with mRNA and UPX as paired
  within-animal sensitivities, and the podocyte/structural result is reported
  under all three (g = −0.673 and −0.686 for the barrier orientation).
- The reporter-tag alias and the phospho-layer block shift are cited only to
  justify excluding OSD-462 phosphoproteomics from the evidence base, which the
  2026-08-11 report already does.
- The residue-indexing and co-modified-peptide corrections are not needed for
  an RNA paper and can simply go unpublished.

An OSDR curation ticket remains available at any time and costs nothing. It is
not a publication and does not conflict with this decision.

---

## 10. Honest risks

- **K1 probably removes the contradiction.** Siew's PODXL statement appears in
  a urine-proteomics paragraph. Plan for the extension framing.
- **The podocyte effect is not large and the prediction interval crosses zero.**
  The claim is coherence across the observed missions, not a guarantee for the
  next one. Say so.
- **The proxy-adjustment result is genuinely limiting.** The barrier-core genes
  do not move beyond the broader podocyte program. Anyone who reads carefully
  will find this; put it in the abstract.
- **The compartment scan followed the targeted result.** This is post-hoc and
  the lock is retrospective. It is a strong hypothesis, not a confirmatory
  endpoint, and the paper must say that in those words.
- **Political care is warranted.** Disagreeing on direction with a ~100-author
  consortium paper is fine and normal science; framing it as a correction to
  their work is not the way to get it read. Frame it as the mission-level
  reproducibility layer their design could not provide, plus a question about
  one gene's direction.

---

## 11. Prior-art search performed, and what it found

Approximately ten targeted searches (PubMed, Consensus, bioRxiv, web) on
2026-08-11. Not a systematic review. Findings:

**No in vivo spaceflight podocyte or glomerular transcript claim was found.**
PubMed `microgravity AND podocyte` returns exactly two records, both in vitro:

- Melica et al. 2024, *Stem Cell Res Ther*, doi:10.1186/s13287-024-03633-3 —
  modeled microgravity (rotary cell culture) **impairs** renal progenitor
  differentiation into podocytes, disrupts F-actin cytoskeleton and focal
  adhesions, reduces nephrin and nestin. This is real prior art on
  microgravity → podocyte and it must be cited. It is in vitro, uses a
  simulation device, and points toward impairment rather than elevation. It
  also independently implicates the **cytoskeleton**, which is uncomfortable
  given §3b.
- Nishimura & Wang 2021, *BBRC*, doi:10.1016/j.bbrc.2021.08.012 — simulated
  microgravity spheroid culture of developing kidney cells. Tangential.

**Astronaut literature points the other way on barrier function.** Urinary
albumin excretion is *reduced* during spaceflight relative to pre- and
post-flight, only partly explained, with hemodynamic and fluid-shift
contributions proposed. This is consistent with the report's refusal to claim
barrier failure, and it is an argument the paper must engage rather than avoid.

**No cross-mission animal-level meta-analysis of spaceflight kidney was found.**
The nearest works are Siew 2024 (pan-omic pooling across 11 mouse missions,
doi:10.1038/s41467-024-49212-1) and da Silveira 2025 (OSD-102/OSD-163 strain
comparison, npj Microgravity, doi:10.1038/s41526-025-00465-0). Neither performs
the mission-level effect-size synthesis with FWER control that this repository
has built.

**Still owed before submission:** forward-citation searches from Siew 2024,
da Silveira 2025, and Melica 2024.

---

## 12. What is actually being asked of the owner

Authorize **K0** (§3b), and authorize drafting the common sections in parallel.

K0 decides the title, not whether there is a paper. Both branches — the podocyte
headline (§5) and the reproducibility headline (§5 alt) — share their
introduction, eligibility gates, methods, four-axis null results, and boundary
statements. Waiting for K0 before writing spends three weeks of a
fourteen-week timeline on a question that only affects one results section.

There is no separate fallback product, by owner decision. The fallback is a
different headline on the same manuscript. If K0 collapses *and* the
reproducibility framing also proves unwritable, that is the point to stop and
re-plan rather than to reach for another asset.
