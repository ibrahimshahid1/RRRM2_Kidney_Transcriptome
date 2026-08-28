# Prior-art verdict: LayerScore, the network bound, and the proposed simulation

**Date:** 2026-07-29
**Scope:** ~12 targeted searches against the claims in
`METHODS_PAPER_ANALYSIS_PLAN_2026-07-29.md` and the LayerScore preregistration.
**Confidence:** moderate-to-high on items marked **covered**; a full systematic
search (Scopus/WoS forward-citation from the anchors below) would take a week and
is unlikely to reverse them.

---

## Headline: I was wrong on three counts

1. Both novelty axes I proposed for the simulation — **prior accuracy** and
   **derivation source** — have substantial prior art.
2. The idea I called "the best new empirical idea in the plan" (S3: use replicate
   measurements to estimate reliability, then bound cross-layer attenuation) is
   **established proteogenomics methodology**, published at least three times.
3. LayerScore's central identity — equivalence testing to separate supported
   absence from insufficient evidence — was **published in proteomics in 2024**.

The honest consequence: **there is no new-methods paper here.** There are still
three publishable products, and one of them does not depend on novelty at all.
Details below.

---

## 1. LayerScore

### Claim: equivalence testing separates "supported absence" from "insufficient evidence"

**Covered.**

- **QuEStVar** (2024) — "Statistical Testing for Protein Equivalence Identifies
  Core Functional Modules Conserved across 360 Cancer Cell Lines and Presents a
  General Approach to Investigating Biological Systems."
  [PMC11166143](https://pmc.ncbi.nlm.nih.gov/articles/PMC11166143/). Combines
  differential *and* equivalence testing to expand the statistical
  classification of analytes when comparing conditions, in proteomics, at scale.
  This is LayerScore's conceptual core, already executed and published.
- **Equivalence in genomics** — p-values, q-values and posterior probabilities
  for equivalence in genomics studies,
  [arXiv:1202.0048](https://arxiv.org/pdf/1202.0048) (2012).
- **TREAT** — McCarthy & Smyth, thresholded null \(H_0: |\beta| \le \tau\), which
  the paper itself connects to equivalence testing with the hypotheses reversed.
  [Bioinformatics 2009](https://academic.oup.com/bioinformatics/article/25/6/765/251641).
- **DESeq2** already ships an equivalence-direction test
  (`altHypothesis="lessAbs"`), so this is in the standard toolchain, not the
  frontier.

**What is left:** the *target*. QuEStVar classifies analytes as changed or
equivalent within one layer. LayerScore's target is **cross-layer localisation** —
using equivalence at the protein layer to license a claim about *where in the
regulatory hierarchy* a response lives — plus observability-matched nulls for
proteomic detection structure. That is a narrower and more defensible
contribution, and it must be framed as an application of an existing decision
framework rather than a new one, with QuEStVar cited prominently.

### Claim: correct cross-layer concordance for measurement-error attenuation

**Covered, repeatedly.**

- **Franks, Airoldi & Slavov (2017)** — "Post-transcriptional regulation across
  human tissues,"
  [PLOS Comp Biol](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005535).
  Applies Spearman's disattenuation to mRNA/protein correlation and estimates the
  true across-tissue correlation at 0.21–0.29. This is exactly the arithmetic in
  proposed S3.
- **Upadhya & Ryan (2022)** — "Experimental reproducibility limits the
  correlation between mRNA and protein abundances in tumour proteomic profiles,"
  [Cell Reports Methods](https://www.cell.com/cell-reports-methods/pdf/S2667-2375(22)00170-9.pdf).
  Uses measurement reproducibility to explain observed mRNA–protein correlation.
  Same logic, same conclusion structure.
- **Antibody reliability influences observed mRNA–protein correlations** (2023),
  [Life Science Alliance](https://www.life-science-alliance.org/content/6/8/e202201885).
- Disattenuation for cross-platform expression comparison goes back to at least
  [2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2211756/).
- **RepExplore** exploits technical-replicate variance for more reliable
  differential statistics,
  [PMC4481852](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4481852/).

**The one real distinction:** those papers disattenuate correlations of
*abundances* across genes or tissues. LayerScore would disattenuate a correlation
of *treatment effects*. That is a genuine difference in estimand — effect
estimates have their own sampling variance on top of assay noise — but the
correction is the same algebra, and a reviewer will say so. It is a sentence in a
methods section, not a paper.

The OSD-462 three-preparation structure is still worth using; it is just a
*measurement*, not a method. And note the disclosure below.

---

## 2. The proposed simulation (S1)

### Axis 1 — power as a function of prior accuracy

**Substantially covered.**

- **"Assessment Method for a Power Analysis to Identify Differentially Expressed
  Pathways"** (2012),
  [PMC3356338](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3356338/). Examines
  the influence of sample size, **detection call — the percentage of genes
  differentially expressed within a pathway** — and pathway size. Detection call
  *is* the set-purity axis.
- **Statistical power of GSEA is a function of gene set correlation structure**,
  [bioRxiv 186288](https://www.biorxiv.org/content/10.1101/186288.full.pdf).
- Noisy/incomplete gene-set simulations exist, including designs that remove
  pathway genes and add controlled noise.
- Gene-set method power comparisons under correlation,
  [PMC3196970](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3196970/); camera
  handles inter-gene correlation explicitly.

### Axis 2 — derivation source (external vs data-derived vs transfer)

**Covered, and by stronger work than a simulation would be.**

- **Curated gene sets versus transcriptomics-derived signatures** (2020),
  [BMC Bioinformatics](https://link.springer.com/article/10.1186/s12859-020-3366-4).
  Direct head-to-head: data-derived AUC 0.67 versus curated 0.59.
- The **double-dipping / post-selection inference** literature treats
  fit-and-test-on-the-same-data rigorously rather than empirically: selective
  inference for hierarchical clustering
  ([arXiv:2012.02936](https://arxiv.org/abs/2012.02936)), count splitting and
  data thinning (Neufeld et al.), FDR control via data splitting after
  clustering ([arXiv:2410.06451](https://arxiv.org/html/2410.06451)), synthetic
  control for double dipping in single-cell and spatial data
  ([PMC10401959](https://pmc.ncbi.nlm.nih.gov/articles/PMC10401959/)).

Prediction P4 ("a \(w\) fitted on the samples you test is anticonservative") is
just **double dipping**, a named and heavily studied phenomenon. Demonstrating it
by simulation in 2026 adds nothing.

**One important technical correction, in the other direction.** Witten's result
that *sample splitting does not fix the problem* applies to testing for a
difference in means **between estimated clusters** — the estimand is defined by
the clustering. It does **not** apply to the transfer case in the plan, where a
*fixed* linear functional \(Xw\) with \(w\) derived elsewhere is tested in new
data; that test is valid. Worth stating precisely, because a reviewer who knows
this literature will otherwise assume the transfer template is also broken.

### Verdict on S1

**Do not run it.** Six to ten weeks to reproduce known results under a new
parameterisation. If any of it survives, it is a single small illustrative figure
inside product C below — not a study.

---

## 3. The network information bound

**The most novel thing in the repository, and it is a synthesis, not a theorem.**

No prior art surfaced for the specific composition: the DPI chain
\(Y \to X \to R \to W \to F \to \Delta\) applied end-to-end to a
co-expression → node2vec → Procrustes → cosine-rewiring pipeline, with the
per-edge variance, BH-floor, rotation non-identifiability and anchor-circularity
costs assembled in one place.

But every *component* is published: the DPI is textbook (Cover & Thomas),
skip-gram = shifted PMI is Levy–Goldberg, the adjacency-polynomial identity is
NetMF, orthogonal Procrustes is classical, Fisher-\(z\) variance is standard.

So the contribution is **pedagogical and diagnostic**, not novel mathematics.
That is a legitimate genre — tutorial, perspective, or "cautionary note" — but it
is not a methods paper, and the title must stop implying a new result. Something
closer to *"Why co-expression rewiring pipelines fail at n = 5: a worked
decomposition"* is honest and still useful, because people keep running these
pipelines.

---

## 4. The OSD-462 data note

### Reporter-tag interference: phenomenon covered, dataset instance not

- **Brenes et al. (2019)** — "Multibatch TMT Reveals False Positives, Batch
  Effects and Missing Values,"
  [PMC6773557](https://pmc.ncbi.nlm.nih.gov/articles/PMC6773557/). Detected
  Y-chromosome peptides in female channels, demonstrating reporter interference
  empirically. This is the canonical demonstration.
- **optTMT** (2026) — "optimizing any experimental design to minimize false
  positives caused by TMT reporter ion interference,"
  [Bioinformatics Advances](https://academic.oup.com/bioinformaticsadvances/article/6/1/vbaf243/8270647).
  A tool exists for designing against exactly this.

So do not present tag interference as a discovery. **Present OSD-462 as an
instance of the failure mode Brenes documented and optTMT was built to prevent** —
condition perfectly aliased with channel block, no cross-plex swap. Citing
established methodology as the standard the dataset fails is stronger than
claiming a finding.

### Library-preparation discordance: covered

Protocol comparisons on the same samples report as little as **17–20 % overlap in
differentially expressed genes** between preparations
([PMC3747248](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3747248/)); total-RNA
versus poly(A) capture differences are well characterised
([PMC4620295](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4620295/)).

The OSD-462 sign flips (−0.09 / +0.98 / −1.16 on the same 20 animals) are
therefore an *expected* magnitude of discordance, not an anomaly. Report them as
a dataset-specific caution with the prior literature cited, and drop the framing
that this is a novel observation about estimand definition. The estimand point is
still correct — it is just a known correct point.

### What is genuinely unaddressed

- **No papers surfaced** documenting metadata errors or reuse limitations in
  GeneLab/OSDR specifically. ISSOP exists as the standardisation consortium
  ([NAR 2025](https://academic.oup.com/nar/article/53/D1/D1697/7903386)) and is
  the natural audience.
- **The Y65/Y68 residue misindexing and the co-modified peptide problem are
  specific factual corrections**, not methodological claims. Nothing in the
  literature covers them because they are properties of one workbook.

**This is the strongest surviving item in the entire portfolio, precisely because
its value does not depend on novelty.** A factual correction to how a reused
public dataset is being interpreted is worth publishing even if every statistical
technique in it is thirty years old.

---

## 5. Revised portfolio

| Product | Genre | Novelty required | Status |
|---|---|---|---|
| **A. OSD-462 data note / matters arising** | Data note, comment, or preprint + OSDR curation ticket | **None** — factual corrections | Recommended. Mostly written. |
| **B. Reanalysis case study** | Reproducibility / reanalysis case study | Low — the contribution is the completeness of the audit trail | Recommended. Assets exist. |
| **C. Rewiring-pipeline cautionary note** | Tutorial / perspective | Low — synthesis of known results | Optional, cheap, genuinely useful. |
| **~~D. New methods paper~~** | — | High, and absent | **Drop.** |
| **~~S1 simulation~~** | — | High, and absent | **Drop**, or reduce to one figure in C. |

**Product B is the one worth thinking hardest about.** Its pitch is not "here is
a new method" but: *here is what happens when every recommended control is
applied to one heavily reused public dataset — the headline claim fails four
independent ways, and here is the complete decision record, including the twelve
prior versions and their dated falsifications.* There is an audience for this;
the small-cohort replicability literature
([PLOS Comp Biol 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1011630))
exists because people want it. And almost nobody can supply what you can: a
versioned trail of claims that were made, tested and retracted by the same
author.

LayerScore survives inside B as the layer-localisation component, citing QuEStVar
as the framework it applies.

## 6. What this saves

Dropping S1 removes 6–10 weeks of simulation whose two headline results were
already published in 2012 and 2020. The surviving portfolio needs writing, not
computation, and the one analysis still worth running is the within-block
reporter-position diagnostic — days, not weeks.

## 7. Caveat on this verdict

Roughly twelve searches, not a systematic review. Items marked **covered** rest
on named papers that clearly do the thing; I would not expect those to reverse.
Before submitting anything, run forward-citation searches from QuEStVar,
Franks 2017, Upadhya 2022, the 2012 pathway-power paper, the 2020
curated-vs-derived comparison, Brenes 2019 and optTMT. If a systematic search
turns up a paper that already does the OSD-462 residue audit, the portfolio
collapses to product C — but that is unlikely, since it requires someone to have
opened that specific workbook with that specific question.
