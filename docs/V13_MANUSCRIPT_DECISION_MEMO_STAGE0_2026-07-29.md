# V13 manuscript decision memo after the OSD-462 Stage 0 audit

**Decision:** v13 must not inherit v12's claim that OSD-462 independently shows
canonical NCC/SPAK regulatory dephosphorylation. Stage 0 found **zero isolated,
assay-qualified canonical NCC/SPAK phosphoforms**. The observable central
features are co-modified peptides indexed to NCC T53 and SPAK S383, and flight
status is perfectly aliased with reporter-tag blocks in both plexes. The only
plausible main-paper novelty is therefore the frozen continuous parent-gene
test of independently defined DCT2/CNT-transition and ASDN programs. Its
machine-readable claim tier must control the title, abstract, Results, figures,
and conclusion.

This memo is anchored to
[`OSD462_STAGE0_PROVENANCE_AUDIT_2026-07-28.md`](OSD462_STAGE0_PROVENANCE_AUDIT_2026-07-28.md)
and the claim gates in
[`v13_continuous_phospho_inference.md`](v13_continuous_phospho_inference.md).

## Final exact-run resolution

The definitive run
`data/results/run_20260729_v13_continuous_phospho_exact_final/` completed all
63,504 balanced within-plex assignments and resolved the provisional decision:

- DCT2/CNT transition: 5 observable genes, below the frozen minimum of 8;
  component status `non_evaluable`.
- ASDN: 16 observable genes, competitive statistic 0.718 and conditional exact
  maxT \(p=0.0291\), but component status `fail` because the prespecified
  podocyte comparator is stronger.
- Overall statistical tier: **`neither`**.
- Publication promotion: **blocked by design/provenance** because condition is
  aliased with reporter-tag blocks and Stage 0 found zero isolated qualified
  canonical anchors.

The compact v13 document must therefore be a technical/null-boundary report.
The title boundaries shown below for passing biological tiers are retained as
the prospective decision rule, not as options permitted by the completed run.

## Required disposition of v12

| V12 item | Action for v13 | Required change |
|---|---|---|
| Title and Abstract | **Delete and rewrite after the claim tier is final** | Remove “dephosphorylates NCC,” “regulatory phosphorylation resolves,” “independent WNK recovery,” and any unconditional DCT2/CNT extension. State that the data are whole-kidney and that condition is aliased with reporter-tag blocks. |
| Figure 1 (`fig:model`) | **Delete in its current form; rebuild last** | It presents the canonical NCC/SPAK result as an assay-validated anchor and the DCT2/CNT extension as settled. A replacement may show only the tested parent-gene sets, exclusions, whole-pipeline permutation, compartment comparators, and the final claim tier. |
| Figure 2 (`fig:axis`) | **Delete as a data figure** | The WNK–SPAK/OSR1–NCC pathway may appear as a small literature-background schematic, but OSD-462 must not be placed on an “active-to-inactive NCC” diagram. The present data do not isolate regulatory-site occupancy or NCC activity. |
| Table 1 (`tab:datasets`) | **Rewrite and promote** | Give exact sample ID, WT genotype, female sex, kidney side, condition, plex, reporter tag, raw-file identifier, and inclusion reason. State that the MS workbooks contain 30 WT samples (5 baseline, 5 flight, and 5 ground per plex) and that the primary contrast is 10 flight versus 10 ground; the p21-null arm is absent. Use the assay nomenclature below. |
| Table 2 (`tab:phospho`) | **Replace with a provenance table** | Report peptide sequence, accession/isoform, modified residues, AScore/localization evidence, reporter completeness, and effect for the T53- and S383-indexed co-modified features. Correct NCC 65/68 to Y65/Y68. Do not divide rows into “regulatory” and “non-regulatory” evidence classes unless a row meets the frozen isolated-site qualification rule. |
| Figures 3–6 (`fig:rnarecurrence`, `fig:pathwayscoresmethods`, `fig:crosslayer`, `fig:propagation`) | **Demote to Supplement or a separate RNA paper** | They are context, not validation of the phosphosite endpoint. Delete “phospho-carry,” “does not propagate,” and other categorical cross-layer claims; use “little evidence of concordant protein-abundance change.” |
| Figure 7 (`fig:osd462`) | **Rebuild or demote** | Total NCC abundance may be shown as having no detected whole-kidney decrease. Remove Panel C's regulatory/non-regulatory cluster and the abundance-versus-activation contrast. If retained, display the co-modified features individually and label them as non-isolating context. |
| Figure 8 (`fig:kinomerecurrence`) | **Demote to Supplement and rewrite** | KSEA uses the same phosphosite effects and is not independent or orthogonal confirmation. Describe only motif-compatible WNK-family substrate enrichment, report WNK1/WNK3 substrate overlap and shared-substrate removal, and do not infer kinase activity or use it as a positive control. |
| Figure 9 and Table 3 (`fig:endoscatter`, `tab:matched`) | **Delete in their current form** | “NCC regulatory phosphorylation” is not an assay-qualified variable, and the pooled correlation is largely condition separation under tag aliasing. A recomputed co-modified-feature correlation, if scientifically useful, belongs in the Supplement as exploratory and condition-stratified. |
| Table 4 and Figure 10 (`tab:dct2compare`, `fig:dctladder`) | **Supersede with the v13 continuous analysis** | The nominal-\(p<0.05\), phosphosite-row Fisher analysis and low-DCT1 decile must not remain the visible primary result. Replace them with one-gene-one-score competitive enrichment, whole-pipeline animal-label permutation, maxT adjustment, independent signatures, compartment comparators, and leave-one-gene-out results. |
| Figure 11 (`fig:composition`) | **Demote and replace** | The M0–M5 reselected-set ladder does not resolve composition or reporter-tag confounding. Retain only prespecified v13 parent-protein, broad-expression, centered/uncentered, and observability sensitivities; do not call mediators “resolved confounders.” |
| Figure 12 (`fig:perturbation`) and Figure 13 (`fig:aldo`) | **Delete from the main paper** | Low-K, IRI spots, and cohort-level MR correlation are hypothesis context, not subtype or phosphosite validation. At most retain compact Supplementary summaries without inferential language. |
| Table 5 (`tab:physpred`) | **Delete from the main paper** | Sodium, potassium, calcium, magnesium, blood pressure, transporter flux, and NCC activation were not measured. Convert any retained content into a short “measurements required next” paragraph. |
| Figure S1 (`fig:tmtqc`) | **Rewrite and promote its design finding** | Centering is a normalization sensitivity, not a solution to tag confounding. Explicitly show the condition-by-tag map and the absence of a cross-plex label swap. |
| Tables S3–S4 (`tab:dct2validation`, `tab:leadingedge`) | **Replace** | Replace the low-DCT1 bin with frozen DCT1, strict DCT2, DCT2/CNT-transition, ASDN, broad-expression, and whole-kidney compartment annotations. Rank leading genes by the continuous permutation-calibrated statistic, not “strongest nominally significant site.” |
| Table S6 (`tab:dct1counts`) | **Delete or retain only as historical sensitivity** | It is a thresholded contingency analysis and cannot carry the conclusion. |
| Introduction, Discussion, and Conclusion | **Rewrite to match the final tier** | Prior antibody work may be cited as prior evidence from the same RR-10 material, not as independent validation of these MS features. Do not retreat to a canonical NCC paper if both new sets fail: Stage 0 also prevents that fallback. |

## Four-section main Results skeleton

### 1. Stage 0 defines the assay, sample contrast, and inferential boundary

Report the exact 30-sample/2-plex map, the 10-flight versus 10-ground primary
contrast, TMTpro provenance, the legacy iTRAQ inconsistency, reporter-tag
aliasing, and peptide-level residue audit. End with the factual result that no
isolated canonical NCC/SPAK feature qualified. The T53- and S383-indexed
co-modified rows are context, not the validation endpoint.

### 2. Continuous parent-gene analysis maps flight-associated phosphorylation across kidney compartments

Present the label-blind site universe, continuous site effects, one score per
parent gene, gene-specific label-permutation calibration, and the prespecified
whole-kidney compartment comparison. This section establishes whether the
pattern is late-distal-selective or part of broader endothelial, stromal,
proximal, immune, or pan-kidney remodeling. Use “flight-associated” throughout
because condition and tag blocks are aliased.

### 3. Frozen DCT2/CNT-transition and ASDN tests determine whether signal extends beyond canonical anchors

Test the independently defined DCT2/CNT-transition and ASDN sets as one maxT
family after the primary and strict canonical-axis exclusions. Show competitive
effect sizes, empirical permutation \(p\)-values, broad-expression exclusion,
independent-reference validation, compartment rank, and leave-one-gene-out
results. Do not use the old low-DCT1 decile or nominal-site membership as the
primary evidence.

### 4. Robustness and the permitted claim tier

Report centered versus uncentered direction, parent-protein subtraction,
alternative gene scores, localization/de-duplication profiles, observability
matching, and tag-bias diagnostics. Close with exactly one tier:

| Tier | Permitted main conclusion | Title boundary |
|---|---|---|
| `DCT2_CNT` | “In whole-kidney tissue, parent proteins assigned independently to a DCT2/CNT-transition program showed disproportionate flight-associated phosphorylation suppression after canonical-axis exclusion.” This does not localize the sites to DCT2/CNT cells. | **Flight-associated whole-kidney phosphosite suppression is enriched among DCT2/CNT-transition parent proteins** |
| `ASDN_only` | “Flight-associated phosphorylation suppression was enriched among ASDN-associated parent proteins, but a distinct DCT2/CNT-transition enrichment was not established.” | **Flight-associated phosphosite suppression is enriched among aldosterone-sensitive distal-nephron parent proteins** |
| `neither` | “Neither predefined late-distal program showed selective enrichment after multiplicity control and robustness gates.” Do not return to “NCC dephosphorylation” as the paper's result. | A null-boundary or methods/resource paper only; no “beyond NCC/DCT1” title. |
| `non_evaluable` | “The available design did not permit a valid tier decision.” State the specific failed gate, including external-reference or reporter-tag/exchangeability failure. | Do not submit a biological subtype-claim paper; frame as a provenance/methods report or stop. |

## Exact language boundaries

### Co-modified phosphoforms

**Use:**

- “The NCC T53-indexed feature was derived from a peptide carrying both T53
  and Y65 modifications.”
- “The SPAK S383-indexed feature was derived from a peptide carrying both S382
  and S383 modifications.”
- “The flight-associated reporter-signal difference for each co-modified
  peptide cannot identify which residue or phosphoform produced the change.”
- “No isolated canonical NCC T53/T58/S71 or SPAK T243/S383 feature satisfied
  the prespecified assay-qualification criteria.”
- “Prior work reported antibody evidence for NCC dephosphorylation in RR-10;
  because it analyzes the same biological experiment, it is external
  contextual evidence rather than independent replication.”

**Do not use:**

- “pT53-NCC decreased,” “pS383-SPAK decreased,” or “canonical regulatory sites
  were suppressed.”
- “NCC was deactivated,” “SPAK activity fell,” “reduced NCC activation,” or
  downstream electrolyte predictions stated as results.
- “The regulatory cluster was validated,” “orthogonal confirmation,” or
  “independent recovery of the WNK–SPAK axis.”
- “Phosphosite occupancy” or “stoichiometry” for reporter-ion abundance.

### Assay and design

**Methods sentence:**

> The detailed OSD-462 protocol specifies TMTpro labeling, a +304.207-Da tag
> mass, and TMTpro isotope-impurity correction; a legacy investigation field
> instead describes iTRAQ. We therefore refer to the measurements as TMTpro
> isobaric-tag proteomics and phosphoproteomics and disclose the inconsistent
> legacy field.

**Design-limitation sentence:**

> Flight status and reporter-tag blocks were perfectly aliased in both plexes,
> with no cross-plex label swap. The observed coefficient therefore cannot
> distinguish a biological flight effect from a systematic reporter-tag-block
> effect; we report it as a flight-associated contrast and evaluate tag-related
> sensitivity without claiming that normalization resolves the aliasing.

Use “TMTpro isobaric-tag proteomics/phosphoproteomics” after the first
definition. Do not alternate among “TMT,” “iTRAQ,” and generic “TMT-like”
terminology.

## Go/no-go rule

Freeze the prose only after `claim_tier.tsv` and all exact-run diagnostics are
final. A DCT2/CNT or ASDN manuscript proceeds only if its declared gates pass
and reporter-tag diagnostics do not make the contrast non-evaluable. Otherwise
the scientifically correct outcome is a null-boundary or technical report, not
a narrower canonical NCC/SPAK paper.
