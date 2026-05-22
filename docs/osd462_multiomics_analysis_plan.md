# OSD-462 Multi-Omics Anchor — Analysis Plan

**Status:** plan / pre-registration draft · May 2026
**Backbone:** the existing RRRM-2 pathway / state / mechanism-axis analysis (`run_20260519_000547_2500g`) is *not* re-derived here; it is the fixed input. This plan only adds an independent protein/phosphoprotein anchor.
**Companion docs:** `docs/statistical_assessment_v5.md`, `docs/osd462_data_manifest.md`.

---

## 0. Purpose in one sentence

Test whether the RNA-level matrix-high / DCT-low remodeling signal — especially the DCT/WNK-SPAK/OSR1-NCC transport axis — corresponds to **protein-level and phosphoprotein-level** evidence of renal transporter remodeling in an independent spaceflight kidney cohort (OSD-462 / RR-10), so that the central claim moves from "a co-expression network shifted" to "a transcript-, protein-, and signaling-concordant remodeling phenotype."

---

## 1. Why this is the right move, stated plainly

An RNA co-expression / embedding shift is weak standalone evidence: it moves with cell composition, batch, stress, sampling time, and mean-expression correlation. "WNK/NCC network rewiring" is not a biological claim until it connects to the actual signaling chain:

> WNK kinases → SPAK/OSR1 phosphorylation → NCC (Slc12a3) phosphorylation → Na-Cl reabsorption in the DCT → electrolyte/volume/Ca handling → (clinically) stone-risk context.

RNA-seq can only see transcript abundance. The anchor adds two independent layers: **proteomics** ("do the transport proteins themselves move?") and **phosphoproteomics** ("did signaling *activity* change?"). This does not make "ECM up / DCT down after spaceflight" a surprising discovery — it is not. It converts a weak-inference result into a defensible, falsifiable, mechanistically specific one, and the DCT/WNK/NCC framing supplies the clinical hook (electrolyte handling, stone risk) that justifies publishing a non-surprising result.

---

## 2. What OSD-462 is — and its hard constraints

OSD-462 = NASA **RR-10** kidney multi-omic study. Verified from the local ISA metadata and data files:

| Property | OSD-462 / RR-10 | RRRM-2 / OSD-771 |
|---|---|---|
| Tissue | Left kidney | Kidney |
| Sex | Female | Female |
| Strain | B6129SF2/J (WT) | C57BL/6NTac |
| Age | 14–15 wk at launch (single age) | Young + Old |
| Arms | Single collection arm (terminal, cryo-preserved) | ISS-T + LAR |
| Conditions | Basal, Ground Control, Space Flight, Vivarium | FLT, HGC, VIV, BSL |
| Modalities | RNA-seq (UPX/mRNA/totRNA), TMT proteomics, phosphoproteomics | RNA-seq only |

**Constraint 1 — it validates the ISS-T side only.** OSD-462 has no LAR-type arm and no old animals. It can corroborate the ISS-T matrix-high/DCT-low *flight* signal. It **cannot** validate the LAR attenuation/divergence finding or the age-stratified observation. Be explicit in the manuscript: the multi-omics anchor strengthens the part of the paper that was already strongest, and the LAR section remains RNA-only and underpowered.

**Constraint 2 — cross-study, cross-strain, cross-modality.** No sample-level pooling. The comparison is **direction of effect** (flight − ground), exactly the cosine/contrast-vector logic already used for OSD-513. This is a *recurrence/concordance* test, not a replication in the strict sense.

**Constraint 3 — bulk whole-kidney mass spec.** DCT is a small fraction of kidney mass, so flight effects on DCT proteins are diluted. Mitigated here (see §3) because NCC is unexpectedly well-quantified, but quantification noise will still cap power.

**Constraint 4 — TMT 2-plex batch structure.** Proteomics is two TMT plexes ("Samp1-5", "Samp6-10"), each containing BL/FL/GC channels. Plex must be handled as a batch factor.

---

## 3. Data inventory — what is confirmed present

All paths under `data/external/osdr/OSD-462/`.

**RNA-seq (downloaded):** VST count matrices and OSDR differential-expression tables for UPX / mRNA / totRNA. mRNA VST = 26,966 genes × 42 samples, **ENSMUSG-indexed** → direct join to RRRM-2 (also ENSMUSG). DE table carries `Log2fc_(Ground Control)v(Space Flight)` and the reverse contrast. Sample table conditions: VIV 12, BSL 10, SF 10, GC 10.

**Proteomics (downloaded):** `GLDS-462_proteomics_..._Protein_WorkUp.xlsx`, sheet `protein_quant_2721`, 7,861 proteins. Columns: `Gene Symbol`, `Number of peptides`, per-channel `*_sn sum` (TMT S/N) and `*_sn scaled` (sum-to-100). Sample channels labeled BL / FL / GC across 2 plexes.

**Verified detection of the target machinery** (peptide counts plex1 / plex2):

- DCT/WNK/NCC chain: **Slc12a3 / NCC 34/36** (well quantified), Oxsr1 12/11, Wnk1 8/11, Calb1 39/56, Kcnj10, Pvalb — all present. Marginal: Wnk4 3/2, Stk39/SPAK 3/9, Klhl3 4/0.
- Matrix/ECM: Fn1 52/58, Col1a1 11/12, Fbn1, Mmp2, Dcn, Lum present; Col3a1 marginal (1/2).
- Broad transport controls: Slc12a1, Aqp2, Slc9a3, Atp1a1, Umod, Lrp2 — all present.
- **Clock proteins: only Arntl detected (1/4).** → the clock/circadian-DCT angle gets *no* protein support; it must stay RNA-only and out of the multi-omics claim.

**Phosphoproteomics (downloaded):** `GLDS-462_phosphproteomics_..._Pho_WorkUp_JM.xlsx`. The main sheets are `siteQuant_360` and `siteQuant_360_compositeSite`; Layer 2 keeps the hard checkpoint that abandons the activity claim if NCC and SPAK/OSR1 phosphosites are not quantified.

---

## 4. Pre-registered hypotheses and falsification rules

Registered **before** running, so a null result is interpretable rather than a failure.

**H1 — protein concordance (matrix).** The matrix/ECM gene set shows positive RNA↔protein flight-effect concordance in OSD-462.
*Predicted:* matrix proteins flight-up; concordance with RRRM-2 ISS-T RNA effect exceeds the random-gene-set null. *Falsified if:* matrix proteins show no directional flight effect or concordance is within the null band.

**H2 — protein concordance (DCT/transport).** DCT/WNK/NCC transport proteins shift in the *same direction* as the RNA signal (transcript-down → protein-down).
*Predicted:* Slc12a3/Calb1/Wnk1/Oxsr1 protein flight effect < 0; gene-set concordance exceeds null. *Falsified if:* DCT proteins flat or opposite, or concordance within null.

**H3 — phospho activity (conditional on the phospho file).** The WNK-SPAK/OSR1-NCC phospho-axis shows reduced phosphorylation in flight, consistent with reduced DCT transport *activity*.
*Predicted:* phospho-NCC (Slc12a3 Thr-cluster), phospho-SPAK/OSR1 sites down in flight. *Falsified if:* those phosphosites are quantified and show no decrease — **report this as a negative result; it would mean the RNA/protein abundance change does not extend to signaling activity.**

**H4 — network-candidate translation.** RRRM-2 LIONESS/node2vec top candidate genes are enriched for OSD-462 protein- or phospho-changing genes more than random genes.
*Predicted:* enrichment > null. *Falsified if:* not enriched — **report honestly: the network layer nominated transcripts without independent protein support.**

**Null-result policy:** every hypothesis above has an explicit "falsified if." A null on H3 or H4 is a publishable, informative outcome and will be reported, not buried.

---

## 5. Analysis layers

### Layer 0 — Harmonization (`scripts/osd462/00_harmonize.py`)

- **ID bridge.** RNA (both cohorts) is ENSMUSG; proteomics is Gene Symbol. Use `data/processed/resources/id_map.tsv` for symbol↔ENSMUSG. Resolve many-to-one collisions explicitly; log unresolved targets.
- **Flight-effect definitions.** OSD-462 RNA flight effect = Space Flight − Ground Control (use the SF-v-GC contrast; keep BSL/VIV only for sensitivity). RRRM-2 reference = the existing ISS-T FLT−GC effect (gene level from `lar_reversal_gene_scatter.tsv`; pathway/mechanism level from the state/mechanism score tables).
- **Proteomics flight effect.** From `*_sn scaled` channels: log2-transform, compute FL−GC **within each plex**, then average the two plex estimates (or fit `protein ~ flight + plex`). Require ≥2 peptides; carry `Number of peptides` as a quality weight. Drop proteins quantified in only one plex from primary tests; keep them in a sensitivity table.
- Output: a tidy `osd462_flight_effects.tsv` with RNA-effect, protein-effect, peptide count, plex coverage per gene.

### Layer 1 — Protein-level concordance (`scripts/osd462/01_protein_concordance.py`)

- **Targeted gene sets** (reuse the existing curated sets in `config/mechanism_gene_sets.yaml`): `dct_ncc_wnk`, `ecm_matrix`, `tlr4_innate`, `s1p`, plus `broad_transport` as a specificity control.
- For each set, compute (a) mean protein flight effect, (b) Spearman concordance between RRRM-2 ISS-T RNA effect and OSD-462 protein effect across set members, (c) sign-agreement rate.
- **Null model (the core statistic).** For each set of size *k*, draw 10,000 random size-*k* gene sets matched on protein abundance/peptide-count decile; recompute the concordance statistic; the empirical p is the fraction of random sets exceeding observed. This answers "do *these specific* pathways move concordantly more than arbitrary genes," which is the honest version of the test.
- Also report the genome-wide RNA↔protein flight-effect correlation as context — expect r ≈ 0.2–0.4; a modest targeted signal *above* that background is the realistic success criterion, not a high r.
- Figure: ISS-T RNA effect (x) vs OSD-462 protein effect (y) scatter, target sets highlighted, with the null band shaded.

### Layer 2 — Phosphoprotein activity (`scripts/osd462/02_phospho_axis.py`) — conditional

- **Hard checkpoint first.** On loading the phospho workbook, grep the phosphosite table for Slc12a3, Wnk1, Wnk4, Stk39, Oxsr1, Klhl3. Print a coverage report. **If phospho-NCC and phospho-SPAK/OSR1 sites are absent, stop the activity claim** and report "phosphoproteomics did not quantify the relevant transporter phosphosites; the activity layer is not testable in OSD-462" — this is an acceptable, honest outcome.
- If present: compute FL−GC phospho effect per site (plex-corrected), focused on the WNK-SPAK/OSR1-NCC chain. Report site-level effects with CIs; do not aggregate away the individual sites (phospho-NCC Thr is the readout that matters).
- Cross-layer integration: state whether transcript-down + protein-down + phospho-down co-occur for the DCT axis. Only if all three align may the manuscript say "evidence for reduced NCC transport-pathway activity."

### Layer 3 — Network-candidate translation check (`scripts/osd462/03_network_translation.py`)

- Take RRRM-2 LIONESS/node2vec top-composite candidates (`gene_axis_priority.tsv`).
- Map to OSD-462 proteins/phosphosites; test whether candidates are enriched among protein- or phospho-changing genes vs a matched random null (same null machinery as Layer 1).
- This is the *only* test that speaks to the network methods. Frame the outcome either way: enrichment = "network layer nominated genes with independent proteomic support"; no enrichment = "network layer did not translate to protein biology." Both are reported.

### Layer 4 — Same-modality RNA sanity check (`scripts/osd462/04_rna_recurrence.py`)

- Before trusting cross-modality results, confirm OSD-462 **RNA** flight effect recurs the RRRM-2 ISS-T matrix-high/DCT-low direction (pathway-vector cosine + leave-one-pathway-out), exactly as done for OSD-513. If OSD-462 RNA does *not* recur the signal, a protein concordance result would be hard to interpret. This is a gate, run first.

---

## 6. Statistics, null models, multiple testing

- **Primary inference = the abundance-matched random-gene-set null** (Layer 1/3). It is robust to the generically weak RNA↔protein correlation and to gene-set size.
- **Plex** is a batch factor in every proteomics/phospho effect estimate; never pool raw channels across plexes without it.
- **Multiple testing:** BH within an explicitly named family — the 4 targeted gene sets in Layer 1 form one family; phosphosites in Layer 2 another. Do not BH across heterogeneous layers.
- **Peptide-count weighting / filtering:** ≥2 peptides for inclusion; report results with and without a stricter ≥4-peptide filter as sensitivity (Wnk4, Stk39, Klhl3, Col3a1 are marginal and will move under this filter).
- **Effect-size honesty:** report concordance as effect size + null band, not just p. A statistically-significant-but-tiny concordance is reported as such.

---

## 7. Decision table — what each outcome means

| Layer 1 (protein) | Layer 2 (phospho) | Interpretation / claim allowed |
|---|---|---|
| Matrix + DCT concordant | Phospho-axis down | Strongest: transcript+protein+signaling concordant transporter remodeling. |
| Matrix + DCT concordant | Phospho null / not quantified | "Protein-level concordant remodeling; activity not assessable" — still a real upgrade. |
| Matrix concordant, DCT not | — | Matrix story protein-anchored; DCT stays transcript-level — narrow the headline to matrix. |
| Neither concordant | — | Honest null: the RNA signal does not extend to protein; reconsider the whole claim. |

The plan commits to whichever row the data land on.

---

## 8. Claim ceiling — what this can and cannot establish

- **Can:** that the RNA remodeling direction has independent protein (and possibly phospho) support in a second flight cohort; that specific transporter machinery moves; that network-nominated genes do or do not translate.
- **Cannot:** physiology. No histology, no urine/serum electrolytes in scope. Even a perfect three-layer concordance yields a **molecular** claim ("evidence for transporter-pathway remodeling"), not a functional one ("kidney electrolyte handling changed"). Histology and urine/serum chemistry remain explicitly future work.
- **Does not rescue "network rewiring" as a method.** A positive result validates the *biology* (matrix↑/DCT↓); LIONESS/node2vec stay a candidate-prioritization layer (Layer 3 is their only test).

---

## 9. Deliverables

- Scripts: `scripts/osd462/00_harmonize.py` … `04_rna_recurrence.py`.
- Module: `src/multiomics/osd462_anchor.py` (effect estimation, plex correction, matched-null sampler) + `tests/test_osd462_anchor.py`.
- Outputs under `data/results/<run>/osd462_anchor/`: `osd462_flight_effects.tsv`, `protein_concordance_summary.tsv`, `phospho_axis_summary.tsv`, `network_translation.tsv`, `osd462_rna_recurrence.tsv`, `manifest.json` (input SHAs, as in the existing pipeline).
- One figure: a 4-panel multi-omics anchor dashboard (RNA-protein scatter; DCT-axis protein bars; phospho-axis sites; network-translation enrichment).
- A short Methods + Results subsection slotted into the manuscript after the mechanism-score section.

## 10. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Phospho file lacks NCC/SPAK sites | Layer 2 checkpoint; fall back to protein-only, report honestly |
| RNA↔protein concordance weak (expected) | Matched-null test detects *relative* enrichment, not absolute r |
| TMT plex confound | Plex as batch factor everywhere |
| DCT dilution in bulk tissue | NCC well-quantified (34/36 pep); still report peptide-count sensitivity |
| Cross-strain divergence (B6129 vs B6N) | Framed as recurrence, not replication; Layer 4 RNA gate |
| Over-claiming physiology | §8 ceiling fixed in advance |

## 11. Work sequence

1. **Layer 0 + Layer 4** — harmonize, confirm OSD-462 RNA recurs the signal (gate).
2. **Layer 1** — protein concordance + matched null (runnable now with downloaded data).
3. **Layer 2** — phospho axis, once the user's download lands and the checkpoint passes.
4. **Layer 3** — network-candidate translation.
5. Figure + manuscript subsection; update claim hierarchy.
