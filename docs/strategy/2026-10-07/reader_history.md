# Strategy reader: Project history and status of every workstream

_Read-only research agent output from workflow wf_1207e53e-984 (2026-10-06), preserved verbatim as structured notes. Input to the next-directions synthesis; not a decision document. Claims cite files/URLs the agent read; verify before acting._

## Summary

The project has run about 28 workstreams since December 2025. Only one is still alive: the cross-mission clinical-axes layer (Workflow D, config/revisions/clinical-axes.yaml). Every OSD-462 phospho line (v7–v13, DCT1/DCT2/ASDN, the 51-set phospho compartment audit) is either retracted or closed. The reason is design: in OSD-462, condition is perfectly aliased with reporter-tag block. Network rewiring, contrast vectors, LAR reversal and Grey60 are retired, each with its negative result written up.

The one live result is bounded. Across 5 terminal missions, a band of structural and compartment gene sets shifts upward. Podocyte high_specificity sits at the top (g=0.689, CI 0.042–1.336), and the structural control is in the same band (0.549; ranks 2–5 lie within 0.02). The direction is robust (G1). It cannot be separated from structural drift (K0: P1=0.112, CI −0.031 to 0.255). Within RRRM-2, persistence after live return is INCONCLUSIVE (G2). The structural shift is smaller in the live-return arm (interaction p=0.009), but recovery is aliased with shorter flight.

Two descriptive hints are more useful than any podocyte headline. In OSD-253, podocyte, structural and Grey60 scores all follow an exposure-duration pattern (day 25 near 0, day 75 near +1.3 to +1.7). In OSD-513 the shift is podocyte-specific (podocyte-minus-structural difference D, g=1.53). The blueprinted manuscript's podocyte title is superseded, so the alternative 'reproducibility' headline (§5 alt) is now the one to use.

Underused assets that matter:
- RRRM-2 heart and brain bulk RNA from the same animals (OSD-580/561/562/613), which can serve as a cross-tissue control.
- OSD-102 and OSD-163 kidney proteomics, never touched.
- OSD-457 kidney.
- The TMS atlas, already on disk.

Constraint: this environment lacks every historical run directory, the DCT, spatial and LINCS inputs, and the v13 manuscript. 11 legacy tests fail here. The README and the 08-25 audit are stale; the structural I² is 0.67%, not 67.3%.

## Items

### Founding network-rewiring pipeline: LIONESS/node2vec 'silent shifters' (Workflow A)

- **Status:** retired / negative-documented
- **Key facts:** Founding ASGSR proposal dated 9 Dec 2025 (node2vec_asgsr_osd771_kidney_transcriptomics.pdf, cited in docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md §1 and §11). Pipeline: MuSiC deconvolution → DESeq2 VST+SVA → 2,500-gene panel → Ledoit-Wolf shrinkage skeleton (top-k=80) → LIONESS → limma edge regression → PecanPy node2vec → Procrustes → cosine rewiring. A 'silent shifter' was high rewiring with |log2FC|<0.3 and DE FDR>0.2. The May 2026 remediation (docs/statistical_assessment_v5.md, docs/annotated_issues_v5.md, commit 35ab40c) removed post-hoc focused permutation and circular anchors. The leakage-safe CV (src/validation/enhanced_cv.py) came back negative on 10 May 2026, and the leave-one-cohort-out run run_20260701_sf_classifier was also negative. The v8 H4 test (network candidates translating to protein/phospho) was null. Registry: config/revisions/network-rewiring.yaml (retired; reference run data/archive/results/run_20260505_remediated_2500g, not present in this checkout). Postmortem: latex_paper/network_journey_retrospective.tex (gitignored, absent here). Reusable code: the fold-safe CV harness, plus src/networks/ (contrast_vectors, stability_test, external_aging_axis), which docs/_reorg_plan.md §3.2 marks load-bearing.
- **Open questions:** Whether the retrospective and the network-information-bound note are ever published (they were sequenced after the primary paper by docs/POST_V13_DIRECTION_DECISION_2026-08-11.md §8).
- **Relevance to next steps:** Do not revive. POST_V13 §9 retires it, and ANALYSIS_TRIAGE says any weighting w built from the tested matrix X is 'dead'. What survives is the 'transfer' template (structure defined in A, tested in B), which the atlas-defined clinical-axes sets already follow. The retrospective and a possible 'why rewiring fails at n=5' cautionary note are post-primary-paper items only.

### Original draft PPAR sex-dimorphism claim and external replication protocol

- **Status:** retired
- **Key facts:** The original draft made a sex-dependent PPAR claim: male OSD-513 enriched (q=0.079), female OSD-102 null. Sex is confounded with mission, strain, duration and processing, and no sex-by-flight interaction was estimated (docs/MANUSCRIPT_V1_V12_CROSS_VERSION_AUDIT_2026-07-29.md). The external replication protocol (docs/external_replication_protocol/: OSD-102 primary, OSD-513 sex checks, OSD-163/253 context, OSD-568 excluded in src/validation/external_replication.py) belongs to this era and is superseded.
- **Relevance to next steps:** None. Listed as retired in NEXT_PAPER_DIRECTION_AFTER_GREY60 'Directions to retire'.

### WGCNA pivot → Grey60 ECM/cell-migration module, plus its adversarial go/no-go

- **Status:** retired (NO_GO_RETIRE_STANDALONE_GREY60; registry grey60.yaml says 'blocked'). A single-study OSD-771 score survives as a bounded statement.
- **Key facts:** Early result: young q=6.1e-5, old q=0.0108, ECM OR 29.8. The Zsummary of 20.6 belongs to a 439-gene GC module; compact Grey60 has Z≈7.8. Adversarial audit of 2026-07-29 (docs/GREY60_GO_NO_GO_REPORT_2026-07-29.md):
- Gate A pass: maxT 0.00235; pooled 1.369, CI 0.956–1.760.
- Gate B fail: 0/27 flight-blind recoveries, Jaccard 0.022, projected p=0.831.
- Gate C pass: 48/48 gene contrasts positive, maximum contribution 3.05%.
- Gate D fail: baseline vs vivarium g=0.563, above the 0.50 ceiling; 89.0% of the effect retained after adjustment.
- Gate E fail: 2/4 positive; meta g=0.163, CI −0.658 to 0.984.
External effects: OSD-102 −0.342; OSD-163 0.625; OSD-253 0.526 (day 25 −0.007, day 75 1.283); OSD-462 −0.112 (UPX −0.091, polyA +0.975, total RNA −1.162); OSD-513 2.456 (live return, excluded). Rescue is retired by the owner (docs/owner_decisions.md, 2026-07-29) and by POST_V13 §9.
- **Open questions:** Whether Grey60, the structural scaffold and the podocyte high-specificity set share a single latent 'exposure-duration' axis. This has never been tested jointly.
- **Relevance to next steps:** Its external profile looks like the current structural drift: positive at OSD-253 day 75 but not day 25, very large in OSD-513, sign flips across OSD-462 preparations. It is probably the same broad ECM/structural phenomenon. Its frozen 48-gene set could serve as a pre-existing comparator in any duration or live-return analysis. It is not a paper.

### Contrast-vector (aging-vector) framework, including the external TMS aging axis

- **Status:** retired (lock 2026-05-18; never carried to a claim)
- **Key facts:** Defined in agents_instruction.md (Guardrails A–E) and config/contrast_vector_framework.yaml. RRRM-2 control aging vectors were bootstrap-unstable in 3/4 arm settings (docs/v11_novelty_extensions_implementation_plan_2026-06-06.md §6). The external TMS female aging axis was built (config/aging_reference.yaml, src/networks/external_aging_axis.py). v3 ISS-T/OSD-513 cosine was 0.641; the v4 OSD-513 fixed-reference permutation gave p=0.127. Reference run: run_20260518_172747_contrast_vectors (absent here). OSD-457 and OSD-913 were listed as 'extension' cohorts in that config and were never analyzed.
- **Relevance to next steps:** Dead as a framing. The TMS aging-axis code could be reused only if someone asks whether the structural drift looks like aging, a claim made in the RRRM-1 multi-organ OSDR descriptions (OSD-918, OSD-913). Low priority.

### Cross-cohort pathway vector / 'matrix-high, DCT-low' state and gene-wise Stouffer recurrence (v3–v4, v12 context)

- **Status:** retired as primary evidence
- **Key facts:** Matrix/ECM signed Stouffer p=7.03e-4 (median I² 63.8%); DCT/NCC-WNK p=0.187 (I² 73.3%) (docs/STATE_OF_PLAY_2026-07-29.md). Both were retired by docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md because they treated correlated genes as inferential units. The animal/mission-level fibrosis estimate is g=0.311, I² 62.5%. The Stage-0 DerSimonian–Laird anchor using config/gene_sets.yaml panels gave fibrosis g=0.798 (docs/STAGE0_GO_NO_GO_2026-07-29.md).
- **Relevance to next steps:** The gap between 0.798 and 0.311 is itself a reportable result under the reproducibility headline: estimator choice changes the field's 'recurrent fibrosis' claim (POST_V13 §5 alt).

### LAR reversal (v5) → attenuation (v6)

- **Status:** retired
- **Key facts:** The young-LAR cosine of about −0.993 was unstable because the LAR effect magnitude was near zero. In v6, young-LAR p=0.361, q=0.542, and the vector interval crosses zero (MANUSCRIPT_V1_V12 audit). Scripts: scripts/run_lar_reversal_analysis.py, scripts/plot_lar_reversal_dashboard.py.
- **Relevance to next steps:** The legitimate version of this question was asked as G2 (preregistered). The new G2 anti-correlation between ISS-T and LAR per-set effects (r=−0.65) must not be reframed as 'reversal'. That is the v5 trap: it is descriptive, unplanned in direction, and confounded with 32 vs 53–56 days of flight.

### OSD-462 multi-omics anchor era (v7–v10)

- **Status:** retired / superseded by Stage 0
- **Key facts:** Covers protein concordance, the phospho anchor, the v8 four-hypothesis record, cross-omic 'decoupling', KSEA/regulator activity and phenotype anchoring.
- v8: H1 (ECM RNA→protein) falsified; H2 (DCT protein) flat/falsified; H3 (phospho axis) 'supported' but now superseded; H4 (network candidates) null.
- KSEA (positive control on the same data, not independent): SPAK/OSR1 z=−6.31, WNK z=−4.12 (docs/regulator_activity_prioritization_plan.md).
- docs/phenotype_anchoring_plan.md scored Slc12a3 pT53/58/65/68, but positions 65 and 68 are tyrosines, so it is invalid.
- docs/mechanism_experiment_design.md (mpkDCT WNK autophosphorylation experiment) rests on Stk39 pS383 and Slc12a3 p65;68, which are co-modified or mis-indexed, so its premise is invalid.
- Decoupler PROGENy/CollecTRI per-cohort activity: src/multiomics/regulator_activity.py.
- **Relevance to next steps:** No biological value now. The per-cohort decoupler machinery could generate hypotheses about the upstream programs behind the structural drift, but only as exploratory context.

### v11 DCT1 phospho mediation, layer-specificity modules, perturbation triangulation, human concordance and CMap

- **Status:** retired / retracted (config/revisions/v11-core.yaml)
- **Key facts:** Layer-specificity modules (docs/v11_layer_specificity_execution_summary_2026-06-07.md):
- Module 2 propagation: ECM protein inversion −0.107 (q 0.091); TLR4 q 0.078.
- Module 3 observability: NCC sites at ≥99.1 percentile, but those features were later shown to be co-modified.
- Module 4 Gate 0 (OSD-105 muscle as a negative-control tissue): conditional, never run.
- Module 1 deconvolution: deferred.
Perturbation triangulation (docs/v11_perturbation_triangulation_analysis.md):
- Low-K GSE228367: cosine −0.286 (CI −0.752 to 0.458, q 0.218) for OSD-462; −0.311 for RRRM-2.
- IRI Visium GSE269622: DCT-adjacent spots −0.0436 at day 14 (p=0.005).
- PXD001729 dDAVP: 60 shared sites, 0 transport targets.
- KLHL3/CUL3: untestable.
Other pieces: Twins Table S8 'Reach D' sign test (config/human_concordance_prereg.yaml); CMap/LINCS GSE92742 (src/v11/cmap_screen.py, appendix tier). Reference run: run_20260526_v11_dct1_phospho_mediation (absent). The tests that depend on it (test_v11_h2_enrichment, test_v11_rna_protein_propagation) fail in this environment because the artifacts are missing.
- **Relevance to next steps:** Retired from all main texts (V13_INDEPENDENT_REVIEW §7). The reusable piece is the abundance/peptide-count matched-null engine (scripts/osd462/01_protein_concordance.py).

### v12 'beyond NCC/DCT1 into DCT2/CNT' headline, and the unsent Casaletto email

- **Status:** retired; the email is unsent and now actively wrong
- **Key facts:** DCT2-leaning odds ratio: 1.77 (abundance), 1.85 (detection-aware), 1.33 (specificity), 1.00 (rank average). Specificity-matched p=0.062; adjusted parent-gene logistic models null. docs/email_to_casaletto_draft.md (dated 2026-06-19) pitches this retracted DCT2/CNT and NCC result.
- **Relevance to next steps:** Do not send the email. Any outreach to Casaletto or the PI should be rewritten around the cross-mission reproducibility result.

### LayerScore method, network information bound, S1 simulation and merged methods paper

- **Status:** dropped
- **Key facts:** Covers the LayerScore preregistration (docs/layerscore_design_preregistration_2026-07-19.tex), docs/METHODS_PAPER_ANALYSIS_PLAN_2026-07-29.md and docs/PUBLICATION_STRATEGY_AFTER_GREY60_NOGO_2026-07-29.md §3. docs/PRIOR_ART_VERDICT_2026-07-29.md found them covered by QuEStVar 2024, Franks 2017, Upadhya & Ryan 2022, pathway-power 2012, curated-vs-derived 2020 and the double-dipping literature. POST_V13 §9 records: 'No S1 simulation. No LayerScore as a new method.'
- **Relevance to next steps:** Not a direction. At most a perspective piece later.

### OSD-462 Stage 0 assay-provenance audit, reporter-position diagnostic and layer-block shift (Workflow B)

- **Status:** complete (lock 2026-07-28); load-bearing only as OSD-462 eligibility justification
- **Key facts:** Findings (config/revisions/osd462-stage0.yaml, docs/STATE_OF_PLAY_2026-07-29.md):
- Protocol says TMTpro (+304.207 Da); a legacy metadata field says iTRAQ.
- Baseline sits in channels 126–128C, flight in 129N–131N, ground in 131C–133C in both plexes, with no swap.
- Zero isolated canonical NCC/SPAK phosphoforms: the T53 feature is T53/Y65 co-modified, the S383 feature is S382/S383, and positions 65/68 are tyrosines.
Diagnostics:
- Within-block reporter slopes: all 6 non-significant; pooled +0.0072 log2/step, which predicts −0.036 against observed −0.179/−0.157.
- Phospho layer shifts −0.179/−0.157 vs protein −0.021/−0.002; 3,356 of 3,869 proteins (86.7%) more negative at the phospho layer.
- Intensity gradient: Spearman −0.633 vs −0.143.
Owner decision 2026-08-11 (POST_V13 §9): no OSD-462 data note, removed from the portfolio. An OSDR curation ticket 'remains available at any time'; there is no record it was filed.
- **Open questions:** Was the curation ticket ever filed? There is no record in docs/.
- **Relevance to next steps:** Cite only in the eligibility and sensitivity subsection of the cross-mission paper. Filing the free GeneLab/OSDR curation ticket (assay label, residue indexing) is a zero-cost loose end.

### Flight-blind subtype and compartment reference build (subtype-reference) and Mouse Kidney Atlas marker tiers

- **Status:** locked; reusable infrastructure
- **Key facts:** GSE228367 paired edgeR pseudobulk, a GSE150338 check, and the Mouse Kidney Atlas (MKA). Frozen membership SHA ed64f0e2…; DCT2/CNT 27 genes, ASDN 29. The MKA was re-acquired 2026-10-05 (Zenodo 10.5281/zenodo.17395591, MD5 match). data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv (sha f080a0d9…) holds 157 podocyte high-specificity genes and 49 evaluable sets. GSE228367 and GSE150338 are not on disk in this environment.
- **Relevance to next steps:** The MKA tier code (scripts/subtype_reference/03_atlas_pseudobulk.py; build_marker_tiers in scripts/v13/compartment_adversarial_audit.py) is the backbone of the live result. It is directly reusable for the TMS second-atlas item 9 and for any new tissue or mission scoring.

### v13 exact continuous phospho inference and the v13 manuscript (Workflow C)

- **Status:** locked; claim_tier = neither; manuscript complete but its product (the OSD-462 note) removed from the portfolio
- **Key facts:** Exact enumeration of 63,504 assignments over 8,021 sites and 3,524 parent genes.
- DCT2/CNT: 5 of 27 observable (minimum 8); negative in every profile (−1.04 to −1.17).
- ASDN: 0.718, p = maxT = 0.0291; fails selectivity because podocyte scores 1.112 (p=2.36e-4, q 0.00213) and fibroblast 0.707. Observability-matched p=0.247; intensity-matched p=0.263 (docs/V13_INDEPENDENT_REVIEW_AND_PATH_FORWARD_2026-07-29.md).
latex_paper/manuscript_v13.tex (8 pp) is referenced in the README but is gitignored and absent from this checkout. The 08-25 technical audit calls it 'publication-ready', which conflicts with POST_V13 §9 (no OSD-462 note, no v14).
- **Open questions:** Whether the owner wants manuscript_v13 shelved permanently or folded in as a supplement of the cross-mission paper.
- **Relevance to next steps:** Closed. The engine (src/v13/continuous_phospho_inference.py) is reusable only for designs without aliasing. Do not cite the phospho podocyte statistic as independent confirmation.

### v13 compartment adversarial audit on phospho (51 sets) and its stop rule

- **Status:** complete; stop rule active since 2026-08-02. rrrm2 doctor reports it BLOCKED because archived inputs are missing.
- **Key facts:** Podocyte all-enriched (89 genes): 0.783, p=0.000819, FWER 0.324; 7-covariate null BH q=0.0455. Mesenchymal 0.577 (FWER 0.703); TAL 0.587 (FWER 0.683); DCT2/CNT 4 genes. No compartment passed every gate. Owner decision (docs/owner_decisions.md, 2026-08-02): no further OSD-462 compartment mining; podocyte retained only as a wet-lab validation target. The audit introduced broad_structural_scaffold_control__all, the frozen control that now drives K0.
- **Relevance to next steps:** Provides the frozen 51-set tier spec used by clinical axes. The rule that phospho may not be mined again remains in force.

### Stage-0 protocol inventory, coverage/confounding gate and axis effect-size anchor (2026-07-29)

- **Status:** complete (design calibration)
- **Key facts:** 5′/3′ gene-body coverage ratios: OSD-102 0.837, OSD-253 0.686 (SD 0.179), OSD-163 0.639, OSD-513 0.533, while RIN stays 8.3–9.4. OSD-253 flight samples have worse coverage than controls (g=−0.556, p=0.0023; flagged CONFOUNDED). OSD-513 alone drives the old distal-transport effect (−2.01). Source: docs/STAGE0_GO_NO_GO_2026-07-29.md; scripts in scripts/stage0/.
- **Open questions:** Does the 5′/3′ coverage ratio, or another quality covariate, explain part of the shared structural drift across missions?
- **Relevance to next steps:** Relevant to the new duration hint: the OSD-253 day-75 effects (podocyte 1.73, structural 1.62) come from the cohort with a measured flight-associated coverage confound. A coverage- or RIN-adjusted sensitivity of the structural drift has never been run. The LAR arm also has a RIN imbalance (g=0.93, p=0.045).

### OSD-462 matched-library robustness paper and cross-mission control-choice sensitivity (the two 07-29 priorities)

- **Status:** never run; retired
- **Key facts:** Recommended in docs/NEXT_PAPER_DIRECTION_AFTER_GREY60_2026-07-29.md and in the owner_decisions 07-29 recommendation, but never authorized. POST_V13 §9 then retired it: 'No matched-library robustness paper' (library-prep discordance is documented in PMC3747248 and PMC4620295). Control choice is partly absorbed into G1 items 5–6 (OSD-462 mRNA/UPX; OSD-253 white-light rerun) and the OSD-253 C3H sensitivity.
- **Relevance to next steps:** Do not reopen as papers. Their sensitivity outputs belong in the cross-mission paper's robustness table.

### Cross-mission clinical renal axes: 4 frozen axes, 5 terminal missions, 83 animals (Workflow D)

- **Status:** complete (retrospective lock 2026-08-11); reproduced exactly on 2026-10-06
- **Key facts:** Pooled results (primary_meta_results.tsv, verified in data/results/clinical-axes/20261006T004825Z_g1g2/):
- Glomerular barrier identity loss: −0.716 (−1.361, −0.072); prediction interval −1.455 to 0.023; maxT 0.0018.
- Tubular injury: 0.570 (−0.095, 1.234).
- Fibrosis: 0.311 (−0.732, 1.355), τ²=0.433.
- Distal transport: −0.153 (−1.319, 1.012), τ²=0.603.
Per-mission barrier (adverse orientation): 0.001, −0.981, −1.288, −0.913, −0.413. Barrier-core adjusted for the disjoint podocyte proxy: 0.507 → 0.070, with I² rising from 0 to 76%. Matched-panel p=0.0062. The CPM threshold was amended from 1.0 to 0.1 before effects were computed (label-blind). The OSD-656 urine context is non-inferential.
- **Open questions:** Whether to report permutation-calibrated intervals alongside the mHK ones, given that K0 simulations show the parametric REML/mHK test is very conservative at k=5 (0.8% rejection at nominal 5%).
- **Relevance to next steps:** This is the results backbone of the only remaining manuscript, under the 'Most reported renal responses to spaceflight do not recur across independent mouse missions' branch (POST_V13 §5 alt). Every number exists, so drafting can start now.

### Podocyte lead: 49-set family, K1–K4, strict matching, disjoint, C3H and Podxl/Nid1 variants

- **Status:** alive but demoted to a bounded observation (blueprint headline superseded by its 2026-10-06 addendum)
- **Key facts:** Podocyte high_specificity: g=0.689 (0.042, 1.336); maxT 0.0189 published / 0.0192 rerun; I²=0; prediction interval −0.053 to 1.431. Mission effects: −0.117, 1.166, 0.676, 0.541, 1.174.
Variants:
- Core-disjoint: 0.675 (FWER 0.0233).
- C3H/HeJ arm: 0.893 (FWER 0.033, I² 42.7%).
- Podxl/Nid1 forced in: 0.686. Podxl alone 0.511 (−0.335, 1.357); Nid1 alone 0.389 (I² 72.7%).
- K4 strict matching (142 targets vs 543 candidates): target-minus-matched 0.542 (−0.165, 1.250); q95 trim 0.500; matched pool is endothelium-heavy (81 of 133).
- K2: 51 papers citing Siew 2024 screened; no competing podocyte-program claim found.
K1: Siew's PODXL/NID1 call is a vote count that includes mouse kidney RNA; it is not urine-only. Per the K1 heatmap read, PODXL was up in two mouse kidney proteome columns.
No .tex exists (docs/PODOCYTE_CROSS_MISSION_PAPER_BLUEPRINT_2026-08-11.md; docs/PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md; docs/PODOCYTE_STRICT_MATCHING_AUDIT_2026-08-11.md).
- **Open questions:** Whether any measurement can separate podocyte number, sampling and per-cell state. The neighboring-glomerular-cell control cannot be evaluated with MKA.
- **Relevance to next steps:** Report it as the top of a structural band, not as the title. The 'podocyte-associated' wording is allowed only as a descriptive row.

### K0: podocyte vs structural scaffold specificity

- **Status:** complete; verdict 'shared structural drift with a podocyte-leaning tail'
- **Key facts:** From podocyte_scaffold_meta.tsv (verified):
- P0 (podocyte unadjusted): 0.198 (0.023, 0.374).
- P1 (podocyte adjusted for structural): 0.112 (−0.031, 0.255).
- S0 (structural unadjusted): 0.119.
- S1 (structural adjusted): 0.031 (−0.130, 0.192).
- D (direct difference): 0.087 (−0.043, 0.218).
P1 family-wise p (podocyte_scaffold_permutation.tsv): Freedman–Lane joint 0.0415 on the baseline map and 0.0563 on the runner-up; joint-label 0.0388 / 0.0536; legacy 0.0605 / 0.0827.
Under the reconstructed map the structural control ranks 2nd, not 4th; ranks 2–5 lie within 0.02. Its I² is 0.67%: the 08-25 technical audit misreports this as 67.3%.
- **Open questions:** Is the drift kidney-specific? No same-animal non-kidney test has ever been run.
- **Relevance to next steps:** K0 decides the title and points to the reproducibility branch. It also reframes the biological question: what drives a shared, low-heterogeneity structural/ECM/cytoskeletal shift (composition, systemic ECM response, collection/handling, RNA quality)?

### G1: preregistered sensitivities of the 157-gene podocyte set (items 1–8)

- **Status:** complete (lock e574feb; code f04f84d); ROBUST_DIRECTION and DEFINITION_ROBUST 8/8; no flips under the runner-up map
- **Key facts:** Item results:
- Median scoring: 0.631 (−0.014, 1.276); family FWER 0.045.
- Common genes: 0.693.
- Leave-one-mission: 0.550–0.846; 4 of 5 intervals cross 0.
- Leave-one-gene minimum: 0.666.
- OSD-462 mRNA / UPX: 0.751 / 0.703.
- OSD-253 rerun control: 0.802.
- OSD-163 mapping-rate adjustment: 0.669.
- Atlas leave-one-study-out (item 8): 0.588–0.752; Miao21 omitted crosses 0.
Gene influence: Mgat5 is the largest contributor at 3.6%; the top 10 genes carry 24%. Removing the top 5 contributors puts the lower bound at −0.044. P1 stays positive in every variant (0.083–0.149), and all its intervals cross 0. Source: docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md; podocyte_sensitivities/*.tsv; atlas_loso/*.tsv.
- **Relevance to next steps:** Supplies the robustness table: direction and magnitude are robust, significance margin is thin.

### G2: recovery persistence (OSD-771 ISS-T vs LAR; OSD-513 and OSD-253 duration descriptive)

- **Status:** complete; primary INCONCLUSIVE (as powered); structural INCONCLUSIVE after Holm
- **Key facts:** Primary and key-secondary results:
- Podocyte: ISS-T g 1.17, LAR g 0.14; interaction −0.281 (−0.622, 0.061); retention ratio ρ̂ 0.10 (Fieller −0.66 to 1.08).
- Structural disjoint: ISS-T 1.36, LAR −0.43; interaction −0.350 (−0.626, −0.073); FL p 0.009; maxT 0.024; Holm 0.053; DECAYS at κ ≥ 0.75.
- Power: P(PERSISTS | full persistence) = 0.117; minimum detectable interaction about 1.8 SD.
Exploratory (descriptive):
- 37 of 48 sets have a negative interaction; 6 pass maxT (4 endothelial tiers, TAL HS, mesenchymal broad).
- ISS-T vs LAR per-set effects anti-correlated, r=−0.648 (one-sided p 0.916).
- Fibrosis axis interaction FL p=0.003.
- OSD-513: podocyte g 1.56, structural −0.18, D 1.53 (maxT 0.009).
- OSD-253 (verified in persistence_duration_context.tsv): podocyte day 25 0.040 → day 75 1.731; structural −0.189 → 1.615; D 0.71 → 0.67 (flat).
Confounds: LAR flight RIN imbalance (p 0.045); winner's curse.
- **Open questions:** Is OSD-513's podocyte-specific D a male/live-return/landing effect? It is k=1, and OSD-457 is the only other live-return kidney RNA.
- **Relevance to next steps:** No recovery claim is allowed. Exposure duration and collection/landing context look like the main source of variation: OSD-253 shows a dose pattern for structural as well as podocyte, and OSD-513 is podocyte-specific. A frozen duration moderator is a more natural next question than podocyte biology.

### G1 item 9: second-atlas validation (Tabula Muris Senis)

- **Status:** deferred by owner 2026-10-06; inputs on disk
- **Key facts:** Plan: docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md. Inputs: data/external/single_cell_atlases/tms_kidney_cellxgene/ (10x object, 21,647 cells; Smart-seq2 object has no podocytes). Of 486 podocytes, 447 carry the ontology label 'plasmatocyte', so free_annotation must be used. Estimated 2–3 days. Gates are already in config/clinical_axes_podocyte_sensitivities.yaml. Un-defer if 'podocyte' stays in the title or a reviewer asks.
- **Open questions:** Whether TMS has enough podocyte donors (at least 3 with at least 25 cells) to pass the evaluability gate.
- **Relevance to next steps:** Optional under the reproducibility branch; mandatory only if 'podocyte' re-enters the title. A cheap supplementary check.

### Human urine and fluid context (OSD-656 Inspiration4 urine, NASA Twins 'Reach D', OSD-575)

- **Status:** complete and non-inferential (OSD-656); Reach D obsolete; OSD-575 never ingested
- **Key facts:** OSD-656 measures only 4 of the 34 axis genes (LCN2, CCL2, EGF, TIMP1). At R+1, LCN2 rose in 3 of 3 crew (+1.124 NPQ) and CCL2 in 3 of 3; most were lower by R+82. There are no barrier markers, no creatinine and no controls. Optional stage osd656-urine. Reach D (config/human_concordance_prereg.yaml) predicted distal/RAAS directions that are now moot. OSD-575 (Inspiration4 serum chemistry) was noted but never downloaded.
- **Relevance to next steps:** Keep out of the main paper (blueprint). There is no human validation route for a structural or podocyte program.

### Spatial and perturbation triangulation (GSE269622 Visium IRI, GSE269719 Xenium, LINCS/CMap GSE92742, PXD001729, kinome atlases)

- **Status:** retired; none of the inputs are on disk in this environment
- **Key facts:** Visium IRI was used only for a DCT-adjacent transport score (−0.0436, p=0.005). Xenium GSE269719 was never central. CMap was appendix only. The Johnson 2023 and Yaron-Barir 2024 kinome atlases fed KSEA only. All were built for the DCT/NCC hypothesis.
- **Relevance to next steps:** Low value for the current question. These are non-spaceflight datasets and cannot separate composition from per-cell state in flight tissue.

### Bulk deconvolution (MuSiC, TMS+Chen hybrid reference)

- **Status:** partial; not the basis of any claim
- **Key facts:** Three generations of sanity-check outputs and no authoritative version (technical audit §3C, §12). v11 Module 1 (formal deconvolution) was deferred. scripts/run_deconvolution.R and scripts/build_hybrid_reference.R (pinned in the frozen config).
- **Relevance to next steps:** Composition is the leading alternative explanation for the shared structural drift. A descriptive composition-adjusted re-score using MKA signatures could bound it, but it cannot resolve composition vs per-cell state.

### Cross-mission manuscript (podocyte blueprint → reproducibility branch) and the retrospective

- **Status:** blueprinted; no .tex; target bioRxiv Dec 2026 / journal Jan 2027
- **Key facts:** Two branches share their introduction, methods, eligibility gates, four-axis nulls and boundary statements (POST_V13 §5, §5 alt, §12). K0 selects the §5-alt title. Figures exist: figures/clinical_renal_axes_cross_mission/figure_1–3 (built for the podocyte framing). network_journey_retrospective is to be published only after the primary paper.
- **Open questions:** Venue choice (npj Microgravity vs a nephrology journal) is still unresolved, as is whether G2 and the duration observations go in the main text or a supplement.
- **Relevance to next steps:** This is the single remaining publication product. Drafting is the critical path, and all results are in hand. Figure 2 and 3 framing needs revising so the structural band, not the podocyte row, is shown.

### External validation and outreach: morphometry ask via PI/Siew consortium, GeneLab ticket, Melica forward search

- **Status:** planned, never done (no record in repo)
- **Key facts:** Morphometry ask (POST_V13 §7): WT1+ nuclei per glomerulus and per tuft area, plus nephrin/podocin/synaptopodin, on existing RR-10 sections (Sarder/Roufosse/Walker-Samuel). A PI one-pager was planned for weeks 2–3. The forward-citation search from Melica 2024 is still owed (technical audit §9).
- **Open questions:** Whether the PI was ever contacted.
- **Relevance to next steps:** Under K0 the ask should be broadened beyond podocyte counts to structural readouts: tuft area, cortical sampling fraction, ECM/cytoskeletal staining. It gates nothing and should go out early.

### Repository governance and verifiability (rrrm2 CLI/registry, reorg plan, freeze tagging)

- **Status:** partially done
- **Key facts:** Done: the rrrm2 CLI and registry (e8d3295), the historical archive move (d9fa6bc), and the G1/G2 preregistration committed before the run (e574feb), which is the first git-verifiable lock. The two earliest G1 stages ran about 5 minutes before the push.
Outstanding:
- latex_paper/* is still gitignored (only v11 tracked) and results directories are gitignored.
- README.md is stale: v13-centric and says '199 passed'.
- docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md is stale or wrong: it says K0 'not run', structural 'rank 4' and I² 67.3%.
- docs/owner_decisions.md has no entries after 2026-08-02. The 08-11 'no data note' decision, the K0 and G1/G2 authorizations and the item-9 deferral are recorded only in other docs.
- docs/_reorg_plan.md (2026-08-28) Batches 1–5 were only partly applied.
- Fresh-environment pytest: 11 failed / 503 passed, all from missing legacy artifacts (scratch log).
- **Relevance to next steps:** Before submission: tag the freeze, commit the run manifests, update the README and owner_decisions, and correct the audit errors. This is needed for every 'frozen before testing' sentence in the paper.

### Underused assets already on disk

- **Status:** underused
- **Key facts:** - OSD-462 proteome and phospho workbooks (data/external/osdr/OSD-462/*.xlsx): excluded by aliasing.
- OSD-462 mRNA and UPX RNA preparations: sensitivities only.
- OSD-771 BSL/VIV groups: used only for Grey60 Gate D (g=0.563); never used to test whether the structural drift tracks habitat/handling.
- OSD-771 STAR counts (data/raw/counts/).
- OSD-253 C3H arm, duration strata and rerun control: descriptive only.
- OSD-513: descriptive only.
- TMS kidney h5ad: deferred.
- MKA pseudobulk: only the podocyte and compartment tiers were used.
- data/processed/multi_study_gene_universe.tsv: 24,059 genes.
- multiqc/QC folders for every mission: coverage covariates never put into a structural-drift sensitivity.
- **Relevance to next steps:** The cheapest new analyses use these: a habitat/handling contrast (VIV vs GC vs BSL) on the structural score, a QC-covariate sensitivity, and a frozen duration moderator.

### Underused assets on OSDR, never ingested (verified via OSDR API on 2026-10-06)

- **Status:** unexplored
- **Key facts:** - RRRM-2 same-animal tissues: OSD-580 heart bulk RNA-seq (GLDS-573; 78 samples, ISS-T and LAR arms; sample names mirror OSD-771, e.g. RRRM2_HRT_BSL_ISS-T_YNG_BY1). Also OSD-561 cerebellum, OSD-562 hippocampus and OSD-613 brain bulk RNA.
- OSD-457 (JAXA MHU-3, live return after 31 days): kidney RNA-seq, WT 3 FLT + 3 GC and Nrf2-KO 3+3. Listed as an 'extension' cohort in config/contrast_vector_framework.yaml and never analyzed.
- OSD-913: kidney miRNA-seq (Grandke et al., 'MiRNAs shape … ECM and developmental pathways').
- OSD-918 and its series: RRRM-1 single-cell data from 28 organs, whose description states spaceflight drives systemic cytoskeleton/ECM remodeling. No public kidney single-cell in that series.
- OSD-102 kidney TMT proteomics (GLDS-102_proteomics_TMT_KDNL.processed.tar.gz, about 3 GB) and OSD-163 kidney proteomics (GLDS-163_proteomics_RR3_KDN_processed.zip, 20.9 MB). Never referenced in the repo; likely Siew's 'mouse kidney proteome' columns.
- OSD-102 and OSD-163 bisulfite methylation.
- OSD-105 muscle multi-omics (v11 Gate 0 candidate).
- **Open questions:** The OSD-102 and OSD-163 proteomics designs (TMT channel layout, possible aliasing) are unverified. The OSD-457 euthanasia timing is still to be confirmed from ISA. It is unknown whether the OSD-580 library prep matches OSD-771.
- **Relevance to next steps:** Three of these bear directly on what K0 left open:
1. OSD-580 and the brain sets give a same-animal cross-tissue negative control. Is the ISS-T structural drift kidney-specific or systemic/handling-related?
2. OSD-102 and OSD-163 proteomes would test whether the podocyte and structural programs move at protein level in independent missions. They need a Stage-0-style design audit first.
3. OSD-457 and OSD-513 give a tiny live-return descriptive replicate.
OSD-918 and OSD-913 are prior art that any 'systemic structural drift' framing must cite.

## Risks and constraints

- Fresh environment: data/archive/results/* (all historical locked runs: network, contrast vectors, v11, Stage 0, v13 exact, Grey60, 2026-08-11 clinical axes), GSE228367/GSE150338, spatial, LINCS, kinome and Twins inputs, and latex_paper/manuscript_v13 and the retrospective are absent. Only OSDR bulk RNA, OSD-462 workbooks, MKA, TMS and OSD-656 are on disk. Grey60, v13-compartment and v11 cannot be re-run here, and 11 legacy tests fail.
- The original id_map is lost. The baseline map was reverse-engineered to reproduce published numbers (agreement on two rows is partly by construction). The K0 P1 family-wise p<0.05 holds only on the baseline map (0.0415 vs 0.0563 on the runner-up).
- Locks of 2026-08-02 and 2026-08-11 are retrospective. The podocyte row was found by scanning a 49-set family after the targeted axes failed. G1/G2 are preregistered but post hoc relative to discovery, and two G1 stages ran about 5 minutes before the lock commit was pushed.
- Small n throughout: k=5 missions, 12–20 animals per mission, 5 per OSD-771 cell. Prediction intervals cross zero. The G2 minimum detectable interaction is about 1.8 SD. Parametric REML/mHK is very conservative at k=5 (0.8% rejection at nominal 5%). G2 carries winner's curse.
- The structural-drift pattern (many sets, τ²=0, narrow band, duration dose in OSD-253, prep sign flips in OSD-462) is the signature of a global compositional, normalization, RNA-quality or collection effect. OSD-253 also has a flight-associated 5′/3′ coverage confound (g=−0.556), and the LAR flight group has a RIN imbalance (p=0.045).
- Recovery is aliased with flight length (32 vs 53–56 days) and with on-Earth vs on-orbit euthanasia. OSD-513 is k=1, male, a separate study, and its ground controls were sampled 3 days apart.
- Prior-art risk for any 'systemic structural/ECM drift' framing: the RRRM-1 multi-organ single-cell study (OSD-918 series) and Grandke et al. (OSD-913 kidney miRNA) report systemic ECM/cytoskeleton remodeling. A forward search is needed, and the Melica 2024 forward search is still owed.
- Conflicting recommendations across docs: the 07-29 documents push an OSD-462 data note, a matched-library paper and a merged methods paper. POST_V13 §9 (2026-08-11, owner decision) retires all three. docs/owner_decisions.md does not record the 08-11, K0, G1/G2 or item-9 decisions.
- Stale or incorrect documentation: README (v13 headline, '199 passed'); the technical audit (structural control I² 67.3%, actually 0.67%; K0 'not run'; structural 'rank 4', now 2nd). docs/email_to_casaletto_draft.md pitches retracted results.
- Political: the direction difference with the Siew consortium on PODXL/NID1 must be framed as a different estimand, not a correction (K1).
- Bulk RNA cannot separate cell representation, cortical sampling and per-cell transcription. Resolving it needs external tissue or morphometry that the project does not control.

## Sources read

- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/_schema.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/network-rewiring.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/contrast-vectors.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/grey60.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/osd462-stage0.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/subtype-reference.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/v11-core.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/v13-compartment.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/v13-phospho.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/revisions/clinical-axes.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/contrast_vector_framework.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/aging_reference.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/human_concordance_prereg.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/human_urine_marker_panel.yaml
- /home/user/RRRM2_Kidney_Transcriptome/docs/STATE_OF_PLAY_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/owner_decisions.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/_reorg_plan.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/plan_changelog.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PUBLICATION_STRATEGY_AFTER_GREY60_NOGO_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/NEXT_PAPER_DIRECTION_AFTER_GREY60_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/POST_V13_DIRECTION_DECISION_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_AXES_G1_G2_PREREGISTRATION_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PODOCYTE_CROSS_MISSION_PAPER_BLUEPRINT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/MANUSCRIPT_V1_V12_CROSS_VERSION_AUDIT_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/ANALYSIS_TRIAGE_AND_NETWORK_REFRAME_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PRIOR_ART_VERDICT_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/METHODS_PAPER_ANALYSIS_PLAN_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/STAGE0_GO_NO_GO_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/V13_INDEPENDENT_REVIEW_AND_PATH_FORWARD_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/GREY60_GO_NO_GO_REPORT_2026-07-29.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/phenotype_anchoring_plan.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/regulator_activity_prioritization_plan.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/mechanism_experiment_design.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/v11_perturbation_triangulation_analysis.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/v11_novelty_extensions_implementation_plan_2026-06-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/v11_layer_specificity_execution_summary_2026-06-07.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/v11_gate0_negative_control_search_2026-06-07.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/v11_execution_results.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/osd462_multiomics_analysis_plan.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/layerscore_design_preregistration_2026-07-19.tex
- /home/user/RRRM2_Kidney_Transcriptome/docs/email_to_casaletto_draft.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/external_replication_protocol/README.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/external_replication_protocol/hypothesis_registry.tsv
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLI.md
- /home/user/RRRM2_Kidney_Transcriptome/README.md
- /home/user/RRRM2_Kidney_Transcriptome/agents_instruction.md
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/primary_meta_results.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/compartment_context_meta_results.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/podocyte_scaffold_meta.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/podocyte_scaffold_permutation.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_contrasts.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_duration_context.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_pattern_concordance.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/processed/multi_study_gene_universe.tsv.manifest.json
- /home/user/RRRM2_Kidney_Transcriptome/data/raw/metadata/s_OSD-771.txt
- https://osdr.nasa.gov/osdr/data/search?term=kidney&from=0&size=100&type=cgene
- https://osdr.nasa.gov/osdr/data/search?term=RRRM-2&from=0&size=100&type=cgene
- https://osdr.nasa.gov/osdr/data/search?term=%22age-stages%22&from=0&size=60&type=cgene
- https://osdr.nasa.gov/osdr/data/osd/files/457
- https://osdr.nasa.gov/osdr/data/osd/files/580
- https://osdr.nasa.gov/osdr/data/osd/files/102
- https://osdr.nasa.gov/osdr/data/osd/files/163
- https://osdr.nasa.gov/geode-py/ws/studies/OSD-457/download?source=datamanager&file=GLDS-457_rna_seq_bulkRNASeq_v2_runsheet.csv
- https://osdr.nasa.gov/geode-py/ws/studies/OSD-580/download?source=datamanager&file=GLDS-573_rna_seq_SampleTable_GLbulkRNAseq.csv
