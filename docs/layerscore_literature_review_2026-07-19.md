# LayerScore and spaceflight-kidney literature review

**Date:** 2026-07-19  
**Scope:** repository-wide review of the v12 biology manuscript, LayerScore preregistration, source/configuration files, tests, and `data/results/` outputs  
**Catalog size:** 111 distinct, DOI-bearing papers: 37 biology, 40 statistics/methodology, and 34 future-direction papers  
**DOI policy:** every DOI below was checked against Crossref metadata; ambiguous records were separately checked against the publisher or PubMed/PMC. Dataset DOIs, preprints, corrections, books, and papers without a DOI were not counted.

## Executive synthesis

The repository now contains two related but scientifically distinct projects.

1. **The v12 biology paper** argues that the reproducible bulk-RNA response to mouse kidney spaceflight is matrix/endothelial-high and DCT/NCC-WNK-low, while matched OSD-462 RNA and protein effects are weakly coupled and the most interpretable distal-nephron signal resolves at regulatory phosphorylation. The DCT2/CNT-leaning extension is a modest parent-gene prioritization from whole-kidney phosphoproteomics, not cell-resolved localization.
2. **LayerScore** is a proposed pathway-level decision procedure for locating a response across RNA, protein, and phosphoprotein layers. Its methodological identity is the distinction between a supported absence of protein carriage and insufficient evidence, implemented through observability-matched nulls, equivalence testing, and explicit abstention.

The literature therefore has to do four jobs simultaneously: establish spaceflight renal physiology; explain mouse distal-nephron, WNK-SPAK/OSR1-NCC, cell-composition, and injury/repair biology; justify every inferential component already used or planned; and identify validation systems that can turn LayerScore from a kidney-specific scoring workflow into a general method.

### Important citation corrections found during this review

- The v12 bibliography assigns da Silveira et al. the DOI `10.1016/j.cell.2020.10.005`. That DOI resolves to an unrelated SARS-CoV-2 rhesus-macaque paper. The correct DOI for *Comprehensive Multi-omics Analysis Reveals Mitochondrial Stress as a Central Biological Hub for Spaceflight Impact* is [10.1016/j.cell.2020.11.002](https://doi.org/10.1016/j.cell.2020.11.002).
- The two entries currently marked “confirm pages/doi” can be completed: Christensen et al., *Renal and sympathoadrenal responses in space*, [10.1053/ajkd.2001.27758](https://doi.org/10.1053/ajkd.2001.27758); Drummer et al., *Water and sodium balance in space*, [10.1053/ajkd.2001.27765](https://doi.org/10.1053/ajkd.2001.27765).
- Pietrzyk et al., *Renal stone formation among astronauts* (2007), and Cockett et al., *Astronautical urolithiasis* (1962), appear not to have registered article DOIs. They remain historically useful but are outside this DOI-only catalog.

## What the repository and results require the literature to cover

| Repository object or result | What the literature must support | Principal catalog blocks |
|---|---|---|
| OSD-771, OSD-513, OSD-253, OSD-102, OSD-163 RNA recurrence | Mission, strain, sex, arm, control-scenario, and tissue-state heterogeneity; contrast-level rather than raw-expression pooling | B01-B13, B30-B33, M31-M33, F05, F33 |
| Five-cohort matrix/ECM recurrence; weak transport recurrence | Fibrosis/remodeling, endothelial/stromal expansion, bulk-composition dilution, and heterogeneity-aware meta-analysis | B02, B05, B33, M18-M20, M31-M33 |
| OSD-462 RNA-protein mismatch | Baseline mRNA-protein coupling limits, assay observability, missingness, matched nulls, and direction-aware pathway effects | M34-M40, F15, F28-F29 |
| NCC/SPAK/WNK regulatory-site suppression with flat parent protein | DCT physiology, WNK-SPAK/OSR1-NCC regulation, potassium/chloride sensing, parent-adjusted phosphoregulation, and functional-site annotation | B14-B29, M28-M30, F16-F20 |
| DCT1/DCT2 and DCT2/CNT prior enrichment | Native subtype atlases, gradient rather than hard-marker biology, parent-gene dependence, abundance/detection effects, and whole-kidney localization limits | B14, B29-B33, F05, F34 |
| Composition-aware and parent-normalized sensitivities | Deconvolution, surrogate/batch effects, paired models, clustered dependence, and missing-data mechanisms | M15-M20, M39-M40 |
| Cross-OSDR contrast-vector framework | Stability gating, shrinkage covariance, bootstrap/permutation uncertainty, cosine/projection geometry, and meta-analysis | M04-M10, M14, M21-M25, M31-M33 |
| LayerScore `rna_only_supported` versus `indeterminate` | TOST equivalence, smallest effect size of interest, regression dilution/errors-in-variables, observability, and abstention | M11-M13, M37-M40 |
| LayerScore simulation and real-data validation | Known-truth simulation, large matched multi-omics benchmarks, positive-control perturbations, tissue specificity, and causal follow-up | F01-F34 |

## Statistical structure observed in the repository

This map is descriptive of the implemented and preregistered analysis, not a claim that every branch has passed its gate.

| Analysis layer | Implemented or planned procedure | Literature anchors |
|---|---|---|
| RNA preprocessing | DESeq2/VST, limma/edgeR resources, within-study standardization, SVA and technical-covariate residualization | M01-M03, M15-M17 |
| Composition | MuSiC plus orthogonal deconvolution/NNLS; agreement gating; residualized pathway effects | M18-M20 |
| Network construction | Ledoit-Wolf shrinkage, frozen top-k backbone, LIONESS sample-specific edges, WGCNA modules | M14, M21-M22 |
| Exploratory geometry | node2vec/PecanPy embeddings, multi-seed stability, generalized Procrustes alignment | M23-M25 |
| Contrast geometry | Aging and flight vectors; weighted/unweighted projection slope, cosine, residual fraction, and resolution hierarchy | M09-M10, M13-M14, M39 |
| Uncertainty | Stratified full-pipeline bootstrap (`B=2000`), stratified permutation (`K=5000`), plus-one empirical p-values, stability gates | M08-M10 |
| Multiplicity | BH/BY FDR; selective-family Benjamini-Bogomolov logic for predeclared families | M04-M07 |
| Gene sets/regulators | Standardized-mean panels, GSEA/camera context, Fisher/logistic enrichment, KSEA, decoupleR/PROGENy | M26-M30 |
| Cross-cohort synthesis | Within-study effects, DerSimonian-Laird random effects, I-squared, signed Stouffer, leave-one-out checks | M31-M33 |
| LayerScore RNA-protein estimands | Sign-oriented carriage and through-origin slope, matched on abundance, peptide count, missingness, and RNA-effect precision | M34-M40 |
| Supported absence | TOST equivalence inside a preregistered assay-scale margin; no “non-significant = absent” shortcut | M11-M12 |
| Phosphoregulation | Paired phosphosite-minus-parent-protein modeling; parent-observable scope; functional-site annotation | M28-M30, F16-F20 |

## Results-folder audit and claim boundaries

The audit covered 4,481 nonempty result files, including 2,155 TSV tables, 1,431 NumPy arrays, 407 PNG figures, 191 JSON records, and 69 PDFs. The following outputs most directly determine the literature priorities and the permissible claim language.

| Result branch | Result observed | Consequence for interpretation and literature |
|---|---|---|
| Cross-OSDR recurrence | Across five cohorts, the 17-gene matrix/ECM set has signed Stouffer p = 0.000703 and median I-squared = 63.8%; the 15-gene DCT/NCC/WNK set has p = 0.187, median I-squared = 73.3%, and unstable leave-one-out sign. | The matrix program is the stronger recurrent RNA result. The transport result is mechanistically important in OSD-462 but should not be described as a uniform pan-mission RNA signature. See `cross_osdr_recurrence/recurrence_meta_verdict.json`; M31-M33. |
| Cross-layer propagation | Among 11 predeclared pathways with 10,000 matched null draws, no pathway receives a calibrated forward RNA-to-protein label. TLR4/innate and ECM are labeled `protein_inverted_calibrated`; DCT/NCC/WNK, broad tubular transport, and S1P are RNA-to-phospho candidates. | This is evidence for layer specificity, not a general claim that flight RNA never propagates to protein. Observability, small pathway counts, and protein missingness are central. See `run_20260606_v11_layer_specificity/propagation/rna_protein_propagation_summary.tsv`; M34-M40. |
| NCC regulatory phosphosites | Several Slc12a3 sites are strongly suppressed (site effects approximately -0.79 to -0.93) while the mapped parent-protein effect is +0.089. | The clearest distal-nephron evidence is regulatory phosphorylation rather than total NCC abundance. Parent subtraction strengthens the contrast but does not estimate occupancy. See `h2_occupancy/h2_occupancy_site_effects.tsv`; B18-B24 and F16-F19. |
| Parent-normalized subtype-prior enrichment | DCT1 top-decile enrichment: OR = 1.522, q = 1.04e-11. DCT2-bottom comparator: OR = 1.324, q = 4.27e-6. | The result supports a distal-nephron prior but not DCT1 exclusivity. Both comparator directions must remain visible. See `h2_occupancy/h2_occupancy_verdict.json`; B14 and B29-B33. |
| Composition-aware sensitivity | Adjusted top-decile enrichment passes, while the full continuous and site-fixed interaction tests fail at FDR 0.10. | Claim a robust enriched subset, not a broad continuous DCT1 gradient. Animal-level bulk composition scores may also lie on the causal path. See `h2_composition_adjusted/h2_composition_aware_verdict_single.json`; M18-M20 and F30-F31. |
| DCT2/CNT-leaning score definitions | Mean-difference OR = 1.769, permutation p = 0.0102; detection-aware OR = 1.847, p = 0.0030; log2-ratio OR = 1.328, p = 0.0624; rank-average OR = 1.000, p = 0.606. | Enrichment is modest and score-definition-sensitive. It is a parent-gene prioritization result, not cell-resolved localization. See `run_20260709_v11_dct2_specificity_enrichment/dct2_score_definition_enrichment.json`; B14, B29-B33, F06-F14. |
| External spatial reference | The best Visium timepoint match is week 6 IRI, cosine = 0.801 across 11 pathways. Xenium supplies 1,374,915 cells but only a 299-gene targeted panel. | This is injury/repair contextualization, not spatial validation of a spaceflight lesion. See `spatial_reference/spatial_reference_projection_verdict.json`; B33 and F06-F12. |

## Category 1 — Biology (37 papers)

### A. Direct spaceflight-kidney, astronaut physiology, and system-wide spaceflight biology

**B01. Siew et al., *Cosmic kidney disease: an integrated pan-omic, physiological and morphological study into spaceflight-induced renal dysfunction*.** The indispensable comparator and external anchor for OSD-462 transporter dephosphorylation, nephron remodeling, stone risk, and radiation injury. [DOI](https://doi.org/10.1038/s41467-024-49212-1)

**B02. Finch et al., *Spaceflight causes strain-dependent gene expression changes in the kidneys of mice*.** Directly supports strain-dependent lipid, ECM, and TGF-beta responses and the need to treat OSD-102/OSD-163-style cohorts as heterogeneous replications. [DOI](https://doi.org/10.1038/s41526-025-00465-0)

**B03. Chouker et al., *Impact of spaceflight on endocrine, metabolic and kidney function: current evidence, open issues, and potential countermeasures*.** Current integrative review for RAAS, metabolism, kidney physiology, countermeasures, and knowledge gaps. [DOI](https://doi.org/10.1186/s12915-025-02471-w)

**B04. Olde Engberink et al., *The kidney, volume homeostasis and osmoregulation in space: current perspective and knowledge gaps*.** Best broad physiological source for fluid redistribution, sodium/water handling, osmoregulation, and renal evidence boundaries. [DOI](https://doi.org/10.1038/s41526-023-00268-1)

**B05. Hammond et al., *Effects of Space Flight on Mouse Liver versus Kidney: Gene Pathway Analyses*.** Early direct mouse-kidney transcriptome study and a useful same-animal tissue-comparison precedent for LayerScore H2. [DOI](https://doi.org/10.3390/ijms19124106)

**B06. Hayashi et al., *Bone mineral loss damages renal tubules in mice*.** Links actual microgravity, bone resorption, FGF23/phosphaturia, and tubular injury, adding a bone-kidney mechanism complementary to the primary renal-transport interpretation. [DOI](https://doi.org/10.1038/s42003-026-09603-0)

**B07. Liakopoulos et al., *The kidney in space*.** Foundational review of renal structural and functional effects of weightlessness; useful for historical framing but should be paired with newer data. [DOI](https://doi.org/10.1007/s11255-012-0289-7)

**B08. Drummer et al., *Water and sodium balance in space*.** Primary source for diminished intake/output and altered sodium-water balance during flight. [DOI](https://doi.org/10.1053/ajkd.2001.27765)

**B09. Grigoriev et al., *Water and electrolyte studies during long-term missions onboard the space stations SALYUT and MIR*.** Long-duration cosmonaut evidence for endocrine, electrolyte, and renal adaptation and postflight readaptation. [DOI](https://doi.org/10.1007/BF00189308)

**B10. Smith et al., *Men and Women in Space: Bone Loss and Kidney Stone Risk After Long-Duration Spaceflight*.** Sex-aware human evidence connecting skeletal demineralization with urinary stone risk. [DOI](https://doi.org/10.1002/jbmr.2185)

**B11. Garrett-Bakelman et al., *The NASA Twins Study*.** Human multi-omic and urinary/physiological context for the repository’s exploratory human-concordance work; n-of-one limitations must remain explicit. [DOI](https://doi.org/10.1126/science.aau8650)

**B12. Afshinnekoo et al., *Fundamental Biological Features of Spaceflight*.** Cross-system synthesis establishing shared stress, immune, mitochondrial, and molecular features of spaceflight. [DOI](https://doi.org/10.1016/j.cell.2020.10.050)

**B13. da Silveira et al., *Comprehensive Multi-omics Analysis Reveals Mitochondrial Stress as a Central Biological Hub for Spaceflight Impact*.** Supports the oxidative-stress/mitochondrial context and demonstrates cross-tissue, cross-species multi-omic integration. [DOI](https://doi.org/10.1016/j.cell.2020.11.002)

**B34. Christensen et al., *Renal and sympathoadrenal responses in space*.** Primary evidence for altered natriuretic responses, renin, sympathetic activity, and renal fluid handling in flight. [DOI](https://doi.org/10.1053/ajkd.2001.27758)

**B35. Drummer et al., *Vasopressin, Hypercalciuria and Aquaporin - The Key Elements for Impaired Renal Water Handling in Astronauts?*.** Mechanistic bridge between water handling, AQP2, vasopressin, and hypercalciuria. [DOI](https://doi.org/10.1159/000064111)

**B36. Pastushkova et al., *Detection of Renal Tissue and Urinary Tract Proteins in the Human Urine after Space Flight*.** Direct human urinary-proteome evidence relevant to translational biomarkers and the v12 urine-concordance discussion. [DOI](https://doi.org/10.1371/journal.pone.0071652)

**B37. Whitson et al., *Effect of Potassium Citrate Therapy on the Risk of Renal Stone Formation During Spaceflight*.** Countermeasure evidence and a reminder that biomolecular mechanisms should be related to operationally testable renal-risk endpoints. [DOI](https://doi.org/10.1016/j.juro.2009.07.010)

### B. Distal-nephron, DCT/CNT, WNK-SPAK/OSR1-NCC, aldosterone, and vasopressin biology

**B14. Su et al., *Enriched Single-Nucleus RNA-Sequencing Reveals Unique Attributes of Distal Convoluted Tubule Cells*.** The direct source for the native DCT1/DCT2 reference prior used in v12; crucial for interpreting subtype gradients and marker limitations. [DOI](https://doi.org/10.1681/ASN.0000000000000297)

**B15. Subramanya and Ellison, *Distal Convoluted Tubule*.** Compact authoritative review of DCT transport, regulation, and disease; a central v12 background citation. [DOI](https://doi.org/10.2215/CJN.05920613)

**B16. McCormick and Ellison, *Distal Convoluted Tubule*.** Deeper physiology review covering DCT segmentation, electrolyte control, and WNK-NCC regulation. [DOI](https://doi.org/10.1002/cphy.c140002)

**B17. Gamba, *Molecular Physiology and Pathophysiology of Electroneutral Cation-Chloride Cotransporters*.** Foundational transporter biology for NCC/NKCC/KCC structure, function, and pathophysiology. [DOI](https://doi.org/10.1152/physrev.00011.2004)

**B18. Moriguchi et al., *WNK1 Regulates Phosphorylation of Cation-Chloride-Coupled Cotransporters via SPAK and OSR1*.** Establishes the core WNK-SPAK/OSR1-to-cotransporter phosphorylation cascade. [DOI](https://doi.org/10.1074/jbc.M510042200)

**B19. Vitari et al., *The WNK1 and WNK4 Protein Kinases ... Phosphorylate and Activate SPAK and OSR1*.** Independent mechanistic definition of the same kinase cascade and Gordon-syndrome context. [DOI](https://doi.org/10.1042/BJ20051180)

**B20. Richardson and Alessi, *The Regulation of Salt Transport and Blood Pressure by the WNK-SPAK/OSR1 Signalling Pathway*.** Integrates the kinase cascade with blood-pressure and salt-transport physiology. [DOI](https://doi.org/10.1242/jcs.029223)

**B21. Pacheco-Alvarez et al., *The Na+:Cl- Cotransporter Is Activated and Phosphorylated ... upon Intracellular Chloride Depletion*.** Direct basis for chloride-sensitive NCC activation and site-level phosphorylation interpretation. [DOI](https://doi.org/10.1074/jbc.M603773200)

**B22. Terker et al., *Potassium Modulates Electrolyte Balance and Blood Pressure through Effects on Distal Cell Voltage and Chloride*.** Key mechanistic source for the potassium-voltage-chloride-WNK-NCC switch used in the proposed upstream model. [DOI](https://doi.org/10.1016/j.cmet.2014.12.006)

**B23. Mukherjee et al., *Roles of WNK4 and SPAK in K+-Mediated Dephosphorylation of the NaCl Cotransporter*.** Especially relevant to v12 because it focuses on rapid NCC/SPAK dephosphorylation and DCT1-specific responses. [DOI](https://doi.org/10.1152/ajprenal.00459.2020)

**B24. Thomson et al., *WNK Bodies Cluster WNK4 and SPAK/OSR1 to Promote NCC Activation in Hypokalemia*.** Provides a subcellular mechanism and a concrete imaging endpoint for future validation. [DOI](https://doi.org/10.1152/ajprenal.00232.2019)

**B25. Boyden et al., *Mutations in Kelch-Like 3 and Cullin 3 Cause Hypertension and Electrolyte Abnormalities*.** Establishes KLHL3/CUL3 as the WNK degradation machinery; supports treating the repository’s current lack of ubiquitin/coverage evidence as a genuine gap. [DOI](https://doi.org/10.1038/nature10814)

**B26. Wilson et al., *Human Hypertension Caused by Mutations in WNK Kinases*.** Foundational genetic evidence tying WNK dysregulation to salt handling and blood pressure. [DOI](https://doi.org/10.1126/science.1062844)

**B27. Arroyo et al., *Aldosterone Paradox: Differential Regulation of Ion Transport in Distal Nephron*.** Necessary for interpreting mineralocorticoid/RAAS directionality without collapsing context-dependent potassium and volume responses. [DOI](https://doi.org/10.1152/physiol.00049.2010)

**B28. Cheng et al., *A Systems Level Analysis of Vasopressin-Mediated Signaling Networks in Kidney Distal Convoluted Tubule Cells*.** Source for PXD001729/mpkDCT phosphoproteomic context and an important warning about incomplete shared-site coverage. [DOI](https://doi.org/10.1038/srep12829)

**B29. Loffing et al., *Distribution of Transcellular Calcium and Sodium Transport Pathways along Mouse Distal Nephron*.** Anatomical basis for DCT1-to-DCT2/CNT transitions and mixed transporter programs. [DOI](https://doi.org/10.1152/ajprenal.0085.2001)

### C. Mouse-kidney cell atlases, aging, and injury/repair context

**B30. Park et al., *Single-Cell Transcriptomics of the Mouse Kidney Reveals Potential Cellular Targets of Kidney Disease*.** Broad mouse-kidney cell reference for deconvolution, marker validation, and cell-of-origin caution. [DOI](https://doi.org/10.1126/science.aar2131)

**B31. Ransick et al., *Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney*.** Supports sex/segment heterogeneity and the need to distinguish proximal, loop, distal, stromal, vascular, and immune programs. [DOI](https://doi.org/10.1016/j.devcel.2019.10.005)

**B32. Tabula Muris Consortium, *A Single-Cell Transcriptomic Atlas Characterizes Ageing Tissues in the Mouse*.** Source for the external kidney aging axis and a core reference for deconvolution and age-generalization. [DOI](https://doi.org/10.1038/s41586-020-2496-1)

**B33. Xuanyuan et al., *Multimodal Spatial Transcriptomic Characterization of Mouse Kidney Injury and Repair*.** Direct source for the IRI Visium/Xenium contextual reference and the boundary between pathway resemblance and true spaceflight localization. [DOI](https://doi.org/10.1038/s41467-025-62599-9)

## Category 2 — Statistics and methodology (40 papers)

### A. Differential expression, multiple testing, resampling, equivalence, and measurement error

**M01. Love et al., *Moderated Estimation of Fold Change and Dispersion for RNA-seq Data with DESeq2*.** Supports count modeling, shrinkage, and external-cohort RNA differential expression. [DOI](https://doi.org/10.1186/s13059-014-0550-8)

**M02. Ritchie et al., *limma Powers Differential Expression Analyses for RNA-Sequencing and Microarray Studies*.** Supports moderated linear-model inference used throughout the RNA and edge-regression history. [DOI](https://doi.org/10.1093/nar/gkv007)

**M03. Robinson et al., *edgeR: A Bioconductor Package for Differential Expression Analysis of Digital Gene Expression Data*.** Alternative negative-binomial framework and a useful corroborating reference for count-based cohort analysis. [DOI](https://doi.org/10.1093/bioinformatics/btp616)

**M04. Benjamini and Hochberg, *Controlling the False Discovery Rate*.** Primary justification for within-family BH correction. [DOI](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

**M05. Benjamini and Yekutieli, *The Control of the False Discovery Rate in Multiple Testing under Dependency*.** Supports the repository’s conservative BY diagnostics where pathway membership or phosphosites are dependent. [DOI](https://doi.org/10.1214/aos/1013699998)

**M06. Storey, *A Direct Approach to False Discovery Rates*.** Useful for explaining q-values and the distinction between estimated discovery proportion and individual p-values. [DOI](https://doi.org/10.1111/1467-9868.00346)

**M07. Benjamini and Bogomolov, *Selective Inference on Multiple Families of Hypotheses*.** Direct anchor for the two-stage family/gene testing logic in `permutation_bootstrap.py`. [DOI](https://doi.org/10.1111/rssb.12028)

**M08. Phipson and Smyth, *Permutation P-values Should Never Be Zero*.** Direct support for the plus-one empirical p-value correction used in the repository and essential for the RR-1 exact-permutation floor discussion. [DOI](https://doi.org/10.2202/1544-6115.1585)

**M09. Efron, *Bootstrap Methods: Another Look at the Jackknife*.** Foundational support for resampling-based uncertainty and stability estimation. [DOI](https://doi.org/10.1214/aos/1176344552)

**M10. Efron, *Better Bootstrap Confidence Intervals*.** Supports BCa sensitivity analyses when percentile intervals show bias or skew; the preregistered percentile interval should remain primary unless amended. [DOI](https://doi.org/10.1080/01621459.1987.10478410)

**M11. Schuirmann, *A Comparison of the Two One-Sided Tests Procedure ... for Assessing Equivalence*.** Foundational TOST source for the `rna_only_supported` gate. [DOI](https://doi.org/10.1007/BF01068419)

**M12. Lakens, *Equivalence Tests: A Practical Primer*.** Clear modern source for smallest-effect-size-of-interest selection, interval interpretation, and the error of treating non-significance as evidence of absence. [DOI](https://doi.org/10.1177/1948550617697177)

**M13. Frost and Thompson, *Correcting for Regression Dilution Bias*.** Directly addresses attenuation of the LayerScore RNA-to-protein slope when RNA effects are measured with error. [DOI](https://doi.org/10.1111/1467-985X.00164)

### B. Shrinkage covariance, unwanted variation, deconvolution, and network geometry

**M14. Ledoit and Wolf, *A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices*.** Justifies shrinkage covariance/precision estimation for the frozen network backbone and noisy high-dimensional vectors. [DOI](https://doi.org/10.1016/S0047-259X(03)00096-4)

**M15. Leek and Storey, *Capturing Heterogeneity in Gene Expression Studies by Surrogate Variable Analysis*.** Supports unknown-confounder adjustment and motivates checking that residualization does not erase the biological contrast. [DOI](https://doi.org/10.1371/journal.pgen.0030161)

**M16. Johnson et al., *Adjusting Batch Effects in Microarray Expression Data Using Empirical Bayes Methods*.** Foundational ComBat reference for study/batch adjustment. [DOI](https://doi.org/10.1093/biostatistics/kxj037)

**M17. Zhang et al., *ComBat-seq: Batch Effect Adjustment for RNA-seq Count Data*.** Count-aware batch adjustment; relevant if any future pipeline moves batch correction before VST. [DOI](https://doi.org/10.1093/nargab/lqaa078)

**M18. Wang et al., *Bulk Tissue Cell Type Deconvolution with Multi-Subject Single-Cell Expression Reference (MuSiC)*.** Primary method citation for the multi-subject reference deconvolution already used/planned. [DOI](https://doi.org/10.1038/s41467-018-08023-x)

**M19. Newman et al., *Robust Enumeration of Cell Subsets from Tissue Expression Profiles (CIBERSORT)*.** Orthogonal reference for signature-based deconvolution and an alternative to MuSiC. [DOI](https://doi.org/10.1038/nmeth.3337)

**M20. Avila Cobos et al., *Benchmarking of Cell Type Deconvolution Pipelines for Transcriptomics Data*.** Supports the two-method agreement gate and emphasizes reference, marker, and mixture-dependent performance. [DOI](https://doi.org/10.1038/s41467-020-19015-1)

**M21. Langfelder and Horvath, *WGCNA: An R Package for Weighted Correlation Network Analysis*.** Method anchor for module construction and eigengene/module summaries. [DOI](https://doi.org/10.1186/1471-2105-9-559)

**M22. Kuijjer et al., *Estimating Sample-Specific Regulatory Networks*.** Primary LIONESS reference for converting an aggregate network into per-sample edge weights. [DOI](https://doi.org/10.1016/j.isci.2019.03.021)

**M23. Grover and Leskovec, *node2vec: Scalable Feature Learning for Networks*.** Supports the older exploratory embedding layer; it should not be allowed to outrank the preregistered contrast-vector analyses. [DOI](https://doi.org/10.1145/2939672.2939754)

**M24. Liu and Krishnan, *PecanPy: A Fast, Efficient and Parallelized Python Implementation of node2vec*.** Software/method reference for the repository’s actual node2vec engine and multi-seed sensitivity. [DOI](https://doi.org/10.1093/bioinformatics/btab202)

**M25. Gower, *Generalized Procrustes Analysis*.** Mathematical anchor for alignment of network embeddings across conditions or random seeds. [DOI](https://doi.org/10.1007/BF02291478)

### C. Gene-set and regulator inference, cross-cohort pooling, and heterogeneity

**M26. Subramanian et al., *Gene Set Enrichment Analysis*.** Foundational source for ranked gene-set interpretation and external recurrence analyses. [DOI](https://doi.org/10.1073/pnas.0506580102)

**M27. Wu and Smyth, *Camera: A Competitive Gene Set Test Accounting for Inter-Gene Correlation*.** Important corrective to naïve gene-set testing when members are correlated; useful for LayerScore simulation benchmarks. [DOI](https://doi.org/10.1093/nar/gks461)

**M28. Casado et al., *Kinase-Substrate Enrichment Analysis*.** Primary KSEA source for inferring kinase output from phosphosite effects. [DOI](https://doi.org/10.1126/scisignal.2003573)

**M29. Badia-i-Mompel et al., *decoupleR*.** Supports consensus regulator/activity inference and provides a benchmark for LayerScore’s pathway-level labels. [DOI](https://doi.org/10.1093/bioadv/vbac016)

**M30. Schubert et al., *Perturbation-Response Genes Reveal Signaling Footprints in Cancer Gene Expression (PROGENy)*.** Shows how perturbation-derived footprints can outperform membership-only pathways; useful for future positive-control benchmarks. [DOI](https://doi.org/10.1038/s41467-017-02391-6)

**M31. DerSimonian and Laird, *Meta-Analysis in Clinical Trials*.** Foundational random-effects estimator used by the recurrence module. [DOI](https://doi.org/10.1016/0197-2456(86)90046-2)

**M32. Higgins and Thompson, *Quantifying Heterogeneity in a Meta-Analysis*.** Primary I-squared reference; critical because the repository’s transport recurrence is heterogeneous even where matrix recurrence is positive. [DOI](https://doi.org/10.1002/sim.1186)

**M33. Whitlock, *Combining Probability from Independent Tests: The Weighted Z-Method Is Superior to Fisher's Approach*.** Supports signed Stouffer combination while preserving direction. [DOI](https://doi.org/10.1111/j.1420-9101.2005.00917.x)

### D. Multi-omics comparators, RNA-protein coupling, Deming regression, and observability

**M34. Paczkowska et al., *Integrative Pathway Enrichment Analysis of Multivariate Omics Data (ActivePathways)*.** Required comparator because it combines layer evidence but does not localize carriage or distinguish supported absence from uncertainty. [DOI](https://doi.org/10.1038/s41467-019-13983-9)

**M35. Meng et al., *A Multivariate Approach to the Integration of Multi-Omics Datasets (MOGSA)*.** Required comparator for per-omic contribution to an integrated gene-set score; LayerScore must demonstrate its distinct decision target. [DOI](https://doi.org/10.1186/1471-2105-15-162)

**M36. Argelaguet et al., *MOFA+: A Statistical Framework for Comprehensive Integration of Multi-Modal Data*.** Required latent-factor comparator and a useful stress test for whether LayerScore adds interpretable decisions beyond variance decomposition. [DOI](https://doi.org/10.1186/s13059-020-02015-1)

**M37. Liu et al., *On the Dependency of Cellular Protein Levels on mRNA Abundance*.** Central biological-statistical reference for why mRNA-protein concordance is neither expected to be perfect nor safely interpreted from raw correlation alone. [DOI](https://doi.org/10.1016/j.cell.2016.03.014)

**M38. Vogel and Marcotte, *Insights into the Regulation of Protein Abundance from Proteomic and Transcriptomic Analyses*.** Broad mechanistic context for transcriptional, translational, and degradation controls that generate buffering or inversion. [DOI](https://doi.org/10.1038/nrg3185)

**M39. Linnet, *Evaluation of Regression Procedures for Methods Comparison Studies*.** Direct methodological anchor for Deming/weighted Deming regression when both RNA and protein effects contain measurement error. [DOI](https://doi.org/10.1093/clinchem/39.3.424)

**M40. Lazar et al., *Accounting for the Multiple Natures of Missing Values in Label-Free Quantitative Proteomics Data Sets*.** Supports treating abundance-dependent missingness as structured observability rather than innocuous random absence. [DOI](https://doi.org/10.1021/acs.jproteome.5b00981)

## Category 3 — Papers that enable new directions (34 papers)

These are not merely additional background citations. Each one supplies an experimental platform, measurement technology, benchmark dataset, perturbational resource, or causal framework that can resolve a limitation visible in the current analyses.

### A. Kidney validation systems and biological generalization

**F01. Chapin et al., *Development of a Kidney Microphysiological System Hardware Platform for Microgravity Studies*.** The most direct experimental bridge from the computational signature to a controlled microgravity-compatible proximal-tubule system; it also provides a model for flight-certified hardware, perfusion, and repeated sampling. [DOI](https://doi.org/10.1038/s41526-024-00398-0)

**F02. Takasato et al., *Kidney Organoids from Human iPS Cells Contain Multiple Lineages and Model Human Nephrogenesis*.** Establishes a multicellular human kidney model in which coordinated RNA, protein, and phosphoprotein responses could be collected under defined perturbations. [DOI](https://doi.org/10.1038/nature15695)

**F03. Schutgens et al., *Tubuloids Derived from Human Adult Kidney and Urine for Personalized Disease Modeling*.** Offers a patient-derived, renewable tubular system and a plausible path from the repository's urinary-protein observations to individualized validation. [DOI](https://doi.org/10.1038/s41587-019-0048-8)

**F04. Wu et al., *Comparative Analysis and Refinement of Human PSC-Derived Kidney Organoid Differentiation with Single-Cell Transcriptomics*.** Documents organoid immaturity and off-target cell states; this is essential quality-control guidance before using organoids as a LayerScore ground truth. [DOI](https://doi.org/10.1016/j.stem.2018.10.010)

**F05. Doke et al., *Multi-Omic and Spatial Analysis of Mouse Kidneys Highlights Sex-Specific Differences in Gene Regulation across the Lifespan*.** Supplies the sex-by-age generalization axis needed to test whether the flight/aging alignment and distal-nephron conclusions persist beyond the present cohorts. [DOI](https://doi.org/10.1038/s41588-025-02161-x)

### B. Spatial, single-cell, proteomic, and phosphosite resolution

**F06. Stahl et al., *Visualization and Analysis of Gene Expression in Tissue Sections by Spatial Transcriptomics*.** Foundational spatial-transcriptomics framework for replacing pathway resemblance to an external injury atlas with direct anatomical localization in flown kidney. [DOI](https://doi.org/10.1126/science.aaf2403)

**F07. Rodriques et al., *Slide-seq: A Scalable Technology for Measuring Genome-Wide Expression at High Spatial Resolution*.** Provides finer spatial resolution for mapping DCT/CNT, vascular, and matrix programs without treating bulk-tissue enrichment as cell localization. [DOI](https://doi.org/10.1126/science.aaw1219)

**F08. Stickels et al., *Highly Sensitive Spatial Transcriptomics at Near-Cellular Resolution with Slide-seqV2*.** Increases sensitivity and near-cellular localization, especially relevant for sparse distal-nephron segments and low-abundance regulators. [DOI](https://doi.org/10.1038/s41587-020-0739-1)

**F09. Liu et al., *High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue*.** A direct technology precedent for measuring spatially matched molecular layers rather than inferring localization separately from bulk RNA and proteomics. [DOI](https://doi.org/10.1016/j.cell.2020.10.026)

**F10. Goltsev et al., *Deep Profiling of Mouse Splenic Architecture with CODEX Multiplexed Imaging*.** Shows how high-plex imaging can preserve cell neighborhoods; an adapted renal panel could jointly localize NCC, phosphorylated NCC/SPAK, endothelial, immune, and matrix markers. [DOI](https://doi.org/10.1016/j.cell.2018.07.010)

**F11. Giesen et al., *Highly Multiplexed Imaging of Tumor Tissues with Subcellular Resolution by Mass Cytometry*.** Imaging mass cytometry offers protein and phosphoprotein localization with morphology retained, directly addressing the current parent-protein and cell-of-origin limitations. [DOI](https://doi.org/10.1038/nmeth.2869)

**F12. Bendall et al., *Single-Cell Mass Cytometry of Differential Immune and Drug Responses across a Human Hematopoietic Continuum*.** Establishes single-cell phosphosignaling measurement and perturbation response, a conceptual template for nephron-segment-specific LayerScore validation. [DOI](https://doi.org/10.1126/science.1198704)

**F13. Specht et al., *Single-Cell Proteomic and Transcriptomic Analysis of Macrophage Heterogeneity Using SCoPE2*.** Demonstrates paired single-cell transcriptomic/proteomic reasoning that could separate true cross-layer buffering from bulk cell-composition shifts. [DOI](https://doi.org/10.1186/s13059-021-02267-5)

**F14. Schoof et al., *Quantitative Single-Cell Proteomics as a Tool to Characterize Cellular Hierarchies*.** Supports direct protein-level resolution of cell states and provides a route to test whether LayerScore labels remain stable after cell-type stratification. [DOI](https://doi.org/10.1038/s41467-021-23667-y)

**F15. Gillet et al., *Targeted Data Extraction of the MS/MS Spectra Generated by Data-Independent Acquisition*.** Foundational SWATH/DIA work; a future DIA design could reduce stochastic peptide missingness and strengthen the observability model underlying supported absence. [DOI](https://doi.org/10.1074/mcp.O111.016717)

**F16. Olsen et al., *Quantitative Phosphoproteomics Reveals Widespread Full Phosphorylation Site Occupancy during Mitosis*.** Establishes occupancy-aware phosphoproteomics and clarifies why phosphosite-minus-parent abundance is regulatory evidence but not, by itself, phosphorylation stoichiometry. [DOI](https://doi.org/10.1126/scisignal.2000475)

**F17. Johnson et al., *An Atlas of Substrate Specificities for the Human Serine/Threonine Kinome*.** Enables sequence-informed kinase attribution for regulated sites and can refine KSEA beyond incomplete literature-curated kinase-substrate lists. [DOI](https://doi.org/10.1038/s41586-022-05575-3)

**F18. Ochoa et al., *The Functional Landscape of the Human Phosphoproteome*.** Supplies functional-priority scores for phosphosites, useful for distinguishing mechanistically informative NCC/SPAK/WNK sites from merely detectable sites. [DOI](https://doi.org/10.1038/s41587-019-0344-3)

**F19. Hornbeck et al., *PhosphoSitePlus: Mutations, PTMs and Recalibrations*.** Core curated resource for site identity, regulatory function, and kinase relationships; mappings should be frozen by version to preserve reproducibility. [DOI](https://doi.org/10.1093/nar/gku1267)

**F20. Turei et al., *OmniPath: Guidelines and Gateway for Literature-Curated Signaling Pathway Resources*.** Provides signed and directed signaling priors that can test whether LayerScore gains accuracy from causal topology rather than pathway membership alone. [DOI](https://doi.org/10.1038/nmeth.4077)

### C. Perturbational positive controls and actionability

**F21. Subramanian et al., *A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles*.** Offers large-scale perturbational signatures for asking which compounds or genetic perturbations reverse, reproduce, or uncouple the repository's RNA programs. [DOI](https://doi.org/10.1016/j.cell.2017.10.049)

**F22. Abelin et al., *Reduced-Representation Phosphosignatures ... Enable Large-Scale Comparison of Drug-Induced Phenotypes*.** The P100 phosphoproteomic resource is a strong positive-control setting for cross-layer state discrimination and drug-induced phosphosignature retrieval. [DOI](https://doi.org/10.1074/mcp.M116.058354)

**F23. Dixit et al., *Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens*.** Establishes a scalable way to perturb candidate regulators and read out heterogeneous cellular consequences. [DOI](https://doi.org/10.1016/j.cell.2016.11.038)

**F24. Datlinger et al., *Pooled CRISPR Screening with Single-Cell Transcriptome Readout*.** CROP-seq supplies an independent pooled-CRISPR architecture for validating predicted regulators and interaction structure. [DOI](https://doi.org/10.1038/nmeth.4177)

**F25. Norman et al., *Exploring Genetic Interaction Manifolds Constructed from Rich Single-Cell Phenotypes*.** Shows how combinatorial perturbations reveal nonlinear interactions; useful for testing whether matrix, endothelial, and transport programs are parallel, hierarchical, or compensatory. [DOI](https://doi.org/10.1126/science.aax4438)

**F26. Replogle et al., *Mapping Information-Rich Genotype-Phenotype Landscapes with Genome-Scale Perturb-seq*.** Demonstrates genome-scale causal mapping and gives a benchmark for prioritizing a tractable regulator set from LayerScore outputs. [DOI](https://doi.org/10.1016/j.cell.2022.05.013)

**F27. Dhainaut et al., *Spatial CRISPR Genomics Identifies Regulators of the Tumor Microenvironment*.** Perturb-map is a model for combining in vivo perturbation with spatial phenotypes; the analogous kidney experiment could distinguish epithelial-autonomous regulation from stromal or vascular mediation. [DOI](https://doi.org/10.1016/j.cell.2022.02.015)

### D. Large matched multi-omics benchmarks, tissue specificity, and causal translation

**F28. Clark et al., *Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma*.** A large, kidney-derived matched RNA/protein/phosphoprotein cohort for testing LayerScore calibration, missingness gates, and label stability at a scale unavailable in OSD-462. Disease context prevents direct biological transfer, but not method benchmarking. [DOI](https://doi.org/10.1016/j.cell.2019.10.007)

**F29. Krug et al., *Proteogenomic Landscape of Breast Cancer Tumorigenesis and Targeted Therapy*.** An independent large matched multi-omics cohort for estimating false classifications, platform portability, and sensitivity to strong molecular subtypes. [DOI](https://doi.org/10.1016/j.cell.2020.10.036)

**F30. Imai et al., *A General Approach to Causal Mediation Analysis*.** Supplies formal estimands and assumptions for future tests of whether phosphoregulation mediates RNA-level or physiological responses; current observational LayerScore labels should not be described as mediation. [DOI](https://doi.org/10.1037/a0020761)

**F31. Greenland et al., *Causal Diagrams for Epidemiologic Research*.** DAGs provide a disciplined way to encode mission, sex, strain, preservation, batch, cell composition, and molecular-layer relationships before covariate adjustment. [DOI](https://doi.org/10.1097/00001648-199901000-00008)

**F32. Uhlen et al., *Tissue-Based Map of the Human Proteome*.** Supplies tissue and protein-expression priors for assessing whether a proposed translational marker is kidney-enriched and genuinely observable at the protein level. [DOI](https://doi.org/10.1126/science.1260419)

**F33. GTEx Consortium, *The GTEx Consortium Atlas of Genetic Regulatory Effects across Human Tissues*.** Provides cross-tissue expression and regulatory context for the preregistered kidney-versus-muscle specificity test and guards against treating a pan-tissue stress response as renal-specific. [DOI](https://doi.org/10.1126/science.aaz1776)

**F34. Lake et al., *An Atlas of Healthy and Injured Cell States and Niches in the Human Kidney*.** The strongest human cell-state and spatial-niche bridge for testing whether mouse matrix/endothelial and distal-nephron signals map to conserved human renal states. [DOI](https://doi.org/10.1038/s41586-023-05769-3)

## Recommended reading order

### For revising the v12 biology paper

1. Start with B01-B05 and B11-B13 for direct spaceflight and multi-omic context.
2. Read B14-B29 as the mechanistic core for the DCT/NCC-WNK interpretation.
3. Use B30-B33 to keep cell-type, age, sex, injury-atlas, and localization claims properly bounded.
4. Use M31-M33 and M37-M40 when explaining cross-cohort heterogeneity and RNA-protein mismatch.

### For developing the LayerScore methods paper

1. Start with M11-M13 and M37-M40: equivalence, measurement error, mRNA-protein coupling, and proteomic missingness define the method's central claim.
2. Compare directly against M34-M36, then benchmark gene-set behavior using M26-M30.
3. Use M04-M10 for multiplicity and resampling, and M14-M20 for covariance, nuisance structure, and composition.
4. Read F15-F20 and F28-F31 before freezing observability, phosphosite, validation, and causal-language rules.

## Concrete research directions derived from the repository

### Priority 1 — Validate the decision rule before adding more score components

Build the preregistered simulation suite with four known-truth states: RNA-to-protein carriage, protein inversion, true protein equivalence, and low-observability indeterminacy. Vary pathway size, RNA uncertainty, protein abundance, peptide count, missingness, correlated genes, and a small number of high-leverage genes. Report label-wise sensitivity, specificity, abstention rate, calibration, and confusion matrices. Then run the identical frozen code on the large matched CPTAC cohorts in F28-F29. The decisive test is whether `rna_only_supported` remains distinct from `indeterminate`, not whether LayerScore produces the most significant pathway ranking.

### Priority 2 — Resolve errors-in-variables and the equivalence margin

The through-origin protein-on-RNA slope is attenuated when the RNA contrast is noisy. Make weighted Deming regression or a hierarchical errors-in-variables model a prespecified sensitivity analysis (M13, M39), and include RNA standard error in the matched-null design. Derive the primary equivalence margin from assay-scale repeatability or negative-control null behavior, then show a complete margin-sensitivity curve (M11-M12). This is the most important statistical vulnerability in the current preregistration.

### Priority 3 — Perform the same-animal tissue-specificity experiment

Use the RR-1 kidney-versus-liver or kidney-versus-muscle comparison as a strict paired test: identical animals, shared permutations, tissue-specific observability models, and no pooling of raw expression across tissues. B05 is the direct biological precedent; F32-F33 supply broad tissue-specificity priors. Predefine whether the target is a different score magnitude, a different label, or both.

### Priority 4 — Turn the DCT2/CNT result into a spatially resolved hypothesis

The current enrichment is score-definition-sensitive, the DCT1 comparator also passes, and the source is whole kidney. A strong next experiment would combine segment markers with NCC, regulatory phospho-NCC, phospho-SPAK/OSR1, WNK bodies, endothelial markers, and matrix proteins in flown kidney. Spatial multi-omics or high-plex phosphoprotein imaging (F06-F12) can test whether suppression is DCT1, DCT2/CNT, broadly distal, or a composition artifact.

### Priority 5 — Establish phosphosite mechanism rather than parent-adjusted association

Measure targeted NCC/WNK/SPAK/OSR1 sites with parallel parent proteins and, where possible, occupancy calibration (F15-F19). Challenge kidney chips, tubuloids, or organoids with potassium/chloride shifts, aldosterone, vasopressin, WNK-SPAK inhibition, radiation, and simulated microgravity (F01-F04; B21-B28). Collect RNA, protein, phosphosite, electrolyte-transport, and injury readouts from the same experimental units.

### Priority 6 — Add perturbational positive controls and rescue hypotheses

Use L1000 and P100 (F21-F22) to select perturbations known to cause RNA-only, concordant RNA-protein, or phospho-dominant responses. These become positive controls for LayerScore's label accuracy. Then prioritize a small set of kidney-relevant regulators for Perturb-seq or spatial CRISPR validation (F23-F27), including tests of whether ECM/endothelial activation drives distal suppression or occurs in parallel.

### Priority 7 — Test demographic and human generalization

Stratify or externally validate by sex and age using F05 and B30-B32. Map the mouse signatures into human healthy/injured kidney cell states (F34) and require kidney-enrichment and protein observability checks using F32-F33 before proposing urinary or clinical biomarkers. Treat conserved direction as evidence of generalization, not proof of the same causal mechanism.

## Citation placement guide

| Draft location or claim | Most useful references |
|---|---|
| v12 Introduction: renal effects and human relevance | B01-B13, B34-B37 |
| v12 Introduction/Discussion: DCT and NCC-WNK mechanism | B14-B29 |
| v12 Methods: RNA, networks, resampling, multiplicity | M01-M10, M14-M25 |
| v12 Results/Discussion: cross-cohort recurrence | M31-M33 |
| v12 Results/Discussion: RNA-protein mismatch and observability | M34-M40 |
| v12 limitations: composition, localization, aging, sex | B30-B33, F05-F14 |
| LayerScore Introduction: methodological gap | M11-M13, M34-M40 |
| LayerScore Methods: equivalence, nulls, phosphosite rules | M08-M13, M28-M30, M39-M40, F16-F20 |
| LayerScore validation and Discussion | F01-F34 |

## Interpretation and verification notes

- The catalog is deliberately DOI-only. Historically important papers without a DOI and data repositories with dataset DOIs should be retained separately rather than forced into this article bibliography.
- A paper's presence here means it is pertinent, not that every statement in it should be adopted. Direct flight data, simulated microgravity, radiation analogs, kidney disease, organoids, and cancer proteogenomics answer different questions and must be labeled accordingly.
- The spatial IRI reference supports contextual resemblance, not localization of the spaceflight signal. The current DCT2/CNT analysis is parent-gene prioritization in bulk phosphoproteomics, not proof of a DCT2 cellular origin.
- Parent-normalized phosphosite change is not phosphorylation occupancy. Occupancy requires additional quantitative assumptions or measurements (F16).
- Random-effects summaries do not erase heterogeneity. The repository's transport recurrence is weaker and more heterogeneous than the matrix recurrence, so cohort-specific estimates and leave-one-out results should stay visible.
- Any future causal or mediation language requires an intervention or defensible identification assumptions. LayerScore is presently a descriptive cross-layer decision framework.
- DOI links and titles were verified on 2026-07-19. Before journal submission, import them into the bibliography manager and retrieve complete author, journal, year, volume, issue, page/article number, and access metadata from Crossref or PubMed rather than transcribing those fields manually.
