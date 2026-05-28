# V11 Reviewer-Ready Revision Packet

## A. Title Options And Chosen Title

Title options considered:

1. Cross-omic decoupling and DCT1-high parent-gene phosphosite enrichment in mouse spaceflight kidney
2. Recurrent RNA remodeling and DCT1-high phosphosite subset enrichment in mouse spaceflight kidney
3. DCT1-high phosphosite subset enrichment within recurrent mouse spaceflight kidney remodeling
4. RNA remodeling and regulatory phosphosite suppression in mouse spaceflight kidney

Chosen title:

**Cross-omic decoupling and DCT1-high parent-gene phosphosite enrichment in mouse spaceflight kidney**

Rationale: this states the two central results and makes clear that the DCT1 result is parent-gene enrichment in whole-kidney phosphoproteomics.

## B. Revised Abstract

The revised abstract is now structured in four moves:

1. Prior work established NCC/SPAK/WNK dephosphorylation; this study asks what public cross-cohort and multi-omic data add around that phenotype.
2. Cross-cohort RNA recurrence and matched OSD-462 protein-abundance decoupling.
3. OSD-462 phosphoproteomic anchor plus DCT1 top-decile parent-gene enrichment, including parent-protein and compartment-score adjustment.
4. Exploratory mediation and spatial-reference analyses are hypothesis-generating, not causal or localizing.

Key numbers preserved:

- OSD-462/RRRM-2 RNA pathway-vector cosine 0.87, 95% CI 0.65-0.90.
- KSEA SPAK/OSR1 z = -6.31 and WNK z = -4.12.
- DCT1 top-decile enrichment OR 1.51, q = 1.13e-11.
- Full parent-protein and compartment-score model top-decile enrichment OR 1.30, q = 0.00158.

## C. Revised Main-Text Outline

1. Introduction
2. Methods
   - Datasets and study roles
   - RNA preprocessing and gene-set scoring
   - Cross-cohort recurrence analysis
   - Matched OSD-462 multi-omic analysis
   - DCT1/DCT2 reference prior and phosphosite enrichment
   - Parent-protein and composition-aware sensitivity models
   - Kinase activity and regulator analyses
   - Exploratory mediation and spatial-reference analyses
   - Statistics, multiple testing, and reproducibility
3. Results
   - A recurrent matrix/endothelial-high and DCT-low RNA response appears across mouse kidney spaceflight cohorts
   - OSD-462 shows RNA-protein decoupling
   - DCT identity and compartment scores suggest that the bulk DCT-low signal partly reflects remodeling context
   - Flight-suppressed phosphosites are enriched among DCT1-high parent genes
   - DCT1 top-decile enrichment persists after adjustment
   - Exploratory mediation links remodeling scores to NCC regulatory phosphorylation but cannot prove causality
   - Spatial, network, live-return, TLR4, and aging analyses provide context and boundaries
4. Discussion
   - Main finding
   - Relationship to prior work
   - Meaning of the DCT1-high parent-gene result
   - Remodeling versus transporter suppression
   - Context and boundary analyses
   - Limitations
   - Future experiments
5. Conclusion

## D-G. Section Revisions

Introduction:

- Removed internal v10/v11 framing.
- Removed framework-heavy prose.
- Established the problem in normal scientific terms: prior phenotype, open context question, bulk RNA ambiguity, matched multi-omic testing, and DCT1-high parent-gene enrichment.

Methods:

- Consolidated into nine conventional subsections.
- Removed state-space equations and de-branded ECM/DCT visualization as ordinary gene-set scoring.
- Replaced “pre-registered” with “specified before computing the protein and phosphoprotein flight effects.”
- Described KSEA as a pathway-level coherence check, not an independent discovery.

Results:

- Reordered around the main biological spine.
- Moved network/LAR/TLR4/aging/spatial inventories into compressed boundary/context language.
- Preserved all key negative results.

Discussion and Limitations:

- Rewritten around reviewer concerns: prior work, DCT1 reference-prior meaning, composition caveats, mediation limits, spatial-reference limits, and future experiments.

## H. Renamed And De-Branded Terms

| Previous wording | Revised wording |
|---|---|
| Tubulointerstitial state-space analysis | ECM and DCT transport score visualization |
| State score Ri | ECM-minus-DCT score or remodeling score |
| Recurrent remodeling-axis construction | Cross-cohort recurrent RNA score |
| Mechanism scoring | Curated gene-set scoring |
| Evidence ladder | Claim hierarchy or evidence hierarchy |
| Phenotype anchoring | Comparison to measured NCC regulatory phosphorylation |
| Causal-skeleton framework language | Exploratory mediation analysis |
| Live-return attenuation/divergence analysis | ISS-terminal versus live-return contrast |
| Gate PASS | Met the recurrence criterion, or report the cosine/CI directly |
| Tier-1 phenotype | Measured phosphoproteomic endpoint |
| DCT1-anchored mechanism | DCT1-high parent-gene subset enrichment |

## I. Claims Audit

| Claim | Evidence | Allowed wording | Forbidden wording | Placement |
|---|---|---|---|---|
| Recurrent RNA remodeling | RRRM-2, OSD-513, OSD-462 RNA; OSD-462/RRRM-2 cosine 0.87 | Recurrent matrix/endothelial-high and DCT-low RNA response | Universal spaceflight kidney signature | Main |
| RNA-protein decoupling | OSD-462 matched RNA/protein; protein concordance null | RNA context does not translate cleanly to protein abundance | RNA proves matrix protein deposition or transporter loss | Main |
| NCC/SPAK/WNK phospho suppression | OSD-462 phosphosites; total NCC flat; KSEA negative | Transporter lesion resolves at regulatory-phosphorylation layer | Newly discovered NCC dephosphorylation | Main |
| DCT1 enrichment | GSE228367 prior mapped to OSD-462 whole-kidney phosphosite parent genes | DCT1-high parent-gene subset enrichment | DCT1-specific phosphoproteomics; DCT1 cell of origin | Main |
| Composition-aware robustness | Full model OR 1.30, q=0.00158 | Persists after parent-protein and bulk-compartment adjustment | Deconfounded; independent of composition; deconvolved | Main |
| Continuous DCT1 model | Continuous/interaction models weak | No broad continuous DCT1-gradient support | Smooth DCT1-reference suppression gradient | Main |
| Mediation | Exploratory matched OSD-462 paths | Consistent with remodeling-linked hypothesis | Endothelial remodeling mediates or drives NCC dephosphorylation | Main short; details supplement |
| Spatial reference | GSE269622/GSE269719 IRI atlas | External spatial contextualization | Spatial validation/localization of spaceflight lesion | Main short; details supplement |
| Network candidates | No OSD-462 protein/phosphoprotein enrichment | Exploratory and negative boundary | Network-derived mechanism | Supplement |
| LAR/TLR4/aging | Secondary analyses | Context and boundaries | Core mechanism | Supplement or brief main |
| PXD001729 | DCT-lineage phosphoproteome; anti-alignment null | Plausibility and coverage boundary | Vasopressin anti-alignment mechanism | Supplement or brief main |
| KLHL3/CUL3 | No ubiquitinomics/coverage | Future mechanism not resolved here | KLHL3/CUL3 turnover shown | Discussion/future |

## J. Analyses/Figures Moved Or Demoted To Supplement

Detailed or secondary material now belongs in the companion results compendium/supplement:

- Full WGCNA module tables and LIONESS/node2vec candidate lists.
- Detailed LAR vector geometry and age-stratified reversal diagnostics.
- OSD-253 TLR4/control-scenario tables.
- External aging-axis projection details.
- Full spatial Visium timepoint/niche tables and Xenium annotation inventory.
- PXD001729 shared-site tables and anti-alignment diagnostics.
- KLHL3/CUL3 coverage boundary checks.
- Artifact/run manifest tables.

## K. Final Checklist

- No causal overclaim: confirmed.
- No DCT1-specific overclaim: confirmed.
- No spatial validation overclaim: confirmed.
- No deconvolution overclaim: confirmed.
- No unnecessary framework language: greatly reduced.
- Key numbers preserved: confirmed.
- Negative results preserved: confirmed.
- Manuscript compiles: confirmed after LaTeX build.
