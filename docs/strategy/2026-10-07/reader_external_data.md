# Strategy reader: External data landscape (OSDR, GEO, human, atlases)

_Read-only research agent output from workflow wf_1207e53e-984 (2026-10-06), preserved verbatim as structured notes. Input to the next-directions synthesis; not a decision document. Claims cite files/URLs the agent read; verify before acting._

## Summary

I checked every OSDR study from OSD-1 to OSD-1150 through the metadata API (627 are public), plus the biodata API, GEO E-utilities, PRIDE, BioStudies, JAXA ibSLS and the NASA S3 mirror. Only three public bulk rodent kidney transcriptome sets exist beyond the six the project already uses. (1) OSD-457 / GSE152382 (JAXA MHU-3): male C57BL/6J, 31 d flight, dissected 2 days after return (R+2), 3 vs 3 per genotype, processed counts available. (2) OSD-913 / GSE294046 (RRRM-1/RR-8): miRNA only, 57 BALB/c females. It has ISS-terminal (22-24 d) and live-return (39.5 d + 2-4 d) arms with two ground-control types, so its flight-length vs recovery aliasing runs opposite to RRRM-2's. (3) OSD-973: gamma-ray microarray only. No public dataset has a same-strain terminal arm matched to the LAR's 32-day flight. The cheapest partial fix is a cross-mission grid: MHU-3 31 d + R+2 against RRRM-2 LAR 32 d + R+24, anchored by terminal RR-10 (28-29 d) and RR-1 (37 d). The only spaceflight single-cell kidney data are RRRM-1 10x, live-return only, 4 vs 4 (GSE295428). Its kidney tissue/age/treatment labels were officially reset to "unknown" in April 2026, and the matching OSDR studies are staged but not public. One new ground analog: GSE235042 (hindlimb unloading kidney, 3 per group, reports glomerular widening). No public kidney data yet from Bion-M2, MHU-1/2/4/5/8, RR-5/6/9/12/18/19, or the GCR datasets in Siew 2024 (OSD-706..712 are not public). No human dataset measures podocyte markers. Good independent podocyte references exist beyond TMS: GSE146912, GSE111107, GSE160048 and GSE39441 (sorted podocytes). Note: repo HEAD is c730ce1, not a550d39.

## Items

### OSD-457 / GSE152382 (JAXA MHU-3) kidney bulk RNA-seq

- **Status:** unexplored; public; processed counts available; NOT on disk
- **Key facts:** URL https://osdr.nasa.gov/bio/repo/data/studies/OSD-457 (resolves). Kidney is 12 of 192 samples: GSM4614065-76 = MHU3_GC02/06/07/10/11/12 and FL01/02/05/08/10/11. Groups: WT flight 3, WT ground 3, Nrf2-KO flight 3, KO ground 3. Male C57BL/6J, 8 wk at launch. SpX-14, launched 2018-04-02, returned 2018-05-05: 31 d on the ISS, live return, kidney dissected at R+2 after behaviour tests. Library: poly(A), 76 bp PE (ibSLS). GeneLab RSEM counts cover all 192 samples in one matrix: GLDS-457_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv, 30.5 MB, header checked. Runsheet copy in scratchpad osd457/runsheet.csv. Paper: Suzuki et al., Kidney Int 2022, 101:92-105. Repo: config/contrast_vector_framework.yaml lines 115-116 (extension only); docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md:224; MHU-3 PODXL direction from Siew Supp Data 3 in docs/PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md:45.
- **Open questions:** Exact hours from splashdown to dissection (R+2 per ibSLS). Whether GeneLab processing is rRNA-removed for kidney. Whether to pool WT+KO, given Nrf2 effects on the kidney.
- **Relevance to next steps:** The cheapest dataset aimed at the G2 aliasing. Its flight length matches RRRM-2 LAR's flight (31 vs 32 d) but it was sampled at R+2 instead of R+24. If podocyte/structural scores are high at R+2 and near zero in LAR, that favours recovery over shorter exposure (but different strain, sex and lab). Preregister the direction first; with 3 vs 3 (WT) it is descriptive only. The KO arm adds a genotype interaction.

### OSD-913 / GSE294046 (RRRM-1 = RR-8) kidney small-RNA/miRNA

- **Status:** unexplored; public; raw FASTQ on OSDR, processed only on GEO; miRNA only
- **Key facts:** URL https://osdr.nasa.gov/bio/repo/data/studies/OSD-913 (resolves). 57 kidney samples, BALB/cAnNTac female, 3 and 8 months. Design: Flight / Habitat GC / Vivarium GC x ISS-terminal (22-24 d; carcasses frozen on the ISS, dissected after thaw) vs LAR (39.5 d on ISS, returned 2019-01-14, collected 2-4 d after splashdown). n = 4-5 per cell, from biodata API factor values. Problems from the OSDR protocol: food-bar mould, terminal arm cut from a planned 8 weeks to about 3, LAR held 7 d in double-density return housing. OSDR has only FASTQ and MultiQC. GEO has GSE294046_miRNA_complete_quantification_raw.tsv.gz (all 686 samples, URL resolves). Paper: Grandke et al., Nat Commun 2026, doi 10.1038/s41467-026-68737-1. Repo: listed as an extension in config/contrast_vector_framework.yaml:116 only.
- **Open questions:** Whether RRRM-1 kidney total RNA or bulk mRNA exists anywhere (none found on OSDR or GEO). Whether a podocyte miRNA signature is defensible.
- **Relevance to next steps:** The only second-mission kidney dataset with a terminal vs live-return contrast. Its aliasing runs opposite to RRRM-2 (live-return mice flew longer and recovered briefly). It cannot score the 157-gene mRNA set. Use it only for a preregistered miRNA-level hypothesis, e.g. podocyte-enriched miRNAs, flight-vs-control within each arm. The terminal-vs-live-return interaction is confounded with dissection method (thawed carcass vs fresh).

### GSE295428 RRRM-1 kidney scRNA-seq (10x)

- **Status:** public on GEO but compromised: kidney labels set to 'unknown' (Apr 2026); the repo already rejects it
- **Key facts:** Kidney samples GSM8947832-39: 'Kidney, HC' mice 1-4 and 'Kidney, LAR' mice 1-4 (2 x 3 mo and 2 x 8 mo per group). BALB/cAnNTac female, 39.5 d flight, dissected about 2-3 d after splashdown. Per-sample barcodes/features/mtx on GEO (mtx URL returns 200). Barcodes per sample: 3290, 3931, 3917, 4016, 3156, 1350, 5689, 3139 (about 28.5k total). Possible batch confound: LAR mice 1-2 were in pool 4, all others in pool 1. The 2026-04-22 GEO update reset tissue/age/treatment to 'unknown' for six tissues: aorta, kidney, mammary, skin, spleen, thymus. On S3 the aorta (OSD-917), mammary (OSD-925) and skin (OSD-934) sets are staged but return no public study through the API. No kidney OSD was found among OSD-600 to OSD-1150. The repo already rejects these labels: docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md:294-297.
- **Open questions:** When corrected labels or a kidney OSD will be released. Whether kidney label swaps can be detected from marker genes, which would itself be a QC exercise.
- **Relevance to next steps:** Do not use it to separate 'more podocytes' from 'more expression per podocyte' until the labels are corrected. Watch for a corrected GEO/OSDR release. Podocyte yield from whole-kidney 10x is probably very low.

### OSD-973 (GLDS-793) kidney microarray, low-dose-rate gamma radiation

- **Status:** public; marginal relevance
- **Key facts:** Male C57BL/6J, 137Cs gamma for 485 d at 0.032, 0.65 or 13 uGy/min (21, 420 or 8000 mGy total). n = 3 per group, 12 in total, Illumina Sentrix arrays. Euthanised right after irradiation. Paper: doi 10.1269/jrr.09011.
- **Open questions:** Array gene coverage of the 157-gene set.
- **Relevance to next steps:** Low-LET chronic radiation only. Could serve as a negative or analog context for whether podocyte/structural scores move with dose. Not a GCR or spaceflight test. Low priority.

### GSE235042 hindlimb-unloading kidney RNA-seq (ER stress / 4-PBA)

- **Status:** new (public 2026-04-15); unexplored
- **Key facts:** 9 samples: ground control 3, HU 3, HU + 4-PBA 3. C57BL/6J, 3 wk of HU. Sex not stated in GEO. Poly(A) RNA-seq, HISAT2, processed FPKM only (GSE235042_Processed_data.txt.gz resolves); raw reads in PRJNA984237. Paper PMID 41946794 reports HU-induced glomerular widening and loss of Bowman's space, partly reversed by 4-PBA.
- **Open questions:** Sex and age of the mice. Whether the glomerular morphometry data are deposited.
- **Relevance to next steps:** The only public ground-analog kidney transcriptome, and its histology phenotype is glomerular. It can test whether unloading and fluid shift alone move the podocyte/structural scores, without radiation, launch or landing. With 3 per group and FPKM it is descriptive at best. Requantify from SRA for consistency.

### Human kidney proximal-tubule microphysiological systems in microgravity (OSD-516 / GSE268236; GSE302560)

- **Status:** public; in vitro; low relevance
- **Key facts:** OSD-516 / GSE268236: ISS (SpaceX-17) PT-MPS with serum and vitamin D, 73 samples, counts.txt.gz. GSE302560 (public 2026-06-23): calcium-oxalate crystals ± potassium citrate in a PT-MPS in microgravity, 144 samples, counts.txt.gz (resolves). Primary human proximal tubule cells only.
- **Relevance to next steps:** Not informative for podocytes or structural drift. Possible context for any kidney-stone/transport discussion.

### Missions/analogs with NO public kidney omics (verified)

- **Status:** blocked / not available
- **Key facts:** OSDR mission searches: RR-5 (skin, quad, oral), RR-6 (skin, thymus, liver, spleen, colon, lung, feces), RR-9 (liver, retina, gastrocnemius, spleen, thymus, brain), RR-12 (heart, OSD-599), RR-18 (brain regions), RR-19 (plasma OSD-342). RR-17 = RRRM-2. RR-20/21/22/25 return 0 hits. In ibSLS, kidney RNA-seq appears only for MHU-3. MHU-1/2/4/5 have none; MHU-8 is OSD-758/759 (eye). JAXA biospecimen lists contain no kidney. Bion-M1: OSD-209/PXD005102 liver proteome, GSE80223/94381 muscle, no kidney, although recovery groups existed. Bion-M2: launched 2025-08-20, landed 2025-09-19 (30 d, 75 male mice, 65 returned alive), no public data found. Siew 2024 cites OSD-706..712, but these have no public study and no S3 keys. NSRL/GCR kidney: no GEO or OSDR deposit found. A 12-day shuttle kidney Affymetrix study (IJMS 2018, PMC6321533) has no accession. NASA LSDA biospecimens (e.g. RR-1 kidney) are available only by request.
- **Open questions:** Whether NASA BSP holds RRRM-1 or RRRM-2 kidney blocks suitable for WT1/nephrin staining (the direct test listed in the decision report, lines 284-297).
- **Relevance to next steps:** No hidden public terminal arm matched to a 32-day flight. The alternatives are cross-mission modelling or a tissue request through NASA LSDA/ALSDA, or watching for Bion-M2 and Siew's GCR deposits.

### Cross-mission duration x recovery grid from available kidney mRNA data

- **Status:** proposal built from verified metadata
- **Key facts:** Terminal arms: RR-7 (OSD-253) 25 d and 75 d; RR-10 (OSD-462) 28-29 d (B6129 female); RR-1 (OSD-102) 37 d (C57BL/6J female); RR-3 (OSD-163) 39-42 d (BALB/c female); RRRM-2 ISS-T (OSD-771) 55-58 d. Live-return arms: OSD-513 38 d + about 1 d (male C57BL/6J); OSD-457 31 d + 2 d (male C57BL/6J); RRRM-2 LAR 32 d + 24 d. RRRM-1 (miRNA or scRNA only) adds 22-24 d terminal and 39.5 d + 2-4 d live-return arms.
- **Open questions:** Whether a duration meta-regression across 5-6 missions is identifiable given strain/sex heterogeneity.
- **Relevance to next steps:** Lets you place RRRM-2 LAR against a duration curve: the predicted terminal effect at about 32 d from RR-10 and RR-1, the short-recovery effect at about 31 d from MHU-3, and the observed LAR effect. Every comparison crosses strain, sex and lab, so any result is hypothesis-generating. Preregister it like G2.

### Human renal-context data (OSD-656, OSD-575, OSD-571, OSD-530, OSD-483, OSD-353/PXD013305, PXD069732, Twins Table S8)

- **Status:** public or controlled; none measures podocyte markers
- **Key facts:** OSD-656 (on disk): Inspiration4 urine NULISAseq, 203 proteins, 22 samples, 4 crew (L-92/-44/-3, R+1/+45/+82). OSD-575: Inspiration4 serum metabolic panel (creatinine, BUN) plus cytokines, 28 samples. OSD-571: Inspiration4 plasma MS proteomics, metabolomics, cfDNA and cfRNA, 117 samples. OSD-530: JAXA CFE plasma cfRNA, 6 astronauts on the ISS for more than 120 d, 64 samples. OSD-483: plasma sEV small RNA, 14 STS astronauts, 42 samples. OSD-353 / PXD013305: plasma proteome from 21-d head-down bed rest and dry immersion, 46 samples. PXD069732 (2026): serum proteome of 8 astronauts on 180-d missions, including in flight. OSD-903/943: blood RNA and methylation from 2 astronauts. NASA Twins Table S8 urine chemistry is configured at config/human_concordance_prereg.yaml:23 (data/external/human_spaceflight/aau8650_table_s8.xlsx) but is NOT on disk. Polaris Dawn, Axiom and Fram2 have no OSDR deposits (search returned 0 relevant hits); the Polaris 'Stone Risk' urine-calcium data are not public. NASA kidney-stone and biochemical-profile urine data are controlled (LSDA request).
- **Open questions:** Whether SOMA or Polaris Dawn will release urine proteomics or ACR.
- **Relevance to next steps:** Context only; this cannot test the podocyte claim (no albumin/ACR, nephrin, podocin or WT1). The cfRNA sets (OSD-530, OSD-571) would allow a long-shot tissue-of-origin check for glomerular transcripts. Restore Table S8 if the human-concordance stage is kept.

### Independent mouse kidney/glomerular references for item 9 (second atlas)

- **Status:** public; TMS already on disk; podocyte-rich glomerular sets unexplored
- **Key facts:** None of these are among the 8 MKA source studies (docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md:32-33). TMS kidney on disk under data/external/single_cell_atlases/tms_kidney_cellxgene/ (486 podocytes via free_annotation). GSE146912 (Chung 2020): glomerular 10x, 17 samples, all male, 3 C57BL/6J 12-wk controls plus doxorubicin, nephritis and BTBR ob/ob models; processed_data and raw_counts tars. GSE111107 (Karaiskos 2018): glomerular Drop-seq, CD1 8-wk males, 4 replicates plus 2 bulk glomeruli. GSE160048 (He 2021): Smart-seq2 glomerular single cells with mouse 'unbiased', CD45-CD31- and Pdgfrb+ RPKM files, 3222 cells. GSE39441 (Boerries 2013): FACS-sorted podocytes vs non-podocytes, 5 vs 5 pools, Affymetrix exon arrays; a cell-sorted bulk reference. GSE129798 (Ransick 2019): 12 samples, 2 male and 2 female C57BL/6J, cortex/outer/inner medulla. GSE141115 (Denisenko 2020): 54 sc/sn/bulk. GSE240375: young vs aged glomeruli, bulk plus scRNA. GSE64959: GUDMAP adult kidney components. GSE109774/GSE132042 Tabula Muris. GSE108097 Mouse Cell Atlas. GSE190887 (Li 2022 sci-RNA-seq, IRI/UUO, more than 300k nuclei). Pippin et al. 2026 podocyte snRNA atlas (bioRxiv 10.64898/2026.02.26.708349; accession not verified). All GEO suppl directories listed return HTTP 200.
- **Open questions:** Podocyte counts per control sample in GSE146912 and GSE111107. Whether glomerulus-only references break the tier-building rule (≥4-fold above every modelled compartment).
- **Relevance to next steps:** TMS may fail the evaluability gate (few podocytes). GSE146912 controls or GSE111107 give far more podocytes and can define podocyte-high-specificity sets against other glomerular cells. Their glomerulus-only scope means the max-other denominator lacks tubules, so pair them with a whole-kidney atlas (Ransick or Denisenko). GSE39441 gives a cleaner check of podocyte specificity. GSE240375 bears on the age axis in OSD-771.

### Kidney spaceflight proteomics/phosphoproteomics

- **Status:** OSD-462 on disk; OSD-102 and OSD-163 processed proteomics available but not downloaded
- **Key facts:** OSD-102 (RR-1): kidney TMT proteomics, 12 samples, GLDS-102_proteomics_TMT_KDNL.processed.tar.gz, 3.04 GB (URL resolves). OSD-163 (RR-3): 19 samples, GLDS-163_proteomics_RR3_KDN_processed.zip, 20.9 MB (resolves). OSD-462 (RR-10): proteomics 30 and phosphoproteomics 30, workbooks under data/external/osdr/OSD-462/. Searches of PRIDE (spaceflight, microgravity, cosmonaut, astronaut, hindlimb, Bion) found no other rodent kidney proteomics. The BNL-3 GCR kidney proteome exists only as Siew 2024 Supplementary Data 3 consensus scores. PXD001729 (mpkDCT) is referenced in docs, but data/external/phosphoproteomics/ is absent.
- **Open questions:** Coverage of podocyte proteins in the 2014-2016 Orbitrap runs.
- **Relevance to next steps:** OSD-102 and OSD-163 can check podocyte/structural protein abundance (NPHS1, NPHS2, PODXL, SYNPO, NID1). This is a second assay layer on the same animals, not independent replication, and p-values must not be combined (same caveat as the OSD-462 rule in the decision report).

## Risks and constraints

- No public dataset gives a same-strain, same-lab terminal arm matched to RRRM-2 LAR's 32-day flight. Any recovery-vs-duration test from external data crosses strain (C57BL/6J, BALB/c, B6129), sex (OSD-457 and OSD-513 are male) and lab, so it can only be descriptive.
- OSD-457 kidney is 3 vs 3 per genotype and dissected after behaviour tests at R+2. OSD-913 is miRNA only. GSE235042 is 3 per group with FPKM only. All are severely underpowered.
- RRRM-1 single-cell kidney labels were officially reset to 'unknown' (GSE295428 update, 2026-04-22). Kidney live-return mice 1-2 sit in a different pool from the other kidney samples. The OSDR counterparts (OSD-917, 925, 934 seen) are staged on S3 but not public.
- RRRM-1 problems: food-bar mould, terminal arm cut to about 3 weeks, 7 days in double-density return housing, and terminal tissue dissected from thawed carcasses versus fresh live-return tissue. The terminal-vs-live-return contrast is therefore confounded with dissection method.
- Assay layers from the same animals (OSD-102, 163 and 462 proteomics) are not independent replication. Do not combine their p-values.
- No human spaceflight dataset measures podocyte or barrier markers (albumin/ACR, nephrin, podocin, WT1). The Twins Table S8 file the config expects is missing from disk.
- Bion-M2 (2025) and Siew 2024's GCR/NSRL kidney data (OSD-706..712) are not public. Any test that depends on them is blocked.
- Glomerulus-only atlases (GSE146912, GSE111107, GSE160048) lack tubular compartments, so they change the max-other denominator used to build tiers. Pair them with a whole-kidney reference.
- Network: www.ncbi.nlm.nih.gov, biorxiv.org and virtualcellmodels via WebFetch were blocked, so GEO, bioRxiv and CZI details came from curl and E-utilities.
- Repo HEAD is c730ce1, not a550d39 as the task brief states. data/external holds only osdr/, human_spaceflight/OSD-656 and single_cell_atlases/ (MKA, TMS). The GSE228367, GSE150338, GSE269622/9 and PXD001729 paths named in configs and docs are not on disk.

## Sources read

- /home/user/RRRM2_Kidney_Transcriptome/config/contrast_vector_framework.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/clinical_renal_axes_cross_mission.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/clinical_axes_recovery_persistence.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/human_concordance_prereg.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/human_urine_marker_panel.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/dct_subtype_reference_freeze_v1.yaml
- /home/user/RRRM2_Kidney_Transcriptome/docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/v11_execution_research_plan.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CODE_COMMENTARY.md
- https://osdr.nasa.gov/osdr/data/search?term=kidney&type=cgene
- https://osdr.nasa.gov/osdr/data/osd/meta/<1..1150> (scanned; saved under scratchpad ext/meta/)
- https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?study.characteristics.material%20type=/[Kk]idney/
- https://visualization.osdr.nasa.gov/biodata/api/v2/query/metadata/?id.accession=OSD-913
- https://osdr.nasa.gov/osdr/data/osd/files/913
- https://osdr.nasa.gov/osdr/data/osd/files/457
- https://nasa-osdr.s3.amazonaws.com/?list-type=2&prefix=OSD-917/ (also OSD-925, OSD-934 ISA zips)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE295428 (SOFT via curl)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE294046
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE152382
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235042
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302560
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268236
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146912
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111107
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160048
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129798
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141115
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE39441
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240375
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324989
- https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi (GEO kidney/spaceflight/HU/GCR/glomerular queries)
- https://ibsls.megabank.tohoku.ac.jp/overview
- https://api.biorxiv.org/details/biorxiv/10.1101/2025.09.08.675003
- https://api.biorxiv.org/details/biorxiv/10.64898/2026.02.26.708349
- https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=spaceflight (and microgravity, cosmonaut, astronaut, hindlimb, Bion)
- https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD069732
- https://www.ebi.ac.uk/biostudies/api/v1/search?query=kidney%20spaceflight
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11167060/ (Siew 2024 data availability via search)
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6321533/fullTextXML
- https://virtualcellmodels.cziscience.com/dataset/mouse-kidney (via curl)
- https://en.wikipedia.org/wiki/Bion-M_No.2
- https://www.nature.com/articles/s41526-025-00465-0
- https://lsda.jsc.nasa.gov/cf/scripts/biospecimens/bio_search_start_adv.cfm (via search result)
