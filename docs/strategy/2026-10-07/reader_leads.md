# Strategy reader: Leads in the newest results (G1/G2 run outputs)

_Read-only research agent output from workflow wf_1207e53e-984 (2026-10-06), preserved verbatim as structured notes. Input to the next-directions synthesis; not a decision document. Claims cite files/URLs the agent read; verify before acting._

## Summary

Read-only review of the newest results (origin/main a550d39, which matches local c730ce1), with my own scratchpad recomputations. The recomputed mission effects match the official tables to 1e-15. Five leads stand out.

(1) The podocyte high-specificity (HS) shift is the most robust signal. It stays rank 1 of 49 after correcting for gene GC content, length and expression (g 0.65–0.71, lower CI above 0). It also stays rank 1 after adjusting for tissue-sampling marker panels (cortex, medulla, blood, hilar smooth muscle): 0.198→0.165 score units, CI 0.009–0.322. The structural, endothelial and mesenchymal shifts mostly disappear under the same adjustment.

(2) OSD-457 (JAXA MHU-3) is a live-return kidney RNA-seq cohort on OSDR that the project has not scored. It shows podocyte HS g 1.16 (wild-type only 1.36), rank 4 of 49, and its 49-set pattern correlates r=0.61 with the pooled terminal pattern. With OSD-513 it forms a descriptive time-since-return series: 0 d about 0.7–1.2, 1–2 d 1.2–1.6, 24 d (LAR) 0.14.

(3) The ISS-T vs LAR anti-correlation (r=-0.65) is not a rebound. Contrasts that use no flight animals (basal minus ground control) show r=-0.65 too, while flight minus vivarium, which uses no ground-control group, shows r=+0.07. The ERCC spike-in mix is confounded with condition within each arm and reversed between arms.

(4) The ISS-T flight group carries on-orbit handling and sampling signatures: immediate-early genes g 2.14, medulla 2.71, hemoglobin 1.73. These explain the exploratory interaction hits (endothelial, mesenchymal, TAL, fibrosis).

(5) Gene GC-content bias is strong in the LAR flight group (rho -0.47). Correcting for it roughly halves the G2 podocyte interaction.

The main threat is control choice: podocyte flight-vs-vivarium meta g is 0.28 (CI -0.45 to 1.01). All of this is post hoc.

## Items

### L1. Podocyte HS shift survives GC/length/expression correction and tissue-sampling covariates (strongest; new robustness, post hoc)

- **Status:** alive - strongest lead; direction robust; still not separable from structural in score units (K0 stands)
- **Key facts:** Recomputation reproduces official per-mission effects exactly (max abs diff 6.7e-16 vs compartment_context_mission_effects.tsv). Per-sample GC/length/expression residualization of gene z-scores, a CQN-like correction (scratchpad cqn_like.py → adj_effects_*.tsv), keeps the 5-mission meta at rank 1 in every mode: GC 0.648 (0.009, 1.287); length 0.675 (0.034, 1.316); GC+length 0.662 (0.024, 1.299); GC+length+expression 0.705 (0.065, 1.344). I² stays 0. Sampling-panel covariates are cortex_PT (Slc34a1/Lrp2/Slc5a2/Slc22a6/Kap/Slc22a8), medulla, blood (Hb, Alas2) and vsmc_hilar. Adjusting with pooled within-cell slopes: podocyte score units 0.198→0.165 (0.009, 0.322); g scale 0.602 (−0.038, 1.243), still rank 1 of 49 (t 2.61 vs 2.96). Structural + sampling gives 0.127 (0.005, 0.249). Structural alone, pooled within-cell slope, gives 0.112, reproducing K0. Under the same sampling adjustment the structural scaffold falls from 0.549 to 0.187 g, endothelial and mesenchymal fall to about 0, and principal_cell broad falls from 0.53 to 0.07. Within-cell animal correlations (n=141): podocyte HS vs structural 0.58; vs vsmc 0.35; vs glomerular endothelium 0.31; vs cortex_PT 0.26; vs blood −0.02; vs medulla −0.08. The structural drift is not circular: its 351 genes with no compartment enrichment still shift (mean gene g 0.093 vs podocyte HS 0.219). The structural set does contain 26 podocyte-, 34 endothelial- and 26 mesenchymal-enriched genes.
- **Open questions:** Covariate panels were chosen after seeing the data, and the covariates may themselves respond to flight. Per-mission (K0-style) covariate slopes with n=12–20 overfit and attenuate to 0.079. The medulla panel contains TAL and principal-cell markers, so those adjustments are circular; the podocyte adjustment is not.
- **Relevance to next steps:** Strongest candidate for a preregistered G3 'technical and sampling robustness' item: CQN/EDASeq GC-length normalization plus prespecified sampling-marker covariates. If locked and passed, the claim becomes 'the podocyte-leaning tail survives technical and sampling correction, while the structural, endothelial and stromal drift mostly does not'. That sharpens K0 without overturning it.

### L2. OSD-457 (JAXA MHU-3) bulk kidney RNA-seq exists on OSDR, has not been scored, and replicates the podocyte shift descriptively

- **Status:** unexplored → quick win; hypothesis-generating (n=3+3 per genotype)
- **Key facts:** Source: OSDR OSD-457 runsheet (GLDS-457_rna_seq_bulkRNASeq_v2_runsheet.csv). It holds 12 kidney samples: wild-type 3 flight / 3 ground control and Nrf2-KO 3/3. Design: male C57BL/6J background, 12 weeks at launch, 31 d on ISS, live return; mice went through behavioural tests, then isoflurane and exsanguination (protocol text from https://osdr.nasa.gov/osdr/data/osd/meta/457; Suzuki et al. 2020, doi:10.1038/s42003-020-01227-2). My standalone scorer (median-of-ratios log2 plus the same z/CPM-eligibility rules) was validated on OSD-513: r=0.995 against the official 49-set g, max |Δ| 0.36, podocyte HS 1.46 vs 1.56. MHU-3 genotype-blocked: podocyte HS g 1.16 (wild-type only 1.36), rank 4 of 49; podocyte all_enriched 1.14; structural 0.73 (wild-type 0.61); endothelial 1.12; PT −0.75; blood −2.17; immediate-early genes 1.50. Its 49-set pattern correlates r=0.61 with the pooled terminal estimates, 0.62 with OSD-513 and 0.55 with OSD-771 ISS-T. The project's audit states 'No sixth mission: none exists' (docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md §12). The decision report mentions OSD-457 only as lacking urine and histology endpoints (docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md L224).
- **Open questions:** VST files from the GeneLab pipeline may be absent, so it must be rerun through the project loaders. Half the animals are Nrf2-KO. n is tiny. The exact time from return to euthanasia needs verifying. It may be a live-return moderator rather than a terminal replicate.
- **Relevance to next steps:** Add OSD-457 as a frozen descriptive moderator cohort, like OSD-513, through the project's loaders (about 1 day of work). Two live-return cohorts showing the shift 1–2 days after landing materially strengthen the recurrence argument and feed the recovery question.

### L3. Time-since-return series: podocyte shift present 1–2 d after landing, gone 24 d after (descriptive recovery lead)

- **Status:** hypothesis-generating; cross-study, confounded
- **Key facts:** Podocyte HS g by days since return: 0 d terminal (pooled 0.69; OSD-771 ISS-T 1.17); about 1 d OSD-513 1.56 (euthanized 14-Jan-2021, splashdown 13-Jan); about 1–2 d MHU-3 1.16; 24 d OSD-771 LAR 0.14. Flight durations of the live-return cohorts (31 d MHU-3, 37–38 d OSD-513) are close to LAR's 32 d. That weakens the doc's alias 'LAR's shorter flight alone could explain a smaller shift'. Counter-evidence: OSD-253 C57BL/6J day 25 g 0.04. The doc's interaction is −0.281 (−0.622, 0.061) and the Fieller ratio 0.10 (−0.66, 1.08) (persistence_contrasts.tsv).
- **Open questions:** Strain differs (B6J vs B6NTac), sex differs (male vs female), and preservation and study differ. Group-level batch effects (L6) could dominate any single cohort.
- **Relevance to next steps:** Use it to motivate, not claim, a recovery interpretation. A cross-study 'days since return' moderator figure is legitimate as descriptive context. The definitive design still needs a terminal arm matched to the 32 d flight, which the doc's §6 already notes.

### L4. The ISS-T vs LAR anti-correlation (r=−0.65) is a ground-control-group artifact, not 'rebound'

- **Status:** strong negative finding: retire the rebound interpretation
- **Key facts:** Recomputed r=−0.648 (Spearman −0.66) over 49 sets (persistence_pattern_concordance.tsv). It is broad, not driven by a few compartments: compartment-level r=−0.68 (11 compartments); leave-one-compartment-out −0.54 to −0.72; dropping endothelial and mesenchymal gives −0.50. The decisive test (scratchpad mixtest.py) uses all four OSD-771 groups per arm (data/processed/vst_normalized/GLDS-674…VST): BSL−GC r=−0.65 and VIV−GC r=−0.48, both without flight animals, while FLT−VIV r=+0.07 and FLT−BSL −0.23. The two ground-control groups therefore deviate in opposite directions along a shared axis. The ERCC spike-in mix is perfectly confounded with condition within arm and reversed between arms (data/raw/metadata/a_OSD-771_…Illumina.txt): ISS-T FLT and VIV are Mix 1, GC and BSL Mix 2; LAR FLT and VIV are Mix 2, GC and BSL Mix 1. Each group was also euthanized on its own day: LAR FLT 19-Sep, GC 18-Sep, VIV 20-Sep. Raw count: 35 of 49 sets have LAR g < ISS-T g.
- **Open questions:** The control groups could also differ biologically, so a processing-group (batch) explanation versus a real ground-control difference cannot be fully separated.
- **Relevance to next steps:** Do not pursue 'rebound'. Report that the G2 interaction is aliased not only with recovery, duration and euthanasia route but also with group-level batch (θ_L−θ_T absorbs any GC-group deviation). This argues for using BSL and VIV as within-arm calibration controls in any persistence claim.

### L5. OSD-771 ISS-T flight group carries on-orbit handling/sampling signatures that explain the exploratory interaction 'hits'

- **Status:** moderate-strong technical explanation
- **Key facts:** Flight g in ISS-T for marker panels: immediate-early genes (Fos/Jun/Egr1/Atf3…) 2.14; medulla (Aqp2/Slc14a2/Umod/Slc12a1/Ptger3/Aqp1) 2.71; mesangial (Itga8/Pdgfrb/Des) 2.10; blood 1.73; vsmc/hilar 1.34; hypoxia 0.86; mitochondrial mRNA −1.12. In LAR all are near 0 (IEG 0.03, medulla −0.41, blood 0.24). Fibrosis-axis collagens in ISS-T: Col1a1 1.87, Col1a2 1.90, Col3a1 2.30 (OSD-462 is opposite, about −1.0). Protocol, from the OSD-771 i_Investigation: ISS-T kidneys were dissected fresh on orbit by crew into RNAlater for 48 h at 4°C, then MELFI, with ketamine/xylazine; LAR tissues were snap-frozen in LN2 after gas anesthesia. Each kidney's RNA comes from a middle coronal slab containing cortex, pyramid, calyx, arteries and veins, so sampling composition varies. TAL HS ISS-T g 1.43 matches medulla over-representation. Endothelial and mesenchymal I² is 72–76%, driven by OSD-462 (−1.0 to −1.4) vs OSD-771 (+1.4 to +1.7). Consistent with Grey60's 'OSD-771 ISS-T vascular/stromal/ECM' result (docs/owner_decisions.md).
- **Open questions:** Whether OSD-513's equally strong mesenchymal (2.3) and blood (1.46) shift reflects post-landing physiology or dissection.
- **Relevance to next steps:** Label the exploratory endothelial, mesenchymal, TAL and fibrosis interaction passes as handling or sampling-confounded. Add the marker panels as a standard QC table in any future cohort.

### L6. Control-choice and group-level variance: the podocyte effect is weaker against vivarium controls (main threat)

- **Status:** moderate; most serious remaining threat to 'spaceflight' attribution
- **Key facts:** Podocyte HS: FLT−GC meta across 6 cohorts 0.81 (0.22, 1.39). FLT−VIV, C57BL/6J only, k=4: 0.28 (−0.45, 1.01), with OSD-253 0.39, OSD-462 0.47, OSD-513 0.46, OSD-771 ISS-T −0.17. FLT−BSL, k=3: 1.96 / −0.14 / 0.47. Structural FLT−VIV is 0.55, stronger than podocyte, so the podocyte > structural ordering depends on control choice. Vivarium mice score higher than habitat ground controls on podocyte in 6 of 8 strata (e.g., ISS-T VIV−GC 1.19; OSD-513 GC−VIV −0.90; OSD-253 C3H 75 d −1.28). Across 35 control-vs-control contrasts (non-independent; 24 from OSD-253), podocyte HS RMS g is 0.88 against a mean sampling variance of 0.42, an excess group-level variance of about 0.36. Inflating the meta variance by that amount gives a CI of (−0.31, 1.68) and z-test p 0.055. In OSD-253 the original GC differs strongly from the rerun white-light GC for B6 (GC−GCRw 1.19 / 1.86), so FLT−GCRw is 0.77 / 3.00. Source: scratchpad nullcontrasts.py → null_contrasts.tsv.
- **Open questions:** Vivarium controls differ in housing hardware, so FLT−VIV mixes flight with hardware effects. The excess-variance estimate is an upper bound because contrast types are heterogeneous.
- **Relevance to next steps:** Run the already-recommended cross-mission control-choice sensitivity (owner_decisions 2026-07-29) as a preregistered item: FLT vs VIV, FLT vs BSL, plus empirical group-variance calibration from control-control contrasts. Also report an empirically calibrated interval. This is the check most likely to change the headline.

### L7. Gene GC-content bias is cohort- and group-specific; it inflates or deflates set effects and roughly halves the G2 podocyte interaction

- **Status:** moderate technical lead; needs a locked correction
- **Key facts:** Spearman correlation between gene-level flight g and GC%: OSD-102 0.07, OSD-163 0.36, OSD-253 B6 −0.22, C3H 0.56 (R² 0.34 with length and expression), OSD-462 0.27, OSD-771 0.17, OSD-513 0.14, LAR −0.48 (R² 0.25). GC came from Ensembl BioMart (scratchpad dl/genes.tsv). At set level, set median GC predicts per-set g: ISS-T r=0.60, LAR −0.49, OSD-462 −0.42, OSD-513 0.47. Endothelial sets are GC-rich (48–50%) and DCT1 GC-poor (40–43%); podocyte HS (45.3%) sits near the background of 45.9%. The LAR bias is specific to the LAR flight group (FLT−VIV ρ −0.46; VIV−GC −0.003) and driven by a few samples: LAR FLT OLD per-sample GC slope SD 0.67 vs about 0.1 in GC. After GC correction: podocyte LAR g 0.14→0.31 (GC only) / 0.45 (GC+length), ISS-T 1.17→0.81 / 0.71. ISS-T vs LAR r goes from −0.65 to −0.48. Within the podocyte set, gene effect correlates with GC (ρ 0.29, p<0.001) but the 5-mission meta is unchanged by correction.
- **Open questions:** Whether the GC bias reflects library or PCR batch or RNA quality; check per-sample GC slope against RIN, 3'/5' coverage ratio and extraction batch.
- **Relevance to next steps:** Add a GC/length-normalized (CQN/EDASeq) sensitivity to G2. The persistence labels are probably technically fragile, so treat any LAR decay estimate cautiously.

### L8. Glomerular barrier 'opposite direction' and podocyte HS up are one coherent ~6% podocyte-transcript increase; compositional vs per-cell cause unresolved

- **Status:** coherent at transcript level; mechanism hypothesis-generating
- **Key facts:** Barrier axis −0.716 (−1.36, −0.07), I² 0, maxT 0.0018, Holm p_mkh 0.147 (primary_meta_results.tsv). The signs mean barrier transcripts are higher in flight. Nphs2 is higher in all 5 missions (g 0.73, 1.29, 1.30, 0.18, 1.81); Nphs1 and Wt1 in 4 of 5 (descriptive_signed_gene_effects.tsv). The barrier genes are a subset of the HS program; the earlier disjoint-proxy adjustment absorbed them (0.507→0.070). Mean VST log2 flight−GC for podocyte HS: OSD-163/253/771/513 about +0.09 (6–7%), C3H +0.24 (18%), OSD-462 +0.04, OSD-102 0. Siew et al. 2024 (Nat Commun 15:4923, Europe PMC PMC11167060) report higher kidney weight relative to body weight with spaceflight or MG+GCR, and in RR-10 lower DCT density with larger retained tubules. A compositional route (relatively less tubular mass, so a larger glomerular share) is plausible. But glomerular endothelial markers (Ehd3 −0.76 to 1.71) and mesangial markers (Itga8 −2.43 to 2.20) are inconsistent across missions. Podocyte rises where blood and vascular panels fall (OSD-462, MHU-3), which argues against blood or hilar content but leaves 'more podocytes' vs 'more per-cell transcription' open. Siew's vote-count reported PODXL down in RR-7, RR-10 and MHU-3; here Podxl g is up (0.74 B6, 1.83 C3H, 0.10 OSD-462), and the GeneLab DGE for OSD-462 gives Space Flight vs Ground Control log2FC +0.017 (total RNA) and +0.33 (mRNA).
- **Open questions:** Whether glomerular number or volume per mg changes; whether kidney weight change is tubular. A glomerular-endothelium HS set was not evaluable in MKA.
- **Relevance to next steps:** Frame as 'bulk podocyte-associated transcript abundance up about 6%; direction opposite to injury'. The decisive next evidence is morphometry: WT1+ nuclei per glomerulus, glomerular density and area on existing RR-10 or RRRM-2 sections (the planned Sarder pipeline).

### L9. OSD-253 exposure pattern across all 49 sets; C3H/HeJ strongly replicates the podocyte shift

- **Status:** descriptive; strain-dependent dose pattern
- **Key facts:** B6 day 75 > day 25 in 36 of 49 sets. Mean |g| rises from 0.32 to 0.82, and r(day 25, day 75) across sets is 0.24. Day 75 B6: podocyte tiers 1.73–2.33; mesenchymal 1.16–1.64; endothelial 0.57–1.20; structural 1.59; immune about 1.0; DCT1 about −0.65; DCT2/CNT −0.26 to −0.45. Day 25 B6 is dominated by TAL down (−1.17 to −1.46). C3H/HeJ, not in the primary pool: podocyte HS 1.43 at day 25 and 2.19 at day 75; podocyte all_enriched 1.37 / 2.25; structural 0.97 / 1.58; mesenchymal broad 0.34 / 2.07. So C3H shows the podocyte shift already at 25 d. Day-75 B6 GC n=4. Source: scratchpad all_scores.tsv, project loader load_osd253_strain_sensitivity.
- **Open questions:** Day-25 and day-75 carcasses had different on-ISS storage times. The original ground controls were blue-light and the rerun controls differ by batch.
- **Relevance to next steps:** The 'shorter flight explains LAR' alias rests mainly on B6 day 25. C3H and the 31–38 d live-return cohorts weaken it. Report strain alongside duration.

### L10. OSD-513 podocyte-specificity vs terminal non-separability

- **Status:** hypothesis-generating
- **Key facts:** OSD-513 (RR-23): male C57BL/6J, 16–17 wk, 38 d flight, live return, euthanized 1 d after splashdown at TAMU, isoflurane plus thoracotomy, RNAlater 24 h at 4°C. Groups were dissected on separate days (FLT 14-Jan, GC 17-Jan, VIV 20-Jan), so condition is fully confounded with date. All terminal cohorts are female (RR-1/3/7/10 and RRRM-2). OSD-513 flight shows PT panel −1.45, mesangial 2.10 and blood 1.46. The structural score tracks PT within cells (r 0.33), which pulls the structural set down. After sampling adjustment, OSD-513 podocyte is 2.03 and structural 0.10, so the specificity survives. MHU-3 (also male live-return) is less specific (podocyte 1.16, structural 0.73). OSD-513's 49-set pattern correlates r=0.65 with ISS-T and 0.37 with the pooled terminal pattern (p 0.18, persistence_pattern_concordance.tsv).
- **Open questions:** Sex vs post-landing reloading vs dissection date; k=1–2.
- **Relevance to next steps:** Not pool-able. Only useful as a moderator together with MHU-3. Do not interpret as recovery.

### L11. Top contributing genes: biology-consistent and only weakly technical

- **Status:** hypothesis-generating (sub-module test possible)
- **Key facts:** From podocyte_gene_contributions.tsv, the top 15 contributors are consistent across missions: Mgat5, Arhgef16, Aplp1, Dcdc2a, Mapt, Nphs2, Col4a3, St3gal3, E2f1, Optn, Rhpn1, Myo1e, Wt1, Dock5, Rab3b; 106 of 145 genes are positive. Podocyte HS genes are long (median span 81 kb vs 26 kb background; Mgat5 283 kb, St3gal3 203 kb, Myo1e 192 kb). Within the set, though, span (ρ −0.16, p 0.06) and transcript length (ρ −0.06) do not predict effect. GC (ρ 0.29, p<0.001) and expression (ρ 0.21, p 0.014) do. Atlas specificity has ρ +0.14 (p 0.10): more-specific markers do not shift less, which argues against a contaminating non-podocyte source. Leave-top-5 drops the lower CI to −0.044 (podocyte_leave_top_k.tsv). Functional clusters: N-glycan branching and sialylation (Mgat5, St3gal3; Podxl rank 35), i.e. the glycocalyx; Rho/actin (Arhgef16, Rhpn1, Myo1e, Dock5, Nck2, Srgap1, Arhgap28, Mapt); slit diaphragm and GBM (Nphs1, Nphs2, Crb2, Col4a3); cell-cycle brakes (Cdkn1c, E2f1). Ambient RNA matters only for atlas derivation and is partly covered by the 8/8 LOSO rebuilds.
- **Open questions:** Whether a glycocalyx or sialylation sub-module carries more than its share (needs a prespecified sub-module test). OSD-462 proteome and phospho data could check Mgat5/St3gal3/Podxl, but are aliased with reporter tag.
- **Relevance to next steps:** Low priority. Optionally preregister 3–4 functional sub-modules before any further gene-level storytelling.

### L12. Technical QC: RIN imbalance negligible; 3'-bias and mapping track structural/endothelial scores

- **Status:** resolved / minor
- **Key facts:** Within-cell correlations (n=121 with QC): RIN vs podocyte −0.05, vs structural −0.14. A RIN-adjusted LAR podocyte goes 0.14→0.17 and structural −0.42→−0.34. ISS-T and OSD-513 change by 0.04–0.07 or less. This closes the doc §6 open item. The 3'/5' gene-body coverage ratio correlates with structural (−0.36), endothelial (−0.44) and mesenchymal (−0.30) scores, and uniquely-mapped % with structural (0.33) and podocyte (0.24). Imbalances from technical_qc.tsv: OSD-163 mapping rate g −1.28 (flagged), 3'/5' −1.06, depth +1.20; OSD-102 depth +1.06; LAR RIN +0.93 (p 0.045). Also the ERCC mix confound (L4) and the group-specific GC bias (L7).
- **Open questions:** Whether the 3'-bias in OSD-163 flight inflates its structural and endothelial (not podocyte) effects.
- **Relevance to next steps:** Report RIN as non-consequential. Include 3'/5' ratio and the ERCC mix as covariates or flags in a technical sensitivity.

### L13. DCT1 down and principal-cell up: secondary compartment signals

- **Status:** weak; DCT1 hypothesis-generating, principal-cell likely sampling
- **Key facts:** DCT1 tiers pooled −0.25 to −0.28, I² 0, 4 of 5 negative. OSD-462 (RR-10) is +0.32 to +0.56; Siew found lower DCT density there, with larger tubules. Other DCT1 values: OSD-513 −0.73 to −1.34, ISS-T −0.53 to −0.69, B6 day 75 about −0.65; mean log2 −4% to −8% in OSD-253/513/771. Principal_cell broad is 0.53 (I² 0, second compartment), but falls to 0.07 after sampling adjustment, which is circular because the medulla panel includes collecting-duct markers. TAL becomes −0.5 after adjustment, also circular.
- **Open questions:** Whether DCT1 down plus podocyte up is a single compositional axis (tubule loss).
- **Relevance to next steps:** Could be tested jointly with L8 morphometry: DCT and glomerular density on the same sections. Not a standalone paper given the project's DCT history.

### L14. Exploratory G2 49-set hits (endothelial×4, TAL HS, mesenchymal broad, fibrosis axis)

- **Status:** retire as biology; explained by L4/L5/L7
- **Key facts:** From persistence_contrasts.tsv: endothelial broad interaction −0.678, FL p 0.0004, maxT 0.005; TAL HS −1.015, maxT 0.015; mesenchymal broad −0.711, maxT 0.031; fibrosis axis −0.846, FL p 0.003. All are driven by large ISS-T flight elevations (1.4–2.1) that co-occur with the IEG, medulla, blood and hilar signatures, plus opposite ground-control-group deviations and GC bias. Orientation flips make 'positive' TAL broad and immune interactions artefacts of sign convention: in raw values they are ISS-T up, LAR down.
- **Open questions:** None worth pursuing without new cohorts.
- **Relevance to next steps:** State explicitly in the write-up that these are not recovery biology.

## Risks and constraints

- Every new number here is post hoc and unregistered. It was computed in the scratchpad (/tmp/claude-0/-home-user-RRRM2-Kidney-Transcriptome/6b4738ea-7428-58b8-9974-e262e0076857/scratchpad/: score_all.py, gene_effects.py, cqn_like.py, nuisance.py, mixtest.py, nullcontrasts.py, standalone_score.py and the output TSVs). Anything used in a paper must be preregistered and rerun through the rrrm2 CLI.
- Repo untouched: git status is clean. origin/main is a550d39; the working branch claude/sleepy-lamport-9b5zqs is at c730ce1.
- The sampling-marker panels were chosen by hand after seeing results. The medulla and cortex panels contain TAL, principal-cell and PT markers, so adjustments of those compartments are circular. Covariates may themselves be flight-affected (mediator or collider).
- MHU-3 (OSD-457) was scored with a standalone median-of-ratios log2 pipeline, not GeneLab VST. It was validated only on OSD-513 (r 0.995). It has n=3+3 per genotype, half the animals are Nrf2-KO, and the exact post-return timing is unverified.
- The control-vs-control 'null' contrasts are non-independent (24 of 35 from OSD-253, with shared groups) and heterogeneous (housing, light, sequencing batch). The excess group-variance estimate (0.36) is an upper bound.
- In the score-unit K0 metric, every adjustment variant still has the podocyte interval near zero. K0's 'not separable' verdict should not be upgraded from these sensitivities.
- In OSD-771, OSD-513 and others, condition is confounded with dissection day and, in OSD-771, with ERCC mix. Group-level batch effects cannot be removed by any within-mission model.
- WebFetch to pmc.ncbi.nlm.nih.gov and sainsburywellcome.org was blocked by egress, so the Siew 2024 text was read through the Europe PMC REST full text.
- Bulk RNA cannot separate cell number from per-cell expression. All mechanistic readings (composition, glycocalyx, recovery) remain hypotheses until histology or morphometry exists.

## Sources read

- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/CLINICAL_RENAL_AXES_DECISION_REPORT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PODOCYTE_STRICT_MATCHING_AUDIT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PODOCYTE_K1_K4_ADVERSARIAL_AUDIT_2026-08-11.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/PROJECT_TECHNICAL_AUDIT_2026-08-25.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/DEFERRED_SECOND_ATLAS_VALIDATION_2026-10-06.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/owner_decisions.md
- /home/user/RRRM2_Kidney_Transcriptome/docs/plan_changelog.md
- /home/user/RRRM2_Kidney_Transcriptome/config/clinical_renal_axes_cross_mission.yaml
- /home/user/RRRM2_Kidney_Transcriptome/config/contrast_vector_framework.yaml
- /home/user/RRRM2_Kidney_Transcriptome/src/clinical_axes/analysis.py
- /home/user/RRRM2_Kidney_Transcriptome/src/clinical_axes/statistics.py
- /home/user/RRRM2_Kidney_Transcriptome/src/clinical_axes/data.py
- /home/user/RRRM2_Kidney_Transcriptome/scripts/clinical_axes/run_compartment_context.py
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_contrasts.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_arm_effects.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_pattern_concordance.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_duration_context.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_adjusted.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_technical_qc.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_cohort_manifest.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_persistence/persistence_terminal_reference.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/podocyte_sensitivities/podocyte_gene_contributions.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/compartment_context_meta_results.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/compartment_context_mission_effects.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/compartment_context_manifest.json
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/podocyte_scaffold_manifest.json
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/primary_meta_results.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/primary_mission_effects.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/descriptive_signed_gene_effects.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/recovery_moderator_effects.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/technical_qc.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/results/clinical-axes/20261006T004825Z_g1g2/sample_manifest.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv
- /home/user/RRRM2_Kidney_Transcriptome/data/raw/metadata/s_OSD-771.txt
- /home/user/RRRM2_Kidney_Transcriptome/data/raw/metadata/a_OSD-771_transcription-profiling_rna-sequencing-(rna-seq)_Illumina.txt
- /home/user/RRRM2_Kidney_Transcriptome/data/external/osdr/OSD-771/metadata/i_Investigation.txt
- /home/user/RRRM2_Kidney_Transcriptome/data/external/osdr/OSD-102/metadata, OSD-163, OSD-253, OSD-462, OSD-513 s_*.txt and i_Investigation.txt
- /home/user/RRRM2_Kidney_Transcriptome/data/external/osdr/OSD-462/GLDS-462_rna_seq_differential_expression_totRNA_GLbulkRNAseq.csv (and _mRNA)
- /home/user/RRRM2_Kidney_Transcriptome/data/external/osdr/*/GLDS-*_runsheet*.csv and VST/RSEM count matrices
- https://osdr.nasa.gov/osdr/data/search?term=kidney&type=cgene
- https://osdr.nasa.gov/osdr/data/osd/files/457
- https://osdr.nasa.gov/osdr/data/osd/meta/457
- https://osdr.nasa.gov/geode-py/ws/studies/OSD-457/download?source=datamanager&file=GLDS-457_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv
- https://www.ensembl.org/biomart/martservice (mouse gene coordinates, GC%, canonical transcript length)
- https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11167060/fullTextXML (Siew et al. 2024, Nat Commun 15:4923)
- https://pmc.ncbi.nlm.nih.gov/articles/PMC11167060/ (search result; fetch blocked)
- https://amaral.northwestern.edu/publications/aging-associated-systemic-length-associated-transcriptome-imbalance/citation (Stoeger et al., length-associated transcriptome imbalance)
- /tmp/claude-0/-home-user-RRRM2-Kidney-Transcriptome/6b4738ea-7428-58b8-9974-e262e0076857/scratchpad/ (all_scores.tsv, my_effects.tsv, gene_effects.tsv, adj_effects_*.tsv, nuisance_scores.tsv, mixtest_effects.tsv, null_contrasts.tsv, sa457.tsv, sa457wt.tsv, sa513.tsv)
