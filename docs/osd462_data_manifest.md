# OSD-462 Local Data Manifest

Local data root: `data/external/osdr/OSD-462/`

OSD-462 is the RR-10 kidney multi-omic study that includes RNA-seq, proteomics, and phosphoproteomics assays. For this project, the useful first-pass files are the processed RNA-seq matrices/tables, the processed proteomics workbook, and the ISA metadata that maps samples and assays. Vendor raw mass-spectrometry archives are not required unless we decide to reprocess proteomics from raw `.raw` files.

## Downloaded And Loaded

| File or folder | Status | Use |
| --- | --- | --- |
| `OSD-462_metadata_OSD-462-ISA.zip` | Downloaded | Original ISA metadata bundle from OSDR. |
| `metadata/` | Extracted | Unpacked ISA investigation, study, RNA-seq, proteomics, and phosphoproteomics assay metadata. |
| `GLDS-462_rna_seq_SampleTable_UPX_GLbulkRNAseq.csv` | Downloaded | UPX RNA-seq sample-condition mapping. |
| `GLDS-462_rna_seq_SampleTable_mRNA_GLbulkRNAseq.csv` | Downloaded | mRNA RNA-seq sample-condition mapping. |
| `GLDS-462_rna_seq_SampleTable_totRNA_GLbulkRNAseq.csv` | Downloaded | total RNA-seq sample-condition mapping. |
| `GLDS-462_rna_seq_VST_Counts_UPX_GLbulkRNAseq.csv` | Downloaded | Variance-stabilized UPX expression matrix for PCA, correlation, enrichment scoring, and network-style analyses. |
| `GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv` | Downloaded | Variance-stabilized mRNA expression matrix; likely the primary transcriptome table for comparison to proteomics. |
| `GLDS-462_rna_seq_VST_Counts_totRNA_GLbulkRNAseq.csv` | Downloaded | Variance-stabilized total RNA expression matrix. |
| `GLDS-462_rna_seq_RSEM_Unnormalized_Counts_UPX_GLbulkRNAseq.csv` | Downloaded | Raw-ish UPX gene count matrix for re-running differential expression if needed. |
| `GLDS-462_rna_seq_RSEM_Unnormalized_Counts_mRNA_GLbulkRNAseq.csv` | Downloaded | Raw-ish mRNA gene count matrix for re-running differential expression if needed. |
| `GLDS-462_rna_seq_RSEM_Unnormalized_Counts_totRNA_GLbulkRNAseq.csv` | Downloaded | Raw-ish total RNA gene count matrix for re-running differential expression if needed. |
| `GLDS-462_rna_seq_differential_expression_UPX_GLbulkRNAseq.csv` | Downloaded | OSDR-provided UPX differential-expression results. |
| `GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv` | Downloaded | OSDR-provided mRNA differential-expression results. |
| `GLDS-462_rna_seq_differential_expression_totRNA_GLbulkRNAseq.csv` | Downloaded | OSDR-provided total-RNA differential-expression results. |
| `GLDS-462_rna_seq_contrasts_UPX_GLbulkRNAseq.csv` | Downloaded | UPX contrast definitions. |
| `GLDS-462_rna_seq_contrasts_mRNA_GLbulkRNAseq.csv` | Downloaded | mRNA contrast definitions. |
| `GLDS-462_rna_seq_contrasts_totRNA_GLbulkRNAseq.csv` | Downloaded | total-RNA contrast definitions. |
| `GLDS-462_proteomics_2021-12-31_tc884-885_Protein_WorkUp.xlsx` | Downloaded | Processed TMT proteomics workbook. Main sheet `protein_quant_2721` has 7,862 protein rows and 66 columns; peptide-level sheets are also included. |
| `GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx` | Downloaded | Processed phosphoproteomics workbook. Main sheets include `siteQuant_360` with 17,251 phosphosite rows and `siteQuant_360_compositeSite` with 3,957 composite-site rows. |
| `GLDS-462_proteomics_proteomics-plex1.tar.gz` | Downloaded but not needed for first pass | Raw vendor mass-spec archive containing `.raw` files. Keep only if we plan to reprocess proteomics from raw files. |

## Still Pending

No additional processed table is currently required for the first-pass RNA/protein/phosphoprotein analysis.

The other missing proteomics/phosphoproteomics archives are large raw-data tarballs. They are not needed unless we deliberately choose to redo the mass-spec processing from vendor `.raw` files.

## First-Pass Analysis Intent

1. Use `GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv` with `GLDS-462_rna_seq_SampleTable_mRNA_GLbulkRNAseq.csv` for expression-level signatures.
2. Use `GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv` for OSDR-provided flight, ground, basal, and vivarium contrasts.
3. Use `GLDS-462_proteomics_2021-12-31_tc884-885_Protein_WorkUp.xlsx`, especially `protein_quant_2721`, for protein-level validation of transcript-level signals.
4. Use `GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx`, especially `siteQuant_360` and `siteQuant_360_compositeSite`, to check whether candidate transporter/signaling changes are visible at the phosphorylation layer.
5. Treat raw mass-spec tarballs as optional archival inputs, not as required inputs for the first analysis pass.

## Sources

- OSDR study page: https://osdr.nasa.gov/bio/repo/data/studies/OSD-462
- OSDR download endpoint pattern: `https://osdr.nasa.gov/geode-py/ws/studies/OSD-462/download?source=datamanager&file=<filename>`
