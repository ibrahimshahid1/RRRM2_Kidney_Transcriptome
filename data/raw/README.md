# Raw Data Directory

## Data Sources

### NASA GeneLab OSD-771
**Study**: Rodent Research Reference Mission 2 (RRRM-2)  
**Accession**: OSD-771  
**URL**: https://genelab.nasa.gov/data/study?acc=OSD-771

### Required Files
Download the following from NASA GeneLab:

1. **Count Matrix** (`counts/`)
   - `GLDS-771_rna_seq_ERCCnorm_counts.csv`
   - Raw read counts normalized to ERCC spike-ins

2. **Metadata** (`metadata/`)
   - `GLDS-771_metadata_OSD-771-samples.csv`
   - Sample annotations (age, flight status, batch info)

3. **Gene Annotations** (`metadata/`)
   - `GLDS-771_annotations.csv`
   - Gene symbol, Ensembl ID, chromosome locations

### Download Instructions

#### Option 1: Manual Download
1. Visit https://genelab.nasa.gov/data/study?acc=OSD-771
2. Navigate to "Processed Data"
3. Download files listed above
4. Place in appropriate subdirectories

#### Option 2: Automated Script
```bash
python scripts/download_genelab.py
```

## Data Structure After Download
```
raw/
├── counts/
│   └── GLDS-771_rna_seq_ERCCnorm_counts.csv
├── metadata/
│   ├── GLDS-771_metadata_OSD-771-samples.csv
│   └── GLDS-771_annotations.csv
└── README.md (this file)
```

## Data Description

- **Organism**: Mus musculus (C57BL/6NTac)
- **Tissue**: Whole kidney
- **Sequencing**: Illumina HiSeq 2500, paired-end
- **Read length**: 2×150 bp
- **Samples**: 80 (4 conditions × 10 biological replicates × 2 technical replicates)
- **Genes**: ~25,000 (before filtering)

## Citation

NASA GeneLab Consortium. (2024). Rodent Research Reference Mission 2 (RRRM-2) Kidney Transcriptome Dataset. OSD-771.