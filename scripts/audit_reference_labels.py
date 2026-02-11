
import sys
import pandas as pd
import re
from pathlib import Path
import os

# We will try to load the converted CSV/TSV files first, as they are faster and guaranteed to be readable
# h5ad loading can be finicky depending on environment (scanpy/anndata versions)
RAW_REF_DIR = Path("/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome/data/processed/deconvolution/raw_ref")
BARCOCES_PATH = RAW_REF_DIR / "barcodes.tsv"
FEATURES_PATH = RAW_REF_DIR / "features.tsv"
MTX_PATH = RAW_REF_DIR / "matrix.mtx"

# Original H5AD
H5AD_PATH = Path("/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome/data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad")

def segment_from_ct(x):
    if pd.isna(x):
        return "NA"
    x2 = str(x).lower().strip()

    # The exact regex logic from src/preprocessing/deconvolution.R
    if re.search(r"podocyte", x2):
        return "Podocyte"
    if re.search(r"mesangial", x2):
        return "Mesangial"
    if re.search(r"pecam|endothelial|capillary|artery", x2):
        return "Endothelial"
    if re.search(r"fibroblast|stroma|stromal", x2):
        return "Fibroblast"
    if re.search(r"^cd45|macrophage|\bt\s*cell\b|\bb\s*cell\b|plasma\s*cell|nk\s*cell|lymph|leukocyte", x2):
        return "Immune"
    
    # Tubules
    if re.search(r"proximal.*tubule|proximal\s+tube|\bpt\b", x2):
        return "PT"
    
    # DCT regex is CRITICAL here
    if re.search(r"distal.*convoluted|distal.*tubule|\bdct\b", x2):
        return "DCT"
        
    if re.search(r"thick.*ascending|loop\s+of\s+henle|henle|\btal\b", x2):
        return "TAL_LOH"
    if re.search(r"collecting\s+duct|principal\s+cell|intercalated", x2):
        return "CD"
    
    if re.search(r"brush\s*cell|tuft\s*cell", x2):
        return "Other"
    
    return "Other (Unmapped)"

def main():
    print(f"Auditing Reference Labels...")
    
    # Prefer H5AD if available for complete metadata
    try:
        import anndata
        print(f"Loading H5AD: {H5AD_PATH}")
        adata = anndata.read_h5ad(H5AD_PATH)
        obs = adata.obs
        print(f"Loaded {len(obs)} cells from H5AD.")
    except Exception as e:
        print(f"Could not load H5AD ({e}). Trying raw_ref TSV export...")
        if BARCOCES_PATH.exists():
            obs = pd.read_csv(BARCOCES_PATH, sep="\t", index_col=0)
            print(f"Loaded {len(obs)} cells from barcodes.tsv.")
        else:
            print("No reference metadata found.")
            return

    # Columns to check
    candidates = ["free_annotation", "cell_type", "subtissue", "cluster"]
    best_col = None
    
    print("\n--- Available Columns ---")
    for col in candidates:
        if col in obs.columns:
            n_unique = obs[col].nunique()
            print(f"'{col}': {n_unique} unique labels")
            if col == "free_annotation": 
                best_col = col

    if not best_col:
        # Fallback to cell_type if free_annotation missing
        if "cell_type" in obs.columns: best_col = "cell_type"
        else: best_col = obs.columns[0]
            
    print(f"\nUsing primary label column: '{best_col}'")
    
    labels = obs[best_col].astype(str).value_counts()
    
    print(f"\n{'='*90}")
    print(f"{'MAPPED SEGMENT':<20} | {'ORIGINAL LABEL (Top 50)':<50} | {'COUNT':>8}")
    print(f"{'-'*20} | {'-'*50} | {'-'*8}")
    
    mapped_counts = {}
    
    for label, count in labels.head(50).items():
        seg = segment_from_ct(label)
        print(f"{seg:<20} | {label:<50} | {count:>8}")
        
    print(f"{'='*90}")

    # Exhaustive search for dangerous unmapped terms
    print("\n--- ⚠️  DANGER AUDIT: Missed Distal/Collecting Terms ---")
    print("Looking for labels containing 'connecting', 'cnt', 'intercalated', 'principal'...")
    
    missed_terms = ["connecting", "cnt", "intercalated", "principal", "distal"]
    
    for label in labels.index:
        seg = segment_from_ct(label)
        label_lower = label.lower()
        
        # Check if mapped to "Other" or generic, but has keywords
        if seg in ["Other (Unmapped)", "Other"]:
            for term in missed_terms:
                if term in label_lower:
                    print(f"❌ MISSED: '{label}' contains '{term}' -> Mapped to '{seg}'")
        
        # Check if mapped to DCT but looks like something else (CNT?)
        if seg == "DCT":
             if "connecting" in label_lower or "cnt" in label_lower:
                 print(f"⚠️  AMBIGUOUS: '{label}' mapped to DCT (might be CNT?)")

    print("\n--- Audit Complete ---")

if __name__ == "__main__":
    main()
