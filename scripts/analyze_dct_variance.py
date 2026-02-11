
import pandas as pd
import numpy as np
import os
import json
from pathlib import Path

# Paths
REPO_ROOT = Path("/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome")
COUNTS_PATH = REPO_ROOT / "data/processed/aligned_outputs/rsem_rRNArm_raw_counts.csv"
DECONV_PATH = REPO_ROOT / "data/processed/deconvolution/music_segment_direct_proportions_CLR.csv"
ID_MAP_PATH = REPO_ROOT / "data/processed/resources/id_map.json"

# Canonical Markers
MARKERS = {
    "DCT1 (NCC)": ["Slc12a3"],
    "DCT2/CNT": ["Pvalb", "Calb1"],
    "Principal (PC)": ["Aqp2"],
    "Intercalated (IC)": ["Atp6v1b1", "Slc4a1", "Slc26a4"],
    "Reference": ["Actb", "Gapdh"]
}

def load_id_map(path):
    print(f"Loading ID map from {path}...")
    with open(path) as f:
        return json.load(f)

def load_counts(path):
    print(f"Loading counts from {path}...")
    # First column is gene ID
    df = pd.read_csv(path, index_col=0)
    # Strip Ensembl versions
    df.index = df.index.str.replace(r"\.\d+$", "", regex=True)
    # Collapse duplicates (if any)
    df = df.groupby(df.index).sum()
    return df

def normalize_cpm(counts):
    print("Normalizing to CPM...")
    # Sum per sample
    libs = counts.sum(axis=0)
    # CPM = (counts / lib_size) * 1e6
    cpm = (counts * 1e6) / libs
    # Log1p(CPM) for variance stability
    return np.log1p(cpm)

def main():
    if not COUNTS_PATH.exists():
        print(f"Error: Counts file not found at {COUNTS_PATH}")
        return

    # 1. Load Data
    counts = load_counts(COUNTS_PATH)
    
    # 2. Map IDs to Symbols
    if ID_MAP_PATH.exists():
        id_map = load_id_map(ID_MAP_PATH)
        
        # Hardcode Slc12a3 if missing from map but present in counts
        slc12a3_id = "ENSMUSG00000031766"
        if slc12a3_id in counts.index:
            print(f"Manually mapping {slc12a3_id} -> Slc12a3")
            id_map[slc12a3_id] = "Slc12a3"
        else:
             print(f"Warning: Slc12a3 ID {slc12a3_id} NOT found in counts index.")

        # Add Symbol column
        counts["Symbol"] = counts.index.map(id_map)
        
        # Aggregate by symbol
        counts_symbol = counts.dropna(subset=["Symbol"]).groupby("Symbol").sum()
        
        print(f"Collapsed genes to {len(counts_symbol)} unique symbols.")
        counts = counts_symbol
    else:
        print("Warning: ID map not found. IDs may be Ensembl.")

    # 3. Normalize
    log_cpm = normalize_cpm(counts)
    
    # 4. Extract Marker Data
    print("\n--- Marker Gene Variance Analysis (Log1p CPM) ---")
    
    stats_rows = []
    
    flat_markers = [g for sublist in MARKERS.values() for g in sublist]
    present_markers = [g for g in flat_markers if g in log_cpm.index]
    
    if not present_markers:
        print(f"No markers found! First 5 count genes: {log_cpm.index[:5].tolist()}")
        return

    # Calculate stats
    marker_data = log_cpm.loc[present_markers].T # Samples x Genes
    
    for group, genes in MARKERS.items():
        for gene in genes:
            if gene not in marker_data.columns:
                continue
            
            vec = marker_data[gene]
            mean_val = vec.mean()
            var_val = vec.std()**2 # safe variance
            cv = vec.std() / mean_val if mean_val > 0 else 0
            
            stats_rows.append({
                "Group": group,
                "Gene": gene,
                "Mean_Expr": mean_val,
                "Variance": var_val,
                "CV": cv
            })
            
    stats_df = pd.DataFrame(stats_rows).sort_values("Variance", ascending=False)
    
    print(stats_df.to_string(index=False, float_format="{:.4f}".format))
    
    # 5. Correlate with 'DCT' Fraction
    if DECONV_PATH.exists():
        print("\n--- Correlation with Deconvolved 'DCT' Fraction ---")
        try:
            deconv = pd.read_csv(DECONV_PATH, index_col=0)
            if "DCT" not in deconv.columns:
                print("Error: 'DCT' column missing in deconvolution results.")
            else:
                dct_frac = deconv["DCT"]
                
                intersect = deconv.index.intersection(marker_data.index)
                if len(intersect) < 10:
                    print(f"Sample mismatch! Intersection: {len(intersect)}")
                    print(f"Deconv: {deconv.index[:3]}")
                    print(f"Counts: {marker_data.index[:3]}")
                    return
                
                dct_vec = dct_frac.loc[intersect]
                expr_aligned = marker_data.loc[intersect]
                
                corrs = []
                for gene in expr_aligned.columns:
                    r = dct_vec.corr(expr_aligned[gene])
                    corrs.append({"Gene": gene, "Corr_with_DCT_Fract": r})
                    
                corr_df = pd.DataFrame(corrs).sort_values("Corr_with_DCT_Fract", ascending=False)
                
                # Add group info safely
                corr_df["Group"] = corr_df["Gene"].apply(lambda g: next((k for k, v in MARKERS.items() if g in v), "Unknown"))
                
                print(corr_df[["Group", "Gene", "Corr_with_DCT_Fract"]].to_string(index=False, float_format="{:.4f}".format))
                
        except Exception as e:
            print(f"Error correlation analysis: {e}")
    else:
        print(f"\nDeconvolution file not found at {DECONV_PATH}. Skipping correlation.")

if __name__ == "__main__":
    main()
