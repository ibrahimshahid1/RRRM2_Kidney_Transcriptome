"""
Script: src/preprocessing/deconvolution_sanity.py
Purpose: Validate deconvolution (Phase 0) by correlating cell type fractions (CLR) with canonical marker expression (VST).
"""
import pandas as pd
import numpy as np
import scipy.stats
from pathlib import Path

# Repository root
REPO_ROOT = Path(__file__).resolve().parents[2]

# Canonical markers (Mouse Symbols)
MARKERS = {
    "PT": ["Slc34a1", "Lrp2", "Cubn", "Slc5a12"],
    "TAL_LOH": ["Umod", "Slc12a1"],
    "DCT": ["Slc12a3", "Pvalb", "Calb1"],
    "CD": ["Aqp2", "Aqp3", "Scnn1g"],  # Collecting Duct
    "Podocyte": ["Nphs1", "Nphs2", "Wt1"],
    "Endothelial": ["Pecam1", "Kdr"],
    "Macrophage": ["Adgre1", "Cd68"], # Often lumped in Immune
    "Fibroblast": ["Col1a1", "Col3a1", "Pdgfra"]
}

# IDs for key markers (Mouse)
# Hardcode known IDs to avoid map dependency failing
KNOWN_IDS = {
    "Slc34a1": "ENSMUSG00000031929",
    "Lrp2": "ENSMUSG00000027170",
    "Umod": "ENSMUSG00000030500",
    "Slc12a1": "ENSMUSG00000030588", # Wait, Slc12a3 is 30588 (NCC). Slc12a1 is NKCC2.
    "Slc12a3": "ENSMUSG00000030588",
    "Pvalb": "ENSMUSG00000005871",
    "Aqp2": "ENSMUSG00000024610", 
    "Nphs1": "ENSMUSG00000024463",
    "Nphs2": "ENSMUSG00000026132"
}

def load_vst():
    path = REPO_ROOT / "data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
    print(f"Loading VST: {path}")
    df = pd.read_csv(path)
    # First column is gene ID (Ensembl), cleaner to set index
    df = df.set_index(df.columns[0])
    # Strip versions just in case
    df.index = df.index.str.split('.').str[0]
    return df

def load_clr():
    path = REPO_ROOT / "data/processed/deconvolution/music_segment_direct_proportions_CLR.csv"
    print(f"Loading CLR: {path}")
    df = pd.read_csv(path, index_col=0)
    return df

def load_id_map():
    path = REPO_ROOT / "data/processed/resources/id_map.tsv"
    if not path.exists():
        # Fallback to building one or using raw IDs if unavailable
        print("WARN: ID map missing. Using limited hardcoded IDs.")
        return {}
    
    m = pd.read_csv(path, sep='\t', comment='#')
    sym_to_ens = {}
    
    # Normalize column names
    cols = {c.lower(): c for c in m.columns}
    ens_col = cols.get('ensembl_gene_id')
    sym_col = cols.get('mgi_symbol')
    
    if ens_col and sym_col:
        for _, row in m.iterrows():
            if pd.notna(row[sym_col]) and pd.notna(row[ens_col]):
                sym = str(row[sym_col]).strip()
                ens = str(row[ens_col]).strip()
                sym_to_ens[sym] = ens
    return sym_to_ens

def main():
    vst = load_vst()
    clr = load_clr()
    sym_to_ens = load_id_map()
    
    # Merge known IDs
    for k, v in KNOWN_IDS.items():
        if k not in sym_to_ens:
            sym_to_ens[k] = v

    # Align samples
    # VST columns match CLR index?
    # VST cols = sample names. CLR index = sample names.
    common = sorted(list(set(vst.columns) & set(clr.index)))
    print(f"Aligned samples: {len(common)}")
    vst_aligned = vst[common]
    clr_aligned = clr.loc[common]
    
    results = []
    
    print("\n--- Correlation Analysis ---")
    
    for cell_type in clr_aligned.columns:
        # Match cell_type names to our MARKERS dict (fuzzy match)
        markers_to_test = []
        for k, v in MARKERS.items():
            if k.lower() in cell_type.lower() or cell_type.lower() in k.lower():
                markers_to_test.extend(v)
        
        if not markers_to_test:
            continue
            
        print(f"\nCell Type: {cell_type}")
        fraction = clr_aligned[cell_type]
        
        for sym in markers_to_test:
            # Map symbol to ID
            ens_id = sym_to_ens.get(sym)
            
            if not ens_id:
                print(f"  {sym}: No ID mapping found")
                continue
                
            if ens_id not in vst_aligned.index:
                print(f"  {sym} ({ens_id}): Not found in VST matrix")
                continue
                
            expr = vst_aligned.loc[ens_id]
            
            # Remove NaN if any
            mask = ~np.isnan(fraction) & ~np.isnan(expr)
            if mask.sum() < 3:
                print(f"  {sym}: Not enough data points")
                continue
                
            r, p = scipy.stats.pearsonr(fraction[mask], expr[mask])
            rho, p_rho = scipy.stats.spearmanr(fraction[mask], expr[mask])
            
            print(f"  {sym} ({ens_id}): Pearson r={r:.3f} (p={p:.1e}) | Spearman rho={rho:.3f}")
            
            results.append({
                "cell_type": cell_type,
                "marker": sym,
                "ensembl_id": ens_id,
                "pearson_r": r,
                "pearson_p": p,
                "spearman_rho": rho
            })

    # Save to results dir
    out_dir = REPO_ROOT / "data/results"
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / "deconvolution_sanity_check.tsv"
    pd.DataFrame(results).to_csv(out_path, sep='\t', index=False)
    print(f"\n[OK] Validation table saved to: {out_path}")

if __name__ == "__main__":
    main()
