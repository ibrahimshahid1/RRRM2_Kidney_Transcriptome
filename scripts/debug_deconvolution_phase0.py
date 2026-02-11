import pandas as pd
import numpy as np
import os
import re

VST_PATH = "/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome/data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
MARKER_MEANS_PATH = "/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome/data/processed/deconvolution/marker_sanity/marker_mean_logcpm_by_segment.csv"
PROPS_PATH = "/Users/ibrahimshahid/Documents/Github/RRRM2_Kidney_Transcriptome/data/processed/deconvolution/music_cluster_proportions.csv"

def load_vst():
    print("\nLoading VST data...")
    try:
        vst = pd.read_csv(VST_PATH, index_col=0)
        # Clean index: remove versions
        vst.index = vst.index.astype(str).str.split(".").str[0]
        print(f"VST loaded: {vst.shape}")
        return vst
    except Exception as e:
        print(f"Failed to load VST: {e}")
        return None

def build_gene_map(marker_file):
    print("Building gene map from marker file...")
    try:
        ref_means = pd.read_csv(marker_file, index_col=0)
        # Index format: "Symbol (EID)"
        gene_map = {}
        for idx in ref_means.index:
            match = re.match(r"(\w+) \((ENS\w+)\)", idx)
            if match:
                sym, eid = match.groups()
                gene_map[sym] = eid
        print(f"Mapped {len(gene_map)} genes.")
        return gene_map, ref_means
    except Exception as e:
        print(f"Failed to build gene map: {e}")
        return {}, None

def step1_check_markers(vst):
    print("\n--- Step 1: Confirm you're not fighting a dumb input bug ---")
    genes = {
      "Slc12a3":"ENSMUSG00000031766",
      "Pvalb":"ENSMUSG00000005716",
      "Calb1":"ENSMUSG00000028222",
      "Slc12a1":"ENSMUSG00000027202",
      "Umod":"ENSMUSG00000030963",
      "Lrp2":"ENSMUSG00000027070",
      "Cubn":"ENSMUSG00000026726"
    }
    
    missing = []
    for sym, eid in genes.items():
        if eid in vst.index:
            print(f"✅ {sym} ({eid}) present")
        else:
            print(f"❌ {sym} ({eid}) MISSING")
            missing.append(sym)
            
    if missing:
        print(f"❌ Step 1 FAILED: Missing markers: {', '.join(missing)}")
        return False
    else:
        print("✅ Step 1 PASSED.")
        return True

def step2_verify_ref(ref_means, gene_map):
    print("\n--- Step 2: Verify what 'DCT' means in your reference ---")
    
    # We need EIDs for Slc12a3 and Pvalb. 
    # Use hardcoded if map fails or doesn't have them (unlikely if they are in ref)
    slc_eid = gene_map.get("Slc12a3", "ENSMUSG00000031766")
    pvalb_eid = gene_map.get("Pvalb", "ENSMUSG00000005716")
    
    print(f"Checking for Slc12a3 ({slc_eid}) and Pvalb ({pvalb_eid})")

    # Find rows by EID in the formatted index
    def get_row(eid):
        if not eid: return None
        for idx in ref_means.index:
            if eid in idx:
                return ref_means.loc[idx]
        return None

    slc_row = get_row(slc_eid)
    pvalb_row = get_row(pvalb_eid)
    
    if slc_row is None:
        print("❌ Slc12a3 not found in reference means.")
    if pvalb_row is None:
        print("❌ Pvalb not found in reference means.")
        
    if slc_row is None or pvalb_row is None:
        return

    # Identify columns
    cols = ref_means.columns
    dct_col = next((c for c in cols if "DCT" in c), None)
    tal_col = next((c for c in cols if "TAL" in c), None)
    
    if not dct_col or not tal_col:
        print(f"Missing DCT or TAL columns. Found: {cols}")
        return

    print(f"Using columns: {dct_col}, {tal_col}")
    
    s_d = slc_row[dct_col]
    s_t = slc_row[tal_col]
    p_d = pvalb_row[dct_col]
    p_t = pvalb_row[tal_col]
    
    print(f"Slc12a3 | DCT: {s_d:.2f} | TAL: {s_t:.2f}")
    print(f"Pvalb   | DCT: {p_d:.2f} | TAL: {p_t:.2f}")
    
    if p_d > s_d:
        print("⚠️  WARNING: Pvalb > Slc12a3 in DCT reference. Label might be DCT2/CNT.")
    else:
        print("✅ Slc12a3 > Pvalb in DCT reference.")

def step3_identifiability(props):
    print("\n--- Step 3: Quantify 'identifiability' ---")
    if "DCT" not in props.columns:
        print("DCT column not in proportions.")
        return

    corrs = props.corr()["DCT"].sort_values()
    print("Correlations with DCT:")
    print(corrs)
    
    # Check for strong correlations excluding self
    others = corrs[corrs.index != "DCT"]
    
    if (others < -0.6).any():
        print("⚠️  Strong negative correlations found (simplex tradeoff).")
    if (others > 0.6).any():
        print("⚠️  Strong positive correlations found (shared program).")

def step4_dynamic_range(props):
    print("\n--- Step 4: Check DCT dynamic range ---")
    if "DCT" not in props.columns: return
    
    dct = props["DCT"]
    desc = dct.describe()[['mean', 'std', 'min', 'max']]
    print(desc)
    
    if dct.std() < 0.01:
        print("⚠️  Very low standard deviation! Signal may be dominated by noise.")

def step5_marker_scores(vst, props, gene_map):
    print("\n--- Step 5: DCT1 vs DCT2 marker scores ---")
    
    dct1_syms = ["Slc12a3", "Wnk1", "Wnk4", "Hsd11b2"]
    dct2_syms = ["Pvalb", "Calb1"]
    
    # Get EIDs from map, fallback to knowing Pvalb/Calb1/Slc12a3
    # Use known EIDs for Wnk1 etc if not in map
    known_eids = {
        "Slc12a3": "ENSMUSG00000031766",
        "Pvalb": "ENSMUSG00000005716",
        "Calb1": "ENSMUSG00000028222"
    }
    
    dct1_eids = []
    for s in dct1_syms:
        eid = gene_map.get(s)
        if not eid and s in known_eids: eid = known_eids[s]
        if eid: dct1_eids.append(eid)
            
    dct2_eids = []
    for s in dct2_syms:
        eid = gene_map.get(s)
        if not eid and s in known_eids: eid = known_eids[s]
        if eid: dct2_eids.append(eid)
    
    # Filter present in VST
    dct1_eids = [e for e in dct1_eids if e in vst.index]
    dct2_eids = [e for e in dct2_eids if e in vst.index]
    
    print(f"DCT1 Genes used: {dct1_eids} (Total {len(dct1_eids)})")
    print(f"DCT2 Genes used: {dct2_eids} (Total {len(dct2_eids)})")
    
    if not dct1_eids or not dct2_eids:
        print("Cannot compute scores, missing genes.")
        return

    # Compute scores (mean of log-counts)
    score1 = vst.loc[dct1_eids].mean(axis=0)
    score2 = vst.loc[dct2_eids].mean(axis=0)
    
    # Intersect samples
    # Check if index matches
    common = props.index.intersection(vst.columns)
    if len(common) == 0:
        print(f"Props index: {props.index[0]}")
        print(f"VST columns: {vst.columns[0]}")
        print("No intersecting samples. Trying to fix index names...")
        # Maybe strip something?
        # If props index has "GLDS-674_..." and vst has "GLDS-674_..." 
        pass
    
    if len(common) > 0:
        y = props.loc[common, "DCT"]
        x1 = score1[common]
        x2 = score2[common]
        
        c1 = y.corr(x1)
        c2 = y.corr(x2)
        
        print(f"Corr(DCT_prop, DCT1_score): {c1:.4f}")
        print(f"Corr(DCT_prop, DCT2_score): {c2:.4f}")
        
        if c2 > c1+0.1 and c2 > 0.4:
             print("⚠️  DCT fraction correlates significantly better with DCT2 markers.")
             print("   Conclusion: Your DCT fraction is functionally 'DCT2/CNT mass'.")
        elif c1 > 0.4:
             print("✅ DCT fraction correlates well with DCT1 markers.")
        else:
             print("⚠️  DCT fraction correlates poorly with both.")

def main():
    vst = load_vst()
    if vst is None: return

    gene_map, ref_means = build_gene_map(MARKER_MEANS_PATH)
    
    if not step1_check_markers(vst):
        return

    if ref_means is not None:
        step2_verify_ref(ref_means, gene_map)
    
    try:
        if os.path.exists(PROPS_PATH):
            props = pd.read_csv(PROPS_PATH, index_col=0)
            step3_identifiability(props)
            step4_dynamic_range(props)
            step5_marker_scores(vst, props, gene_map)
        else:
            print(f"Proportions file not found at {PROPS_PATH}")
    except Exception as e:
        print(f"Failed to load proportions: {e}")

if __name__ == "__main__":
    main()
