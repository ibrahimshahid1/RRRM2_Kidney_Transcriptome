# src/preprocessing/normalization.py
import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

def filter_low_expression(counts_df, threshold_tpm=1, min_fraction=0.2):
    """
    Filters genes with TPM < 1 in > 80% of samples (keeps genes present in > 20%).
    Note: If lengths are unavailable, we approximate TPM with CPM (Counts Per Million).
    
    Args:
        counts_df (pd.DataFrame): Raw counts (genes x samples).
    """
    # Calculate CPM (approximation for TPM if gene lengths unknown)
    cpm = counts_df * 1e6 / counts_df.sum(axis=0)
    
    # Logic: Keep gene if (TPM >= threshold) in at least (min_fraction * n_samples)
    mask = (cpm >= threshold_tpm).sum(axis=1) >= (min_fraction * counts_df.shape[1])
    
    print(f"Filtering: Dropped {np.sum(~mask)} genes. Keeping {np.sum(mask)} genes.")
    return counts_df[mask]

def apply_vst(counts_df, metadata_df, design_factors="Age + Flight"):
    """
    Applies Variance Stabilizing Transformation (VST) using PyDESeq2.
    Matches PDF Requirement: [cite: 43]
    """
    print("Initializing DESeq2 Object...")
    # PyDESeq2 expects Samples x Genes, so we transpose
    dds = DeseqDataSet(
        counts=counts_df.T,
        metadata=metadata_df,
        design_factors=design_factors
    )
    
    print("Running VST normalization...")
    dds.vst()
    
    # Return Genes x Samples
    return pd.DataFrame(
        dds.layers['vst_counts'].T, 
        index=counts_df.index, 
        columns=counts_df.columns
    )