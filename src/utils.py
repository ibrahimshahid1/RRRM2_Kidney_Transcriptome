import pandas as pd
import os

def load_raw_data(data_dir):
    """
    Loads the raw counts and metadata from the specified directory.
    Assumes standard NASA OSDR naming conventions or your specific filenames.
    
    Args:
        data_dir (str): Path to 'data/raw'
        
    Returns:
        tuple: (counts_df, metadata_df)
    """
    # 1. Define paths (Adjust filenames to match your exact download)
    counts_path = os.path.join(data_dir, "counts", "OSD-771_counts.csv") 
    meta_path = os.path.join(data_dir, "metadata", "OSD-771_metadata.csv")
    
    # 2. Check if files exist
    if not os.path.exists(counts_path):
        raise FileNotFoundError(f"Counts file not found at: {counts_path}")
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Metadata file not found at: {meta_path}")

    # 3. Load Data
    # Index_col=0 assumes the first column is Gene ID / Sample Name
    print(f"Loading counts from {counts_path}...")
    counts = pd.read_csv(counts_path, index_col=0)
    
    print(f"Loading metadata from {meta_path}...")
    meta = pd.read_csv(meta_path, index_col=0)
    
    return counts, meta