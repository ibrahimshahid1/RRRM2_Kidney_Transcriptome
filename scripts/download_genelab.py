# scripts/download_genelab.py
import os
import zipfile
import pandas as pd
from urllib.request import urlretrieve
import ssl

# CONFIGURATION
OSD_ID = "771"
GLDS_ID = "674"  # OSD-771 files are prefixed with GLDS-674
# Adjusted path to be relative to where the script is run (inside scripts/)
RAW_DIR = "../data/raw"

# Bypass SSL certification verification (often needed for NASA servers)
ssl._create_default_https_context = ssl._create_unverified_context

def download_file(study_id, filename, output_path):
    """Helper to download a single file from GeneLab."""
    url = f'https://osdr.nasa.gov/geode-py/ws/studies/OSD-{study_id}/download?source=datamanager&file={filename}'
    print(f"Fetching {filename}...")
    try:
        urlretrieve(url, output_path)
        print(f" -> Success: Saved to {output_path}")
    except Exception as e:
        print(f" -> Error downloading {filename}: {e}")

def download_osd_771():
    # Ensure directories exist
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(f"{RAW_DIR}/metadata", exist_ok=True)
    os.makedirs(f"{RAW_DIR}/counts", exist_ok=True)
    
    print(f"Downloading Data for Study OSD-{OSD_ID} (GLDS-{GLDS_ID})...")

    # ---------------------------------------------------------
    # 1. Download Metadata
    # ---------------------------------------------------------
    meta_filename = f"OSD-{OSD_ID}_metadata_OSD-{OSD_ID}-ISA.zip"
    meta_zip_path = f"{RAW_DIR}/metadata/{meta_filename}"
    
    download_file(OSD_ID, meta_filename, meta_zip_path)
    
    # Extract Metadata
    try:
        with zipfile.ZipFile(meta_zip_path, 'r') as zip_ref:
            zip_ref.extractall(f"{RAW_DIR}/metadata")
        print(" -> Metadata extracted successfully.")
    except zipfile.BadZipFile:
        print(" -> Warning: Metadata ZIP appears corrupted or failed to download.")

    # ---------------------------------------------------------
    # 2. Download Unnormalized Counts (For Deconvolution)
    # ---------------------------------------------------------
    # Use RSEM Unnormalized rRNArm counts for accurate cell-type deconvolution
    raw_filename = f"GLDS-{GLDS_ID}_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv"
    raw_path = f"{RAW_DIR}/counts/unnormalized_counts.csv"
    
    download_file(OSD_ID, raw_filename, raw_path)

    # ---------------------------------------------------------
    # 3. Download VST Normalized Counts (For Network Construction)
    # ---------------------------------------------------------
    # Use NASA's VST counts for the main node2vec pipeline
    vst_filename = f"GLDS-{GLDS_ID}_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
    vst_path = f"{RAW_DIR}/counts/vst_counts.csv"
    
    download_file(OSD_ID, vst_filename, vst_path)

    print("\nDownload process complete.")
    print(f"Unnormalized data (for deconvolution): {os.path.abspath(raw_path)}")
    print(f"VST data (for networks): {os.path.abspath(vst_path)}")

if __name__ == "__main__":
    download_osd_771()