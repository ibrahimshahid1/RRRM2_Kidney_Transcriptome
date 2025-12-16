# scripts/download_genelab.py
import os
import zipfile
import pandas as pd
from urllib.request import urlretrieve

# CONFIGURATION
OSD_ID = "771"
RAW_DIR = "../data/raw"

def download_osd_771():
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(f"{RAW_DIR}/metadata", exist_ok=True)
    os.makedirs(f"{RAW_DIR}/counts", exist_ok=True)
    
    print(f"Downloading OSD-{OSD_ID} Data...")

    # 1. Download Metadata (Adapted from read_meta_data)
    meta_url = f'https://osdr.nasa.gov/geode-py/ws/studies/OSD-{OSD_ID}/download?source=datamanager&file=OSD-{OSD_ID}_metadata_OSD-{OSD_ID}-ISA.zip'
    meta_zip_path = f"{RAW_DIR}/metadata/OSD-{OSD_ID}_meta.zip"
    
    print(f"Fetching Metadata from: {meta_url}")
    urlretrieve(meta_url, meta_zip_path)
    
    # Extract
    with zipfile.ZipFile(meta_zip_path, 'r') as zip_ref:
        zip_ref.extractall(f"{RAW_DIR}/metadata")
    print("Metadata downloaded and extracted.")

    # 2. Download Counts (Adapted from read_rnaseq_data)
    # Note: We need the specific filename. Usually for OSD-771 it follows this pattern:
    # Check the repo manually or assume standard naming. 
    # If this fails, you might need to grab the exact filename from the website.
    counts_filename = f"GLDS-{OSD_ID}_rna_seq_Unnormalized_Counts.csv"
    counts_url = f'https://osdr.nasa.gov/geode-py/ws/studies/OSD-{OSD_ID}/download?source=datamanager&file={counts_filename}'
    
    counts_path = f"{RAW_DIR}/counts/OSD-{OSD_ID}_counts.csv"
    
    print(f"Fetching Counts from: {counts_url}")
    try:
        urlretrieve(counts_url, counts_path)
        print("Counts downloaded successfully.")
    except Exception as e:
        print(f"Error downloading counts: {e}")
        print("You may need to check the exact filename on the GeneLab OSD-771 page.")

if __name__ == "__main__":
    download_osd_771()