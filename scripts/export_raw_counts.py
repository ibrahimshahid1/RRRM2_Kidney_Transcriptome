
import scanpy as sc
import scipy.io
import pandas as pd
import os

H5AD_PATH = "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"
OUT_DIR = "data/processed/deconvolution/raw_ref"

def main():
    print(f"Loading {H5AD_PATH}...")
    adata = sc.read_h5ad(H5AD_PATH)
    
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
        
    # Check if raw exists
    if adata.raw is None:
        print("Error: No .raw found in H5AD.")
        return

    print("Extracting raw counts...")
    raw_adata = adata.raw.to_adata()
    
    # Save as mtx
    print(f"Saving to {OUT_DIR}...")
    scipy.io.mmwrite(os.path.join(OUT_DIR, "matrix.mtx"), raw_adata.X)
    
    # Save features and barcodes
    raw_adata.var.to_csv(os.path.join(OUT_DIR, "features.tsv"), sep="\t", header=True)
    raw_adata.obs.to_csv(os.path.join(OUT_DIR, "barcodes.tsv"), sep="\t", header=True)
    
    # Also overwrite the metadata csv with the raw indices to be sure
    # (Though barcodes should match if no filtering happened, but raw might have diff # of genes)
    print("Done exporting raw counts.")

if __name__ == "__main__":
    main()
