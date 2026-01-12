
import scanpy as sc
import numpy as np
import scipy.sparse

H5AD_PATH = "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"

def check_counts():
    print(f"Loading {H5AD_PATH}...")
    try:
        adata = sc.read_h5ad(H5AD_PATH)
    except Exception as e:
        print(f"Error loading: {e}")
        return

    print("\n--- Data Matrix Inspection (adata.X) ---")
    if scipy.sparse.issparse(adata.X):
        data_sample = adata.X.data[:1000]
    else:
        data_sample = adata.X.flatten()[:1000]

    print(f"First 10 values: {data_sample[:10]}")
    print(f"Min: {np.min(data_sample)}")
    print(f"Max: {np.max(data_sample)}")
    
    is_integer = np.all(np.equal(np.mod(data_sample, 1), 0))
    print(f"Are values integers? {is_integer}")
    
    if not is_integer:
        print("Values appear to be floats (likely normalized or log-transformed).")
    else:
        print("Values appear to be integers (likely raw counts).")

    print(f"\n--- Layers available: {list(adata.layers.keys())} ---")
    
    # Check raw if present
    if adata.raw is not None:
        print("\n--- adata.raw Inspection ---")
        try:
            raw_data = adata.raw.X
            if scipy.sparse.issparse(raw_data):
                raw_sample = raw_data.data[:1000]
            else:
                raw_sample = raw_data.flatten()[:1000]
            
            print(f"First 10 values (raw): {raw_sample[:10]}")
            is_raw_integer = np.all(np.equal(np.mod(raw_sample, 1), 0))
            print(f"Are raw values integers? {is_raw_integer}")
        except Exception as e:
            print(f"Could not inspect adata.raw: {e}")
    else:
        print("\nNo adata.raw found.")

if __name__ == "__main__":
    check_counts()
