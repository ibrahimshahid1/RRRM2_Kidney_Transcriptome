
import scanpy as sc
import pandas as pd
import numpy as np

# Path to the 10x female kidney H5AD found in the workspace
H5AD_PATH = "data/external/single_cell_atlases/kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad"

def inspect_reference(path):
    print(f"Loading {path}...")
    try:
        adata = sc.read_h5ad(path)
    except Exception as e:
        print(f"Error loading file: {e}")
        return

    print("\n--- Obs Columns ---")
    print(adata.obs.columns.tolist())

    # Look for likely cell type columns
    obs_cols = adata.obs.columns
    cell_type_cols = [c for c in obs_cols if any(x in c.lower() for x in ['type', 'annot', 'label', 'class'])]
    
    print("\n--- Likely Cell Type Columns ---")
    print(cell_type_cols)

    for col in cell_type_cols:
        print(f"\nValues in '{col}':")
        unique_vals = adata.obs[col].unique().tolist()
        # Convert to string to safely check
        unique_vals_str = [str(x) for x in unique_vals]
        
        # Check for podocyte/glomerular
        podo_matches = [x for x in unique_vals_str if 'podo' in x.lower()]
        glom_matches = [x for x in unique_vals_str if 'glom' in x.lower()]
        
        print(f"  Total unique values: {len(unique_vals)}")
        print(f"  Examples: {unique_vals[:10]}")
        print(f"  Contains 'podo'?: {podo_matches}")
        print(f"  Contains 'glom'?: {glom_matches}")

    print("\n--- Obs Keys with 'podo' or 'glom' in value ---")
    # Brute force search in all columns if explicitly requested, but cell_type_cols should cover it.
    
if __name__ == "__main__":
    inspect_reference(H5AD_PATH)
