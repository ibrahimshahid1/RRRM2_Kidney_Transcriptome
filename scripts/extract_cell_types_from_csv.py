
import pandas as pd

CSV_PATH = "data/external/single_cell_atlases/obs_metadata_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.csv"

def main():
    print(f"Reading {CSV_PATH}...")
    df = pd.read_csv(CSV_PATH)
    
    print("\nColumns:", df.columns.tolist())
    
    for col in ['cell_type', 'free_annotation', 'subtissue']:
        if col in df.columns:
            print(f"\n--- Unique values in '{col}' ---")
            vals = sorted(df[col].dropna().unique().astype(str))
            for v in vals:
                print(f"  {v}")
            
            podo = [v for v in vals if 'podo' in v.lower()]
            glom = [v for v in vals if 'glom' in v.lower()]
            print(f"  Includes 'podo'?: {podo}")
            print(f"  Includes 'glom'?: {glom}")

if __name__ == "__main__":
    main()
