
import pandas as pd

CSV_PATH = "data/external/single_cell_atlases/obs_metadata_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.csv"

def main():
    df = pd.read_csv(CSV_PATH)
    col = 'free_annotation' 
    if col in df.columns:
        vals = sorted(df[col].dropna().unique().astype(str))
        with open("unique_cell_types.txt", "w") as f:
            for v in vals:
                f.write(v + "\n")
        print(f"Wrote {len(vals)} unique values to unique_cell_types.txt")
    else:
        print(f"Column '{col}' not found.")

if __name__ == "__main__":
    main()
