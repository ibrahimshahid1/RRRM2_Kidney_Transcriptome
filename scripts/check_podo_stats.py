
import pandas as pd

try:
    df = pd.read_csv("data/processed/deconvolution/music_segment_proportions.csv")
    print("Columns:", df.columns.tolist())
    if "Podocyte" in df.columns:
        print("Max Podocyte prop:", df["Podocyte"].max())
        print("Mean Podocyte prop:", df["Podocyte"].mean())
        print("Non-zero Podocyte samples:", (df["Podocyte"] > 0).sum())
    else:
        print("Podocyte column missing.")
        
    if "Glomerular" in df.columns: # If aggregated differently
        print("Max Glom:", df["Glomerular"].max())

except Exception as e:
    print(e)
