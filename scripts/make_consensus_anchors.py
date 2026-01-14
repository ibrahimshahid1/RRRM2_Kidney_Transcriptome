from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def load_one(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns or "rewiring_abs" not in df.columns:
        raise ValueError(f"{path} must contain columns: gene, rewiring_abs")
    df = df[["gene", "rewiring_abs"]].copy()
    return df

def consensus(files: list[str], k: int, out: str):
    dfs = []
    for f in files:
        df = load_one(f)
        # rank within this file (1=best)
        df["rank"] = df["rewiring_abs"].rank(method="average", ascending=True)
        df["pct"]  = df["rank"] / df["rank"].max()
        tag = Path(f).stem
        df = df.rename(columns={
            "rewiring_abs": f"abs__{tag}",
            "rank":        f"rank__{tag}",
            "pct":         f"pct__{tag}",
        })
        dfs.append(df)

    # outer-merge on gene (in case sets differ)
    m = dfs[0]
    for df in dfs[1:]:
        m = m.merge(df, on="gene", how="outer")

    pct_cols  = [c for c in m.columns if c.startswith("pct__")]
    rank_cols = [c for c in m.columns if c.startswith("rank__")]
    abs_cols  = [c for c in m.columns if c.startswith("abs__")]

    # robustness stats across comparisons
    m["median_pct"] = m[pct_cols].median(axis=1, skipna=True)
    m["max_pct"]    = m[pct_cols].max(axis=1, skipna=True)
    m["iqr_pct"]    = (m[pct_cols].quantile(0.75, axis=1) - m[pct_cols].quantile(0.25, axis=1))

    # if a gene is missing in some file, penalize it mildly
    m["n_present"]  = m[pct_cols].notna().sum(axis=1)
    m["missing"]    = len(files) - m["n_present"]

    # score: lower is better
    # - median drives stability
    # - max punishes “good usually, bad once”
    # - iqr punishes inconsistency
    # - missing penalizes absent genes
    m["score"] = (
        m["median_pct"]
        + 0.50 * m["max_pct"]
        + 0.25 * m["iqr_pct"]
        + 0.10 * (m["missing"] / len(files))
    )

    out_df = m.sort_values("score", ascending=True).head(k).copy()

    # add a few helpful summary columns in raw units too
    out_df["median_abs"] = m[abs_cols].median(axis=1, skipna=True)
    out_df["max_abs"]    = m[abs_cols].max(axis=1, skipna=True)

    outpath = Path(out)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(outpath, sep="\t", index=False)

    # quick overlap stats
    sets = []
    for f in files:
        sets.append(set(load_one(f)["gene"]))
    inter = set.intersection(*sets)
    union = set.union(*sets)
    print(f"Files: {len(files)}")
    print(f"Intersection size: {len(inter)}")
    print(f"Union size: {len(union)}")
    print(f"Jaccard: {len(inter)/len(union):.4f}")
    print(f"Saved consensus anchors (k={k}) -> {outpath}")
    print(out_df[["gene","score","median_pct","max_pct","iqr_pct","n_present","median_abs","max_abs"]].head(15).to_string(index=False))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--files", nargs="+", required=True)
    ap.add_argument("--k", type=int, default=150)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    consensus(args.files, args.k, args.out)

if __name__ == "__main__":
    main()
