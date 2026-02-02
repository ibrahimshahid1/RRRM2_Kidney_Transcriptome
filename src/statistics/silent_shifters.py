# src/statistics/silent_shifters.py
"""
Phase 5: Silent shifters (STRICT prereg definition).

Silent shifter = high rewiring (top decile) AND low DE:
  |log2FC| < 0.3 AND FDR > 0.2

Optionally adds Phase 6 "supported" tag using permutation q-values.

Inputs:
  - rewiring agg tables (Phase 3.3)
  - optional: gene DE table for the matching contrast
  - optional: phase6 perm pval table

Outputs:
  data/results/phase5_silent_shifters_strict/<contrast>_silent_shifters.tsv
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

# Repository root (2 levels up from src/statistics/)
REPO_ROOT = Path(__file__).resolve().parents[2]


def load_rewiring(path: Path) -> pd.DataFrame:
    """Load rewiring aggregate table, normalizing column names."""
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        raise ValueError(f"{path} missing 'gene' column")
    
    # Normalize mean/std columns
    if "rewiring_mean" not in df.columns:
        for c in ["delta_mean", "mean", "rewiring"]:
            if c in df.columns:
                df = df.rename(columns={c: "rewiring_mean"})
                break
    if "rewiring_std" not in df.columns:
        for c in ["delta_std", "std", "sd"]:
            if c in df.columns:
                df = df.rename(columns={c: "rewiring_std"})
                break
    return df


def main():
    ap = argparse.ArgumentParser(description="Build strict silent shifters")
    ap.add_argument("--rewiring", required=True, 
                    help="*_rewiring_agg.tsv from phase3")
    ap.add_argument("--gene_de", default="", 
                    help="Optional gene-level DE TSV with gene, log2FC/logFC, FDR/adj.P.Val")
    ap.add_argument("--perm", default="", 
                    help="Optional phase6 perm table with gene, q_BH")
    ap.add_argument("--outdir", 
                    default=str(REPO_ROOT / "data/results/phase5_silent_shifters_strict"),
                    help="Output directory")
    ap.add_argument("--top_quantile", type=float, default=0.9, 
                    help="Top decile = 0.9")
    ap.add_argument("--fdr_min", type=float, default=0.2,
                    help="Minimum FDR for 'not DE' (default: 0.2)")
    ap.add_argument("--lfc_max", type=float, default=0.3,
                    help="Maximum |logFC| for 'not DE' (default: 0.3)")
    ap.add_argument("--support_q", type=float, default=0.1,
                    help="q-value threshold for 'supported' (default: 0.1)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load rewiring
    rw_path = Path(args.rewiring)
    rw = load_rewiring(rw_path)
    print(f"Loaded {len(rw)} genes from {rw_path}")

    # High rewiring threshold
    thr = rw["rewiring_mean"].quantile(args.top_quantile)
    rw["high_rewiring"] = rw["rewiring_mean"] >= thr
    print(f"High rewiring threshold (q={args.top_quantile}): {thr:.4f}")
    print(f"  High rewiring genes: {rw['high_rewiring'].sum()}")

    # Load DE if provided
    if args.gene_de and Path(args.gene_de).exists():
        de = pd.read_csv(args.gene_de, sep=None, engine="python")
        de = de.rename(columns={
            "adj.P.Val": "FDR", "padj": "FDR", "qvalue": "FDR",
            "log2FC": "logFC"
        })
        
        if "FDR" not in de.columns or "logFC" not in de.columns or "gene" not in de.columns:
            print(f"[WARN] gene_de missing required columns. Found: {list(de.columns)}")
            de = None
        else:
            de = de[["gene", "logFC", "FDR"]].drop_duplicates("gene")
            rw = rw.merge(de, on="gene", how="left")
            rw["low_DE"] = (rw["FDR"] > args.fdr_min) & (rw["logFC"].abs() < args.lfc_max)
            print(f"Loaded DE for {de['gene'].nunique()} genes")
            print(f"  Low DE genes: {rw['low_DE'].sum()}")
    else:
        # No DE available - use rewiring only
        rw["logFC"] = np.nan
        rw["FDR"] = np.nan
        rw["low_DE"] = True  # All genes pass DE filter when not available
        print("[INFO] No gene_de provided, using rewiring-only filter")

    # Load permutation support if provided
    if args.perm and Path(args.perm).exists():
        perm = pd.read_csv(args.perm, sep="\t")
        if "q_BH" in perm.columns:
            perm = perm[["gene", "q_BH"]].drop_duplicates("gene")
            rw = rw.merge(perm, on="gene", how="left")
            rw["supported"] = rw["q_BH"].fillna(1.0) < args.support_q
            print(f"Loaded permutation q-values for {perm['gene'].nunique()} genes")
            print(f"  Supported genes (q<{args.support_q}): {rw['supported'].sum()}")
        else:
            rw["q_BH"] = np.nan
            rw["supported"] = np.nan
    else:
        rw["q_BH"] = np.nan
        rw["supported"] = np.nan
        print("[INFO] No perm file provided, 'supported' column will be NA")

    # Filter to silent shifters
    silent = rw[rw["high_rewiring"] & rw["low_DE"]].copy()
    silent = silent.sort_values("rewiring_mean", ascending=False).reset_index(drop=True)

    # Output filename based on input
    out_name = rw_path.stem.replace("_rewiring_agg", "") + "_silent_shifters.tsv"
    out_path = outdir / out_name
    silent.to_csv(out_path, sep="\t", index=False)

    print(f"\nSilent shifter summary:")
    print(f"  Total genes: {len(rw)}")
    print(f"  High rewiring: {rw['high_rewiring'].sum()}")
    print(f"  Low DE: {rw['low_DE'].sum()}")
    print(f"  Silent shifters: {len(silent)}")
    if silent["supported"].notna().any():
        print(f"  Silent + supported: {(silent['supported'] == True).sum()}")

    print(f"\n[OK] Wrote: {out_path}")


if __name__ == "__main__":
    main()
