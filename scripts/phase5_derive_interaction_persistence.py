# scripts/phase5_derive_interaction_persistence.py
"""
Phase 5 (derived): interaction + persistence/recovery metrics from Phase 3.3 outputs.

Inputs (from Phase 3.3):
  data/results/phase3_rewiring/*_rewiring_agg.tsv

Outputs:
  data/results/phase5_derived/
    ISS_T_interaction.tsv
    LAR_interaction.tsv
    ISS_minus_LAR_YNG_persistence.tsv
    ISS_minus_LAR_OLD_persistence.tsv
"""

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_rewiring_agg(path: Path) -> pd.DataFrame:
    """Load rewiring aggregate table, normalizing column names."""
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        raise ValueError(f"{path} missing 'gene' column. Columns: {list(df.columns)}")

    # Normalize column names across possible variants
    col_map = {}
    # mean
    for c in ["rewiring_mean", "delta_mean", "mean", "rewiring"]:
        if c in df.columns:
            col_map[c] = "rewiring_mean"
            break
    # std
    for c in ["rewiring_std", "delta_std", "std", "sd"]:
        if c in df.columns:
            col_map[c] = "rewiring_std"
            break
    # rank stability (optional)
    for c in ["rank_std", "rank_sd"]:
        if c in df.columns:
            col_map[c] = "rank_std"
            break

    df = df.rename(columns=col_map).copy()
    keep = ["gene", "rewiring_mean"]
    if "rewiring_std" in df.columns:
        keep.append("rewiring_std")
    if "rank_std" in df.columns:
        keep.append("rank_std")
    return df[keep]


def merge_on_gene(a: pd.DataFrame, b: pd.DataFrame, suffixes=("_A", "_B")) -> pd.DataFrame:
    """Merge two dataframes on 'gene' column."""
    return a.merge(b, on="gene", how="inner", suffixes=suffixes)


def main():
    ap = argparse.ArgumentParser(description="Compute interaction + persistence metrics")
    ap.add_argument("--in_dir", default=str(REPO_ROOT / "data/results/phase3_rewiring"),
                    help="Directory containing Phase 3.3 rewiring agg outputs")
    ap.add_argument("--out_dir", default=str(REPO_ROOT / "data/results/phase5_derived"),
                    help="Output directory for derived metrics")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Required Phase 3.3 agg outputs
    files = {
        "ISS_YNG": in_dir / "ISS_T_YNG_FLT_minus_GC_rewiring_agg.tsv",
        "ISS_OLD": in_dir / "ISS_T_OLD_FLT_minus_GC_rewiring_agg.tsv",
        "LAR_YNG": in_dir / "LAR_YNG_FLT_minus_GC_rewiring_agg.tsv",
        "LAR_OLD": in_dir / "LAR_OLD_FLT_minus_GC_rewiring_agg.tsv",
    }
    
    for k, p in files.items():
        if not p.exists():
            raise FileNotFoundError(f"Missing {k}: {p}")

    print(f"Loading rewiring tables from {in_dir}...")
    iss_y = load_rewiring_agg(files["ISS_YNG"])
    iss_o = load_rewiring_agg(files["ISS_OLD"])
    lar_y = load_rewiring_agg(files["LAR_YNG"])
    lar_o = load_rewiring_agg(files["LAR_OLD"])

    # Interaction: |Δ_old - Δ_yng|
    print("Computing ISS-T interaction (age-dependent rewiring)...")
    iss_io = merge_on_gene(iss_o, iss_y, suffixes=("_OLD", "_YNG"))
    iss_inter = pd.DataFrame({
        "gene": iss_io["gene"],
        "delta_interaction": (iss_io["rewiring_mean_OLD"] - iss_io["rewiring_mean_YNG"]).abs(),
        "delta_old": iss_io["rewiring_mean_OLD"],
        "delta_yng": iss_io["rewiring_mean_YNG"],
    })
    iss_inter = iss_inter.sort_values("delta_interaction", ascending=False).reset_index(drop=True)
    iss_inter.to_csv(out_dir / "ISS_T_interaction.tsv", sep="\t", index=False)
    print(f"  Wrote ISS_T_interaction.tsv ({len(iss_inter)} genes)")

    print("Computing LAR interaction (age-dependent rewiring)...")
    lar_io = merge_on_gene(lar_o, lar_y, suffixes=("_OLD", "_YNG"))
    lar_inter = pd.DataFrame({
        "gene": lar_io["gene"],
        "delta_interaction": (lar_io["rewiring_mean_OLD"] - lar_io["rewiring_mean_YNG"]).abs(),
        "delta_old": lar_io["rewiring_mean_OLD"],
        "delta_yng": lar_io["rewiring_mean_YNG"],
    })
    lar_inter = lar_inter.sort_values("delta_interaction", ascending=False).reset_index(drop=True)
    lar_inter.to_csv(out_dir / "LAR_interaction.tsv", sep="\t", index=False)
    print(f"  Wrote LAR_interaction.tsv ({len(lar_inter)} genes)")

    # Persistence/recovery: (ISS - LAR) within age
    print("Computing ISS-LAR persistence for Young...")
    yl = merge_on_gene(iss_y, lar_y, suffixes=("_ISS", "_LAR"))
    persist_y = pd.DataFrame({
        "gene": yl["gene"],
        "ISS_minus_LAR_delta": yl["rewiring_mean_ISS"] - yl["rewiring_mean_LAR"],
        "ISS_delta": yl["rewiring_mean_ISS"],
        "LAR_delta": yl["rewiring_mean_LAR"],
    }).sort_values("ISS_minus_LAR_delta", ascending=False).reset_index(drop=True)
    persist_y.to_csv(out_dir / "ISS_minus_LAR_YNG_persistence.tsv", sep="\t", index=False)
    print(f"  Wrote ISS_minus_LAR_YNG_persistence.tsv ({len(persist_y)} genes)")

    print("Computing ISS-LAR persistence for Old...")
    ol = merge_on_gene(iss_o, lar_o, suffixes=("_ISS", "_LAR"))
    persist_o = pd.DataFrame({
        "gene": ol["gene"],
        "ISS_minus_LAR_delta": ol["rewiring_mean_ISS"] - ol["rewiring_mean_LAR"],
        "ISS_delta": ol["rewiring_mean_ISS"],
        "LAR_delta": ol["rewiring_mean_LAR"],
    }).sort_values("ISS_minus_LAR_delta", ascending=False).reset_index(drop=True)
    persist_o.to_csv(out_dir / "ISS_minus_LAR_OLD_persistence.tsv", sep="\t", index=False)
    print(f"  Wrote ISS_minus_LAR_OLD_persistence.tsv ({len(persist_o)} genes)")

    print(f"\n[OK] All outputs written to: {out_dir}")


if __name__ == "__main__":
    main()
