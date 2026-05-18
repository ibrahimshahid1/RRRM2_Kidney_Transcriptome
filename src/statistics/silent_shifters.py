# src/statistics/silent_shifters.py
"""
Phase 5: DE-aware silent shifter generation.

Definitions:
  * candidate rewired genes: top rewiring quantile, regardless of DE.
  * DE-supported genes: high rewiring with differential-expression support.
  * strict silent shifters: high rewiring and bounded/small mean-expression
    change (shrunken |log2FC| < threshold, most of the 95% CI lies inside a
    small-effect interval, and DE FDR is not significant).
  * supported strict subset: strict silent shifters with Phase 6 support.

Missing DE is an error by default. Older rewiring-only behavior can be requested
only with --allow_missing_de and is labelled exploratory.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT


def load_rewiring(path: Path) -> pd.DataFrame:
    """Load rewiring aggregate table, normalizing column names."""
    df = pd.read_csv(path, sep="\t")
    if "gene" not in df.columns:
        raise ValueError(f"{path} missing 'gene' column")
    if "rewiring_mean" not in df.columns:
        for c in ["delta_mean", "mean", "rewiring", "rewiring_abs_obs", "edge_sum_node_rewiring_obs"]:
            if c in df.columns:
                df = df.rename(columns={c: "rewiring_mean"})
                break
    if "rewiring_mean" not in df.columns:
        raise ValueError(f"{path} missing a rewiring mean/statistic column")
    if "rewiring_std" not in df.columns:
        for c in ["delta_std", "std", "sd"]:
            if c in df.columns:
                df = df.rename(columns={c: "rewiring_std"})
                break
    return df


def contrast_from_rewiring_path(path: Path) -> str:
    return path.stem.replace("_rewiring_agg", "").replace("_node_rewiring", "")


def de_path_for_contrast(contrast: str, de_dir: Path) -> Path:
    if contrast.endswith("_minus_GC"):
        stem = contrast.replace("_minus_GC", "_vs_GC")
    elif contrast.endswith("_minus_GND"):
        stem = contrast.replace("_minus_GND", "_vs_GC")
    else:
        stem = contrast.replace("_minus_", "_vs_")
    return de_dir / f"{stem}_gene_DE.tsv"


def load_de_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Gene-level DE table not found: {path}")
    de = pd.read_csv(path, sep=None, engine="python")
    de = de.rename(columns={
        "adj.P.Val": "FDR",
        "padj": "FDR",
        "qvalue": "FDR",
        "log2FC_shrunk": "log2FC_shrunken",
        "log2FoldChange": "log2FC",
        "logFC": "log2FC",
        "lfcSE_shrunken": "lfcSE",
        "ci_low": "lfc_ci_low",
        "ci_high": "lfc_ci_high",
    })
    if "log2FC_shrunken" not in de.columns and "log2FC" in de.columns:
        de["log2FC_shrunken"] = de["log2FC"]
    required = {"gene", "log2FC_shrunken", "FDR", "lfc_ci_low", "lfc_ci_high"}
    missing = required - set(de.columns)
    if missing:
        raise ValueError(
            f"{path} missing required DE columns for bounded silent shifters: {sorted(missing)}"
        )
    cols = ["gene", "log2FC_shrunken", "FDR", "lfc_ci_low", "lfc_ci_high"]
    if "lfcSE" in de.columns:
        cols.append("lfcSE")
    out = de[cols].drop_duplicates("gene").copy()
    out["log2FC_shrunken"] = pd.to_numeric(out["log2FC_shrunken"], errors="coerce")
    out["FDR"] = pd.to_numeric(out["FDR"], errors="coerce")
    out["lfc_ci_low"] = pd.to_numeric(out["lfc_ci_low"], errors="coerce")
    out["lfc_ci_high"] = pd.to_numeric(out["lfc_ci_high"], errors="coerce")
    if "lfcSE" in out.columns:
        out["lfcSE"] = pd.to_numeric(out["lfcSE"], errors="coerce")
    out["log2FC"] = out["log2FC_shrunken"]
    return out


def ci_inside_fraction(
    low: pd.Series | np.ndarray,
    high: pd.Series | np.ndarray,
    interval_low: float,
    interval_high: float,
) -> np.ndarray:
    """Fraction of each CI width that lies inside the small-effect interval."""
    low = np.asarray(low, dtype=float)
    high = np.asarray(high, dtype=float)
    lo = np.minimum(low, high)
    hi = np.maximum(low, high)
    width = hi - lo
    overlap = np.maximum(0.0, np.minimum(hi, interval_high) - np.maximum(lo, interval_low))
    out = np.zeros_like(width, dtype=float)
    nonzero = width > 0
    out[nonzero] = overlap[nonzero] / width[nonzero]
    degenerate = ~nonzero & np.isfinite(lo) & np.isfinite(hi)
    out[degenerate] = ((lo[degenerate] >= interval_low) & (hi[degenerate] <= interval_high)).astype(float)
    out[~np.isfinite(out)] = 0.0
    return out


def attach_support(rw: pd.DataFrame, support_path: Path | None, support_q: float) -> pd.DataFrame:
    rw = rw.copy()
    if support_path is None or not support_path.exists():
        rw["support_q_value"] = np.nan
        rw["supported"] = np.nan
        rw["support_source"] = ""
        return rw

    support = pd.read_csv(support_path, sep="\t")
    q_col = None
    for candidate in [
        "q_BH_empirical_signed",
        "q_BH_edge_sum",
        "q_BH_candidate",
        "q_BB_two_stage",
        "q_BH",
    ]:
        if candidate in support.columns:
            q_col = candidate
            break
    if q_col is None or "gene" not in support.columns:
        raise ValueError(f"Support table {support_path} must contain gene and a recognized q-value column")

    support = support[["gene", q_col]].drop_duplicates("gene").rename(columns={q_col: "support_q_value"})
    rw = rw.merge(support, on="gene", how="left")
    rw["supported"] = rw["support_q_value"].fillna(1.0) < support_q
    rw["support_source"] = support_path.name
    return rw


def hypergeom_upper_tail(population: int, successes: int, draws: int, observed: int) -> float:
    """P[X >= observed] for Hypergeometric(population, successes, draws)."""
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(observed - 1, population, successes, draws))
    except Exception:
        return float("nan")


def build_silent_shifter_tables(
    rw: pd.DataFrame,
    de_fdr_alpha: float = 0.05,
    lfc_max: float = 0.3,
    ci_inside_min: float = 0.80,
    top_quantile: float = 0.9,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Return annotated all genes, high-rewiring candidates, strict silent table, summary."""
    rw = rw.copy()
    thr = rw["rewiring_mean"].quantile(top_quantile)
    rw["high_rewiring"] = rw["rewiring_mean"] >= thr
    rw["lfc_ci_inside_fraction"] = ci_inside_fraction(
        rw["lfc_ci_low"], rw["lfc_ci_high"], -lfc_max, lfc_max
    )
    rw["bounded_lfc"] = (
        (rw["log2FC_shrunken"].abs() < lfc_max) &
        (rw["lfc_ci_inside_fraction"] >= ci_inside_min)
    )
    rw["not_de_significant"] = rw["FDR"] >= de_fdr_alpha
    rw["de_null"] = rw["bounded_lfc"] & rw["not_de_significant"]
    rw["de_supported"] = (rw["FDR"] < de_fdr_alpha) | (~rw["bounded_lfc"])
    rw["strict_silent_shifter"] = rw["high_rewiring"] & rw["de_null"]

    candidates = rw[rw["high_rewiring"]].sort_values("rewiring_mean", ascending=False).reset_index(drop=True)
    strict = rw[rw["strict_silent_shifter"]].sort_values("rewiring_mean", ascending=False).reset_index(drop=True)
    de_supported = rw[rw["high_rewiring"] & rw["de_supported"]].sort_values(
        "rewiring_mean", ascending=False
    ).reset_index(drop=True)

    total = len(rw)
    high = int(rw["high_rewiring"].sum())
    low_de = int(rw["de_null"].sum())
    observed = len(strict)
    expected = high * low_de / total if total else np.nan
    calibration_p = hypergeom_upper_tail(total, low_de, high, observed) if total else np.nan

    supported_strict = np.nan
    if "supported" in strict.columns and strict["supported"].notna().any():
        supported_strict = int((strict["supported"] == True).sum())

    summary = {
        "n_total_genes": total,
        "top_quantile": top_quantile,
        "rewiring_threshold": float(thr),
        "small_effect_interval_low": -lfc_max,
        "small_effect_interval_high": lfc_max,
        "ci_inside_min": ci_inside_min,
        "de_fdr_alpha": de_fdr_alpha,
        "n_candidate_rewired": high,
        "n_bounded_small_expression_change": low_de,
        "n_de_supported_high_rewiring": len(de_supported),
        "n_strict_silent_shifters": observed,
        "n_supported_strict_silent_shifters": supported_strict,
        "expected_strict_under_independence": float(expected),
        "observed_over_expected": float(observed / expected) if expected and expected > 0 else np.nan,
        "hypergeom_p_overrepresentation": calibration_p,
    }
    return rw, candidates, strict, summary


def main() -> None:
    ap = argparse.ArgumentParser(description="Build DE-aware strict silent shifters")
    ap.add_argument("--rewiring", required=True, help="*_rewiring_agg.tsv from Phase 3")
    ap.add_argument("--gene_de", default="", help="Contrast-matched gene-level DE TSV")
    ap.add_argument("--gene_de_dir", default=str(REPO_ROOT / "data/processed/gene_level_DE"))
    ap.add_argument("--perm", default="", help="Optional Phase 6 permutation support table")
    ap.add_argument("--regression_results", default="", help="Optional Phase 6 regression support table")
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/phase5_silent_shifters_strict"))
    ap.add_argument("--top_quantile", type=float, default=0.9)
    ap.add_argument("--de_fdr_alpha", type=float, default=0.05)
    ap.add_argument("--fdr_min", type=float, default=None,
                    help="Deprecated alias for --de_fdr_alpha")
    ap.add_argument("--lfc_max", type=float, default=0.3)
    ap.add_argument("--ci_inside_min", type=float, default=0.80)
    ap.add_argument("--support_q", type=float, default=0.1)
    ap.add_argument("--allow_missing_de", action="store_true",
                    help="Exploratory rewiring-only mode; strict silent shifters are not emitted.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rw_path = Path(args.rewiring)
    contrast = contrast_from_rewiring_path(rw_path)

    rw = load_rewiring(rw_path)
    print(f"Loaded {len(rw)} genes from {rw_path}")

    if args.gene_de:
        de_path = Path(args.gene_de)
    else:
        de_path = de_path_for_contrast(contrast, Path(args.gene_de_dir))

    if de_path.exists():
        de = load_de_table(de_path)
        rw = rw.merge(de, on="gene", how="left")
        missing_de = int(rw["FDR"].isna().sum())
        if missing_de:
            print(f"[WARN] {missing_de} rewiring genes lack DE rows in {de_path.name}; they cannot be DE-null.")
        print(f"Loaded DE for {de['gene'].nunique()} genes from {de_path}")
    elif args.allow_missing_de:
        rw["log2FC_shrunken"] = np.nan
        rw["lfc_ci_low"] = np.nan
        rw["lfc_ci_high"] = np.nan
        rw["FDR"] = np.nan
        rw["de_null"] = False
        print("[WARN] Missing DE allowed: output is exploratory and no strict silent shifters will be called.")
    else:
        raise FileNotFoundError(
            f"DE-aware silent shifters require contrast-matched DE. Expected {de_path}. "
            f"Use --gene_de to supply it, or --allow_missing_de only for exploratory rewiring-only output."
        )

    support_path = None
    if args.regression_results:
        support_path = Path(args.regression_results)
    elif args.perm:
        support_path = Path(args.perm)
    rw = attach_support(rw, support_path, args.support_q)

    de_fdr_alpha = args.de_fdr_alpha if args.fdr_min is None else args.fdr_min
    annotated, candidates, strict, summary = build_silent_shifter_tables(
        rw,
        de_fdr_alpha=de_fdr_alpha,
        lfc_max=args.lfc_max,
        ci_inside_min=args.ci_inside_min,
        top_quantile=args.top_quantile,
    )

    if args.allow_missing_de and not de_path.exists():
        strict = strict.iloc[0:0].copy()
        summary["n_strict_silent_shifters"] = 0
        summary["n_supported_strict_silent_shifters"] = np.nan

    prefix = contrast
    annotated.to_csv(outdir / f"{prefix}_all_genes_DE_rewiring.tsv", sep="\t", index=False)
    candidates.to_csv(outdir / f"{prefix}_candidate_rewired_genes.tsv", sep="\t", index=False)
    candidates[candidates["de_supported"]].to_csv(
        outdir / f"{prefix}_DE_supported_rewired_genes.tsv", sep="\t", index=False
    )
    strict.to_csv(outdir / f"{prefix}_silent_shifters.tsv", sep="\t", index=False)
    if "supported" in strict.columns:
        strict[strict["supported"] == True].to_csv(
            outdir / f"{prefix}_supported_strict_silent_shifters.tsv", sep="\t", index=False
        )
    pd.DataFrame([summary]).to_csv(outdir / f"{prefix}_silent_shifter_calibration.tsv", sep="\t", index=False)

    print("\nSilent shifter summary:")
    for key, value in summary.items():
        print(f"  {key}: {value}")
    print(f"\n[OK] Wrote DE-aware silent-shifter outputs to: {outdir}")


if __name__ == "__main__":
    main()
