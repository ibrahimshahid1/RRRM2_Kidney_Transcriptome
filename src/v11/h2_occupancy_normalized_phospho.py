#!/usr/bin/env python3
"""Parent-protein-normalized OSD-462 phosphosite robustness analysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_ROOT = REPO_ROOT / "data/results/run_20260526_v11_dct1_phospho_mediation"
DCT_PRIOR = RUN_ROOT / "dct_prior/osd462_phosphosite_dct1_prior.tsv"

TARGET_GENES = {
    "Slc12a3",
    "Stk39",
    "Oxsr1",
    "Wnk1",
    "Wnk4",
    "Klhl3",
    "Cul3",
    "Kcnj10",
    "Kcnj16",
    "Clcnkb",
    "Bsnd",
    "Ppp1ca",
    "Ppp1r1a",
    "Calb1",
}


def bh(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.ones_like(p)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(vals)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[ok] = restored
    return out


def odds_ci(a: float, b: float, c: float, d: float) -> tuple[float, float, float]:
    aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    log_or = np.log((aa * dd) / (bb * cc))
    se = np.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return float(np.exp(log_or)), float(np.exp(log_or - 1.96 * se)), float(np.exp(log_or + 1.96 * se))


def fisher_rows(df: pd.DataFrame, suppressed_col: str, analysis: str) -> list[dict]:
    rows = []
    for flag in ["dct1_top_decile", "dct1_top_quartile"]:
        table_df = df[df[flag].notna()].copy()
        sup = table_df[suppressed_col].astype(bool)
        in_flag = table_df[flag].astype(bool)
        a = int((sup & in_flag).sum())
        b = int((sup & ~in_flag).sum())
        c = int((~sup & in_flag).sum())
        d = int((~sup & ~in_flag).sum())
        odds, ci_low, ci_high = odds_ci(a, b, c, d)
        _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "analysis": analysis,
                "flag": flag,
                "n_background": int(len(table_df)),
                "n_suppressed": int(sup.sum()),
                "table_suppressed_in_flag": a,
                "table_suppressed_not_flag": b,
                "table_background_in_flag": c,
                "table_background_not_flag": d,
                "odds_ratio": odds,
                "ci_low": ci_low,
                "ci_high": ci_high,
                "p_value": float(p),
            }
        )
    return rows


def single_site_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    sort_cols = ["gene_symbol", "phospho_p_value", "occupancy_effect", "site_id"]
    return (
        df.sort_values(sort_cols, ascending=[True, True, True, True])
        .groupby("gene_symbol", as_index=False)
        .head(1)
        .copy()
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=RUN_ROOT)
    args = ap.parse_args()

    out_dir = args.run_root / "h2_occupancy"
    out_dir.mkdir(parents=True, exist_ok=True)
    prior_path = args.run_root / "dct_prior/osd462_phosphosite_dct1_prior.tsv"
    df = pd.read_csv(prior_path, sep="\t")
    df["parent_protein_effect"] = pd.to_numeric(df["flight_effect"], errors="coerce")
    df["occupancy_effect"] = pd.to_numeric(df["phospho_effect"], errors="coerce") - df["parent_protein_effect"]
    df["has_parent_protein_effect"] = df["parent_protein_effect"].notna()
    df["is_occupancy_suppressed_p05"] = (
        df["has_parent_protein_effect"]
        & (df["occupancy_effect"] < 0)
        & (pd.to_numeric(df["phospho_p_value"], errors="coerce") < 0.05)
    )
    df["is_occupancy_suppressed_q10"] = (
        df["has_parent_protein_effect"]
        & (df["occupancy_effect"] < 0)
        & (pd.to_numeric(df["phospho_q_value_all_sites"], errors="coerce") < 0.10)
    )
    df.to_csv(out_dir / "h2_occupancy_site_effects.tsv", sep="\t", index=False)

    rows = []
    base = df[df["has_parent_protein_effect"]].copy()
    rows.extend(fisher_rows(base, "is_occupancy_suppressed_p05", "occupancy_p05"))
    rows.extend(fisher_rows(base, "is_occupancy_suppressed_q10", "occupancy_q10"))
    rows.extend(fisher_rows(base[~base["is_anchor_gene"].astype(bool)].copy(), "is_occupancy_suppressed_p05", "occupancy_exclude_anchor_genes"))
    rows.extend(fisher_rows(base[~base["is_ncc_site"].astype(bool)].copy(), "is_occupancy_suppressed_p05", "occupancy_exclude_ncc_sites"))
    single = single_site_per_gene(base)
    rows.extend(fisher_rows(single, "is_occupancy_suppressed_p05", "occupancy_single_site_p05"))
    enrich = pd.DataFrame(rows)
    enrich["q_value"] = bh(enrich["p_value"])
    enrich.to_csv(out_dir / "h2_occupancy_dct1_enrichment.tsv", sep="\t", index=False)

    target = base[base["gene_symbol"].isin(TARGET_GENES) | base["is_anchor_gene"].astype(bool) | base["is_ncc_site"].astype(bool)].copy()
    target = target.sort_values(["gene_symbol", "site_position_int", "site_position_str"])
    target.to_csv(out_dir / "h2_occupancy_target_sites.tsv", sep="\t", index=False)

    primary = enrich[
        enrich["analysis"].eq("occupancy_p05") & enrich["flag"].eq("dct1_top_decile")
    ]
    primary_record = primary.iloc[0].to_dict() if not primary.empty else {}
    verdict = {
        "status": "complete",
        "primary_record": primary_record,
        "supports_parent_normalized_dct1_enrichment": bool(
            primary_record and primary_record.get("odds_ratio", 0) > 1 and primary_record.get("q_value", 1) < 0.05
        ),
        "boundary": (
            "occupancy_effect = phosphosite flight effect - parent protein flight effect. "
            "The p-value threshold still comes from the phosphosite model; this is a parent-protein-normalized "
            "enrichment robustness test, not a direct phospho-stoichiometry measurement."
        ),
    }
    (out_dir / "h2_occupancy_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(f"[h2-occupancy] complete: {out_dir}")


if __name__ == "__main__":
    main()
