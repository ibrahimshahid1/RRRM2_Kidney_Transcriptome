#!/usr/bin/env python3
"""TMT channel-centering QC for the v11 phosphosite enrichment analysis.

This module answers two audit questions around the OSD-462 TMT phosphosite
layer:

1. Does the raw phosphosite channel pattern show lower flight-channel medians?
2. Does the DCT-subtype-prior enrichment persist when phosphosite effects are
   recomputed without the within-plex channel-centering step?

The primary analysis uses channel-centered phosphosite effects written by
``scripts/osd462/02_phospho_axis.py``.  The uncentered recomputation here is a
sensitivity check, not a replacement estimator.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.multiomics.osd462_anchor import compute_site_flight_effect_lm, parse_tmt_sheet


RUN_ROOT = Path("data/results/run_20260526_v11_dct1_phospho_mediation")
OSD462 = Path("data/results/run_20260522_osd462_anchor/osd462_anchor")
PHOSPHO_XLSX = Path("data/external/osdr/OSD-462/GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx")
FISHER_ALTERNATIVE = "greater"

SHEET_CONFIG = [
    ("siteQuant_360", "gene_symbol", "Site Position", "Protein Id", "single"),
    ("siteQuant_360_compositeSite", "geneSymbol", "sitePosStr", "proteinID", "composite"),
]

PRIOR_COLS = [
    "gene_symbol",
    "mean_expr_dct1",
    "pct_detected_dct1",
    "mean_expr_dct2",
    "pct_detected_dct2",
    "log2_mean_ratio_dct1_vs_dct2",
    "dct1_enrichment_score",
    "p_value",
    "fdr",
    "dct_expression_class",
    "dct1_top_quartile",
    "dct2_bottom_quartile",
    "dct1_top_decile",
    "dct2_bottom_decile",
    "dct1_core_fdr",
    "dct2_core_fdr",
    "dct1_score_percentile",
    "dct2_leaning_percentile",
]


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


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


def load_dct_prior(run_root: Path) -> pd.DataFrame:
    de = read_tsv(run_root / "dct_prior" / "gse228367_dct1_vs_dct2_de.tsv")
    expressed = (
        (de["pct_detected_dct1"] >= 0.05)
        | (de["pct_detected_dct2"] >= 0.05)
        | (de["mean_expr_dct1"] >= 0.05)
        | (de["mean_expr_dct2"] >= 0.05)
    )
    score = pd.to_numeric(de["dct1_enrichment_score"], errors="coerce")
    q75 = score.loc[expressed].quantile(0.75)
    q25 = score.loc[expressed].quantile(0.25)
    q90 = score.loc[expressed].quantile(0.90)
    q10 = score.loc[expressed].quantile(0.10)
    de = de.copy()
    de["dct1_top_quartile"] = expressed & (score >= q75)
    de["dct2_bottom_quartile"] = expressed & (score <= q25)
    de["dct1_top_decile"] = expressed & (score >= q90)
    de["dct2_bottom_decile"] = expressed & (score <= q10)
    de["dct1_core_fdr"] = de["dct_expression_class"].eq("DCT1_core")
    de["dct2_core_fdr"] = de["dct_expression_class"].eq("DCT2_core")
    de["dct1_score_percentile"] = np.nan
    de.loc[expressed, "dct1_score_percentile"] = score.loc[expressed].rank(pct=True)
    de["dct2_leaning_percentile"] = 1 - de["dct1_score_percentile"]
    return de[PRIOR_COLS]


def load_phosphosite_effects(channel_center: bool) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for sheet, gene_col, site_col, protein_col, kind in SHEET_CONFIG:
        tab = parse_tmt_sheet(
            PHOSPHO_XLSX,
            sheet,
            gene_col=gene_col,
            peptide_cols={"Samp1-5": "Samp1-5~num_quant", "Samp6-10": "Samp6-10~num_quant"},
            id_col=protein_col,
            extra_meta_cols=[site_col],
        )
        eff = compute_site_flight_effect_lm(tab, min_per_condition=3, channel_center=channel_center)
        eff = eff.rename(columns={gene_col: "gene_symbol", site_col: "site_position"})
        eff["source_sheet"] = sheet
        eff["source_row"] = np.arange(len(eff), dtype=int)
        eff["site_kind"] = kind
        frames.append(eff)
    sites = pd.concat(frames, ignore_index=True)
    keep = [
        "source_sheet",
        "source_row",
        "gene_symbol",
        "site_position",
        "site_kind",
        "phospho_effect",
        "phospho_se",
        "phospho_p_value",
        "n_fl",
        "n_gc",
    ]
    return sites[np.isfinite(sites["phospho_effect"])][keep].copy()


def annotate_phosphosites(effects: pd.DataFrame, prior: pd.DataFrame) -> pd.DataFrame:
    phospho = effects.copy()
    phospho["phospho_q_value_all_sites"] = bh(phospho["phospho_p_value"])
    phospho["site_position_str"] = phospho["site_position"].astype(str)
    phospho["is_single_site"] = phospho["site_position_str"].str.fullmatch(r"\d+")
    phospho["site_id"] = (
        phospho["source_sheet"].astype(str)
        + ":"
        + phospho["source_row"].astype(str)
        + ":"
        + phospho["gene_symbol"].astype(str)
        + ":"
        + phospho["site_position_str"]
    )
    phospho["is_suppressed_p05"] = (phospho["phospho_effect"] < 0) & (phospho["phospho_p_value"] < 0.05)
    phospho["is_suppressed_q10"] = (phospho["phospho_effect"] < 0) & (phospho["phospho_q_value_all_sites"] < 0.10)
    effect_q25 = pd.to_numeric(phospho["phospho_effect"], errors="coerce").quantile(0.25)
    phospho["is_effect_bottom_quartile"] = pd.to_numeric(phospho["phospho_effect"], errors="coerce") <= effect_q25
    return phospho.merge(prior, on="gene_symbol", how="left", indicator="dct_prior_merge")


def fisher_table(df: pd.DataFrame, flag: str, suppressed_col: str = "is_suppressed_p05") -> tuple[float, float, float, np.ndarray]:
    sub = df[df["dct1_enrichment_score"].notna()].copy()
    row_flag = sub[flag].fillna(False).astype(bool)
    row_suppressed = sub[suppressed_col].fillna(False).astype(bool)
    arr = np.array(
        [
            [int((row_suppressed & row_flag).sum()), int((row_suppressed & ~row_flag).sum())],
            [int((~row_suppressed & row_flag).sum()), int((~row_suppressed & ~row_flag).sum())],
        ]
    )
    odds, p_value = stats.fisher_exact(arr, alternative=FISHER_ALTERNATIVE)
    frac_suppressed = arr[0, 0] / max(arr[0].sum(), 1)
    frac_background = arr[1, 0] / max(arr[1].sum(), 1)
    fold = frac_suppressed / frac_background if frac_background > 0 else np.inf
    return float(odds), float(p_value), float(fold), arr


def one_representative_site_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["gene_symbol", "phospho_p_value", "phospho_effect", "site_id"]
    sub = df.dropna(subset=["gene_symbol"]).copy()
    return sub.sort_values(cols, ascending=[True, True, True, True]).groupby("gene_symbol", as_index=False).head(1).copy()


def one_single_position_representative_site_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    return one_representative_site_per_gene(df[df["is_single_site"].fillna(False).astype(bool)].copy())


def enrichment_sensitivity(annotated_by_setting: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    analyses = [
        ("primary_p05", "phosphosite_row", lambda d: d),
        ("single_position_one_site_per_parent_gene", "single_position_parent_gene_representative_site",
         one_single_position_representative_site_per_gene),
    ]
    for setting, df in annotated_by_setting.items():
        for analysis, unit, make_df in analyses:
            sub = make_df(df)
            for flag in ["dct1_top_decile", "dct2_bottom_decile", "dct1_top_quartile", "dct2_bottom_quartile"]:
                odds, p_value, fold, arr = fisher_table(sub, flag)
                scored = sub[sub["dct1_enrichment_score"].notna()]
                rows.append(
                    {
                        "effect_estimator": setting,
                        "analysis": analysis,
                        "test": f"fisher_{flag}",
                        "unit": unit,
                        "n_background": int(len(scored)),
                        "n_suppressed": int(scored["is_suppressed_p05"].sum()),
                        "n_parent_genes": int(scored["gene_symbol"].nunique()),
                        "odds_ratio": odds,
                        "p_value": p_value,
                        "fisher_alternative": FISHER_ALTERNATIVE,
                        "fold_enrichment": fold,
                        "table_suppressed_in_flag": int(arr[0, 0]),
                        "table_suppressed_not_flag": int(arr[0, 1]),
                        "table_background_in_flag": int(arr[1, 0]),
                        "table_background_not_flag": int(arr[1, 1]),
                    }
                )
    out = pd.DataFrame(rows)
    out["q_value"] = bh(out["p_value"])
    return out


def effect_comparison(centered: pd.DataFrame, uncentered: pd.DataFrame) -> pd.DataFrame:
    keys = ["source_sheet", "source_row"]
    merged = centered.merge(
        uncentered,
        on=keys,
        suffixes=("_centered", "_uncentered"),
        how="inner",
    )
    x = pd.to_numeric(merged["phospho_effect_centered"], errors="coerce")
    y = pd.to_numeric(merged["phospho_effect_uncentered"], errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() >= 3:
        spearman_rho, spearman_p = stats.spearmanr(x[ok], y[ok])
        pearson_r = float(np.corrcoef(x[ok], y[ok])[0, 1])
    else:
        spearman_rho, spearman_p, pearson_r = np.nan, np.nan, np.nan
    return pd.DataFrame(
        [
            {
                "n_sites_compared": int(ok.sum()),
                "effect_spearman_rho": float(spearman_rho),
                "effect_spearman_p": float(spearman_p),
                "effect_pearson_r": pearson_r,
                "median_effect_centered": float(x[ok].median()) if ok.any() else np.nan,
                "median_effect_uncentered": float(y[ok].median()) if ok.any() else np.nan,
                "median_delta_centered_minus_uncentered": float((x[ok] - y[ok]).median()) if ok.any() else np.nan,
                "median_abs_delta": float((x[ok] - y[ok]).abs().median()) if ok.any() else np.nan,
                "max_abs_delta": float((x[ok] - y[ok]).abs().max()) if ok.any() else np.nan,
                "n_suppressed_p05_centered": int(centered["is_suppressed_p05"].sum()),
                "n_suppressed_p05_uncentered": int(uncentered["is_suppressed_p05"].sum()),
            }
        ]
    )


def raw_channel_pattern_summary(run_root: Path) -> pd.DataFrame:
    qc_path = run_root / "baseline" / "osd462_tmt_channel_qc.tsv"
    if not qc_path.exists():
        return pd.DataFrame()
    qc = read_tsv(qc_path)
    sub = qc[qc["layer"].astype(str).str.startswith("phosphosite")].copy()
    grouped = (
        sub.groupby(["layer", "metric", "plex", "condition"], as_index=False)
        .agg(
            n_channels=("sample_label", "nunique"),
            median_channel_median=("median", "median"),
            median_missing_fraction=("missing_fraction", "median"),
        )
    )
    wide = grouped.pivot(index=["layer", "metric", "plex"], columns="condition", values="median_channel_median").reset_index()
    wide.columns.name = None
    for cond in ["BL", "FL", "GC"]:
        if cond not in wide.columns:
            wide[cond] = np.nan
    wide = wide.rename(columns={"BL": "median_bl", "FL": "median_fl", "GC": "median_gc"})
    wide["fl_minus_gc_channel_median"] = wide["median_fl"] - wide["median_gc"]
    wide["fl_minus_bl_channel_median"] = wide["median_fl"] - wide["median_bl"]
    wide["fl_over_gc_channel_median"] = wide["median_fl"] / wide["median_gc"]
    return wide.sort_values(["layer", "metric", "plex"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default=str(RUN_ROOT))
    args = parser.parse_args()
    run_root = Path(args.run_root)
    (run_root / "baseline").mkdir(parents=True, exist_ok=True)
    (run_root / "h2_enrichment").mkdir(parents=True, exist_ok=True)

    raw_pattern = raw_channel_pattern_summary(run_root)
    raw_pattern.to_csv(run_root / "baseline" / "osd462_tmt_raw_channel_pattern_summary.tsv", sep="\t", index=False)

    prior = load_dct_prior(run_root)
    centered = annotate_phosphosites(load_phosphosite_effects(channel_center=True), prior)
    uncentered = annotate_phosphosites(load_phosphosite_effects(channel_center=False), prior)
    annotated = {
        "channel_centered_primary": centered,
        "uncentered_sensitivity": uncentered,
    }

    sens = enrichment_sensitivity(annotated)
    sens.to_csv(run_root / "h2_enrichment" / "h2_dct_channel_centering_sensitivity.tsv", sep="\t", index=False)
    comp = effect_comparison(centered, uncentered)
    comp.to_csv(run_root / "h2_enrichment" / "h2_dct_channel_centering_effect_comparison.tsv", sep="\t", index=False)

    verdict = {
        "raw_channel_pattern_table": str(run_root / "baseline" / "osd462_tmt_raw_channel_pattern_summary.tsv"),
        "enrichment_sensitivity_table": str(run_root / "h2_enrichment" / "h2_dct_channel_centering_sensitivity.tsv"),
        "effect_comparison_table": str(run_root / "h2_enrichment" / "h2_dct_channel_centering_effect_comparison.tsv"),
        "interpretation": (
            "Raw phosphosite channel medians can differ by condition, but the primary phosphosite effects are "
            "estimated after within-plex channel median centering; uncentered recomputation is reported as a "
            "sensitivity check."
        ),
    }
    (run_root / "h2_enrichment" / "h2_dct_channel_centering_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(f"channel-centering QC complete: {run_root}")


if __name__ == "__main__":
    main()
