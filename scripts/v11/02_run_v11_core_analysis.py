#!/usr/bin/env python3
"""Execute the core v11 DCT1/phosphoproteome/mediation analyses."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats


RUN_ROOT = Path("data/results/run_20260526_v11_dct1_phospho_mediation")
OSD462 = Path("data/results/run_20260522_osd462_anchor/osd462_anchor")
PHENO = Path("data/results/run_20260522_phenotype_anchor")
CELLTYPE = Path("data/results/run_20260522_celltype_decomposition")
REGULATOR = Path("data/results/run_20260522_regulator_activity")
DCT_PRIOR_DIR = RUN_ROOT / "dct_prior"
PXD_DIR = Path("data/external/phosphoproteomics/PXD001729")


ANCHOR_GENES = {"Slc12a3", "Stk39", "Oxsr1", "Wnk1", "Wnk4"}
TRANSPORT_TARGETS = {
    "Slc12a3",
    "Stk39",
    "Oxsr1",
    "Wnk1",
    "Wnk4",
    "Klhl3",
    "Cul3",
    "Nedd4l",
    "Sgk1",
    "Kcnj10",
    "Kcnj16",
    "Trpm6",
    "Pvalb",
    "Calb1",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(pvals: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvals), dtype=float)
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


def safe_float(s):
    return pd.to_numeric(s, errors="coerce")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() == 0:
        return float("nan")
    a = a[keep]
    b = b[keep]
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return float("nan")
    return float(np.dot(a, b) / denom)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def ensure_dirs(root: Path):
    for name in [
        "baseline",
        "external_qc",
        "dct_prior",
        "h2_enrichment",
        "h2_pxd",
        "h2_klhl3",
        "h3_mediation",
        "manifests",
    ]:
        (root / name).mkdir(parents=True, exist_ok=True)


def baseline_lock(root: Path):
    summary = json.loads((OSD462 / "results_summary.json").read_text())
    pheno_summary = read_tsv(PHENO / "phenotype_anchor_summary.tsv")
    cell_corr = read_tsv(CELLTYPE / "osd462_compartment_vs_ncc_phospho.tsv")
    ksea = read_tsv(REGULATOR / "osd462_kinase_activity_summary.tsv")

    rows = [
        {
            "component": "OSD462_RNA_recurrence",
            "metric": "pathway_cosine",
            "value": summary["layer4_rna_gate"]["point_cosine"],
            "status": "PASS" if summary["layer4_rna_gate"]["recurrence_pass"] else "FAIL",
            "interpretation": "OSD-462 RNA recurs RRRM-2 ISS-T matrix-high/DCT-low direction",
        },
        {
            "component": "OSD462_protein",
            "metric": "any_target_set_concordant",
            "value": summary["layer1_protein"]["any_set_concordant_in_predicted_direction"],
            "status": "NULL",
            "interpretation": "targeted protein abundance concordance remains null",
        },
        {
            "component": "OSD462_phospho",
            "metric": "ncc_regulatory_mean_phospho",
            "value": summary["layer2_phospho"]["ncc_regulatory_mean_phospho"],
            "status": "SUPPORTED",
            "interpretation": "NCC regulatory phosphosites suppressed with total NCC protein flat",
        },
    ]
    for _, row in ksea.iterrows():
        rows.append(
            {
                "component": "KSEA",
                "metric": row["kinase"],
                "value": row["z_score"],
                "status": row["direction"],
                "interpretation": "positive-control kinase activity anchor",
            }
        )
    for _, row in pheno_summary.iterrows():
        rows.append(
            {
                "component": "phenotype_anchor",
                "metric": row["comparison"],
                "value": row["spearman_condition_adjusted"],
                "status": row["interpretation"],
                "interpretation": "animal-matched RNA-phospho link remains suggestive/underpowered except non-regulatory wrinkle",
            }
        )
    for _, row in cell_corr.iterrows():
        rows.append(
            {
                "component": "celltype_vs_ncc_phospho",
                "metric": row["panel"],
                "value": row["spearman_vs_ncc_phospho"],
                "status": "correlation",
                "interpretation": "compartment score correlation with NCC regulatory phosphorylation",
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(root / "baseline" / "v11_baseline_lock_summary.tsv", sep="\t", index=False)

    inputs = [
        Path("docs/v11_execution_research_plan.md"),
        OSD462 / "results_summary.json",
        PHENO / "phenotype_anchor_summary.tsv",
        PHENO / "phenotype_anchor_per_animal.tsv",
        CELLTYPE / "osd462_compartment_scores_per_sample.tsv",
        CELLTYPE / "osd462_compartment_vs_ncc_phospho.tsv",
        REGULATOR / "osd462_kinase_activity_summary.tsv",
    ]
    manifest = {
        "analysis": "v11 baseline lock",
        "run_root": str(root),
        "inputs": [{"path": str(p), "sha256": sha256(p)} for p in inputs if p.exists()],
    }
    (root / "baseline" / "v11_baseline_input_manifest.json").write_text(json.dumps(manifest, indent=2))


def parse_site(site) -> list[tuple[int, str, str]]:
    if pd.isna(site):
        return []
    text = str(site)
    found = []
    for match in re.finditer(r"(\d+)_([STY])([*?]?)", text):
        found.append((int(match.group(1)), match.group(2), match.group(3) or ""))
    return found


def load_pxd_tables(root: Path):
    xls = PXD_DIR / "41598_2015_BFsrep12829_MOESM6_ESM.xls"
    inventory_rows = []
    for sheet in pd.ExcelFile(xls, engine="xlrd").sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet, engine="xlrd", header=None)
        inventory_rows.append(
            {
                "file": xls.name,
                "sheet": sheet,
                "n_rows_raw": len(df),
                "n_cols_raw": df.shape[1],
            }
        )
    pd.DataFrame(inventory_rows).to_csv(root / "external_qc" / "pxd001729_table_inventory.tsv", sep="\t", index=False)

    total = pd.read_excel(xls, sheet_name="Total phosphosites", engine="xlrd", header=13)
    total = total.rename(columns=lambda c: str(c).strip())
    total["ddavp_effect_log2"] = np.nan
    total["p_neg_log10"] = np.nan
    total["ddavp_direction"] = "identified_total"

    changed_frames = []
    for sheet, direction in [
        ("phos_sites_increased_with_dDAVP", "increased"),
        ("phos_sites_decreased_with_dDAVP", "decreased"),
    ]:
        df = pd.read_excel(xls, sheet_name=sheet, engine="xlrd", header=2)
        df = df.rename(columns=lambda c: str(c).strip())
        df["ddavp_direction"] = direction
        changed_frames.append(df)
    changed = pd.concat(changed_frames, ignore_index=True)

    def normalize(df: pd.DataFrame) -> pd.DataFrame:
        colmap = {
            "Peptide Sequence": "peptide_sequence",
            "Phospho. Sites": "phospho_sites",
            "Phospho. sites": "phospho_sites",
            "Accession No.": "protein_accession",
            "Protein Accession": "protein_accession",
            "Protein Name": "protein_name",
            "Gene Symbol": "gene_symbol",
            "Mean (log2)": "ddavp_effect_log2",
            "p (-log10)": "p_neg_log10",
            "Site Conservation Score": "site_conservation_score",
            "Motif Conservation Score": "motif_conservation_score",
        }
        keep = [c for c in df.columns if c in colmap]
        out = df[keep].rename(columns=colmap).copy()
        for c in ["gene_symbol", "phospho_sites", "peptide_sequence", "protein_accession", "protein_name"]:
            if c in out.columns:
                out[c] = out[c].astype(str).str.strip()
        if "ddavp_effect_log2" not in out.columns:
            out["ddavp_effect_log2"] = np.nan
        if "p_neg_log10" not in out.columns:
            out["p_neg_log10"] = np.nan
        return out

    total_n = normalize(total)
    total_n["source_sheet"] = "Total phosphosites"
    total_n["ddavp_direction"] = "identified_total"
    changed_n = normalize(changed)
    changed_n["source_sheet"] = changed["ddavp_direction"].map(
        {
            "increased": "phos_sites_increased_with_dDAVP",
            "decreased": "phos_sites_decreased_with_dDAVP",
        }
    )
    changed_n["ddavp_direction"] = changed["ddavp_direction"].to_numpy()

    rows = []
    for _, row in pd.concat([total_n, changed_n], ignore_index=True).iterrows():
        for pos, residue, confidence in parse_site(row.get("phospho_sites")):
            d = row.to_dict()
            d["site_position"] = pos
            d["site_residue"] = residue
            d["site_assignment"] = confidence
            d["site_id"] = f"{d.get('gene_symbol')}:{residue}{pos}"
            rows.append(d)
    parsed = pd.DataFrame(rows)
    parsed["ddavp_effect_log2"] = safe_float(parsed["ddavp_effect_log2"])
    parsed["p_neg_log10"] = safe_float(parsed["p_neg_log10"])
    parsed["p_value"] = 10 ** (-parsed["p_neg_log10"])
    parsed.loc[parsed["p_neg_log10"].isna(), "p_value"] = np.nan
    parsed.to_csv(root / "external_qc" / "pxd001729_phosphosite_effects.tsv", sep="\t", index=False)

    target = parsed[parsed["gene_symbol"].isin(TRANSPORT_TARGETS)].copy()
    cov = (
        target.groupby("gene_symbol", dropna=False)
        .agg(
            n_sites=("site_id", "count"),
            n_changed=("ddavp_effect_log2", lambda s: int(s.notna().sum())),
            n_increased=("ddavp_direction", lambda s: int((s == "increased").sum())),
            n_decreased=("ddavp_direction", lambda s: int((s == "decreased").sum())),
        )
        .reset_index()
    )
    cov.to_csv(root / "external_qc" / "pxd001729_target_site_coverage.tsv", sep="\t", index=False)
    changed_n.to_csv(root / "h2_pxd" / "pxd001729_ddavp_direction.tsv", sep="\t", index=False)
    target.to_csv(root / "h2_pxd" / "pxd001729_dct_target_direction.tsv", sep="\t", index=False)
    return parsed


def build_dct_prior_mapping(root: Path):
    de = read_tsv(DCT_PRIOR_DIR / "gse228367_dct1_vs_dct2_de.tsv")
    expressed = (
        (de["pct_detected_dct1"] >= 0.05)
        | (de["pct_detected_dct2"] >= 0.05)
        | (de["mean_expr_dct1"] >= 0.05)
        | (de["mean_expr_dct2"] >= 0.05)
    )
    q75 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.75)
    q25 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.25)
    q90 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.90)
    de["dct1_top_quartile"] = expressed & (de["dct1_enrichment_score"] >= q75)
    de["dct2_bottom_quartile"] = expressed & (de["dct1_enrichment_score"] <= q25)
    de["dct1_top_decile"] = expressed & (de["dct1_enrichment_score"] >= q90)
    de["dct1_core_fdr"] = de["dct_expression_class"].eq("DCT1_core")
    de["dct2_core_fdr"] = de["dct_expression_class"].eq("DCT2_core")

    proteins = read_tsv(OSD462 / "protein_effects_gene_anyplex.tsv")
    phospho = read_tsv(OSD462 / "phospho_all_sites.tsv")
    phospho["phospho_q_value_all_sites"] = bh(phospho["phospho_p_value"])
    phospho["site_position_str"] = phospho["site_position"].astype(str)
    phospho["is_single_site"] = phospho["site_position_str"].str.fullmatch(r"\d+")
    phospho["site_position_int"] = pd.to_numeric(phospho["site_position_str"], errors="coerce")
    phospho["site_id"] = phospho["gene_symbol"].astype(str) + ":" + phospho["site_position_str"]
    phospho["is_anchor_gene"] = phospho["gene_symbol"].isin(ANCHOR_GENES)
    phospho["is_ncc_site"] = phospho["gene_symbol"].eq("Slc12a3")
    phospho["is_suppressed_p05"] = (phospho["phospho_effect"] < 0) & (phospho["phospho_p_value"] < 0.05)
    phospho["is_suppressed_q10"] = (phospho["phospho_effect"] < 0) & (phospho["phospho_q_value_all_sites"] < 0.10)

    prior_cols = [
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
        "dct1_core_fdr",
        "dct2_core_fdr",
    ]
    proteins_prior = proteins.merge(de[prior_cols], on="gene_symbol", how="left", indicator="dct_prior_merge")
    phospho_prior = phospho.merge(de[prior_cols], on="gene_symbol", how="left", indicator="dct_prior_merge")
    phospho_prior = phospho_prior.merge(
        proteins[["gene_symbol", "flight_effect", "n_peptides", "abundance_log2", "plex_coverage"]],
        on="gene_symbol",
        how="left",
        suffixes=("", "_protein"),
    )
    phospho_prior["abundance_bin"] = pd.qcut(phospho_prior["abundance_log2"].rank(method="first"), 4, labels=False, duplicates="drop")
    phospho_prior["peptide_bin"] = pd.qcut(phospho_prior["n_peptides"].rank(method="first"), 4, labels=False, duplicates="drop")
    phospho_prior["match_stratum"] = phospho_prior["abundance_bin"].astype("Int64").astype(str) + "_p" + phospho_prior[
        "peptide_bin"
    ].astype("Int64").astype(str)

    proteins_prior.to_csv(root / "dct_prior" / "osd462_protein_dct1_prior.tsv", sep="\t", index=False)
    phospho_prior.to_csv(root / "dct_prior" / "osd462_phosphosite_dct1_prior.tsv", sep="\t", index=False)
    de.to_csv(root / "dct_prior" / "dct1_enrichment_prior_v1.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(
        [
            {
                "table": "osd462_proteins",
                "n_rows": len(proteins_prior),
                "n_mapped": int(proteins_prior["dct1_enrichment_score"].notna().sum()),
                "fraction_mapped": proteins_prior["dct1_enrichment_score"].notna().mean(),
            },
            {
                "table": "osd462_phosphosites",
                "n_rows": len(phospho_prior),
                "n_mapped": int(phospho_prior["dct1_enrichment_score"].notna().sum()),
                "fraction_mapped": phospho_prior["dct1_enrichment_score"].notna().mean(),
            },
        ]
    )
    target_cov = phospho_prior[phospho_prior["gene_symbol"].isin(TRANSPORT_TARGETS)].groupby("gene_symbol").agg(
        n_phosphosites=("site_id", "count"),
        mapped=("dct1_enrichment_score", lambda s: bool(s.notna().any())),
        dct1_enrichment_score=("dct1_enrichment_score", "first"),
        dct_expression_class=("dct_expression_class", "first"),
        dct1_top_quartile=("dct1_top_quartile", "first"),
    )
    target_cov = target_cov.reset_index()
    coverage.to_csv(root / "dct_prior" / "osd462_dct1_prior_coverage.tsv", sep="\t", index=False)
    target_cov.to_csv(root / "dct_prior" / "osd462_dct1_prior_target_gene_coverage.tsv", sep="\t", index=False)
    return de, proteins_prior, phospho_prior


def fisher_table(df: pd.DataFrame, flag: str, suppressed_col: str):
    sub = df[df["dct1_enrichment_score"].notna()].copy()
    tab = pd.crosstab(sub[suppressed_col], sub[flag].fillna(False))
    for r in [False, True]:
        if r not in tab.index:
            tab.loc[r] = 0
    for c in [False, True]:
        if c not in tab.columns:
            tab[c] = 0
    tab = tab.sort_index().sort_index(axis=1)
    arr = np.array([[tab.loc[True, True], tab.loc[True, False]], [tab.loc[False, True], tab.loc[False, False]]])
    odds, p = stats.fisher_exact(arr)
    frac_sup = arr[0, 0] / max(arr[0].sum(), 1)
    frac_bg = arr[1, 0] / max(arr[1].sum(), 1)
    fold = frac_sup / frac_bg if frac_bg > 0 else np.inf
    return odds, p, fold, arr


def matched_null_mean(df: pd.DataFrame, suppressed_col: str, n_draws=5000, seed=20260526):
    rng = np.random.default_rng(seed)
    sub = df[df["dct1_enrichment_score"].notna()].copy()
    sup = sub[sub[suppressed_col]]
    bg = sub[~sub[suppressed_col]]
    obs = sup["dct1_enrichment_score"].mean()
    if len(sup) == 0 or len(bg) == 0:
        return obs, np.nan, np.nan, np.nan, 0
    draws = []
    strata_counts = sup["match_stratum"].fillna("missing").value_counts().to_dict()
    bg_strata = {k: g for k, g in bg.assign(match_stratum=bg["match_stratum"].fillna("missing")).groupby("match_stratum")}
    all_bg = bg
    for _ in range(n_draws):
        sampled = []
        for stratum, n in strata_counts.items():
            pool = bg_strata.get(stratum, all_bg)
            if len(pool) == 0:
                pool = all_bg
            sampled.append(pool.sample(n=n, replace=len(pool) < n, random_state=int(rng.integers(0, 2**31 - 1))))
        draw = pd.concat(sampled, ignore_index=True)
        draws.append(draw["dct1_enrichment_score"].mean())
    draws = np.asarray(draws)
    p_emp = (np.sum(draws >= obs) + 1) / (len(draws) + 1)
    return obs, float(np.median(draws)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975)), float(p_emp)


def run_h2_enrichment(root: Path, phospho_prior: pd.DataFrame):
    rows = []
    families = [
        ("primary_p05", "is_suppressed_p05", phospho_prior),
        ("strict_q10", "is_suppressed_q10", phospho_prior),
        ("exclude_anchor_genes", "is_suppressed_p05", phospho_prior[~phospho_prior["is_anchor_gene"]]),
        ("exclude_ncc_sites", "is_suppressed_p05", phospho_prior[~phospho_prior["is_ncc_site"]]),
        ("single_sites_only", "is_suppressed_p05", phospho_prior[phospho_prior["is_single_site"]]),
    ]
    for label, suppressed_col, df in families:
        sub = df[df["dct1_enrichment_score"].notna()].copy()
        sup = sub[sub[suppressed_col]]
        nonsup = sub[~sub[suppressed_col]]
        mw = stats.mannwhitneyu(
            sup["dct1_enrichment_score"],
            nonsup["dct1_enrichment_score"],
            alternative="greater",
        ) if len(sup) and len(nonsup) else (np.nan, np.nan)
        if hasattr(mw, "statistic"):
            mw_stat, mw_p = mw.statistic, mw.pvalue
        else:
            mw_stat, mw_p = mw
        obs, null_med, null_lo, null_hi, null_p = matched_null_mean(sub, suppressed_col)
        rows.append(
            {
                "analysis": label,
                "test": "mann_whitney_continuous_dct1_score",
                "n_background": len(sub),
                "n_suppressed": len(sup),
                "observed_mean_dct1_score_suppressed": obs,
                "background_mean_dct1_score": nonsup["dct1_enrichment_score"].mean(),
                "statistic": mw_stat,
                "p_value": mw_p,
                "fold_enrichment": np.nan,
                "null_median": np.nan,
                "null_lo": np.nan,
                "null_hi": np.nan,
            }
        )
        rows.append(
            {
                "analysis": label,
                "test": "matched_null_mean_dct1_score",
                "n_background": len(sub),
                "n_suppressed": len(sup),
                "observed_mean_dct1_score_suppressed": obs,
                "background_mean_dct1_score": nonsup["dct1_enrichment_score"].mean(),
                "statistic": obs,
                "p_value": null_p,
                "fold_enrichment": np.nan,
                "null_median": null_med,
                "null_lo": null_lo,
                "null_hi": null_hi,
            }
        )
        for flag in ["dct1_core_fdr", "dct1_top_quartile", "dct1_top_decile", "dct2_bottom_quartile"]:
            odds, p, fold, arr = fisher_table(sub, flag, suppressed_col)
            rows.append(
                {
                    "analysis": label,
                    "test": f"fisher_{flag}",
                    "n_background": len(sub),
                    "n_suppressed": len(sup),
                    "observed_mean_dct1_score_suppressed": obs,
                    "background_mean_dct1_score": nonsup["dct1_enrichment_score"].mean(),
                    "statistic": odds,
                    "p_value": p,
                    "fold_enrichment": fold,
                    "null_median": np.nan,
                    "null_lo": np.nan,
                    "null_hi": np.nan,
                    "table_suppressed_in_flag": int(arr[0, 0]),
                    "table_suppressed_not_flag": int(arr[0, 1]),
                    "table_background_in_flag": int(arr[1, 0]),
                    "table_background_not_flag": int(arr[1, 1]),
                }
            )
    summary = pd.DataFrame(rows)
    summary["q_value"] = bh(summary["p_value"])
    summary.to_csv(root / "h2_enrichment" / "h2_dct1_phosphosite_enrichment_summary.tsv", sep="\t", index=False)
    phospho_prior.to_csv(root / "h2_enrichment" / "h2_dct1_phosphosite_enrichment_background.tsv", sep="\t", index=False)

    primary_continuous = summary[
        (summary["analysis"] == "primary_p05")
        & (summary["test"].isin(["mann_whitney_continuous_dct1_score", "matched_null_mean_dct1_score"]))
    ]
    primary_binary = summary[
        (summary["analysis"] == "primary_p05")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    anchor_continuous = summary[
        (summary["analysis"] == "exclude_anchor_genes")
        & (summary["test"].isin(["mann_whitney_continuous_dct1_score", "matched_null_mean_dct1_score"]))
    ]
    anchor_binary = summary[
        (summary["analysis"] == "exclude_anchor_genes")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    ncc_binary = summary[
        (summary["analysis"] == "exclude_ncc_sites")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    primary_continuous_pass = bool((primary_continuous["q_value"] <= 0.10).any())
    primary_binary_pass = bool((primary_binary["q_value"] <= 0.10).any())
    survives_anchor_continuous = bool((anchor_continuous["q_value"] <= 0.10).any())
    survives_anchor_binary = bool((anchor_binary["q_value"] <= 0.10).any())
    survives_ncc_binary = bool((ncc_binary["q_value"] <= 0.10).any())
    passes = primary_continuous_pass or primary_binary_pass
    survives_anchor_exclusion = survives_anchor_continuous or survives_anchor_binary
    if primary_continuous_pass and survives_anchor_continuous:
        interpretation = "Broad continuous DCT1-prior enrichment of suppressed phosphosites is supported."
    elif primary_binary_pass and survives_anchor_binary and survives_ncc_binary:
        interpretation = (
            "DCT1-high percentile enrichment is supported and survives anchor/NCC exclusion, "
            "but the continuous DCT1-score shift is weak; claim as DCT1-prioritized subset enrichment only."
        )
    elif passes:
        interpretation = (
            "A limited DCT1-prior signal is present, but sensitivity analyses do not support a broad DCT1-specific claim."
        )
    else:
        interpretation = "DCT1 enrichment is not supported under the pre-specified tests."
    verdict = {
        "primary_h2_supported_at_fdr_0_10": passes,
        "survives_anchor_gene_exclusion_at_fdr_0_10": survives_anchor_exclusion,
        "primary_continuous_pass": primary_continuous_pass,
        "primary_binary_percentile_pass": primary_binary_pass,
        "survives_anchor_continuous": survives_anchor_continuous,
        "survives_anchor_binary_percentile": survives_anchor_binary,
        "survives_ncc_binary_percentile": survives_ncc_binary,
        "interpretation": interpretation,
        "claim_caution": "OSD-462 remains whole-kidney phosphoproteomics; DCT1 specificity is reference-prior inference only.",
    }
    (root / "h2_enrichment" / "h2_dct1_enrichment_verdict.json").write_text(json.dumps(verdict, indent=2))
    summary.to_csv(root / "h2_enrichment" / "h2_dct1_sensitivity_summary.tsv", sep="\t", index=False)
    return summary, verdict


def run_pxd_antialignment(root: Path, pxd: pd.DataFrame, phospho_prior: pd.DataFrame):
    changed = pxd[pxd["ddavp_effect_log2"].notna()].copy()
    changed = changed[changed["gene_symbol"].notna() & changed["site_position"].notna()]
    changed["site_position_int"] = changed["site_position"].astype(int)
    changed_agg = (
        changed.groupby(["gene_symbol", "site_position_int"], as_index=False)
        .agg(
            ddavp_effect_log2=("ddavp_effect_log2", "mean"),
            p_neg_log10=("p_neg_log10", "max"),
            n_pxd_rows=("site_id", "count"),
        )
    )
    osd = phospho_prior[phospho_prior["is_single_site"]].copy()
    osd["site_position_int"] = osd["site_position_int"].astype("Int64")
    shared = osd.merge(changed_agg, on=["gene_symbol", "site_position_int"], how="inner")
    shared["site_id_shared"] = shared["gene_symbol"] + ":" + shared["site_position_int"].astype(str)
    shared.to_csv(root / "h2_pxd" / "pxd001729_osd462_shared_phosphosites.tsv", sep="\t", index=False)
    target_shared = shared[shared["gene_symbol"].isin(TRANSPORT_TARGETS)].copy()

    def summarize(label, df):
        c = cosine(df["ddavp_effect_log2"], df["phospho_effect"])
        rng = np.random.default_rng(20260526)
        boots = []
        if len(df) >= 3:
            for _ in range(2000):
                idx = rng.integers(0, len(df), len(df))
                boots.append(cosine(df["ddavp_effect_log2"].to_numpy()[idx], df["phospho_effect"].to_numpy()[idx]))
        boots = np.asarray(boots, dtype=float)
        return {
            "comparison": label,
            "n_shared_sites": len(df),
            "n_shared_genes": df["gene_symbol"].nunique() if len(df) else 0,
            "cosine_ddavp_vs_spaceflight": c,
            "ci_low": np.nanquantile(boots, 0.025) if len(boots) else np.nan,
            "ci_high": np.nanquantile(boots, 0.975) if len(boots) else np.nan,
            "interpretation": (
                "anti_aligned_with_dct_activation"
                if len(df) >= 3 and c < 0 and np.nanquantile(boots, 0.975) < 0
                else "suggestive_or_not_stable_reference_signal"
            ),
        }

    rows = [
        summarize("all_shared_single_sites", shared),
        summarize("transport_target_shared_sites", target_shared),
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(root / "h2_pxd" / "h2_pxd001729_ddavp_antialignment_summary.tsv", sep="\t", index=False)

    missing_targets = sorted(TRANSPORT_TARGETS - set(target_shared["gene_symbol"]))
    verdict = {
        "n_all_shared_single_sites": int(len(shared)),
        "n_transport_target_shared_sites": int(len(target_shared)),
        "transport_targets_missing_from_shared_map": missing_targets,
        "claim_caution": "PXD001729 is mpkDCT cell-line dDAVP phosphoproteomics, not spaceflight tissue.",
        "verdict": summary.loc[0, "interpretation"] if len(summary) else "not_testable",
    }
    (root / "h2_pxd" / "h2_pxd001729_ddavp_antialignment_verdict.json").write_text(json.dumps(verdict, indent=2))
    failures = phospho_prior[phospho_prior["gene_symbol"].isin(TRANSPORT_TARGETS) & phospho_prior["is_single_site"]].copy()
    failures["present_in_pxd_changed_same_site"] = failures.set_index(["gene_symbol", "site_position_int"]).index.isin(
        changed_agg.set_index(["gene_symbol", "site_position_int"]).index
    )
    failures.to_csv(root / "h2_pxd" / "pxd001729_osd462_mapping_failures.tsv", sep="\t", index=False)
    return summary, verdict


def run_klhl3_check(root: Path, pxd: pd.DataFrame, phospho_prior: pd.DataFrame, proteins_prior: pd.DataFrame):
    genes = ["Klhl3", "Cul3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Slc12a3", "Nedd4l", "Sgk1"]
    osd_sites = phospho_prior[phospho_prior["gene_symbol"].isin(genes)].copy()
    pxd_sites = pxd[pxd["gene_symbol"].isin(genes)].copy()
    prot = proteins_prior[proteins_prior["gene_symbol"].isin(genes)].copy()
    osd_sites.to_csv(root / "h2_klhl3" / "h2_klhl3_cul3_site_coverage.tsv", sep="\t", index=False)
    combined = []
    for g in genes:
        combined.append(
            {
                "gene_symbol": g,
                "osd462_n_phosphosites": int((osd_sites["gene_symbol"] == g).sum()),
                "osd462_min_phospho_p": osd_sites.loc[osd_sites["gene_symbol"] == g, "phospho_p_value"].min(),
                "osd462_mean_phospho_effect": osd_sites.loc[osd_sites["gene_symbol"] == g, "phospho_effect"].mean(),
                "pxd001729_n_phosphosites": int((pxd_sites["gene_symbol"] == g).sum()),
                "pxd001729_n_ddavp_changed": int(pxd_sites.loc[pxd_sites["gene_symbol"] == g, "ddavp_effect_log2"].notna().sum()),
                "osd462_protein_effect": prot.loc[prot["gene_symbol"] == g, "flight_effect"].mean(),
                "osd462_rna_effect": prot.loc[prot["gene_symbol"] == g, "osd462_rna_effect"].mean()
                if "osd462_rna_effect" in prot.columns
                else np.nan,
            }
        )
    effects = pd.DataFrame(combined)
    effects.to_csv(root / "h2_klhl3" / "h2_klhl3_cul3_effects.tsv", sep="\t", index=False)
    has_klhl3_433 = bool(((osd_sites["gene_symbol"] == "Klhl3") & (osd_sites["site_position"].astype(str) == "433")).any())
    interpretation = [
        "# KLHL3/CUL3 exploratory check",
        "",
        f"KLHL3 Ser433 detected in OSD-462 targeted phosphosite table: {has_klhl3_433}.",
        "",
        "This check is exploratory because public OSD-462 does not include ubiquitinomics.",
        "A KLHL3/WNK turnover mechanism cannot be distinguished from ionic/osmotic WNK-SPAK suppression without perturbation or ubiquitin-remnant data.",
    ]
    (root / "h2_klhl3" / "h2_klhl3_cul3_interpretation.md").write_text("\n".join(interpretation) + "\n")


def sample_to_animal(sample: str) -> tuple[str, str] | None:
    m = re.search(r"RR10_KDN_WT_(FLT|GC|BSL|VIV)_([A-Z])(\d+)", sample)
    if not m:
        return None
    cond_code, _, number = m.groups()
    cond = {
        "FLT": "Space Flight",
        "GC": "Ground Control",
        "BSL": "Basal",
        "VIV": "Vivarium",
    }[cond_code]
    return cond, f"{cond}|{int(number)}"


def approximate_bayes_linear(y, X, n_draws=10000, seed=20260526):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    beta_hat = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta_hat
    df = max(n - p, 1)
    sse = float(np.sum(resid**2))
    xtx_inv = np.linalg.pinv(X.T @ X)
    sigma2 = sse / rng.chisquare(df, size=n_draws)
    draws = np.empty((n_draws, p))
    for k in range(n_draws):
        cov = sigma2[k] * xtx_inv
        draws[k] = rng.multivariate_normal(beta_hat, cov)
    return draws, np.sqrt(sigma2), beta_hat, resid


def run_mediation(root: Path):
    pheno = read_tsv(PHENO / "phenotype_anchor_per_animal.tsv")
    scores = pd.read_csv(CELLTYPE / "osd462_compartment_scores_per_sample.tsv", sep="\t", index_col=0).T
    rows = []
    for sample, vals in scores.iterrows():
        parsed = sample_to_animal(sample)
        if parsed is None:
            continue
        condition, animal = parsed
        d = {"sample": sample, "condition": condition, "animal": animal}
        d.update(vals.to_dict())
        rows.append(d)
    comp = pd.DataFrame(rows)
    comp = comp.groupby(["animal", "condition"], as_index=False).mean(numeric_only=True)
    merged = pheno.merge(comp, on=["animal", "condition"], how="inner")
    merged = merged[merged["condition"].isin(["Space Flight", "Ground Control"])].copy()
    merged["flight"] = (merged["condition"] == "Space Flight").astype(float)
    merged["matrix_endothelial_composite"] = merged[[c for c in ["endothelial", "stromal_fibroblast"] if c in merged]].mean(axis=1)
    merged.to_csv(root / "h3_mediation" / "h3_mediation_input_scores.tsv", sep="\t", index=False)

    mediators = [
        "endothelial",
        "stromal_fibroblast",
        "dct_identity",
        "matrix_endothelial_composite",
    ]
    summary_rows = []
    posterior_frames = []
    power_rows = []
    for idx, mediator in enumerate(mediators):
        df = merged[["flight", mediator, "ncc_activity_score_regulatory"]].dropna().copy()
        if len(df) < 10:
            continue
        M = (df[mediator] - df[mediator].mean()) / df[mediator].std(ddof=1)
        Y = (df["ncc_activity_score_regulatory"] - df["ncc_activity_score_regulatory"].mean()) / df[
            "ncc_activity_score_regulatory"
        ].std(ddof=1)
        X_med = np.column_stack([np.ones(len(df)), df["flight"].to_numpy()])
        a_draws, sigma_m, a_hat, resid_m = approximate_bayes_linear(M, X_med, seed=20260526 + idx)
        X_out = np.column_stack([np.ones(len(df)), df["flight"].to_numpy(), M.to_numpy()])
        b_draws, sigma_y, b_hat, resid_y = approximate_bayes_linear(Y, X_out, seed=20260626 + idx)
        a = a_draws[:, 1]
        c_prime = b_draws[:, 1]
        b = b_draws[:, 2]
        indirect = a * b
        total = c_prime + indirect
        for name, draws in [("a", a), ("b", b), ("c_prime", c_prime), ("indirect", indirect), ("total", total)]:
            summary_rows.append(
                {
                    "mediator": mediator,
                    "parameter": name,
                    "posterior_median": np.median(draws),
                    "ci_low": np.quantile(draws, 0.025),
                    "ci_high": np.quantile(draws, 0.975),
                    "p_less_than_zero": np.mean(draws < 0),
                    "p_greater_than_zero": np.mean(draws > 0),
                    "n_animals": len(df),
                    "model": "approximate_bayesian_ols_weak_prior_fallback",
                }
            )
        posterior_frames.append(
            pd.DataFrame(
                {
                    "mediator": mediator,
                    "draw": np.arange(len(indirect)),
                    "a": a,
                    "b": b,
                    "c_prime": c_prime,
                    "indirect": indirect,
                    "total": total,
                }
            )
        )

        # Simple future-n simulation using posterior median paths and residual SDs.
        path_a = float(np.median(a))
        path_b = float(np.median(b))
        path_c = float(np.median(c_prime))
        sd_m = float(np.median(sigma_m))
        sd_y = float(np.median(sigma_y))
        rng = np.random.default_rng(20260726 + idx)
        for n_total in [20, 30, 40, 60, 80, 100, 140]:
            sign_hits = 0
            ci_hits = 0
            sims = 300
            n1 = n_total // 2
            n0 = n_total - n1
            xsim = np.r_[np.zeros(n0), np.ones(n1)]
            for _ in range(sims):
                msim = path_a * xsim + rng.normal(0, sd_m, size=n_total)
                ysim = path_c * xsim + path_b * msim + rng.normal(0, sd_y, size=n_total)
                xm = np.column_stack([np.ones(n_total), xsim])
                xy = np.column_stack([np.ones(n_total), xsim, msim])
                beta_m = np.linalg.lstsq(xm, msim, rcond=None)[0]
                beta_y = np.linalg.lstsq(xy, ysim, rcond=None)[0]
                resid_m_sim = msim - xm @ beta_m
                resid_y_sim = ysim - xy @ beta_y
                cov_m = np.sum(resid_m_sim**2) / max(n_total - xm.shape[1], 1) * np.linalg.pinv(xm.T @ xm)
                cov_y = np.sum(resid_y_sim**2) / max(n_total - xy.shape[1], 1) * np.linalg.pinv(xy.T @ xy)
                ad = beta_m[1]
                bd = beta_y[2]
                se_a = math.sqrt(max(cov_m[1, 1], 0))
                se_b = math.sqrt(max(cov_y[2, 2], 0))
                ind = ad * bd
                se_ind = math.sqrt((bd**2) * (se_a**2) + (ad**2) * (se_b**2))
                lo = ind - 1.96 * se_ind
                hi = ind + 1.96 * se_ind
                sign_hits += int(ind < 0)
                ci_hits += int(ind < 0 and hi < 0)
            power_rows.append(
                {
                    "mediator": mediator,
                    "n_total": n_total,
                    "directional_sign_recovery": sign_hits / sims,
                    "ci_exclusion_power_estimate": ci_hits / sims,
                    "simulation_note": "Sobel-style interval exclusion for posterior-median effect; still approximate and not causal proof",
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(root / "h3_mediation" / "h3_mediation_model_summary.tsv", sep="\t", index=False)
    if posterior_frames:
        pd.concat(posterior_frames, ignore_index=True).to_csv(
            root / "h3_mediation" / "h3_mediation_posterior_draws.tsv.gz", sep="\t", index=False, compression="gzip"
        )
    pd.DataFrame(power_rows).to_csv(root / "h3_mediation" / "h3_mediation_power_simulation.tsv", sep="\t", index=False)
    indirect_rows = summary[summary["parameter"] == "indirect"].copy()
    verdict = {
        "model_caveat": "Approximate Bayesian linear-model fallback; cross-sectional n=20 bulk tissue, not causal proof.",
        "mediators": indirect_rows.to_dict(orient="records"),
        "overall_interpretation": "Use as mediation-specified causal-structure bound. Wide intervals or composition controls prevent mechanistic claims.",
    }
    (root / "h3_mediation" / "h3_mediation_verdict.json").write_text(json.dumps(verdict, indent=2))


def write_manifest(root: Path):
    input_paths = [
        Path("docs/v11_execution_research_plan.md"),
        DCT_PRIOR_DIR / "gse228367_dct1_vs_dct2_de.tsv",
        OSD462 / "phospho_all_sites.tsv",
        OSD462 / "protein_effects_gene_anyplex.tsv",
        PHENO / "phenotype_anchor_per_animal.tsv",
        CELLTYPE / "osd462_compartment_scores_per_sample.tsv",
        PXD_DIR / "41598_2015_BFsrep12829_MOESM6_ESM.xls",
    ]
    manifest = {
        "analysis": "v11 core DCT1/phosphoproteome/mediation analysis",
        "run_root": str(root),
        "claim_discipline": "GSE228367 and PXD001729 are reference priors/scaffolds, not spaceflight evidence.",
        "inputs": [{"path": str(p), "sha256": sha256(p)} for p in input_paths if p.exists()],
    }
    (root / "manifests" / "v11_core_manifest.json").write_text(json.dumps(manifest, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default=str(RUN_ROOT))
    args = parser.parse_args()
    root = Path(args.run_root)
    ensure_dirs(root)

    baseline_lock(root)
    pxd = load_pxd_tables(root)
    _, proteins_prior, phospho_prior = build_dct_prior_mapping(root)
    run_h2_enrichment(root, phospho_prior)
    run_pxd_antialignment(root, pxd, phospho_prior)
    run_klhl3_check(root, pxd, phospho_prior, proteins_prior)
    run_mediation(root)
    write_manifest(root)
    print(f"v11 core analysis complete: {root}")


if __name__ == "__main__":
    main()
