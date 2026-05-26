#!/usr/bin/env python3
"""Composition-aware H2 robustness tests for OSD-462 phosphoproteomics.

This script tests whether the v11 DCT1-prior phosphosite signal survives a
conservative adjustment ladder using per-animal TMT phosphosite intensities.

Claim discipline:
  * This is not DCT-specific phosphoproteomics.
  * Composition scores are bulk RNA marker-panel estimates and may sit on the
    biological path from flight to phosphosite dilution.
  * The output should be described as composition-aware robustness, not
    deconvolution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.multiomics.osd462_anchor import parse_tmt_sheet  # noqa: E402


RUN_ROOT = REPO_ROOT / "data/results/run_20260526_v11_dct1_phospho_mediation"
OSD462_DIR = REPO_ROOT / "data/external/osdr/OSD-462"
PHOSPHO_XLSX = OSD462_DIR / "GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx"
PROTEIN_XLSX = OSD462_DIR / "GLDS-462_proteomics_2021-12-31_tc884-885_Protein_WorkUp.xlsx"
CELLTYPE = REPO_ROOT / "data/results/run_20260522_celltype_decomposition/osd462_compartment_scores_per_sample.tsv"
DCT_PRIOR = RUN_ROOT / "dct_prior/osd462_phosphosite_dct1_prior.tsv"
OUT_DIR = RUN_ROOT / "h2_composition_adjusted"

SEED = 20260526
EPS = 1e-12


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def zscore(s: pd.Series) -> pd.Series:
    x = pd.to_numeric(s, errors="coerce")
    mu = x.mean(skipna=True)
    sd = x.std(skipna=True, ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return x * np.nan
    return (x - mu) / sd


def tmt_key(condition: str, sample_label: str) -> str | None:
    m = re.match(r"^(BL|FL|GC)-0?(\d+)$", str(sample_label).strip())
    if not m:
        return None
    cond = {"BL": "BSL", "FL": "FLT", "GC": "GC"}[m.group(1)]
    return f"{cond}|{int(m.group(2))}"


def rna_sample_key(sample_name: str) -> str | None:
    m = re.search(r"RR10_KDN_WT_(FLT|GC|BSL|VIV)_([A-Z])(\d+)", str(sample_name))
    if not m:
        return None
    return f"{m.group(1)}|{int(m.group(3))}"


def channel_metadata(table) -> pd.DataFrame:
    meta = table.channels.copy()
    meta["sample_key"] = [tmt_key(c, s) for c, s in zip(meta["condition"], meta["sample"])]
    meta["flight"] = (meta["condition"] == "FL").astype(float)
    meta["plex2"] = (meta["plex"] == "Samp6-10").astype(float)
    return meta


def centered_log2_matrix(table) -> pd.DataFrame:
    """log2 scaled TMT values with per-channel median centering inside plex."""
    cols = list(table.channels["column"])
    mat = table.scaled[cols].to_numpy(dtype=float)
    mat = np.where(mat > 0, mat, np.nan)
    log2 = np.log2(mat)
    centers = np.zeros(log2.shape[1], dtype=float)
    for plex in sorted(table.channels["plex"].unique()):
        idx = np.flatnonzero(table.channels["plex"].to_numpy() == plex)
        centers[idx] = np.nanmedian(log2[:, idx], axis=0)
    log2 = log2 - centers[None, :]
    return pd.DataFrame(log2, columns=cols)


def load_compartment_scores() -> pd.DataFrame:
    scores = pd.read_csv(CELLTYPE, sep="\t", index_col=0).T
    rows = []
    for sample, vals in scores.iterrows():
        key = rna_sample_key(sample)
        if key is None:
            continue
        rows.append(
            {
                "sample_key": key,
                "dct_identity_score": vals.get("dct_identity", np.nan),
                "endothelial_score": vals.get("endothelial", np.nan),
                "stromal_score": vals.get("stromal_fibroblast", np.nan),
            }
        )
    out = pd.DataFrame(rows)
    out = out.groupby("sample_key", as_index=False).mean(numeric_only=True)
    comp = out[["dct_identity_score", "endothelial_score", "stromal_score"]].copy()
    comp = comp.apply(zscore)
    comp = comp.fillna(0.0)
    u, s, _ = np.linalg.svd(comp.to_numpy(dtype=float), full_matrices=False)
    out["composition_pc1"] = u[:, 0] * s[0]
    for c in ["dct_identity_score", "endothelial_score", "stromal_score", "composition_pc1"]:
        out[f"{c}_z"] = zscore(out[c])
    return out


def load_protein_long() -> pd.DataFrame:
    tab = parse_tmt_sheet(
        PROTEIN_XLSX,
        "protein_quant_2721",
        gene_col="Gene Symbol",
        peptide_cols={"Samp1-5": "Samp1-5 Peptides", "Samp6-10": "Samp6-10 Peptides"},
        id_col="Protein Id",
    )
    values = centered_log2_matrix(tab)
    ch = channel_metadata(tab)
    sample_cols = ch[ch["condition"].isin(["FL", "GC"])].copy()
    sample_cols = sample_cols[sample_cols["sample_key"].notna()].copy()

    meta = tab.meta.copy()
    meta["row_weight"] = pd.to_numeric(meta.get("n_peptides", np.nan), errors="coerce")
    if "n_peptides" not in meta.columns:
        pep1 = pd.to_numeric(meta.get("n_pep_Samp1-5", 0), errors="coerce").fillna(0)
        pep2 = pd.to_numeric(meta.get("n_pep_Samp6-10", 0), errors="coerce").fillna(0)
        meta["row_weight"] = pep1 + pep2
    meta["row_weight"] = meta["row_weight"].fillna(1).clip(lower=1)

    blocks = []
    for _, c in sample_cols.iterrows():
        d = pd.DataFrame(
            {
                "gene_symbol": meta["gene_symbol"].astype(str),
                "sample_key": c["sample_key"],
                "protein_abundance": values[c["column"]].to_numpy(dtype=float),
                "weight": meta["row_weight"].to_numpy(dtype=float),
            }
        )
        blocks.append(d)
    long = pd.concat(blocks, ignore_index=True).dropna(subset=["protein_abundance"])

    def wavg(g):
        w = g["weight"].to_numpy(dtype=float)
        x = g["protein_abundance"].to_numpy(dtype=float)
        return pd.Series({"parent_protein_abundance": float(np.average(x, weights=w))})

    return long.groupby(["gene_symbol", "sample_key"], sort=False).apply(wavg).reset_index()


def load_phosphosite_long(site_scope: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    configs = [("siteQuant_360", "gene_symbol", "Site Position", "Protein Id", "single")]
    if site_scope == "all":
        configs.append(("siteQuant_360_compositeSite", "geneSymbol", "sitePosStr", "proteinID", "composite"))

    all_long = []
    all_meta = []
    for sheet, gcol, scol, idcol, kind in configs:
        tab = parse_tmt_sheet(
            PHOSPHO_XLSX,
            sheet,
            gene_col=gcol,
            peptide_cols={"Samp1-5": "Samp1-5~num_quant", "Samp6-10": "Samp6-10~num_quant"},
            id_col=idcol,
            extra_meta_cols=[scol],
        )
        values = centered_log2_matrix(tab)
        ch = channel_metadata(tab)
        sample_cols = ch[ch["condition"].isin(["FL", "GC"])].copy()
        sample_cols = sample_cols[sample_cols["sample_key"].notna()].copy()
        meta = tab.meta.copy().rename(columns={scol: "site_position"})
        meta["site_kind"] = kind
        meta["sheet"] = sheet
        meta["site_position_str"] = meta["site_position"].astype(str)
        meta["site_row_id"] = [
            f"{sheet}|{i}|{g}|{p}"
            for i, (g, p) in enumerate(zip(meta["gene_symbol"], meta["site_position_str"]))
        ]
        meta["mean_intensity"] = values[list(sample_cols["column"])].mean(axis=1, skipna=True)
        meta["n_valid_samples"] = values[list(sample_cols["column"])].notna().sum(axis=1)
        all_meta.append(
            meta[
                [
                    "site_row_id",
                    "gene_symbol",
                    "site_position_str",
                    "site_kind",
                    "sheet",
                    "protein_id",
                    "mean_intensity",
                    "n_valid_samples",
                ]
            ]
        )
        for _, c in sample_cols.iterrows():
            d = pd.DataFrame(
                {
                    "site_row_id": meta["site_row_id"],
                    "gene_symbol": meta["gene_symbol"].astype(str),
                    "site_position_str": meta["site_position_str"],
                    "site_kind": kind,
                    "sample_key": c["sample_key"],
                    "flight": c["flight"],
                    "plex2": c["plex2"],
                    "phosphosite_abundance": values[c["column"]].to_numpy(dtype=float),
                }
            )
            all_long.append(d)
    long = pd.concat(all_long, ignore_index=True).dropna(subset=["phosphosite_abundance"])
    meta = pd.concat(all_meta, ignore_index=True)
    return long, meta


def attach_site_priors(site_meta: pd.DataFrame) -> pd.DataFrame:
    prior = pd.read_csv(DCT_PRIOR, sep="\t")
    gene_prior = (
        prior.sort_values(["gene_symbol", "site_kind"])
        .groupby("gene_symbol", as_index=False)
        .agg(
            dct1_enrichment_score=("dct1_enrichment_score", "first"),
            dct1_top_quartile=("dct1_top_quartile", "max"),
            dct1_top_decile=("dct1_top_decile", "max"),
            dct2_bottom_quartile=("dct2_bottom_quartile", "max"),
            parent_protein_flight_effect=("flight_effect", "first"),
            parent_abundance_log2=("abundance_log2", "first"),
            parent_n_peptides=("n_peptides", "first"),
        )
    )
    out = site_meta.merge(gene_prior, on="gene_symbol", how="left")
    out["n_sites_on_parent_gene"] = out.groupby("gene_symbol")["site_row_id"].transform("count")
    out["dct1_enrichment_score_z"] = zscore(out["dct1_enrichment_score"])
    out["mean_intensity_z"] = zscore(out["mean_intensity"])
    out["missingness"] = 20 - pd.to_numeric(out["n_valid_samples"], errors="coerce")
    out["missingness_z"] = zscore(out["missingness"])
    out["n_sites_on_parent_gene_log_z"] = zscore(np.log1p(out["n_sites_on_parent_gene"]))
    out["parent_protein_flight_effect_z"] = zscore(out["parent_protein_flight_effect"])
    return out


def first_stage_models(long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    model_covars = {
        "M0_raw": [],
        "M1_parent_protein_abundance": ["parent_protein_abundance_z"],
        "M2_dct_score": ["dct_identity_score_z"],
        "M3_endothelial_stromal": ["endothelial_score_z", "stromal_score_z"],
        "M4_full_parent_dct_endothelial_stromal": [
            "parent_protein_abundance_z",
            "dct_identity_score_z",
            "endothelial_score_z",
            "stromal_score_z",
        ],
        "M5_parent_composition_pc1": ["parent_protein_abundance_z", "composition_pc1_z"],
    }
    rows = []
    design_rows = []
    grouped = long.groupby("site_row_id", sort=False)
    for site_id, g in grouped:
        for model, covars in model_covars.items():
            use_cols = ["phosphosite_abundance", "flight", "plex2"] + covars
            d = g.dropna(subset=use_cols)
            n = len(d)
            p = 3 + len(covars)
            if n <= p + 1 or d["flight"].nunique() < 2:
                continue
            y = d["phosphosite_abundance"].to_numpy(dtype=float)
            X = np.column_stack(
                [
                    np.ones(n),
                    d["flight"].to_numpy(dtype=float),
                    d["plex2"].to_numpy(dtype=float),
                ]
                + [d[c].to_numpy(dtype=float) for c in covars]
            )
            try:
                beta, *_ = np.linalg.lstsq(X, y, rcond=None)
                resid = y - X @ beta
                df = n - X.shape[1]
                if df <= 0:
                    continue
                sigma2 = float((resid @ resid) / df)
                xtxi = np.linalg.pinv(X.T @ X)
                se = math.sqrt(max(sigma2 * xtxi[1, 1], 0))
                tval = beta[1] / se if se > 0 else np.nan
                pval = float(2 * stats.t.sf(abs(tval), df)) if np.isfinite(tval) else np.nan
            except np.linalg.LinAlgError:
                continue
            rows.append(
                {
                    "site_row_id": site_id,
                    "model": model,
                    "adjusted_flight_effect": float(beta[1]),
                    "adjusted_flight_se": float(se),
                    "adjusted_flight_p": pval,
                    "n_obs": int(n),
                    "n_parameters": int(X.shape[1]),
                }
            )
            design_rows.append({"model": model, "covariates": ",".join(covars) if covars else "none"})
    effects = pd.DataFrame(rows)
    design = pd.DataFrame(design_rows).drop_duplicates().sort_values("model")
    if len(effects):
        effects["adjusted_flight_q_all_sites"] = effects.groupby("model")["adjusted_flight_p"].transform(bh)
    return effects, design


def effect_level_regression(effects: pd.DataFrame, site_meta: pd.DataFrame) -> pd.DataFrame:
    df = effects.merge(site_meta, on="site_row_id", how="left")
    covars_full = [
        "dct1_enrichment_score_z",
        "parent_protein_flight_effect_z",
        "mean_intensity_z",
        "missingness_z",
        "n_sites_on_parent_gene_log_z",
    ]
    covars_dct_only = ["dct1_enrichment_score_z"]
    rows = []
    for model, d0 in df.groupby("model", sort=False):
        for adjustment, covars in [
            ("dct1_only_second_stage", covars_dct_only),
            ("full_second_stage", covars_full),
        ]:
            d = d0.dropna(subset=["adjusted_flight_effect", "gene_symbol"] + covars).copy()
            if len(d) < 50:
                continue
            X = sm.add_constant(d[covars], has_constant="add")
            y = d["adjusted_flight_effect"].astype(float)
            fit = sm.OLS(y, X).fit()
            try:
                fit = fit.get_robustcov_results(cov_type="cluster", groups=d["gene_symbol"].astype(str))
                names = list(X.columns)
                params = pd.Series(fit.params, index=names)
                pvals = pd.Series(fit.pvalues, index=names)
                ses = pd.Series(fit.bse, index=names)
            except Exception:
                fit = sm.OLS(y, X).fit(cov_type="HC3")
                params = fit.params
                pvals = fit.pvalues
                ses = fit.bse
            rows.append(
                {
                    "model": model,
                    "second_stage_adjustment": adjustment,
                    "term": "dct1_enrichment_score_z",
                    "coefficient": float(params["dct1_enrichment_score_z"]),
                    "se": float(ses["dct1_enrichment_score_z"]),
                    "p_value": float(pvals["dct1_enrichment_score_z"]),
                    "n_sites": int(len(d)),
                    "n_parent_genes": int(d["gene_symbol"].nunique()),
                    "interpretation": "negative means higher DCT1 prior predicts stronger flight suppression",
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value"] = bh(out["p_value"])
    return out


def adjusted_enrichment(effects: pd.DataFrame, site_meta: pd.DataFrame) -> pd.DataFrame:
    df = effects.merge(site_meta, on="site_row_id", how="left")
    rows = []
    for model, d0 in df.groupby("model", sort=False):
        d = d0.dropna(subset=["adjusted_flight_effect", "adjusted_flight_p", "dct1_top_decile", "dct1_top_quartile"])
        if len(d) < 50:
            continue
        suppressed = (d["adjusted_flight_effect"] < 0) & (d["adjusted_flight_p"] < 0.05)
        for flag in ["dct1_top_quartile", "dct1_top_decile"]:
            in_flag = d[flag].astype(bool)
            table = [
                [int((suppressed & in_flag).sum()), int((suppressed & ~in_flag).sum())],
                [int((~suppressed & in_flag).sum()), int((~suppressed & ~in_flag).sum())],
            ]
            odds, p = stats.fisher_exact(table, alternative="greater")
            rows.append(
                {
                    "model": model,
                    "flag": flag,
                    "n_sites": int(len(d)),
                    "n_suppressed_nominal": int(suppressed.sum()),
                    "suppressed_in_flag": table[0][0],
                    "suppressed_not_flag": table[0][1],
                    "background_in_flag": table[1][0],
                    "background_not_flag": table[1][1],
                    "odds_ratio": float(odds),
                    "p_value": float(p),
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value"] = bh(out["p_value"])
    return out


def site_fixed_long_models(long: pd.DataFrame, site_meta: pd.DataFrame) -> pd.DataFrame:
    df = long.copy()
    if "dct1_enrichment_score_z" not in df.columns:
        meta = site_meta[["site_row_id", "gene_symbol", "dct1_enrichment_score_z"]].copy()
        df = df.merge(meta, on=["site_row_id", "gene_symbol"], how="left")
    df["flight_x_dct1_prior"] = df["flight"] * df["dct1_enrichment_score_z"]
    model_covars = {
        "LM0_raw_site_fixed": ["flight", "flight_x_dct1_prior", "plex2"],
        "LM1_parent_site_fixed": ["flight", "flight_x_dct1_prior", "plex2", "parent_protein_abundance_z"],
        "LM2_dct_site_fixed": ["flight", "flight_x_dct1_prior", "plex2", "dct_identity_score_z"],
        "LM3_endothelial_stromal_site_fixed": [
            "flight",
            "flight_x_dct1_prior",
            "plex2",
            "endothelial_score_z",
            "stromal_score_z",
        ],
        "LM4_full_site_fixed": [
            "flight",
            "flight_x_dct1_prior",
            "plex2",
            "parent_protein_abundance_z",
            "dct_identity_score_z",
            "endothelial_score_z",
            "stromal_score_z",
        ],
        "LM5_parent_composition_pc1_site_fixed": [
            "flight",
            "flight_x_dct1_prior",
            "plex2",
            "parent_protein_abundance_z",
            "composition_pc1_z",
        ],
    }
    rows = []
    for model, covars in model_covars.items():
        d = df.dropna(subset=["phosphosite_abundance", "gene_symbol"] + covars).copy()
        if len(d) < 1000:
            continue
        cols = ["phosphosite_abundance"] + covars
        # Absorb phosphosite baseline by within-site demeaning. This is the
        # scalable fixed-effect analogue of a random phosphosite intercept.
        means = d.groupby("site_row_id")[cols].transform("mean")
        yd = d["phosphosite_abundance"] - means["phosphosite_abundance"]
        Xd = d[covars] - means[covars]
        nonzero = Xd.std(axis=0) > 1e-10
        Xd = Xd.loc[:, nonzero]
        fit = sm.OLS(yd, Xd).fit()
        try:
            fit = fit.get_robustcov_results(cov_type="cluster", groups=d["gene_symbol"].astype(str))
            names = list(Xd.columns)
            params = pd.Series(fit.params, index=names)
            pvals = pd.Series(fit.pvalues, index=names)
            ses = pd.Series(fit.bse, index=names)
        except Exception:
            fit = sm.OLS(yd, Xd).fit(cov_type="HC3")
            params = fit.params
            pvals = fit.pvalues
            ses = fit.bse
        term = "flight_x_dct1_prior"
        rows.append(
            {
                "model": model,
                "term": term,
                "coefficient": float(params.get(term, np.nan)),
                "se": float(ses.get(term, np.nan)),
                "p_value": float(pvals.get(term, np.nan)),
                "n_observations": int(len(d)),
                "n_sites": int(d["site_row_id"].nunique()),
                "n_parent_genes": int(d["gene_symbol"].nunique()),
                "interpretation": "negative interaction means flight suppresses higher-DCT1-prior sites more strongly",
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value"] = bh(out["p_value"])
    return out


def build_verdict(effect_summary: pd.DataFrame, long_summary: pd.DataFrame, enrichment: pd.DataFrame) -> dict:
    full_eff = effect_summary[
        (effect_summary["model"] == "M4_full_parent_dct_endothelial_stromal")
        & (effect_summary["second_stage_adjustment"] == "full_second_stage")
    ]
    full_long = long_summary[long_summary["model"] == "LM4_full_site_fixed"]
    eff_pass = bool(
        len(full_eff)
        and full_eff["coefficient"].iloc[0] < 0
        and full_eff["q_value"].iloc[0] <= 0.10
    )
    long_pass = bool(
        len(full_long)
        and full_long["coefficient"].iloc[0] < 0
        and full_long["q_value"].iloc[0] <= 0.10
    )
    full_top_decile = enrichment[
        (enrichment["model"] == "M4_full_parent_dct_endothelial_stromal")
        & (enrichment["flag"] == "dct1_top_decile")
    ]
    top_decile_pass = bool(
        len(full_top_decile)
        and full_top_decile["odds_ratio"].iloc[0] > 1
        and full_top_decile["q_value"].iloc[0] <= 0.10
    )
    if eff_pass and long_pass:
        interp = (
            "DCT1-prior association remains negative after parent-protein and conservative composition adjustment "
            "in both two-stage and site-fixed models."
        )
    elif top_decile_pass and not (eff_pass or long_pass):
        interp = (
            "DCT1 top-decile adjusted-suppression enrichment survives the full conservative model, "
            "but continuous DCT1-prior coefficients and the site-fixed interaction are weak. Claim a robust "
            "DCT1-high subset enrichment, not a broad continuous DCT1-gradient effect."
        )
    elif (len(full_eff) and full_eff["coefficient"].iloc[0] < 0) or (
        len(full_long) and full_long["coefficient"].iloc[0] < 0
    ):
        interp = (
            "DCT1-prior association remains directionally negative in at least one conservative model, "
            "but evidence is attenuated or model-dependent."
        )
    else:
        interp = "Conservative adjustment does not preserve a negative DCT1-prior association."
    return {
        "claim_label": "composition-aware robustness, not DCT-specific deconvolution",
        "full_two_stage_pass_fdr_0_10": eff_pass,
        "full_site_fixed_pass_fdr_0_10": long_pass,
        "full_adjusted_top_decile_enrichment_pass_fdr_0_10": top_decile_pass,
        "interpretation": interp,
        "caveat": (
            "DCT/endothelial/stromal scores are animal-level bulk RNA estimates; adjusting for them may remove "
            "biological signal if tissue remodeling is on the causal path."
        ),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=RUN_ROOT)
    ap.add_argument("--site-scope", choices=["single", "all"], default="single")
    args = ap.parse_args()

    out_dir = args.run_root / "h2_composition_adjusted"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("[h2-composition] loading per-sample compartment scores")
    comp = load_compartment_scores()
    print("[h2-composition] loading parent protein abundance matrix")
    prot_long = load_protein_long()
    print("[h2-composition] loading phosphosite abundance matrix")
    phospho_long, site_meta = load_phosphosite_long(args.site_scope)
    site_meta = attach_site_priors(site_meta)

    long = phospho_long.merge(comp, on="sample_key", how="left")
    long = long.merge(prot_long, on=["gene_symbol", "sample_key"], how="left")
    long["parent_protein_abundance_z"] = zscore(long["parent_protein_abundance"])

    keep_meta = site_meta[["site_row_id", "dct1_enrichment_score_z"]]
    long = long.merge(keep_meta, on="site_row_id", how="left")
    long = long.dropna(subset=["dct1_enrichment_score_z", "phosphosite_abundance"])

    print(f"[h2-composition] long rows={len(long):,}; sites={long['site_row_id'].nunique():,}")
    site_meta.to_csv(out_dir / f"h2_composition_site_metadata_{args.site_scope}.tsv", sep="\t", index=False)

    print("[h2-composition] fitting per-site model ladder")
    effects, design = first_stage_models(long)
    effects.to_csv(
        out_dir / f"h2_composition_adjusted_site_effects_{args.site_scope}.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    design.to_csv(out_dir / f"h2_composition_model_ladder_{args.site_scope}.tsv", sep="\t", index=False)

    print("[h2-composition] fitting effect-level DCT1 regressions")
    effect_summary = effect_level_regression(effects, site_meta)
    effect_summary.to_csv(
        out_dir / f"h2_composition_effect_level_dct1_ladder_{args.site_scope}.tsv",
        sep="\t",
        index=False,
    )

    print("[h2-composition] computing adjusted suppression enrichment")
    enrichment = adjusted_enrichment(effects, site_meta)
    enrichment.to_csv(
        out_dir / f"h2_composition_adjusted_suppression_enrichment_{args.site_scope}.tsv",
        sep="\t",
        index=False,
    )

    print("[h2-composition] fitting one-shot site-fixed interaction models")
    long_summary = site_fixed_long_models(long, site_meta)
    long_summary.to_csv(
        out_dir / f"h2_composition_site_fixed_interaction_ladder_{args.site_scope}.tsv",
        sep="\t",
        index=False,
    )

    verdict = build_verdict(effect_summary, long_summary, enrichment)
    verdict["site_scope"] = args.site_scope
    (out_dir / f"h2_composition_aware_verdict_{args.site_scope}.json").write_text(json.dumps(verdict, indent=2))

    manifest = {
        "analysis": "H2 composition-aware phosphosite robustness",
        "site_scope": args.site_scope,
        "inputs": [
            {"path": str(PHOSPHO_XLSX), "sha256": sha256(PHOSPHO_XLSX)},
            {"path": str(PROTEIN_XLSX), "sha256": sha256(PROTEIN_XLSX)},
            {"path": str(CELLTYPE), "sha256": sha256(CELLTYPE)},
            {"path": str(DCT_PRIOR), "sha256": sha256(DCT_PRIOR)},
        ],
        "outputs": {
            "site_metadata": str(out_dir / f"h2_composition_site_metadata_{args.site_scope}.tsv"),
            "site_effects": str(out_dir / f"h2_composition_adjusted_site_effects_{args.site_scope}.tsv.gz"),
            "effect_level": str(out_dir / f"h2_composition_effect_level_dct1_ladder_{args.site_scope}.tsv"),
            "site_fixed": str(out_dir / f"h2_composition_site_fixed_interaction_ladder_{args.site_scope}.tsv"),
            "verdict": str(out_dir / f"h2_composition_aware_verdict_{args.site_scope}.json"),
        },
        "claim_discipline": "composition-aware robustness; not deconvolved or DCT-isolated phosphoproteomics",
    }
    (out_dir / f"h2_composition_manifest_{args.site_scope}.json").write_text(json.dumps(manifest, indent=2))
    print(f"[h2-composition] complete: {out_dir}")


if __name__ == "__main__":
    main()
