#!/usr/bin/env python3
"""Repair C -- DCT1<->DCT2 *continuous* contrast for flight-suppressed phosphosites."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.v11.core_analysis import (
    ANCHOR_GENES,
    RUN_ROOT,
    bh,
    build_parent_gene_table,
    read_tsv,
)

# Parent-gene-level covariates (already z-scored by build_parent_gene_table).
COVARIATE_COLS = [
    "n_quantified_phosphosites_z",
    "parent_n_peptides_z",
    "parent_abundance_log2_z",
    "mean_missing_samples_z",
]
COORD_COL = "dct1_enrichment_score"
MIN_PARENT_GENES = 30


# #
def _aligned(coord, outcome, covariates=None):
    """Return a clean numeric design frame (y, coord, covars) with NA rows dropped."""
    y = pd.to_numeric(pd.Series(outcome).reset_index(drop=True), errors="coerce")
    x = pd.to_numeric(pd.Series(coord).reset_index(drop=True), errors="coerce")
    frame = pd.DataFrame({"y": y, "coord": x})
    cov_cols: list[str] = []
    if covariates is not None:
        cov = pd.DataFrame(covariates).reset_index(drop=True)
        for c in cov.columns:
            frame[c] = pd.to_numeric(cov[c], errors="coerce")
            cov_cols.append(str(c))
    return frame.dropna(), cov_cols


def fit_linear_gradient(coord, outcome, covariates=None, min_n: int = MIN_PARENT_GENES) -> dict:
    """OLS slope of ``outcome`` on ``coord`` (+ optional covariates).

    The slope is reported in the units of ``coord`` as passed in: standardize
    ``coord`` before calling for a per-SD slope. A planted positive trend is
    recovered as a positive slope (sign-faithful), which is what the unit test
    checks.
    """
    try:
        import statsmodels.api as sm
    except Exception as exc:  # pragma: no cover - optional dependency
        return {"model_status": f"not_fit: {exc}"}

    d, cov_cols = _aligned(coord, outcome, covariates)
    if len(d) < min_n or d["coord"].nunique() < 3 or d["y"].nunique() < 3:
        return {"model_status": "not_fit: insufficient variation", "n": int(len(d))}
    X = sm.add_constant(d[["coord"] + cov_cols].astype(float), has_constant="add")
    try:
        fit = sm.OLS(d["y"].astype(float), X).fit()
        slope = float(fit.params["coord"])
        se = float(fit.bse["coord"])
        return {
            "model_status": "fit",
            "slope": slope,
            "se": se,
            "t_stat": float(fit.tvalues["coord"]),
            "p_value": float(fit.pvalues["coord"]),
            "ci_low": float(fit.conf_int().loc["coord", 0]),
            "ci_high": float(fit.conf_int().loc["coord", 1]),
            "r_squared": float(fit.rsquared),
            "n": int(len(d)),
            "covariate_adjusted": bool(cov_cols),
        }
    except Exception as exc:  # pragma: no cover
        return {"model_status": f"not_fit: {exc}", "n": int(len(d))}


def fit_logistic_gradient(coord, indicator, covariates=None, min_n: int = MIN_PARENT_GENES) -> dict:
    """GLM-binomial log-odds slope of a binary ``indicator`` on ``coord``."""
    try:
        import statsmodels.api as sm
    except Exception as exc:  # pragma: no cover
        return {"model_status": f"not_fit: {exc}"}

    y = pd.Series(indicator).reset_index(drop=True)
    y = y.map({True: 1.0, False: 0.0}).where(y.isin([True, False]), pd.to_numeric(y, errors="coerce"))
    d, cov_cols = _aligned(coord, y, covariates)
    if len(d) < min_n or d["coord"].nunique() < 3 or d["y"].nunique() < 2:
        return {"model_status": "not_fit: insufficient variation", "n": int(len(d))}
    X = sm.add_constant(d[["coord"] + cov_cols].astype(float), has_constant="add")
    try:
        fit = sm.GLM(d["y"].astype(float), X, family=sm.families.Binomial()).fit()
        coef = float(fit.params["coord"])
        se = float(fit.bse["coord"])
        return {
            "model_status": "fit",
            "log_odds": coef,
            "se": se,
            "odds_ratio": float(np.exp(coef)),
            "z_stat": float(coef / se) if se > 0 else np.nan,
            "p_value": float(fit.pvalues["coord"]),
            "ci_low": float(np.exp(coef - 1.96 * se)),
            "ci_high": float(np.exp(coef + 1.96 * se)),
            "n": int(len(d)),
            "n_positive": int(d["y"].sum()),
            "covariate_adjusted": bool(cov_cols),
        }
    except Exception as exc:  # pragma: no cover
        return {"model_status": f"not_fit: {exc}", "n": int(len(d))}


def spearman_gradient(coord, outcome, min_n: int = MIN_PARENT_GENES) -> dict:
    """Rank-based (Spearman) association of ``coord`` with ``outcome``."""
    d, _ = _aligned(coord, outcome)
    if len(d) < min_n or d["coord"].nunique() < 3 or d["y"].nunique() < 3:
        return {"model_status": "not_fit: insufficient variation", "n": int(len(d))}
    rho, p = stats.spearmanr(d["coord"], d["y"])
    return {
        "model_status": "fit",
        "spearman_rho": float(rho),
        "p_value": float(p),
        "n": int(len(d)),
    }


def fit_spline_gradient(
    coord, outcome, covariates=None, df_spline: int = 4, min_n: int = MIN_PARENT_GENES
) -> dict:
    """Natural-cubic-spline non-linearity test (nested F vs the linear contrast).

    Fits a reduced model (const + coord [+ covars]) and a full model
    (const + natural-cubic-spline(coord, df) [+ covars]). The natural cubic
    spline space contains all linear functions, so the reduced model is nested
    in the full one and a standard partial F-test is valid. A small ``f_p``
    means the gradient departs from a straight line (non-monotone or curved).
    """
    try:
        import patsy
        import statsmodels.api as sm
    except Exception as exc:  # pragma: no cover
        return {"model_status": f"not_fit: {exc}"}

    d, cov_cols = _aligned(coord, outcome, covariates)
    if len(d) < min_n or d["coord"].nunique() < (df_spline + 1) or d["y"].nunique() < 3:
        return {"model_status": "not_fit: insufficient variation", "n": int(len(d))}
    try:
        cov = d[cov_cols].astype(float) if cov_cols else None
        # Reduced: const + coord (+ covars)
        Xr = sm.add_constant(
            pd.concat([d[["coord"]].astype(float).reset_index(drop=True),
                       cov.reset_index(drop=True)] if cov is not None
                      else [d[["coord"]].astype(float).reset_index(drop=True)], axis=1),
            has_constant="add",
        )
        reduced = sm.OLS(d["y"].astype(float).reset_index(drop=True), Xr).fit()
        # Full: const + cr(coord, df) (+ covars)
        basis = patsy.dmatrix(
            f"cr(coord, df={df_spline})",
            {"coord": d["coord"].astype(float).values},
            return_type="dataframe",
        ).reset_index(drop=True)
        basis.columns = [f"spline_{i}" for i in range(basis.shape[1])]
        parts = [basis] + ([cov.reset_index(drop=True)] if cov is not None else [])
        Xf = sm.add_constant(pd.concat(parts, axis=1), has_constant="add")
        full = sm.OLS(d["y"].astype(float).reset_index(drop=True), Xf).fit()

        df_num = reduced.df_resid - full.df_resid
        ssr_r, ssr_f = float(reduced.ssr), float(full.ssr)
        if df_num <= 0 or full.df_resid <= 0 or ssr_f <= 0 or ssr_r < ssr_f:
            return {"model_status": "not_fit: spline not nested/degenerate", "n": int(len(d))}
        f_stat = ((ssr_r - ssr_f) / df_num) / (ssr_f / full.df_resid)
        f_p = float(stats.f.sf(f_stat, df_num, full.df_resid))
        return {
            "model_status": "fit",
            "f_stat": float(f_stat),
            "df_num": int(df_num),
            "df_den": int(full.df_resid),
            "p_value": f_p,
            "aic_linear": float(reduced.aic),
            "aic_spline": float(full.aic),
            "aic_linear_minus_spline": float(reduced.aic - full.aic),
            "n": int(len(d)),
            "covariate_adjusted": bool(cov_cols),
        }
    except Exception as exc:  # pragma: no cover
        return {"model_status": f"not_fit: {exc}", "n": int(len(d))}


# #
def _zscore(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    sd = v.std(ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return v * np.nan
    return (v - v.mean()) / sd


# (outcome column, label, predicted supportive sign, family)
_OUTCOMES = [
    ("mean_phospho_effect", "mean_phospho_effect", "negative", "linear"),
    ("most_negative_phospho_effect", "most_negative_phospho_effect", "negative", "linear"),
    ("any_suppressed_p05", "any_suppressed_p05", "positive", "logistic"),
]


def _variant_parent_tables(phospho_prior: pd.DataFrame) -> dict[str, pd.DataFrame]:
    """Build the parent-gene table for each robustness variant (z-scores recomputed within each)."""
    sites = phospho_prior.copy()
    is_anchor = sites["gene_symbol"].isin(ANCHOR_GENES)
    is_single = sites["is_single_site"].astype(bool) if "is_single_site" in sites.columns else pd.Series(True, index=sites.index)
    return {
        "all_parent_genes": build_parent_gene_table(sites),
        "exclude_anchor_genes": build_parent_gene_table(sites[~is_anchor]),
        "single_position_sites_only": build_parent_gene_table(sites[is_single]),
    }


def _row(**kw) -> dict:
    base = {
        "variant": None, "outcome": None, "method": None, "term": None,
        "covariate_adjusted": None, "coord_units": None, "supportive_sign": None,
        "n_parent_genes": np.nan, "n_used": np.nan, "estimate": np.nan,
        "std_error": np.nan, "test_stat": np.nan, "p_value": np.nan,
        "ci_low": np.nan, "ci_high": np.nan, "extra_stat_name": None,
        "extra_stat": np.nan, "model_status": None, "provenance_key": None,
    }
    base.update(kw)
    return base


def run_dct_continuous_gradient(root: Path) -> pd.DataFrame:
    root = Path(root)
    prior_path = root / "dct_prior" / "osd462_phosphosite_dct1_prior.tsv"
    if not prior_path.exists():
        raise FileNotFoundError(
            f"Missing DCT phospho prior: {prior_path}. Run core_analysis (phase 11) first."
        )
    phospho_prior = read_tsv(prior_path)
    variants = _variant_parent_tables(phospho_prior)

    rows: list[dict] = []
    for vname, parent in variants.items():
        n_genes = int(len(parent))
        if not n_genes or COORD_COL not in parent.columns:
            rows.append(_row(variant=vname, method="skipped", n_parent_genes=n_genes,
                             model_status="not_fit: empty variant"))
            continue
        coord_z = _zscore(parent[COORD_COL])
        covars = parent[COVARIATE_COLS] if all(c in parent.columns for c in COVARIATE_COLS) else None
        is_primary = vname == "all_parent_genes"

        for col, label, sign, family in _OUTCOMES:
            if col not in parent.columns:
                continue
            y = parent[col]

            if family == "linear":
                lin = fit_linear_gradient(coord_z, y, covariates=covars)
                prov = "dct_gradient_slope" if (is_primary and col == "mean_phospho_effect") else None
                rows.append(_row(
                    variant=vname, outcome=label, method="ols_linear_per_sd",
                    term=COORD_COL, covariate_adjusted=lin.get("covariate_adjusted"),
                    coord_units="per_sd_dct1_score", supportive_sign=sign,
                    n_parent_genes=n_genes, n_used=lin.get("n", np.nan),
                    estimate=lin.get("slope", np.nan), std_error=lin.get("se", np.nan),
                    test_stat=lin.get("t_stat", np.nan), p_value=lin.get("p_value", np.nan),
                    ci_low=lin.get("ci_low", np.nan), ci_high=lin.get("ci_high", np.nan),
                    extra_stat_name="r_squared", extra_stat=lin.get("r_squared", np.nan),
                    model_status=lin.get("model_status"), provenance_key=prov,
                ))
                # Raw-unit slope for the primary mean-effect contrast (interpretability).
                if is_primary and col == "mean_phospho_effect":
                    lin_raw = fit_linear_gradient(parent[COORD_COL], y, covariates=covars)
                    rows.append(_row(
                        variant=vname, outcome=label, method="ols_linear_raw_unit",
                        term=COORD_COL, covariate_adjusted=lin_raw.get("covariate_adjusted"),
                        coord_units="per_raw_dct1_score_unit", supportive_sign=sign,
                        n_parent_genes=n_genes, n_used=lin_raw.get("n", np.nan),
                        estimate=lin_raw.get("slope", np.nan), std_error=lin_raw.get("se", np.nan),
                        test_stat=lin_raw.get("t_stat", np.nan), p_value=lin_raw.get("p_value", np.nan),
                        ci_low=lin_raw.get("ci_low", np.nan), ci_high=lin_raw.get("ci_high", np.nan),
                        model_status=lin_raw.get("model_status"),
                    ))

                spl = fit_spline_gradient(coord_z, y, covariates=covars)
                rows.append(_row(
                    variant=vname, outcome=label, method="natural_cubic_spline",
                    term="nonlinearity_F_vs_linear", covariate_adjusted=spl.get("covariate_adjusted"),
                    coord_units="per_sd_dct1_score", supportive_sign="n/a (shape test)",
                    n_parent_genes=n_genes, n_used=spl.get("n", np.nan),
                    test_stat=spl.get("f_stat", np.nan), p_value=spl.get("p_value", np.nan),
                    extra_stat_name="aic_linear_minus_spline",
                    extra_stat=spl.get("aic_linear_minus_spline", np.nan),
                    model_status=spl.get("model_status"),
                ))

                sp = spearman_gradient(parent[COORD_COL], y)
                prov_sp = "dct_gradient_spearman" if (is_primary and col == "mean_phospho_effect") else None
                rows.append(_row(
                    variant=vname, outcome=label, method="spearman",
                    term=COORD_COL, covariate_adjusted=False,
                    coord_units="rank", supportive_sign=sign,
                    n_parent_genes=n_genes, n_used=sp.get("n", np.nan),
                    estimate=sp.get("spearman_rho", np.nan), p_value=sp.get("p_value", np.nan),
                    model_status=sp.get("model_status"), provenance_key=prov_sp,
                ))

            elif family == "logistic":
                log = fit_logistic_gradient(coord_z, y, covariates=covars)
                rows.append(_row(
                    variant=vname, outcome=label, method="glm_binomial_per_sd",
                    term=COORD_COL, covariate_adjusted=log.get("covariate_adjusted"),
                    coord_units="per_sd_dct1_score", supportive_sign=sign,
                    n_parent_genes=n_genes, n_used=log.get("n", np.nan),
                    estimate=log.get("log_odds", np.nan), std_error=log.get("se", np.nan),
                    test_stat=log.get("z_stat", np.nan), p_value=log.get("p_value", np.nan),
                    ci_low=log.get("ci_low", np.nan), ci_high=log.get("ci_high", np.nan),
                    extra_stat_name="odds_ratio_per_sd", extra_stat=log.get("odds_ratio", np.nan),
                    model_status=log.get("model_status"),
                ))

    summary = pd.DataFrame(rows)

    # BH across the inferential family (point-estimate tests, adjusted models only).
    infer = summary["method"].isin(["ols_linear_per_sd", "glm_binomial_per_sd", "spearman"])
    summary["q_value"] = np.nan
    if infer.any():
        summary.loc[infer, "q_value"] = bh(summary.loc[infer, "p_value"])

    # Provenance p for dct_gradient_slope -> dct_gradient_p (same primary linear row).
    prim_slope = summary[summary["provenance_key"] == "dct_gradient_slope"]
    if len(prim_slope):
        pr = prim_slope.iloc[0]
        summary = pd.concat([summary, pd.DataFrame([_row(
            variant="all_parent_genes", outcome="mean_phospho_effect",
            method="ols_linear_per_sd", term=COORD_COL, covariate_adjusted=bool(pr["covariate_adjusted"]),
            coord_units="per_sd_dct1_score", supportive_sign="negative",
            n_parent_genes=pr["n_parent_genes"], n_used=pr["n_used"],
            estimate=pr["estimate"], p_value=pr["p_value"], model_status="fit",
            provenance_key="dct_gradient_p", extra_stat_name="duplicate_of_slope_row_pvalue",
        )])], ignore_index=True)

    out_dir = root / "h2_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "h2_dct_continuous_gradient.tsv"
    summary.to_csv(out_path, sep="\t", index=False)

    _write_verdict(root, summary)
    print(f"Repair C continuous DCT gradient written: {out_path} ({len(summary)} rows)")
    return summary


def _supportive(row: pd.Series) -> bool:
    """Is a single inferential row supportive of the DCT1-prior direction at q<=0.10?"""
    if row.get("model_status") != "fit":
        return False
    est, q = row.get("estimate"), row.get("q_value")
    if not np.isfinite(est) or not np.isfinite(q) or q > 0.10:
        return False
    return est < 0 if row.get("supportive_sign") == "negative" else est > 0


def _write_verdict(root: Path, summary: pd.DataFrame) -> None:
    def pick(variant, outcome, method):
        m = summary[(summary["variant"] == variant) & (summary["outcome"] == outcome)
                    & (summary["method"] == method)]
        return m.iloc[0] if len(m) else None

    primary_linear = pick("all_parent_genes", "mean_phospho_effect", "ols_linear_per_sd")
    primary_spear = pick("all_parent_genes", "mean_phospho_effect", "spearman")
    primary_logit = pick("all_parent_genes", "any_suppressed_p05", "glm_binomial_per_sd")
    anchor_linear = pick("exclude_anchor_genes", "mean_phospho_effect", "ols_linear_per_sd")
    anchor_spear = pick("exclude_anchor_genes", "mean_phospho_effect", "spearman")
    spline_primary = pick("all_parent_genes", "mean_phospho_effect", "natural_cubic_spline")

    linear_ok = primary_linear is not None and _supportive(primary_linear)
    spear_ok = primary_spear is not None and _supportive(primary_spear)
    logit_ok = primary_logit is not None and _supportive(primary_logit)
    survives_anchor = ((anchor_linear is not None and _supportive(anchor_linear))
                       or (anchor_spear is not None and _supportive(anchor_spear)))
    nonlinear = bool(spline_primary is not None
                     and np.isfinite(spline_primary.get("p_value", np.nan))
                     and spline_primary["p_value"] <= 0.05)

    primary_pass = linear_ok or spear_ok or logit_ok
    if primary_pass and survives_anchor:
        interp = ("Continuous DCT-subtype-prior gradient of flight phospho suppression is "
                  "supported and survives NCC/SPAK/WNK anchor exclusion. State as an "
                  "exploratory subtype-prior gradient, not a mechanistic DCT1 claim.")
    elif primary_pass:
        interp = ("Continuous DCT-prior gradient is supported in the full parent-gene set but "
                  "weakens when anchor genes are removed; report as anchor-influenced.")
    else:
        interp = ("No continuous DCT-prior gradient detected; the binned H2 enrichment, if any, "
                  "is not corroborated on the continuous coordinate.")
    if nonlinear:
        interp += (" Spline F-test flags non-monotone structure on the coordinate; the linear "
                   "slope is an incomplete summary.")

    verdict = {
        "analysis": "Repair C -- continuous DCT1<->DCT2 contrast (parent-gene gradient)",
        "coordinate": COORD_COL,
        "claim_discipline": "GSE228367 DCT-subtype RNA prior; exploratory, not spaceflight evidence.",
        "primary_linear_slope_per_sd": None if primary_linear is None else float(primary_linear["estimate"]),
        "primary_linear_q": None if primary_linear is None else float(primary_linear.get("q_value", np.nan)),
        "primary_spearman_rho": None if primary_spear is None else float(primary_spear["estimate"]),
        "primary_logistic_log_odds_per_sd": None if primary_logit is None else float(primary_logit["estimate"]),
        "linear_supportive": bool(linear_ok),
        "spearman_supportive": bool(spear_ok),
        "logistic_supportive": bool(logit_ok),
        "survives_anchor_exclusion": bool(survives_anchor),
        "nonlinearity_flagged_by_spline": nonlinear,
        "interpretation": interp,
    }
    (root / "h2_enrichment" / "h2_dct_continuous_gradient_verdict.json").write_text(
        json.dumps(verdict, indent=2)
    )


def main():
    parser = argparse.ArgumentParser(description="Repair C: continuous DCT1<->DCT2 phospho gradient.")
    parser.add_argument("--run-root", default=str(RUN_ROOT))
    args = parser.parse_args()
    run_dct_continuous_gradient(Path(args.run_root))


if __name__ == "__main__":
    main()
