"""Shared matched-null helpers for v11 layer-specificity modules."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from src.multiomics.osd462_anchor import (
    MatchedNullResult,
    assign_match_strata,
    matched_null_test,
)

EPS = 1e-12


def finite_mask(*arrays: object) -> np.ndarray:
    """Rows where all supplied arrays are finite numeric values."""
    if not arrays:
        return np.asarray([], dtype=bool)
    masks = [np.isfinite(pd.to_numeric(pd.Series(a), errors="coerce").to_numpy(dtype=float)) for a in arrays]
    out = masks[0].copy()
    for mask in masks[1:]:
        out &= mask
    return out


def ols_slope(x: object, y: object, min_n: int = 3) -> float:
    """Intercept-inclusive least-squares slope for ``y ~ x`` with NaN handling."""
    xx = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    yy = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(xx) & np.isfinite(yy)
    if int(keep.sum()) < min_n:
        return float("nan")
    xx = xx[keep]
    yy = yy[keep]
    xc = xx - xx.mean()
    denom = float(np.dot(xc, xc))
    if denom <= EPS:
        return float("nan")
    return float(np.dot(xc, yy - yy.mean()) / denom)


def sign_agreement_rate(x: object, y: object, min_n: int = 1) -> float:
    """Fraction of finite, non-zero pairs with matching signs."""
    xx = pd.to_numeric(pd.Series(x), errors="coerce").to_numpy(dtype=float)
    yy = pd.to_numeric(pd.Series(y), errors="coerce").to_numpy(dtype=float)
    keep = np.isfinite(xx) & np.isfinite(yy) & (xx != 0) & (yy != 0)
    if int(keep.sum()) < min_n:
        return float("nan")
    return float(np.mean(np.sign(xx[keep]) == np.sign(yy[keep])))


def prepare_matched_pool(
    df: pd.DataFrame,
    required_cols: list[str],
    abundance_col: str = "abundance_log2",
    peptide_col: str = "n_peptides",
    peptide_filter: int = 2,
    n_abundance_bins: int = 5,
    n_peptide_bins: int = 4,
) -> pd.DataFrame:
    """Filter a dataframe to a finite analysis pool and assign fresh strata."""
    cols = list(dict.fromkeys(required_cols + [abundance_col, peptide_col]))
    missing = [col for col in cols if col not in df.columns]
    if missing:
        raise KeyError(f"matched-null pool missing columns: {missing}")

    pool = df.copy()
    for col in cols:
        pool[col] = pd.to_numeric(pool[col], errors="coerce")
    keep = np.ones(len(pool), dtype=bool)
    for col in cols:
        keep &= np.isfinite(pool[col].to_numpy(dtype=float))
    keep &= pool[peptide_col].to_numpy(dtype=float) >= float(peptide_filter)
    pool = pool.loc[keep].reset_index(drop=True)
    if pool.empty:
        pool["match_stratum"] = pd.Series(dtype=object)
        return pool
    pool["match_stratum"] = assign_match_strata(
        pool,
        abundance_col=abundance_col,
        peptide_col=peptide_col,
        n_abundance_bins=n_abundance_bins,
        n_peptide_bins=n_peptide_bins,
    ).values
    return pool


def run_matched_null(
    pool: pd.DataFrame,
    target_mask: np.ndarray,
    stat_fn: Callable[[pd.DataFrame], float],
    statistic_name: str,
    n_null: int,
    rng: np.random.Generator,
) -> MatchedNullResult:
    """Thin v11 wrapper around the OSD-462 anchor matched-null implementation."""
    return matched_null_test(
        pool=pool,
        target_mask=np.asarray(target_mask, dtype=bool),
        stat_fn=stat_fn,
        strata=pool["match_stratum"],
        statistic_name=statistic_name,
        n_null=n_null,
        rng=rng,
    )


__all__ = [
    "MatchedNullResult",
    "finite_mask",
    "ols_slope",
    "prepare_matched_pool",
    "run_matched_null",
    "sign_agreement_rate",
]
