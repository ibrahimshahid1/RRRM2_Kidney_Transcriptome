"""Phenotype-anchoring layer: animal-matched RNA-state vs NCC-activity.

OSD-462 / RR-10 ran kidney RNA-seq and phosphoproteomics on the *same animals*
(channel ``FL-01`` <-> RNA sample ``FLT_F1`` etc.). This module builds a
per-animal NCC/SPAK regulatory-phosphorylation activity score and a per-animal
DCT/NCC-low RNA state score, matches them by animal, and compares them at three
levels of stringency:

1. group level (Space Flight vs Ground Control) -- the primary, robust claim;
2. all-sample correlation -- descriptive, inflated by the group means;
3. condition-adjusted correlation -- the real test of an animal-level link,
   on within-condition-centered residuals; underpowered at n~10/condition.

Controls: a non-regulatory NCC phosphosite score (should not track the RNA
state) and an RNA score with the Slc12a3 transcript removed (the link must not
be one gene covarying with its own phospho).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, asdict

import numpy as np
import pandas as pd
from scipy import stats

EPS = 1e-12

# Pre-specified NCC/SPAK regulatory (SPAK/OSR1-target activating) phosphosites
NCC_REGULATORY_SITES = (
    ("Slc12a3", "53"), ("Slc12a3", "58"), ("Slc12a3", "65"), ("Slc12a3", "68"),
    ("Stk39", "382"), ("Stk39", "383"),
)
NCC_REGULATORY_SITES_SENS = NCC_REGULATORY_SITES + (("Slc12a3", "89"),)
# Non-regulatory NCC sites -- negative control (should not track the RNA state)
NCC_NONREGULATORY_SITES = (
    ("Slc12a3", "96"), ("Slc12a3", "120"), ("Slc12a3", "122"), ("Slc12a3", "124"),
)


@dataclass(frozen=True)
class ComparisonResult:
    n_animals: int
    group_rna_flt_minus_gc: float
    group_phospho_flt_minus_gc: float
    spearman_all: float
    spearman_all_p: float
    spearman_condition_adjusted: float
    spearman_condition_adjusted_p: float
    interpretation: str


def channel_to_animal(label: str) -> tuple[str, str] | None:
    """Map a phospho channel label (``FL-01``/``GC-12``/``BL-09``) to
    (condition, animal-number). Returns None for unrecognized labels."""
    m = re.fullmatch(r"\s*(BL|FL|GC)-0*(\d+)\s*", str(label))
    if not m:
        return None
    cond = {"BL": "Basal", "FL": "Space Flight", "GC": "Ground Control"}[m.group(1)]
    return cond, m.group(2)


def rna_sample_animal(sample_name: str) -> tuple[str, str] | None:
    """Map an OSD-462 RNA sample name to (condition, animal-number)."""
    m = re.search(r"_(BSL|FLT|GC)_([BFG])(\d+)", str(sample_name))
    if not m:
        return None
    cond = {"BSL": "Basal", "FLT": "Space Flight", "GC": "Ground Control"}[m.group(1)]
    return cond, m.group(3)


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    """Z-score each row (feature) across samples; drop zero-variance rows."""
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=1)
    keep = sd > EPS
    return mat.loc[keep].sub(mu[keep], axis=0).div(sd[keep], axis=0)


def per_sample_score(values: pd.DataFrame, feature_keys: list[str]) -> pd.Series:
    """Mean z-scored value across the requested features, per sample (column).

    ``values``: features x samples. Missing features are skipped (logged by
    caller via the returned coverage in n)."""
    present = [k for k in feature_keys if k in values.index]
    if not present:
        return pd.Series(dtype=float)
    z = zscore_rows(values.loc[present])
    return z.mean(axis=0)


def condition_adjusted_residuals(score: pd.Series, condition: pd.Series) -> pd.Series:
    """Within-condition-centered residuals (removes the group-mean structure)."""
    df = pd.DataFrame({"s": score, "c": condition}).dropna()
    return df["s"] - df.groupby("c")["s"].transform("mean")


def compare_scores(
    rna_score: pd.Series,
    phospho_score: pd.Series,
    condition: pd.Series,
    *,
    flight_label: str = "Space Flight",
    ground_label: str = "Ground Control",
) -> ComparisonResult:
    """Three-level comparison of two animal-matched scores."""
    df = pd.DataFrame({"rna": rna_score, "phospho": phospho_score, "cond": condition}).dropna()
    flt = df[df["cond"].eq(flight_label)]
    gc = df[df["cond"].eq(ground_label)]
    g_rna = float(flt["rna"].mean() - gc["rna"].mean()) if len(flt) and len(gc) else np.nan
    g_pho = float(flt["phospho"].mean() - gc["phospho"].mean()) if len(flt) and len(gc) else np.nan

    paired = df[df["cond"].isin([flight_label, ground_label])]
    if len(paired) >= 4:
        rho_all, p_all = stats.spearmanr(paired["rna"], paired["phospho"])
        r_rna = condition_adjusted_residuals(paired["rna"], paired["cond"])
        r_pho = condition_adjusted_residuals(paired["phospho"], paired["cond"])
        if len(r_rna) >= 4 and r_rna.std() > EPS and r_pho.std() > EPS:
            rho_adj, p_adj = stats.spearmanr(r_rna, r_pho)
        else:
            rho_adj, p_adj = np.nan, np.nan
    else:
        rho_all = p_all = rho_adj = p_adj = np.nan

    same_group_dir = np.isfinite(g_rna) and np.isfinite(g_pho) and np.sign(g_rna) == np.sign(g_pho)
    if same_group_dir and np.isfinite(p_adj) and p_adj < 0.05:
        interp = "animal_level_link_supported"
    elif same_group_dir:
        interp = "group_level_concordant_animal_level_underpowered"
    else:
        interp = "no_consistent_link"

    return ComparisonResult(
        n_animals=int(len(paired)),
        group_rna_flt_minus_gc=g_rna,
        group_phospho_flt_minus_gc=g_pho,
        spearman_all=float(rho_all) if np.isfinite(rho_all) else np.nan,
        spearman_all_p=float(p_all) if np.isfinite(p_all) else np.nan,
        spearman_condition_adjusted=float(rho_adj) if np.isfinite(rho_adj) else np.nan,
        spearman_condition_adjusted_p=float(p_adj) if np.isfinite(p_adj) else np.nan,
        interpretation=interp,
    )


def result_to_dict(r: ComparisonResult) -> dict:
    return asdict(r)
