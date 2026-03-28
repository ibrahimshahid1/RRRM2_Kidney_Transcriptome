# src/common.py
"""
Shared utilities for the RRRM-2 pipeline.

Centralizes functions that were previously copy-pasted across multiple modules:
  - REPO_ROOT: repository root path
  - find_sample_col: detect sample identifier column in metadata
  - normalize_labels: canonical Age/Arm/EnvGroup label normalization
  - bh_fdr: Benjamini-Hochberg FDR correction
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


# Repository root — single source of truth
REPO_ROOT = Path(__file__).resolve().parents[1]


def find_sample_col(meta: pd.DataFrame) -> str:
    """Find the sample identifier column in metadata.

    Checks for common column names used across OSD-771 metadata files.
    Falls back to the first column if none match.
    """
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            return col
    return meta.columns[0]


def normalize_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Normalize Age, Arm, and EnvGroup labels to canonical forms.

    Canonical forms:
        Age:      YNG, OLD
        Arm:      ISS-T, LAR
        EnvGroup: FLT, GC, VIV, BSL

    This is the single authoritative normalization used across all pipeline
    phases.  Previous versions in full_regression.py applied .str.upper()
    before replacement, which caused HGC to remain unmapped — that bug is
    fixed here.
    """
    meta = meta.copy()

    if "Age" in meta.columns:
        meta["Age"] = meta["Age"].astype(str).replace({
            "Young": "YNG", "Yng": "YNG", "young": "YNG",
            "Old": "OLD", "old": "OLD",
            "YOUNG": "YNG", "OLD": "OLD",
        })

    if "Arm" in meta.columns:
        meta["Arm"] = meta["Arm"].astype(str).replace({
            "ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T", "ISS T": "ISS-T",
            "LAR_T": "LAR", "LAR-T": "LAR", "LAR T": "LAR",
        })

    if "EnvGroup" in meta.columns:
        meta["EnvGroup"] = meta["EnvGroup"].astype(str).replace({
            "HGC": "GC", "VGC": "VIV",
            "HGC/GC": "GC", "VIV/VGC": "VIV",
            "FLIGHT": "FLT", "GROUND CONTROL": "GC",
        })

    return meta


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p : array-like of raw p-values

    Returns
    -------
    q : ndarray of adjusted p-values (same shape as input), clipped to [0, 1]
    """
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out
