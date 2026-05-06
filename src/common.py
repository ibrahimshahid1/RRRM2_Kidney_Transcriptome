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
import re
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


def by_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Yekutieli FDR correction for arbitrarily dependent tests."""
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return np.asarray([], dtype=float)
    harmonic = np.sum(1.0 / np.arange(1, n + 1))
    return bh_fdr(np.clip(p * harmonic, 0, 1))


def load_id_map(path: str | Path) -> pd.DataFrame:
    """Load the Ensembl-to-symbol map with canonical column names."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"ID map not found: {path}")
    df = pd.read_csv(path, sep="\t", comment="#")
    renames: dict[str, str] = {}
    for col in df.columns:
        cl = col.lower().strip()
        if "ensembl" in cl and "gene" in cl:
            renames[col] = "ensembl_gene_id"
        elif "symbol" in cl or cl == "mgi":
            renames[col] = "mgi_symbol"
    df = df.rename(columns=renames)
    required = {"ensembl_gene_id", "mgi_symbol"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing required columns: {sorted(missing)}")
    df["ensembl_gene_id"] = df["ensembl_gene_id"].astype(str).str.strip()
    df["mgi_symbol"] = df["mgi_symbol"].astype(str).str.strip()
    df = df[(df["ensembl_gene_id"] != "") & (df["mgi_symbol"] != "") & (df["mgi_symbol"] != "nan")]
    return df


def id_map_lookup(path: str | Path) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Return Ensembl→symbol and lowercase symbol→Ensembl lookup maps."""
    df = load_id_map(path)
    ens_to_symbol = dict(zip(df["ensembl_gene_id"], df["mgi_symbol"]))
    symbol_to_ens: dict[str, set[str]] = {}
    for _, row in df.iterrows():
        symbol_to_ens.setdefault(str(row["mgi_symbol"]).lower(), set()).add(str(row["ensembl_gene_id"]))
    return ens_to_symbol, symbol_to_ens


def flatten_yaml_gene_groups(config: dict) -> list[dict[str, str]]:
    """Flatten config['anchors']-style groups into records with group/symbol."""
    records: list[dict[str, str]] = []
    groups = config.get("anchors", {}) if isinstance(config, dict) else {}
    for group_name, group_data in groups.items():
        if isinstance(group_data, dict):
            genes = group_data.get("genes", [])
            role = group_data.get("role", group_name)
        else:
            genes = group_data
            role = group_name
        for gene in genes or []:
            gene = str(gene).strip()
            if gene:
                records.append({"group": str(group_name), "role": str(role), "symbol": gene})
    return records


_ENSEMBL_RE = re.compile(r"^ENS[A-Z]*G\d+")


def resolve_configured_genes(
    symbols_or_ids: list[str],
    id_map_path: str | Path,
    panel_genes: set[str] | None = None,
) -> pd.DataFrame:
    """Resolve configured gene symbols/IDs to Ensembl IDs and panel status.

    No observed rewiring statistic is consulted here.  Ambiguous symbols are kept
    as separate resolved rows; callers can decide whether ambiguity is acceptable.
    """
    ens_to_symbol, symbol_to_ens = id_map_lookup(id_map_path)
    panel = set(panel_genes) if panel_genes is not None else None
    rows: list[dict[str, object]] = []

    for raw in symbols_or_ids:
        query = str(raw).strip()
        if not query:
            continue
        resolved: set[str]
        if _ENSEMBL_RE.match(query):
            resolved = {query}
            symbol = ens_to_symbol.get(query, query)
            query_type = "ensembl"
        else:
            resolved = set(symbol_to_ens.get(query.lower(), set()))
            symbol = query
            query_type = "symbol"

        if not resolved:
            rows.append({
                "query": query,
                "query_type": query_type,
                "symbol": symbol,
                "ensembl_gene_id": "",
                "status": "unmapped",
                "in_panel": False if panel is not None else np.nan,
            })
            continue

        for eid in sorted(resolved):
            in_panel = eid in panel if panel is not None else np.nan
            rows.append({
                "query": query,
                "query_type": query_type,
                "symbol": ens_to_symbol.get(eid, symbol),
                "ensembl_gene_id": eid,
                "status": "mapped",
                "in_panel": in_panel,
            })
    return pd.DataFrame(rows)


def load_anchor_config(anchor_config: str | Path) -> tuple[dict, list[dict[str, str]]]:
    """Load anchor_genes.yaml and return raw config plus flattened records."""
    try:
        import yaml
    except Exception as exc:  # pragma: no cover - dependency guard
        raise ImportError("PyYAML is required to load configured anchors") from exc
    anchor_config = Path(anchor_config)
    if not anchor_config.exists():
        raise FileNotFoundError(f"Anchor config not found: {anchor_config}")
    cfg = yaml.safe_load(anchor_config.read_text()) or {}
    records = flatten_yaml_gene_groups(cfg)
    if not records:
        raise ValueError(f"No anchor genes found in {anchor_config}")
    return cfg, records
