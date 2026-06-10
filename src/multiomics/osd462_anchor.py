"""OSD-462 / RR-10 multi-omics anchor.

Dataset-agnostic estimation and inference helpers used by the
``scripts/osd462/*`` layer scripts to test whether the RRRM-2 RNA-level
matrix-high / DCT-low remodeling signal recurs at the protein and
phosphoprotein level in the independent OSD-462 spaceflight kidney cohort.

Design notes
------------
* Cross-study, cross-strain, cross-modality: no sample-level pooling.  The
  unit of comparison is the *direction of the flight effect* (Space Flight
  minus Ground Control), exactly the contrast-vector logic already used for
  OSD-513 in ``src/networks/cross_osdr_projection.py``.
* TMT 2-plex batch structure ("Samp1-5", "Samp6-10").  Plex is handled as a
  batch factor everywhere: flight effects are estimated *within* each plex and
  then averaged, which removes any per-plex (and, after channel centering,
  per-channel) loading constant.
* Primary inference is an abundance/peptide-matched random-gene-set null,
  which is robust to the generically weak RNA<->protein correlation and to
  gene-set size.

The TMT workbook layout (both proteomics and phosphoproteomics) is::

    row 1 : group banners ("... scaled ..." marks the scaled blocks)
    row 2 : per-column sample labels  (BL-01, FL-03, GC-05, ...)
    row 3 : machine column headers     (Samp1-5~rq_129n_sn scaled, ...)
    row 4+: one protein / phosphosite per row
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

PLEX1 = "Samp1-5"
PLEX2 = "Samp6-10"
CONDITIONS = ("BL", "FL", "GC")
EPS = 1e-12


# ─────────────────────────────────────────────────────────────────────────────
# TMT workbook parsing
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class TmtTable:
    """Parsed TMT workbook sheet (proteomics or phosphoproteomics)."""

    meta: pd.DataFrame          # row metadata (gene symbol, ids, peptide counts, ...)
    scaled: pd.DataFrame        # rows x value-columns of scaled S/N
    channels: pd.DataFrame      # one row per value-column: column, plex, condition, sample
    sheet: str
    n_rows: int = field(init=False)

    def __post_init__(self) -> None:
        self.n_rows = int(self.scaled.shape[0])

    def condition_columns(self, plex: str, condition: str) -> list[str]:
        sel = self.channels[(self.channels["plex"] == plex)
                            & (self.channels["condition"] == condition)]
        return list(sel["column"])


def _classify_header(header: str, sample_label: str) -> dict | None:
    """Map a row-3 machine header + row-2 sample label to channel metadata.

    Returns ``None`` for non-quantitative columns (metadata columns).
    """
    if not isinstance(header, str) or "~" not in header:
        return None
    low = header.lower()
    if "_sn" not in low and "sn " not in low and "sn_" not in low:
        return None
    plex = header.split("~", 1)[0].strip()
    if plex not in (PLEX1, PLEX2):
        return None
    value_type = "scaled" if "scaled" in low else "sum"
    cond = None
    sample = None
    if isinstance(sample_label, str) and "-" in sample_label:
        prefix = sample_label.split("-", 1)[0].strip().upper()
        if prefix in CONDITIONS:
            cond = prefix
            sample = sample_label.strip()
    if cond is None:
        return None
    return {"plex": plex, "value_type": value_type, "condition": cond, "sample": sample}


def parse_tmt_sheet(
    path: str | Path,
    sheet: str,
    gene_col: str,
    peptide_cols: dict[str, str],
    id_col: str | None = None,
    extra_meta_cols: Sequence[str] | None = None,
    header_row: int = 3,
) -> TmtTable:
    """Parse one TMT workbook sheet into a :class:`TmtTable`.

    Parameters
    ----------
    path, sheet
        Workbook path and sheet name.
    gene_col
        Name of the gene-symbol column in the row-``header_row`` header.
    peptide_cols
        Mapping ``{PLEX1: <peptide col name>, PLEX2: <peptide col name>}``.
    id_col
        Optional protein-id column to retain.
    extra_meta_cols
        Additional metadata column names to retain (e.g. site position).
    header_row
        1-based index of the machine-header row (3 for these workbooks).
    """
    import openpyxl

    path = Path(path)
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[sheet]
    rows_iter = ws.iter_rows(values_only=True)
    pre: list[tuple] = []
    for _ in range(header_row):
        pre.append(next(rows_iter))
    sample_labels = pre[header_row - 2]   # row 2
    headers = pre[header_row - 1]         # row 3 (machine headers)

    # Locate metadata columns by header name.
    name_to_idx: dict[str, int] = {}
    for j, h in enumerate(headers):
        if isinstance(h, str):
            name_to_idx[h.strip()] = j

    def need(col: str) -> int:
        if col not in name_to_idx:
            raise KeyError(f"column {col!r} not found in sheet {sheet!r}")
        return name_to_idx[col]

    gene_idx = need(gene_col)
    id_idx = name_to_idx.get(id_col) if id_col else None
    pep_idx = {p: need(c) for p, c in peptide_cols.items()}
    extra_idx = {c: name_to_idx[c] for c in (extra_meta_cols or []) if c in name_to_idx}

    # Build channel metadata for scaled value columns only.
    channels: list[dict] = []
    for j, h in enumerate(headers):
        info = _classify_header(h, sample_labels[j] if j < len(sample_labels) else None)
        if info is None or info["value_type"] != "scaled":
            continue
        col_name = f"{info['plex']}|{info['condition']}|{info['sample']}"
        channels.append({"column": col_name, "idx": j, **info})
    channel_df = pd.DataFrame(channels)
    if channel_df.empty:
        raise ValueError(f"no scaled TMT channels detected in sheet {sheet!r}")

    meta_records: list[dict] = []
    value_records: list[list[float]] = []
    value_idx = list(channel_df["idx"])
    value_names = list(channel_df["column"])
    for row in rows_iter:
        if row is None:
            continue
        gene = row[gene_idx]
        if gene is None or (isinstance(gene, str) and gene.strip() == ""):
            # skip rows without a gene symbol
            continue
        rec = {"gene_symbol": str(gene).strip()}
        if id_idx is not None:
            rec["protein_id"] = row[id_idx]
        for p, idx in pep_idx.items():
            rec[f"n_pep_{p}"] = _to_float(row[idx])
        for c, idx in extra_idx.items():
            rec[c] = row[idx]
        meta_records.append(rec)
        value_records.append([_to_float(row[j]) for j in value_idx])

    meta = pd.DataFrame(meta_records).reset_index(drop=True)
    scaled = pd.DataFrame(value_records, columns=value_names).reset_index(drop=True)
    channel_df = channel_df.drop(columns=["idx"]).reset_index(drop=True)
    wb.close()
    return TmtTable(meta=meta, scaled=scaled, channels=channel_df, sheet=sheet)


def _to_float(v) -> float:
    try:
        if v is None:
            return np.nan
        return float(v)
    except (TypeError, ValueError):
        return np.nan


# ─────────────────────────────────────────────────────────────────────────────
# Flight-effect estimation (within-plex FL - GC, plex-averaged)
# ─────────────────────────────────────────────────────────────────────────────

def _log2_channel_matrix(
    table: TmtTable, plex: str, condition: str, center: pd.Series | None
) -> np.ndarray:
    cols = table.condition_columns(plex, condition)
    if not cols:
        return np.empty((table.n_rows, 0))
    block = table.scaled[cols].to_numpy(dtype=float)
    block = np.where(block > 0, block, np.nan)
    log2 = np.log2(block)
    if center is not None:
        log2 = log2 - center[cols].to_numpy(dtype=float)[None, :]
    return log2


def _plex_channel_centers(table: TmtTable, plex: str) -> pd.Series:
    """Per-channel (sample-loading) median of log2 scaled values within a plex."""
    cols = list(table.channels[table.channels["plex"] == plex]["column"])
    block = table.scaled[cols].to_numpy(dtype=float)
    block = np.where(block > 0, block, np.nan)
    with np.errstate(invalid="ignore"):
        med = np.nanmedian(np.log2(block), axis=0)
    return pd.Series(med, index=cols)


def compute_flight_effect(
    table: TmtTable,
    min_channels_per_condition: int = 2,
    channel_center: bool = True,
) -> pd.DataFrame:
    """Per-row Flight - Ground flight effect, estimated within each plex.

    For each plex the effect is ``mean(log2 FL) - mean(log2 GC)`` over the
    available channels; the two plex estimates are then averaged.  Optional
    per-channel median centering (``channel_center``) performs standard TMT
    sample-loading normalization within each plex; because the differencing is
    within plex, both the per-plex normalization constant and (after centering)
    per-channel loading cancel.

    Returns one row per input row with the flight effect, per-plex effects,
    plex coverage, peptide counts, and a mean-abundance column used for
    matched-null binning.
    """
    centers = {p: (_plex_channel_centers(table, p) if channel_center else None)
               for p in (PLEX1, PLEX2)}

    def plex_effect(plex: str) -> tuple[np.ndarray, np.ndarray]:
        fl = _log2_channel_matrix(table, plex, "FL", centers[plex])
        gc = _log2_channel_matrix(table, plex, "GC", centers[plex])
        n_fl = np.sum(np.isfinite(fl), axis=1)
        n_gc = np.sum(np.isfinite(gc), axis=1)
        ok = (n_fl >= min_channels_per_condition) & (n_gc >= min_channels_per_condition)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            eff = np.nanmean(fl, axis=1) - np.nanmean(gc, axis=1)
        eff = np.where(ok, eff, np.nan)
        return eff, ok

    eff1, ok1 = plex_effect(PLEX1)
    eff2, ok2 = plex_effect(PLEX2)
    stacked = np.vstack([eff1, eff2])
    with np.errstate(invalid="ignore"):
        flight_effect = np.nanmean(stacked, axis=0)
    plex_coverage = (np.isfinite(eff1).astype(int) + np.isfinite(eff2).astype(int))

    # Mean abundance across all present scaled channels (for matched-null bins).
    all_cols = list(table.channels["column"])
    allmat = table.scaled[all_cols].to_numpy(dtype=float)
    allmat = np.where(allmat > 0, allmat, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        abundance = np.nanmean(np.log2(allmat), axis=1)
    n_channels_used = np.sum(np.isfinite(allmat), axis=1)

    out = table.meta.copy()
    out["effect_plex1"] = eff1
    out["effect_plex2"] = eff2
    out["flight_effect"] = flight_effect
    out["plex_coverage"] = plex_coverage
    out["abundance_log2"] = abundance
    out["n_channels_used"] = n_channels_used
    pep1 = out.get(f"n_pep_{PLEX1}", pd.Series(np.nan, index=out.index)).fillna(0)
    pep2 = out.get(f"n_pep_{PLEX2}", pd.Series(np.nan, index=out.index)).fillna(0)
    out["n_peptides"] = (pep1 + pep2).to_numpy()
    return out


def compute_site_flight_effect_lm(
    table: TmtTable,
    min_per_condition: int = 3,
    channel_center: bool = True,
) -> pd.DataFrame:
    """Per-row FL - GC effect with CI from a plex-adjusted linear model.

    For each row, log2 scaled FL/GC channels (both plexes) are fit with
    ``y ~ flight + plex``; the flight coefficient is the plex-corrected
    FL - GC effect.  By default, per-channel median centering performs the same
    within-plex TMT loading normalization used for protein effects.  Returns
    effect, SE, 95% CI, two-sided p, and channel counts.  Suitable for the
    small phosphosite tables where a per-site interval is wanted.
    """
    from scipy.stats import t as t_dist

    centers = {p: (_plex_channel_centers(table, p) if channel_center else None)
               for p in (PLEX1, PLEX2)}
    # Build the long design once: rows are (plex, condition) channel columns.
    obs_cols: list[tuple[str, int, int]] = []  # (column, flight, plex2)
    for _, ch in table.channels.iterrows():
        if ch["condition"] not in ("FL", "GC"):
            continue
        flight = 1 if ch["condition"] == "FL" else 0
        plex2 = 1 if ch["plex"] == PLEX2 else 0
        obs_cols.append((ch["column"], flight, plex2))
    cols = [c for c, _, _ in obs_cols]
    flight_vec = np.array([f for _, f, _ in obs_cols], dtype=float)
    plex2_vec = np.array([p for _, _, p in obs_cols], dtype=float)

    block = table.scaled[cols].to_numpy(dtype=float)
    block = np.where(block > 0, block, np.nan)
    log2 = np.log2(block)
    if channel_center:
        for p in (PLEX1, PLEX2):
            pcols = [i for i, (c, _, _) in enumerate(obs_cols)
                     if c in set(table.channels[table.channels["plex"] == p]["column"])]
            if pcols and centers[p] is not None:
                log2[:, pcols] = log2[:, pcols] - centers[p][[obs_cols[i][0] for i in pcols]].to_numpy()[None, :]

    n_rows = log2.shape[0]
    out = {k: np.full(n_rows, np.nan) for k in
           ("effect", "se", "ci_low", "ci_high", "p_value", "n_fl", "n_gc")}
    crit_cache: dict[int, float] = {}
    for i in range(n_rows):
        y = log2[i]
        m = np.isfinite(y)
        n_fl = int(((flight_vec == 1) & m).sum())
        n_gc = int(((flight_vec == 0) & m).sum())
        out["n_fl"][i] = n_fl
        out["n_gc"][i] = n_gc
        if n_fl < min_per_condition or n_gc < min_per_condition:
            continue
        yy = y[m]
        # design: intercept, flight, plex2 (drop plex2 column if only one plex present)
        cols_d = [np.ones(m.sum()), flight_vec[m]]
        if len(np.unique(plex2_vec[m])) > 1:
            cols_d.append(plex2_vec[m])
        X = np.column_stack(cols_d)
        df_resid = X.shape[0] - X.shape[1]
        if df_resid <= 0:
            continue
        beta, _, _, _ = np.linalg.lstsq(X, yy, rcond=None)
        resid = yy - X @ beta
        sigma2 = float(resid @ resid) / df_resid
        try:
            xtx_inv = np.linalg.inv(X.T @ X)
        except np.linalg.LinAlgError:
            continue
        se = float(np.sqrt(sigma2 * xtx_inv[1, 1]))
        eff = float(beta[1])
        out["effect"][i] = eff
        out["se"][i] = se
        if se > 0:
            if df_resid not in crit_cache:
                crit_cache[df_resid] = float(t_dist.ppf(0.975, df_resid))
            tc = crit_cache[df_resid]
            out["ci_low"][i] = eff - tc * se
            out["ci_high"][i] = eff + tc * se
            out["p_value"][i] = float(2 * t_dist.sf(abs(eff / se), df_resid))

    res = table.meta.copy()
    for k, v in out.items():
        res[f"phospho_{k}" if k in ("effect", "se", "ci_low", "ci_high", "p_value") else k] = v
    return res


def collapse_to_gene(
    effects: pd.DataFrame,
    effect_col: str = "flight_effect",
    require_both_plex: bool = True,
) -> pd.DataFrame:
    """Collapse multiple protein rows per gene symbol into one gene-level row.

    Many-to-one collisions (isoforms / shared symbols) are resolved by a
    peptide-weighted mean of the flight effect; peptide counts are summed and
    the number of collapsed rows is logged in ``n_protein_rows``.
    """
    df = effects.copy()
    if require_both_plex:
        df = df[df["plex_coverage"] == 2]
    df = df[np.isfinite(df[effect_col])]
    df = df[df["gene_symbol"].notna() & (df["gene_symbol"] != "")]

    def agg(group: pd.DataFrame) -> pd.Series:
        w = group["n_peptides"].to_numpy(dtype=float)
        if not np.isfinite(w).any() or np.nansum(w) <= 0:
            w = np.ones(len(group))
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() <= 0:
            w = np.ones(len(group))
        eff = np.average(group[effect_col].to_numpy(dtype=float), weights=w)
        return pd.Series({
            effect_col: eff,
            "n_peptides": float(np.nansum(group["n_peptides"].to_numpy(dtype=float))),
            "abundance_log2": float(np.nanmean(group["abundance_log2"].to_numpy(dtype=float))),
            "plex_coverage": int(group["plex_coverage"].max()),
            "n_protein_rows": int(len(group)),
        })

    collapsed = df.groupby("gene_symbol", sort=False).apply(agg).reset_index()
    return collapsed


# ─────────────────────────────────────────────────────────────────────────────
# Gene-set statistics
# ─────────────────────────────────────────────────────────────────────────────

def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rho with NaN handling and tie-aware ranking."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return float("nan")
    from scipy.stats import rankdata
    ra = rankdata(a[mask])
    rb = rankdata(b[mask])
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt(np.sum(ra ** 2) * np.sum(rb ** 2))
    if denom <= EPS:
        return float("nan")
    return float(np.sum(ra * rb) / denom)


def sign_agreement(rna: np.ndarray, prot: np.ndarray) -> float:
    rna = np.asarray(rna, dtype=float)
    prot = np.asarray(prot, dtype=float)
    mask = np.isfinite(rna) & np.isfinite(prot) & (rna != 0) & (prot != 0)
    if mask.sum() == 0:
        return float("nan")
    return float(np.mean(np.sign(rna[mask]) == np.sign(prot[mask])))


# ─────────────────────────────────────────────────────────────────────────────
# Abundance / peptide-matched random-gene-set null
# ─────────────────────────────────────────────────────────────────────────────

def assign_match_strata(
    pool: pd.DataFrame,
    abundance_col: str = "abundance_log2",
    peptide_col: str = "n_peptides",
    n_abundance_bins: int = 5,
    n_peptide_bins: int = 4,
) -> pd.Series:
    """Joint abundance x peptide-count stratum label for each pool gene."""
    def safe_qcut(series: pd.Series, q: int, prefix: str) -> pd.Series:
        s = series.astype(float)
        try:
            codes = pd.qcut(s.rank(method="first"), q=q, labels=False, duplicates="drop")
        except ValueError:
            codes = pd.Series(0, index=s.index)
        return prefix + codes.astype("Int64").astype(str)

    ab = safe_qcut(pool[abundance_col], n_abundance_bins, "a")
    pep = safe_qcut(pool[peptide_col], n_peptide_bins, "p")
    return (ab.astype(str) + "_" + pep.astype(str)).reset_index(drop=True)


@dataclass
class MatchedNullResult:
    statistic: str
    observed: float
    null_median: float
    null_ci_low: float
    null_ci_high: float
    p_greater: float
    p_two_sided: float
    n_null_valid: int
    n_target: int


def matched_null_test(
    pool: pd.DataFrame,
    target_mask: np.ndarray,
    stat_fn: Callable[[pd.DataFrame], float],
    strata: pd.Series,
    statistic_name: str,
    n_null: int = 10000,
    rng: np.random.Generator | None = None,
    alpha: float = 0.05,
) -> MatchedNullResult:
    """Stratum-matched random-gene-set null for a gene-set statistic.

    ``pool`` is the background of all genes (with effects + match columns);
    ``target_mask`` selects the gene set; ``stat_fn`` maps a sub-frame of
    ``pool`` to a scalar; ``strata`` is the matched-sampling stratum per pool
    row.  Each null draw samples, within each stratum, the same number of genes
    the target set has there, then recomputes ``stat_fn``.
    """
    if rng is None:
        rng = np.random.default_rng()
    pool = pool.reset_index(drop=True)
    strata = pd.Series(strata).reset_index(drop=True)
    target_mask = np.asarray(target_mask, dtype=bool)

    observed = float(stat_fn(pool[target_mask]))

    # Build per-stratum index pools and the target's per-stratum counts.
    stratum_to_idx: dict[str, np.ndarray] = {
        s: np.where(strata.values == s)[0] for s in pd.unique(strata)
    }
    target_strata = strata[target_mask]
    counts = target_strata.value_counts().to_dict()

    null_vals = np.full(n_null, np.nan)
    for b in range(n_null):
        picks: list[int] = []
        ok = True
        for s, c in counts.items():
            avail = stratum_to_idx.get(s, np.array([], dtype=int))
            if avail.size == 0:
                ok = False
                break
            replace = avail.size < c
            picks.extend(rng.choice(avail, size=c, replace=replace).tolist())
        if not ok or not picks:
            continue
        null_vals[b] = stat_fn(pool.iloc[picks])

    finite = null_vals[np.isfinite(null_vals)]
    if finite.size == 0 or not np.isfinite(observed):
        return MatchedNullResult(statistic_name, observed, float("nan"), float("nan"),
                                 float("nan"), float("nan"), float("nan"), 0,
                                 int(target_mask.sum()))
    med = float(np.median(finite))
    lo = float(np.percentile(finite, 100 * alpha / 2))
    hi = float(np.percentile(finite, 100 * (1 - alpha / 2)))
    p_greater = float((np.sum(finite >= observed) + 1) / (finite.size + 1))
    # two-sided around the null median
    dev = abs(observed - med)
    p_two = float((np.sum(np.abs(finite - med) >= dev) + 1) / (finite.size + 1))
    return MatchedNullResult(statistic_name, observed, med, lo, hi, p_greater, p_two,
                             int(finite.size), int(target_mask.sum()))


# ─────────────────────────────────────────────────────────────────────────────
# Pathway-vector cosine recurrence (Layer 4)
# ─────────────────────────────────────────────────────────────────────────────

def pathway_effect_vector(
    gene_effect: pd.Series,
    gene_sets: dict[str, list[str]],
    min_genes: int = 3,
) -> tuple[pd.Series, pd.DataFrame]:
    """Mean gene effect per pathway -> ordered pathway vector + coverage table.

    ``gene_effect`` is indexed by the same id space as the gene-set members
    (e.g. ENSMUSG).  Pathways with fewer than ``min_genes`` mapped, finite
    members are dropped from the vector.
    """
    records = []
    values = {}
    for name, members in gene_sets.items():
        present = [g for g in members if g in gene_effect.index]
        vals = gene_effect.reindex(present).astype(float)
        vals = vals[np.isfinite(vals)]
        records.append({"pathway": name, "n_members": len(members),
                        "n_mapped": int(vals.size), "mean_effect": float(vals.mean()) if vals.size else np.nan})
        if vals.size >= min_genes:
            values[name] = float(vals.mean())
    vec = pd.Series(values)
    return vec, pd.DataFrame(records)


def aligned_pathway_cosine(vec_a: pd.Series, vec_b: pd.Series) -> tuple[float, list[str]]:
    """Cosine between two pathway vectors on their shared, finite pathways."""
    from src.networks.contrast_vectors import cosine
    shared = [p for p in vec_a.index if p in vec_b.index]
    a = vec_a.reindex(shared).to_numpy(dtype=float)
    b = vec_b.reindex(shared).to_numpy(dtype=float)
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 2:
        return float("nan"), [p for p, m in zip(shared, mask) if m]
    return float(cosine(a[mask], b[mask])), [p for p, m in zip(shared, mask) if m]


__all__ = [
    "TmtTable",
    "parse_tmt_sheet",
    "compute_flight_effect",
    "compute_site_flight_effect_lm",
    "collapse_to_gene",
    "spearman",
    "sign_agreement",
    "assign_match_strata",
    "MatchedNullResult",
    "matched_null_test",
    "pathway_effect_vector",
    "aligned_pathway_cosine",
]
