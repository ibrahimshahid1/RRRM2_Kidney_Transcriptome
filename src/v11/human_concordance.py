#!/usr/bin/env python3
"""Reach D -- human urine/fluid-axis concordance."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from src.common import REPO_ROOT
from src.v11.core_analysis import RUN_ROOT, sha256


HUMAN_DIR = Path("data/external/human_spaceflight")
TWINS_TABLE_S8 = HUMAN_DIR / "aau8650_table_s8.xlsx"
TWINS_SM_PDF = HUMAN_DIR / "aau8650_garrett-bakelman_sm.pdf"
OSD656_DIR = HUMAN_DIR / "OSD-656"
OSD656_SUBMITTED = OSD656_DIR / "LSDS-64_Multiplex_urine.immune.AlamarPanel_SUBMITTED.xlsx"
OSD656_TRANSFORMED = OSD656_DIR / "LSDS-64_Multiplex_urine.immune.AlamarPanel_TRANSFORMED.csv"

TABLE_S8_TIMEPOINT_ROW = 4
TABLE_S8_DATA_START_ROW = 5
TABLE_S8_ANALYTE_COL = 0
TABLE_S8_TW_FIRST_COL = 1
TABLE_S8_TW_LAST_COL = 15
TABLE_S8_HR_SUMMARY_COL = 16

#  Externalized, version-controlled prediction / scoring spec
CONFIG_DIR = REPO_ROOT / "config"
PREREG_YAML = CONFIG_DIR / "human_concordance_prereg.yaml"
MARKER_PANEL_YAML = CONFIG_DIR / "human_urine_marker_panel.yaml"


def _load_yaml(path: str | Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def _clean(value: object) -> object:
    """Collapse whitespace in YAML folded/multiline strings to keep TSV cells single-line."""
    if isinstance(value, str):
        return re.sub(r"\s+", " ", value).strip()
    return value


def _resolve_gene_set(ref: Mapping[str, object]) -> list[str]:
    """Resolve a {file, set} reference into a flat symbol list."""
    data = _load_yaml(REPO_ROOT / str(ref["file"]))
    block = data.get(str(ref["set"]), {}) or {}
    genes = list(block.get("genes", []) or [])
    genes += list(block.get("protected_genes", []) or [])
    return [str(g) for g in genes]


def load_prereg(path: str | Path = PREREG_YAML) -> tuple[dict[str, dict[str, object]], list[dict[str, object]], str]:
    """Load the table-analyte spec and digitized figure rows from config."""
    cfg = _load_yaml(path)
    table_specs: dict[str, dict[str, object]] = {}
    for row in cfg.get("table_analytes", []) or []:
        label = str(row["label"])
        table_specs[label] = {k: _clean(v) for k, v in row.items() if k != "label"}
    fig = cfg.get("figure_evidence", {}) or {}
    panel = str(fig.get("source_panel", ""))
    figure_rows: list[dict[str, object]] = []
    for r in fig.get("rows", []) or []:
        figure_rows.append(
            {
                "analyte": r["analyte"],
                "human_symbol": r.get("human_symbol", r["analyte"]),
                "axis": r.get("axis", ""),
                "expected_direction": r.get("expected_direction", "context"),
                "observed_direction": r.get("observed_direction", "requires_digitization"),
                "evidence_type": "figure_level",
                "source_panel": panel,
                "scored": bool(r.get("scored", False)),
                "confounded": bool(r.get("confounded", False)),
                "excluded_reason": _clean(r.get("excluded_reason", "")),
                "caveat": _clean(r.get("caveat", "")),
            }
        )
    return table_specs, figure_rows, panel


def load_marker_panel(
    path: str | Path = MARKER_PANEL_YAML,
) -> tuple[set[str], set[str], tuple[str, ...], str]:
    """Load OSD-656 context categories (composed from sourced gene sets)."""
    cfg = _load_yaml(path)
    cats = cfg.get("categories", {}) or {}

    def build(name: str) -> set[str]:
        block = cats.get(name, {}) or {}
        genes: list[str] = []
        for ref in block.get("from_gene_sets", []) or []:
            genes += _resolve_gene_set(ref)
        for m in block.get("extra_markers", []) or []:
            genes.append(str(m["symbol"]) if isinstance(m, Mapping) else str(m))
        return set(genes)

    direct = build("direct_distal_raas_marker")
    recovery = build("recovery_inflammation_matrix_context")
    fb = cfg.get("family_prefix_fallback", {}) or {}
    prefixes = tuple(str(p) for p in fb.get("prefixes", [])) if fb.get("enabled", False) else tuple()
    default_cat = str(cfg.get("default_category", "panel_marker"))
    return direct, recovery, prefixes, default_cat


TABLE_ANALYTE_SPECS, FIGURE_EVIDENCE_ROWS, FIGURE_PANEL = load_prereg()
(
    OSD656_DIRECT_AXIS_MARKERS,
    OSD656_RECOVERY_CONTEXT_MARKERS,
    OSD656_RECOVERY_PREFIXES,
    OSD656_DEFAULT_CATEGORY,
) = load_marker_panel()


def parse_timepoint(label: object) -> dict[str, object]:
    """Classify a Twins Table S8 timepoint label."""
    s = "" if label is None or (isinstance(label, float) and np.isnan(label)) else str(label).strip()
    qc_excluded = "tube frozen before aliquotting" in s.lower()
    if s.startswith("L-"):
        phase = "preflight"
    elif s.startswith("FD"):
        phase = "inflight"
    elif s.startswith("R+"):
        phase = "recovery"
    else:
        phase = "other"
    day = np.nan
    m = re.search(r"(?:L-|FD|R\+)(\d+)", s)
    if m:
        day = int(m.group(1))
        if s.startswith("L-"):
            day = -day
    return {"timepoint": s, "phase": phase, "day": day, "qc_excluded": qc_excluded}


def parse_osd656_timepoint(timepoint: object, timepoint2: object = None) -> dict[str, object]:
    """Classify an OSD-656 Inspiration4 urine collection timepoint."""
    primary = "" if timepoint is None or (isinstance(timepoint, float) and np.isnan(timepoint)) else str(timepoint).strip()
    secondary = "" if timepoint2 is None or (isinstance(timepoint2, float) and np.isnan(timepoint2)) else str(timepoint2).strip()
    label = primary or secondary
    phase_label = secondary.lower()
    if phase_label == "preflight" or label.startswith("L-"):
        phase = "preflight"
    elif label.startswith("R+") or secondary.startswith("R+"):
        phase = "recovery"
    elif label.startswith("FD") or secondary.startswith("FD"):
        phase = "inflight"
    else:
        phase = "other"
    day = np.nan
    m = re.search(r"(?:L-|FD|R\+)(\d+)", label or secondary)
    if m:
        day = int(m.group(1))
        if (label or secondary).startswith("L-"):
            day = -day
    return {"timepoint": label, "timepoint2": secondary, "phase": phase, "day": day}


def coerce_numeric(value: object) -> float:
    """Coerce a table cell to float, treating below-limit strings as missing."""
    if value is None:
        return float("nan")
    s = str(value).strip()
    if not s or s.lower() in {"nan", "na"} or s.startswith("<"):
        return float("nan")
    try:
        return float(s)
    except ValueError:
        return float("nan")


def split_analyte_unit(label: object) -> tuple[str, str]:
    """Best-effort split of a Table S8 analyte label into name and unit."""
    text = str(label).strip()
    if not text:
        return "", ""
    unit_patterns = [
        r"(.+?)\s+(mmol/day|umol/day|mg/day|mL|mmol/L|umol/L|mg/dL|g/dL|U/L|uIU/mL|pg/ml)$",
    ]
    for pat in unit_patterns:
        m = re.match(pat, text)
        if m:
            return m.group(1).strip(), m.group(2).strip()
    return text, ""


def parse_twins_table_s8(path: str | Path = TWINS_TABLE_S8) -> pd.DataFrame:
    """Return long-format TW flight-subject Table S8 values."""
    raw = pd.read_excel(path, header=None)
    timepoints = {
        col: parse_timepoint(raw.iloc[TABLE_S8_TIMEPOINT_ROW, col])
        for col in range(TABLE_S8_TW_FIRST_COL, TABLE_S8_TW_LAST_COL + 1)
    }
    rows: list[dict[str, object]] = []
    category = ""
    for i in range(TABLE_S8_DATA_START_ROW, len(raw)):
        label = raw.iloc[i, TABLE_S8_ANALYTE_COL]
        if label is None or (isinstance(label, float) and np.isnan(label)):
            continue
        label_s = str(label).strip()
        numeric_values = [coerce_numeric(raw.iloc[i, col]) for col in timepoints]
        if not any(np.isfinite(v) for v in numeric_values):
            category = label_s
            continue
        analyte_name, unit = split_analyte_unit(label_s)
        for col, meta in timepoints.items():
            value = coerce_numeric(raw.iloc[i, col])
            rows.append(
                {
                    "source": "Twins Table S8",
                    "row_index": i,
                    "category": category,
                    "analyte_label": label_s,
                    "analyte_name": analyte_name,
                    "unit": unit,
                    "subject": "TW",
                    "subject_role": "flight_subject",
                    "timepoint": meta["timepoint"],
                    "phase": meta["phase"],
                    "mission_day": meta["day"],
                    "qc_excluded": bool(meta["qc_excluded"]),
                    "value": value,
                }
            )
    return pd.DataFrame(rows)


def observed_direction(delta: float, *, eps: float = 1e-9) -> str:
    if not np.isfinite(delta) or abs(delta) <= eps:
        return "flat"
    return "up" if delta > 0 else "down"


def normalize_marker(label: object) -> str:
    return re.sub(r"[^A-Z0-9]", "", str(label).upper())


def osd656_relevance_category(analyte: object) -> str:
    marker = normalize_marker(analyte)
    if marker in {normalize_marker(m) for m in OSD656_DIRECT_AXIS_MARKERS}:
        return "direct_distal_raas_marker"
    if marker in {normalize_marker(m) for m in OSD656_RECOVERY_CONTEXT_MARKERS}:
        return "recovery_inflammation_matrix_context"
    if any(marker.startswith(prefix) for prefix in OSD656_RECOVERY_PREFIXES):
        return "recovery_inflammation_matrix_context"
    return OSD656_DEFAULT_CATEGORY


def concordant(observed: str, expected: str) -> bool | None:
    if expected not in {"up", "down"} or observed not in {"up", "down"}:
        return None
    return observed == expected


def summarize_table_analytes(
    long: pd.DataFrame,
    specs: Mapping[str, Mapping[str, object]] = TABLE_ANALYTE_SPECS,
) -> pd.DataFrame:
    """Summarize preflight vs inflight direction for selected Table S8 analytes."""
    rows: list[dict[str, object]] = []
    use = long[~long["qc_excluded"].astype(bool)].copy()
    for label, spec in specs.items():
        sub = use[use["analyte_label"].astype(str).eq(label)].copy()
        pre = sub[sub["phase"] == "preflight"]["value"].dropna()
        inflight = sub[sub["phase"] == "inflight"]["value"].dropna()
        recovery = sub[sub["phase"] == "recovery"]["value"].dropna()
        if len(pre) == 0 or len(inflight) == 0:
            delta = float("nan")
            obs = "missing"
        else:
            delta = float(inflight.mean() - pre.mean())
            obs = observed_direction(delta)
        expected = str(spec.get("expected_direction", "context"))
        rows.append(
            {
                "source": "Twins Table S8",
                "analyte": spec.get("analyte", label),
                "analyte_label": label,
                "axis": spec.get("axis", ""),
                "evidence_type": spec.get("evidence_type", "machine_readable_table"),
                "expected_direction": expected,
                "observed_direction": obs,
                "concordant": concordant(obs, expected),
                "scored": bool(spec.get("scored", False)) and expected in {"up", "down"} and obs in {"up", "down"},
                "preflight_mean": float(pre.mean()) if len(pre) else np.nan,
                "inflight_mean": float(inflight.mean()) if len(inflight) else np.nan,
                "recovery_mean": float(recovery.mean()) if len(recovery) else np.nan,
                "delta_inflight_minus_preflight": delta,
                "n_preflight": int(len(pre)),
                "n_inflight": int(len(inflight)),
                "n_recovery": int(len(recovery)),
                "confounded": bool(spec.get("confounded", False)),
                "excluded_reason": str(spec.get("excluded_reason", "")),
                "interpretation": spec.get("interpretation", ""),
                "caveat": "Single flight subject; sign-style context only.",
            }
        )
    return pd.DataFrame(rows)


def figure_evidence_table(rows: Iterable[Mapping[str, object]] = FIGURE_EVIDENCE_ROWS) -> pd.DataFrame:
    out = pd.DataFrame(list(rows))
    out["source"] = "Twins SM PDF"
    out["source_file"] = str(TWINS_SM_PDF)
    out["concordant"] = [concordant(str(o), str(e)) for o, e in zip(out["observed_direction"], out["expected_direction"])]
    if "excluded_reason" not in out.columns:
        out["excluded_reason"] = ""
    out["excluded_reason"] = out["excluded_reason"].fillna("")
    if "confounded" not in out.columns:
        out["confounded"] = False
    out["confounded"] = out["confounded"].fillna(False)
    return out


def build_axis_concordance(table_summary: pd.DataFrame, figure_summary: pd.DataFrame) -> pd.DataFrame:
    fig = figure_summary.copy()
    for col in ["preflight_mean", "inflight_mean", "recovery_mean", "delta_inflight_minus_preflight",
                "n_preflight", "n_inflight", "n_recovery"]:
        if col not in fig.columns:
            fig[col] = np.nan
    if "analyte_label" not in fig.columns:
        fig["analyte_label"] = fig["analyte"]
    cols = [
        "source", "analyte", "analyte_label", "axis", "evidence_type",
        "expected_direction", "observed_direction", "concordant", "scored",
        "confounded", "excluded_reason",
        "preflight_mean", "inflight_mean", "recovery_mean", "delta_inflight_minus_preflight",
        "n_preflight", "n_inflight", "n_recovery", "interpretation", "caveat",
    ]
    for col in cols:
        if col not in fig.columns:
            fig[col] = ""
        if col not in table_summary.columns:
            table_summary[col] = ""
    return pd.concat([table_summary[cols], fig[cols]], ignore_index=True)


def axis_level_concordance(axis: pd.DataFrame) -> pd.DataFrame:
    """Collapse scored analytes to independent physiological axes.

    The scored machine-readable + figure analytes are not statistically
    independent: several index the same physiology (e.g. 24 h urine volume and
    AQP2 both report water balance).  Treating each as its own Bernoulli trial
    inflates the sign test, so we collapse scored rows to their axis and treat
    each axis as a single trial.  An axis is concordant only if *every* scored
    analyte on it is concordant.
    """
    scored = axis[axis["scored"].astype(bool)].copy()
    rows: list[dict[str, object]] = []
    for axis_name, sub in scored.groupby("axis", sort=True):
        conc = sub["concordant"].astype(bool)
        analytes = sorted(sub["analyte"].astype(str).tolist())
        rows.append(
            {
                "axis": str(axis_name),
                "n_scored_analytes": int(len(sub)),
                "analytes": ";".join(analytes),
                "all_analytes_concordant": bool(conc.all()),
                "any_analyte_concordant": bool(conc.any()),
                "axis_concordant": bool(conc.all()),
            }
        )
    return pd.DataFrame(rows)


def concordance_verdict(axis: pd.DataFrame, *, osd656_status: str = "not_inspected") -> dict[str, object]:
    scored = axis[axis["scored"].astype(bool)].copy()
    n_analytes = int(len(scored))
    k_analytes = int(scored["concordant"].astype(bool).sum()) if n_analytes else 0
    axis_summary = axis_level_concordance(axis)
    n_axes = int(len(axis_summary))
    k_axes = int(axis_summary["axis_concordant"].astype(bool).sum()) if n_axes else 0
    p_two_sided = (
        float(stats.binomtest(k_axes, n_axes, 0.5, alternative="two-sided").pvalue)
        if n_axes else np.nan
    )
    frac_axes = float(k_axes / n_axes) if n_axes else np.nan
    if n_axes == 0:
        status = "no_scored_human_axes"
    elif k_axes == n_axes and n_axes >= 3:
        status = "directionally_concordant_all_axes_underpowered"
    elif frac_axes >= 0.6:
        status = "directionally_concordant_partial"
    else:
        status = "mixed_or_nonconcordant"
    return {
        "analysis": "Reach D -- human urine/fluid-axis concordance",
        "status": status,
        "claim_discipline": (
            "Human urine / clinical-chemistry sign concordance only; not human kidney omics "
            "validation and not a regression. The sign test is taken over independent "
            "physiological axes (water balance, aldosterone/fluid, DCT transport), not over "
            "individual analytes, and a single Twins flight trajectory cannot reach significance: "
            "all three axes concordant give an exact two-sided sign-test p = 0.25. The per-analyte "
            "directional trajectories, not the p-value, are the substantive result. AGT (predicted "
            "direction read from the same Fig. S4A panel that supplies the observation) and urinary "
            "calcium (opposing spaceflight bone-loss hypercalciuria confounder) are reported but not "
            "scored."
        ),
        "n_axes_scored": n_axes,
        "n_axes_concordant": k_axes,
        "fraction_axes_concordant": None if not np.isfinite(frac_axes) else frac_axes,
        "sign_test_p_two_sided": None if not np.isfinite(p_two_sided) else p_two_sided,
        "n_analytes_scored": n_analytes,
        "n_analytes_concordant": k_analytes,
        "scored_axes": axis_summary.to_dict(orient="records"),
        "evidence_types": sorted(axis["evidence_type"].dropna().astype(str).unique().tolist()),
        "osd656_status": osd656_status,
        "interpretation": (
            "Report as directional (not significant) human fluid-axis concordance across three "
            "independent axes. AQP2's direction is digitized from Twins Fig. S4A (TW flight twin, "
            "in-flight vs pre-flight; see config/human_concordance_prereg.yaml and the figS4A "
            "provenance crops) and enters the water-balance axis; the single Twins trajectory "
            "remains underpowered (all three axes concordant give exact two-sided p = 0.25)."
        ),
    }


def parse_osd656_submitted(path: str | Path = OSD656_SUBMITTED) -> pd.DataFrame:
    """Return long-format OSD-656 urine inflammation panel data.

    The submitted workbook is an Inspiration4 urine Multiplex/NULISAseq-style
    inflammation panel with preflight and recovery samples.  It is not treated as
    inflight kidney proteomics.
    """
    raw = pd.read_excel(path)
    required = {
        "Analyte",
        "Concentration",
        "Timepoint",
        "ID",
        "Unit",
        "Timepoint2",
        "Type",
        "Percent_normalized_value",
    }
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"{path} missing required OSD-656 columns: {sorted(missing)}")
    rows: list[dict[str, object]] = []
    for _, row in raw.iterrows():
        meta = parse_osd656_timepoint(row["Timepoint"], row["Timepoint2"])
        analyte = str(row["Analyte"]).strip()
        rows.append(
            {
                "source": "OSD-656 urine Multiplex/NULISAseq inflammation panel",
                "analyte": analyte,
                "normalized_marker": normalize_marker(analyte),
                "subject": str(row["ID"]).strip(),
                "timepoint": meta["timepoint"],
                "timepoint2": meta["timepoint2"],
                "phase": meta["phase"],
                "mission_day": meta["day"],
                "unit": row["Unit"],
                "assay_type": row["Type"],
                "concentration_npq": pd.to_numeric(row["Concentration"], errors="coerce"),
                "percent_normalized_value": pd.to_numeric(row["Percent_normalized_value"], errors="coerce"),
                "relevance_category": osd656_relevance_category(analyte),
            }
        )
    return pd.DataFrame(rows)


def _mean_or_nan(values: pd.Series) -> float:
    clean = pd.to_numeric(values, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(clean.mean()) if len(clean) else np.nan


def summarize_osd656_prepost(long: pd.DataFrame) -> pd.DataFrame:
    """Summarize OSD-656 preflight vs recovery directions by marker."""
    if long.empty:
        return pd.DataFrame(
            [{
                "status": "parsed_no_rows",
                "context_use": "optional_recovery_inflammation_context",
                "caveat": "OSD-656 is not used as primary Reach D evidence.",
            }]
        )
    rows: list[dict[str, object]] = []
    for analyte, sub in long.groupby("analyte", sort=True):
        pre = sub[sub["phase"].eq("preflight")]
        rec = sub[sub["phase"].eq("recovery")]
        pre_conc = _mean_or_nan(pre["concentration_npq"])
        rec_conc = _mean_or_nan(rec["concentration_npq"])
        delta_conc = rec_conc - pre_conc if np.isfinite(pre_conc) and np.isfinite(rec_conc) else np.nan
        pre_pct = _mean_or_nan(pre["percent_normalized_value"])
        rec_pct = _mean_or_nan(rec["percent_normalized_value"])
        delta_pct = rec_pct - pre_pct if np.isfinite(pre_pct) and np.isfinite(rec_pct) else np.nan
        rec_days = pd.to_numeric(rec["mission_day"], errors="coerce").dropna()
        category = str(sub["relevance_category"].dropna().iloc[0]) if sub["relevance_category"].notna().any() else "panel_marker"
        rows.append(
            {
                "status": "summarized",
                "source": "OSD-656 urine Multiplex/NULISAseq inflammation panel",
                "analyte": analyte,
                "normalized_marker": normalize_marker(analyte),
                "unit": str(sub["unit"].dropna().iloc[0]) if sub["unit"].notna().any() else "",
                "assay_type": str(sub["assay_type"].dropna().iloc[0]) if sub["assay_type"].notna().any() else "",
                "relevance_category": category,
                "context_use": (
                    "potential_direct_marker_if_present"
                    if category == "direct_distal_raas_marker"
                    else "optional_recovery_inflammation_context"
                ),
                "n_preflight_measurements": int(pre["concentration_npq"].notna().sum()),
                "n_recovery_measurements": int(rec["concentration_npq"].notna().sum()),
                "n_preflight_subjects": int(pre["subject"].dropna().nunique()),
                "n_recovery_subjects": int(rec["subject"].dropna().nunique()),
                "preflight_mean_concentration_npq": pre_conc,
                "recovery_mean_concentration_npq": rec_conc,
                "delta_recovery_minus_preflight_npq": delta_conc,
                "observed_direction_npq": observed_direction(delta_conc) if np.isfinite(delta_conc) else "missing",
                "preflight_mean_percent_normalized": pre_pct,
                "recovery_mean_percent_normalized": rec_pct,
                "delta_percent_normalized": delta_pct,
                "observed_direction_percent_normalized": observed_direction(delta_pct) if np.isfinite(delta_pct) else "missing",
                "earliest_recovery_day": float(rec_days.min()) if len(rec_days) else np.nan,
                "latest_recovery_day": float(rec_days.max()) if len(rec_days) else np.nan,
                "caveat": "I4 urine recovery/inflammation context only; not independent inflight AQP2 validation.",
            }
        )
    out = pd.DataFrame(rows)
    category_order = {
        "direct_distal_raas_marker": 0,
        "recovery_inflammation_matrix_context": 1,
        "panel_marker": 2,
    }
    out["_category_order"] = out["relevance_category"].map(category_order).fillna(9)
    out["_abs_delta"] = out["delta_recovery_minus_preflight_npq"].abs().fillna(-1)
    out = out.sort_values(["_category_order", "_abs_delta", "analyte"], ascending=[True, False, True])
    return out.drop(columns=["_category_order", "_abs_delta"]).reset_index(drop=True)


def _find_osd656_submitted_file(root: Path) -> Path | None:
    preferred = root / OSD656_SUBMITTED.name
    if preferred.exists():
        return preferred
    hits = sorted(root.rglob("*SUBMITTED.xlsx"))
    return hits[0] if hits else None


def inspect_osd656(osd656_dir: str | Path = OSD656_DIR) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Catalog optional OSD-656 files and summarize detectable pre/recovery markers.

    The submitted result workbook is summarized as recovery/inflammation context
    only.  It never enters the primary Twins concordance sign test.
    """
    root = Path(osd656_dir)
    if not root.exists():
        return pd.DataFrame(), pd.DataFrame(), "not_downloaded"
    files = [p for p in root.rglob("*") if p.is_file()]
    index = pd.DataFrame(
        {
            "file_name": [p.name for p in files],
            "relative_path": [str(p.relative_to(root)) for p in files],
            "size_bytes": [p.stat().st_size for p in files],
            "suffix": [p.suffix.lower() for p in files],
        }
    )
    candidate = index[
        index["file_name"].str.contains("nulisa|multiplex|inflamm|urine|protein", case=False, na=False)
        & index["suffix"].isin([".tsv", ".txt", ".csv", ".xlsx"])
    ].copy()
    if candidate.empty:
        summary = pd.DataFrame(
            [{
                "status": "catalog_only_no_parseable_relevant_table",
                "n_candidate_files": 0,
                "summary": "OSD-656 is optional I4 urine inflammation/recovery context; no table parsed.",
            }]
        )
        return index, summary, "catalog_only"
    submitted = _find_osd656_submitted_file(root)
    if submitted is None:
        summary = pd.DataFrame(
            [{
                "status": "candidate_files_present_not_summarized",
                "n_candidate_files": int(len(candidate)),
                "summary": "Relevant-looking OSD-656 files are present, but the submitted urine panel workbook is absent.",
            }]
        )
        return index, summary, "candidate_files_present"
    long = parse_osd656_submitted(submitted)
    summary = summarize_osd656_prepost(long)
    return index, summary, "parsed_prepost_recovery_context"


def osd656_verdict_summary(summary: pd.DataFrame) -> dict[str, object]:
    if summary.empty or "status" not in summary.columns or not summary["status"].astype(str).eq("summarized").any():
        return {}
    direct = summary[summary["relevance_category"].eq("direct_distal_raas_marker")]
    context = summary[summary["relevance_category"].eq("recovery_inflammation_matrix_context")]
    return {
        "n_markers_summarized": int(len(summary)),
        "n_direct_distal_raas_markers": int(len(direct)),
        "n_recovery_inflammation_matrix_markers": int(len(context)),
        "context_framing": "I4 urine short-duration post-flight recovery/inflammation context only.",
    }


def write_human_summary(out_dir: Path, verdict: Mapping[str, object]) -> Path:
    osd = verdict.get("osd656_summary", {})
    osd = osd if isinstance(osd, Mapping) else {}
    summary = pd.DataFrame(
        [{
            "analysis": "human_concordance",
            "status": verdict.get("status", ""),
            "n_axes_scored": verdict.get("n_axes_scored", np.nan),
            "n_axes_concordant": verdict.get("n_axes_concordant", np.nan),
            "fraction_axes_concordant": verdict.get("fraction_axes_concordant", np.nan),
            "sign_test_p_two_sided": verdict.get("sign_test_p_two_sided", np.nan),
            "n_analytes_scored": verdict.get("n_analytes_scored", np.nan),
            "n_analytes_concordant": verdict.get("n_analytes_concordant", np.nan),
            "osd656_status": verdict.get("osd656_status", ""),
            "osd656_n_markers_summarized": osd.get("n_markers_summarized", np.nan),
            "osd656_n_direct_distal_raas_markers": osd.get("n_direct_distal_raas_markers", np.nan),
            "osd656_n_recovery_inflammation_matrix_markers": osd.get(
                "n_recovery_inflammation_matrix_markers", np.nan
            ),
            "claim_discipline": verdict.get("claim_discipline", ""),
        }]
    )
    path = out_dir / "human_concordance_summary.tsv"
    summary.to_csv(path, sep="\t", index=False)
    return path


def run_human_concordance(
    root: str | Path,
    *,
    twins_table_s8: str | Path = TWINS_TABLE_S8,
    twins_sm_pdf: str | Path = TWINS_SM_PDF,
    osd656_dir: str | Path = OSD656_DIR,
) -> dict[str, object]:
    root = Path(root)
    out_dir = root / "human_concordance"
    out_dir.mkdir(parents=True, exist_ok=True)

    long = parse_twins_table_s8(twins_table_s8)
    table_summary = summarize_table_analytes(long)
    fig = figure_evidence_table()
    axis = build_axis_concordance(table_summary, fig)

    long.to_csv(out_dir / "twins_table_s8_long.tsv", sep="\t", index=False)
    fig.to_csv(out_dir / "twins_aqp2_figure_evidence.tsv", sep="\t", index=False)
    axis.to_csv(out_dir / "twins_axis_concordance.tsv", sep="\t", index=False)
    axis_level = axis_level_concordance(axis)
    axis_level.to_csv(out_dir / "twins_axis_level_concordance.tsv", sep="\t", index=False)

    osd_index, osd_summary, osd_status = inspect_osd656(osd656_dir)
    if not osd_index.empty:
        osd_index.to_csv(out_dir / "osd656_file_index.tsv", sep="\t", index=False)
    if not osd_summary.empty:
        osd_summary.to_csv(out_dir / "osd656_prepost_summary.tsv", sep="\t", index=False)

    verdict = concordance_verdict(axis, osd656_status=osd_status)
    verdict["osd656_summary"] = osd656_verdict_summary(osd_summary)
    summary_path = write_human_summary(out_dir, verdict)
    verdict["outputs"] = {
        "twins_axis_concordance": str(out_dir / "twins_axis_concordance.tsv"),
        "twins_axis_level_concordance": str(out_dir / "twins_axis_level_concordance.tsv"),
        "twins_aqp2_figure_evidence": str(out_dir / "twins_aqp2_figure_evidence.tsv"),
        "twins_table_s8_long": str(out_dir / "twins_table_s8_long.tsv"),
        "human_concordance_summary": str(summary_path),
    }
    if not osd_index.empty:
        verdict["outputs"]["osd656_file_index"] = str(out_dir / "osd656_file_index.tsv")
    if not osd_summary.empty:
        verdict["outputs"]["osd656_prepost_summary"] = str(out_dir / "osd656_prepost_summary.tsv")
    manifest_inputs = [Path(twins_table_s8), Path(twins_sm_pdf)]
    verdict["input_manifest"] = [
        {"path": str(p), "sha256": sha256(p)} for p in manifest_inputs if p.exists()
    ]
    (out_dir / "human_concordance_verdict.json").write_text(json.dumps(verdict, indent=2))
    return verdict


def main() -> None:
    ap = argparse.ArgumentParser(description="Reach D: human urine/fluid-axis concordance.")
    ap.add_argument("--run-root", default=str(RUN_ROOT))
    ap.add_argument("--twins-table-s8", default=str(TWINS_TABLE_S8))
    ap.add_argument("--twins-sm-pdf", default=str(TWINS_SM_PDF))
    ap.add_argument("--osd656-dir", default=str(OSD656_DIR))
    args = ap.parse_args()
    verdict = run_human_concordance(
        args.run_root,
        twins_table_s8=args.twins_table_s8,
        twins_sm_pdf=args.twins_sm_pdf,
        osd656_dir=args.osd656_dir,
    )
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
