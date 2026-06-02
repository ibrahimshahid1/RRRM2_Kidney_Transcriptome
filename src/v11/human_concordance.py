#!/usr/bin/env python3
"""Reach D -- human urine/fluid-axis concordance.

This module deliberately does *not* validate the mouse kidney omics signature in
human kidney tissue.  The human evidence is urine / clinical chemistry from the
NASA Twins Study, so the analysis is an analyte-level sign concordance check
against the mouse-derived fluid, mineralocorticoid, and distal-nephron transport
predictions.

Inputs on disk:
  * Twins Table S8: urinary electrolyte / biochemistry timepoints.
  * Twins SM PDF Fig. S4A: AQP2 / AGT / RENR urine-protein figure evidence.
  * Optional OSD-656 files: Inspiration4 urine inflammation/NULISAseq recovery
    context, never the main Reach D source.

Outputs:
  * ``human_concordance/twins_axis_concordance.tsv``
  * ``human_concordance/twins_aqp2_figure_evidence.tsv``
  * ``human_concordance/human_concordance_verdict.json``
  * optional OSD-656 catalog/summary files if an OSD-656 directory is present.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy import stats

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

FIGURE_PANEL = "Twins Supplement Fig. S4A Fluid Regulation"

OSD656_DIRECT_AXIS_MARKERS = {
    "AQP2",
    "AGT",
    "ATP6AP2",
    "RENR",
    "SLC12A3",
    "STK39",
    "OXSR1",
    "WNK1",
    "WNK3",
    "WNK4",
    "KLHL3",
}
OSD656_RECOVERY_CONTEXT_MARKERS = {
    "AGER",
    "ANXA1",
    "CCL2",
    "CCL3",
    "CCL4",
    "CCL5",
    "CCL7",
    "CCL8",
    "CCL13",
    "CCL14",
    "CCL15",
    "CCL16",
    "CCL17",
    "CCL19",
    "CCL20",
    "CCL22",
    "CCL23",
    "CCL24",
    "CCL25",
    "CCL28",
    "CHI3L1",
    "CX3CL1",
    "CXCL1",
    "CXCL2",
    "CXCL5",
    "CXCL6",
    "CXCL8",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CXCL12",
    "IL1A",
    "IL1B",
    "IL6",
    "IL10",
    "IL18",
    "LCN2",
    "MMP1",
    "MMP2",
    "MMP3",
    "MMP7",
    "MMP9",
    "MMP10",
    "MMP12",
    "SPP1",
    "TGFA",
    "TGFB1",
    "TIMP1",
    "TIMP2",
    "TNF",
    "VCAM1",
    "VEGFA",
}
OSD656_RECOVERY_PREFIXES = ("CCL", "CXCL", "IL", "MMP", "TIMP")


# Explicitly scored machine-readable analytes.  Rows not listed can still be
# emitted as context if selected, but they do not enter the sign test.
TABLE_ANALYTE_SPECS: dict[str, dict[str, object]] = {
    "Urinary sodium mmol/day": {
        "analyte": "urinary_sodium",
        "axis": "aldosterone_fluid",
        "expected_direction": "up",
        "evidence_type": "machine_readable_table",
        "interpretation": "RAAS/MR down-shift context predicts natriuresis.",
        "scored": True,
    },
    "24 h volume mL": {
        "analyte": "urine_volume_24h",
        "axis": "water_balance",
        "expected_direction": "down",
        "evidence_type": "machine_readable_table",
        "interpretation": "AQP2-up water-retention context predicts lower 24 h urine volume.",
        "scored": True,
    },
    "Urinary magnesium mg/day": {
        "analyte": "urinary_magnesium",
        "axis": "dct_transport_context",
        "expected_direction": "up",
        "evidence_type": "machine_readable_table",
        "interpretation": "DCT transport suppression context predicts reduced Mg reabsorption.",
        "scored": True,
    },
    "Urinary potassium mmol/day": {
        "analyte": "urinary_potassium",
        "axis": "electrolyte_context",
        "expected_direction": "context",
        "evidence_type": "machine_readable_table",
        "interpretation": "Potassium excretion is downstream of mixed MR and distal-delivery effects.",
        "scored": False,
    },
    "Urinary phosphorus mg/day": {
        "analyte": "urinary_phosphorus",
        "axis": "stone_risk_context",
        "expected_direction": "context",
        "evidence_type": "machine_readable_table",
        "interpretation": "Stone-risk / mineral context; not a direct DCT/MR prediction.",
        "scored": False,
    },
    "Creatinine mmol/day": {
        "analyte": "urinary_creatinine",
        "axis": "collection_qc_context",
        "expected_direction": "context",
        "evidence_type": "machine_readable_table",
        "interpretation": "Collection normalization context.",
        "scored": False,
    },
}


FIGURE_EVIDENCE_ROWS: list[dict[str, object]] = [
    {
        "analyte": "AQP2",
        "human_symbol": "AQP2",
        "axis": "water_balance",
        "expected_direction": "up",
        "observed_direction": "up",
        "evidence_type": "figure_level",
        "source_panel": FIGURE_PANEL,
        "scored": True,
        "caveat": "Direction encoded from figure-level evidence; no numeric digitization.",
    },
    {
        "analyte": "AGT",
        "human_symbol": "AGT",
        "axis": "aldosterone_fluid",
        "expected_direction": "down",
        "observed_direction": "down",
        "evidence_type": "figure_level",
        "source_panel": FIGURE_PANEL,
        "scored": True,
        "caveat": "Direction encoded from figure-level evidence; no numeric digitization.",
    },
    {
        "analyte": "RENR",
        "human_symbol": "ATP6AP2",
        "axis": "raas_context",
        "expected_direction": "context",
        "observed_direction": "requires_digitization",
        "evidence_type": "figure_level",
        "source_panel": FIGURE_PANEL,
        "scored": False,
        "caveat": "Panel is recorded as figure-level context; direction left unscored until Fig. S4A is digitized.",
    },
]


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
    return "panel_marker"


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
        "preflight_mean", "inflight_mean", "recovery_mean", "delta_inflight_minus_preflight",
        "n_preflight", "n_inflight", "n_recovery", "interpretation", "caveat",
    ]
    for col in cols:
        if col not in fig.columns:
            fig[col] = ""
        if col not in table_summary.columns:
            table_summary[col] = ""
    return pd.concat([table_summary[cols], fig[cols]], ignore_index=True)


def concordance_verdict(axis: pd.DataFrame, *, osd656_status: str = "not_inspected") -> dict[str, object]:
    scored = axis[axis["scored"].astype(bool)].copy()
    n = int(len(scored))
    k = int(scored["concordant"].astype(bool).sum()) if n else 0
    p = float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue) if n else np.nan
    frac = float(k / n) if n else np.nan
    if n == 0:
        status = "no_scored_human_analytes"
    elif k == n and n >= 3:
        status = "directionally_concordant_all_scored"
    elif frac >= 0.6:
        status = "directionally_concordant_partial"
    else:
        status = "mixed_or_nonconcordant"
    return {
        "analysis": "Reach D -- human urine/fluid-axis concordance",
        "status": status,
        "claim_discipline": (
            "Human urine/clinical chemistry sign concordance only; not human kidney omics validation, "
            "not a regression, and not powered beyond one Twins flight trajectory."
        ),
        "n_scored_analytes": n,
        "n_concordant": k,
        "fraction_concordant": None if not np.isfinite(frac) else frac,
        "binomial_p_greater_than_chance": None if not np.isfinite(p) else p,
        "evidence_types": sorted(axis["evidence_type"].dropna().astype(str).unique().tolist()),
        "osd656_status": osd656_status,
        "interpretation": (
            "Report as figure/table-level human fluid-axis concordance. AQP2 remains figure-level "
            "directional evidence unless Fig. S4A is digitized."
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
            "n_scored_analytes": verdict.get("n_scored_analytes", np.nan),
            "n_concordant": verdict.get("n_concordant", np.nan),
            "fraction_concordant": verdict.get("fraction_concordant", np.nan),
            "binomial_p_greater_than_chance": verdict.get("binomial_p_greater_than_chance", np.nan),
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
