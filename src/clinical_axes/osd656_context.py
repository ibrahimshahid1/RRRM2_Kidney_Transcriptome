"""Descriptive OSD-656 urine context for the frozen renal tissue axes.

OSD-656 contains post-flight urine measurements from four Inspiration4 crew
members.  This module intersects that assay with the *already frozen* mouse
kidney tissue-axis genes and preserves subjects and recovery timepoints.  It
does not perform a hypothesis test and its outputs must not be used as
validation or added to the mouse cross-mission meta-analysis.
"""

from __future__ import annotations

import re
from typing import Mapping

import numpy as np
import pandas as pd


def normalize_analyte(value: object) -> str:
    """Normalize gene/analyte labels for an exact symbol-style intersection."""

    return re.sub(r"[^A-Z0-9]", "", str(value).upper())


def frozen_axis_gene_catalog(config: Mapping[str, object]) -> pd.DataFrame:
    """Flatten only the frozen primary-family subdomain genes.

    Sensitivity additions and secondary markers are intentionally excluded.
    """

    family = config.get("primary_family", {})
    if not isinstance(family, Mapping) or not family:
        raise ValueError("config contains no primary_family")
    rows: list[dict[str, object]] = []
    for axis, axis_spec in family.items():
        if not isinstance(axis_spec, Mapping):
            continue
        subdomains = axis_spec.get("subdomains", {})
        if not isinstance(subdomains, Mapping):
            continue
        for subdomain, subdomain_spec in subdomains.items():
            if not isinstance(subdomain_spec, Mapping):
                continue
            genes = subdomain_spec.get("genes", {})
            if not isinstance(genes, Mapping):
                continue
            for gene, sign in genes.items():
                direction = int(sign)
                if direction not in (-1, 1):
                    raise ValueError(f"invalid frozen direction for {gene}: {sign}")
                rows.append(
                    {
                        "axis": str(axis),
                        "subdomain": str(subdomain),
                        "gene_symbol": str(gene),
                        "normalized_symbol": normalize_analyte(gene),
                        "frozen_tissue_axis_sign": direction,
                        "frozen_tissue_direction": "higher" if direction > 0 else "lower",
                    }
                )
    catalog = pd.DataFrame(rows)
    if catalog.empty:
        raise ValueError("no genes found in primary_family subdomains")
    if catalog["normalized_symbol"].duplicated().any():
        duplicates = catalog.loc[
            catalog["normalized_symbol"].duplicated(keep=False), "gene_symbol"
        ].tolist()
        raise ValueError("frozen genes occur more than once: " + ", ".join(duplicates))
    return catalog


def axis_panel_coverage(
    long: pd.DataFrame,
    catalog: pd.DataFrame,
) -> pd.DataFrame:
    """Audit which frozen genes were or were not on the OSD-656 panel."""

    required = {"normalized_marker", "analyte", "subject", "phase", "timepoint"}
    missing = required - set(long.columns)
    if missing:
        raise ValueError("OSD-656 long table is missing: " + ", ".join(sorted(missing)))
    panel = long.copy()
    panel["normalized_marker"] = panel["normalized_marker"].map(normalize_analyte)
    panel_lookup = (
        panel.groupby("normalized_marker", sort=False)
        .agg(
            osd656_analyte=("analyte", "first"),
            n_panel_rows=("analyte", "size"),
            n_subjects=("subject", "nunique"),
            n_preflight_rows=("phase", lambda x: int(x.eq("preflight").sum())),
            n_recovery_rows=("phase", lambda x: int(x.eq("recovery").sum())),
            observed_timepoints=(
                "timepoint",
                lambda x: ";".join(sorted(set(x.dropna().astype(str)))),
            ),
        )
        .reset_index()
    )
    coverage = catalog.merge(
        panel_lookup,
        left_on="normalized_symbol",
        right_on="normalized_marker",
        how="left",
    )
    coverage["panel_status"] = np.where(
        coverage["osd656_analyte"].notna(), "measured_overlap", "absent_from_panel"
    )
    for column in (
        "n_panel_rows",
        "n_subjects",
        "n_preflight_rows",
        "n_recovery_rows",
    ):
        coverage[column] = coverage[column].fillna(0).astype(int)
    coverage["scope_note"] = np.where(
        coverage["panel_status"].eq("measured_overlap"),
        "postflight human urine context only; not a tissue or clinical concentration measurement",
        "no OSD-656 urine measurement for this frozen tissue-axis gene",
    )
    return coverage.drop(columns=["normalized_marker"]).sort_values(
        ["axis", "subdomain", "gene_symbol"]
    ).reset_index(drop=True)


def paired_recovery_changes(
    long: pd.DataFrame,
    coverage: pd.DataFrame,
) -> pd.DataFrame:
    """Pair each recovery value to that subject's mean preflight baseline.

    Multiple preflight collections are first averaged within subject and
    analyte.  Recovery timepoints remain separate; they are never pooled or
    treated as independent subjects.
    """

    required = {
        "analyte",
        "normalized_marker",
        "subject",
        "timepoint",
        "phase",
        "mission_day",
        "concentration_npq",
        "percent_normalized_value",
        "unit",
        "assay_type",
    }
    missing = required - set(long.columns)
    if missing:
        raise ValueError("OSD-656 long table is missing: " + ", ".join(sorted(missing)))
    measured = coverage.loc[coverage["panel_status"].eq("measured_overlap")].copy()
    if measured.empty:
        return pd.DataFrame()
    measured_symbols = set(measured["normalized_symbol"])
    selected = long.copy()
    selected["normalized_marker"] = selected["normalized_marker"].map(normalize_analyte)
    selected = selected[selected["normalized_marker"].isin(measured_symbols)].copy()
    selected["concentration_npq"] = pd.to_numeric(
        selected["concentration_npq"], errors="coerce"
    )
    selected["percent_normalized_value"] = pd.to_numeric(
        selected["percent_normalized_value"], errors="coerce"
    )
    selected["mission_day"] = pd.to_numeric(selected["mission_day"], errors="coerce")

    duplicate = selected.duplicated(
        ["normalized_marker", "subject", "timepoint"], keep=False
    )
    if duplicate.any():
        keys = selected.loc[
            duplicate, ["normalized_marker", "subject", "timepoint"]
        ].drop_duplicates()
        raise ValueError("duplicate subject/analyte/timepoint rows: " + keys.to_json())

    pre = selected[selected["phase"].eq("preflight")].copy()
    baseline = (
        pre.groupby(["normalized_marker", "subject"], sort=False)
        .agg(
            n_preflight_timepoints=("concentration_npq", "count"),
            preflight_timepoints=(
                "timepoint",
                lambda x: ";".join(
                    value
                    for _, value in sorted(
                        zip(
                            pd.to_numeric(
                                x.astype(str).str.extract(r"(\d+)")[0], errors="coerce"
                            ).fillna(-1),
                            x.astype(str),
                        ),
                        reverse=True,
                    )
                ),
            ),
            subject_preflight_mean_npq=("concentration_npq", "mean"),
            subject_preflight_mean_percent_normalized=(
                "percent_normalized_value",
                "mean",
            ),
        )
        .reset_index()
    )
    recovery = selected[selected["phase"].eq("recovery")].copy()
    paired = recovery.merge(
        baseline,
        on=["normalized_marker", "subject"],
        how="left",
        validate="many_to_one",
    )
    paired = paired.merge(
        measured[
            [
                "axis",
                "subdomain",
                "gene_symbol",
                "normalized_symbol",
                "frozen_tissue_axis_sign",
                "frozen_tissue_direction",
            ]
        ],
        left_on="normalized_marker",
        right_on="normalized_symbol",
        how="left",
        validate="many_to_one",
    )
    paired["delta_recovery_minus_subject_preflight_npq"] = (
        paired["concentration_npq"] - paired["subject_preflight_mean_npq"]
    )
    paired["delta_recovery_minus_subject_preflight_percent_points"] = (
        paired["percent_normalized_value"]
        - paired["subject_preflight_mean_percent_normalized"]
    )
    delta = paired["delta_recovery_minus_subject_preflight_npq"]
    paired["observed_urine_direction"] = np.select(
        [delta > 1e-12, delta < -1e-12], ["higher", "lower"], default="flat"
    )
    paired["interpretation_boundary"] = (
        "descriptive postflight urine context only; no validation or inferential weight"
    )
    columns = [
        "axis",
        "subdomain",
        "gene_symbol",
        "analyte",
        "subject",
        "timepoint",
        "mission_day",
        "n_preflight_timepoints",
        "preflight_timepoints",
        "subject_preflight_mean_npq",
        "concentration_npq",
        "delta_recovery_minus_subject_preflight_npq",
        "subject_preflight_mean_percent_normalized",
        "percent_normalized_value",
        "delta_recovery_minus_subject_preflight_percent_points",
        "observed_urine_direction",
        "frozen_tissue_axis_sign",
        "frozen_tissue_direction",
        "unit",
        "assay_type",
        "interpretation_boundary",
    ]
    return paired[columns].sort_values(
        ["axis", "gene_symbol", "mission_day", "subject"]
    ).reset_index(drop=True)


def summarize_paired_changes(paired: pd.DataFrame) -> pd.DataFrame:
    """Describe paired changes by analyte and recovery day without p-values."""

    if paired.empty:
        return pd.DataFrame()
    rows: list[dict[str, object]] = []
    group_columns = [
        "axis",
        "subdomain",
        "gene_symbol",
        "analyte",
        "timepoint",
        "mission_day",
        "frozen_tissue_axis_sign",
        "frozen_tissue_direction",
        "unit",
        "assay_type",
    ]
    for keys, sub in paired.groupby(group_columns, sort=False, dropna=False):
        values = sub["delta_recovery_minus_subject_preflight_npq"].dropna()
        percent = sub[
            "delta_recovery_minus_subject_preflight_percent_points"
        ].dropna()
        n_higher = int((values > 1e-12).sum())
        n_lower = int((values < -1e-12).sum())
        n_flat = int((values.abs() <= 1e-12).sum())
        if n_higher == len(values):
            pattern = "all_higher"
        elif n_lower == len(values):
            pattern = "all_lower"
        else:
            pattern = "mixed"
        row = dict(zip(group_columns, keys, strict=True))
        row.update(
            {
                "n_paired_subjects": int(len(values)),
                "paired_subjects": ";".join(sorted(sub.loc[values.index, "subject"])),
                "mean_delta_npq": float(values.mean()),
                "median_delta_npq": float(values.median()),
                "minimum_delta_npq": float(values.min()),
                "maximum_delta_npq": float(values.max()),
                "n_higher": n_higher,
                "n_lower": n_lower,
                "n_flat": n_flat,
                "subject_pattern": pattern,
                "mean_delta_percent_points": float(percent.mean()),
                "median_delta_percent_points": float(percent.median()),
                "inference_status": "not_performed_by_design",
                "interpretation_boundary": (
                    "descriptive postflight urine context only; no validation or "
                    "inferential weight"
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["axis", "gene_symbol", "mission_day"]
    ).reset_index(drop=True)


def build_osd656_context(
    long: pd.DataFrame,
    config: Mapping[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build coverage, paired-detail, and paired-summary context tables."""

    catalog = frozen_axis_gene_catalog(config)
    coverage = axis_panel_coverage(long, catalog)
    paired = paired_recovery_changes(long, coverage)
    summary = summarize_paired_changes(paired)
    return coverage, paired, summary
