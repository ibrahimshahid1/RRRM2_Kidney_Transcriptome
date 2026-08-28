"""Analysis orchestration helpers for frozen renal tissue axes."""

from __future__ import annotations

from dataclasses import asdict
from typing import Mapping

import numpy as np
import pandas as pd
from scipy import stats

from .data import MissionData, cpm_eligible_genes
from .statistics import (
    AxisScoreResult,
    combine_fixed_effects,
    hedges_g,
    random_effects_reml_mkh,
    score_signed_axis,
)


def axis_directions(axis_spec: Mapping[str, object]) -> dict[str, float]:
    directions: dict[str, float] = {}
    for subdomain in axis_spec["subdomains"].values():
        for gene, sign in subdomain["genes"].items():
            if gene in directions:
                raise ValueError(f"Axis contains duplicate gene {gene}")
            directions[str(gene)] = float(sign)
    return directions


def eligible_axis_spec(
    data: MissionData,
    axis_name: str,
    axis_spec: Mapping[str, object],
    *,
    cpm_threshold: float,
) -> tuple[dict[str, float], dict[str, list[str]], list[dict[str, object]]]:
    """Apply the frozen label-blind observability and subdomain coverage rules."""
    eligible = cpm_eligible_genes(data, cpm_threshold)
    directions: dict[str, float] = {}
    subdomains: dict[str, list[str]] = {}
    audit = []
    for subdomain_name, subdomain in axis_spec["subdomains"].items():
        requested = list(subdomain["genes"])
        used = [
            gene
            for gene in requested
            if gene in eligible and gene in data.expression.index
        ]
        minimum = int(subdomain["minimum_present"])
        audit.append(
            {
                "mission": data.mission,
                "axis": axis_name,
                "subdomain": subdomain_name,
                "n_requested": len(requested),
                "n_used": len(used),
                "minimum_required": minimum,
                "coverage_pass": len(used) >= minimum,
                "genes_used": "|".join(used),
                "genes_missing_or_ineligible": "|".join(
                    gene for gene in requested if gene not in used
                ),
            }
        )
        if len(used) < minimum:
            raise ValueError(
                f"{data.mission}/{axis_name}/{subdomain_name}: "
                f"{len(used)} genes pass observability; need {minimum}"
            )
        subdomains[str(subdomain_name)] = used
        directions.update(
            {gene: float(subdomain["genes"][gene]) for gene in used}
        )
    return directions, subdomains, audit


def score_mission_axes(
    data: MissionData,
    primary_family: Mapping[str, object],
    *,
    cpm_threshold: float,
    method: str = "mean",
) -> tuple[pd.DataFrame, dict[str, AxisScoreResult], pd.DataFrame]:
    scores: dict[str, pd.Series] = {}
    details: dict[str, AxisScoreResult] = {}
    coverage_rows: list[dict[str, object]] = []
    for axis_name, axis_spec in primary_family.items():
        directions, subdomains, audit = eligible_axis_spec(
            data,
            axis_name,
            axis_spec,
            cpm_threshold=cpm_threshold,
        )
        result = score_signed_axis(
            data.expression,
            directions,
            method=method,
            subdomains=subdomains,
            min_genes_per_subdomain=1,
            require_all_genes=True,
        )
        if result.scores.isna().any():
            raise ValueError(f"{data.mission}/{axis_name}: nonfinite animal score")
        scores[axis_name] = result.scores
        details[axis_name] = result
        coverage_rows.extend(audit)
    return pd.DataFrame(scores), details, pd.DataFrame(coverage_rows)


def combined_score_design(
    missions: Mapping[str, MissionData],
    primary_family: Mapping[str, object],
    *,
    cpm_threshold: float,
    method: str = "mean",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, dict[str, AxisScoreResult]]]:
    score_frames = []
    design_frames = []
    coverage = []
    all_details: dict[str, dict[str, AxisScoreResult]] = {}
    for mission, data in missions.items():
        scores, details, audit = score_mission_axes(
            data,
            primary_family,
            cpm_threshold=cpm_threshold,
            method=method,
        )
        all_details[mission] = details
        unique_index = pd.Index(
            [f"{mission}::{sample}" for sample in scores.index], name="sample"
        )
        scores.index = unique_index
        design = pd.DataFrame(
            {
                "mission": mission,
                "stratum": data.metadata["block"].to_numpy(),
                "group": data.metadata["condition"].map(
                    {"FLT": "flight", "GC": "control"}
                ).to_numpy(),
                "animal": data.metadata.index.to_numpy(),
                "context": data.context,
            },
            index=unique_index,
        )
        score_frames.append(scores)
        design_frames.append(design)
        coverage.append(audit)
    return (
        pd.concat(score_frames, axis=0),
        pd.concat(design_frames, axis=0),
        pd.concat(coverage, axis=0, ignore_index=True),
        all_details,
    )


def mission_effect_from_score(
    score: pd.Series, metadata: pd.DataFrame
) -> tuple[dict[str, object], pd.DataFrame]:
    stratum_rows = []
    for block, sub in metadata.groupby("block", sort=False):
        flight = score.loc[sub.index[sub["condition"] == "FLT"]]
        control = score.loc[sub.index[sub["condition"] == "GC"]]
        result = hedges_g(flight, control)
        row = {"stratum": block, **asdict(result)}
        stratum_rows.append(row)
    table = pd.DataFrame(stratum_rows)
    pooled = combine_fixed_effects(table["estimate"], table["variance"])
    summary = {
        "estimate": pooled.estimate,
        "variance": pooled.variance,
        "standard_error": pooled.standard_error,
        "ci_low": pooled.ci_low,
        "ci_high": pooled.ci_high,
        "p_normal": pooled.p,
        "q_within_mission": pooled.q,
        "q_p_within_mission": pooled.q_p,
        "n_strata": pooled.k,
        "n_flight": int((metadata["condition"] == "FLT").sum()),
        "n_control": int((metadata["condition"] == "GC").sum()),
    }
    return summary, table


def descriptive_gene_effects(
    missions: Mapping[str, MissionData],
    details: Mapping[str, Mapping[str, AxisScoreResult]],
) -> pd.DataFrame:
    """Mission effects for signed gene-z values; descriptive, not a test family."""
    rows = []
    for mission, axis_details in details.items():
        data = missions[mission]
        for axis, result in axis_details.items():
            for gene in result.signed_gene_z.index:
                summary, _ = mission_effect_from_score(
                    result.signed_gene_z.loc[gene], data.metadata
                )
                rows.append(
                    {
                        "mission": mission,
                        "axis": axis,
                        "gene": gene,
                        **summary,
                    }
                )
    return pd.DataFrame(rows)


def meta_from_mission_effects(mission_effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for axis, sub in mission_effects.groupby("axis", sort=False):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        rows.append(
            {
                "axis": axis,
                "estimate": fit.estimate,
                "tau2": fit.tau2,
                "standard_error_mkh": fit.standard_error,
                "t_mkh": fit.t,
                "df": fit.df,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh": fit.p,
                "prediction_low": fit.prediction_low,
                "prediction_high": fit.prediction_high,
                "q": fit.q,
                "q_p": fit.q_p,
                "i_squared": fit.i_squared,
                "maximum_weight": float(np.max(fit.weights)),
                "k_missions": fit.k,
            }
        )
    return pd.DataFrame(rows)


def technical_qc_audit(missions: Mapping[str, MissionData]) -> pd.DataFrame:
    """Compare exact FLT/GC analysis cells for frozen technical metrics."""
    rows = []
    metrics = [
        "ratio_genebody_cov_3_to_5",
        "rin",
        "read_depth",
        "uniquely_mapped_percent",
    ]
    for mission, data in missions.items():
        for metric in metrics:
            if metric not in data.qc or data.qc[metric].isna().all():
                rows.append(
                    {
                        "mission": mission,
                        "metric": metric,
                        "evaluable": False,
                    }
                )
                continue
            values = data.qc[metric]
            flight = values.loc[data.metadata.index[data.metadata["condition"] == "FLT"]]
            control = values.loc[data.metadata.index[data.metadata["condition"] == "GC"]]
            if flight.isna().any() or control.isna().any():
                rows.append(
                    {
                        "mission": mission,
                        "metric": metric,
                        "evaluable": False,
                        "n_flight": int(flight.notna().sum()),
                        "n_control": int(control.notna().sum()),
                    }
                )
                continue
            welch = stats.ttest_ind(flight, control, equal_var=False)
            effect = hedges_g(flight, control)
            rows.append(
                {
                    "mission": mission,
                    "metric": metric,
                    "evaluable": True,
                    "n_flight": len(flight),
                    "n_control": len(control),
                    "mean_flight": float(flight.mean()),
                    "mean_control": float(control.mean()),
                    "mean_difference": float(flight.mean() - control.mean()),
                    "hedges_g": effect.estimate,
                    "welch_t": float(welch.statistic),
                    "welch_p": float(welch.pvalue),
                    "imbalance_flag": bool(
                        welch.pvalue < 0.05 and abs(effect.estimate) >= 0.50
                    ),
                }
            )
    return pd.DataFrame(rows)


def sample_manifest(missions: Mapping[str, MissionData]) -> pd.DataFrame:
    rows = []
    for mission, data in missions.items():
        for animal, row in data.metadata.iterrows():
            out = {
                "mission": mission,
                "context": data.context,
                "animal": animal,
                "condition": row["condition"],
                "block": row["block"],
                "source_sample": row["source_sample"],
            }
            for metric in (
                "ratio_genebody_cov_3_to_5",
                "rin",
                "read_depth",
                "uniquely_mapped_percent",
            ):
                if metric in data.qc:
                    out[metric] = data.qc.loc[animal, metric]
            rows.append(out)
    return pd.DataFrame(rows)

