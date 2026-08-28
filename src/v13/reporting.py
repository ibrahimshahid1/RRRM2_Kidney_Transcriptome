"""Read-only reporting for v13 continuous phosphoproteomic inference outputs.

This module deliberately consumes the frozen inference artifacts instead of
recomputing any site, parent-gene, permutation, multiplicity, or claim-gate
quantity.  Its only derived values are descriptive joins, observable-member
counts for sets omitted as non-evaluable, labels, and display ordering.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "rrrm2-mpl-cache")
)
os.environ.setdefault(
    "XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "rrrm2-xdg-cache")
)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

plt.rcParams["axes.unicode_minus"] = False


PRIMARY_PROFILE = "primary"
PRIMARY_EXCLUSION = "primary"
DEFAULT_PRIMARY_SCORE = "median_negative_site_effect"

DISPLAY_NAMES = {
    "proximal_tubule": "Proximal tubule",
    "thick_ascending_limb": "Thick ascending limb",
    "DCT1": "DCT1",
    "DCT2_CNT_transition": "DCT2/CNT transition",
    "DCT2_CNT_external_validation": "DCT2/CNT external",
    "ASDN": "ASDN",
    "principal_cell": "Principal cell",
    "intercalated_cell": "Intercalated cell",
    "podocyte": "Podocyte",
    "endothelial": "Endothelial",
    "fibroblast": "Fibroblast",
    "immune": "Immune",
}

FOREST_ORDER = [
    "proximal_tubule",
    "thick_ascending_limb",
    "DCT1",
    "DCT2_CNT_transition",
    "ASDN",
    "principal_cell",
    "intercalated_cell",
    "podocyte",
    "endothelial",
    "fibroblast",
    "immune",
]

SUBTYPE_ANNOTATION_SETS = [
    "DCT1_core",
    "DCT1",
    "strict_DCT2_peaked",
    "DCT2_CNT_transition",
    "DCT2_CNT_external_validation",
    "DCT2_CNT_atlas_context",
]

SUBTYPE_DISPLAY = {
    "DCT1_core": "DCT1 core",
    "DCT1": "DCT1",
    "strict_DCT2_peaked": "strict DCT2",
    "DCT2_CNT_transition": "transition",
    "DCT2_CNT_external_validation": "external",
    "DCT2_CNT_atlas_context": "atlas context",
}

REQUIRED_FILES = {
    "set_results": "set_level_permutation_inference.tsv",
    "gene_results": "parent_gene_null_calibration.tsv",
    "membership": "gene_set_membership_frozen_copy.tsv",
    "claim_gates": "claim_gates.tsv",
    "claim_tier": "claim_tier.tsv",
    "manifest": "manifest.json",
}


def _read_tsv(path: Path, required: Iterable[str] = ()) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required inference artifact is missing: {path}")
    frame = pd.read_csv(path, sep="\t")
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return frame


def _truthy(values: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)
    return (
        values.astype(str)
        .str.strip()
        .str.lower()
        .isin({"true", "1", "yes", "y"})
    )


def _scalar_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(str(value) for value in values))


def _format_number(value: Any, digits: int = 3) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(numeric):
        return "NA"
    if numeric == 0:
        return "0"
    if abs(numeric) < 0.001:
        return f"{numeric:.2e}"
    return f"{numeric:.{digits}f}"


def _format_probability(value: Any) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return "NA"
    if not np.isfinite(numeric):
        return "NA"
    if numeric < 0.001:
        return f"{numeric:.2e}"
    return f"{numeric:.3f}"


def _format_gate(value: Any) -> str:
    if value is None or pd.isna(value) or str(value).strip() == "":
        return "NA"
    return str(value)


def _markdown_escape(value: Any) -> str:
    return str(value).replace("|", r"\|").replace("\n", " ")


def _load_manifest(run_dir: Path) -> dict[str, Any]:
    path = run_dir / REQUIRED_FILES["manifest"]
    if not path.exists():
        raise FileNotFoundError(f"Required inference artifact is missing: {path}")
    with path.open() as handle:
        return json.load(handle)


def _resolve_recorded_path(
    recorded: str | Path | None, run_dir: Path
) -> Path | None:
    if not recorded:
        return None
    path = Path(recorded).expanduser()
    candidates = [path] if path.is_absolute() else [Path.cwd() / path]
    if not path.is_absolute():
        candidates.extend(parent / path for parent in [run_dir, *run_dir.parents])
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def _load_config(manifest: Mapping[str, Any], run_dir: Path) -> dict[str, Any]:
    path = _resolve_recorded_path(manifest.get("config"), run_dir)
    if path is None:
        return {}
    with path.open() as handle:
        loaded = yaml.safe_load(handle)
    return loaded if isinstance(loaded, dict) else {}


def load_inference_artifacts(run_dir: str | Path) -> dict[str, Any]:
    """Load and minimally validate a completed smoke or exact inference run."""

    root = Path(run_dir).expanduser().resolve()
    if not root.is_dir():
        raise NotADirectoryError(f"Inference output directory not found: {root}")
    manifest = _load_manifest(root)
    artifacts: dict[str, Any] = {
        "run_dir": root,
        "manifest": manifest,
        "config": _load_config(manifest, root),
    }
    artifacts["set_results"] = _read_tsv(
        root / REQUIRED_FILES["set_results"],
        {
            "profile",
            "exclusion",
            "gene_score",
            "set_test",
            "gene_set",
            "n_observable_genes",
            "observed_statistic",
            "empirical_p_greater",
            "maxT_fwer",
            "comparator_bh_q",
        },
    )
    artifacts["gene_results"] = _read_tsv(
        root / REQUIRED_FILES["gene_results"],
        {
            "profile",
            "gene_symbol",
            "gene_score",
            "observed_raw_score",
            "observed_gene_z",
            "n_sites",
            "fixed_null_universe_eligible",
        },
    )
    artifacts["membership"] = _read_tsv(
        root / REQUIRED_FILES["membership"],
        {"gene_symbol", "gene_set", "final_for_testing"},
    )
    artifacts["claim_gates"] = _read_tsv(
        root / REQUIRED_FILES["claim_gates"],
        {"gene_set", "claim_gate_status", "boundary"},
    )
    artifacts["claim_tier"] = _read_tsv(
        root / REQUIRED_FILES["claim_tier"], {"claim_tier", "rationale"}
    )
    artifacts["run_summary"] = (
        json.loads((root / "run_summary.json").read_text())
        if (root / "run_summary.json").exists()
        else {}
    )
    return artifacts


def _primary_score(artifacts: Mapping[str, Any]) -> str:
    return str(
        artifacts["manifest"].get("primary_gene_score", DEFAULT_PRIMARY_SCORE)
    )


def _minimum_observable(artifacts: Mapping[str, Any]) -> int:
    config = artifacts.get("config", {})
    return int(config.get("set_test", {}).get("minimum_observable_genes", 8))


def _primary_sets(artifacts: Mapping[str, Any]) -> list[str]:
    configured = (
        artifacts.get("config", {}).get("set_test", {}).get("primary_family", [])
    )
    if configured:
        return _unique(configured)
    return _unique(artifacts["claim_gates"]["gene_set"].dropna())


def _comparator_sets(artifacts: Mapping[str, Any]) -> list[str]:
    configured = artifacts.get("config", {}).get("comparator_sets", [])
    if configured:
        return _unique(configured)
    rows = artifacts["set_results"]
    if "set_role" in rows:
        rows = rows[rows["set_role"].eq("kidney_comparator")]
    return _unique(rows["gene_set"].dropna())


def _exclusion_genes(
    artifacts: Mapping[str, Any], exclusion: str
) -> set[str]:
    if exclusion == "none":
        return set()
    values = (
        artifacts.get("config", {})
        .get("canonical_exclusions", {})
        .get(exclusion, [])
    )
    return {str(value) for value in values}


def _final_membership(artifacts: Mapping[str, Any]) -> pd.DataFrame:
    membership = artifacts["membership"].copy()
    membership = membership[_truthy(membership["final_for_testing"])].copy()
    membership["gene_symbol"] = membership["gene_symbol"].astype(str)
    membership["gene_set"] = membership["gene_set"].astype(str)
    return membership


def observable_members(
    artifacts: Mapping[str, Any],
    gene_set: str,
    *,
    profile: str,
    exclusion: str,
    gene_score: str,
) -> list[str]:
    """Return descriptive observable members without recalculating inference."""

    membership = _final_membership(artifacts)
    members = set(
        membership.loc[membership["gene_set"].eq(gene_set), "gene_symbol"]
    )
    gene_results = artifacts["gene_results"]
    eligible = gene_results[
        gene_results["profile"].eq(profile)
        & gene_results["gene_score"].eq(gene_score)
        & _truthy(gene_results["fixed_null_universe_eligible"])
        & pd.to_numeric(
            gene_results["observed_gene_z"], errors="coerce"
        ).notna()
    ]
    observed = set(eligible["gene_symbol"].astype(str))
    members -= _exclusion_genes(artifacts, exclusion)
    return sorted(members & observed)


def _select_set_result(
    artifacts: Mapping[str, Any],
    gene_set: str,
    *,
    profile: str,
    exclusion: str,
    gene_score: str,
) -> pd.Series | None:
    rows = artifacts["set_results"]
    selected = rows[
        rows["profile"].eq(profile)
        & rows["exclusion"].eq(exclusion)
        & rows["gene_score"].eq(gene_score)
        & rows["set_test"].eq("competitive")
        & rows["gene_set"].eq(gene_set)
    ]
    if len(selected) > 1:
        raise ValueError(
            "Inference output contains duplicate competitive result rows for "
            f"{profile}/{exclusion}/{gene_score}/{gene_set}"
        )
    return selected.iloc[0] if len(selected) else None


def _adjusted_probability(
    row: pd.Series | None, gene_set: str, primary_sets: set[str]
) -> tuple[float, str]:
    if row is None:
        return float("nan"), "not_evaluable"
    if gene_set in primary_sets:
        return float(row.get("maxT_fwer", np.nan)), "maxT_FWER"
    comparator_q = float(row.get("comparator_bh_q", np.nan))
    if np.isfinite(comparator_q):
        return comparator_q, "BH_q"
    return float(row.get("empirical_p_greater", np.nan)), "unadjusted"


def _is_smoke(artifacts: Mapping[str, Any]) -> bool:
    mode = (
        artifacts["manifest"].get("permutation", {}).get("mode", "")
    )
    return bool(artifacts.get("run_summary", {}).get("smoke", False)) or (
        str(mode).lower() != "exact"
    )


def build_claim_decision_summary(
    artifacts: Mapping[str, Any],
) -> pd.DataFrame:
    """Create one compact reporting row per frozen primary claim set."""

    manifest = artifacts["manifest"]
    claim_tier = artifacts["claim_tier"].iloc[0]
    gates = artifacts["claim_gates"]
    score = _primary_score(artifacts)
    primary_sets = set(_primary_sets(artifacts))
    permutation = manifest.get("permutation", {})
    smoke = _is_smoke(artifacts)
    gate_columns = [
        column for column in gates.columns if column.startswith("gate_")
    ]
    rows: list[dict[str, Any]] = []
    for _, gate in gates.iterrows():
        gene_set = str(gate["gene_set"])
        result = _select_set_result(
            artifacts,
            gene_set,
            profile=PRIMARY_PROFILE,
            exclusion=PRIMARY_EXCLUSION,
            gene_score=score,
        )
        members = observable_members(
            artifacts,
            gene_set,
            profile=PRIMARY_PROFILE,
            exclusion=PRIMARY_EXCLUSION,
            gene_score=score,
        )
        adjusted, adjustment = _adjusted_probability(
            result, gene_set, primary_sets
        )
        failed = [
            column.removeprefix("gate_")
            for column in gate_columns
            if str(gate.get(column, "")).lower() == "fail"
        ]
        non_evaluable = [
            column.removeprefix("gate_")
            for column in gate_columns
            if str(gate.get(column, "")).lower() == "non_evaluable"
        ]
        if gene_set == "DCT2_CNT_transition":
            statistical_eligible = _scalar_bool(
                claim_tier.get("statistical_dct2_title_eligible", False)
            )
            permitted = _scalar_bool(
                claim_tier.get("dct2_title_permitted", False)
            )
        elif gene_set == "ASDN":
            statistical_eligible = _scalar_bool(
                claim_tier.get("statistical_asdn_claim_eligible", False)
            )
            permitted = _scalar_bool(
                claim_tier.get("asdn_beyond_axis_claim_permitted", False)
            )
        else:
            statistical_eligible = (
                str(gate["claim_gate_status"]) == "pass"
            )
            permitted = str(gate["claim_gate_status"]) == "pass"
        coverage_gate = str(
            gate.get("gate_minimum_observable_genes", "")
        ).lower()
        claim_gate_status = str(gate["claim_gate_status"]).lower()
        if coverage_gate == "fail" or (
            claim_gate_status == "non_evaluable"
            and len(members) < _minimum_observable(artifacts)
        ):
            evaluation_status = "coverage_non_evaluable"
            interpretation = (
                "Below the frozen observable-gene threshold; this is "
                "non-evaluability, not evidence of a null effect."
            )
        elif claim_gate_status == "fail":
            evaluation_status = "evaluated_gate_failure"
            interpretation = (
                "Statistically evaluable, but at least one frozen claim gate "
                "failed."
            )
        elif claim_gate_status == "pass":
            evaluation_status = "statistically_eligible"
            interpretation = "Every applicable frozen statistical gate passed."
        else:
            evaluation_status = "non_evaluable_other"
            interpretation = "One or more required statistical gates were unavailable."
        rows.append(
            {
                "analysis_id": manifest.get("analysis_id", ""),
                "permutation_mode": permutation.get("mode", ""),
                "n_assignments_run": permutation.get("n_assignments_run", np.nan),
                "inferential_status": (
                    "smoke_only_not_citable" if smoke else "exact_run"
                ),
                "claim_tier": claim_tier["claim_tier"],
                "gene_set": gene_set,
                "claim_gate_status": gate["claim_gate_status"],
                "statistical_evaluation_status": evaluation_status,
                "statistical_interpretation": interpretation,
                "statistical_claim_eligible": statistical_eligible,
                "publication_claim_permitted": permitted,
                "claim_permitted": permitted,
                "publication_promotion_status": claim_tier.get(
                    "publication_promotion_status", ""
                ),
                "publication_design_and_provenance_gate": claim_tier.get(
                    "publication_design_and_provenance_gate", ""
                ),
                "publication_blockers": claim_tier.get(
                    "publication_blockers", ""
                ),
                "n_observable_genes": len(members),
                "minimum_observable_genes": _minimum_observable(artifacts),
                "observed_competitive_statistic": (
                    result["observed_statistic"] if result is not None else np.nan
                ),
                "empirical_p_greater": (
                    result["empirical_p_greater"] if result is not None else np.nan
                ),
                "family_adjusted_p": adjusted,
                "multiplicity_adjustment": adjustment,
                "n_leave_one_gene_out_rows": gate.get(
                    "n_leave_one_gene_out_rows", np.nan
                ),
                "gate_leave_one_gene_out_completeness": gate.get(
                    "gate_leave_one_gene_out_completeness", ""
                ),
                "gate_all_leave_one_gene_out_positive": gate.get(
                    "gate_all_leave_one_gene_out_positive", ""
                ),
                "gate_centered_vs_uncentered_direction": gate.get(
                    "gate_centered_vs_uncentered_direction", ""
                ),
                "gate_scaled_vs_summed_direction": gate.get(
                    "gate_scaled_vs_summed_direction", ""
                ),
                "gate_normalization_direction_concordance": gate.get(
                    "gate_normalization_direction_concordance", ""
                ),
                "gate_multicompartment_broad_expression_exclusion_direction": (
                    gate.get(
                        "gate_multicompartment_broad_expression_exclusion_direction",
                        "",
                    )
                ),
                "gate_no_unrelated_compartment_equal_or_stronger": gate.get(
                    "gate_no_unrelated_compartment_equal_or_stronger", ""
                ),
                "strongest_unrelated_comparator": gate.get(
                    "strongest_unrelated_comparator", ""
                ),
                "strongest_unrelated_competitive_statistic": gate.get(
                    "strongest_unrelated_competitive_statistic", np.nan
                ),
                "failed_gates": ";".join(failed),
                "non_evaluable_gates": ";".join(non_evaluable),
                "boundary": gate.get("boundary", ""),
                "design_warning": (
                    "Condition is aligned with reporter-tag position; label "
                    "permutation does not remove systematic tag effects."
                    if manifest.get("condition_reporter_position_confounded")
                    else ""
                ),
                "overall_rationale": claim_tier.get("rationale", ""),
            }
        )
    return pd.DataFrame(rows)


def write_claim_decision_markdown(
    summary: pd.DataFrame, path: str | Path
) -> None:
    path = Path(path)
    if summary.empty:
        path.write_text("# V13 claim decision\n\nNo claim rows were available.\n")
        return
    first = summary.iloc[0]
    lines = [
        "# V13 claim-decision summary",
        "",
        f"- Claim tier: **{_markdown_escape(first['claim_tier'])}**",
        f"- Run status: **{_markdown_escape(first['inferential_status'])}**",
        "- Publication promotion: "
        f"**{_markdown_escape(first['publication_promotion_status'])}**",
        f"- Overall rationale: {_markdown_escape(first['overall_rationale'])}",
    ]
    blockers = str(first.get("publication_blockers", ""))
    if blockers:
        lines.append(f"- Publication blockers: {_markdown_escape(blockers)}")
    warning = str(first.get("design_warning", ""))
    if warning:
        lines.append(f"- Design limitation: {warning}")
    lines.extend(
        [
            "",
            "| Program | Observable n | Statistic | Empirical p | Adjusted p | Statistical evaluation | Eligible | Publication claim |",
            "|---|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for _, row in summary.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    _markdown_escape(DISPLAY_NAMES.get(row["gene_set"], row["gene_set"])),
                    str(int(row["n_observable_genes"])),
                    _format_number(row["observed_competitive_statistic"]),
                    _format_probability(row["empirical_p_greater"]),
                    _format_probability(row["family_adjusted_p"]),
                    _markdown_escape(row["statistical_evaluation_status"]),
                    (
                        "yes"
                        if bool(row["statistical_claim_eligible"])
                        else "no"
                    ),
                    "yes" if bool(row["claim_permitted"]) else "no",
                ]
            )
            + " |"
        )
    lines.extend(["", "## Decisive gates", ""])
    for _, row in summary.iterrows():
        failures = row["failed_gates"] or "none"
        non_eval = row["non_evaluable_gates"] or "none"
        lines.append(
            f"- **{DISPLAY_NAMES.get(row['gene_set'], row['gene_set'])}:** "
            f"{_markdown_escape(row['statistical_interpretation'])} "
            f"failed = {_markdown_escape(failures)}; "
            f"non-evaluable = {_markdown_escape(non_eval)}."
        )
    lines.extend(
        [
            "",
            "## Selected robustness gates",
            "",
            "| Program | LOO rows | LOO complete | All LOO positive | Centered vs uncentered | Scaled vs summed | Multi-compartment broad exclusion | No unrelated compartment as strong |",
            "|---|---:|---|---|---|---|---|---|",
        ]
    )
    for _, row in summary.iterrows():
        loo_rows = row["n_leave_one_gene_out_rows"]
        loo_rows_display = (
            str(int(loo_rows)) if pd.notna(loo_rows) else "NA"
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    _markdown_escape(
                        DISPLAY_NAMES.get(row["gene_set"], row["gene_set"])
                    ),
                    loo_rows_display,
                    _markdown_escape(
                        _format_gate(
                            row["gate_leave_one_gene_out_completeness"]
                        )
                    ),
                    _markdown_escape(
                        _format_gate(row["gate_all_leave_one_gene_out_positive"])
                    ),
                    _markdown_escape(
                        _format_gate(
                            row["gate_centered_vs_uncentered_direction"]
                        )
                    ),
                    _markdown_escape(
                        _format_gate(row["gate_scaled_vs_summed_direction"])
                    ),
                    _markdown_escape(
                        _format_gate(
                            row[
                                "gate_multicompartment_broad_expression_exclusion_direction"
                            ]
                        )
                    ),
                    _markdown_escape(
                        _format_gate(
                            row[
                                "gate_no_unrelated_compartment_equal_or_stronger"
                            ]
                        )
                    ),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "Boundary: these are parent-protein annotation enrichments in "
            "whole-kidney phosphoproteomics, not cell-of-origin localization.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def build_primary_compartment_table(
    artifacts: Mapping[str, Any],
) -> pd.DataFrame:
    """Build the exact rows displayed by the primary compartment dot plot."""

    score = _primary_score(artifacts)
    primary_sets = set(_primary_sets(artifacts))
    comparators = _comparator_sets(artifacts)
    expected = _unique(
        FOREST_ORDER
        + [name for name in comparators if name not in FOREST_ORDER]
        + [name for name in primary_sets if name not in FOREST_ORDER]
    )
    minimum = _minimum_observable(artifacts)
    rows: list[dict[str, Any]] = []
    for order, gene_set in enumerate(expected):
        result = _select_set_result(
            artifacts,
            gene_set,
            profile=PRIMARY_PROFILE,
            exclusion=PRIMARY_EXCLUSION,
            gene_score=score,
        )
        members = observable_members(
            artifacts,
            gene_set,
            profile=PRIMARY_PROFILE,
            exclusion=PRIMARY_EXCLUSION,
            gene_score=score,
        )
        adjusted, adjustment = _adjusted_probability(
            result, gene_set, primary_sets
        )
        evaluated = (
            result is not None
            and len(members) >= minimum
            and np.isfinite(float(result["observed_statistic"]))
        )
        rows.append(
            {
                "display_order": order,
                "gene_set": gene_set,
                "display_name": DISPLAY_NAMES.get(gene_set, gene_set),
                "set_role": (
                    "primary"
                    if gene_set in primary_sets
                    else "kidney_comparator"
                ),
                "evaluation_status": (
                    "evaluated" if evaluated else "non_evaluable"
                ),
                "n_observable_genes": len(members),
                "minimum_observable_genes": minimum,
                "observed_statistic": (
                    result["observed_statistic"] if result is not None else np.nan
                ),
                "null_ci_low": (
                    result.get("null_ci_low", np.nan)
                    if result is not None
                    else np.nan
                ),
                "null_ci_high": (
                    result.get("null_ci_high", np.nan)
                    if result is not None
                    else np.nan
                ),
                "empirical_p_greater": (
                    result["empirical_p_greater"] if result is not None else np.nan
                ),
                "adjusted_p": adjusted,
                "multiplicity_adjustment": adjustment,
                "n_null_valid": (
                    result.get("n_null_valid", np.nan)
                    if result is not None
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("display_order").reset_index(drop=True)


def plot_primary_compartment_table(
    table: pd.DataFrame,
    output_base: str | Path,
    *,
    smoke: bool,
) -> list[Path]:
    """Draw observed competitive statistics against permutation-null intervals."""

    output_base = Path(output_base)
    n = max(len(table), 1)
    fig, ax = plt.subplots(figsize=(9.8, max(5.2, 0.48 * n + 1.8)))
    y = np.arange(len(table))
    colors = {
        "primary": "#b45309",
        "kidney_comparator": "#2563a6",
    }
    for index, row in table.iterrows():
        ypos = y[index]
        if row["evaluation_status"] != "evaluated":
            ax.scatter(
                [0],
                [ypos],
                marker="x",
                s=52,
                color="#9ca3af",
                linewidth=1.5,
                zorder=4,
            )
            right_label = (
                f"NE (n={int(row['n_observable_genes'])}; "
                f"minimum {int(row['minimum_observable_genes'])})"
            )
        else:
            low = float(row["null_ci_low"])
            high = float(row["null_ci_high"])
            if np.isfinite(low) and np.isfinite(high):
                ax.hlines(
                    ypos,
                    low,
                    high,
                    color="#a8afb8",
                    linewidth=2.2,
                    zorder=1,
                )
            adjusted = float(row["adjusted_p"])
            significant = np.isfinite(adjusted) and adjusted <= 0.05
            face = colors[row["set_role"]] if significant else "white"
            edge = colors[row["set_role"]]
            ax.scatter(
                [float(row["observed_statistic"])],
                [ypos],
                s=68,
                facecolor=face,
                edgecolor=edge,
                linewidth=1.8,
                zorder=4,
            )
            right_label = (
                f"p={_format_probability(row['empirical_p_greater'])}; "
                f"{row['multiplicity_adjustment']}="
                f"{_format_probability(row['adjusted_p'])}; "
                f"n={int(row['n_observable_genes'])}"
            )
        ax.text(
            1.01,
            ypos,
            right_label,
            transform=ax.get_yaxis_transform(),
            va="center",
            ha="left",
            fontsize=8.4,
            color="#374151",
            clip_on=False,
        )
    ax.axvline(0, color="#4b5563", linewidth=1, linestyle="--", zorder=0)
    ax.set_yticks(y, table["display_name"].tolist())
    ax.invert_yaxis()
    ax.set_xlabel(
        "Observed competitive statistic\n"
        "(mean parent-gene Z in set - eligible background)"
    )
    ax.set_title(
        "Primary competitive enrichment across kidney programs",
        loc="left",
        fontsize=13,
        fontweight="bold",
        pad=16,
    )
    ax.text(
        0,
        1.012,
        "Bars: 95% permutation-null intervals "
        "(not effect-size confidence intervals).",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=8.5,
        color="#4b5563",
    )
    if smoke:
        ax.text(
            0.5,
            -0.16,
            "SMOKE RUN - p/q values are pipeline checks, not inferential results",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=9,
            color="#b91c1c",
            fontweight="bold",
        )
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e5e7eb", linewidth=0.7)
    fig.subplots_adjust(left=0.24, right=0.69, top=0.86, bottom=0.16)
    outputs = [
        output_base.with_suffix(".png"),
        output_base.with_suffix(".svg"),
        output_base.with_suffix(".pdf"),
    ]
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[1], bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[2], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return outputs


def _broad_expression_flags(
    artifacts: Mapping[str, Any],
) -> tuple[pd.Series, str]:
    # Prefer the immutable run-internal profile difference. A reference table
    # named in an old manifest can be rebuilt in place after the run.
    gene_results = artifacts["gene_results"]
    score = _primary_score(artifacts)
    primary_rows = gene_results[
        gene_results["profile"].eq(PRIMARY_PROFILE)
        & gene_results["gene_score"].eq(score)
        & _truthy(gene_results["fixed_null_universe_eligible"])
    ]
    broad_profile_rows = gene_results[
        gene_results["profile"].eq(
            "exclude_multicompartment_broad_expression"
        )
        & gene_results["gene_score"].eq(score)
        & _truthy(gene_results["fixed_null_universe_eligible"])
    ]
    if not primary_rows.empty and not broad_profile_rows.empty:
        primary = set(primary_rows["gene_symbol"].astype(str))
        retained = set(broad_profile_rows["gene_symbol"].astype(str))
        return (
            pd.Series(
                {gene: gene not in retained for gene in primary},
                dtype="boolean",
            ),
            "inference_profile_difference",
        )

    manifest = artifacts["manifest"]
    frozen_path = _resolve_recorded_path(
        manifest.get("inputs", {}).get("frozen_gene_sets"),
        artifacts["run_dir"],
    )
    if frozen_path is not None:
        flags_path = frozen_path.parent / "broad_expression_flags.tsv"
        if flags_path.exists():
            flags = pd.read_csv(flags_path, sep="\t")
            if {"gene_symbol", "broadly_expressed"}.issubset(flags.columns):
                values = _truthy(flags["broadly_expressed"])
                return (
                    pd.Series(
                        values.to_numpy(),
                        index=flags["gene_symbol"].astype(str),
                        dtype="boolean",
                    ),
                    "recorded_reference_table",
                )

    return pd.Series(dtype="boolean"), "unavailable"


def build_leading_parent_gene_matrix(
    artifacts: Mapping[str, Any], *, top_n_per_set: int = 20
) -> pd.DataFrame:
    """Join observed primary genes to sensitivity and frozen annotations."""

    score = _primary_score(artifacts)
    primary_sets = [
        value
        for value in _primary_sets(artifacts)
        if value in {"ASDN", "DCT2_CNT_transition"}
    ]
    membership = _final_membership(artifacts)
    gene_results = artifacts["gene_results"].copy()
    primary = gene_results[
        gene_results["profile"].eq(PRIMARY_PROFILE)
        & gene_results["gene_score"].eq(score)
        & _truthy(gene_results["fixed_null_universe_eligible"])
    ].copy()
    primary = primary.drop_duplicates("gene_symbol").set_index("gene_symbol")
    protein = gene_results[
        gene_results["profile"].eq("parent_protein_subtracted")
        & gene_results["gene_score"].eq(score)
        & _truthy(gene_results["fixed_null_universe_eligible"])
    ].copy()
    protein = protein.drop_duplicates("gene_symbol").set_index("gene_symbol")
    membership_by_gene = (
        membership.groupby("gene_symbol")["gene_set"]
        .agg(lambda values: sorted(set(map(str, values))))
        .to_dict()
    )
    category = (
        membership.assign(
            category=membership.get("category", pd.Series("", index=membership.index))
            .fillna("")
            .astype(str)
        )
        .groupby(["gene_symbol", "gene_set"])["category"]
        .agg(lambda values: ";".join(sorted({value for value in values if value})))
        .to_dict()
    )
    broad, broad_source = _broad_expression_flags(artifacts)
    rows: list[dict[str, Any]] = []
    for focal_set in primary_sets:
        observable = observable_members(
            artifacts,
            focal_set,
            profile=PRIMARY_PROFILE,
            exclusion=PRIMARY_EXCLUSION,
            gene_score=score,
        )
        ranked = sorted(
            observable,
            key=lambda gene: float(primary.loc[gene, "observed_gene_z"]),
            reverse=True,
        )
        if top_n_per_set > 0:
            ranked = ranked[:top_n_per_set]
        for rank, gene in enumerate(ranked, start=1):
            row = primary.loc[gene]
            protein_row = protein.loc[gene] if gene in protein.index else None
            annotations = membership_by_gene.get(gene, [])
            subtype = [
                value for value in SUBTYPE_ANNOTATION_SETS if value in annotations
            ]
            primary_z = float(row["observed_gene_z"])
            protein_z = (
                float(protein_row["observed_gene_z"])
                if protein_row is not None
                else np.nan
            )
            rows.append(
                {
                    "focal_gene_set": focal_set,
                    "rank_within_set": rank,
                    "n_observable_set_members_total": len(observable),
                    "gene_symbol": gene,
                    "observed_gene_z": primary_z,
                    "observed_raw_score": row["observed_raw_score"],
                    "gene_empirical_p_greater": row.get(
                        "empirical_p_greater", np.nan
                    ),
                    "n_sites": int(row["n_sites"]),
                    "median_log2_signal": row.get(
                        "median_log2_signal", np.nan
                    ),
                    "mean_missing_fraction": row.get(
                        "mean_missing_fraction", np.nan
                    ),
                    "parent_protein_sensitivity_available": protein_row
                    is not None,
                    "parent_protein_observed_gene_z": protein_z,
                    "parent_protein_observed_raw_score": (
                        protein_row["observed_raw_score"]
                        if protein_row is not None
                        else np.nan
                    ),
                    "parent_protein_n_sites": (
                        int(protein_row["n_sites"])
                        if protein_row is not None
                        else np.nan
                    ),
                    "parent_protein_direction_concordant": (
                        bool(
                            np.sign(primary_z) == np.sign(protein_z)
                            and np.sign(primary_z) != 0
                        )
                        if np.isfinite(protein_z)
                        else pd.NA
                    ),
                    "ASDN_category": category.get((gene, "ASDN"), ""),
                    "subtype_annotations": ";".join(subtype),
                    "all_frozen_annotations": ";".join(annotations),
                    "in_DCT2_CNT_external_validation": (
                        "DCT2_CNT_external_validation" in annotations
                    ),
                    "in_strict_DCT2_peaked": (
                        "strict_DCT2_peaked" in annotations
                    ),
                    "broadly_expressed": (
                        bool(broad.loc[gene]) if gene in broad.index else pd.NA
                    ),
                    "broad_flag_source": broad_source,
                }
            )
    return pd.DataFrame(rows)


def plot_leading_parent_gene_matrix(
    table: pd.DataFrame,
    output_base: str | Path,
    *,
    smoke: bool,
) -> list[Path]:
    output_base = Path(output_base)
    if table.empty:
        fig, ax = plt.subplots(figsize=(8, 2.5))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No ASDN or DCT2/CNT parent genes were observable after the primary exclusion.",
            ha="center",
            va="center",
        )
    else:
        values = table[
            ["observed_gene_z", "parent_protein_observed_gene_z"]
        ].to_numpy(dtype=float)
        finite = np.abs(values[np.isfinite(values)])
        limit = max(2.0, float(np.quantile(finite, 0.95))) if finite.size else 2.0
        norm = mcolors.TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
        cmap = plt.get_cmap("RdBu_r").copy()
        cmap.set_bad("#e5e7eb")
        n = len(table)
        fig = plt.figure(figsize=(13.8, max(4.8, 0.42 * n + 2.2)))
        grid = fig.add_gridspec(1, 2, width_ratios=[2.5, 7.5], wspace=0.03)
        ax = fig.add_subplot(grid[0, 0])
        detail = fig.add_subplot(grid[0, 1], sharey=ax)
        ax.imshow(
            np.ma.masked_invalid(values),
            aspect="auto",
            cmap=cmap,
            norm=norm,
            interpolation="none",
        )
        labels = [
            f"{DISPLAY_NAMES.get(row.focal_gene_set, row.focal_gene_set)} | {row.gene_symbol}"
            for row in table.itertuples()
        ]
        ax.set_yticks(np.arange(n), labels, fontsize=8.2)
        ax.set_xticks(
            [0, 1], ["Primary\nparent-gene Z", "Parent-protein-\nadjusted Z"]
        )
        ax.xaxis.tick_top()
        ax.tick_params(axis="x", labelsize=8.5, length=0)
        ax.tick_params(axis="y", length=0)
        for i in range(n):
            for j in range(2):
                if np.isfinite(values[i, j]):
                    ax.text(
                        j,
                        i,
                        f"{values[i, j]:.2f}",
                        ha="center",
                        va="center",
                        fontsize=7.6,
                        color=(
                            "white"
                            if abs(values[i, j]) > limit * 0.52
                            else "#111827"
                        ),
                    )
                else:
                    ax.text(j, i, "NA", ha="center", va="center", fontsize=7)
        ax.spines[:].set_visible(False)
        detail.set_xlim(0, 1)
        detail.set_ylim(n - 0.5, -0.5)
        detail.axis("off")
        headers = [
            (0.02, "Raw score"),
            (0.15, "Sites"),
            (0.24, "Broad"),
            (0.34, "ASDN category / subtype annotation"),
        ]
        for x, header in headers:
            detail.text(
                x,
                1.015,
                header,
                transform=detail.transAxes,
                fontsize=8.3,
                fontweight="bold",
                ha="left",
                va="bottom",
            )
        for i, row in table.reset_index(drop=True).iterrows():
            broad = row["broadly_expressed"]
            broad_text = (
                "yes"
                if pd.notna(broad) and bool(broad)
                else "no"
                if pd.notna(broad)
                else "NA"
            )
            annotation = " / ".join(
                value
                for value in [
                    str(row["ASDN_category"] or "").replace("_", " "),
                    "; ".join(
                        SUBTYPE_DISPLAY.get(value, value)
                        for value in str(
                            row["subtype_annotations"] or ""
                        ).split(";")
                        if value
                    ),
                ]
                if value
            )
            detail.text(
                0.02,
                i,
                _format_number(row["observed_raw_score"]),
                fontsize=7.7,
                va="center",
            )
            detail.text(0.15, i, str(int(row["n_sites"])), fontsize=7.7, va="center")
            detail.text(0.24, i, broad_text, fontsize=7.7, va="center")
            detail.text(
                0.34,
                i,
                annotation or "-",
                fontsize=7.4,
                va="center",
                color="#374151",
            )
            detail.hlines(
                i + 0.5, 0, 1, color="#f0f1f3", linewidth=0.6, clip_on=False
            )
        title = "Leading observable parent genes in frozen ASDN and DCT2/CNT sets"
        fig.suptitle(title, x=0.02, y=0.995, ha="left", fontsize=13, fontweight="bold")
        fig.text(
            0.02,
            0.965,
            "Cell color and value are standardized parent-gene Z "
            "(blue = negative; red = stronger flight-associated suppression).",
            ha="left",
            va="top",
            fontsize=8.5,
            color="#4b5563",
        )
        if smoke:
            fig.text(
                0.98,
                0.99,
                "SMOKE RUN",
                ha="right",
                va="top",
                fontsize=9,
                color="#b91c1c",
                fontweight="bold",
            )
        fig.subplots_adjust(left=0.22, right=0.985, top=0.87, bottom=0.04)
    outputs = [
        output_base.with_suffix(".png"),
        output_base.with_suffix(".svg"),
        output_base.with_suffix(".pdf"),
    ]
    fig.savefig(outputs[0], dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[1], bbox_inches="tight", facecolor="white")
    fig.savefig(outputs[2], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return outputs


def build_robustness_summary(
    artifacts: Mapping[str, Any],
) -> pd.DataFrame:
    """Describe primary-score competitive results across profiles/exclusions."""

    score = _primary_score(artifacts)
    set_results = artifacts["set_results"]
    manifest_profiles = [
        str(value.get("name"))
        for value in artifacts["manifest"].get("profiles", [])
        if isinstance(value, dict) and value.get("name")
    ]
    profiles = _unique(manifest_profiles or set_results["profile"].dropna())
    exclusions = _unique(set_results["exclusion"].dropna())
    focal_sets = _primary_sets(artifacts)
    if "DCT2_CNT_external_validation" in set(
        _final_membership(artifacts)["gene_set"]
    ):
        focal_sets = _unique(focal_sets + ["DCT2_CNT_external_validation"])
    minimum = _minimum_observable(artifacts)
    primary_set_names = set(_primary_sets(artifacts))
    primary_direction: dict[str, float] = {}
    for gene_set in focal_sets:
        result = _select_set_result(
            artifacts,
            gene_set,
            profile=PRIMARY_PROFILE,
            exclusion=PRIMARY_EXCLUSION,
            gene_score=score,
        )
        primary_direction[gene_set] = (
            float(result["observed_statistic"]) if result is not None else np.nan
        )
    nonindependent = bool(
        artifacts["manifest"]
        .get("normalization_equivalence_audit", {})
        .get("independent_robustness_evidence")
        is False
    )
    rows: list[dict[str, Any]] = []
    for gene_set in focal_sets:
        for profile in profiles:
            for exclusion in exclusions:
                result = _select_set_result(
                    artifacts,
                    gene_set,
                    profile=profile,
                    exclusion=exclusion,
                    gene_score=score,
                )
                members = observable_members(
                    artifacts,
                    gene_set,
                    profile=profile,
                    exclusion=exclusion,
                    gene_score=score,
                )
                adjusted, adjustment = _adjusted_probability(
                    result, gene_set, primary_set_names
                )
                statistic = (
                    float(result["observed_statistic"])
                    if result is not None
                    else np.nan
                )
                evaluated = (
                    result is not None
                    and len(members) >= minimum
                    and np.isfinite(statistic)
                )
                reference = primary_direction.get(gene_set, np.nan)
                rows.append(
                    {
                        "gene_set": gene_set,
                        "profile": profile,
                        "exclusion": exclusion,
                        "gene_score": score,
                        "evaluation_status": (
                            "evaluated" if evaluated else "non_evaluable"
                        ),
                        "n_observable_genes": len(members),
                        "minimum_observable_genes": minimum,
                        "observed_statistic": statistic,
                        "direction": (
                            "positive"
                            if statistic > 0
                            else "negative"
                            if statistic < 0
                            else "zero"
                            if np.isfinite(statistic)
                            else "not_evaluable"
                        ),
                        "empirical_p_greater": (
                            result["empirical_p_greater"]
                            if result is not None
                            else np.nan
                        ),
                        "adjusted_p": adjusted,
                        "multiplicity_adjustment": adjustment,
                        "direction_concordant_with_primary": (
                            bool(
                                np.sign(statistic) == np.sign(reference)
                                and np.sign(reference) != 0
                            )
                            if np.isfinite(statistic) and np.isfinite(reference)
                            else pd.NA
                        ),
                        "is_declared_primary_result": (
                            profile == PRIMARY_PROFILE
                            and exclusion == PRIMARY_EXCLUSION
                        ),
                        "independence_note": (
                            "Scaled-versus-summed normalization is an algebraic "
                            "audit, not independent robustness evidence."
                            if nonindependent
                            and profile == "signal_to_noise_sum_centered"
                            else ""
                        ),
                    }
                )
    return pd.DataFrame(rows)


def write_robustness_markdown(
    robustness: pd.DataFrame, path: str | Path, *, smoke: bool
) -> None:
    path = Path(path)
    lines = ["# V13 robustness summary", ""]
    if smoke:
        lines.extend(
            [
                "> This is a smoke-run pipeline report. Its p/q values are not "
                "inferential results.",
                "",
            ]
        )
    lines.extend(
        [
            "The table below keeps the frozen primary parent-gene score and "
            "competitive test fixed while varying analysis profile and "
            "canonical-axis exclusion.",
            "",
            "| Program | Profile | Exclusion | n | Statistic | Direction | p | Adjusted | Status |",
            "|---|---|---|---:|---:|---|---:|---:|---|",
        ]
    )
    for _, row in robustness.iterrows():
        lines.append(
            "| "
            + " | ".join(
                [
                    _markdown_escape(DISPLAY_NAMES.get(row["gene_set"], row["gene_set"])),
                    _markdown_escape(row["profile"]),
                    _markdown_escape(row["exclusion"]),
                    str(int(row["n_observable_genes"])),
                    _format_number(row["observed_statistic"]),
                    _markdown_escape(row["direction"]),
                    _format_probability(row["empirical_p_greater"]),
                    _format_probability(row["adjusted_p"]),
                    _markdown_escape(row["evaluation_status"]),
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "The summed-S/N profile is retained as a normalization provenance "
            "audit. Where the run manifest marks scaled and summed contrasts as "
            "algebraically equivalent, it is not counted as independent evidence.",
            "",
        ]
    )
    path.write_text("\n".join(lines))


def build_report(
    run_dir: str | Path,
    output_dir: str | Path | None = None,
    *,
    top_n_per_set: int = 20,
) -> dict[str, Path]:
    """Build all reporting artifacts without modifying the inference run."""

    artifacts = load_inference_artifacts(run_dir)
    destination = (
        Path(output_dir).expanduser().resolve()
        if output_dir is not None
        else artifacts["run_dir"] / "reporting"
    )
    destination.mkdir(parents=True, exist_ok=True)
    smoke = _is_smoke(artifacts)

    claim = build_claim_decision_summary(artifacts)
    claim_tsv = destination / "claim_decision_summary.tsv"
    claim_md = destination / "claim_decision_summary.md"
    claim.to_csv(claim_tsv, sep="\t", index=False)
    write_claim_decision_markdown(claim, claim_md)

    compartment = build_primary_compartment_table(artifacts)
    compartment_tsv = destination / "primary_compartment_enrichment.tsv"
    compartment.to_csv(compartment_tsv, sep="\t", index=False)
    forest_outputs = plot_primary_compartment_table(
        compartment,
        destination / "primary_compartment_enrichment",
        smoke=smoke,
    )

    leading = build_leading_parent_gene_matrix(
        artifacts, top_n_per_set=top_n_per_set
    )
    leading_tsv = destination / "leading_parent_gene_matrix.tsv"
    leading.to_csv(leading_tsv, sep="\t", index=False)
    leading_outputs = plot_leading_parent_gene_matrix(
        leading,
        destination / "leading_parent_gene_matrix",
        smoke=smoke,
    )

    robustness = build_robustness_summary(artifacts)
    robustness_tsv = destination / "robustness_profiles_exclusions.tsv"
    robustness_md = destination / "robustness_profiles_exclusions.md"
    robustness.to_csv(robustness_tsv, sep="\t", index=False)
    write_robustness_markdown(robustness, robustness_md, smoke=smoke)

    outputs = {
        "claim_summary_tsv": claim_tsv,
        "claim_summary_md": claim_md,
        "compartment_table_tsv": compartment_tsv,
        "compartment_plot_png": forest_outputs[0],
        "compartment_plot_svg": forest_outputs[1],
        "compartment_plot_pdf": forest_outputs[2],
        "leading_matrix_tsv": leading_tsv,
        "leading_matrix_png": leading_outputs[0],
        "leading_matrix_svg": leading_outputs[1],
        "leading_matrix_pdf": leading_outputs[2],
        "robustness_tsv": robustness_tsv,
        "robustness_md": robustness_md,
    }
    report_manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_inference_directory": str(artifacts["run_dir"]),
        "source_analysis_id": artifacts["manifest"].get("analysis_id", ""),
        "source_permutation_mode": artifacts["manifest"]
        .get("permutation", {})
        .get("mode", ""),
        "smoke_only": smoke,
        "analysis_recomputed": False,
        "reporting_boundary": (
            "Read-only rendering of frozen inference outputs; no thresholds, "
            "memberships, statistics, p-values, multiplicity corrections, or "
            "claim gates were changed."
        ),
        "outputs": {key: str(value) for key, value in outputs.items()},
    }
    manifest_path = destination / "report_manifest.json"
    manifest_path.write_text(json.dumps(report_manifest, indent=2) + "\n")
    outputs["report_manifest"] = manifest_path
    return outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Render read-only reports from a v13 inference directory."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Completed smoke or exact v13 inference output directory.",
    )
    parser.add_argument(
        "--output-dir",
        help="Reporting destination (default: <input-dir>/reporting).",
    )
    parser.add_argument(
        "--leading-top-n",
        type=int,
        default=20,
        help="Maximum observable genes displayed per primary set; <=0 keeps all.",
    )
    args = parser.parse_args(argv)
    outputs = build_report(
        args.input_dir,
        args.output_dir,
        top_n_per_set=args.leading_top_n,
    )
    print(json.dumps({key: str(value) for key, value in outputs.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
