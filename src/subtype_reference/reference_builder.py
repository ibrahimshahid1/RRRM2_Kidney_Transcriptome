#!/usr/bin/env python3
"""Build frozen distal-nephron signatures without consulting flight results.

The module consumes reference-only differential-expression and whole-kidney
expression summaries.  It deliberately has no dependency on the OSD-462
phosphoproteomic code.  Final data-derived signatures are emitted only when
both the independent distal-nephron validation and whole-kidney specificity
inputs are present; otherwise the output records a non-evaluable gate rather
than silently promoting discovery-only genes.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

FORBIDDEN_INPUT_TOKENS = (
    "flight",
    "osd462",
    "phospho",
    "spaceflight",
    "suppressed_site",
)

DISCOVERY_COLUMNS = {
    "gene_symbol",
    "log2_fc_dct2_vs_dct1",
    "fdr",
    "pct_detected_dct1",
    "pct_detected_dct2",
    "n_consistent_pairs",
    "n_pairs",
}

VALIDATION_COLUMNS = {
    "gene_symbol",
    "log2_fc_dct2_vs_dct1",
    "fdr_dct2_vs_dct1",
    "log2_fc_dct2_vs_cnt",
    "fdr_dct2_vs_cnt",
    "log2_fc_dct1_vs_cnt",
    "fdr_dct1_vs_cnt",
    "pct_detected_dct1",
    "pct_detected_dct2",
    "pct_detected_cnt",
    "n_consistent_dct2_vs_dct1",
    "n_pairs",
}

ATLAS_COLUMNS = {
    "gene_symbol",
    "compartment",
    "mean_cpm",
    "source_study_detection_fraction",
}

SEGMENT_VALIDATION_COLUMNS = {
    "gene_symbol",
    "log2_fc_cnt_vs_dct",
    "ci_low_cnt_vs_dct",
    "p_noninferiority_cnt_vs_dct",
    "fdr_noninferiority_cnt_vs_dct",
    "pct_replicates_detected_cnt",
    "n_consistent_cnt_retained",
    "n_pairs",
}


def load_config(path: str | Path) -> dict[str, Any]:
    """Load and validate the dedicated subtype-reference freeze config."""
    path = Path(path)
    cfg = yaml.safe_load(path.read_text()) or {}
    if not cfg.get("flight_blind", False):
        raise ValueError("Subtype-reference config must declare flight_blind: true")
    if "reference_builder" not in cfg:
        raise ValueError("Config is missing the reference_builder section")
    return cfg


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_commit() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL
        ).strip()
    except Exception:
        return None


def _assert_reference_only(path: Path, frame: pd.DataFrame | None = None) -> None:
    """Refuse obvious outcome leakage into marker construction."""
    lower_path = str(path).lower()
    bad_path = [token for token in FORBIDDEN_INPUT_TOKENS if token in lower_path]
    bad_cols: list[str] = []
    if frame is not None:
        lower_cols = [str(col).lower() for col in frame.columns]
        bad_cols = [
            col
            for col in lower_cols
            if any(token in col for token in FORBIDDEN_INPUT_TOKENS)
        ]
    if bad_path or bad_cols:
        raise ValueError(
            "Reference construction refused a potential outcome-leaking input: "
            f"path_tokens={bad_path}, columns={bad_cols}"
        )


def _read_table(path: Path | None, required: set[str]) -> pd.DataFrame | None:
    if path is None or not path.exists():
        return None
    frame = pd.read_csv(path, sep="\t")
    _assert_reference_only(path, frame)
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    frame = frame.copy()
    frame["gene_symbol"] = frame["gene_symbol"].astype(str)
    return frame


def _broad_expression_table(
    atlas: pd.DataFrame | None, config: dict[str, Any]
) -> pd.DataFrame:
    """Separate true expression breadth from the signature-specificity filter.

    ``broadly_expressed`` requires expression in the distal target and in the
    configured number of unrelated compartments. The compatibility column
    ``non_distal_specific_or_broad`` additionally captures a single unrelated
    compartment whose expression is comparable with the distal target; frozen
    DCT signature construction continues to use that stricter flag.
    """

    columns = [
        "gene_symbol",
        "target_max_cpm",
        "target_expressed",
        "unrelated_max_cpm",
        "unrelated_to_target_ratio",
        "n_unrelated_detected",
        "broadly_expressed",
        "non_distal_specific_or_broad",
    ]
    if atlas is None or atlas.empty:
        return pd.DataFrame(columns=columns)

    broad_cfg = config["reference_builder"]["broad_expression"]
    target = set(broad_cfg["target_compartments"])
    min_cpm = float(broad_cfg["min_cpm"])
    ratio_cutoff = float(broad_cfg["unrelated_to_target_ratio"])
    compartment_ratio_floor = float(
        broad_cfg["unrelated_compartment_ratio_floor"]
    )
    n_cutoff = int(broad_cfg["min_unrelated_compartments"])
    source_study_fraction = float(
        broad_cfg["min_source_study_detection_fraction"]
    )

    rows: list[dict[str, Any]] = []
    for gene, part in atlas.groupby("gene_symbol", sort=False):
        eligible = part[
            pd.to_numeric(
                part["source_study_detection_fraction"], errors="coerce"
            )
            >= source_study_fraction
        ].copy()
        target_values = pd.to_numeric(
            eligible.loc[eligible["compartment"].isin(target), "mean_cpm"],
            errors="coerce",
        )
        unrelated_values = pd.to_numeric(
            eligible.loc[~eligible["compartment"].isin(target), "mean_cpm"],
            errors="coerce",
        )
        target_max = float(target_values.max()) if len(target_values) else 0.0
        target_expressed = target_max >= min_cpm
        unrelated_max = float(unrelated_values.max()) if len(unrelated_values) else 0.0
        ratio = unrelated_max / target_max if target_max > 0 else np.inf
        compartment_floor = max(min_cpm, compartment_ratio_floor * target_max)
        n_unrelated = int((unrelated_values >= compartment_floor).sum())
        broadly_expressed = bool(target_expressed and n_unrelated >= n_cutoff)
        rows.append(
            {
                "gene_symbol": gene,
                "target_max_cpm": target_max,
                "target_expressed": bool(target_expressed),
                "unrelated_max_cpm": unrelated_max,
                "unrelated_to_target_ratio": ratio,
                "n_unrelated_detected": n_unrelated,
                "broadly_expressed": broadly_expressed,
                "non_distal_specific_or_broad": bool(
                    target_expressed
                    and (broadly_expressed or ratio >= ratio_cutoff)
                ),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _curated_membership(config: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    rb = config["reference_builder"]
    for set_name, categories in rb.get("curated_sets", {}).items():
        for category, genes in categories.items():
            for gene in genes:
                rows.append(
                    {
                        "gene_symbol": str(gene),
                        "gene_set": str(set_name),
                        "category": str(category),
                        "definition_source": "frozen_config",
                        "status": "defined",
                    }
                )
    for set_name, genes in rb.get("major_kidney_comparator_sets", {}).items():
        for gene in genes:
            rows.append(
                {
                    "gene_symbol": str(gene),
                    "gene_set": str(set_name),
                    "category": "major_kidney_comparator",
                    "definition_source": "frozen_config",
                    "status": "defined",
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def _atlas_comparator_membership(
    atlas: pd.DataFrame | None, config: dict[str, Any]
) -> pd.DataFrame:
    """Derive major-compartment marker sets from the whole-kidney atlas."""
    columns = [
        "gene_symbol",
        "gene_set",
        "category",
        "definition_source",
        "status",
        "target_cpm",
        "max_other_cpm",
        "log2_target_to_max_other",
        "target_source_study_fraction",
    ]
    if atlas is None or atlas.empty:
        return pd.DataFrame(columns=columns)
    rb = config["reference_builder"]
    targets = rb["atlas_comparator_targets"]
    rule = rb["atlas_comparator_rule"]
    min_cpm = float(rule["min_target_cpm"])
    min_detect = float(rule["min_target_source_study_fraction"])
    min_specificity = float(rule["min_log2_target_to_max_other"])
    max_genes = int(rule["max_genes_per_set"])
    rows: list[pd.DataFrame] = []
    for set_name, target_compartments in targets.items():
        target_mask = atlas["compartment"].isin(target_compartments)
        target = (
            atlas.loc[target_mask]
            .groupby("gene_symbol", as_index=False)
            .agg(
                target_cpm=("mean_cpm", "max"),
                target_source_study_fraction=(
                    "source_study_detection_fraction",
                    "max",
                ),
            )
        )
        other = (
            atlas.loc[~target_mask]
            .groupby("gene_symbol", as_index=False)
            .agg(max_other_cpm=("mean_cpm", "max"))
        )
        table = target.merge(other, on="gene_symbol", how="left")
        table["max_other_cpm"] = table["max_other_cpm"].fillna(0.0)
        table["log2_target_to_max_other"] = np.log2(
            (table["target_cpm"] + 0.1) / (table["max_other_cpm"] + 0.1)
        )
        table = table[
            (table["target_cpm"] >= min_cpm)
            & (table["target_source_study_fraction"] >= min_detect)
            & (table["log2_target_to_max_other"] >= min_specificity)
        ].copy()
        table = table.sort_values(
            ["log2_target_to_max_other", "target_cpm"], ascending=False
        ).head(max_genes)
        table["gene_set"] = set_name
        table["category"] = "major_kidney_comparator"
        table["definition_source"] = "whole_kidney_atlas"
        table["status"] = "defined"
        rows.append(table.reindex(columns=columns))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=columns)


def _assert_asdn_consistency(config: dict[str, Any]) -> None:
    """Require the reference and production configs to freeze the same ASDN set."""
    consistency = config["reference_builder"].get("production_config_consistency", {})
    production_path = consistency.get("production_analysis_config")
    if not production_path:
        return
    production = yaml.safe_load(Path(production_path).read_text()) or {}
    production_genes = set(production.get("asdn_gene_set", {}).get("genes", []))
    categories = config["reference_builder"].get("curated_sets", {}).get("ASDN", {})
    reference_genes = {gene for genes in categories.values() for gene in genes}
    if production_genes != reference_genes:
        raise ValueError(
            "ASDN mismatch between subtype-reference and production configs: "
            f"only_reference={sorted(reference_genes - production_genes)}, "
            f"only_production={sorted(production_genes - reference_genes)}"
        )


def _classify_data_driven(
    discovery: pd.DataFrame | None,
    validation: pd.DataFrame | None,
    segment_validation: pd.DataFrame | None,
    broad: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    rb = config["reference_builder"]
    thresholds = rb["thresholds"]
    empty_cols = [
        "gene_symbol",
        "discovery_dct2_pass",
        "discovery_dct1_pass",
        "discovery_transition_pass",
        "validation_dct2_pass",
        "validation_dct1_pass",
        "validation_dct2_peak_pass",
        "validation_transition_pass",
        "broadly_expressed",
        "non_distal_specific_or_broad",
        "DCT1_core",
        "strict_DCT2_peaked",
        "DCT2_CNT_transition",
        "classification_status",
    ]
    if discovery is None:
        return pd.DataFrame(columns=empty_cols)

    d = discovery.copy()
    numeric = DISCOVERY_COLUMNS - {"gene_symbol"}
    for col in numeric:
        d[col] = pd.to_numeric(d[col], errors="coerce")

    fc = float(thresholds["discovery_abs_log2_fc"])
    fdr = float(thresholds["discovery_fdr"])
    detect = float(thresholds["strict_min_target_detection_fraction"])
    transition_detect = float(
        thresholds["transition_min_target_detection_fraction"]
    )
    transition_fc = float(thresholds["transition_discovery_log2_fc"])
    consistent = int(thresholds["min_consistent_pairs"])
    d["discovery_dct2_pass"] = (
        (d["log2_fc_dct2_vs_dct1"] >= fc)
        & (d["fdr"] <= fdr)
        & (d["pct_detected_dct2"] >= detect)
        & (d["n_consistent_pairs"] >= consistent)
    )
    d["discovery_dct1_pass"] = (
        (d["log2_fc_dct2_vs_dct1"] <= -fc)
        & (d["fdr"] <= fdr)
        & (d["pct_detected_dct1"] >= detect)
        & (d["n_consistent_pairs"] >= consistent)
    )
    d["discovery_transition_pass"] = (
        (d["log2_fc_dct2_vs_dct1"] >= transition_fc)
        & (d["pct_detected_dct2"] >= transition_detect)
        & (d["n_consistent_pairs"] >= consistent)
    )

    if validation is not None:
        v = validation.copy()
        for col in VALIDATION_COLUMNS - {"gene_symbol"}:
            v[col] = pd.to_numeric(v[col], errors="coerce")
        d = d.merge(v, on="gene_symbol", how="left", suffixes=("", "_validation"))
        vfc = float(thresholds["validation_direction_min_abs_log2_fc"])
        vfdr = float(thresholds["validation_fdr"])
        peak_fc = float(thresholds["dct2_peak_vs_cnt_log2_fc"])
        d["validation_dct2_pass"] = (
            (d["log2_fc_dct2_vs_dct1_validation"] >= vfc)
            & (d["fdr_dct2_vs_dct1"] <= vfdr)
            & (d["pct_detected_dct2_validation"] >= detect)
            & (d["n_consistent_dct2_vs_dct1"] >= consistent)
        )
        d["validation_dct1_pass"] = (
            (d["log2_fc_dct2_vs_dct1_validation"] <= -vfc)
            & (d["fdr_dct2_vs_dct1"] <= vfdr)
            & (d["pct_detected_dct1_validation"] >= detect)
        )
        d["validation_dct2_peak_pass"] = (
            (d["log2_fc_dct2_vs_cnt"] >= peak_fc)
            & (d["fdr_dct2_vs_cnt"] <= vfdr)
        )
    else:
        for col in (
            "validation_dct2_pass",
            "validation_dct1_pass",
            "validation_dct2_peak_pass",
        ):
            d[col] = False

    if segment_validation is not None:
        segment = segment_validation.copy()
        for col in SEGMENT_VALIDATION_COLUMNS - {"gene_symbol"}:
            segment[col] = pd.to_numeric(segment[col], errors="coerce")
        d = d.merge(segment, on="gene_symbol", how="left")
        retention_margin = float(thresholds["external_cnt_retention_log2_margin"])
        retention_fdr = float(thresholds["external_retention_fdr"])
        min_detect_reps = float(
            thresholds["external_min_replicate_detection_fraction"]
        )
        d["validation_transition_pass"] = (
            (d["log2_fc_cnt_vs_dct"] >= float(
                thresholds["cnt_retention_log2_fc_floor"]
            ))
            &
            (d["ci_low_cnt_vs_dct"] > retention_margin)
            & (d["fdr_noninferiority_cnt_vs_dct"] <= retention_fdr)
            & (d["pct_replicates_detected_cnt"] >= min_detect_reps)
            & (d["n_consistent_cnt_retained"] >= consistent)
        )
    else:
        d["validation_transition_pass"] = False

    if not broad.empty:
        d = d.merge(
            broad[
                [
                    "gene_symbol",
                    "broadly_expressed",
                    "non_distal_specific_or_broad",
                ]
            ],
            on="gene_symbol",
            how="left",
        )
        d["broadly_expressed"] = d["broadly_expressed"].fillna(True).astype(bool)
        d["non_distal_specific_or_broad"] = (
            d["non_distal_specific_or_broad"].fillna(True).astype(bool)
        )
    else:
        d["broadly_expressed"] = True
        d["non_distal_specific_or_broad"] = True

    has_fine_inputs = validation is not None and not broad.empty
    has_transition_inputs = segment_validation is not None and not broad.empty
    d["DCT1_core"] = (
        d["discovery_dct1_pass"]
        & d["validation_dct1_pass"]
        & ~d["non_distal_specific_or_broad"]
        & has_fine_inputs
    )
    d["strict_DCT2_peaked"] = (
        d["discovery_dct2_pass"]
        & d["validation_dct2_pass"]
        & ~d["non_distal_specific_or_broad"]
        & has_fine_inputs
    )
    d["DCT2_CNT_transition"] = (
        d["discovery_transition_pass"]
        & d["validation_transition_pass"]
        & ~d["non_distal_specific_or_broad"]
        & has_transition_inputs
    )
    d["classification_status"] = np.where(
        has_fine_inputs,
        "complete_fine_and_transition",
        np.where(
            has_transition_inputs,
            "transition_evaluable_fine_subtypes_missing",
            "not_evaluable_missing_validation_or_atlas",
        ),
    )
    return d


def _external_validation_signature(
    validation: pd.DataFrame | None,
    segment_validation: pd.DataFrame | None,
    broad: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    """Derive a GSE150338-only DCT2/CNT signature without GSE228367 membership.

    The fine-subtype and microdissected-segment inputs are both from
    GSE150338. The whole-kidney atlas is used only as the predeclared combined
    distal-specificity/breadth exclusion. No discovery-table membership is
    consulted.
    """

    columns = [
        "gene_symbol",
        "fine_dct2_enrichment_pass",
        "fine_cnt_retention_pass",
        "segment_cnt_retention_pass",
        "broadly_expressed",
        "non_distal_specific_or_broad",
        "DCT2_CNT_external_validation",
    ]
    if (
        validation is None
        or segment_validation is None
        or broad.empty
    ):
        return pd.DataFrame(columns=columns)

    thresholds = config["reference_builder"]["thresholds"]
    v = validation.copy()
    for col in VALIDATION_COLUMNS - {"gene_symbol"}:
        v[col] = pd.to_numeric(v[col], errors="coerce")
    segment = segment_validation.copy()
    for col in SEGMENT_VALIDATION_COLUMNS - {"gene_symbol"}:
        segment[col] = pd.to_numeric(segment[col], errors="coerce")

    table = v.merge(segment, on="gene_symbol", how="inner")
    table = table.merge(
        broad[
            [
                "gene_symbol",
                "broadly_expressed",
                "non_distal_specific_or_broad",
            ]
        ],
        on="gene_symbol",
        how="left",
    )
    table["broadly_expressed"] = (
        table["broadly_expressed"].fillna(True).astype(bool)
    )
    table["non_distal_specific_or_broad"] = (
        table["non_distal_specific_or_broad"].fillna(True).astype(bool)
    )

    transition_fc = float(thresholds["transition_discovery_log2_fc"])
    validation_fdr = float(thresholds["validation_fdr"])
    detection = float(thresholds["transition_min_target_detection_fraction"])
    consistent = int(thresholds["min_consistent_pairs"])
    retention_floor = float(thresholds["cnt_retention_log2_fc_floor"])
    retention_margin = float(thresholds["external_cnt_retention_log2_margin"])
    retention_fdr = float(thresholds["external_retention_fdr"])
    min_detect_reps = float(
        thresholds["external_min_replicate_detection_fraction"]
    )

    table["fine_dct2_enrichment_pass"] = (
        (table["log2_fc_dct2_vs_dct1"] >= transition_fc)
        & (table["fdr_dct2_vs_dct1"] <= validation_fdr)
        & (table["pct_detected_dct2"] >= detection)
        & (table["n_consistent_dct2_vs_dct1"] >= consistent)
    )
    # log2_fc_dct2_vs_cnt <= 0.5 is equivalent to
    # log2_fc_cnt_vs_dct2 >= -0.5 under the frozen retention floor.
    table["fine_cnt_retention_pass"] = (
        (table["log2_fc_dct2_vs_cnt"] <= -retention_floor)
        & (table["pct_detected_cnt"] >= detection)
    )
    table["segment_cnt_retention_pass"] = (
        (table["log2_fc_cnt_vs_dct"] >= retention_floor)
        & (table["ci_low_cnt_vs_dct"] > retention_margin)
        & (table["fdr_noninferiority_cnt_vs_dct"] <= retention_fdr)
        & (table["pct_replicates_detected_cnt"] >= min_detect_reps)
        & (table["n_consistent_cnt_retained"] >= consistent)
    )
    table["DCT2_CNT_external_validation"] = (
        table["fine_dct2_enrichment_pass"]
        & table["fine_cnt_retention_pass"]
        & table["segment_cnt_retention_pass"]
        & ~table["non_distal_specific_or_broad"]
    )
    return table.reindex(columns=columns)


def build_signatures(
    config_path: str | Path,
    output_dir: str | Path,
    discovery_path: str | Path | None = None,
    validation_path: str | Path | None = None,
    segment_validation_path: str | Path | None = None,
    atlas_path: str | Path | None = None,
) -> dict[str, Any]:
    """Build curated and data-derived signature tables plus a gate report."""
    config_path = Path(config_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(config_path)
    _assert_asdn_consistency(cfg)

    dpath = Path(discovery_path) if discovery_path else None
    vpath = Path(validation_path) if validation_path else None
    spath = Path(segment_validation_path) if segment_validation_path else None
    apath = Path(atlas_path) if atlas_path else None
    discovery = _read_table(dpath, DISCOVERY_COLUMNS)
    validation = _read_table(vpath, VALIDATION_COLUMNS)
    segment_validation = _read_table(spath, SEGMENT_VALIDATION_COLUMNS)
    atlas = _read_table(apath, ATLAS_COLUMNS)

    broad = _broad_expression_table(atlas, cfg)
    classified = _classify_data_driven(
        discovery, validation, segment_validation, broad, cfg
    )
    external_validation = _external_validation_signature(
        validation, segment_validation, broad, cfg
    )
    curated = _curated_membership(cfg)
    atlas_comparators = _atlas_comparator_membership(atlas, cfg)

    broad.to_csv(output_dir / "broad_expression_flags.tsv", sep="\t", index=False)
    classified.to_csv(
        output_dir / "data_driven_signature_gene_table.tsv", sep="\t", index=False
    )
    external_validation.to_csv(
        output_dir / "external_validation_signature_gene_table.tsv",
        sep="\t",
        index=False,
    )
    curated.to_csv(output_dir / "curated_gene_set_membership.tsv", sep="\t", index=False)
    atlas_comparators.to_csv(
        output_dir / "atlas_comparator_gene_sets.tsv", sep="\t", index=False
    )

    membership_rows: list[dict[str, Any]] = []
    if not classified.empty:
        for set_name in ("DCT1_core", "strict_DCT2_peaked", "DCT2_CNT_transition"):
            for gene in classified.loc[classified[set_name], "gene_symbol"]:
                membership_rows.append(
                    {
                        "gene_symbol": gene,
                        "gene_set": set_name,
                        "definition_source": "reference_data",
                        "status": "defined",
                    }
                )
    if not external_validation.empty:
        for gene in external_validation.loc[
            external_validation["DCT2_CNT_external_validation"], "gene_symbol"
        ]:
            membership_rows.append(
                {
                    "gene_symbol": gene,
                    "gene_set": "DCT2_CNT_external_validation",
                    "definition_source": (
                        "GSE150338_fine_and_segment_with_atlas_broad_filter"
                    ),
                    "status": "defined",
                }
            )
    data_membership = pd.DataFrame(
        membership_rows,
        columns=["gene_symbol", "gene_set", "definition_source", "status"],
    )
    data_membership.to_csv(
        output_dir / "data_driven_gene_set_membership.tsv", sep="\t", index=False
    )

    frozen_curated = curated.copy()
    frozen_curated["final_for_testing"] = frozen_curated["gene_set"].eq("ASDN")
    frozen_curated = frozen_curated[frozen_curated["gene_set"].eq("ASDN")].copy()
    frozen_atlas = atlas_comparators[
        ["gene_symbol", "gene_set", "category", "definition_source", "status"]
    ].copy()
    frozen_atlas["final_for_testing"] = True
    frozen_data = data_membership.copy()
    if not frozen_data.empty:
        frozen_data["category"] = np.where(
            frozen_data["gene_set"].eq("DCT2_CNT_external_validation"),
            "external_validation_signature",
            "data_driven_reference",
        )
        frozen_data["analysis_role"] = np.where(
            frozen_data["gene_set"].eq("DCT2_CNT_external_validation"),
            "external_validation",
            "primary_or_supporting_reference",
        )
        frozen_data["final_for_testing"] = True
    frozen_curated["analysis_role"] = "primary_or_supporting_reference"
    frozen_atlas["analysis_role"] = "kidney_comparator"
    frozen = pd.concat(
        [
            frozen_curated[
                [
                    "gene_symbol",
                    "gene_set",
                    "category",
                    "definition_source",
                    "status",
                    "final_for_testing",
                    "analysis_role",
                ]
            ],
            frozen_atlas,
            frozen_data.reindex(
                columns=[
                    "gene_symbol",
                    "gene_set",
                    "category",
                    "definition_source",
                    "status",
                    "final_for_testing",
                    "analysis_role",
                ]
            ),
        ],
        ignore_index=True,
    )
    frozen = frozen.sort_values(["gene_set", "gene_symbol"]).reset_index(drop=True)
    frozen.to_csv(output_dir / "frozen_gene_sets.tsv", sep="\t", index=False)

    gate = {
        "freeze_id": cfg.get("freeze_id"),
        "flight_blind": True,
        "discovery_present": discovery is not None,
        "fine_subtype_validation_present": validation is not None,
        "independent_segment_validation_present": segment_validation is not None,
        "whole_kidney_specificity_present": atlas is not None,
        "transition_set_final": bool(
            discovery is not None
            and segment_validation is not None
            and atlas is not None
        ),
        "fine_subtype_sets_final": bool(
            discovery is not None and validation is not None and atlas is not None
        ),
        "external_validation_set_final": bool(
            validation is not None
            and segment_validation is not None
            and atlas is not None
        ),
        "n_members": {
            name: int(len(data_membership[data_membership["gene_set"] == name]))
            for name in (
                "DCT1_core",
                "strict_DCT2_peaked",
                "DCT2_CNT_transition",
                "DCT2_CNT_external_validation",
            )
        },
        "n_atlas_comparator_members": {
            name: int(len(atlas_comparators[atlas_comparators["gene_set"] == name]))
            for name in cfg["reference_builder"].get("atlas_comparator_targets", {})
        },
        "missing_requirements": [
            label
            for label, present in (
                ("GSE228367 count-pseudobulk DCT1-vs-DCT2 table", discovery is not None),
                ("fine-subtype DCT1/DCT2/CNT validation table", validation is not None),
                (
                    "GSE150338 segment-level transition validation table",
                    segment_validation is not None,
                ),
                ("whole-kidney atlas expression summary", atlas is not None),
            )
            if not present
        ],
    }
    (output_dir / "signature_gate.json").write_text(json.dumps(gate, indent=2) + "\n")

    inputs = [
        path
        for path in (dpath, vpath, spath, apath)
        if path is not None and path.exists()
    ]
    provenance = {
        "config": str(config_path),
        "config_sha256": sha256_file(config_path),
        "git_commit": _git_commit(),
        "flight_result_inputs_used": [],
        "reference_inputs": [
            {"path": str(path), "sha256": sha256_file(path)} for path in inputs
        ],
        "rules": cfg["reference_builder"]["thresholds"],
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    return gate
