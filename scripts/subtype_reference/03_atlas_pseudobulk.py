#!/usr/bin/env python3
"""Create source-study whole-kidney atlas pseudobulks and mapping QC."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


def _sha256(path: Path) -> str:
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


def _assign_compartments(labels: pd.Series, aliases: dict[str, list[str]]) -> pd.Series:
    """Assign one compartment using config order; unmatched cells remain NA."""
    output = pd.Series(pd.NA, index=labels.index, dtype="object")
    lower = labels.astype("string").fillna("").str.lower()
    for compartment, patterns in aliases.items():
        available = output.isna()
        if not available.any():
            break
        expression = "|".join(f"(?:{pattern})" for pattern in patterns)
        hit = lower.str.contains(expression, regex=True, flags=re.IGNORECASE, na=False)
        output.loc[available & hit] = compartment
    return output


def _integer_like(matrix: sp.spmatrix) -> bool:
    if matrix.nnz == 0:
        return False
    for start in range(0, matrix.nnz, 1_000_000):
        values = matrix.data[start : start + 1_000_000]
        if not np.all(np.isfinite(values)):
            return False
        if np.max(np.abs(values - np.round(values))) >= 1e-6:
            return False
    return True


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--source-study-column", default="Origin")
    parser.add_argument("--cell-type-column", default="cell_type")
    parser.add_argument("--fallback-cell-type-column")
    parser.add_argument(
        "--reference-label",
        default="whole_kidney_atlas",
        help="Source label written to output provenance.",
    )
    args = parser.parse_args()

    config_path = Path(args.config)
    input_path = Path(args.input)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    cfg = yaml.safe_load(config_path.read_text()) or {}
    rb = cfg["reference_builder"]
    aliases = rb["whole_kidney_compartment_aliases"]
    min_cells = int(rb["thresholds"]["min_cells_per_pseudobulk"])
    min_cpm = float(rb["broad_expression"]["min_cpm"])

    adata = ad.read_h5ad(input_path)
    required_columns = [args.source_study_column, args.cell_type_column]
    if args.fallback_cell_type_column:
        required_columns.append(args.fallback_cell_type_column)
    for column in required_columns:
        if column not in adata.obs:
            raise ValueError(f"{input_path} lacks required obs column {column!r}")
    # The atlas layer is restricted to 3,000 variable genes. Full raw.X is
    # required for broad-expression exclusion and comparator construction.
    if adata.raw is not None:
        counts = adata.raw.X
        count_source = "raw.X"
        genes = pd.Index(adata.raw.var_names.astype(str), name="gene_symbol")
    elif args.counts_layer in adata.layers:
        counts = adata.layers[args.counts_layer]
        count_source = f"layers[{args.counts_layer!r}]"
        genes = pd.Index(adata.var_names.astype(str), name="gene_symbol")
    else:
        counts = adata.X
        count_source = "X"
        genes = pd.Index(adata.var_names.astype(str), name="gene_symbol")
    counts = sp.csr_matrix(counts)
    expected_shape = (adata.n_obs, len(genes))
    if counts.shape != expected_shape:
        raise ValueError(
            f"Selected count matrix has shape {counts.shape}, expected {expected_shape}"
        )
    if not _integer_like(counts):
        raise ValueError("Selected matrix is not integer-like raw counts")

    obs = adata.obs[required_columns].copy()
    obs["compartment"] = _assign_compartments(obs[args.cell_type_column], aliases)
    if args.fallback_cell_type_column:
        fallback = _assign_compartments(
            obs[args.fallback_cell_type_column], aliases
        )
        obs.loc[obs["compartment"].isna(), "compartment"] = fallback
    obs["source_study_id"] = obs[args.source_study_column].astype(str)
    mapped = obs["compartment"].notna()
    mapping_fraction = float(mapped.mean())
    qc_cfg = rb["atlas_qc"]
    minimum_mapping_fraction = float(qc_cfg["minimum_mapping_fraction"])
    if mapping_fraction < minimum_mapping_fraction:
        raise ValueError(
            f"Atlas mapping fraction {mapping_fraction:.4f} is below frozen "
            f"minimum {minimum_mapping_fraction:.4f}"
        )
    group = (
        obs.loc[mapped, "source_study_id"]
        + "::"
        + obs.loc[mapped, "compartment"].astype(str)
    )
    levels = sorted(group.unique())
    codes = pd.Categorical(group, categories=levels).codes
    membership = sp.csr_matrix(
        (np.ones(len(codes)), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(levels)),
    )
    mapped_counts = counts[mapped.to_numpy(), :]
    pseudobulk = membership.T @ mapped_counts
    pseudobulk = sp.csr_matrix(pseudobulk)
    cell_counts = np.asarray(membership.sum(axis=0)).ravel().astype(int)

    sample = pd.DataFrame({"sample_id": levels, "n_cells": cell_counts})
    sample[["source_study_id", "compartment"]] = sample["sample_id"].str.split(
        "::", n=1, expand=True
    )
    keep = sample["n_cells"] >= min_cells
    sample["included"] = keep
    sample.to_csv(out / "atlas_pseudobulk_sample_metadata.tsv", sep="\t", index=False)
    pseudobulk = pseudobulk[keep.to_numpy(), :]
    kept_sample = sample.loc[keep].reset_index(drop=True)

    dense = pseudobulk.toarray().astype(np.float64)
    library = dense.sum(axis=1)
    cpm = np.divide(
        dense * 1e6,
        library[:, None],
        out=np.zeros_like(dense),
        where=library[:, None] > 0,
    )
    count_frame = pd.DataFrame(dense.T.astype(np.int64), index=genes)
    count_frame.columns = kept_sample["sample_id"]
    count_frame.to_csv(out / "atlas_pseudobulk_counts.tsv.gz", sep="\t")

    summary_rows: list[pd.DataFrame] = []
    for compartment, indexes in kept_sample.groupby(
        "compartment", observed=True
    ).groups.items():
        idx = np.asarray(list(indexes), dtype=int)
        summary_rows.append(
            pd.DataFrame(
                {
                    "gene_symbol": genes,
                    "compartment": compartment,
                    "mean_cpm": cpm[idx, :].mean(axis=0),
                    "median_cpm": np.median(cpm[idx, :], axis=0),
                    "source_study_detection_fraction": (
                        cpm[idx, :] >= min_cpm
                    ).mean(axis=0),
                    "n_source_studies": len(idx),
                    "reference_label": args.reference_label,
                }
            )
        )
    summary = pd.concat(summary_rows, ignore_index=True)
    summary.to_csv(out / "atlas_compartment_expression.tsv.gz", sep="\t", index=False)

    studies_by_compartment = (
        kept_sample.groupby("compartment", observed=True)["source_study_id"]
        .nunique()
        .to_dict()
    )
    minimum_studies = int(qc_cfg["minimum_source_studies_per_compartment"])
    required_compartments = list(map(str, qc_cfg["required_compartments"]))
    failed_compartments = {
        compartment: int(studies_by_compartment.get(compartment, 0))
        for compartment in required_compartments
        if int(studies_by_compartment.get(compartment, 0)) < minimum_studies
    }
    if failed_compartments:
        raise ValueError(
            "Atlas compartments failed frozen source-study coverage QC: "
            f"{failed_compartments}; minimum={minimum_studies}"
        )

    mapping_columns = [args.cell_type_column]
    if args.fallback_cell_type_column:
        mapping_columns.append(args.fallback_cell_type_column)
    cell_mapping = (
        obs.groupby(
            mapping_columns + ["compartment"], dropna=False, observed=True
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    cell_mapping.to_csv(out / "atlas_cell_type_mapping_qc.tsv", sep="\t", index=False)

    qc = {
        "reference_label": args.reference_label,
        "n_cells": int(adata.n_obs),
        "n_genes_selected_count_matrix": int(len(genes)),
        "selected_count_matrix_shape": [int(value) for value in counts.shape],
        "adata_shape": [int(adata.n_obs), int(adata.n_vars)],
        "adata_raw_shape": (
            [int(adata.raw.n_obs), int(adata.raw.n_vars)]
            if adata.raw is not None
            else None
        ),
        "n_cells_mapped": int(mapped.sum()),
        "mapping_fraction": mapping_fraction,
        "minimum_mapping_fraction_required": minimum_mapping_fraction,
        "n_pseudobulks_included": int(keep.sum()),
        "source_studies_by_compartment": studies_by_compartment,
        "minimum_source_studies_per_compartment_required": minimum_studies,
        "required_compartments": required_compartments,
        "source_study_column": args.source_study_column,
        "count_source": count_source,
        "integer_like_counts": True,
        "limitations": (
            "This atlas summary supports broad-expression and major-compartment QC. "
            "It does not independently validate DCT1/DCT2/CNT identity unless those "
            "author labels are present in the source metadata. Origin identifies a "
            "source study, not an animal or donor."
        ),
    }
    (out / "atlas_pseudobulk_qc.json").write_text(json.dumps(qc, indent=2) + "\n")
    provenance = {
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "input": str(input_path),
        "input_sha256": _sha256(input_path),
        "git_commit": _git_commit(),
        "flight_result_inputs_used": [],
        "command_parameters": vars(args),
    }
    (out / "atlas_pseudobulk_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
