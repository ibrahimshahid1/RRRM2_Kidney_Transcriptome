#!/usr/bin/env python3
"""Reconstruct GSE228367 subtype pseudobulks from official raw 10x matrices.

The author-separated RDS objects supply cell membership only. This helper
matches those cells to the NK1--NK3 filtered 10x H5 matrices in the official
GEO raw archive, verifies complete barcode coverage and integer counts, then
emits gene-symbol-collapsed count and detection pseudobulks for edgeR.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
import tempfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _decode(values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            value.decode("utf-8") if isinstance(value, (bytes, np.bytes_)) else str(value)
            for value in values
        ],
        dtype=object,
    )


def _read_10x_h5(
    archive: tarfile.TarFile, member: tarfile.TarInfo
) -> tuple[sp.csc_matrix, np.ndarray, np.ndarray, np.ndarray]:
    source = archive.extractfile(member)
    if source is None:
        raise ValueError(f"Could not read archive member {member.name}")
    with tempfile.NamedTemporaryFile(suffix=".h5") as temporary:
        while chunk := source.read(1024 * 1024):
            temporary.write(chunk)
        temporary.flush()
        with h5py.File(temporary.name, "r") as handle:
            matrix = handle["matrix"]
            shape = tuple(int(value) for value in matrix["shape"][:])
            counts = sp.csc_matrix(
                (
                    matrix["data"][:],
                    matrix["indices"][:],
                    matrix["indptr"][:],
                ),
                shape=shape,
            )
            barcodes = _decode(matrix["barcodes"][:])
            genes = _decode(matrix["features"]["name"][:])
            feature_ids = _decode(matrix["features"]["id"][:])
            feature_types = _decode(matrix["features"]["feature_type"][:])
    gene_expression = feature_types == "Gene Expression"
    if not gene_expression.any():
        raise ValueError(f"{member.name} contains no Gene Expression features")
    return (
        counts[gene_expression, :].tocsc(),
        genes[gene_expression],
        feature_ids[gene_expression],
        barcodes,
    )


def _symbol_collapse_matrix(
    genes: np.ndarray,
) -> tuple[sp.csr_matrix, np.ndarray]:
    symbols = pd.Index(genes.astype(str))
    unique_symbols = pd.Index(pd.unique(symbols))
    codes = pd.Categorical(symbols, categories=unique_symbols).codes
    mapping = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.int8), (codes, np.arange(len(codes)))),
        shape=(len(unique_symbols), len(codes)),
    )
    return mapping, unique_symbols.to_numpy(dtype=object)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True)
    parser.add_argument("--membership", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    archive_path = Path(args.archive)
    membership_path = Path(args.membership)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    membership = pd.read_csv(membership_path, sep="\t")
    required = {"cell_id", "raw_barcode", "replicate", "subtype"}
    missing = required - set(membership)
    if missing:
        raise ValueError(f"Membership table lacks columns: {sorted(missing)}")
    membership = membership.copy()
    for column in required:
        membership[column] = membership[column].astype(str)
    if membership.duplicated(["replicate", "raw_barcode"]).any():
        duplicated = membership.loc[
            membership.duplicated(["replicate", "raw_barcode"], keep=False),
            ["replicate", "raw_barcode", "subtype"],
        ]
        raise ValueError(
            "A raw barcode is assigned more than once within a replicate: "
            f"{duplicated.head().to_dict(orient='records')}"
        )

    expected_subtypes = ["DCT1", "DCT2"]
    count_columns: dict[str, np.ndarray] = {}
    detection_columns: dict[str, np.ndarray] = {}
    mapping_rows: list[dict[str, object]] = []
    reference_genes: np.ndarray | None = None
    reference_feature_ids: np.ndarray | None = None
    symbol_mapping: sp.csr_matrix | None = None
    collapsed_genes: np.ndarray | None = None

    with tarfile.open(archive_path, "r") as archive:
        members = archive.getmembers()
        for replicate in sorted(membership["replicate"].unique()):
            suffix = f"_{replicate}_filtered_feature_bc_matrix.h5"
            matches = [member for member in members if member.name.endswith(suffix)]
            if len(matches) != 1:
                raise ValueError(
                    f"Expected one filtered 10x matrix for {replicate}; found "
                    f"{[member.name for member in matches]}"
                )
            matrix, genes, feature_ids, barcodes = _read_10x_h5(
                archive, matches[0]
            )
            if matrix.shape != (len(genes), len(barcodes)):
                raise ValueError(f"Invalid 10x matrix dimensions for {replicate}")
            if matrix.nnz == 0 or not np.all(np.isfinite(matrix.data)):
                raise ValueError(f"Non-finite or empty count matrix for {replicate}")
            if np.max(np.abs(matrix.data - np.round(matrix.data))) > 1e-8:
                raise ValueError(
                    f"Official filtered 10x counts are not integer-like for {replicate}"
                )
            if (matrix.data < 0).any():
                raise ValueError(f"Negative official counts found for {replicate}")

            if reference_genes is None:
                reference_genes = genes
                reference_feature_ids = feature_ids
                symbol_mapping, collapsed_genes = _symbol_collapse_matrix(genes)
            elif not (
                np.array_equal(genes, reference_genes)
                and np.array_equal(feature_ids, reference_feature_ids)
            ):
                raise ValueError(
                    "Gene order/identifiers differ across official NK matrices"
                )

            barcode_index = pd.Index(barcodes)
            for subtype in expected_subtypes:
                selected = membership[
                    membership["replicate"].eq(replicate)
                    & membership["subtype"].eq(subtype)
                ]
                if selected.empty:
                    raise ValueError(f"No {subtype} cells assigned in {replicate}")
                positions = barcode_index.get_indexer(selected["raw_barcode"])
                matched = positions >= 0
                mapping_fraction = float(matched.mean())
                if not matched.all():
                    missing_barcodes = selected.loc[
                        ~matched, ["cell_id", "raw_barcode"]
                    ]
                    raise ValueError(
                        f"Incomplete {replicate}/{subtype} barcode mapping "
                        f"({mapping_fraction:.6f}); examples="
                        f"{missing_barcodes.head().to_dict(orient='records')}"
                    )
                if symbol_mapping is None:
                    raise RuntimeError("Gene-symbol collapse matrix was not initialized")
                selected_counts = symbol_mapping @ matrix[:, positions]
                column = f"{subtype}::{replicate}"
                count_columns[column] = np.asarray(
                    selected_counts.sum(axis=1)
                ).ravel()
                detection_columns[column] = np.asarray(
                    (selected_counts > 0).sum(axis=1)
                ).ravel() / len(positions)
                mapping_rows.append(
                    {
                        "sample_id": column,
                        "replicate": replicate,
                        "subtype": subtype,
                        "archive_member": matches[0].name,
                        "n_cells_membership": int(len(selected)),
                        "n_barcodes_matched": int(matched.sum()),
                        "mapping_fraction": mapping_fraction,
                        "raw_matrix_n_cells": int(matrix.shape[1]),
                        "raw_matrix_n_features": int(matrix.shape[0]),
                        "integer_like_counts": True,
                        "pseudobulk_library_size": int(
                            count_columns[column].sum()
                        ),
                    }
                )

    if reference_genes is None or collapsed_genes is None:
        raise ValueError("No official 10x matrices were processed")
    count_matrix = np.column_stack(
        [count_columns[name] for name in sorted(count_columns)]
    )
    detection_matrix = np.column_stack(
        [detection_columns[name] for name in sorted(detection_columns)]
    )
    if np.max(np.abs(count_matrix - np.round(count_matrix))) > 1e-8:
        raise ValueError("Collapsed pseudobulk counts are not integer-like")

    ordered_columns = sorted(count_columns)
    counts_frame = pd.DataFrame(
        count_matrix.astype(np.int64),
        columns=ordered_columns,
    )
    counts_frame.insert(0, "gene_symbol", collapsed_genes)
    counts_frame.to_csv(
        output_dir / "gse228367_raw_pseudobulk_counts.tsv.gz",
        sep="\t",
        index=False,
    )
    detection_frame = pd.DataFrame(
        detection_matrix,
        columns=ordered_columns,
    )
    detection_frame.insert(0, "gene_symbol", collapsed_genes)
    detection_frame.to_csv(
        output_dir / "gse228367_raw_pseudobulk_detection.tsv.gz",
        sep="\t",
        index=False,
    )
    mapping = pd.DataFrame(mapping_rows).sort_values(["subtype", "replicate"])
    mapping.to_csv(
        output_dir / "gse228367_barcode_mapping_qc.tsv",
        sep="\t",
        index=False,
    )
    qc = {
        "archive": str(archive_path),
        "archive_sha256": _sha256(archive_path),
        "membership": str(membership_path),
        "membership_sha256": _sha256(membership_path),
        "n_membership_cells": int(len(membership)),
        "n_pseudobulks": int(len(mapping)),
        "all_barcodes_matched": bool(mapping["mapping_fraction"].eq(1.0).all()),
        "minimum_mapping_fraction": float(mapping["mapping_fraction"].min()),
        "all_counts_integer_like": bool(mapping["integer_like_counts"].all()),
        "n_gene_expression_features": int(len(reference_genes)),
        "n_unique_gene_symbols": int(len(collapsed_genes)),
        "source_count_matrix": "official filtered_feature_bc_matrix.h5",
        "subtype_membership_source": "author-separated RDS metadata only",
        "flight_result_inputs_used": [],
    }
    (output_dir / "gse228367_raw_mapping_qc.json").write_text(
        json.dumps(qc, indent=2) + "\n"
    )
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
