"""Bootstrap orchestration helpers for contrast-vector decomposition.

This module is the statistics-facing wrapper around
``src.networks.contrast_vectors``. It keeps file/summary conventions in one
place so Phase 3-7 runners emit consistent artifacts.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from src.networks.contrast_vectors import (
    compute_decomposition,
    bootstrap_decomposition,
    permutation_decomposition,
    precision_weights,
    summarize_bootstrap_decomposition,
    stratified_bootstrap_indices,
)


def bootstrap_gc_precision_weights(
    agc_builder: Callable[[np.ndarray], np.ndarray],
    strata: Sequence,
    n_iterations: int,
    rng: np.random.Generator,
    floor_percentile: float = 10.0,
) -> np.ndarray:
    """Estimate Guardrail C precision weights from bootstrapped A_GC vectors."""
    vectors: list[np.ndarray] = []
    for _ in range(int(n_iterations)):
        idx = stratified_bootstrap_indices(strata, rng)
        vectors.append(np.asarray(agc_builder(idx), dtype=float).ravel())
    return precision_weights(np.vstack(vectors), floor_percentile=floor_percentile)


def run_bootstrap_decomposition(
    vector_builder: Callable[[np.ndarray | None], tuple[np.ndarray, np.ndarray]],
    strata: Sequence,
    n_iterations: int,
    rng: np.random.Generator,
    weights: np.ndarray | None = None,
    permutation: pd.DataFrame | None = None,
    alpha: float = 0.05,
    category_config: Mapping[str, float] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute point estimate, bootstrap distribution, and summary table."""
    a_gc, a_flt = vector_builder(None)
    point = compute_decomposition(a_flt, a_gc, weights=weights)

    def builder_nonnull(idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        return vector_builder(idx)

    boot = bootstrap_decomposition(
        builder_nonnull,
        strata,
        n_iterations=n_iterations,
        weights=weights,
        rng=rng,
    )
    summary = summarize_bootstrap_decomposition(
        point,
        boot,
        permutation=permutation,
        alpha=alpha,
        category_config=category_config,
    )
    return boot, summary


def write_decomposition_artifacts(
    outdir: str | Path,
    bootstrap: pd.DataFrame,
    summary: pd.DataFrame,
    permutation: pd.DataFrame | None = None,
    bootstrap_filename: str = "bootstrap.tsv",
    summary_filename: str = "bootstrap_decomposition.tsv",
    permutation_filename: str = "permutation.tsv",
) -> dict[str, Path]:
    """Persist per-arm/per-resolution decomposition artifacts."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {
        "bootstrap": outdir / bootstrap_filename,
        "summary": outdir / summary_filename,
    }
    bootstrap.to_csv(paths["bootstrap"], sep="\t", index=False)
    summary.to_csv(paths["summary"], sep="\t", index=False)
    if permutation is not None:
        paths["permutation"] = outdir / permutation_filename
        permutation.to_csv(paths["permutation"], sep="\t", index=False)
    return paths


__all__ = [
    "bootstrap_gc_precision_weights",
    "permutation_decomposition",
    "run_bootstrap_decomposition",
    "write_decomposition_artifacts",
]
