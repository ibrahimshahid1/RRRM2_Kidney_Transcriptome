#!/usr/bin/env python3
"""Run the Cross-OSDR Network-Contrast Framework.

This script is intentionally separate from the legacy phase runner and is
invoked from ``src/run_all_phases.py`` only when the new opt-in flag is used.
It implements the RRRM-2 stability gate and within-cohort decomposition path,
with explicit external-axis and decision artifacts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd
import yaml

SCRIPT_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_REPO_ROOT))

from src.common import REPO_ROOT, find_sample_col, normalize_labels, load_id_map
from src.networks.contrast_vectors import (
    build_aging_vectors,
    compute_decomposition,
    empirical_pvalue,
    permutation_decomposition,
    precision_weights,
    stratified_bootstrap_indices,
    summarize_bootstrap_decomposition,
)
from src.networks.lioness import pearson_from_sums
from src.networks.stability_test import (
    StabilityDecision,
    StabilityReport,
    apply_stability_gate,
    assert_within_projection_allowed,
    gates_from_config,
    load_stability_decision,
    on_fail_from_config,
    resolution_role_from_config,
    write_stability_artifacts,
)
from src.networks.external_aging_axis import (
    bootstrap_external_projection,
    load_or_build_aging_axis,
    load_gene_feature_matrix,
    load_tms_donor_axis_source,
    build_tms_kidney_female_aging_axis,
)
from src.networks.cross_osdr_projection import bootstrap_cosine_alignment, permutation_cosine_alignment
from src.networks.manuscript_decision import write_manuscript_decision


ARMS = ("ISS-T", "LAR")
AGE_OLD = "OLD"
AGE_YOUNG = "YNG"
ENV_GC = "GC"
ENV_FLT = "FLT"


def log(message: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"{ts} [CVF] {message}")


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
    except Exception:
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def load_config(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def latest_existing_dir(pattern: str) -> Path | None:
    matches = sorted(REPO_ROOT.glob(pattern), reverse=True)
    for match in matches:
        if match.is_dir():
            return match
    return None


def existing_dirs(pattern: str) -> list[Path]:
    """Return matching directories newest-first."""
    return [p for p in sorted(REPO_ROOT.glob(pattern), reverse=True) if p.is_dir()]


def resolve_phase2_dir(requested: str | None, results_dir: Path) -> Path | None:
    candidates: list[Path] = []
    if requested:
        candidates.append(Path(requested))
    candidates.append(results_dir / "networks")
    candidates.extend(existing_dirs("data/results/run_*/networks"))
    candidates.extend(existing_dirs("data/processed/networks/run_*"))
    for candidate in candidates:
        candidate = candidate if candidate.is_absolute() else REPO_ROOT / candidate
        if (candidate / "skeleton_edges.tsv").exists():
            return candidate
    return None


def resolve_wgcna_dir(requested: str | None, results_dir: Path) -> Path | None:
    candidates: list[Path] = []
    if requested:
        candidates.append(Path(requested))
    candidates.extend([
        results_dir / "wgcna",
        latest_existing_dir("data/results/run_*/wgcna") or Path("__missing__"),
    ])
    for candidate in candidates:
        candidate = candidate if candidate.is_absolute() else REPO_ROOT / candidate
        if (candidate / "module_eigengenes.csv").exists():
            return candidate
    return None


def load_rrrm2_inputs(rtech_path: Path, meta_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Return residualized expression (genes x samples), normalized metadata, sample column."""
    rtech = pd.read_csv(rtech_path, sep="\t", compression="gzip", index_col=0)
    rtech.index = rtech.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    meta = pd.read_csv(meta_path, sep="\t", compression="gzip")
    meta = normalize_labels(meta)
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    common = [s for s in rtech.columns if s in meta.index]
    if not common:
        raise ValueError("No overlapping samples between Rtech and metadata")
    return rtech[common], meta.loc[common].copy(), sample_col


def edge_correlations(rtech: pd.DataFrame, samples: Sequence[str], edges: pd.DataFrame) -> np.ndarray:
    """Vectorized Pearson correlations for skeleton edges in one sample set."""
    if len(samples) < 3:
        return np.zeros(len(edges), dtype=float)
    edge_genes = pd.Index(pd.unique(edges[["gene_i", "gene_j"]].to_numpy().ravel()))
    missing = sorted(set(edge_genes) - set(rtech.index))
    if missing:
        raise ValueError(f"Skeleton references {len(missing)} genes missing from Rtech; first={missing[:5]}")
    gene_to_pos = {g: i for i, g in enumerate(edge_genes)}
    ii = edges["gene_i"].map(gene_to_pos).to_numpy(dtype=np.int32)
    jj = edges["gene_j"].map(gene_to_pos).to_numpy(dtype=np.int32)
    x = rtech.loc[edge_genes, list(samples)].to_numpy(dtype=float)
    n = x.shape[1]
    sx = x.sum(axis=1)
    sxx = (x**2).sum(axis=1)
    sxy = (x[ii, :] * x[jj, :]).sum(axis=1)
    return pearson_from_sums(n, sx[ii], sx[jj], sxx[ii], sxx[jj], sxy)


def edge_age_delta(
    rtech: pd.DataFrame,
    meta_rows: pd.DataFrame,
    sample_col: str,
    edges: pd.DataFrame,
    env: str,
    age_values: Sequence[str] | None = None,
) -> np.ndarray:
    """Per-edge Old-Young Fisher-z correlation delta for one EnvGroup."""
    m = meta_rows.copy()
    age = np.asarray(age_values) if age_values is not None else m["Age"].astype(str).to_numpy()
    if len(age) != len(m):
        raise ValueError("age_values length must match metadata rows")
    env_values = m["EnvGroup"].astype(str).to_numpy()
    old_samples = m.loc[(env_values == env) & (age == AGE_OLD), sample_col].tolist()
    young_samples = m.loc[(env_values == env) & (age == AGE_YOUNG), sample_col].tolist()
    if len(old_samples) < 3 or len(young_samples) < 3:
        raise ValueError(
            f"Need at least 3 old and 3 young samples for edge correlations in env={env}; "
            f"found old={len(old_samples)}, young={len(young_samples)}"
        )
    r_old = np.clip(edge_correlations(rtech, old_samples, edges), -0.9995, 0.9995)
    r_young = np.clip(edge_correlations(rtech, young_samples, edges), -0.9995, 0.9995)
    return np.arctanh(r_old) - np.arctanh(r_young)


def aggregate_edge_delta_to_genes(edge_delta: np.ndarray, edges: pd.DataFrame, genes: Sequence[str]) -> np.ndarray:
    """Aggregate incident |delta-z| edge changes to a gene-level node vector."""
    pos = {g: i for i, g in enumerate(genes)}
    out = np.zeros(len(genes), dtype=float)
    counts = np.zeros(len(genes), dtype=float)
    for value, gi, gj in zip(np.abs(edge_delta), edges["gene_i"], edges["gene_j"]):
        if gi in pos:
            out[pos[gi]] += value
            counts[pos[gi]] += 1.0
        if gj in pos:
            out[pos[gj]] += value
            counts[pos[gj]] += 1.0
    return np.divide(out, counts, out=np.zeros_like(out), where=counts > 0)


def load_pathway_scores(
    rtech: pd.DataFrame,
    gene_sets_path: Path,
    id_map_path: Path,
    min_genes: int = 3,
) -> tuple[pd.DataFrame, list[str]]:
    """Build sample x pathway matrix from curated YAML gene sets."""
    with open(gene_sets_path) as fh:
        raw = yaml.safe_load(fh) or {}
    id_map = load_id_map(id_map_path)
    symbol_to_ids: dict[str, set[str]] = {}
    for row in id_map.itertuples():
        symbol_to_ids.setdefault(str(row.mgi_symbol).casefold(), set()).add(str(row.ensembl_gene_id))
    expression_genes = set(rtech.index)

    scores: dict[str, pd.Series] = {}
    for name, values in raw.items():
        if not isinstance(values, dict) or "genes" not in values:
            continue
        flat: list[str] = []
        for item in values.get("genes", []):
            if isinstance(item, str):
                flat.append(item)
            elif isinstance(item, list):
                flat.extend(str(x) for x in item)
            elif isinstance(item, dict):
                for sub in item.values():
                    if isinstance(sub, list):
                        flat.extend(str(x) for x in sub)
        genes: set[str] = set()
        for item in flat:
            item = item.strip()
            if item in expression_genes:
                genes.add(item)
            genes |= symbol_to_ids.get(item.casefold(), set()) & expression_genes
        if len(genes) >= min_genes:
            scores[name] = rtech.loc[sorted(genes)].mean(axis=0)
    if not scores:
        raise ValueError(f"No curated pathway sets reached min_genes={min_genes}")
    df = pd.DataFrame(scores)
    return df, list(df.columns)


def load_module_scores(wgcna_dir: Path, samples: Sequence[str]) -> tuple[pd.DataFrame, list[str]]:
    eig = pd.read_csv(wgcna_dir / "module_eigengenes.csv", index_col=0)
    common = [s for s in samples if s in eig.index]
    if not common:
        raise ValueError(f"No overlapping samples in {wgcna_dir / 'module_eigengenes.csv'}")
    eig = eig.loc[common]
    return eig, list(eig.columns)


def resolve_lioness_weights_path(phase2_dir: Path) -> Path:
    """Return the preferred LIONESS sample-specific edge-weight matrix."""
    for name in ("lioness_edges.npy", "lioness_z_edges.npy"):
        candidate = phase2_dir / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No LIONESS edge-weight matrix found in {phase2_dir}; expected lioness_edges.npy "
        "or legacy lioness_z_edges.npy"
    )


def load_lioness_inputs(phase2_dir: Path, meta: pd.DataFrame, sample_col: str) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load LIONESS edge features as a sample x edge matrix aligned to metadata."""
    edge_path = phase2_dir / "edge_index.tsv"
    if not edge_path.exists():
        edge_path = phase2_dir / "skeleton_edges.tsv"
    if not edge_path.exists():
        raise FileNotFoundError(f"No edge_index.tsv or skeleton_edges.tsv found in {phase2_dir}")
    samples_path = phase2_dir / "lioness_samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(f"No lioness_samples.txt found in {phase2_dir}")

    edges = pd.read_csv(edge_path, sep="\t")
    required = {"gene_i", "gene_j"}
    missing = required - set(edges.columns)
    if missing:
        raise ValueError(f"{edge_path} missing required columns: {sorted(missing)}")
    samples = [line.strip() for line in samples_path.read_text().splitlines() if line.strip()]
    weights_path = resolve_lioness_weights_path(phase2_dir)
    weights = np.load(weights_path, mmap_mode="r")
    if weights.shape[0] != len(samples) or weights.shape[1] != len(edges):
        raise ValueError(
            f"LIONESS matrix shape {weights.shape} does not match samples={len(samples)} edges={len(edges)}"
        )
    common = [s for s in samples if s in set(meta[sample_col].astype(str))]
    if not common:
        raise ValueError("No overlap between LIONESS samples and metadata samples")
    sample_pos = {s: i for i, s in enumerate(samples)}
    row_idx = [sample_pos[s] for s in common]
    feature_names = [f"{r.gene_i}__{r.gene_j}" for r in edges.itertuples()]
    edge_features = pd.DataFrame(
        np.asarray(weights[row_idx, :], dtype=np.float32),
        index=common,
        columns=feature_names,
    )
    return edges, edge_features, common


def aggregate_lioness_edges_to_nodes(
    edge_features: pd.DataFrame,
    edges: pd.DataFrame,
    genes: Sequence[str],
) -> pd.DataFrame:
    """Aggregate sample-specific LIONESS edges to mean absolute incident node strength."""
    genes = list(dict.fromkeys(str(g) for g in genes))
    pos = {g: i for i, g in enumerate(genes)}
    values = np.abs(edge_features.to_numpy(dtype=np.float32, copy=False))
    out = np.zeros((values.shape[0], len(genes)), dtype=np.float32)
    counts = np.zeros(len(genes), dtype=np.float32)
    gi = edges["gene_i"].astype(str).to_numpy()
    gj = edges["gene_j"].astype(str).to_numpy()
    for e, (a, b) in enumerate(zip(gi, gj)):
        ia = pos.get(a)
        if ia is not None:
            out[:, ia] += values[:, e]
            counts[ia] += 1.0
        ib = pos.get(b)
        if ib is not None:
            out[:, ib] += values[:, e]
            counts[ib] += 1.0
    out = np.divide(out, counts, out=np.zeros_like(out), where=counts > 0)
    return pd.DataFrame(out, index=edge_features.index, columns=genes)


def sample_level_age_delta(
    feature_matrix: pd.DataFrame,
    meta_rows: pd.DataFrame,
    sample_col: str,
    env: str,
    age_values: Sequence[str] | None = None,
) -> np.ndarray:
    m = meta_rows.copy()
    age = np.asarray(age_values) if age_values is not None else m["Age"].astype(str).to_numpy()
    env_values = m["EnvGroup"].astype(str).to_numpy()
    old_samples = m.loc[(env_values == env) & (age == AGE_OLD), sample_col].tolist()
    young_samples = m.loc[(env_values == env) & (age == AGE_YOUNG), sample_col].tolist()
    if not old_samples or not young_samples:
        raise ValueError(f"Missing Old/Young cells for env={env}")
    return (
        feature_matrix.loc[old_samples].mean(axis=0).to_numpy(dtype=float)
        - feature_matrix.loc[young_samples].mean(axis=0).to_numpy(dtype=float)
    )


@dataclass
class ResolutionAdapter:
    name: str
    feature_names: list[str]
    vector_builder: Callable[[str, np.ndarray | None, Sequence[str] | None], tuple[np.ndarray, np.ndarray]]
    agc_builder: Callable[[str, np.ndarray | None], np.ndarray]
    stability_strata: Callable[[str], pd.Series]
    decomposition_strata: Callable[[str], pd.Series]
    age_permutation: Callable[[str, np.random.Generator], np.ndarray]


def make_sample_adapter(
    name: str,
    feature_matrix: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
) -> ResolutionAdapter:
    feature_matrix = feature_matrix.copy()

    def arm_rows(arm: str, envs: set[str]) -> pd.DataFrame:
        return meta[(meta["Arm"].astype(str) == arm) & (meta["EnvGroup"].astype(str).isin(envs))].copy()

    def vector_builder(arm: str, indices: np.ndarray | None, age_values: Sequence[str] | None) -> tuple[np.ndarray, np.ndarray]:
        rows = arm_rows(arm, {ENV_GC, ENV_FLT})
        if indices is not None:
            rows = rows.iloc[indices].copy()
        return (
            sample_level_age_delta(feature_matrix, rows, sample_col, ENV_GC, age_values=age_values),
            sample_level_age_delta(feature_matrix, rows, sample_col, ENV_FLT, age_values=age_values),
        )

    def agc_builder(arm: str, indices: np.ndarray | None) -> np.ndarray:
        rows = arm_rows(arm, {ENV_GC})
        if indices is not None:
            rows = rows.iloc[indices].copy()
        return sample_level_age_delta(feature_matrix, rows, sample_col, ENV_GC)

    def stability_strata(arm: str) -> pd.Series:
        rows = arm_rows(arm, {ENV_GC})
        return rows[["Age", "Arm", "EnvGroup"]].astype(str).agg("|".join, axis=1).reset_index(drop=True)

    def decomposition_strata(arm: str) -> pd.Series:
        rows = arm_rows(arm, {ENV_GC, ENV_FLT})
        return rows[["Age", "Arm", "EnvGroup"]].astype(str).agg("|".join, axis=1).reset_index(drop=True)

    def age_permutation(arm: str, rng: np.random.Generator) -> np.ndarray:
        rows = arm_rows(arm, {ENV_GC, ENV_FLT})
        strata = rows[["Arm", "EnvGroup"]].astype(str).agg("|".join, axis=1)
        return stratified_permute(rows["Age"].astype(str).to_numpy(), strata, rng)

    return ResolutionAdapter(
        name=name,
        feature_names=list(feature_matrix.columns.astype(str)),
        vector_builder=vector_builder,
        agc_builder=agc_builder,
        stability_strata=stability_strata,
        decomposition_strata=decomposition_strata,
        age_permutation=age_permutation,
    )


def make_edge_or_gene_adapter(
    name: str,
    rtech: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
    edges: pd.DataFrame,
    genes: list[str],
) -> ResolutionAdapter:
    def arm_rows(arm: str, envs: set[str]) -> pd.DataFrame:
        return meta[(meta["Arm"].astype(str) == arm) & (meta["EnvGroup"].astype(str).isin(envs))].copy()

    def transform(edge_delta: np.ndarray) -> np.ndarray:
        if name == "edge":
            return edge_delta
        return aggregate_edge_delta_to_genes(edge_delta, edges, genes)

    def vector_builder(arm: str, indices: np.ndarray | None, age_values: Sequence[str] | None) -> tuple[np.ndarray, np.ndarray]:
        rows = arm_rows(arm, {ENV_GC, ENV_FLT})
        if indices is not None:
            rows = rows.iloc[indices].copy()
        return (
            transform(edge_age_delta(rtech, rows, sample_col, edges, ENV_GC, age_values=age_values)),
            transform(edge_age_delta(rtech, rows, sample_col, edges, ENV_FLT, age_values=age_values)),
        )

    def agc_builder(arm: str, indices: np.ndarray | None) -> np.ndarray:
        rows = arm_rows(arm, {ENV_GC})
        if indices is not None:
            rows = rows.iloc[indices].copy()
        return transform(edge_age_delta(rtech, rows, sample_col, edges, ENV_GC))

    def stability_strata(arm: str) -> pd.Series:
        rows = arm_rows(arm, {ENV_GC})
        return rows[["Age", "Arm", "EnvGroup"]].astype(str).agg("|".join, axis=1).reset_index(drop=True)

    def decomposition_strata(arm: str) -> pd.Series:
        rows = arm_rows(arm, {ENV_GC, ENV_FLT})
        return rows[["Age", "Arm", "EnvGroup"]].astype(str).agg("|".join, axis=1).reset_index(drop=True)

    def age_permutation(arm: str, rng: np.random.Generator) -> np.ndarray:
        rows = arm_rows(arm, {ENV_GC, ENV_FLT})
        strata = rows[["Arm", "EnvGroup"]].astype(str).agg("|".join, axis=1)
        return stratified_permute(rows["Age"].astype(str).to_numpy(), strata, rng)

    feature_names = (
        [f"{r.gene_i}__{r.gene_j}" for r in edges.itertuples()]
        if name == "edge"
        else list(genes)
    )
    return ResolutionAdapter(
        name=name,
        feature_names=feature_names,
        vector_builder=vector_builder,
        agc_builder=agc_builder,
        stability_strata=stability_strata,
        decomposition_strata=decomposition_strata,
        age_permutation=age_permutation,
    )


def stratified_permute(labels: Sequence[str], strata: Sequence[str], rng: np.random.Generator) -> np.ndarray:
    labels = np.asarray(labels).copy()
    strata = pd.Series(strata).reset_index(drop=True)
    out = labels.copy()
    for value in pd.unique(strata):
        idx = np.where(strata.values == value)[0]
        out[idx] = labels[idx][rng.permutation(len(idx))]
    return out


def build_adapters(args: argparse.Namespace, config: dict, results_dir: Path) -> tuple[dict[str, ResolutionAdapter], dict]:
    rtech_path = Path(args.rtech)
    meta_path = Path(args.meta)
    if not rtech_path.is_absolute():
        rtech_path = REPO_ROOT / rtech_path
    if not meta_path.is_absolute():
        meta_path = REPO_ROOT / meta_path
    rtech, meta, sample_col = load_rrrm2_inputs(rtech_path, meta_path)

    adapters: dict[str, ResolutionAdapter] = {}
    manifest: dict[str, object] = {
        "rtech": str(rtech_path),
        "meta": str(meta_path),
        "rtech_sha256": file_sha256(rtech_path),
        "meta_sha256": file_sha256(meta_path),
    }

    requested = set(args.resolutions)
    if "pathway" in requested:
        pathway_scores, pathway_names = load_pathway_scores(
            rtech,
            REPO_ROOT / args.gene_sets,
            REPO_ROOT / args.id_map,
            min_genes=args.min_pathway_genes,
        )
        adapters["pathway"] = make_sample_adapter("pathway", pathway_scores, meta, sample_col)
        manifest["pathways"] = {"n": len(pathway_names), "source": args.gene_sets}

    if "module" in requested:
        wgcna_dir = resolve_wgcna_dir(args.wgcna_dir, results_dir)
        if wgcna_dir is None:
            raise FileNotFoundError("module resolution requested but no WGCNA module_eigengenes.csv was found")
        module_scores, module_names = load_module_scores(wgcna_dir, rtech.columns)
        adapters["module"] = make_sample_adapter("module", module_scores, meta.loc[module_scores.index], sample_col)
        manifest["wgcna_dir"] = str(wgcna_dir)
        manifest["modules"] = {"n": len(module_names)}

    needs_edges = bool({"edge", "gene", "lioness_node"} & requested)
    if needs_edges:
        phase2_dir = resolve_phase2_dir(args.phase2_dir, results_dir)
        if phase2_dir is None:
            raise FileNotFoundError("edge/gene resolution requested but no skeleton_edges.tsv was found")
        edges, edge_features, lioness_samples = load_lioness_inputs(phase2_dir, meta, sample_col)
        genes_path = phase2_dir / "phase2_genes.txt"
        genes = genes_path.read_text().splitlines() if genes_path.exists() else sorted(pd.unique(edges[["gene_i", "gene_j"]].to_numpy().ravel()))
        genes = [g for g in genes if g in set(edges["gene_i"].astype(str)) or g in set(edges["gene_j"].astype(str))]
        node_features = aggregate_lioness_edges_to_nodes(edge_features, edges, genes)
        manifest["phase2_dir"] = str(phase2_dir)
        manifest["skeleton_edges"] = int(len(edges))
        manifest["node_genes"] = int(len(genes))
        manifest["lioness_samples"] = int(len(lioness_samples))
        manifest["lioness_edge_weights"] = str(resolve_lioness_weights_path(phase2_dir))
        if "edge" in requested:
            adapters["edge"] = make_sample_adapter("edge", edge_features, meta.loc[edge_features.index], sample_col)
        if "gene" in requested:
            adapters["gene"] = make_sample_adapter("gene", node_features, meta.loc[node_features.index], sample_col)
        if "lioness_node" in requested:
            adapters["lioness_node"] = make_sample_adapter("lioness_node", node_features, meta.loc[node_features.index], sample_col)

    return adapters, manifest


def feature_flight_vector(
    feature_matrix: pd.DataFrame,
    sample_rows: pd.DataFrame,
    sample_col: str,
    flt_label: str = ENV_FLT,
    gc_label: str = ENV_GC,
) -> pd.Series:
    """Return FLT-minus-GC mean feature vector from sample-level features."""
    flt_samples = sample_rows.loc[sample_rows["condition"].eq(flt_label), sample_col].tolist()
    gc_samples = sample_rows.loc[sample_rows["condition"].eq(gc_label), sample_col].tolist()
    if len(flt_samples) < 2 or len(gc_samples) < 2:
        raise ValueError(f"Need at least two FLT and GC samples; found FLT={len(flt_samples)}, GC={len(gc_samples)}")
    return feature_matrix.loc[flt_samples].mean(axis=0) - feature_matrix.loc[gc_samples].mean(axis=0)


def rrrm2_arm_flight_reference(
    feature_matrix: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
    arm: str,
) -> pd.Series:
    """Age-pooled RRRM-2 FLT-minus-GC vector for one arm."""
    rows = meta[
        meta["Arm"].astype(str).eq(arm)
        & meta["EnvGroup"].astype(str).isin([ENV_FLT, ENV_GC])
    ].copy()
    rows = rows.assign(condition=rows["EnvGroup"].astype(str))
    return feature_flight_vector(feature_matrix, rows, sample_col)


def load_external_vst(path: Path) -> pd.DataFrame:
    """Load external processed VST matrix as genes x samples."""
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if df.index.duplicated().any():
        df = df.groupby(df.index).mean()
    return df


def resolve_external_file(study_dir: Path, pattern: str) -> Path:
    matches = sorted(study_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {pattern!r} in {study_dir}")
    return matches[0]


def external_sample_rows(study: str, study_dir: Path, vst_columns: Sequence[str]) -> list[tuple[str, pd.DataFrame, str]]:
    """Return named external FLT/GC sample tables for recurrence/context tests."""
    columns = set(vst_columns)
    if study == "OSD-513":
        flt = sorted(c for c in columns if "_FLT_" in c)
        gc = sorted(c for c in columns if "_GC_" in c)
        rows = pd.DataFrame({
            "Sample Name": flt + gc,
            "condition": [ENV_FLT] * len(flt) + [ENV_GC] * len(gc),
        })
        return [("powered_hgc", rows, "directional_recurrence")]

    if study == "OSD-253":
        runsheet = pd.read_csv(resolve_external_file(study_dir, "*runsheet.csv"))
        runsheet = runsheet[runsheet["Sample Name"].isin(columns)].copy()
        strain = runsheet["Factor Value[Strain]"].astype(str).eq("C57BL/6J")
        space = runsheet["Factor Value[Spaceflight]"].astype(str)
        treatment = runsheet["Factor Value[Treatment]"].astype(str).str.casefold()
        flt_mask = strain & space.eq("Space Flight") & treatment.eq("white light")
        original_gc = strain & space.eq("Ground Control") & treatment.eq("blue light")
        rerun_white_gc = strain & space.eq("Ground Control Rerun") & treatment.eq("white light")

        scenarios = []
        for name, gc_mask in [
            ("original_gc_blue_light", original_gc),
            ("rerun_gc_white_light", rerun_white_gc),
        ]:
            flt = sorted(runsheet.loc[flt_mask, "Sample Name"].astype(str))
            gc = sorted(runsheet.loc[gc_mask, "Sample Name"].astype(str))
            rows = pd.DataFrame({
                "Sample Name": flt + gc,
                "condition": [ENV_FLT] * len(flt) + [ENV_GC] * len(gc),
            })
            scenarios.append((name, rows, "context_sensitivity"))
        return scenarios

    raise ValueError(f"Unsupported cross-OSDR study: {study}")


def run_cross_osdr(
    args: argparse.Namespace,
    config: dict,
    outdir: Path,
    n_bootstrap: int,
    n_permutation: int,
    rng: np.random.Generator,
) -> None:
    """Run pathway-level OSD-513/OSD-253 recurrence/context alignment."""
    rtech, meta, sample_col = load_rrrm2_inputs(Path(args.rtech), Path(args.meta))
    rrrm2_pathways, pathway_names = load_pathway_scores(
        rtech,
        REPO_ROOT / args.gene_sets,
        REPO_ROOT / args.id_map,
        min_genes=args.min_pathway_genes,
    )
    references = {
        arm: rrrm2_arm_flight_reference(rrrm2_pathways, meta, sample_col, arm)
        for arm in ARMS
    }
    external_root = Path(args.external_root)
    if not external_root.is_absolute():
        external_root = REPO_ROOT / external_root
    studies = [s.strip() for s in args.cross_osdr_studies.split(",") if s.strip()]
    thresholds = config.get("recurrence", {})
    cross_dir = outdir / "cross_osdr_recurrence"
    cross_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for study in studies:
        study_dir = external_root / study
        if not study_dir.exists():
            log(f"Cross-OSDR skipped missing study folder: {study_dir}")
            continue
        vst = load_external_vst(resolve_external_file(study_dir, "*VST_Counts*.csv"))
        ext_pathways, _ = load_pathway_scores(
            vst,
            REPO_ROOT / args.gene_sets,
            REPO_ROOT / args.id_map,
            min_genes=args.min_pathway_genes,
        )
        for scenario, sample_rows, claim_boundary in external_sample_rows(study, study_dir, ext_pathways.index):
            sample_rows = sample_rows[sample_rows["Sample Name"].isin(ext_pathways.index)].copy()
            if sample_rows["condition"].eq(ENV_FLT).sum() < 2 or sample_rows["condition"].eq(ENV_GC).sum() < 2:
                log(f"Cross-OSDR skipped {study}/{scenario}: insufficient FLT/GC samples")
                continue
            sample_col_ext = "Sample Name"
            strata = sample_rows["condition"].reset_index(drop=True)

            def builder(idx: np.ndarray | None, rows_in: pd.DataFrame = sample_rows) -> pd.Series:
                rows_b = rows_in if idx is None else rows_in.iloc[idx].copy()
                return feature_flight_vector(ext_pathways, rows_b, sample_col_ext)

            def perm_builder(labels: np.ndarray | None, rows_in: pd.DataFrame = sample_rows) -> pd.Series:
                rows_p = rows_in.copy()
                if labels is not None:
                    if len(labels) != len(rows_p):
                        raise ValueError("Permuted labels length must match sample rows")
                    rows_p["condition"] = np.asarray(labels, dtype=object)
                return feature_flight_vector(ext_pathways, rows_p, sample_col_ext)

            for arm, reference in references.items():
                contrast = f"{study}:{scenario}:{arm}:pathway"
                boot, result = bootstrap_cosine_alignment(
                    builder,
                    reference,
                    strata=strata,
                    n_iterations=n_bootstrap,
                    rng=rng,
                    contrast=contrast,
                    thresholds=thresholds,
                )
                boot_path = cross_dir / f"{study}_{scenario}_{arm}_pathway_bootstrap.tsv".replace("/", "_")
                boot.to_csv(boot_path, sep="\t", index=False)
                perm, perm_result = permutation_cosine_alignment(
                    perm_builder,
                    reference,
                    labels=strata,
                    n_iterations=n_permutation,
                    rng=rng,
                    contrast=contrast,
                )
                perm_path = cross_dir / f"{study}_{scenario}_{arm}_pathway_permutation.tsv".replace("/", "_")
                perm.to_csv(perm_path, sep="\t", index=False)
                row = result.as_dict()
                perm_row = perm_result.as_dict()
                perm_row.pop("contrast", None)
                perm_row.pop("point_estimate", None)
                row.update({f"permutation_{k}": v for k, v in perm_row.items()})
                row.update({
                    "study": study,
                    "scenario": scenario,
                    "arm": arm,
                    "resolution": "pathway",
                    "claim_boundary": claim_boundary,
                    "claimable_recurrence_pass": bool(row["recurrence_pass"] and claim_boundary == "directional_recurrence"),
                    "n_flt": int(sample_rows["condition"].eq(ENV_FLT).sum()),
                    "n_gc": int(sample_rows["condition"].eq(ENV_GC).sum()),
                    "n_features": int(len(pathway_names)),
                    "bootstrap_path": str(boot_path),
                    "permutation_path": str(perm_path),
                })
                rows.append(row)

    summary = pd.DataFrame(rows)
    summary.to_csv(cross_dir / "cross_osdr_alignment_summary.tsv", sep="\t", index=False)
    manifest = {
        "studies": studies,
        "resolution": "pathway",
        "n_bootstrap": int(n_bootstrap),
        "n_permutation": int(n_permutation),
        "thresholds": thresholds,
        "permutation_null": "External FLT/GC labels are permuted within each study/scenario; RRRM-2 arm reference vector is held fixed.",
        "claim_rule": "Only directional_recurrence rows contribute to cross_osdr_recurrence; OSD-253 rows are context_sensitivity.",
    }
    (cross_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def run_stability(
    adapters: dict[str, ResolutionAdapter],
    config: dict,
    outdir: Path,
    n_bootstrap: int,
    rng: np.random.Generator,
    manifest: dict,
) -> None:
    reports: list[StabilityReport] = []
    decisions: list[StabilityDecision] = []
    gates = gates_from_config(config)
    for arm in ARMS:
        for resolution, adapter in adapters.items():
            log(f"Stability gate: arm={arm}, resolution={resolution}, B={n_bootstrap}")
            strata = adapter.stability_strata(arm)

            def builder(idx: np.ndarray | None, arm: str = arm, adapter: ResolutionAdapter = adapter) -> np.ndarray:
                return adapter.agc_builder(arm, idx)

            from src.networks.stability_test import estimate_agc_stability

            report = estimate_agc_stability(
                builder,
                strata,
                arm=arm,
                resolution=resolution,
                n_iterations=n_bootstrap,
                rng=rng,
            )
            decision = apply_stability_gate(
                report,
                gates=gates,
                role=resolution_role_from_config(config, resolution),
                on_fail_action=on_fail_from_config(config, resolution),
            )
            reports.append(report)
            decisions.append(decision)
    write_stability_artifacts(reports, decisions, outdir, manifest=manifest)


def run_decomposition(
    adapters: dict[str, ResolutionAdapter],
    config: dict,
    outdir: Path,
    n_bootstrap: int,
    n_permutation: int,
    rng: np.random.Generator,
    bypass_stability: bool,
) -> None:
    decision_path = outdir / "agc_stability_decision.json"
    decision = load_stability_decision(decision_path)
    assert_within_projection_allowed(decision, bypass_stability=bypass_stability)

    decisions = decision.get("decisions", [])
    pass_map = {
        (row["arm"], row["resolution"]): bool(row["passed"])
        for row in decisions
    }
    category_config = config.get("interpretation_categories", {})
    shrink_cfg = config.get("shrinkage", {})
    floor_pct = float(shrink_cfg.get("variance_floor_percentile", 10.0))

    category_rows: list[dict[str, object]] = []
    for arm in ARMS:
        for resolution, adapter in adapters.items():
            if not pass_map.get((arm, resolution), False) and not bypass_stability:
                log(f"Skipping decomposition after failed gate: arm={arm}, resolution={resolution}")
                continue
            log(f"Decomposition: arm={arm}, resolution={resolution}, B={n_bootstrap}, K={n_permutation}")
            strata = adapter.decomposition_strata(arm)

            def vector_builder(idx: np.ndarray | None, age_override: Sequence[str] | None = None,
                               arm: str = arm, adapter: ResolutionAdapter = adapter) -> tuple[np.ndarray, np.ndarray]:
                return adapter.vector_builder(arm, idx, age_override)

            point_gc, point_flt = vector_builder(None)
            point = compute_decomposition(point_flt, point_gc)

            gc_vectors = []
            gc_strata = adapter.stability_strata(arm)
            for _ in range(n_bootstrap):
                gc_idx = stratified_bootstrap_indices(gc_strata, rng)
                gc_vectors.append(adapter.agc_builder(arm, gc_idx))
            weights = precision_weights(np.vstack(gc_vectors), floor_percentile=floor_pct)
            point_w = compute_decomposition(point_flt, point_gc, weights=weights)

            boot_rows: list[dict[str, float | int | bool]] = []
            boot_rows_w: list[dict[str, float | int | bool]] = []
            for b in range(n_bootstrap):
                idx = stratified_bootstrap_indices(strata, rng)
                a_gc, a_flt = vector_builder(idx)
                dec = compute_decomposition(a_flt, a_gc)
                dec_w = compute_decomposition(a_flt, a_gc, weights=weights)
                boot_rows.append({"iteration": b, **dec.as_dict()})
                boot_rows_w.append({"iteration": b, **dec_w.as_dict()})
            boot = pd.DataFrame(boot_rows)
            boot_w = pd.DataFrame(boot_rows_w)

            perm_rows: list[dict[str, float | int]] = []
            perm_rows_w: list[dict[str, float | int]] = []
            for k in range(n_permutation):
                age_perm = adapter.age_permutation(arm, rng)
                a_gc, a_flt = vector_builder(None, age_perm)
                dec = compute_decomposition(a_flt, a_gc)
                dec_w = compute_decomposition(a_flt, a_gc, weights=weights)
                perm_rows.append({"iteration": k, "beta": dec.beta, "cos": dec.cos, "rho": dec.rho})
                perm_rows_w.append({"iteration": k, "beta": dec_w.beta, "cos": dec_w.cos, "rho": dec_w.rho})
            perm = pd.DataFrame(perm_rows)
            perm_w = pd.DataFrame(perm_rows_w)

            summary = summarize_bootstrap_decomposition(
                point,
                boot,
                permutation=perm,
                alpha=float(config.get("bootstrap", {}).get("ci_alpha", 0.05)),
                category_config=category_config,
            )
            summary_w = summarize_bootstrap_decomposition(
                point_w,
                boot_w,
                permutation=perm_w,
                alpha=float(config.get("bootstrap", {}).get("ci_alpha", 0.05)),
                category_config=category_config,
            )
            summary_all = pd.concat([summary, summary_w], ignore_index=True)
            for stat in ("beta", "cos", "rho"):
                summary_all.loc[summary_all["statistic"] == stat, "empirical_p"] = [
                    empirical_pvalue(
                        float(point_w.as_dict()[stat] if weighted else point.as_dict()[stat]),
                        (perm_w if weighted else perm)[stat],
                    )
                    for weighted in summary_all.loc[summary_all["statistic"] == stat, "weighted"]
                ]

            arm_dir = outdir / arm / resolution
            arm_dir.mkdir(parents=True, exist_ok=True)
            boot.to_csv(arm_dir / "bootstrap.tsv", sep="\t", index=False)
            boot_w.to_csv(arm_dir / "bootstrap_weighted.tsv", sep="\t", index=False)
            perm.to_csv(arm_dir / "permutation.tsv", sep="\t", index=False)
            perm_w.to_csv(arm_dir / "permutation_weighted.tsv", sep="\t", index=False)
            summary_all.to_csv(arm_dir / "bootstrap_decomposition.tsv", sep="\t", index=False)
            with open(arm_dir / "point_estimates.json", "w") as fh:
                json.dump({"unweighted": point.as_dict(), "weighted": point_w.as_dict()}, fh, indent=2)

            beta_row = summary_w[summary_w["statistic"] == "beta"].iloc[0].to_dict()
            category_rows.append({
                "arm": arm,
                "resolution": resolution,
                "weighted": True,
                "interpretation": beta_row["interpretation"],
                "frac_amplify": beta_row["frac_amplify"],
                "frac_dampen": beta_row["frac_dampen"],
                "frac_reverse": beta_row["frac_reverse"],
                "frac_redirect": beta_row["frac_redirect"],
            })
    if category_rows:
        pd.DataFrame(category_rows).to_csv(outdir / "interpretation_categories.tsv", sep="\t", index=False)


def run_external_axis(
    args: argparse.Namespace,
    config: dict,
    outdir: Path,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> None:
    axis_cfg = load_config(REPO_ROOT / args.aging_config)
    tms_cfg = axis_cfg.get("tms_kidney_female", {})
    axis_path = Path(args.external_aging_axis or tms_cfg.get("path", ""))
    if not axis_path.is_absolute():
        axis_path = REPO_ROOT / axis_path
    gene_col = tms_cfg.get("columns", {}).get("gene_id", "ensembl_gene_id")
    effect_col = tms_cfg.get("columns", {}).get("effect", "old_vs_young_logfc")
    source_h5ad = tms_cfg.get("source_h5ad", "")
    source_h5ad_path = Path(source_h5ad) if source_h5ad else None
    if source_h5ad_path is not None and not source_h5ad_path.is_absolute():
        source_h5ad_path = REPO_ROOT / source_h5ad_path
    young_months = tms_cfg.get("age_contrast", {}).get("young_months", [3])
    old_months = tms_cfg.get("age_contrast", {}).get("old_months", [18])
    if not axis_path.exists():
        log(f"External aging-axis TSV missing; building from local TMS H5AD: {source_h5ad_path or 'auto-discover'}")
    axis = load_or_build_aging_axis(
        axis_path,
        gene_col,
        effect_col,
        source_h5ad=source_h5ad_path,
        young_months=young_months,
        old_months=old_months,
    )
    tms_donor_expr = None
    tms_donor_meta = None
    if bool(axis_cfg.get("projection", {}).get("bootstrap_with_tms_samples", False)):
        try:
            tms_donor_expr, tms_donor_meta = load_tms_donor_axis_source(axis_path)
        except FileNotFoundError:
            if source_h5ad_path is None:
                raise
            log("TMS donor bootstrap files missing; rebuilding aging axis from H5AD")
            build_tms_kidney_female_aging_axis(
                source_h5ad_path,
                axis_path,
                young_months=young_months,
                old_months=old_months,
            )
            axis = load_or_build_aging_axis(
                axis_path,
                gene_col,
                effect_col,
                source_h5ad=source_h5ad_path,
                young_months=young_months,
                old_months=old_months,
            )
            tms_donor_expr, tms_donor_meta = load_tms_donor_axis_source(axis_path)
    feature_matrix = load_gene_feature_matrix(args.rtech)
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    ext_dir = outdir / "external_aging_axis"
    ext_dir.mkdir(parents=True, exist_ok=True)
    for arm in ARMS:
        log(f"External aging-axis projection: arm={arm}")
        point, boot, summary = bootstrap_external_projection(
            feature_matrix,
            meta,
            axis,
            arm=arm,
            n_iterations=n_bootstrap,
            rng=rng,
            alpha=float(config.get("bootstrap", {}).get("ci_alpha", 0.05)),
            tms_donor_expr=tms_donor_expr,
            tms_donor_meta=tms_donor_meta,
        )
        arm_dir = ext_dir / arm
        arm_dir.mkdir(parents=True, exist_ok=True)
        boot.to_csv(arm_dir / "bootstrap.tsv", sep="\t", index=False)
        summary.to_csv(arm_dir / f"tms_projection_{arm}.tsv", sep="\t", index=False)
        with open(arm_dir / "point_estimates.json", "w") as fh:
            json.dump(point.as_dict(), fh, indent=2)


def run_mechanism_axis(
    args: argparse.Namespace,
    outdir: Path,
    n_bootstrap: int,
) -> None:
    """Run secondary recurrent remodeling-axis mechanism prioritization."""
    cmd = [
        sys.executable,
        str(REPO_ROOT / "scripts" / "run_mechanism_axis_prioritization.py"),
        f"--run-id={args.run_id}",
        f"--contrast-dir={outdir}",
        f"--outdir={outdir / 'mechanism_axis'}",
        f"--rtech={args.rtech}",
        f"--meta={args.meta}",
        f"--external-root={args.external_root}",
        f"--id-map={args.id_map}",
        f"--gene-sets={args.gene_sets}",
        f"--mechanism-gene-sets={args.mechanism_gene_sets}",
        f"--axis-source={args.axis_source}",
        f"--n-bootstrap={args.mechanism_axis_bootstrap or n_bootstrap}",
        f"--seed={args.seed or 20260101}",
        f"--primary-transport-set={args.primary_transport_set}",
        f"--sensitivity-transport-set={args.sensitivity_transport_set}",
    ]
    if args.phase2_dir:
        cmd.append(f"--phase2-dir={args.phase2_dir}")
    if args.rewiring_dir:
        cmd.append(f"--rewiring-dir={args.rewiring_dir}")
    if args.include_osd253_strain:
        cmd.append("--include-osd253-strain")
    log("Mechanism-axis prioritization")
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise RuntimeError("Mechanism-axis prioritization failed")


def run_decision(outdir: Path) -> None:
    decision_path = outdir / "agc_stability_decision.json"
    stability = load_stability_decision(decision_path) if decision_path.exists() else {"decisions": []}
    any_primary_passed = any(
        row.get("passed") and row.get("resolution") in {"module", "pathway"}
        for row in stability.get("decisions", [])
    )
    within_stable = any_primary_passed and not bool(stability.get("fallback_to_external_axis_only", False))
    summaries = list(outdir.glob("*/*/bootstrap_decomposition.tsv"))
    within_signal = False
    module_signal = False
    nonmodule_signal = False
    any_projection = False
    for path in summaries:
        df = pd.read_csv(path, sep="\t")
        beta = df[(df["statistic"] == "beta") & (df["weighted"] == True)]
        if beta.empty:
            continue
        row = beta.iloc[0]
        excludes_zero = row["ci_low"] > 0 or row["ci_high"] < 0
        excludes_one = row["ci_low"] > 1 or row["ci_high"] < 1
        signal = bool(excludes_zero and excludes_one)
        any_projection = any_projection or signal
        within_signal = within_signal or signal
        if signal and "/module/" in str(path):
            module_signal = True
        elif signal:
            nonmodule_signal = True
    cross_osdr_recurrence = False
    cross_summary_path = outdir / "cross_osdr_recurrence" / "cross_osdr_alignment_summary.tsv"
    if cross_summary_path.exists():
        cross = pd.read_csv(cross_summary_path, sep="\t")
        if "claimable_recurrence_pass" in cross.columns:
            cross_osdr_recurrence = bool(cross["claimable_recurrence_pass"].fillna(False).astype(bool).any())
    external_axis_signal = False
    for path in (outdir / "external_aging_axis").glob("*/tms_projection_*.tsv"):
        df = pd.read_csv(path, sep="\t")
        beta = df[df["statistic"] == "beta"]
        if beta.empty:
            continue
        row = beta.iloc[0]
        external_axis_signal = external_axis_signal or bool(row["ci_low"] > 0 or row["ci_high"] < 0)
    evidence = {
        "within_rrrm2_stability_passed": within_stable,
        "within_rrrm2_any_primary_resolution_passed": any_primary_passed,
        "within_rrrm2_fallback_to_external_axis_only": bool(stability.get("fallback_to_external_axis_only", False)),
        "within_rrrm2_projection_significant": within_signal,
        "cross_osdr_recurrence": cross_osdr_recurrence,
        "external_axis_significant": external_axis_signal,
        "module_activity_only_positive": module_signal and not nonmodule_signal and not external_axis_signal and not cross_osdr_recurrence,
        "any_projection_signal": any_projection,
    }
    decision = write_manuscript_decision(evidence, outdir / "manuscript_decision.json")
    log(f"Decision branch: {decision.branch}")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Cross-OSDR Network-Contrast Framework runner")
    ap.add_argument("--config", default="config/contrast_vector_framework.yaml")
    ap.add_argument("--aging-config", default="config/aging_reference.yaml")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--phase2-dir", default="")
    ap.add_argument("--wgcna-dir", default="")
    ap.add_argument("--id-map", default="data/processed/resources/id_map.tsv")
    ap.add_argument("--gene-sets", default="config/gene_sets.yaml")
    ap.add_argument("--mechanism-gene-sets", default="config/mechanism_gene_sets.yaml")
    ap.add_argument("--external-aging-axis", default="")
    ap.add_argument("--outdir", default="")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--phases", default="stability,decomposition,cross-osdr,external-axis,mechanism-axis,decision")
    ap.add_argument("--resolutions", default="module,pathway,gene,lioness_node")
    ap.add_argument("--external-root", default="data/external/osdr")
    ap.add_argument("--cross-osdr-studies", default="OSD-513,OSD-253")
    ap.add_argument("--min-pathway-genes", type=int, default=3)
    ap.add_argument("--n-bootstrap", type=int, default=None)
    ap.add_argument("--n-permutation", type=int, default=None)
    ap.add_argument("--mechanism-axis-bootstrap", type=int, default=None)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--bypass-stability", action="store_true")
    ap.add_argument("--axis-source", default="rrrm2_iss_osd513", choices=["rrrm2_iss_osd513"])
    ap.add_argument("--include-osd253-strain", action="store_true", default=True)
    ap.add_argument("--primary-transport-set", default="dct_ncc_wnk")
    ap.add_argument("--sensitivity-transport-set", default="ion_transport")
    ap.add_argument("--rewiring-dir", default="")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()
    args.config = str(REPO_ROOT / args.config) if not Path(args.config).is_absolute() else args.config
    args.aging_config = str(REPO_ROOT / args.aging_config) if not Path(args.aging_config).is_absolute() else args.aging_config
    args.rtech = str(REPO_ROOT / args.rtech) if not Path(args.rtech).is_absolute() else args.rtech
    args.meta = str(REPO_ROOT / args.meta) if not Path(args.meta).is_absolute() else args.meta
    args.id_map = str(REPO_ROOT / args.id_map) if not Path(args.id_map).is_absolute() else args.id_map
    args.gene_sets = str(REPO_ROOT / args.gene_sets) if not Path(args.gene_sets).is_absolute() else args.gene_sets
    args.mechanism_gene_sets = str(REPO_ROOT / args.mechanism_gene_sets) if not Path(args.mechanism_gene_sets).is_absolute() else args.mechanism_gene_sets
    args.phases = [p.strip() for p in args.phases.split(",") if p.strip()]
    args.resolutions = [r.strip() for r in args.resolutions.split(",") if r.strip()]
    return args


def main() -> int:
    args = parse_args()
    config = load_config(Path(args.config))
    run_id = args.run_id or datetime.now().strftime("run_%Y%m%d_%H%M%S_contrast_vectors")
    results_dir = Path(args.outdir) if args.outdir else REPO_ROOT / "data" / "results" / run_id / "contrast_vectors"
    if not results_dir.is_absolute():
        results_dir = REPO_ROOT / results_dir

    n_bootstrap = int(args.n_bootstrap or config.get("bootstrap", {}).get("n_iterations", 2000))
    n_permutation = int(args.n_permutation or config.get("permutation", {}).get("n_iterations", 5000))
    seed = int(args.seed or config.get("bootstrap", {}).get("random_seed", 20260101))
    rng = np.random.default_rng(seed)

    log(f"Run ID: {run_id}")
    log(f"Output: {results_dir}")
    log(f"Phases: {args.phases}")
    log(f"Resolutions: {args.resolutions}")
    if args.dry_run:
        log("Dry run complete; no files written.")
        return 0

    results_dir.mkdir(parents=True, exist_ok=True)
    adapter_phases = {"stability", "decomposition"} & set(args.phases)
    if adapter_phases:
        adapters, manifest_inputs = build_adapters(args, config, results_dir.parent)
    else:
        adapters = {}
        manifest_inputs = {
            "rtech": str(args.rtech),
            "meta": str(args.meta),
            "rtech_sha256": file_sha256(Path(args.rtech)) if Path(args.rtech).exists() else "",
            "meta_sha256": file_sha256(Path(args.meta)) if Path(args.meta).exists() else "",
        }
    manifest = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "git_commit": git_commit(),
        "config": str(args.config),
        "config_sha256": file_sha256(Path(args.config)),
        "parameters": {
            "n_bootstrap": n_bootstrap,
            "n_permutation": n_permutation,
            "seed": seed,
            "resolutions": args.resolutions,
            "phases": args.phases,
            "bypass_stability": args.bypass_stability,
        },
        "inputs": manifest_inputs,
    }
    with open(results_dir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)

    if "stability" in args.phases:
        run_stability(adapters, config, results_dir, n_bootstrap, rng, manifest)
    if "decomposition" in args.phases:
        try:
            run_decomposition(adapters, config, results_dir, n_bootstrap, n_permutation, rng, args.bypass_stability)
        except RuntimeError as exc:
            external_requested = "external-axis" in args.phases or "external_axis" in args.phases
            if external_requested and "Guardrail A failed" in str(exc):
                log(f"Within-RRRM-2 decomposition skipped: {exc}")
                log("Continuing with external aging-axis projection per Guardrail E.")
            else:
                raise
    if "cross-osdr" in args.phases or "cross_osdr" in args.phases:
        run_cross_osdr(args, config, results_dir, n_bootstrap, n_permutation, rng)
    if "external-axis" in args.phases or "external_axis" in args.phases:
        run_external_axis(args, config, results_dir, rng, n_bootstrap)
    if "mechanism-axis" in args.phases or "mechanism_axis" in args.phases:
        run_mechanism_axis(args, results_dir, n_bootstrap)
    if "decision" in args.phases:
        run_decision(results_dir)

    log("Contrast-vector framework completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
