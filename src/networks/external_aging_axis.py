"""Guardrail E: project RRRM-2 flight vectors onto an external aging axis.

The implementation is deliberately feature-agnostic: callers provide a flight
effect vector and an independently estimated aging-axis vector indexed by the
same feature IDs. The CLI provides the gene-level RRRM-2 path used by the
contrast-vector orchestrator.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind

from src.common import REPO_ROOT, bh_fdr, find_sample_col, normalize_labels
from src.networks.contrast_vectors import (
    Decomposition,
    compute_decomposition,
    stratified_bootstrap_indices,
    summarize_bootstrap_decomposition,
)


def load_aging_axis(path: str | Path, gene_col: str, effect_col: str) -> pd.Series:
    """Load a per-feature external aging axis as a numeric Series."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"External aging axis not found: {path}")
    axis = pd.read_csv(path, sep="\t")
    missing = {gene_col, effect_col} - set(axis.columns)
    if missing:
        raise ValueError(f"{path} missing columns required for aging axis: {sorted(missing)}")
    out = axis[[gene_col, effect_col]].dropna().copy()
    out[gene_col] = out[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    out = out.groupby(gene_col)[effect_col].mean()
    out = out.astype(float)
    if out.empty:
        raise ValueError(f"External aging axis is empty after loading {path}")
    return out


def parse_age_months(value: object) -> float:
    """Parse TMS age labels like ``3m`` / ``18-month-old stage`` to months."""
    text = str(value).strip().lower()
    match = re.search(r"(\d+(?:\.\d+)?)\s*(m|mo|month)", text)
    if not match:
        match = re.search(r"^(\d+(?:\.\d+)?)$", text)
    if not match:
        return float("nan")
    return float(match.group(1))


def discover_tms_h5ad(base_dir: str | Path | None = None) -> Path:
    """Find the local TMS kidney female H5AD used by deconvolution."""
    base = Path(base_dir) if base_dir is not None else REPO_ROOT / "data/external/single_cell_atlases"
    candidates = [
        base / "tms_kidney_female_ALLDATASETS_counts_innerGenes.h5ad",
        base / "tms_kidney_female_ALLDATASETS_innerGenes.h5ad",
        base / "kidney_female_b8c618e5-4b3d-4566-8a3f-7e40047f5c54.h5ad",
        base / "kidney_female_524179b0-b406-4723-9c46-293ffa77ca81.h5ad",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        f"No local TMS kidney female H5AD found under {base}. "
        "Expected one of: " + ", ".join(c.name for c in candidates)
    )


def _axis_from_matrix(
    x,
    obs: pd.DataFrame,
    var_names: Sequence[str],
    donor_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build donor-pseudobulk logCPM expression from an AnnData matrix."""
    donors = sorted(obs[donor_col].astype(str).unique())
    rows: list[np.ndarray] = []
    meta_rows: list[dict[str, object]] = []
    for donor in donors:
        mask = obs[donor_col].astype(str).to_numpy() == donor
        donor_x = x[mask, :]
        if sparse.issparse(donor_x):
            counts = np.asarray(donor_x.sum(axis=0)).ravel()
        else:
            counts = np.asarray(donor_x).sum(axis=0)
        lib_size = float(counts.sum())
        if lib_size <= 0:
            continue
        rows.append(np.log2((counts / lib_size) * 1_000_000.0 + 1.0))
        donor_obs = obs.loc[mask]
        meta_rows.append({
            "donor_id": donor,
            "age_months": float(donor_obs["age_months"].median()),
            "age_group": donor_obs["age_group"].mode().iloc[0],
            "n_cells": int(mask.sum()),
            "library_size": lib_size,
        })
    if not rows:
        raise ValueError("No donor pseudobulk rows could be built from the TMS H5AD")
    expr = pd.DataFrame(np.vstack(rows), index=[r["donor_id"] for r in meta_rows], columns=list(var_names))
    donor_meta = pd.DataFrame(meta_rows).set_index("donor_id")
    return expr, donor_meta


def build_tms_kidney_female_aging_axis(
    h5ad_path: str | Path,
    outpath: str | Path,
    young_months: Sequence[float] = (3,),
    old_months: Sequence[float] = (18,),
    donor_col: str = "donor_id",
    age_col: str = "age",
    sex_col: str = "sex",
    sex_value: str = "female",
) -> pd.DataFrame:
    """Create the TMS kidney female old-vs-young aging-axis TSV.

    Cells are first summed to donor-level pseudobulk counts, then converted to
    log2(CPM + 1). The exported effect is the old-donor mean minus the young-
    donor mean for each gene. Age ranges are inclusive when a two-value range is
    supplied; a single-value list such as ``[3]`` selects exactly 3-month donors.
    """
    try:
        import scanpy as sc
    except ImportError as exc:  # pragma: no cover - environment guard
        raise ImportError("scanpy is required to build the TMS aging axis from H5AD") from exc

    h5ad_path = Path(h5ad_path)
    if not h5ad_path.exists():
        raise FileNotFoundError(f"TMS H5AD not found: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path)
    obs = adata.obs.copy()
    if donor_col not in obs.columns or age_col not in obs.columns:
        raise ValueError(f"{h5ad_path} must contain obs columns {donor_col!r} and {age_col!r}")
    obs["age_months"] = obs[age_col].map(parse_age_months).astype(float)
    young_lo, young_hi = min(young_months), max(young_months)
    old_lo, old_hi = min(old_months), max(old_months)
    young_mask = obs["age_months"].between(young_lo, young_hi, inclusive="both")
    old_mask = obs["age_months"].between(old_lo, old_hi, inclusive="both")
    keep = young_mask | old_mask
    if sex_col in obs.columns and sex_value:
        keep &= obs[sex_col].astype(str).str.casefold().eq(str(sex_value).casefold())
    if int(keep.sum()) == 0:
        raise ValueError("No TMS cells matched the configured age/sex filters")

    obs = obs.loc[keep].copy()
    obs["age_group"] = np.where(obs["age_months"].between(young_lo, young_hi), "young", "old")
    x = adata.X[keep.to_numpy(), :]
    expr, donor_meta = _axis_from_matrix(x, obs, adata.var_names, donor_col)

    young_donors = donor_meta.index[donor_meta["age_group"] == "young"].tolist()
    old_donors = donor_meta.index[donor_meta["age_group"] == "old"].tolist()
    if len(young_donors) < 2 or len(old_donors) < 2:
        raise ValueError(
            "Need at least two young and two old TMS donors for a donor-level aging axis; "
            f"found young={len(young_donors)}, old={len(old_donors)}"
        )

    young_expr = expr.loc[young_donors].to_numpy(dtype=float)
    old_expr = expr.loc[old_donors].to_numpy(dtype=float)
    effect = old_expr.mean(axis=0) - young_expr.mean(axis=0)
    se = np.sqrt(old_expr.var(axis=0, ddof=1) / len(old_donors)
                 + young_expr.var(axis=0, ddof=1) / len(young_donors))
    test = ttest_ind(old_expr, young_expr, axis=0, equal_var=False, nan_policy="omit")
    p_value = np.nan_to_num(test.pvalue, nan=1.0)
    fdr = bh_fdr(p_value)

    axis = pd.DataFrame({
        "ensembl_gene_id": pd.Index(adata.var_names).astype(str).str.replace(r"\.\d+$", "", regex=True),
        "old_vs_young_logfc": effect,
        "se": se,
        "p_value": p_value,
        "fdr": fdr,
        "mean_old_logcpm": old_expr.mean(axis=0),
        "mean_young_logcpm": young_expr.mean(axis=0),
        "n_old_donors": len(old_donors),
        "n_young_donors": len(young_donors),
        "n_old_cells": int(donor_meta.loc[old_donors, "n_cells"].sum()),
        "n_young_cells": int(donor_meta.loc[young_donors, "n_cells"].sum()),
        "old_months_min": old_lo,
        "old_months_max": old_hi,
        "young_months_min": young_lo,
        "young_months_max": young_hi,
        "source_h5ad": str(h5ad_path),
    })
    axis = axis.groupby("ensembl_gene_id", as_index=False).mean(numeric_only=True)

    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    axis.to_csv(outpath, sep="\t", index=False)
    donor_expr_path = outpath.with_suffix(".donor_logcpm.tsv.gz")
    donor_meta_path = outpath.with_suffix(".donor_meta.tsv")
    expr.to_csv(donor_expr_path, sep="\t", compression="gzip")
    donor_meta.to_csv(donor_meta_path, sep="\t")
    manifest = {
        "source_h5ad": str(h5ad_path),
        "young_months": [young_lo, young_hi],
        "old_months": [old_lo, old_hi],
        "n_young_donors": len(young_donors),
        "n_old_donors": len(old_donors),
        "young_donors": young_donors,
        "old_donors": old_donors,
        "n_genes": int(len(axis)),
        "donor_logcpm": str(donor_expr_path),
        "donor_meta": str(donor_meta_path),
    }
    with open(outpath.with_suffix(".manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
    return axis


def load_or_build_aging_axis(
    path: str | Path,
    gene_col: str,
    effect_col: str,
    source_h5ad: str | Path | None = None,
    young_months: Sequence[float] = (3,),
    old_months: Sequence[float] = (18,),
) -> pd.Series:
    """Load an aging-axis TSV, building it from local TMS H5AD if absent."""
    path = Path(path)
    if not path.exists():
        h5ad_path = Path(source_h5ad) if source_h5ad else discover_tms_h5ad()
        build_tms_kidney_female_aging_axis(
            h5ad_path,
            path,
            young_months=young_months,
            old_months=old_months,
        )
    return load_aging_axis(path, gene_col, effect_col)


def load_tms_donor_axis_source(axis_path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load donor-level pseudobulk files paired with a built TMS aging axis."""
    axis_path = Path(axis_path)
    donor_expr_path = axis_path.with_suffix(".donor_logcpm.tsv.gz")
    donor_meta_path = axis_path.with_suffix(".donor_meta.tsv")
    if not donor_expr_path.exists() or not donor_meta_path.exists():
        raise FileNotFoundError(
            "TMS donor bootstrap files are missing. Rebuild the aging axis from "
            f"H5AD so these files exist: {donor_expr_path}, {donor_meta_path}"
        )
    expr = pd.read_csv(donor_expr_path, sep="\t", compression="gzip", index_col=0)
    expr.columns = expr.columns.astype(str).str.replace(r"\.\d+$", "", regex=True)
    meta = pd.read_csv(donor_meta_path, sep="\t", index_col=0)
    if "age_group" not in meta.columns:
        raise ValueError(f"{donor_meta_path} missing required age_group column")
    common = [d for d in expr.index if d in meta.index]
    if not common:
        raise ValueError("No overlapping donor IDs between TMS donor expression and metadata")
    return expr.loc[common], meta.loc[common]


def bootstrap_tms_aging_axis(
    donor_expr: pd.DataFrame,
    donor_meta: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.Series:
    """Bootstrap the TMS donor-level old-minus-young aging axis."""
    young = donor_meta.index[donor_meta["age_group"].astype(str).eq("young")].to_numpy()
    old = donor_meta.index[donor_meta["age_group"].astype(str).eq("old")].to_numpy()
    if len(young) < 2 or len(old) < 2:
        raise ValueError(f"Need at least two young and two old TMS donors; found young={len(young)}, old={len(old)}")
    young_draw = rng.choice(young, size=len(young), replace=True)
    old_draw = rng.choice(old, size=len(old), replace=True)
    effect = donor_expr.loc[old_draw].mean(axis=0) - donor_expr.loc[young_draw].mean(axis=0)
    effect.index = effect.index.astype(str)
    return effect.astype(float)


def intersect_axis_features(
    flight_vector: pd.Series,
    aging_axis: pd.Series,
    min_features: int = 3,
) -> tuple[pd.Series, pd.Series]:
    """Align flight and aging vectors on a shared feature universe."""
    f = flight_vector.copy()
    a = aging_axis.copy()
    f.index = f.index.astype(str)
    a.index = a.index.astype(str)
    common = f.index.intersection(a.index)
    if len(common) < min_features:
        raise ValueError(
            f"External-axis projection needs at least {min_features} shared features; "
            f"found {len(common)}"
        )
    return f.loc[common].astype(float), a.loc[common].astype(float)


def project_vector_onto_axis(
    flight_vector: pd.Series,
    aging_axis: pd.Series,
    weights: pd.Series | None = None,
    min_features: int = 3,
) -> tuple[Decomposition, pd.DataFrame]:
    """Project a flight vector onto an external aging axis."""
    f, a = intersect_axis_features(flight_vector, aging_axis, min_features=min_features)
    w_arr = None
    if weights is not None:
        w = weights.copy()
        w.index = w.index.astype(str)
        w_arr = w.reindex(f.index).to_numpy(dtype=float)
        if np.isnan(w_arr).any():
            raise ValueError("weights must cover every shared external-axis feature")
    dec = compute_decomposition(f.to_numpy(), a.to_numpy(), weights=w_arr)
    per_feature = pd.DataFrame({
        "feature": f.index,
        "flight_effect": f.to_numpy(dtype=float),
        "external_aging_effect": a.to_numpy(dtype=float),
    })
    if w_arr is not None:
        per_feature["weight"] = w_arr
    return dec, per_feature


def build_rrrm2_flight_vector(
    feature_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    arm: str,
    env_flight: str = "FLT",
    env_control: str = "GC",
) -> pd.Series:
    """Return age-pooled FLT-GC feature effect for one RRRM-2 arm."""
    meta = normalize_labels(metadata)
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    common = [s for s in feature_matrix.index if s in meta.index]
    if not common:
        raise ValueError("No overlapping samples between feature matrix and metadata")
    x = feature_matrix.loc[common]
    m = meta.loc[common]
    fmask = (m["Arm"].astype(str) == arm) & (m["EnvGroup"].astype(str) == env_flight)
    cmask = (m["Arm"].astype(str) == arm) & (m["EnvGroup"].astype(str) == env_control)
    if fmask.sum() < 1 or cmask.sum() < 1:
        raise ValueError(
            f"Need both FLT and GC samples for arm={arm}; found FLT={int(fmask.sum())}, GC={int(cmask.sum())}"
        )
    return x.loc[fmask].mean(axis=0) - x.loc[cmask].mean(axis=0)


def bootstrap_external_projection(
    feature_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    aging_axis: pd.Series,
    arm: str,
    n_iterations: int,
    rng: np.random.Generator,
    alpha: float = 0.05,
    tms_donor_expr: pd.DataFrame | None = None,
    tms_donor_meta: pd.DataFrame | None = None,
) -> tuple[Decomposition, pd.DataFrame, pd.DataFrame]:
    """Bootstrap the RRRM-2 side of the external-axis projection."""
    meta = normalize_labels(metadata)
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    common = [s for s in feature_matrix.index if s in meta.index]
    x = feature_matrix.loc[common]
    m = meta.loc[common].reset_index(drop=True)
    x = x.reset_index(drop=True)

    arm_mask = m["Arm"].astype(str).eq(arm)
    x_arm = x.loc[arm_mask].reset_index(drop=True)
    m_arm = m.loc[arm_mask].reset_index(drop=True)
    strata = m_arm[["Age", "EnvGroup"]].astype(str).agg("|".join, axis=1)

    point_vector = build_rrrm2_flight_vector(
        x_arm.set_index(pd.Index(m_arm[sample_col])),
        m_arm,
        arm=arm,
    )
    point, _ = project_vector_onto_axis(point_vector, aging_axis)

    rows: list[dict[str, float | int]] = []
    for b in range(int(n_iterations)):
        idx = stratified_bootstrap_indices(strata, rng)
        xb = x_arm.iloc[idx].reset_index(drop=True)
        mb = m_arm.iloc[idx].reset_index(drop=True).copy()
        mb[sample_col] = [f"boot_{b}_{i}" for i in range(len(mb))]
        xb.index = mb[sample_col]
        fv = build_rrrm2_flight_vector(xb, mb, arm=arm)
        axis_b = (
            bootstrap_tms_aging_axis(tms_donor_expr, tms_donor_meta, rng)
            if tms_donor_expr is not None and tms_donor_meta is not None
            else aging_axis
        )
        dec, _ = project_vector_onto_axis(fv, axis_b)
        rows.append({
            "iteration": b,
            "beta": dec.beta,
            "cos": dec.cos,
            "rho": dec.rho,
            "a_flt_norm": dec.a_flt_norm,
            "a_gc_norm": dec.a_gc_norm,
            "residual_norm": dec.residual_norm,
        })
    boot = pd.DataFrame(rows)
    summary = summarize_bootstrap_decomposition(point, boot, alpha=alpha)
    return point, boot, summary


def load_gene_feature_matrix(rtech_path: str | Path) -> pd.DataFrame:
    """Load RRRM-2 residualized expression as samples x genes."""
    rtech = pd.read_csv(rtech_path, sep="\t", compression="gzip", index_col=0)
    rtech.index = rtech.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    return rtech.T


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="Project RRRM-2 flight effects onto an external aging axis")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--axis", default="data/external/aging_reference/tms_kidney_female_aging_axis.tsv")
    ap.add_argument("--source-h5ad", default="",
                    help="Local TMS kidney female H5AD to use when --axis is absent")
    ap.add_argument("--gene-col", default="ensembl_gene_id")
    ap.add_argument("--effect-col", default="old_vs_young_logfc")
    ap.add_argument("--outdir", default="data/results/contrast_vectors/external_aging_axis")
    ap.add_argument("--arms", default="ISS-T,LAR")
    ap.add_argument("--n-bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260101)
    args = ap.parse_args(argv)

    feature_matrix = load_gene_feature_matrix(REPO_ROOT / args.rtech)
    meta = pd.read_csv(REPO_ROOT / args.meta, sep="\t", compression="gzip")
    axis = load_or_build_aging_axis(
        REPO_ROOT / args.axis,
        args.gene_col,
        args.effect_col,
        source_h5ad=(REPO_ROOT / args.source_h5ad) if args.source_h5ad else None,
    )
    outdir = REPO_ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    manifest: dict[str, object] = {
        "axis": args.axis,
        "gene_col": args.gene_col,
        "effect_col": args.effect_col,
        "n_bootstrap": args.n_bootstrap,
        "arms": [],
    }
    for arm in [a.strip() for a in args.arms.split(",") if a.strip()]:
        arm_dir = outdir / arm
        arm_dir.mkdir(parents=True, exist_ok=True)
        point, boot, summary = bootstrap_external_projection(
            feature_matrix, meta, axis, arm, args.n_bootstrap, rng
        )
        boot.to_csv(arm_dir / "bootstrap.tsv", sep="\t", index=False)
        summary.to_csv(arm_dir / f"tms_projection_{arm}.tsv", sep="\t", index=False)
        manifest["arms"].append({"arm": arm, "point": point.as_dict()})

    with open(outdir / "manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
