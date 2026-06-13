#!/usr/bin/env python3
"""Run recurrent remodeling-axis and network-neighborhood prioritization."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from scipy import stats

SCRIPT_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_REPO_ROOT))

from src.common import REPO_ROOT, bh_fdr, find_sample_col, id_map_lookup, normalize_labels
from src.networks.mechanism_axis import (
    age_balanced_feature_effect,
    aggregate_lioness_nodes,
    composite_gene_priority,
    define_recurrent_axis,
    fisher_topk_enrichment,
    flight_minus_control,
    load_gene_set_config,
    osd253_discriminator_table,
    preranked_enrichment,
    proportional_top_ks,
    remodeling_score,
    resolve_gene_sets,
    rrrm2_arm_rows,
    safe_z,
    score_gene_sets,
    select_osd253_scenario,
    stratified_arm_flight_effect,
)


ARMS = ("ISS-T", "LAR")
DEFAULT_OSD253_SCENARIOS = ("original_gc_blue_light", "rerun_gc_white_light")


def log(message: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    print(f"{ts} [MECH] {message}")


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


def read_rrrm2_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    return expr


def read_vst(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, index_col=0)
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    return expr


def load_rrrm2_inputs(rtech_path: Path, meta_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    expr = read_rrrm2_expression(rtech_path)
    meta = pd.read_csv(meta_path, sep="\t", compression="gzip")
    meta = normalize_labels(meta)
    sample_col = find_sample_col(meta)
    meta[sample_col] = meta[sample_col].astype(str)
    meta = meta.set_index(sample_col, drop=False)
    common = [s for s in expr.columns.astype(str) if s in meta.index]
    if not common:
        raise ValueError("No overlapping samples between RRRM-2 expression and metadata")
    return expr[common], meta.loc[common].copy(), sample_col


def resolve_external_file(study_dir: Path, pattern: str) -> Path:
    matches = sorted(study_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No file matching {pattern!r} in {study_dir}")
    return matches[0]


def existing_dirs(pattern: str) -> list[Path]:
    return [p for p in sorted(REPO_ROOT.glob(pattern), reverse=True) if p.is_dir()]


def resolve_phase2_dir(requested: str | None, run_root: Path) -> Path:
    candidates: list[Path] = []
    if requested:
        candidates.append(Path(requested))
    candidates.append(run_root / "networks")
    candidates.extend(existing_dirs("data/results/run_*/networks"))
    candidates.extend(existing_dirs("data/processed/networks/run_*"))
    candidates.append(REPO_ROOT / "data/processed/networks/phase2")
    for candidate in candidates:
        candidate = candidate if candidate.is_absolute() else REPO_ROOT / candidate
        if (candidate / "edge_index.tsv").exists() and (candidate / "lioness_samples.txt").exists():
            return candidate
    raise FileNotFoundError("Could not resolve Phase 2/LIONESS network directory")


def resolve_rewiring_dir(requested: str | None, run_root: Path) -> Path | None:
    candidates: list[Path] = []
    if requested:
        candidates.append(Path(requested))
    candidates.append(run_root / "phase3_rewiring")
    candidates.extend(existing_dirs("data/results/run_*/phase3_rewiring"))
    for candidate in candidates:
        candidate = candidate if candidate.is_absolute() else REPO_ROOT / candidate
        if (candidate / "ISS_T_YNG_FLT_minus_GC_rewiring_agg.tsv").exists():
            return candidate
    return None


def osd513_rows_from_columns(columns: Sequence[str]) -> pd.DataFrame:
    flt = sorted(c for c in columns if "_FLT_" in c)
    gc = sorted(c for c in columns if "_GC_" in c)
    viv = sorted(c for c in columns if "_VIV_" in c)
    return pd.DataFrame({
        "Sample Name": flt + gc + viv,
        "condition": ["FLT"] * len(flt) + ["GC"] * len(gc) + ["VIV"] * len(viv),
        "study": ["OSD-513"] * (len(flt) + len(gc) + len(viv)),
        "scenario": ["powered_hgc"] * (len(flt) + len(gc) + len(viv)),
    })


def rrrm2_metadata_for_output(meta: pd.DataFrame, sample_col: str) -> pd.DataFrame:
    cols = [c for c in [sample_col, "Age", "Arm", "EnvGroup"] if c in meta.columns]
    out = meta[cols].copy()
    out = out.rename(columns={sample_col: "Sample Name"})
    out["study"] = "RRRM-2"
    out["scenario"] = "primary"
    out["condition"] = out["EnvGroup"].astype(str)
    return out


def osd253_metadata_for_output(meta: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "Sample Name",
        "Factor Value[Strain]",
        "Factor Value[Spaceflight]",
        "Factor Value[Duration]",
        "Factor Value[Treatment]",
    ]
    out = meta[[c for c in cols if c in meta.columns]].copy()
    out["study"] = "OSD-253"
    out["scenario"] = "all_samples"
    out["condition"] = out["Factor Value[Spaceflight]"].astype(str)
    out["strain"] = out.get("Factor Value[Strain]", pd.Series(index=out.index, dtype=str)).astype(str)
    out["duration"] = out.get("Factor Value[Duration]", pd.Series(index=out.index, dtype=str)).astype(str)
    return out


def build_score_output(
    metadata: pd.DataFrame,
    scores: pd.Series,
    *,
    score_name: str = "remodeling_score",
) -> pd.DataFrame:
    out = metadata.copy()
    out[score_name] = out["Sample Name"].astype(str).map(scores)
    return out


def load_phase2_lioness(phase2_dir: Path, meta: pd.DataFrame, sample_col: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    edges = pd.read_csv(phase2_dir / "edge_index.tsv", sep="\t")
    samples = [x.strip() for x in (phase2_dir / "lioness_samples.txt").read_text().splitlines() if x.strip()]
    genes_path = phase2_dir / "phase2_genes.txt"
    if genes_path.exists():
        genes = [x.strip() for x in genes_path.read_text().splitlines() if x.strip()]
    else:
        genes = sorted(pd.unique(edges[["gene_i", "gene_j"]].to_numpy().ravel()))
    weights_path = phase2_dir / "lioness_edges.npy"
    if not weights_path.exists():
        weights_path = phase2_dir / "lioness_z_edges.npy"
    weights = np.load(weights_path, mmap_mode="r")
    if weights.shape[0] != len(samples) or weights.shape[1] != len(edges):
        raise ValueError(f"LIONESS shape {weights.shape} does not match samples={len(samples)} edges={len(edges)}")
    sample_pos = {s: i for i, s in enumerate(samples)}
    common = [s for s in samples if s in set(meta[sample_col].astype(str))]
    if not common:
        raise ValueError("No overlapping LIONESS samples and RRRM-2 metadata")
    row_idx = [sample_pos[s] for s in common]
    node_features = aggregate_lioness_nodes(np.asarray(weights[row_idx, :]), edges, genes)
    node_features.index = common
    return node_features, edges


def load_embedding_displacement(rewiring_dir: Path | None, genes: Sequence[str]) -> tuple[pd.Series, pd.Series]:
    """Load age-averaged node2vec displacement from Phase 3 rewiring aggregates."""
    gene_index = pd.Index([str(g) for g in genes])
    if rewiring_dir is None:
        log("No phase3_rewiring directory found; embedding displacement terms will be zero.")
        zero = pd.Series(0.0, index=gene_index)
        return zero, zero

    def load_one(name: str) -> pd.Series:
        path = rewiring_dir / f"{name}_FLT_minus_GC_rewiring_agg.tsv"
        if not path.exists():
            return pd.Series(np.nan, index=gene_index)
        df = pd.read_csv(path, sep="\t")
        return df.set_index("gene")["rewiring_mean"].reindex(gene_index)

    iss = pd.concat([load_one("ISS_T_YNG"), load_one("ISS_T_OLD")], axis=1).mean(axis=1)
    lar = pd.concat([load_one("LAR_YNG"), load_one("LAR_OLD")], axis=1).mean(axis=1)
    return iss.fillna(iss.median()), lar.fillna(lar.median())


def load_osd513_gene_support(vst: pd.DataFrame, rows: pd.DataFrame) -> pd.Series:
    flt = rows.loc[rows["condition"].eq("FLT"), "Sample Name"].astype(str).tolist()
    gc = rows.loc[rows["condition"].eq("GC"), "Sample Name"].astype(str).tolist()
    flt = [s for s in flt if s in vst.columns]
    gc = [s for s in gc if s in vst.columns]
    if not flt or not gc:
        return pd.Series(dtype=float)
    effect = vst[flt].mean(axis=1) - vst[gc].mean(axis=1)
    return pd.Series(safe_z(effect.abs()), index=effect.index)


def load_iss_t_abs_logfc(run_root: Path) -> pd.Series:
    candidates = [
        run_root / "phase5_silent_shifters_strict/ISS_T_YNG_FLT_minus_GC_all_genes_DE_rewiring.tsv",
        run_root / "phase5_silent_shifters_strict/ISS_T_OLD_FLT_minus_GC_all_genes_DE_rewiring.tsv",
        REPO_ROOT / "data/processed/gene_level_DE/ISS_T_YNG_FLT_vs_GC_gene_DE.tsv",
        REPO_ROOT / "data/processed/gene_level_DE/ISS_T_OLD_FLT_vs_GC_gene_DE.tsv",
    ]
    series: list[pd.Series] = []
    for path in candidates:
        if not path.exists():
            continue
        df = pd.read_csv(path, sep="\t")
        col = "log2FC_shrunken" if "log2FC_shrunken" in df.columns else "log2FC"
        if "gene" in df.columns and col in df.columns:
            series.append(df.set_index("gene")[col].abs())
    if not series:
        return pd.Series(dtype=float)
    return pd.concat(series, axis=1).mean(axis=1)


def annotate_gene_priority(
    ranked: pd.DataFrame,
    id_map_path: Path,
    resolved_sets: dict,
) -> pd.DataFrame:
    ens_to_symbol, _ = id_map_lookup(id_map_path)
    membership: dict[str, list[str]] = {}
    roles: dict[str, list[str]] = {}
    for name, res in resolved_sets.items():
        for gid in res.scored_gene_ids:
            membership.setdefault(gid, []).append(name)
            roles.setdefault(gid, []).append(res.role)
    out = ranked.copy()
    out["mgi_symbol"] = out["gene"].astype(str).map(ens_to_symbol).fillna("")
    out["mechanism_sets"] = out["gene"].astype(str).map(lambda g: ",".join(sorted(membership.get(g, []))))
    out["mechanism_roles"] = out["gene"].astype(str).map(lambda g: ",".join(sorted(set(roles.get(g, [])))))
    return out


def mechanism_associations(
    sample_scores: pd.DataFrame,
    mechanism_cols: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (study, scenario), sub in sample_scores.groupby(["study", "scenario"], dropna=False):
        if "remodeling_score" not in sub.columns:
            continue
        for mech in mechanism_cols:
            if mech not in sub.columns:
                continue
            vals = sub[["remodeling_score", mech]].replace([np.inf, -np.inf], np.nan).dropna()
            if len(vals) >= 4 and vals["remodeling_score"].nunique() > 1 and vals[mech].nunique() > 1:
                rho, p = stats.spearmanr(vals["remodeling_score"], vals[mech])
            else:
                rho, p = np.nan, np.nan
            rows.append({
                "study": study,
                "scenario": scenario,
                "comparison": "score_vs_remodeling_spearman",
                "mechanism": mech,
                "estimate": float(rho) if np.isfinite(rho) else np.nan,
                "p_value": float(p) if np.isfinite(p) else np.nan,
                "n_samples": int(len(vals)),
            })
    out = pd.DataFrame(rows)
    if not out.empty and out["p_value"].notna().any():
        mask = out["p_value"].notna()
        out.loc[mask, "q_bh"] = bh_fdr(out.loc[mask, "p_value"].to_numpy())
    return out


def append_flight_effect_associations(
    assoc: pd.DataFrame,
    rrrm2_mech: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
    osd513_mech: pd.DataFrame,
    osd513_rows: pd.DataFrame,
    mechanism_cols: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for mech in mechanism_cols:
        for arm in ARMS:
            effect = stratified_arm_flight_effect(rrrm2_mech, meta, sample_col, arm, mech)
            rows.append({
                "study": "RRRM-2",
                "scenario": "primary",
                "comparison": f"{arm}_age_balanced_flight_minus_gc",
                "mechanism": mech,
                "estimate": effect,
                "p_value": np.nan,
                "q_bh": np.nan,
                "n_samples": int(len(rrrm2_arm_rows(meta, arm))),
            })
        flt = osd513_rows.loc[osd513_rows["condition"].eq("FLT"), "Sample Name"].astype(str).tolist()
        gc = osd513_rows.loc[osd513_rows["condition"].eq("GC"), "Sample Name"].astype(str).tolist()
        flt = [s for s in flt if s in osd513_mech.index]
        gc = [s for s in gc if s in osd513_mech.index]
        if flt and gc and mech in osd513_mech.columns:
            effect = float(osd513_mech.loc[flt, mech].mean() - osd513_mech.loc[gc, mech].mean())
            _, p = stats.ttest_ind(osd513_mech.loc[flt, mech], osd513_mech.loc[gc, mech], equal_var=False)
        else:
            effect, p = np.nan, np.nan
        rows.append({
            "study": "OSD-513",
            "scenario": "powered_hgc",
            "comparison": "flight_minus_gc",
            "mechanism": mech,
            "estimate": effect,
            "p_value": float(p) if np.isfinite(p) else np.nan,
            "q_bh": np.nan,
            "n_samples": int(len(flt) + len(gc)),
        })
    extra = pd.DataFrame(rows)
    if not extra.empty and extra["p_value"].notna().any():
        mask = extra["p_value"].notna()
        extra.loc[mask, "q_bh"] = bh_fdr(extra.loc[mask, "p_value"].to_numpy())
    return pd.concat([assoc, extra], ignore_index=True)


def convergence_summary(
    coverage: pd.DataFrame,
    associations: pd.DataFrame,
    discriminator: pd.DataFrame,
    topk: pd.DataFrame,
    preranked: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    mechanisms = sorted(set(coverage["gene_set"]) | set(associations.get("mechanism", [])))
    for mech in mechanisms:
        cov = coverage[coverage["gene_set"].eq(mech)]
        assoc = associations[associations["mechanism"].eq(mech)] if "mechanism" in associations else pd.DataFrame()
        disc = discriminator[discriminator["score_name"].eq(mech)] if "score_name" in discriminator else pd.DataFrame()
        tk = topk[(topk["gene_set"].eq(mech)) & (topk["top_k"].eq(50))] if not topk.empty else pd.DataFrame()
        pr = preranked[preranked["gene_set"].eq(mech)] if not preranked.empty else pd.DataFrame()
        best_assoc_q = float(assoc["q_bh"].min()) if not assoc.empty and assoc["q_bh"].notna().any() else np.nan
        best_top50_q = float(tk["q_bh"].min()) if not tk.empty and tk["q_bh"].notna().any() else np.nan
        best_ranked_q = float(pr["q_bh"].min()) if not pr.empty and pr["q_bh"].notna().any() else np.nan
        interpretation = ""
        if not disc.empty:
            interpretation = ",".join(sorted(set(disc["interpretation"].dropna().astype(str))))
        rows.append({
            "mechanism": mech,
            "role": cov["role"].iloc[0] if not cov.empty else "",
            "n_scored": int(cov["n_scored"].iloc[0]) if not cov.empty else np.nan,
            "best_score_association_q": best_assoc_q,
            "top50_gene_priority_q": best_top50_q,
            "preranked_gene_priority_q": best_ranked_q,
            "osd253_interpretation": interpretation,
            "summary_label": (
                "network_priority_and_sample_score"
                if np.isfinite(best_top50_q) and best_top50_q <= 0.10 and np.isfinite(best_assoc_q) and best_assoc_q <= 0.10
                else "context"
            ),
        })
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run-id", default="")
    ap.add_argument("--outdir", default="")
    ap.add_argument("--contrast-dir", default="")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--phase2-dir", default="")
    ap.add_argument("--rewiring-dir", default="")
    ap.add_argument("--external-root", default="data/external/osdr")
    ap.add_argument("--id-map", default="data/processed/resources/id_map.tsv")
    ap.add_argument("--gene-sets", default="config/gene_sets.yaml")
    ap.add_argument("--mechanism-gene-sets", default="config/mechanism_gene_sets.yaml")
    ap.add_argument("--axis-source", default="rrrm2_iss_osd513", choices=["rrrm2_iss_osd513"])
    ap.add_argument("--include-osd253-strain", action="store_true", default=True)
    ap.add_argument("--n-bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260519)
    ap.add_argument("--primary-transport-set", default="dct_ncc_wnk")
    ap.add_argument("--sensitivity-transport-set", default="ion_transport")
    ap.add_argument("--min-genes", type=int, default=3)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    run_id = args.run_id or datetime.now().strftime("run_%Y%m%d_%H%M%S_mechanism_axis")

    contrast_dir = Path(args.contrast_dir) if args.contrast_dir else Path(args.outdir).parent if args.outdir else None
    if contrast_dir is not None and not contrast_dir.is_absolute():
        contrast_dir = REPO_ROOT / contrast_dir
    outdir = Path(args.outdir) if args.outdir else (contrast_dir / "mechanism_axis" if contrast_dir else REPO_ROOT / "data/results" / run_id / "contrast_vectors/mechanism_axis")
    if not outdir.is_absolute():
        outdir = REPO_ROOT / outdir
    outdir.mkdir(parents=True, exist_ok=True)

    rtech_path = REPO_ROOT / args.rtech if not Path(args.rtech).is_absolute() else Path(args.rtech)
    meta_path = REPO_ROOT / args.meta if not Path(args.meta).is_absolute() else Path(args.meta)
    id_map_path = REPO_ROOT / args.id_map if not Path(args.id_map).is_absolute() else Path(args.id_map)
    gene_sets_path = REPO_ROOT / args.gene_sets if not Path(args.gene_sets).is_absolute() else Path(args.gene_sets)
    mechanism_sets_path = REPO_ROOT / args.mechanism_gene_sets if not Path(args.mechanism_gene_sets).is_absolute() else Path(args.mechanism_gene_sets)
    external_root = REPO_ROOT / args.external_root if not Path(args.external_root).is_absolute() else Path(args.external_root)
    run_root = outdir.parents[1] if outdir.name == "mechanism_axis" else outdir.parent

    log(f"Run ID: {run_id}")
    log(f"Output: {outdir}")
    log("Loading RRRM-2, OSD-513, and OSD-253 expression matrices")

    rrrm2_expr, rrrm2_meta, rrrm2_sample_col = load_rrrm2_inputs(rtech_path, meta_path)
    osd513_dir = external_root / "OSD-513"
    osd253_dir = external_root / "OSD-253"
    osd513_vst = read_vst(resolve_external_file(osd513_dir, "*VST_Counts*.csv"))
    osd253_vst = read_vst(resolve_external_file(osd253_dir, "*VST_Counts*.csv"))
    osd253_meta = pd.read_csv(resolve_external_file(osd253_dir, "*runsheet.csv"))

    axis_gene_sets = load_gene_set_config(gene_sets_path)
    axis_resolved, axis_coverage = resolve_gene_sets(
        axis_gene_sets,
        id_map_path,
        set(rrrm2_expr.index) | set(osd513_vst.index) | set(osd253_vst.index),
        min_genes=args.min_genes,
    )
    rrrm2_pathways = score_gene_sets(rrrm2_expr, axis_resolved)
    osd513_pathways = score_gene_sets(osd513_vst, axis_resolved)
    osd253_pathways = score_gene_sets(osd253_vst, axis_resolved)

    osd513_rows = osd513_rows_from_columns(osd513_pathways.index)
    rrrm2_iss_effect = flight_minus_control(
        rrrm2_pathways,
        rrrm2_arm_rows(rrrm2_meta, "ISS-T"),
        sample_col=rrrm2_sample_col,
    )
    osd513_effect = flight_minus_control(osd513_pathways, osd513_rows, sample_col="Sample Name")
    axis_table, axis_weights = define_recurrent_axis(
        rrrm2_iss_effect,
        osd513_effect,
        primary_transport_set=args.primary_transport_set,
        sensitivity_transport_set=args.sensitivity_transport_set,
    )
    axis_table.to_csv(outdir / "recurrent_axis_weights.tsv", sep="\t", index=False)
    axis_coverage.to_csv(outdir / "axis_gene_set_coverage.tsv", sep="\t", index=False)

    rrrm2_remodeling = remodeling_score(rrrm2_pathways, axis_weights)
    osd513_remodeling = remodeling_score(osd513_pathways, axis_weights)
    osd253_remodeling = remodeling_score(osd253_pathways, axis_weights)

    axis_gene_ids = set()
    for name, res in axis_resolved.items():
        if name in set(axis_weights.index) | {args.primary_transport_set, args.sensitivity_transport_set}:
            axis_gene_ids |= set(res.matched_gene_ids)

    mechanism_config = load_gene_set_config(mechanism_sets_path)
    mechanism_resolved, mechanism_coverage = resolve_gene_sets(
        mechanism_config,
        id_map_path,
        set(rrrm2_expr.index) | set(osd513_vst.index) | set(osd253_vst.index),
        min_genes=args.min_genes,
        deoverlap_gene_ids=axis_gene_ids,
    )
    mechanism_coverage.to_csv(outdir / "mechanism_gene_set_coverage.tsv", sep="\t", index=False)

    rrrm2_mech = score_gene_sets(rrrm2_expr, mechanism_resolved)
    osd513_mech = score_gene_sets(osd513_vst, mechanism_resolved)
    osd253_mech = score_gene_sets(osd253_vst, mechanism_resolved)

    sample_score_frames = [
        build_score_output(rrrm2_metadata_for_output(rrrm2_meta, rrrm2_sample_col), rrrm2_remodeling),
        build_score_output(osd513_rows, osd513_remodeling),
        build_score_output(osd253_metadata_for_output(osd253_meta), osd253_remodeling),
    ]
    sample_scores = pd.concat(sample_score_frames, ignore_index=True, sort=False)
    sample_scores.to_csv(outdir / "sample_remodeling_scores.tsv", sep="\t", index=False)

    mechanism_frames = []
    for study, scenario, meta_out, mech, rem in [
        ("RRRM-2", "primary", rrrm2_metadata_for_output(rrrm2_meta, rrrm2_sample_col), rrrm2_mech, rrrm2_remodeling),
        ("OSD-513", "powered_hgc", osd513_rows, osd513_mech, osd513_remodeling),
        ("OSD-253", "all_samples", osd253_metadata_for_output(osd253_meta), osd253_mech, osd253_remodeling),
    ]:
        wide = meta_out.copy()
        wide["remodeling_score"] = wide["Sample Name"].astype(str).map(rem)
        wide = wide.merge(mech, left_on="Sample Name", right_index=True, how="left")
        wide["study"] = study
        wide["scenario"] = scenario
        mechanism_frames.append(wide)
    mechanism_matrix = pd.concat(mechanism_frames, ignore_index=True, sort=False)
    mechanism_matrix.to_csv(outdir / "mechanism_score_matrix.tsv", sep="\t", index=False)

    mechanism_cols = list(mechanism_resolved.keys())
    assoc = mechanism_associations(mechanism_matrix, mechanism_cols)
    assoc = append_flight_effect_associations(
        assoc,
        rrrm2_mech,
        rrrm2_meta,
        rrrm2_sample_col,
        osd513_mech,
        osd513_rows,
        mechanism_cols,
    )
    assoc.to_csv(outdir / "mechanism_score_associations.tsv", sep="\t", index=False)

    osd253_score_table = osd253_mech.copy()
    osd253_score_table.insert(0, "remodeling_score", osd253_remodeling)
    discriminator = pd.DataFrame()
    discriminator_boot = pd.DataFrame()
    if args.include_osd253_strain:
        score_cols = [
            c for c in [
                "remodeling_score",
                "tlr4_innate",
                "fibrosis_tgfb_emt",
                "ecm_organization",
                "integrin_cell_adhesion",
                "mmp_adam_proteolysis",
                "dct_ncc_wnk_transport",
            ] if c in osd253_score_table.columns
        ]
        discriminator, discriminator_boot = osd253_discriminator_table(
            osd253_score_table,
            osd253_meta,
            scenarios=DEFAULT_OSD253_SCENARIOS,
            score_cols=score_cols,
            n_bootstrap=args.n_bootstrap,
            rng=rng,
        )
        discriminator.to_csv(outdir / "osd253_tlr4_discriminator.tsv", sep="\t", index=False)
        discriminator_boot.to_csv(outdir / "osd253_tlr4_discriminator_bootstrap.tsv", sep="\t", index=False)

    log("Computing LIONESS/node2vec gene prioritization")
    phase2_dir = resolve_phase2_dir(args.phase2_dir, run_root)
    rewiring_dir = resolve_rewiring_dir(args.rewiring_dir, run_root)
    node_features, _edges = load_phase2_lioness(phase2_dir, rrrm2_meta, rrrm2_sample_col)
    node_iss = age_balanced_feature_effect(node_features, rrrm2_meta, rrrm2_sample_col, "ISS-T")
    node_lar = age_balanced_feature_effect(node_features, rrrm2_meta, rrrm2_sample_col, "LAR")
    embedding_iss, embedding_lar = load_embedding_displacement(rewiring_dir, node_features.columns)
    osd513_support = load_osd513_gene_support(osd513_vst, osd513_rows)
    iss_abs_logfc = load_iss_t_abs_logfc(run_root)
    ranked = composite_gene_priority(
        embedding_iss,
        embedding_lar,
        node_iss,
        node_lar,
        osd513_support=osd513_support,
        iss_abs_logfc=iss_abs_logfc,
    )
    ranked = annotate_gene_priority(ranked, id_map_path, mechanism_resolved)
    ranked.to_csv(outdir / "gene_axis_priority.tsv", sep="\t", index=False)
    for k in (25, 50, 100):
        ranked.head(k).to_csv(outdir / f"gene_axis_priority_top{k}.tsv", sep="\t", index=False)

    top_ks = proportional_top_ks(len(ranked))
    topk_enrich = fisher_topk_enrichment(
        ranked,
        mechanism_resolved,
        score_col="composite_score",
        top_ks=top_ks,
    )
    prerank = preranked_enrichment(ranked, mechanism_resolved, score_col="composite_score")
    if ranked["silent_composite_score"].notna().any():
        topk_silent = fisher_topk_enrichment(
            ranked,
            mechanism_resolved,
            score_col="silent_composite_score",
            top_ks=top_ks,
        )
        prerank_silent = preranked_enrichment(ranked, mechanism_resolved, score_col="silent_composite_score")
        topk_enrich = pd.concat([topk_enrich, topk_silent], ignore_index=True)
        prerank = pd.concat([prerank, prerank_silent], ignore_index=True)
    topk_enrich.to_csv(outdir / "gene_priority_enrichment_topk.tsv", sep="\t", index=False)
    prerank.to_csv(outdir / "gene_priority_enrichment_preranked.tsv", sep="\t", index=False)

    summary = convergence_summary(mechanism_coverage, assoc, discriminator, topk_enrich, prerank)
    summary.to_csv(outdir / "mechanism_convergence_summary.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "recurrent ISS-T remodeling-axis mechanism scoring and exploratory LIONESS/node2vec gene prioritization",
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "git_commit": git_commit(),
        "parameters": {
            "axis_source": args.axis_source,
            "n_bootstrap": int(args.n_bootstrap),
            "seed": int(args.seed),
            "primary_transport_set": args.primary_transport_set,
            "sensitivity_transport_set": args.sensitivity_transport_set,
            "top_k_thresholds": top_ks,
            "top_k_rule": "proportional 1%/2%/5%/10% of ranked network universe",
            "gene_level_claim_boundary": "exploratory prioritization only; no individual-gene rewiring confirmation",
        },
        "inputs": {
            "rtech": str(rtech_path),
            "meta": str(meta_path),
            "gene_sets": str(gene_sets_path),
            "mechanism_gene_sets": str(mechanism_sets_path),
            "id_map": str(id_map_path),
            "external_root": str(external_root),
            "phase2_dir": str(phase2_dir),
            "rewiring_dir": str(rewiring_dir) if rewiring_dir else "",
            "rtech_sha256": file_sha256(rtech_path),
            "meta_sha256": file_sha256(meta_path),
            "gene_sets_sha256": file_sha256(gene_sets_path),
            "mechanism_gene_sets_sha256": file_sha256(mechanism_sets_path),
        },
        "outputs": [
            "recurrent_axis_weights.tsv",
            "sample_remodeling_scores.tsv",
            "mechanism_score_matrix.tsv",
            "mechanism_score_associations.tsv",
            "osd253_tlr4_discriminator.tsv",
            "gene_axis_priority.tsv",
            "gene_priority_enrichment_topk.tsv",
            "gene_priority_enrichment_preranked.tsv",
            "mechanism_convergence_summary.tsv",
        ],
    }
    (outdir / "mechanism_axis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log("Mechanism-axis prioritization completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
