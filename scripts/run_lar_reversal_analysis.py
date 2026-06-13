#!/usr/bin/env python3
"""Run LAR reversal and mechanism-switch analysis."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPT_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_REPO_ROOT))

from src.common import REPO_ROOT, bh_fdr, find_sample_col, normalize_labels
from src.networks.lar_reversal import (
    MECHANISM_FEATURES,
    STATE_COMPONENT_FEATURES,
    STATE_FEATURES,
    TARGETED_PANELS,
    arm_effect,
    component_effect_rows,
    cosine_similarity,
    gene_effect_scatter,
    interaction_table,
    leave_one_feature_out,
    mechanism_switch_decomposition,
    projection_beta,
    reversal_summary_for_features,
    rrrm2_flt_gc_rows,
    sample_feature_frame,
    score_symbol_panels,
    spearman_pair,
    vector_stats,
)


DEFAULT_MECH = REPO_ROOT / "data/results/run_20260519_000547_2500g/contrast_vectors/mechanism_axis"
DEFAULT_TSTATE = DEFAULT_MECH / "tubulointerstitial_state"
DEFAULT_WGCNA = REPO_ROOT / "data/results/run_20260505_remediated_2500g/wgcna/module_eigengenes.csv"
DEFAULT_CROSS = REPO_ROOT / "data/results/run_20260519_cosine_perm_null/contrast_vectors/cross_osdr_recurrence/cross_osdr_alignment_summary.tsv"


def log(message: str) -> None:
    print(f"{datetime.now().strftime('%H:%M:%S')} [LAR] {message}")


def resolve(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else REPO_ROOT / p


def file_sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def read_expression(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    expr = pd.read_csv(path, sep="\t", compression=compression, index_col=0)
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    return expr


def read_metadata(path: Path) -> tuple[pd.DataFrame, str]:
    meta = pd.read_csv(path, sep="\t", compression="gzip" if path.suffix == ".gz" else None)
    meta = normalize_labels(meta)
    sample_col = find_sample_col(meta)
    meta[sample_col] = meta[sample_col].astype(str)
    if "condition" not in meta.columns and "EnvGroup" in meta.columns:
        meta["condition"] = meta["EnvGroup"].astype(str)
    return meta, sample_col


def read_wgcna(path: Path, state_scores: pd.DataFrame) -> pd.DataFrame:
    eig = pd.read_csv(path)
    sample_col = "Sample Name" if "Sample Name" in eig.columns else eig.columns[0]
    eig = eig.rename(columns={sample_col: "Sample Name"})
    meta_cols = ["Sample Name", "Age", "Arm", "EnvGroup", "condition", "study", "scenario"]
    meta = state_scores[state_scores["study"].eq("RRRM-2")][[c for c in meta_cols if c in state_scores.columns]].drop_duplicates("Sample Name")
    out = meta.merge(eig, on="Sample Name", how="inner")
    out["study"] = "RRRM-2"
    out["scenario"] = "wgcna"
    return out


def node_summary_rows(priority: pd.DataFrame, *, feature_set: str, iss_col: str, lar_col: str) -> pd.DataFrame:
    if priority.empty or iss_col not in priority.columns or lar_col not in priority.columns:
        return pd.DataFrame()
    df = priority.dropna(subset=[iss_col, lar_col]).copy()
    if df.empty:
        return pd.DataFrame()
    iss = df.set_index("gene")[iss_col].astype(float)
    lar = df.set_index("gene")[lar_col].astype(float)
    stats = vector_stats(feature_set, "pooled", iss, lar).__dict__
    stats.update({
        "bootstrap_median": np.nan,
        "bootstrap_ci_low": np.nan,
        "bootstrap_ci_high": np.nan,
        "bootstrap_beta_ci_low": np.nan,
        "bootstrap_beta_ci_high": np.nan,
        "permutation_p_less": np.nan,
        "permutation_p_two_sided_abs": np.nan,
        "n_bootstrap": 0,
        "n_permutation": 0,
        "interpretation": "exploratory_precomputed_network_summary",
    })
    return pd.DataFrame([stats])


def build_target_gene_frame(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    id_map_path: Path,
    sample_col: str,
    max_gene_features: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    special_scores, coverage = score_symbol_panels(expression, str(id_map_path), TARGETED_PANELS)
    special_frame = sample_feature_frame(metadata, special_scores, sample_col=sample_col)
    common = [s for s in expression.columns if s in set(metadata[sample_col].astype(str))]
    expr_t = expression.loc[:, common].T
    if max_gene_features > 0 and expr_t.shape[1] > max_gene_features:
        variances = expr_t.var(axis=0, skipna=True).sort_values(ascending=False)
        expr_t = expr_t.loc[:, variances.head(max_gene_features).index]
    gene_frame = sample_feature_frame(metadata, expr_t, sample_col=sample_col)
    return special_frame, gene_frame.merge(special_frame[["Sample Name"] + list(special_scores.columns)], on="Sample Name", how="left"), coverage


def add_q_values(summary: pd.DataFrame) -> pd.DataFrame:
    out = summary.copy()
    for col in ["permutation_p_less", "permutation_p_two_sided_abs"]:
        if col in out.columns and out[col].notna().any():
            mask = out[col].notna()
            out[f"q_bh_{col}"] = np.nan
            out.loc[mask, f"q_bh_{col}"] = bh_fdr(out.loc[mask, col].to_numpy(dtype=float))
    return out


def child_rng(seed: int, label: str) -> np.random.Generator:
    """Return a stable per-analysis RNG so reordering blocks does not change draws."""
    digest = hashlib.blake2b(f"{seed}:{label}".encode("utf-8"), digest_size=8).digest()
    return np.random.default_rng(int.from_bytes(digest, byteorder="little", signed=False))


def comparison_external_from_summary(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t")
    keep = df[df["contrast"].astype(str).str.contains("OSD-513:powered_hgc", na=False)].copy()
    if keep.empty:
        return pd.DataFrame()
    keep["comparison"] = keep["arm"].astype(str) + "_vs_OSD513_cross_osdr_pathway"
    return keep


def model_call_table(vector_summary: pd.DataFrame, interactions: pd.DataFrame, switch: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    primary = vector_summary[
        vector_summary["feature_set"].isin(["state_space", "mechanism_scores", "special_clock_s1p_rbm3"])
        & vector_summary["age_scope"].eq("pooled")
    ].copy()
    for _, row in primary.iterrows():
        rows.append({
            "evidence_layer": row["feature_set"],
            "model_A_attenuation": bool(row["interpretation"] == "model_A_attenuation_candidate"),
            "model_B_reversal": bool(row["interpretation"] == "model_B_reversal_candidate"),
            "model_C_switch_or_partial": bool(row["interpretation"] in {"model_C_partial_reversal_or_switch", "mixed_or_inconclusive"}),
            "cos_lar_iss": row["cos_lar_iss"],
            "beta_lar_to_iss": row["beta_lar_to_iss"],
            "interpretation": row["interpretation"],
        })
    key = interactions[interactions["feature"].isin(["matrix_component", "dct_transport_component", "matrix_minus_dct", "s1p_s1pr3", "clock_core", "rbm3_preservation_expression"])]
    for _, row in key.iterrows():
        rows.append({
            "evidence_layer": f"component:{row['feature']}",
            "model_A_attenuation": bool(abs(row["lar_effect"]) < 0.5 * max(abs(row["iss_t_effect"]), 1e-12)),
            "model_B_reversal": bool(row["lar_opposes_iss"]),
            "model_C_switch_or_partial": bool(row["lar_opposes_iss"] or row["q_bh_within_interactions"] <= 0.10),
            "cos_lar_iss": np.nan,
            "beta_lar_to_iss": np.nan,
            "interpretation": "component_arm_divergence",
        })
    if not switch.empty:
        piv = switch.pivot_table(index="axis", columns="arm", values="multi_axis_regression_coefficient", aggfunc="first")
        for axis, row in piv.iterrows():
            iss = float(row.get("ISS-T", np.nan))
            lar = float(row.get("LAR", np.nan))
            rows.append({
                "evidence_layer": f"mechanism_switch_axis:{axis}",
                "model_A_attenuation": bool(np.isfinite(iss) and np.isfinite(lar) and abs(lar) < 0.5 * max(abs(iss), 1e-12)),
                "model_B_reversal": bool(np.isfinite(iss) and np.isfinite(lar) and np.sign(iss) != np.sign(lar) and abs(lar) >= 0.25 * max(abs(iss), 1e-12)),
                "model_C_switch_or_partial": bool(np.isfinite(lar) and abs(lar) > 0.15),
                "cos_lar_iss": np.nan,
                "beta_lar_to_iss": np.nan,
                "interpretation": "axis_decomposition",
            })
    return pd.DataFrame(rows)


def targeted_correlation_table(combined: pd.DataFrame, pairs: list[tuple[str, str]], family: str) -> pd.DataFrame:
    rows = []
    for x, y in pairs:
        if x not in combined.columns or y not in combined.columns:
            continue
        n, rho, p = spearman_pair(combined, x, y)
        rows.append({
            "analysis_family": family,
            "row_type": "correlation",
            "x": x,
            "y": y,
            "n": n,
            "spearman_rho": rho,
            "p_value": p,
        })
    out = pd.DataFrame(rows)
    if not out.empty and out["p_value"].notna().any():
        mask = out["p_value"].notna()
        out.loc[mask, "q_bh"] = bh_fdr(out.loc[mask, "p_value"].to_numpy(dtype=float))
    return out


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state-space-scores", default=str(DEFAULT_TSTATE / "state_space_scores.tsv"))
    ap.add_argument("--mechanism-score-matrix", default=str(DEFAULT_MECH / "mechanism_score_matrix.tsv"))
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--id-map", default="data/processed/resources/id_map.tsv")
    ap.add_argument("--wgcna-module-eigengenes", default=str(DEFAULT_WGCNA))
    ap.add_argument("--gene-priority", default=str(DEFAULT_MECH / "gene_axis_priority.tsv"))
    ap.add_argument("--cross-osdr-summary", default=str(DEFAULT_CROSS))
    ap.add_argument("--outdir", default="")
    ap.add_argument("--n-bootstrap", type=int, default=2000)
    ap.add_argument("--n-permutation", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260520)
    ap.add_argument("--max-gene-features", type=int, default=5000)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    state_path = resolve(args.state_space_scores)
    mechanism_path = resolve(args.mechanism_score_matrix)
    rtech_path = resolve(args.rtech)
    meta_path = resolve(args.meta)
    id_map_path = resolve(args.id_map)
    wgcna_path = resolve(args.wgcna_module_eigengenes)
    priority_path = resolve(args.gene_priority)
    cross_path = resolve(args.cross_osdr_summary)
    outdir = resolve(args.outdir) if args.outdir else state_path.parent / "lar_reversal"
    outdir.mkdir(parents=True, exist_ok=True)

    log("Loading state, expression, WGCNA, and network-priority inputs")
    state = pd.read_csv(state_path, sep="\t")
    mechanism = pd.read_csv(mechanism_path, sep="\t") if mechanism_path.exists() else state.copy()
    expression = read_expression(rtech_path)
    expression = expression[expression.index.astype(str).str.startswith("ENSMUSG")].copy()
    meta, sample_col = read_metadata(meta_path)
    priority = pd.read_csv(priority_path, sep="\t") if priority_path.exists() else pd.DataFrame()

    special_frame, gene_frame, special_coverage = build_target_gene_frame(
        expression,
        meta,
        id_map_path,
        sample_col,
        args.max_gene_features,
    )
    special_cols = [
        c for c in special_frame.columns
        if c not in {"Sample Name", "Age", "Arm", "EnvGroup", "condition", "study", "scenario"}
    ]
    state_vector_cols = [c for c in STATE_FEATURES if c in state.columns]
    state_component_cols = [c for c in STATE_COMPONENT_FEATURES if c in state.columns]
    mech_cols = [c for c in MECHANISM_FEATURES if c in state.columns]

    combined = state[state["study"].eq("RRRM-2")].copy()
    combined = combined.merge(
        special_frame[["Sample Name"] + special_cols],
        on="Sample Name",
        how="left",
    )
    combined.to_csv(outdir / "lar_reversal_sample_scores.tsv", sep="\t", index=False)

    vector_tables: list[pd.DataFrame] = []
    external_tables: list[pd.DataFrame] = []
    loo_tables: list[pd.DataFrame] = []
    feature_frames = [
        ("state_space", state, state_vector_cols, True),
        ("mechanism_scores", state, mech_cols, True),
        ("special_clock_s1p_rbm3", special_frame, special_cols, False),
        ("all_gene_expression", gene_frame, [c for c in gene_frame.columns if str(c).startswith("ENSMUSG")], False),
    ]

    if wgcna_path.exists():
        wgcna = read_wgcna(wgcna_path, state)
        wgcna_cols = [c for c in wgcna.columns if c.startswith("ME")]
        feature_frames.append(("wgcna_module_eigengenes", wgcna, wgcna_cols, False))

    log("Testing LAR attenuation/reversal/switch geometry")
    for feature_set, frame, cols, include_osd513 in feature_frames:
        if not cols:
            continue
        for age_scope in ("pooled", "YNG", "OLD"):
            summary, ext = reversal_summary_for_features(
                frame,
                cols,
                feature_set,
                age_scope=age_scope,
                n_bootstrap=args.n_bootstrap,
                n_permutation=args.n_permutation,
                rng=child_rng(args.seed, f"vector:{feature_set}:{age_scope}"),
                include_osd513=include_osd513 and age_scope == "pooled",
            )
            vector_tables.append(summary)
            if not ext.empty:
                external_tables.append(ext)
            if feature_set in {"state_space", "mechanism_scores", "special_clock_s1p_rbm3", "wgcna_module_eigengenes"}:
                rows = rrrm2_flt_gc_rows(frame, age_scope=age_scope)
                iss = arm_effect(rows, cols, "ISS-T", age_scope=age_scope)
                lar = arm_effect(rows, cols, "LAR", age_scope=age_scope)
                loo_tables.append(leave_one_feature_out(iss, lar, feature_set=feature_set, age_scope=age_scope))

    vector_tables.append(node_summary_rows(priority, feature_set="lioness_node_strength_summary", iss_col="N_ISS", lar_col="N_LAR"))
    node2vec_magnitude = node_summary_rows(priority, feature_set="node2vec_displacement_magnitude_summary", iss_col="D_ISS", lar_col="D_LAR")
    if not node2vec_magnitude.empty:
        node2vec_magnitude.to_csv(outdir / "lar_reversal_node2vec_magnitude_context.tsv", sep="\t", index=False)
    vector_summary = pd.concat([x for x in vector_tables if not x.empty], ignore_index=True, sort=False)
    vector_summary = add_q_values(vector_summary)
    vector_summary.to_csv(outdir / "lar_reversal_vector_summary.tsv", sep="\t", index=False)

    external = pd.concat(external_tables, ignore_index=True, sort=False) if external_tables else pd.DataFrame()
    cross_external = comparison_external_from_summary(cross_path)
    if not cross_external.empty:
        cross_external.to_csv(outdir / "lar_reversal_cross_osdr_context.tsv", sep="\t", index=False)
    external.to_csv(outdir / "lar_reversal_external_alignment.tsv", sep="\t", index=False)

    loo = pd.concat(loo_tables, ignore_index=True, sort=False) if loo_tables else pd.DataFrame()
    loo.to_csv(outdir / "lar_reversal_leave_one_feature_out.tsv", sep="\t", index=False)

    log("Computing component effects and arm-by-flight interactions")
    effect_tables = [
        component_effect_rows(state, state_component_cols, feature_family="state_space"),
        component_effect_rows(state, mech_cols, feature_family="mechanism_scores"),
        component_effect_rows(special_frame, special_cols, feature_family="special_clock_s1p_rbm3"),
    ]
    component_effects = pd.concat(effect_tables, ignore_index=True, sort=False)
    component_effects.to_csv(outdir / "lar_reversal_component_effects.tsv", sep="\t", index=False)

    interaction_tables = [
        interaction_table(
            state,
            state_component_cols,
            feature_family="state_space",
            n_bootstrap=args.n_bootstrap,
            n_permutation=args.n_permutation,
            rng=child_rng(args.seed, "interaction:state_space"),
        ),
        interaction_table(
            state,
            mech_cols,
            feature_family="mechanism_scores",
            n_bootstrap=args.n_bootstrap,
            n_permutation=args.n_permutation,
            rng=child_rng(args.seed, "interaction:mechanism_scores"),
        ),
        interaction_table(
            special_frame,
            special_cols,
            feature_family="special_clock_s1p_rbm3",
            n_bootstrap=args.n_bootstrap,
            n_permutation=args.n_permutation,
            rng=child_rng(args.seed, "interaction:special_clock_s1p_rbm3"),
        ),
    ]
    interactions = pd.concat(interaction_tables, ignore_index=True, sort=False)
    interactions.to_csv(outdir / "lar_reversal_pathway_interactions.tsv", sep="\t", index=False)

    log("Writing gene-level scatter and targeted mechanism tables")
    gene_scatter = gene_effect_scatter(
        expression,
        meta,
        str(id_map_path),
        sample_col=sample_col,
        priority=priority,
        panel_symbols=TARGETED_PANELS,
    )
    gene_scatter.to_csv(outdir / "lar_reversal_gene_scatter.tsv", sep="\t", index=False)

    switch_cols = sorted(set(state_vector_cols + mech_cols + special_cols))
    switch = mechanism_switch_decomposition(combined, switch_cols)
    switch.to_csv(outdir / "lar_mechanism_switch_scores.tsv", sep="\t", index=False)

    circ_cols = [c for c in ["clock_core", "per_cry", "bmal_clock", "dct_wnk_expression", "aldosterone_enac_expression", "Per2_expression", "dct_transport_component", "matrix_minus_dct"] if c in combined.columns]
    circ_effects = component_effect_rows(combined, circ_cols, feature_family="circadian_dct")
    circ_effects["analysis_family"] = "circadian_dct"
    circ_effects["row_type"] = "effect"
    circ_corrs = targeted_correlation_table(
        combined,
        [
            ("clock_core", "dct_transport_component"),
            ("per_cry", "dct_transport_component"),
            ("Per2_expression", "dct_transport_component"),
            ("Per2_expression", "matrix_minus_dct"),
            ("clock_core", "matrix_minus_dct"),
        ],
        "circadian_dct",
    )
    pd.concat([circ_effects, circ_corrs], ignore_index=True, sort=False).to_csv(
        outdir / "lar_circadian_dct_analysis.tsv",
        sep="\t",
        index=False,
    )

    s1p_cols = [c for c in ["s1p_s1pr3", "s1p_axis_expression", "S1pr3_expression", "matrix_component", "immune_context_component", "matrix_minus_dct"] if c in combined.columns]
    s1p_effects = component_effect_rows(combined, s1p_cols, feature_family="s1p_axis")
    s1p_effects["analysis_family"] = "s1p_axis"
    s1p_effects["row_type"] = "effect"
    s1p_corrs = targeted_correlation_table(
        combined,
        [
            ("s1p_axis_expression", "matrix_component"),
            ("S1pr3_expression", "matrix_component"),
            ("s1p_s1pr3", "immune_context_component"),
            ("S1pr3_expression", "matrix_minus_dct"),
        ],
        "s1p_axis",
    )
    pd.concat([s1p_effects, s1p_corrs], ignore_index=True, sort=False).to_csv(
        outdir / "lar_s1p_axis_analysis.tsv",
        sep="\t",
        index=False,
    )

    rbm3_cols = [c for c in ["Rbm3_expression", "rbm3_preservation_expression", "preservation_stress_component", "preservation_stress_response", "matrix_minus_dct", "dct_transport_component"] if c in combined.columns]
    rbm3_effects = component_effect_rows(combined, rbm3_cols, feature_family="rbm3_preservation")
    rbm3_effects["analysis_family"] = "rbm3_preservation"
    rbm3_effects["row_type"] = "effect"
    rbm3_corrs = targeted_correlation_table(
        combined,
        [
            ("Rbm3_expression", "preservation_stress_component"),
            ("Rbm3_expression", "preservation_stress_response"),
            ("Rbm3_expression", "matrix_minus_dct"),
            ("Rbm3_expression", "dct_transport_component"),
            ("rbm3_preservation_expression", "matrix_minus_dct"),
        ],
        "rbm3_preservation",
    )
    pd.concat([rbm3_effects, rbm3_corrs], ignore_index=True, sort=False).to_csv(
        outdir / "lar_rbm3_preservation_analysis.tsv",
        sep="\t",
        index=False,
    )

    special_coverage.to_csv(outdir / "lar_targeted_gene_set_coverage.tsv", sep="\t", index=False)
    model_calls = model_call_table(vector_summary, interactions, switch)
    model_calls.to_csv(outdir / "lar_reversal_model_calls.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "LAR reversal and mechanism-switch analysis",
        "timestamp": datetime.now().isoformat(),
        "git_commit": git_commit(),
        "parameters": {
            "n_bootstrap": int(args.n_bootstrap),
            "n_permutation": int(args.n_permutation),
            "seed": int(args.seed),
            "max_gene_features": int(args.max_gene_features),
            "models": {
                "A": "attenuation/recovery-to-baseline: LAR FLT-minus-GC near zero relative to ISS-T",
                "B": "true reversal: LAR vector is stably opposite to ISS-T",
                "C": "mechanism switch: LAR loads onto alternative clock/DCT, S1P, immune, or preservation-context axes",
            },
            "gene_level_claim_boundary": "gene and LIONESS/node2vec outputs are exploratory prioritization/context",
        },
        "inputs": {
            "state_space_scores": str(state_path),
            "state_space_scores_sha256": file_sha256(state_path),
            "mechanism_score_matrix": str(mechanism_path),
            "mechanism_score_matrix_sha256": file_sha256(mechanism_path),
            "rtech": str(rtech_path),
            "rtech_sha256": file_sha256(rtech_path),
            "meta": str(meta_path),
            "meta_sha256": file_sha256(meta_path),
            "id_map": str(id_map_path),
            "id_map_sha256": file_sha256(id_map_path),
            "wgcna_module_eigengenes": str(wgcna_path),
            "wgcna_module_eigengenes_sha256": file_sha256(wgcna_path),
            "gene_priority": str(priority_path),
            "gene_priority_sha256": file_sha256(priority_path),
            "cross_osdr_summary": str(cross_path),
            "cross_osdr_summary_sha256": file_sha256(cross_path),
        },
        "outputs": {
            "lar_reversal_vector_summary": str(outdir / "lar_reversal_vector_summary.tsv"),
            "lar_reversal_node2vec_magnitude_context": str(outdir / "lar_reversal_node2vec_magnitude_context.tsv"),
            "lar_reversal_component_effects": str(outdir / "lar_reversal_component_effects.tsv"),
            "lar_reversal_gene_scatter": str(outdir / "lar_reversal_gene_scatter.tsv"),
            "lar_reversal_pathway_interactions": str(outdir / "lar_reversal_pathway_interactions.tsv"),
            "lar_mechanism_switch_scores": str(outdir / "lar_mechanism_switch_scores.tsv"),
            "lar_circadian_dct_analysis": str(outdir / "lar_circadian_dct_analysis.tsv"),
            "lar_s1p_axis_analysis": str(outdir / "lar_s1p_axis_analysis.tsv"),
            "lar_rbm3_preservation_analysis": str(outdir / "lar_rbm3_preservation_analysis.tsv"),
            "lar_reversal_sample_scores": str(outdir / "lar_reversal_sample_scores.tsv"),
            "lar_reversal_manifest": str(outdir / "lar_reversal_manifest.json"),
        },
        "counts": {
            "state_rows": int(len(state)),
            "mechanism_rows": int(len(mechanism)),
            "rrrm2_special_rows": int(len(special_frame)),
            "gene_scatter_rows": int(len(gene_scatter)),
            "vector_summary_rows": int(len(vector_summary)),
        },
    }
    (outdir / "lar_reversal_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"Wrote LAR reversal outputs to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
