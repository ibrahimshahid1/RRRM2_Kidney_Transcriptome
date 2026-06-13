"""Mechanism-axis scoring and exploratory network-neighborhood prioritization."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from src.common import bh_fdr, load_id_map


EPS = 1e-12


@dataclass(frozen=True)
class GeneSetResolution:
    """Resolved gene-set membership against an expression universe."""

    name: str
    role: str
    description: str
    requested_symbols: tuple[str, ...]
    matched_gene_ids: tuple[str, ...]
    matched_symbols: tuple[str, ...]
    scored_gene_ids: tuple[str, ...]
    scored_symbols: tuple[str, ...]
    dropped_overlap_symbols: tuple[str, ...]

    @property
    def n_requested(self) -> int:
        return len(self.requested_symbols)

    @property
    def n_matched(self) -> int:
        return len(self.matched_gene_ids)

    @property
    def n_scored(self) -> int:
        return len(self.scored_gene_ids)


def safe_z(values: Sequence[float]) -> np.ndarray:
    """Return finite z-scores; constant vectors become zeros."""
    arr = np.asarray(values, dtype=float)
    out = np.zeros_like(arr, dtype=float)
    finite = np.isfinite(arr)
    if finite.sum() < 2:
        return out
    mu = float(arr[finite].mean())
    sd = float(arr[finite].std(ddof=1))
    if sd < EPS:
        return out
    out[finite] = (arr[finite] - mu) / sd
    return out


def unit_vector(series: pd.Series) -> pd.Series:
    """Normalize a vector to unit Euclidean length while preserving labels."""
    vals = series.astype(float).replace([np.inf, -np.inf], np.nan)
    norm = float(np.linalg.norm(vals.dropna().to_numpy(dtype=float)))
    if norm < EPS:
        return vals * 0.0
    return vals / norm


def cosine_distance_rows(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """1 - cosine similarity row-wise."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    denom = (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)) + EPS
    sim = np.sum(a * b, axis=1) / denom
    return 1.0 - np.clip(sim, -1.0, 1.0)


def flatten_gene_symbols(values: object) -> list[str]:
    """Flatten common YAML gene-list shapes into a list of symbols/IDs."""
    out: list[str] = []
    if values is None:
        return out
    if isinstance(values, str):
        return [values]
    if isinstance(values, list):
        for item in values:
            out.extend(flatten_gene_symbols(item))
    elif isinstance(values, dict):
        if "genes" in values:
            out.extend(flatten_gene_symbols(values.get("genes")))
        else:
            for item in values.values():
                out.extend(flatten_gene_symbols(item))
    return [x.strip() for x in out if str(x).strip()]


def load_gene_set_config(path: Path) -> dict[str, dict[str, object]]:
    """Load gene-set YAML, dropping non-set parameter blocks."""
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    out: dict[str, dict[str, object]] = {}
    for name, payload in raw.items():
        if isinstance(payload, dict) and "genes" in payload:
            out[str(name)] = payload
    return out


def resolve_gene_sets(
    gene_sets: Mapping[str, Mapping[str, object]],
    id_map_path: Path,
    expression_genes: set[str],
    *,
    min_genes: int = 3,
    deoverlap_gene_ids: set[str] | None = None,
    protected_symbols: set[str] | None = None,
) -> tuple[dict[str, GeneSetResolution], pd.DataFrame]:
    """Resolve symbol/Ensembl YAML sets to scored Ensembl IDs.

    ``deoverlap_gene_ids`` is used for candidate-upstream mechanism scores so
    axis-defining genes can be reported and removed.  Symbols listed in a set's
    ``protected_genes`` field remain scored even if they overlap the axis.
    """
    id_map = load_id_map(id_map_path)
    symbol_to_ids: dict[str, set[str]] = {}
    id_to_symbol: dict[str, str] = {}
    for row in id_map.itertuples():
        gid = str(row.ensembl_gene_id)
        sym = str(row.mgi_symbol)
        id_to_symbol[gid] = sym
        symbol_to_ids.setdefault(sym.casefold(), set()).add(gid)

    deoverlap_gene_ids = deoverlap_gene_ids or set()
    protected_symbols = {s.casefold() for s in (protected_symbols or set())}
    resolved: dict[str, GeneSetResolution] = {}
    coverage_rows: list[dict[str, object]] = []

    for name, payload in gene_sets.items():
        symbols = tuple(dict.fromkeys(flatten_gene_symbols(payload.get("genes"))))
        role = str(payload.get("role", "unspecified"))
        description = str(payload.get("description", ""))
        local_protected = set(protected_symbols)
        local_protected |= {str(x).casefold() for x in flatten_gene_symbols(payload.get("protected_genes"))}

        matched_ids: set[str] = set()
        requested_to_ids: dict[str, set[str]] = {}
        for sym in symbols:
            ids = set()
            if sym in expression_genes:
                ids.add(sym)
            ids |= symbol_to_ids.get(sym.casefold(), set()) & expression_genes
            requested_to_ids[sym] = ids
            matched_ids |= ids

        dropped_ids: set[str] = set()
        scored_ids: set[str] = set()
        for gid in matched_ids:
            sym = id_to_symbol.get(gid, gid)
            protected = sym.casefold() in local_protected
            if gid in deoverlap_gene_ids and not protected:
                dropped_ids.add(gid)
            else:
                scored_ids.add(gid)

        if len(scored_ids) < min_genes:
            scored_ids = set(matched_ids)
            dropped_ids = set()

        matched_ids_sorted = tuple(sorted(matched_ids))
        scored_ids_sorted = tuple(sorted(scored_ids))
        dropped_symbols = tuple(sorted({id_to_symbol.get(g, g) for g in dropped_ids}))
        res = GeneSetResolution(
            name=name,
            role=role,
            description=description,
            requested_symbols=symbols,
            matched_gene_ids=matched_ids_sorted,
            matched_symbols=tuple(sorted({id_to_symbol.get(g, g) for g in matched_ids_sorted})),
            scored_gene_ids=scored_ids_sorted,
            scored_symbols=tuple(sorted({id_to_symbol.get(g, g) for g in scored_ids_sorted})),
            dropped_overlap_symbols=dropped_symbols,
        )
        if res.n_scored >= min_genes:
            resolved[name] = res
        coverage_rows.append({
            "gene_set": name,
            "role": role,
            "description": description,
            "n_requested": res.n_requested,
            "n_matched": res.n_matched,
            "n_scored": res.n_scored,
            "coverage_requested": res.n_matched / max(res.n_requested, 1),
            "n_dropped_axis_overlap": len(dropped_symbols),
            "dropped_axis_overlap_symbols": ",".join(dropped_symbols),
            "scored_symbols": ",".join(res.scored_symbols),
            "dropped_for_min_genes": bool(res.n_scored < min_genes),
        })

    return resolved, pd.DataFrame(coverage_rows)


def score_gene_sets(
    expression: pd.DataFrame,
    resolved_sets: Mapping[str, GeneSetResolution],
) -> pd.DataFrame:
    """Score gene sets as sample x set mean gene-wise z-expression."""
    scores: dict[str, np.ndarray] = {}
    expr = expression.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    for name, res in resolved_sets.items():
        genes = [g for g in res.scored_gene_ids if g in expr.index]
        if len(genes) < 3:
            continue
        mat = expr.loc[genes].to_numpy(dtype=float)
        mu = np.nanmean(mat, axis=1, keepdims=True)
        sd = np.nanstd(mat, axis=1, ddof=1, keepdims=True)
        sd[~np.isfinite(sd) | (sd < EPS)] = 1.0
        z = (mat - mu) / sd
        scores[name] = np.nanmean(z, axis=0)
    return pd.DataFrame(scores, index=expr.columns)


def flight_minus_control(
    scores: pd.DataFrame,
    sample_rows: pd.DataFrame,
    *,
    sample_col: str = "Sample Name",
    condition_col: str = "condition",
    flt_label: str = "FLT",
    gc_label: str = "GC",
) -> pd.Series:
    """Mean FLT-minus-control vector for sample-level feature scores."""
    rows = sample_rows[sample_rows[sample_col].isin(scores.index)].copy()
    flt = rows.loc[rows[condition_col].eq(flt_label), sample_col].astype(str).tolist()
    gc = rows.loc[rows[condition_col].eq(gc_label), sample_col].astype(str).tolist()
    if len(flt) < 2 or len(gc) < 2:
        raise ValueError(f"Need at least two FLT and GC samples; found FLT={len(flt)}, GC={len(gc)}")
    return scores.loc[flt].mean(axis=0) - scores.loc[gc].mean(axis=0)


def rrrm2_arm_rows(meta: pd.DataFrame, arm: str) -> pd.DataFrame:
    """Return RRRM-2 FLT/GC rows for one arm with a generic condition column."""
    rows = meta[
        meta["Arm"].astype(str).eq(arm)
        & meta["EnvGroup"].astype(str).isin(["FLT", "GC"])
    ].copy()
    rows["condition"] = rows["EnvGroup"].astype(str)
    return rows


def stratified_arm_flight_effect(
    scores: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
    arm: str,
    feature: str,
    *,
    strata_col: str = "Age",
) -> float:
    """Balanced FLT-minus-GC effect averaged across age strata."""
    rows = rrrm2_arm_rows(meta, arm)
    diffs: list[float] = []
    for _, sub in rows.groupby(strata_col):
        flt = sub.loc[sub["condition"].eq("FLT"), sample_col].astype(str).tolist()
        gc = sub.loc[sub["condition"].eq("GC"), sample_col].astype(str).tolist()
        flt = [s for s in flt if s in scores.index]
        gc = [s for s in gc if s in scores.index]
        if flt and gc:
            diffs.append(float(scores.loc[flt, feature].mean() - scores.loc[gc, feature].mean()))
    return float(np.mean(diffs)) if diffs else float("nan")


def define_recurrent_axis(
    rrrm2_iss: pd.Series,
    osd513: pd.Series,
    *,
    primary_transport_set: str = "dct_ncc_wnk",
    sensitivity_transport_set: str = "ion_transport",
    positive_anchor_sets: Sequence[str] = ("ecm_remodeling", "fibrosis"),
) -> tuple[pd.DataFrame, pd.Series]:
    """Define a unit-normalized recurrent ISS-T/OSD-513 pathway axis."""
    features = sorted(set(rrrm2_iss.index) & set(osd513.index))
    a = unit_vector(rrrm2_iss.loc[features])
    b = unit_vector(osd513.loc[features])
    raw = (a + b) / 2.0
    concordant = (rrrm2_iss.loc[features] * osd513.loc[features]) > 0

    orientation = 1.0
    anchor_vals = [raw.get(x, np.nan) for x in positive_anchor_sets if x in raw.index]
    anchor_vals = [float(x) for x in anchor_vals if np.isfinite(x)]
    if anchor_vals and float(np.nanmean(anchor_vals)) < 0:
        orientation = -1.0
    raw_oriented = raw * orientation

    rows = []
    included_weights: dict[str, float] = {}
    for feat in features:
        include = bool(concordant.loc[feat])
        rationale = "concordant_rrrm2_osd513"
        if feat == sensitivity_transport_set:
            include = False
            rationale = "sensitivity_only_broad_transport"
        if not np.isfinite(raw_oriented.loc[feat]) or abs(float(raw_oriented.loc[feat])) < EPS:
            include = False
            rationale = "zero_or_nonfinite_weight"
        if include:
            included_weights[feat] = float(raw_oriented.loc[feat])
        rows.append({
            "pathway": feat,
            "rrrm2_iss_effect": float(rrrm2_iss.loc[feat]),
            "osd513_effect": float(osd513.loc[feat]),
            "rrrm2_unit": float(a.loc[feat]),
            "osd513_unit": float(b.loc[feat]),
            "concordant": bool(concordant.loc[feat]),
            "raw_weight_oriented": float(raw_oriented.loc[feat]),
            "included_primary_axis": include,
            "rationale": rationale,
        })

    if not included_weights:
        raise ValueError("No concordant RRRM-2/OSD-513 pathways available for recurrent axis")
    weights = pd.Series(included_weights, dtype=float)
    weights = weights / (float(np.linalg.norm(weights.to_numpy())) + EPS)
    table = pd.DataFrame(rows)
    table["weight"] = table["pathway"].map(weights).fillna(0.0)
    table["axis_orientation"] = orientation
    table["primary_transport_set"] = primary_transport_set
    table["sensitivity_transport_set"] = sensitivity_transport_set
    return table.sort_values("pathway").reset_index(drop=True), weights.sort_index()


def remodeling_score(score_matrix: pd.DataFrame, weights: pd.Series) -> pd.Series:
    """Project sample x pathway scores onto a recurrent pathway axis."""
    common = [p for p in weights.index if p in score_matrix.columns]
    if not common:
        raise ValueError("No overlapping pathway scores for remodeling projection")
    z = score_matrix[common].apply(safe_z, axis=0, result_type="broadcast")
    return z.dot(weights.loc[common])


def select_osd253_scenario(meta: pd.DataFrame, scenario: str) -> pd.DataFrame:
    """Metadata-aware OSD-253 FLT/control selection for both strains."""
    scenario_norm = scenario.replace("-", "_").casefold()
    aliases = {
        "original_gc_blue": "original_gc_blue_light",
        "original_gc_blue_light": "original_gc_blue_light",
        "rerun_white_gc": "rerun_gc_white_light",
        "rerun_gc_white": "rerun_gc_white_light",
        "rerun_gc_white_light": "rerun_gc_white_light",
    }
    scenario_key = aliases.get(scenario_norm)
    if scenario_key is None:
        raise ValueError(f"Unknown OSD-253 scenario: {scenario}")

    base = meta[
        meta["Factor Value[Strain]"].astype(str).isin(["C57BL/6J", "C3H/HeJ"])
        & meta["Factor Value[Spaceflight]"].astype(str).isin(["Space Flight", "Ground Control", "Ground Control Rerun"])
        & meta["Factor Value[Duration]"].astype(str).isin(["~25 day", "~75 day", "~25", "~75"])
    ].copy()
    treatment = base["Factor Value[Treatment]"].astype(str).str.casefold()
    space = base["Factor Value[Spaceflight]"].astype(str)
    if scenario_key == "original_gc_blue_light":
        keep = (
            (space.eq("Space Flight") & treatment.eq("white light"))
            | (space.eq("Ground Control") & treatment.eq("blue light"))
        )
        base = base[keep].copy()
        base["condition"] = np.where(base["Factor Value[Spaceflight]"].eq("Space Flight"), "FLT", "GC")
        base["control_definition"] = "original_gc_blue_light"
    else:
        keep = (
            (space.eq("Space Flight") & treatment.eq("white light"))
            | (space.eq("Ground Control Rerun") & treatment.eq("white light"))
        )
        base = base[keep].copy()
        base["condition"] = np.where(base["Factor Value[Spaceflight]"].eq("Space Flight"), "FLT", "GC")
        base["control_definition"] = "rerun_gc_white_light"
    base["strain"] = base["Factor Value[Strain]"].astype(str)
    base["duration"] = base["Factor Value[Duration]"].astype(str)
    return base


def strain_duration_effect(df: pd.DataFrame, score_col: str, strain: str) -> float:
    """Average FLT-control effect across OSD-253 durations for one strain."""
    diffs: list[float] = []
    for _, sub in df[df["strain"].eq(strain)].groupby("duration"):
        flt = sub.loc[sub["condition"].eq("FLT"), score_col]
        gc = sub.loc[sub["condition"].eq("GC"), score_col]
        if len(flt) and len(gc):
            diffs.append(float(flt.mean() - gc.mean()))
    return float(np.mean(diffs)) if diffs else float("nan")


def osd253_observed_stats(df: pd.DataFrame, score_col: str) -> dict[str, float]:
    """C57, C3H, and strain-delta effects for one OSD-253 score."""
    c57 = strain_duration_effect(df, score_col, "C57BL/6J")
    c3h = strain_duration_effect(df, score_col, "C3H/HeJ")
    return {
        "c57_effect": c57,
        "c3h_effect": c3h,
        "c3h_minus_c57": c3h - c57,
    }


def permute_osd253_conditions(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Permute FLT/control labels within Strain x Duration strata."""
    out = df.copy()
    for _, idx in out.groupby(["strain", "duration"]).groups.items():
        labels = out.loc[idx, "condition"].to_numpy(copy=True)
        rng.shuffle(labels)
        out.loc[idx, "condition"] = labels
    return out


def bootstrap_osd253_stats(
    df: pd.DataFrame,
    score_col: str,
    rng: np.random.Generator,
    n_bootstrap: int,
) -> pd.DataFrame:
    """Bootstrap OSD-253 score effects within Strain x Duration x Condition cells."""
    rows: list[dict[str, float | int]] = []
    groups = list(df.groupby(["strain", "duration", "condition"]).groups.items())
    for b in range(n_bootstrap):
        parts = []
        for _, idx in groups:
            idx = list(idx)
            chosen = rng.choice(idx, size=len(idx), replace=True)
            parts.append(df.loc[chosen])
        sampled = pd.concat(parts, ignore_index=True)
        rows.append({"iteration": b, **osd253_observed_stats(sampled, score_col)})
    return pd.DataFrame(rows)


def classify_tlr4_discriminator(row: Mapping[str, object]) -> str:
    """Translate OSD-253 strain attenuation into an interpretation label."""
    c57 = float(row.get("c57_effect", np.nan))
    c3h = float(row.get("c3h_effect", np.nan))
    q_delta = float(row.get("q_strain_delta", np.nan))
    if not np.isfinite(c57) or abs(c57) < EPS or not np.isfinite(c3h):
        return "inconclusive"
    attenuation = 1.0 - (c3h / c57)
    if np.isfinite(q_delta) and q_delta <= 0.10 and attenuation >= 0.50:
        return "tlr4_supported"
    if np.sign(c57) == np.sign(c3h) and abs(c3h) >= 0.70 * abs(c57):
        return "tlr4_not_required"
    return "inconclusive"


def osd253_discriminator_table(
    score_table: pd.DataFrame,
    meta: pd.DataFrame,
    *,
    scenarios: Sequence[str],
    score_cols: Sequence[str],
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run OSD-253 strain-discriminator tests for selected score columns."""
    summary_rows: list[dict[str, object]] = []
    boot_rows: list[pd.DataFrame] = []
    for scenario in scenarios:
        selected = select_osd253_scenario(meta, scenario)
        merged = selected.merge(score_table, left_on="Sample Name", right_index=True, how="inner")
        for score_col in score_cols:
            if score_col not in merged.columns:
                continue
            obs = osd253_observed_stats(merged, score_col)
            boot = bootstrap_osd253_stats(merged, score_col, rng, n_bootstrap)
            boot.insert(0, "score_name", score_col)
            boot.insert(0, "scenario", scenario)
            boot_rows.append(boot)
            perm_counts = {"c57_effect": 0, "c3h_effect": 0, "c3h_minus_c57": 0}
            for _ in range(n_bootstrap):
                stat = osd253_observed_stats(permute_osd253_conditions(merged, rng), score_col)
                for key in perm_counts:
                    perm_counts[key] += int(abs(stat[key]) >= abs(obs[key]))
            row: dict[str, object] = {
                "scenario": scenario,
                "score_name": score_col,
                **obs,
                "c57_ci_low": float(boot["c57_effect"].quantile(0.025)),
                "c57_ci_high": float(boot["c57_effect"].quantile(0.975)),
                "c3h_ci_low": float(boot["c3h_effect"].quantile(0.025)),
                "c3h_ci_high": float(boot["c3h_effect"].quantile(0.975)),
                "delta_ci_low": float(boot["c3h_minus_c57"].quantile(0.025)),
                "delta_ci_high": float(boot["c3h_minus_c57"].quantile(0.975)),
                "p_c57_effect": (perm_counts["c57_effect"] + 1) / (n_bootstrap + 1),
                "p_c3h_effect": (perm_counts["c3h_effect"] + 1) / (n_bootstrap + 1),
                "p_strain_delta": (perm_counts["c3h_minus_c57"] + 1) / (n_bootstrap + 1),
                "n_samples": int(len(merged)),
                "n_bootstrap": int(n_bootstrap),
            }
            row["attenuation_fraction_vs_c57"] = (
                1.0 - (obs["c3h_effect"] / obs["c57_effect"])
                if np.isfinite(obs["c57_effect"]) and abs(obs["c57_effect"]) > EPS
                else np.nan
            )
            summary_rows.append(row)
    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        for scenario, idx in summary.groupby("scenario").groups.items():
            loc = list(idx)
            summary.loc[loc, "q_c57_effect"] = bh_fdr(summary.loc[loc, "p_c57_effect"].to_numpy())
            summary.loc[loc, "q_c3h_effect"] = bh_fdr(summary.loc[loc, "p_c3h_effect"].to_numpy())
            summary.loc[loc, "q_strain_delta"] = bh_fdr(summary.loc[loc, "p_strain_delta"].to_numpy())
        summary["interpretation"] = summary.apply(classify_tlr4_discriminator, axis=1)
    boots = pd.concat(boot_rows, ignore_index=True) if boot_rows else pd.DataFrame()
    return summary, boots


def aggregate_lioness_nodes(
    weights: np.ndarray,
    edges: pd.DataFrame,
    genes: Sequence[str],
) -> pd.DataFrame:
    """Aggregate sample x edge LIONESS weights to absolute incident node strength."""
    genes = list(dict.fromkeys(str(g) for g in genes))
    pos = {g: i for i, g in enumerate(genes)}
    values = np.abs(np.asarray(weights, dtype=np.float32))
    out = np.zeros((values.shape[0], len(genes)), dtype=np.float32)
    counts = np.zeros(len(genes), dtype=np.float32)
    gi = edges["gene_i"].astype(str).to_numpy()
    gj = edges["gene_j"].astype(str).to_numpy()
    for e, (a, b) in enumerate(zip(gi, gj)):
        ia = pos.get(a)
        ib = pos.get(b)
        if ia is not None:
            out[:, ia] += values[:, e]
            counts[ia] += 1.0
        if ib is not None:
            out[:, ib] += values[:, e]
            counts[ib] += 1.0
    out = np.divide(out, counts, out=np.zeros_like(out), where=counts > 0)
    return pd.DataFrame(out, columns=genes)


def age_balanced_feature_effect(
    feature_matrix: pd.DataFrame,
    meta: pd.DataFrame,
    sample_col: str,
    arm: str,
) -> pd.Series:
    """Age-balanced arm-specific FLT-minus-GC effect for sample x gene features."""
    rows = rrrm2_arm_rows(meta, arm)
    diffs: list[pd.Series] = []
    for _, sub in rows.groupby("Age"):
        flt = [s for s in sub.loc[sub["condition"].eq("FLT"), sample_col].astype(str) if s in feature_matrix.index]
        gc = [s for s in sub.loc[sub["condition"].eq("GC"), sample_col].astype(str) if s in feature_matrix.index]
        if flt and gc:
            diffs.append(feature_matrix.loc[flt].mean(axis=0) - feature_matrix.loc[gc].mean(axis=0))
    if not diffs:
        return pd.Series(np.nan, index=feature_matrix.columns)
    return pd.concat(diffs, axis=1).mean(axis=1)


def specificity_score(iss_values: pd.Series, lar_values: pd.Series, *, absolute: bool = False) -> pd.Series:
    """ISS-T-specific score: z(ISS-T) - z(LAR), aligned by feature."""
    common = sorted(set(iss_values.index) & set(lar_values.index))
    iss = iss_values.loc[common].astype(float)
    lar = lar_values.loc[common].astype(float)
    if absolute:
        iss = iss.abs()
        lar = lar.abs()
    return pd.Series(safe_z(iss) - safe_z(lar), index=common)


def composite_gene_priority(
    embedding_iss: pd.Series,
    embedding_lar: pd.Series,
    node_iss: pd.Series,
    node_lar: pd.Series,
    osd513_support: pd.Series | None = None,
    iss_abs_logfc: pd.Series | None = None,
) -> pd.DataFrame:
    """Build exploratory composite and silent-shifter gene priority scores."""
    genes = sorted(set(embedding_iss.index) | set(embedding_lar.index) | set(node_iss.index) | set(node_lar.index))
    df = pd.DataFrame(index=genes)
    df["D_ISS"] = embedding_iss.reindex(genes)
    df["D_LAR"] = embedding_lar.reindex(genes)
    df["embedding_specificity"] = specificity_score(
        df["D_ISS"].fillna(df["D_ISS"].median()),
        df["D_LAR"].fillna(df["D_LAR"].median()),
    ).reindex(genes)
    df["N_ISS"] = node_iss.reindex(genes)
    df["N_LAR"] = node_lar.reindex(genes)
    df["lioness_node_specificity"] = specificity_score(
        df["N_ISS"].fillna(0.0),
        df["N_LAR"].fillna(0.0),
        absolute=True,
    ).reindex(genes)
    if osd513_support is not None and not osd513_support.empty:
        df["OSD513_support"] = osd513_support.reindex(genes)
        osd_term = 0.5 * safe_z(df["OSD513_support"].fillna(0.0))
    else:
        df["OSD513_support"] = np.nan
        osd_term = np.zeros(len(df), dtype=float)
    df["composite_score"] = (
        df["embedding_specificity"].fillna(0.0).to_numpy(dtype=float)
        + df["lioness_node_specificity"].fillna(0.0).to_numpy(dtype=float)
        + osd_term
    )
    if iss_abs_logfc is not None and not iss_abs_logfc.empty:
        df["iss_t_abs_logfc"] = iss_abs_logfc.reindex(genes)
        df["silent_composite_score"] = df["composite_score"] - safe_z(df["iss_t_abs_logfc"].fillna(0.0))
    else:
        df["iss_t_abs_logfc"] = np.nan
        df["silent_composite_score"] = np.nan
    return df.reset_index(names="gene").sort_values("composite_score", ascending=False, kind="mergesort")


def fisher_topk_enrichment(
    ranked: pd.DataFrame,
    resolved_sets: Mapping[str, GeneSetResolution],
    *,
    score_col: str,
    top_ks: Sequence[int] = (25, 50, 100),
) -> pd.DataFrame:
    """Top-k Fisher enrichment for ranked candidate genes."""
    universe = set(ranked["gene"].astype(str))
    rows: list[dict[str, object]] = []
    ordered = ranked.sort_values(score_col, ascending=False, kind="mergesort")
    for k in top_ks:
        top = set(ordered.head(k)["gene"].astype(str))
        for name, res in resolved_sets.items():
            gs = set(res.scored_gene_ids) & universe
            if not gs:
                continue
            a = len(top & gs)
            b = len(top - gs)
            c = len(gs - top)
            d = len(universe - top - gs)
            _, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "score_col": score_col,
                "top_k": int(k),
                "gene_set": name,
                "role": res.role,
                "overlap": int(a),
                "top_k_size": int(len(top)),
                "set_size": int(len(gs)),
                "universe_size": int(len(universe)),
                "fold_enrichment": (a / max(len(top), 1)) / (len(gs) / max(len(universe), 1)),
                "p_value": float(p),
                "overlap_genes": ",".join(sorted(top & gs)),
            })
    out = pd.DataFrame(rows)
    if not out.empty:
        for (score_col_value, k), idx in out.groupby(["score_col", "top_k"]).groups.items():
            loc = list(idx)
            out.loc[loc, "q_bh"] = bh_fdr(out.loc[loc, "p_value"].to_numpy())
    return out


def proportional_top_ks(
    universe_size: int,
    *,
    fractions: Sequence[float] = (0.01, 0.02, 0.05, 0.10),
) -> list[int]:
    """Return proportional top-k thresholds for a ranked universe."""
    ks: set[int] = set()
    for frac in fractions:
        if frac <= 0:
            continue
        ks.add(max(1, int(round(universe_size * float(frac)))))
    return sorted(k for k in ks if k <= universe_size)


def preranked_enrichment(
    ranked: pd.DataFrame,
    resolved_sets: Mapping[str, GeneSetResolution],
    *,
    score_col: str,
) -> pd.DataFrame:
    """Full ranked-list enrichment using one-sided Mann-Whitney tests."""
    scores = ranked.set_index("gene")[score_col].astype(float).replace([np.inf, -np.inf], np.nan).dropna()
    universe = set(scores.index.astype(str))
    rows: list[dict[str, object]] = []
    for name, res in resolved_sets.items():
        gs = set(res.scored_gene_ids) & universe
        if len(gs) < 3 or len(universe - gs) < 3:
            continue
        in_scores = scores.loc[sorted(gs)]
        out_scores = scores.loc[sorted(universe - gs)]
        stat, p = stats.mannwhitneyu(in_scores, out_scores, alternative="greater", method="asymptotic")
        rows.append({
            "score_col": score_col,
            "gene_set": name,
            "role": res.role,
            "set_size": int(len(gs)),
            "universe_size": int(len(universe)),
            "median_in_set": float(in_scores.median()),
            "median_out_set": float(out_scores.median()),
            "u_statistic": float(stat),
            "p_value": float(p),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_bh"] = bh_fdr(out["p_value"].to_numpy())
    return out
