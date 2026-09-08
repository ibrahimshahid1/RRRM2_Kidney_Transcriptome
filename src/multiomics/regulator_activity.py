"""Regulator-activity prioritization (v10 trimmed plan)."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from src.common import bh_fdr

EPS = 1e-12


# Layer A -- KSEA (self-contained, no external dependency)

@dataclass(frozen=True)
class KseaResult:
    kinase: str
    n_substrates_quantified: int
    mean_substrate_effect: float
    z_score: float
    p_value: float


def _site_key(gene: str, position: object) -> str:
    return f"{str(gene).strip()}|{str(position).strip()}"


def load_kinase_substrate_net(path: str) -> pd.DataFrame:
    """Load a kinase-substrate table; requires kinase, substrate_gene, substrate_site columns."""
    net = pd.read_csv(path, sep="\t")
    required = {"kinase", "substrate_gene", "substrate_site"}
    missing = required - set(net.columns)
    if missing:
        raise ValueError(f"kinase-substrate net missing columns: {sorted(missing)}")
    net = net.copy()
    net["substrate_gene"] = net["substrate_gene"].astype(str).str.strip()
    net["substrate_site"] = net["substrate_site"].astype(str).str.strip()
    net["kinase"] = net["kinase"].astype(str).str.strip()
    return net


def ksea(
    sites: pd.DataFrame,
    ks_net: pd.DataFrame,
    *,
    effect_col: str = "phospho_effect",
    gene_col: str = "gene_symbol",
    site_col: str = "site_position",
    min_substrates: int = 3,
) -> pd.DataFrame:
    """Casado (2013) KSEA z-score; positive z = substrates more phosphorylated in flight."""
    s = sites[[gene_col, site_col, effect_col]].copy()
    s = s.replace([np.inf, -np.inf], np.nan).dropna(subset=[effect_col])
    s["__key__"] = [_site_key(g, p) for g, p in zip(s[gene_col], s[site_col])]
    # collapse duplicate site rows (e.g. composite vs single) by mean
    site_effect = s.groupby("__key__")[effect_col].mean()
    all_effects = site_effect.to_numpy(dtype=float)
    bg_mean = float(np.mean(all_effects))
    bg_sd = float(np.std(all_effects, ddof=1))
    if bg_sd < EPS:
        raise ValueError("background phosphosite effect standard deviation is ~0")

    rows: list[dict[str, object]] = []
    for kinase, grp in ks_net.groupby("kinase"):
        keys = {_site_key(g, p) for g, p in zip(grp["substrate_gene"], grp["substrate_site"])}
        matched = site_effect.reindex(sorted(keys)).dropna()
        n = int(len(matched))
        if n < min_substrates:
            rows.append({
                "kinase": kinase,
                "n_substrates_quantified": n,
                "mean_substrate_effect": float(matched.mean()) if n else np.nan,
                "z_score": np.nan,
                "p_value": np.nan,
                "status": "insufficient_substrates",
            })
            continue
        m = float(matched.mean())
        z = (m - bg_mean) * np.sqrt(n) / bg_sd
        p = float(2.0 * stats.norm.sf(abs(z)))
        rows.append({
            "kinase": kinase,
            "n_substrates_quantified": n,
            "mean_substrate_effect": m,
            "z_score": float(z),
            "p_value": p,
            "status": "scored",
        })
    out = pd.DataFrame(rows)
    scored = out["p_value"].notna()
    out["q_value"] = np.nan
    if scored.any():
        out.loc[scored, "q_value"] = bh_fdr(out.loc[scored, "p_value"].to_numpy(dtype=float))
    out["direction"] = np.where(
        out["z_score"] > 0, "inferred_activity_up",
        np.where(out["z_score"] < 0, "inferred_activity_down", "flat"),
    )
    out["background_n_sites"] = int(len(site_effect))
    return out.sort_values("z_score", na_position="last").reset_index(drop=True)


def ksea_positive_control_passes(ksea_table: pd.DataFrame, *, control_kinases: Sequence[str]) -> bool:
    """Evaluate a predeclared down-direction control; the caller must first enforce provenance."""
    sub = ksea_table[ksea_table["kinase"].isin(list(control_kinases))]
    sub = sub[sub["status"].eq("scored")]
    if sub.empty:
        return False
    return bool((sub["z_score"] < 0).all() and (sub["p_value"] < 0.05).all())


# Layer B -- TF / pathway activity inference (decoupler ULM)

def run_ulm_activity(
    gene_effects: pd.DataFrame,
    net: pd.DataFrame,
    *,
    tmin: int = 5,
) -> pd.DataFrame:
    """decoupler ULM activity from a contrasts x genes matrix and source/target/weight prior."""
    import decoupler as dc  # lazy import

    res = dc.mt.ulm(data=gene_effects, net=net, tmin=tmin)
    # decoupler 2.x returns (estimate, pvalue) DataFrames or stores on obj
    if isinstance(res, tuple):
        estimate, pvalue = res[0], res[1]
    else:  # pragma: no cover - depends on decoupler return convention
        estimate, pvalue = res, None
    long = estimate.reset_index().melt(
        id_vars=estimate.index.name or "index",
        var_name="source",
        value_name="activity_score",
    )
    long = long.rename(columns={estimate.index.name or "index": "contrast"})
    if pvalue is not None:
        pv = pvalue.reset_index().melt(
            id_vars=pvalue.index.name or "index",
            var_name="source",
            value_name="p_value",
        ).rename(columns={pvalue.index.name or "index": "contrast"})
        long = long.merge(pv, on=["contrast", "source"], how="left")
        mask = long["p_value"].notna()
        long["q_value"] = np.nan
        if mask.any():
            long.loc[mask, "q_value"] = bh_fdr(long.loc[mask, "p_value"].to_numpy(dtype=float))
    return long


def recurrence_class(activity_by_cohort: Mapping[str, float], *, threshold: float = 1.0) -> str:
    """Classify cross-cohort recurrence from per-cohort activity; requires a consistent sign."""
    vals = np.array([v for v in activity_by_cohort.values() if np.isfinite(v)], dtype=float)
    if len(vals) < 2:
        return "insufficient_cohorts"
    strong = vals[np.abs(vals) >= threshold]
    if len(strong) < 2:
        return "not_recurrent"
    if np.all(strong > 0):
        return "recurrent_up"
    if np.all(strong < 0):
        return "recurrent_down"
    return "mixed_direction"


# Integration -- evidence grading

EVIDENCE_GRADES = (
    "activity_anchor",            # phospho-supported, positive control
    "candidate_upstream_organizer",
    "candidate_context_axis",
    "technical_or_biological_caution",
    "negative_boundary",          # tested and failed (e.g. network translation)
)


def grade_axis(
    *,
    rna_recurrence: str,
    phospho_support: bool,
    protein_boundary_null: bool,
    is_positive_control: bool = False,
) -> str:
    """Assign a conservative evidence grade to a candidate axis."""
    if is_positive_control and phospho_support:
        return "activity_anchor"
    if phospho_support and rna_recurrence in {"recurrent_up", "recurrent_down"}:
        return "candidate_upstream_organizer"
    if rna_recurrence in {"recurrent_up", "recurrent_down"}:
        return "candidate_context_axis"
    if rna_recurrence == "mixed_direction":
        return "technical_or_biological_caution"
    return "negative_boundary"
