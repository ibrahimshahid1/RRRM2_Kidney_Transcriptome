"""Repair B -- cross-cohort recurrence meta-analysis.

Turns the descriptive cosine recurrence (0.87 / 0.64 / -0.51 across three
cohorts) into a precision-weighted random-effects meta-analysis with per-gene
FDR and heterogeneity, across *five* on-disk mouse-kidney spaceflight cohorts
(OSD-102, OSD-163, OSD-253, OSD-462, OSD-513; two strains).

The contract the manuscript leans on is *sign-faithful, precision-weighted
pooling*: a gene coherently suppressed across cohorts must yield a negative
pooled effect with a small meta-p and (when the cohorts agree) a low I^2; the
per-cohort effects feeding the pool are honest flight-vs-ground contrasts on the
VST scale, computed in Python (Welch t; no R). The DCT/NCC-WNK transport program
and the ECM/matrix program are then re-scored with the *pooled* statistic in
place of cosine, and leave-one-cohort-out reports stability.

Design provenance
-----------------
Per-cohort flight/ground membership is read from the GeneLab VST column names:
``_FLT_`` -> flight, ``_GC_`` -> ground. Vivarium (``VIV``), basal (``BSL``) and
the OSD-253 ground re-run batch (``GCrerun``) are excluded so every cohort
contributes a clean hardware-ground-control contrast. OSD-462 uses the totRNA
library (rRNA-depleted total RNA), matching the rRNA-removed library type of the
other four cohorts.

Existing ``contrast_vectors/cross_osdr_recurrence/`` artifacts are PROGENy
*pathway*-level bootstraps, not per-gene effect+SE, so they cannot substitute for
the per-gene meta inputs; per-cohort effects are computed from VST here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from src.common import id_map_lookup
from src.v11.core_analysis import bh
from src.networks.cross_osdr_projection import signed_stouffer_z

# --------------------------------------------------------------------------- #
# Configuration / on-disk locations
# --------------------------------------------------------------------------- #

ID_MAP_PATH = Path("data/processed/resources/id_map.tsv")
MECHANISM_GENE_SETS = Path("config/mechanism_gene_sets.yaml")

#: Outcome-anchor gene sets re-scored with the pooled statistic (config keys).
MATRIX_SET_KEY = "ecm_organization"          # matrix / ECM-high program
TRANSPORT_SET_KEY = "dct_ncc_wnk_transport"  # DCT / NCC-WNK transport-low program

#: Five on-disk GeneLab mouse-kidney cohorts. OSD-462 -> totRNA library to match
#: the rRNA-depleted library type of the others.
COHORT_VST: dict[str, str] = {
    "OSD-102": "data/external/osdr/OSD-102/GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-163": "data/external/osdr/OSD-163/GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-253": "data/external/osdr/OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-462": "data/external/osdr/OSD-462/GLDS-462_rna_seq_VST_Counts_totRNA_GLbulkRNAseq.csv",
    "OSD-513": "data/external/osdr/OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
}

#: OSD-163 is the BAL-TAL strain; the rest are C57-derived. Recorded for the
#: heterogeneity narrative (a strain difference *adds* robustness if the signal
#: survives it).
COHORT_STRAIN: dict[str, str] = {
    "OSD-102": "C57-6J",
    "OSD-163": "BAL-TAL",
    "OSD-253": "C3H-HeJ + C57-6J (pooled)",
    "OSD-462": "C57 WT (RR-10)",
    "OSD-513": "C57 (RR-23)",
}

MIN_REPLICATES = 2   # per arm, to admit a cohort's effect for a gene
MIN_COHORTS = 2      # genes must appear in >= this many cohorts to pool


# --------------------------------------------------------------------------- #
# Design resolution
# --------------------------------------------------------------------------- #

def resolve_design(columns: Sequence[str]) -> dict[str, str]:
    """Map VST sample columns to ``flight`` / ``ground`` from GeneLab names.

    Flight := ``_FLT_``; ground := hardware ground control ``_GC_``. Vivarium,
    basal and the OSD-253 ``GCrerun`` batch are excluded (returned for neither
    arm) so each cohort contributes one clean hardware-ground contrast.
    """
    design: dict[str, str] = {}
    for col in columns:
        c = str(col).upper()
        if "GCRERUN" in c:
            continue
        if "_FLT_" in c or c.endswith("_FLT"):
            design[col] = "flight"
        elif "_GC_" in c or c.endswith("_GC"):
            design[col] = "ground"
    return design


# --------------------------------------------------------------------------- #
# (1) per-cohort effect + SE  (pure)
# --------------------------------------------------------------------------- #

def per_cohort_effect(
    vst: pd.DataFrame,
    design: Mapping[str, str],
    *,
    min_replicates: int = MIN_REPLICATES,
    moderate: bool = False,
) -> pd.DataFrame:
    """Per-gene flight-vs-ground effect and SE on the VST scale (Welch).

    ``vst`` is genes x samples (index = gene id). ``design`` maps a subset of
    columns to ``flight`` / ``ground``. The effect is ``mean(flight) -
    mean(ground)`` (positive = up in flight; sign-faithful), with Welch SE
    ``sqrt(var_f/n_f + var_g/n_g)`` and a Welch-Satterthwaite df.

    ``moderate=True`` adds an empirical-Bayes variance floor (the 10th-percentile
    sampling variance across genes) to stabilise low-variance genes, a light
    stand-in for limma's variance shrinkage.
    """
    flight_cols = [c for c, g in design.items() if g == "flight" and c in vst.columns]
    ground_cols = [c for c, g in design.items() if g == "ground" and c in vst.columns]
    if len(flight_cols) < min_replicates or len(ground_cols) < min_replicates:
        raise ValueError(
            f"insufficient replicates: flight={len(flight_cols)} ground={len(ground_cols)} "
            f"(need >= {min_replicates} each)"
        )

    f = vst[flight_cols].to_numpy(dtype=float)
    g = vst[ground_cols].to_numpy(dtype=float)
    n_f, n_g = f.shape[1], g.shape[1]

    mean_f = np.nanmean(f, axis=1)
    mean_g = np.nanmean(g, axis=1)
    var_f = np.nanvar(f, axis=1, ddof=1)
    var_g = np.nanvar(g, axis=1, ddof=1)

    effect = mean_f - mean_g
    samp_var = var_f / n_f + var_g / n_g

    if moderate:
        finite = samp_var[np.isfinite(samp_var) & (samp_var > 0)]
        floor = float(np.percentile(finite, 10)) if finite.size else 0.0
        samp_var = samp_var + floor

    se = np.sqrt(samp_var)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, effect / se, np.nan)
        # Welch-Satterthwaite degrees of freedom
        num = samp_var ** 2
        den = (var_f / n_f) ** 2 / max(n_f - 1, 1) + (var_g / n_g) ** 2 / max(n_g - 1, 1)
        df = np.where(den > 0, num / den, np.nan)

    out = pd.DataFrame(
        {
            "effect": effect,
            "se": se,
            "t": t,
            "df": df,
            "n_flight": n_f,
            "n_ground": n_g,
        },
        index=vst.index,
    )
    # genes with no usable variance cannot be weighted -> drop
    out = out[np.isfinite(out["effect"]) & np.isfinite(out["se"]) & (out["se"] > 0)]
    return out


# --------------------------------------------------------------------------- #
# (2) random-effects meta-analysis  (pure)
# --------------------------------------------------------------------------- #

def _dersimonian_laird(y: np.ndarray, v: np.ndarray) -> dict[str, float]:
    """DerSimonian-Laird random-effects pool for one gene.

    ``y`` = per-cohort effects, ``v`` = per-cohort variances (se^2). Returns the
    pooled effect, its SE, tau^2, Cochran's Q, I^2 and a two-sided Wald p.
    """
    k = y.size
    w = 1.0 / v
    sw = w.sum()
    ybar_fe = float((w * y).sum() / sw)
    Q = float((w * (y - ybar_fe) ** 2).sum())
    df = k - 1
    C = sw - (w ** 2).sum() / sw
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0

    wstar = 1.0 / (v + tau2)
    swstar = wstar.sum()
    pooled = float((wstar * y).sum() / swstar)
    se_pooled = float(np.sqrt(1.0 / swstar))
    z = pooled / se_pooled if se_pooled > 0 else np.nan
    p = float(2.0 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    i2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0

    # Stouffer cross-check on per-cohort z = y / sqrt(v)
    stouffer = signed_stouffer_z((y / np.sqrt(v)).tolist())

    return {
        "k_cohorts": k,
        "pooled_effect": pooled,
        "se_pooled": se_pooled,
        "ci_low": pooled - 1.96 * se_pooled,
        "ci_high": pooled + 1.96 * se_pooled,
        "z": z,
        "p_value": p,
        "tau2": tau2,
        "Q": Q,
        "df": df,
        "I2": i2,
        "stouffer_z": float(stouffer["z"]),
        "stouffer_p": float(stouffer["p_two_sided"]),
    }


def meta_random_effects(
    effects: Mapping[str, pd.DataFrame],
    *,
    min_cohorts: int = MIN_COHORTS,
) -> pd.DataFrame:
    """Per-gene DerSimonian-Laird meta across cohorts, with BH-FDR.

    ``effects`` maps cohort -> per-gene table (output of :func:`per_cohort_effect`,
    index = gene id, columns include ``effect`` and ``se``). Genes present in at
    least ``min_cohorts`` cohorts are pooled. Adds a BH-FDR over the pooled
    p-values and a Stouffer's-Z cross-check column.
    """
    eff = pd.concat({c: d["effect"] for c, d in effects.items()}, axis=1)
    sev = pd.concat({c: d["se"] for c, d in effects.items()}, axis=1)

    rows: list[dict[str, float]] = []
    index: list[str] = []
    for gene in eff.index:
        y = eff.loc[gene].to_numpy(dtype=float)
        s = sev.loc[gene].to_numpy(dtype=float)
        ok = np.isfinite(y) & np.isfinite(s) & (s > 0)
        if ok.sum() < min_cohorts:
            continue
        res = _dersimonian_laird(y[ok], s[ok] ** 2)
        rows.append(res)
        index.append(gene)

    out = pd.DataFrame(rows, index=pd.Index(index, name="gene_id"))
    if not out.empty:
        out["fdr"] = bh(out["p_value"].to_numpy())
    return out


# --------------------------------------------------------------------------- #
# Gene-set scoring with the pooled statistic
# --------------------------------------------------------------------------- #

def _load_gene_set(key: str) -> list[str]:
    cfg = yaml.safe_load(MECHANISM_GENE_SETS.read_text())
    return [str(g) for g in cfg[key]["genes"]]


def score_gene_set(
    meta: pd.DataFrame,
    genes: Sequence[str],
    symbol_to_ens: Mapping[str, set[str]],
) -> dict[str, object]:
    """Aggregate the pooled per-gene statistic over a curated gene set.

    Maps set symbols (case-insensitive) to Ensembl ids, intersects with the meta
    table, and returns a precision-weighted set effect, a Stouffer combination of
    the per-gene meta z, median I^2 and median per-gene FDR.
    """
    ens_ids: set[str] = set()
    for sym in genes:
        ens_ids |= set(symbol_to_ens.get(str(sym).lower(), set()))
    present = meta.index.intersection(sorted(ens_ids))
    if len(present) == 0:
        return {
            "n_genes_in_set": len(genes),
            "n_genes_present": 0,
            "set_effect": np.nan,
            "set_effect_precision_wt": np.nan,
            "stouffer_z": np.nan,
            "stouffer_p": np.nan,
            "median_fdr": np.nan,
            "min_fdr": np.nan,
            "median_I2": np.nan,
        }
    sub = meta.loc[present]
    w = 1.0 / (sub["se_pooled"].to_numpy() ** 2)
    eff = sub["pooled_effect"].to_numpy()
    prec_wt = float((w * eff).sum() / w.sum())
    stou = signed_stouffer_z(sub["z"].tolist())
    return {
        "n_genes_in_set": len(genes),
        "n_genes_present": int(len(present)),
        "set_effect": float(np.mean(eff)),
        "set_effect_precision_wt": prec_wt,
        "stouffer_z": float(stou["z"]),
        "stouffer_p": float(stou["p_two_sided"]),
        "median_fdr": float(sub["fdr"].median()),
        "min_fdr": float(sub["fdr"].min()),
        "median_I2": float(sub["I2"].median()),
    }


# --------------------------------------------------------------------------- #
# (3) leave-one-cohort-out  (pure)
# --------------------------------------------------------------------------- #

def leave_one_cohort_out(
    effects: Mapping[str, pd.DataFrame],
    gene_sets: Mapping[str, Sequence[str]],
    symbol_to_ens: Mapping[str, set[str]],
    *,
    min_cohorts: int = MIN_COHORTS,
) -> pd.DataFrame:
    """Re-estimate each gene set's pooled score dropping one cohort at a time.

    Returns one row per (dropped_cohort, gene_set) plus the ``__none__`` full-set
    baseline, so stability of the set-level effect/FDR/I^2 is auditable.
    """
    cohorts = list(effects.keys())
    rows: list[dict[str, object]] = []

    def _emit(drop_label: str, subset: Mapping[str, pd.DataFrame]) -> None:
        if len(subset) < min_cohorts:
            return
        meta = meta_random_effects(subset, min_cohorts=min_cohorts)
        for set_name, genes in gene_sets.items():
            sc = score_gene_set(meta, genes, symbol_to_ens)
            rows.append({"dropped_cohort": drop_label, "gene_set": set_name, **sc})

    _emit("__none__", effects)
    for c in cohorts:
        _emit(c, {k: v for k, v in effects.items() if k != c})

    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #

def _load_cohort_effects(
    ens_to_symbol: Mapping[str, str],
    *,
    moderate: bool = False,
) -> dict[str, pd.DataFrame]:
    effects: dict[str, pd.DataFrame] = {}
    for cohort, path in COHORT_VST.items():
        p = Path(path)
        if not p.exists():
            continue
        vst = pd.read_csv(p, index_col=0)
        design = resolve_design(vst.columns)
        eff = per_cohort_effect(vst, design, moderate=moderate)
        eff["gene_symbol"] = [ens_to_symbol.get(g, "") for g in eff.index]
        eff["strain"] = COHORT_STRAIN.get(cohort, "")
        effects[cohort] = eff
    return effects


def run_recurrence_meta(root: str | Path, *, moderate: bool = False) -> dict[str, object]:
    """End-to-end Repair B: per-cohort effects -> meta -> gene-set re-score -> LOO."""
    root = Path(root)
    out_dir = root / "cross_osdr_recurrence"
    out_dir.mkdir(parents=True, exist_ok=True)

    ens_to_symbol, symbol_to_ens = id_map_lookup(ID_MAP_PATH)
    effects = _load_cohort_effects(ens_to_symbol, moderate=moderate)
    if len(effects) < MIN_COHORTS:
        raise RuntimeError(f"only {len(effects)} cohorts loaded; need >= {MIN_COHORTS}")

    meta = meta_random_effects(effects)
    meta = meta.copy()
    meta["gene_symbol"] = [ens_to_symbol.get(g, "") for g in meta.index]

    matrix_genes = _load_gene_set(MATRIX_SET_KEY)
    transport_genes = _load_gene_set(TRANSPORT_SET_KEY)
    gene_sets = {"matrix_ecm_high": matrix_genes, "dct_ncc_wnk_transport": transport_genes}

    set_scores = {
        name: score_gene_set(meta, genes, symbol_to_ens)
        for name, genes in gene_sets.items()
    }
    loo = leave_one_cohort_out(effects, gene_sets, symbol_to_ens)

    # ---- write artifacts -------------------------------------------------- #
    summary_path = out_dir / "recurrence_meta_summary.tsv"
    loo_path = out_dir / "recurrence_meta_leave_one_out.tsv"
    set_scores_path = out_dir / "recurrence_meta_gene_set_scores.tsv"
    meta_out = meta.reset_index().sort_values("p_value")
    meta_out.to_csv(summary_path, sep="\t", index=False)
    loo.to_csv(loo_path, sep="\t", index=False)
    # index-compatible set-level scores: one row per gene set, selectable by name
    set_scores_df = pd.DataFrame(
        [{"gene_set": name, **sc} for name, sc in set_scores.items()]
    )
    set_scores_df.to_csv(set_scores_path, sep="\t", index=False)

    transport = set_scores["dct_ncc_wnk_transport"]
    matrix = set_scores["matrix_ecm_high"]

    # leave-one-out stability for the transport set
    t_loo = loo[(loo["gene_set"] == "dct_ncc_wnk_transport") & (loo["dropped_cohort"] != "__none__")]
    transport_loo_effects = t_loo["set_effect"].to_numpy(dtype=float)
    transport_loo_sign_stable = bool(
        np.all(np.sign(transport_loo_effects) == np.sign(transport["set_effect"]))
    ) if transport_loo_effects.size else False

    verdict = {
        "analysis": "repair_b_cross_cohort_recurrence_meta",
        "n_cohorts": len(effects),
        "cohorts": list(effects.keys()),
        "strains": {c: COHORT_STRAIN.get(c, "") for c in effects},
        "variance_moderated": bool(moderate),
        "n_genes_meta": int(len(meta)),
        "n_genes_fdr05": int((meta["fdr"] < 0.05).sum()) if not meta.empty else 0,
        "transport_set": transport,
        "matrix_set": matrix,
        "transport_loo_sign_stable": transport_loo_sign_stable,
        # ---- headline-index provenance keys ----
        "recurrence_meta_transport_effect": transport["set_effect"],
        "recurrence_meta_transport_fdr": transport["median_fdr"],
        "recurrence_meta_transport_i2": transport["median_I2"],
        "outputs": {
            "summary": str(summary_path),
            "leave_one_out": str(loo_path),
            "gene_set_scores": str(set_scores_path),
        },
    }
    (out_dir / "recurrence_meta_verdict.json").write_text(json.dumps(verdict, indent=2))
    return verdict


def main() -> None:
    ap = argparse.ArgumentParser(description="Repair B -- cross-cohort recurrence meta-analysis")
    ap.add_argument("--run-root", required=True, help="results run directory to write under")
    ap.add_argument("--moderate", action="store_true", help="apply empirical-Bayes variance floor")
    args = ap.parse_args()
    verdict = run_recurrence_meta(args.run_root, moderate=args.moderate)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
