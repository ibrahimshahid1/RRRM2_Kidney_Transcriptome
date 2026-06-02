#!/usr/bin/env python3
"""Reach F -- explicit aldosterone / mineralocorticoid-axis test.

Motivation (manuscript glue). NCC is aldosterone-regulated and spaceflight's
cephalad fluid shift perturbs the renin-angiotensin-aldosterone system, yet the
manuscript never tests the aldosterone axis -- even though the regulator run
flagged a steroid-receptor ("Androgen") PROGENy pathway as up (MR and androgen
receptors share hormone-response elements). This module scores a curated
mineralocorticoid panel across the same five cohorts used by the recurrence
analysis, cross-checks it against the already-computed PROGENy activities, and
asks whether the aldosterone-axis direction is *coherent* with the predicted
distal-nephron (NCC/WNK-SPAK) suppression -- reported side by side with the
transport-anchor score.

Design notes.
  * The MR panel (``data/external/gene_sets/mineralocorticoid_axis_core.tsv``)
    is directional: ``up`` genes are aldosterone-induced (Sgk1, Tsc22d3, Per1,
    Klf15, Fxyd2, Scnn1a/b/g, Slc12a3); ``context`` genes (Hsd11b2, Nr3c2,
    Wnk1/4, Stk39, Oxsr1) gate or organize the axis and are reported but
    excluded from the directional score.
  * Scoring is decoupleR-free on purpose: the panel is small and curated, the
    sandbox has no decoupler, and a direction-corrected mean with a gene-label
    permutation null is fully reproducible. Cross-cohort pooling reuses the
    project's ``signed_stouffer_z`` and ``recurrence_class``.
  * Sign convention: a *positive* axis effect means the aldosterone-responsive
    program moves up; the manuscript's distal-nephron suppression prediction is
    a *negative* axis effect that tracks NCC/WNK-SPAK transport suppression.

Discipline (failure mode is informative). Bulk RNA cannot separate systemic
hormone exposure from intrinsic DCT response; a flat or incoherent result
argues against an endocrine-driven story and narrows the mechanism. Nothing
here is a mechanistic claim.

Inputs (on disk): the five cohort gene-stat tables under the regulator run's
``rna_effects/`` (gene, stat), the curated MR panel, and the existing
``rna_progeny_pathway_activity.tsv``.
Output: ``regulator_activity/aldosterone_axis_summary.tsv`` (+ verdict JSON).
Provenance key for the headline index: ``aldo_axis_meta_effect``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.networks.cross_osdr_projection import signed_stouffer_z
from src.multiomics.regulator_activity import recurrence_class
from src.v11.core_analysis import ANCHOR_GENES, RUN_ROOT, bh, read_tsv

MR_PANEL_PATH = Path("data/external/gene_sets/mineralocorticoid_axis_core.tsv")
REGULATOR_RUN = Path("data/results/run_20260522_regulator_activity")
# NCC/WNK-SPAK distal-transport anchor (same five as core_analysis.ANCHOR_GENES).
NCC_WNK_TRANSPORT = sorted(ANCHOR_GENES)
PROGENY_STEROID_SOURCE = "Androgen"  # the puzzling steroid-receptor hit to reinterpret
N_PERM = 2000
PERM_SEED = 20260601
RECURRENCE_THRESHOLD = 0.5  # on the z-like gene-stat scale; sign-consistency reported separately


# --------------------------------------------------------------------------- #
# Pure scoring (unit-tested; no I/O).                                          #
# --------------------------------------------------------------------------- #
def panel_signs(panel: pd.DataFrame) -> pd.Series:
    """Map each panel gene to a directional sign: ``up`` -> +1, else 0 (context)."""
    direction = panel.set_index("gene_symbol")["expected_aldosterone_direction"].astype(str).str.lower()
    return direction.map({"up": 1.0, "down": -1.0}).fillna(0.0)


def score_axis_in_cohort(gene_stats: pd.Series, signs: pd.Series) -> dict:
    """Direction-corrected mean of a directional panel within one cohort.

    ``gene_stats``: Series gene -> flight-effect statistic (z-like).
    ``signs``: Series gene -> {+1, -1, 0}; only nonzero-sign genes contribute.
    Returns the signed axis effect (positive = aldosterone program up), the
    unsigned mean for transparency, a within-panel Wilcoxon p vs 0, and the
    list of contributing genes.
    """
    directional = signs[signs != 0.0]
    present = [g for g in directional.index if g in gene_stats.index and np.isfinite(gene_stats.get(g, np.nan))]
    if not present:
        return {"effect": np.nan, "effect_unsigned": np.nan, "n_present": 0,
                "wilcoxon_p": np.nan, "present_genes": []}
    raw = gene_stats.loc[present].astype(float)
    corrected = raw.values * directional.loc[present].values
    effect = float(np.mean(corrected))
    wilcoxon_p = np.nan
    nz = corrected[corrected != 0.0]
    if nz.size >= 3 and not np.allclose(nz, nz[0]):
        try:
            wilcoxon_p = float(stats.wilcoxon(nz, alternative="two-sided").pvalue)
        except Exception:
            wilcoxon_p = np.nan
    return {
        "effect": effect,
        "effect_unsigned": float(np.mean(raw.values)),
        "n_present": int(len(present)),
        "wilcoxon_p": wilcoxon_p,
        "present_genes": present,
    }


def competitive_permutation(
    gene_stats: pd.Series, present_genes: list[str], present_signs: np.ndarray,
    n_perm: int = N_PERM, seed: int = PERM_SEED,
) -> dict:
    """Gene-label permutation null for the direction-corrected panel mean.

    Draws ``n_perm`` random gene sets of the same size from the cohort's gene
    universe, assigns them the observed panel sign vector, and recomputes the
    direction-corrected mean. Returns a two-sided competitive p and a one-sided
    *suppression* p (observed <= null), the relevant tail for the manuscript's
    distal-nephron suppression prediction.
    """
    universe = gene_stats.dropna().astype(float)
    k = len(present_genes)
    if k == 0 or len(universe) <= k:
        return {"p_two_sided": np.nan, "p_suppression": np.nan, "null_mean": np.nan, "null_sd": np.nan, "obs": np.nan}
    obs = float(np.mean(universe.loc[present_genes].values * present_signs))
    rng = np.random.default_rng(seed)
    vals = universe.values
    signs = np.asarray(present_signs, dtype=float)
    null = np.empty(n_perm, dtype=float)
    for i in range(n_perm):
        idx = rng.choice(len(vals), size=k, replace=False)
        null[i] = np.mean(vals[idx] * signs)
    null_mean = float(np.mean(null))
    p_two = float((np.sum(np.abs(null - null_mean) >= abs(obs - null_mean)) + 1) / (n_perm + 1))
    p_supp = float((np.sum(null <= obs) + 1) / (n_perm + 1))
    return {"p_two_sided": p_two, "p_suppression": p_supp, "null_mean": null_mean,
            "null_sd": float(np.std(null)), "obs": obs}


def meta_axis(per_cohort: dict[str, float], threshold: float = RECURRENCE_THRESHOLD) -> dict:
    """Pool per-cohort axis effects into a cross-cohort verdict.

    Reuses ``signed_stouffer_z`` for a combined Z/p and ``recurrence_class`` for
    the up/down/mixed label, and adds a directional sign-test (how many cohorts
    are negative, i.e. consistent with predicted suppression).
    """
    items = {k: float(v) for k, v in per_cohort.items() if np.isfinite(v)}
    n = len(items)
    if n == 0:
        return {"meta_effect": np.nan, "meta_z": np.nan, "p_two_sided": np.nan, "n_cohorts": 0,
                "n_negative": 0, "n_positive": 0, "sign_consistent": False, "binomial_p": np.nan,
                "recurrence_class": "insufficient_cohorts"}
    effects = np.array(list(items.values()), dtype=float)
    stouffer = signed_stouffer_z(effects)
    n_neg = int(np.sum(effects < 0))
    n_pos = int(np.sum(effects > 0))
    k = max(n_neg, n_pos)
    binom_p = float(stats.binomtest(k, n, 0.5).pvalue) if n else np.nan
    return {
        "meta_effect": float(np.mean(effects)),
        "meta_z": float(stouffer["z"]),
        "p_two_sided": float(stouffer["p_two_sided"]),
        "n_cohorts": n,
        "n_negative": n_neg,
        "n_positive": n_pos,
        "sign_consistent": bool(n_neg == n or n_pos == n),
        "binomial_p": binom_p,
        "recurrence_class": recurrence_class(items, threshold=threshold),
    }


# --------------------------------------------------------------------------- #
# Orchestration (I/O).                                                         #
# --------------------------------------------------------------------------- #
def _load_cohorts(regulator_run: Path) -> dict[str, pd.Series]:
    spec_path = regulator_run / "rna_effects" / "rna_effects_spec.json"
    spec = json.loads(spec_path.read_text())
    cohorts: dict[str, pd.Series] = {}
    for name, rel in spec.items():
        p = Path(rel)
        if not p.exists():
            p = regulator_run / "rna_effects" / Path(rel).name
        if not p.exists():
            continue
        df = read_tsv(p)
        gcol = "gene" if "gene" in df.columns else df.columns[0]
        scol = "stat" if "stat" in df.columns else df.columns[1]
        cohorts[name] = pd.to_numeric(df.set_index(gcol)[scol], errors="coerce").dropna()
    return cohorts


def _load_progeny_steroid(regulator_run: Path) -> dict[str, float]:
    path = regulator_run / "rna_progeny_pathway_activity.tsv"
    if not path.exists():
        return {}
    prog = read_tsv(path)
    sub = prog[prog["source"] == PROGENY_STEROID_SOURCE]
    return dict(zip(sub["contrast"].astype(str), pd.to_numeric(sub["activity_score"], errors="coerce")))


def _row(**kw) -> dict:
    base = {"panel": "mineralocorticoid_core", "level": None, "cohort": None,
            "statistic": None, "value": np.nan, "p_value": np.nan, "n_genes": np.nan,
            "detail": None, "provenance_key": None}
    base.update(kw)
    return base


def run_aldosterone_axis(root: Path, regulator_run: Path = REGULATOR_RUN) -> pd.DataFrame:
    root = Path(root)
    panel = read_tsv(MR_PANEL_PATH)
    signs = panel_signs(panel)
    cohorts = _load_cohorts(regulator_run)
    if not cohorts:
        raise FileNotFoundError(f"No cohort gene-stat tables under {regulator_run}/rna_effects/.")
    androgen = _load_progeny_steroid(regulator_run)

    rows: list[dict] = []
    per_cohort_axis: dict[str, float] = {}
    per_cohort_transport: dict[str, float] = {}
    for name, gs in cohorts.items():
        sc = score_axis_in_cohort(gs, signs)
        present = sc["present_genes"]
        present_signs = signs.loc[present].values if present else np.array([])
        perm = competitive_permutation(gs, present, present_signs)
        # NCC/WNK-SPAK transport anchor score (unsigned mean; negative = suppressed).
        tgenes = [g for g in NCC_WNK_TRANSPORT if g in gs.index and np.isfinite(gs.get(g, np.nan))]
        transport = float(np.mean(gs.loc[tgenes].values)) if tgenes else np.nan
        per_cohort_axis[name] = sc["effect"]
        per_cohort_transport[name] = transport

        rows.append(_row(level="per_cohort", cohort=name, statistic="mr_axis_effect",
                         value=sc["effect"], p_value=perm["p_two_sided"], n_genes=sc["n_present"],
                         detail=f"perm_suppression_p={perm['p_suppression']:.4g}; wilcoxon_p={sc['wilcoxon_p']}"))
        rows.append(_row(level="per_cohort", cohort=name, statistic="mr_axis_effect_unsigned",
                         value=sc["effect_unsigned"], n_genes=sc["n_present"],
                         detail="mean stat over core up-genes (no direction correction)"))
        rows.append(_row(level="per_cohort", cohort=name, statistic="ncc_wnk_transport_effect",
                         value=transport, n_genes=len(tgenes),
                         detail=f"anchor genes: {','.join(tgenes)}"))
        if name in androgen:
            rows.append(_row(level="per_cohort", cohort=name, statistic="androgen_progeny_activity",
                             value=float(androgen[name]), detail="from rna_progeny_pathway_activity.tsv"))

    meta = meta_axis(per_cohort_axis)

    # Cross-checks across cohorts (low n; report r and sign).
    common = [c for c in per_cohort_axis if np.isfinite(per_cohort_axis[c])]
    def _corr(a: dict, b: dict) -> tuple[float, float, int]:
        keys = [k for k in common if k in b and np.isfinite(a.get(k, np.nan)) and np.isfinite(b.get(k, np.nan))]
        if len(keys) < 3:
            return np.nan, np.nan, len(keys)
        r, p = stats.spearmanr([a[k] for k in keys], [b[k] for k in keys])
        return float(r), float(p), len(keys)

    r_transport, p_transport, n_tr = _corr(per_cohort_axis, per_cohort_transport)
    r_andro, p_andro, n_an = _corr(per_cohort_axis, androgen)
    both_suppressed = int(sum(
        1 for c in common if per_cohort_axis[c] < 0 and np.isfinite(per_cohort_transport.get(c, np.nan))
        and per_cohort_transport[c] < 0
    ))

    rows.append(_row(level="meta", cohort="ALL", statistic="aldo_axis_meta_effect",
                     value=meta["meta_effect"], p_value=meta["p_two_sided"], n_genes=meta["n_cohorts"],
                     detail=(f"meta_z={meta['meta_z']:.4g}; recurrence={meta['recurrence_class']}; "
                             f"n_negative={meta['n_negative']}/{meta['n_cohorts']}; "
                             f"sign_consistent={meta['sign_consistent']}; binomial_p={meta['binomial_p']:.4g}"),
                     provenance_key="aldo_axis_meta_effect"))
    rows.append(_row(level="meta", cohort="ALL", statistic="mr_axis_vs_ncc_wnk_transport_spearman",
                     value=r_transport, p_value=p_transport, n_genes=n_tr,
                     detail=f"both_suppressed_cohorts={both_suppressed}/{len(common)}"))
    rows.append(_row(level="meta", cohort="ALL", statistic="mr_axis_vs_androgen_progeny_spearman",
                     value=r_andro, p_value=p_andro, n_genes=n_an,
                     detail="reinterprets the PROGENy Androgen hit in a kidney MR frame"))

    summary = pd.DataFrame(rows)
    # BH across the per-cohort axis permutation tests.
    axis_mask = (summary["level"] == "per_cohort") & (summary["statistic"] == "mr_axis_effect")
    summary["q_value"] = np.nan
    if axis_mask.any():
        summary.loc[axis_mask, "q_value"] = bh(summary.loc[axis_mask, "p_value"])

    out_dir = root / "regulator_activity"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "aldosterone_axis_summary.tsv"
    summary.to_csv(out_path, sep="\t", index=False)
    _write_verdict(root, meta, per_cohort_axis, per_cohort_transport,
                   (r_transport, both_suppressed, len(common)), (r_andro, androgen))
    print(f"Reach F aldosterone axis written: {out_path} ({len(summary)} rows)")
    return summary


def _write_verdict(root, meta, axis, transport, transport_cc, andro_cc):
    r_transport, both_suppressed, n_common = transport_cc
    r_andro, androgen = andro_cc
    suppressed = bool(np.isfinite(meta["meta_effect"]) and meta["meta_effect"] < 0)
    coherent = bool(suppressed and meta["sign_consistent"]
                    and (not np.isfinite(r_transport) or r_transport > 0)
                    and both_suppressed >= max(1, n_common - 1))
    andro_mean = float(np.nanmean(list(androgen.values()))) if androgen else np.nan

    if coherent:
        interp = ("The mineralocorticoid axis moves down coherently with NCC/WNK-SPAK transport "
                  "suppression across cohorts -- consistent with a distal-nephron aldosterone-axis "
                  "down-shift. Report as an exploratory endocrine-axis corroboration, and pre-register "
                  "the same direction for the human (Reach D) test.")
    elif suppressed:
        interp = ("The mineralocorticoid axis trends down but is not fully consistent across cohorts or "
                  "with the transport anchor; report as a weak/partial aldosterone-axis suppression.")
    elif np.isfinite(andro_mean) and andro_mean > 0:
        interp = ("The curated MR panel does NOT move down; the PROGENy Androgen-up hit is therefore not "
                  "an aldosterone-suppression signal in this panel -- likely a shared hormone-response-element "
                  "artifact or systemic/vascular steroid signal, not intrinsic DCT MR suppression. This "
                  "narrows the mechanism and argues against an endocrine-driven distal story.")
    else:
        interp = ("No coherent mineralocorticoid-axis signal; bulk RNA cannot resolve an endocrine driver here.")

    verdict = {
        "analysis": "Reach F -- aldosterone / mineralocorticoid-axis test",
        "panel": "mineralocorticoid_core",
        "claim_discipline": "Curated panel over bulk RNA; cannot separate systemic hormone exposure from intrinsic DCT response.",
        "meta_effect": None if not np.isfinite(meta["meta_effect"]) else float(meta["meta_effect"]),
        "meta_z": None if not np.isfinite(meta["meta_z"]) else float(meta["meta_z"]),
        "meta_p_two_sided": None if not np.isfinite(meta["p_two_sided"]) else float(meta["p_two_sided"]),
        "recurrence_class": meta["recurrence_class"],
        "n_cohorts": meta["n_cohorts"],
        "n_negative": meta["n_negative"],
        "sign_consistent": meta["sign_consistent"],
        "binomial_sign_p": None if not np.isfinite(meta["binomial_p"]) else float(meta["binomial_p"]),
        "mr_axis_vs_transport_spearman": None if not np.isfinite(r_transport) else float(r_transport),
        "both_suppressed_cohorts": f"{both_suppressed}/{n_common}",
        "androgen_progeny_mean": None if not np.isfinite(andro_mean) else andro_mean,
        "coherent_with_distal_suppression": coherent,
        "interpretation": interp,
    }
    (root / "regulator_activity" / "aldosterone_axis_verdict.json").write_text(json.dumps(verdict, indent=2))


def main():
    parser = argparse.ArgumentParser(description="Reach F: aldosterone/mineralocorticoid-axis test.")
    parser.add_argument("--run-root", default=str(RUN_ROOT))
    parser.add_argument("--regulator-run", default=str(REGULATOR_RUN),
                        help="Run dir holding rna_effects/ and rna_progeny_pathway_activity.tsv.")
    args = parser.parse_args()
    run_aldosterone_axis(Path(args.run_root), Path(args.regulator_run))


if __name__ == "__main__":
    main()
