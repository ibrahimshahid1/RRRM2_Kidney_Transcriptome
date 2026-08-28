#!/usr/bin/env python3
"""Test barrier-core expression beyond a disjoint podocyte-marker proxy."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions  # noqa: E402
from src.clinical_axes.statistics import random_effects_reml_mkh, score_signed_axis  # noqa: E402


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
DEFAULT_RESULTS = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"


def regression_rows(mission, data, outcome, proxy, adjusted):
    design = pd.DataFrame(
        {
            "flight": (data.metadata["condition"] == "FLT").astype(float),
        },
        index=data.metadata.index,
    )
    block = pd.get_dummies(data.metadata["block"], drop_first=True, dtype=float)
    design = design.join(block)
    if adjusted:
        design["disjoint_podocyte_proxy"] = proxy
    design = sm.add_constant(design, has_constant="add")
    ordinary = sm.OLS(outcome, design).fit()
    robust = ordinary.get_robustcov_results(cov_type="HC3")
    index = list(design.columns).index("flight")
    rows = []
    for variance_type, fit in (("model_based", ordinary), ("HC3", robust)):
        estimate = float(fit.params[index] if variance_type == "HC3" else fit.params["flight"])
        se = float(fit.bse[index] if variance_type == "HC3" else fit.bse["flight"])
        rows.append(
            {
                "mission": mission,
                "model": "adjusted_disjoint_podocyte_proxy" if adjusted else "unadjusted_blocked",
                "variance_type": variance_type,
                "estimate": estimate,
                "standard_error": se,
                "variance": se**2,
                "p": float(fit.pvalues[index] if variance_type == "HC3" else fit.pvalues["flight"]),
                "n": len(outcome),
                "outcome_proxy_pearson_r": float(outcome.corr(proxy)),
            }
        )
    return rows


def run(args):
    config = yaml.safe_load(args.config.read_text())
    missions = load_primary_missions(config, REPO)
    threshold = float(config["eligibility"]["cpm_threshold"])
    barrier = list(
        config["primary_family"]["glomerular_barrier_identity_loss"]["subdomains"]
        ["barrier_core"]["genes"]
    )
    tiers = pd.read_csv(args.tiers, sep="\t")
    podocyte = list(
        dict.fromkeys(
            tiers.loc[
                (tiers["gene_set"] == "podocyte__high_specificity")
                & tiers["final_for_testing"].astype(bool),
                "gene_symbol",
            ].astype(str)
        )
    )
    # The proxy is disjoint by construction, so adjustment cannot mechanically
    # subtract the same genes used in the outcome.
    proxy_requested = [gene for gene in podocyte if gene not in set(barrier)]
    rows = []
    coverage = []
    score_rows = []
    for mission, data in missions.items():
        eligible = cpm_eligible_genes(data, threshold)
        barrier_used = [g for g in barrier if g in eligible and g in data.expression.index]
        proxy_used = [
            g for g in proxy_requested if g in eligible and g in data.expression.index
        ]
        if len(barrier_used) != len(barrier) or len(proxy_used) < 8:
            raise RuntimeError(
                f"{mission}: barrier={len(barrier_used)}, proxy={len(proxy_used)}"
            )
        outcome = score_signed_axis(
            data.expression, {gene: 1 for gene in barrier_used}
        ).scores
        proxy = score_signed_axis(
            data.expression, {gene: 1 for gene in proxy_used}
        ).scores
        coverage.append(
            {
                "mission": mission,
                "n_barrier_genes": len(barrier_used),
                "n_disjoint_proxy_genes": len(proxy_used),
                "barrier_genes": "|".join(barrier_used),
                "proxy_genes": "|".join(proxy_used),
            }
        )
        rows.extend(regression_rows(mission, data, outcome, proxy, adjusted=False))
        rows.extend(regression_rows(mission, data, outcome, proxy, adjusted=True))
        for animal in data.metadata.index:
            score_rows.append(
                {
                    "mission": mission,
                    "animal": animal,
                    "condition": data.metadata.loc[animal, "condition"],
                    "block": data.metadata.loc[animal, "block"],
                    "barrier_expression_score": outcome.loc[animal],
                    "disjoint_podocyte_proxy": proxy.loc[animal],
                }
            )
    effects = pd.DataFrame(rows)
    meta_rows = []
    for (model, variance_type), sub in effects.groupby(
        ["model", "variance_type"], sort=False
    ):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        meta_rows.append(
            {
                "model": model,
                "variance_type": variance_type,
                "estimate": fit.estimate,
                "standard_error_mkh": fit.standard_error,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh": fit.p,
                "tau2": fit.tau2,
                "i_squared": fit.i_squared,
                "maximum_weight": float(fit.weights.max()),
            }
        )
    meta = pd.DataFrame(meta_rows)
    args.results.mkdir(parents=True, exist_ok=True)
    effects.to_csv(
        args.results / "barrier_proxy_adjustment_mission_effects.tsv",
        sep="\t",
        index=False,
    )
    meta.to_csv(
        args.results / "barrier_proxy_adjustment_meta.tsv", sep="\t", index=False
    )
    pd.DataFrame(coverage).to_csv(
        args.results / "barrier_proxy_adjustment_coverage.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(score_rows).to_csv(
        args.results / "barrier_proxy_adjustment_scores.tsv",
        sep="\t",
        index=False,
    )
    print(meta.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

