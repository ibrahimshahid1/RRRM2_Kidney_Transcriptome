#!/usr/bin/env python3
"""Run influence and technical sensitivities for the frozen renal-axis study."""

from __future__ import annotations

import argparse
from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.analysis import (  # noqa: E402
    combined_score_design,
    mission_effect_from_score,
)
from src.clinical_axes.data import (  # noqa: E402
    cpm_eligible_genes,
    load_osd253_control_sensitivity,
    load_osd462_preparation,
    load_primary_missions,
)
from src.clinical_axes.statistics import random_effects_reml_mkh  # noqa: E402


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_RESULTS = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"


def _meta_for_family(missions, family, threshold, method="mean"):
    scores, design, coverage, _ = combined_score_design(
        missions, family, cpm_threshold=threshold, method=method
    )
    rows = []
    for mission, data in missions.items():
        idx = design.index[design["mission"] == mission]
        for axis in scores.columns:
            local = pd.Series(
                scores.loc[idx, axis].to_numpy(), index=data.metadata.index
            )
            summary, _ = mission_effect_from_score(local, data.metadata)
            rows.append({"axis": axis, "mission": mission, **summary})
    mission_effects = pd.DataFrame(rows)
    meta_rows = []
    for axis, sub in mission_effects.groupby("axis", sort=False):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        meta_rows.append(
            {
                "axis": axis,
                "estimate": fit.estimate,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh": fit.p,
                "tau2": fit.tau2,
                "i_squared": fit.i_squared,
                "maximum_weight": float(fit.weights.max()),
                "k_missions": fit.k,
            }
        )
    return pd.DataFrame(meta_rows), mission_effects, coverage


def _all_genes(axis_spec):
    return [
        gene
        for subdomain in axis_spec["subdomains"].values()
        for gene in subdomain["genes"]
    ]


def _omit_gene(axis_spec, omitted):
    modified = deepcopy(axis_spec)
    for subdomain in modified["subdomains"].values():
        if omitted in subdomain["genes"]:
            del subdomain["genes"][omitted]
            subdomain["minimum_present"] = max(
                1, int(subdomain["minimum_present"]) - 1
            )
    return modified


def leave_one_gene(missions, family, threshold, full_meta):
    rows = []
    full = full_meta.set_index("axis")["estimate"]
    for axis, axis_spec in family.items():
        for gene in _all_genes(axis_spec):
            if "gene_level_requirement" in axis_spec:
                rows.append(
                    {
                        "axis": axis,
                        "omitted_gene": gene,
                        "evaluable": False,
                        "reason": "primary endpoint requires both genes",
                    }
                )
                continue
            reduced = {axis: _omit_gene(axis_spec, gene)}
            meta, _, _ = _meta_for_family(missions, reduced, threshold)
            estimate = float(meta.loc[0, "estimate"])
            baseline = float(full[axis])
            rows.append(
                {
                    "axis": axis,
                    "omitted_gene": gene,
                    "evaluable": True,
                    "estimate": estimate,
                    "ci_low_mkh": meta.loc[0, "ci_low_mkh"],
                    "ci_high_mkh": meta.loc[0, "ci_high_mkh"],
                    "direction_retained": bool(np.sign(estimate) == np.sign(baseline)),
                    "magnitude_fraction": abs(estimate) / abs(baseline)
                    if abs(baseline) > 1e-12
                    else np.nan,
                }
            )
    return pd.DataFrame(rows)


def common_intersection_family(missions, family, threshold):
    eligible = {
        mission: cpm_eligible_genes(data, threshold)
        for mission, data in missions.items()
    }
    common = set.intersection(*eligible.values())
    modified = deepcopy(family)
    for axis_spec in modified.values():
        for subdomain in axis_spec["subdomains"].values():
            subdomain["genes"] = {
                gene: sign
                for gene, sign in subdomain["genes"].items()
                if gene in common
            }
            subdomain["minimum_present"] = len(subdomain["genes"])
            if not subdomain["genes"]:
                raise ValueError("Common intersection emptied a subdomain")
    return modified


def expanded_barrier_family(family):
    spec = deepcopy(family["glomerular_barrier_identity_loss"])
    additions = spec.pop("sensitivity_additions")["barrier_expanded"]
    core = spec["subdomains"]["barrier_core"]
    core["genes"].update(additions)
    core["minimum_present"] = 6
    return {"glomerular_barrier_identity_loss_expanded": spec}


def subdomain_family(family):
    out = {}
    for axis, spec in family.items():
        if len(spec["subdomains"]) < 2:
            continue
        for name, subdomain in spec["subdomains"].items():
            out[f"{axis}__{name}"] = {
                "role": "subdomain_sensitivity",
                "subdomains": {name: deepcopy(subdomain)},
            }
    return out


def gene_contributions(gene_effect_path: Path):
    gene_effects = pd.read_csv(gene_effect_path, sep="\t")
    rows = []
    for (axis, gene), sub in gene_effects.groupby(["axis", "gene"], sort=False):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        rows.append(
            {
                "axis": axis,
                "gene": gene,
                "pooled_signed_g": fit.estimate,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh_descriptive": fit.p,
            }
        )
    table = pd.DataFrame(rows)
    table["absolute_contribution_share"] = table.groupby("axis")[
        "pooled_signed_g"
    ].transform(lambda values: values.abs() / values.abs().sum())
    return table


def secondary_qc_covariate_sensitivity(missions, family, threshold):
    """Residualize the flagged OSD-163 mapping metric, then recompute Hedges g.

    This is deliberately a sensitivity rather than the primary estimand.  The
    residualization model is label-blind (score ~ QC metric), and the resulting
    residuals enter the same standardized flight-control effect calculation as
    the other mission scores.
    """
    scores, design, _, _ = combined_score_design(
        missions, family, cpm_threshold=threshold
    )
    _, mission_effects, _ = _meta_for_family(missions, family, threshold)
    target = missions["OSD-163"]
    metric = "uniquely_mapped_percent"
    qc = target.qc[metric].astype(float)
    if qc.isna().any() or qc.std(ddof=1) <= 0:
        raise RuntimeError(f"OSD-163 {metric} is unavailable for sensitivity")
    covariate = (qc - qc.mean()) / qc.std(ddof=1)
    idx = design.index[design["mission"] == "OSD-163"]
    rows = []
    for axis in scores.columns:
        outcome = pd.Series(
            scores.loc[idx, axis].to_numpy(), index=target.metadata.index
        )
        matrix = np.column_stack([np.ones(len(covariate)), covariate.to_numpy()])
        coefficients = np.linalg.lstsq(matrix, outcome.to_numpy(), rcond=None)[0]
        residual = pd.Series(
            outcome.to_numpy() - matrix @ coefficients,
            index=outcome.index,
        )
        adjusted, _ = mission_effect_from_score(residual, target.metadata)
        raw_row = mission_effects[
            (mission_effects["axis"] == axis)
            & (mission_effects["mission"] == "OSD-163")
        ].iloc[0]
        replacement = mission_effects[mission_effects["axis"] == axis].copy()
        mask = replacement["mission"] == "OSD-163"
        replacement.loc[mask, "estimate"] = adjusted["estimate"]
        replacement.loc[mask, "variance"] = adjusted["variance"]
        fit = random_effects_reml_mkh(
            replacement["estimate"], replacement["variance"]
        )
        rows.append(
            {
                "sensitivity": "OSD163_residualized_uniquely_mapped_percent",
                "axis": axis,
                "qc_metric": metric,
                "outcome_qc_pearson_r": float(outcome.corr(qc)),
                "osd163_raw_estimate": raw_row["estimate"],
                "osd163_adjusted_estimate": adjusted["estimate"],
                "meta_estimate": fit.estimate,
                "meta_ci_low_mkh": fit.ci_low,
                "meta_ci_high_mkh": fit.ci_high,
                "meta_p_mkh": fit.p,
                "meta_tau2": fit.tau2,
                "meta_i_squared": fit.i_squared,
            }
        )
    return pd.DataFrame(rows)


def run(args):
    config = yaml.safe_load(args.config.read_text())
    family = config["primary_family"]
    threshold = float(config["eligibility"]["cpm_threshold"])
    missions = load_primary_missions(config, REPO)
    outdir = args.results
    outdir.mkdir(parents=True, exist_ok=True)

    full_meta, _, _ = _meta_for_family(missions, family, threshold)
    loo = leave_one_gene(missions, family, threshold, full_meta)
    loo.to_csv(outdir / "leave_one_gene.tsv", sep="\t", index=False)

    common_family = common_intersection_family(missions, family, threshold)
    common_meta, common_missions, common_coverage = _meta_for_family(
        missions, common_family, threshold
    )
    common_meta.insert(0, "sensitivity", "common_gene_intersection")
    common_meta.to_csv(
        outdir / "common_intersection_meta_results.tsv", sep="\t", index=False
    )
    common_missions.to_csv(
        outdir / "common_intersection_mission_effects.tsv", sep="\t", index=False
    )

    expanded_meta, expanded_missions, expanded_coverage = _meta_for_family(
        missions, expanded_barrier_family(family), threshold
    )
    expanded_meta.insert(0, "sensitivity", "barrier_plus_Podxl_Cd2ap")
    expanded_meta.to_csv(
        outdir / "expanded_barrier_meta_results.tsv", sep="\t", index=False
    )
    expanded_missions.to_csv(
        outdir / "expanded_barrier_mission_effects.tsv", sep="\t", index=False
    )

    domains = subdomain_family(family)
    sub_meta, sub_missions, sub_coverage = _meta_for_family(
        missions, domains, threshold
    )
    sub_meta.to_csv(outdir / "subdomain_meta_results.tsv", sep="\t", index=False)
    sub_missions.to_csv(
        outdir / "subdomain_mission_effects.tsv", sep="\t", index=False
    )

    technical_rows = []
    for preparation in ("mRNA", "UPX"):
        alternate = dict(missions)
        alternate["OSD-462"] = load_osd462_preparation(
            config, REPO, preparation
        )
        meta, _, _ = _meta_for_family(alternate, family, threshold)
        meta.insert(0, "sensitivity", f"OSD462_{preparation}")
        technical_rows.append(meta)
    alternate = dict(missions)
    alternate["OSD-253"] = load_osd253_control_sensitivity(config, REPO)
    meta, _, _ = _meta_for_family(alternate, family, threshold)
    meta.insert(0, "sensitivity", "OSD253_white_light_rerun_control")
    technical_rows.append(meta)
    pd.concat(technical_rows, ignore_index=True).to_csv(
        outdir / "technical_design_sensitivities.tsv", sep="\t", index=False
    )

    secondary_qc_covariate_sensitivity(
        missions, family, threshold
    ).to_csv(
        outdir / "technical_qc_covariate_sensitivity.tsv",
        sep="\t",
        index=False,
    )

    contributions = gene_contributions(
        outdir / "descriptive_signed_gene_effects.tsv"
    )
    contributions.to_csv(
        outdir / "gene_contribution_summary.tsv", sep="\t", index=False
    )

    coverage_tables = []
    for label, table in (
        ("common_intersection", common_coverage),
        ("expanded_barrier", expanded_coverage),
        ("subdomains", sub_coverage),
    ):
        coverage_tables.append(table.assign(sensitivity=label))
    pd.concat(coverage_tables, ignore_index=True).to_csv(
        outdir / "sensitivity_gene_coverage.tsv", sep="\t", index=False
    )

    print("Leave-one-gene barrier:")
    print(
        loo[loo["axis"] == "glomerular_barrier_identity_loss"].to_string(
            index=False
        )
    )
    print("\nTechnical/design sensitivities:")
    print(
        pd.concat(technical_rows, ignore_index=True)[
            ["sensitivity", "axis", "estimate", "ci_low_mkh", "ci_high_mkh"]
        ].to_string(index=False)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
