#!/usr/bin/env python3
"""Strict atlas-selection-matched audit of the podocyte RNA program.

This is a post-hoc adversarial sensitivity analysis. It replaces the original
all-gene nearest-neighbour null with candidates selected by the *same frozen
high-specificity atlas rule* as the podocyte target. Matching uses no flight
labels or flight-effect estimates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import sparse
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.clinical_axes.run_matched_panel_null import (  # noqa: E402
    combined_gene_z,
    expression_covariates,
)
from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions  # noqa: E402
from src.clinical_axes.matching import (  # noqa: E402
    cross_group_nearest_distances,
    draw_balanced_unique_panels,
    euclidean_distance_matrix,
    optimal_unique_assignment,
    robust_standardize,
)
from src.clinical_axes.statistics import (  # noqa: E402
    _batch_reml_mkh,
    _mission_effect_batch,
    _prepare_permutation_blocks,
    blocked_meta_permutation,
    random_effects_reml_mkh,
)


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_GENE_SETS = (
    REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
)
DEFAULT_RESULTS = (
    REPO
    / "data/results/run_20260811_clinical_renal_axes_cross_mission"
    / "strict_podocyte_matching"
)
MATCHING_COLUMNS = (
    "median_vst_mean",
    "median_log_vst_variance",
    "log1p_target_cpm",
    "log1p_max_other_cpm",
    "log2_target_to_max_other",
    "target_source_study_fraction",
    "n_non_target_compartments_detected",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _reference_covariates(missions, common_genes, gene_sets_path: Path):
    table = pd.read_csv(gene_sets_path, sep="\t")
    table = table[
        table["final_for_testing"].astype(bool)
        & table["tier"].eq("high_specificity")
    ].copy()
    if table["gene_symbol"].duplicated().any():
        raise RuntimeError("a gene occurs in multiple frozen high-specificity sets")
    table = table.set_index("gene_symbol")
    expression = expression_covariates(missions, sorted(common_genes))
    covariates = table.join(expression, how="inner")
    covariates["log1p_target_cpm"] = np.log1p(covariates["target_cpm"])
    covariates["log1p_max_other_cpm"] = np.log1p(covariates["max_other_cpm"])
    covariates = covariates.dropna(subset=list(MATCHING_COLUMNS))
    target = sorted(
        covariates.index[covariates["report_compartment"].eq("podocyte")]
    )
    candidates = sorted(
        covariates.index[~covariates["report_compartment"].eq("podocyte")]
    )
    if len(target) < 8 or len(candidates) < len(target):
        raise RuntimeError("strict target/candidate reference is not evaluable")
    return covariates, target, candidates


def _panel_meta_statistics(gene_z, design, panels):
    gene_index = {gene: index for index, gene in enumerate(gene_z.columns)}
    rows: list[int] = []
    columns: list[int] = []
    weights: list[float] = []
    for panel_index, genes in enumerate(panels):
        weight = 1.0 / len(genes)
        for gene in genes:
            rows.append(gene_index[gene])
            columns.append(panel_index)
            weights.append(weight)
    membership = sparse.csc_matrix(
        (weights, (rows, columns)),
        shape=(len(gene_z.columns), len(panels)),
    )
    score_array = gene_z.to_numpy(dtype=float) @ membership
    score_frame = pd.DataFrame(
        np.asarray(score_array),
        index=gene_z.index,
        columns=[f"panel_{index:05d}" for index in range(len(panels))],
    )
    blocks, mission_labels, _ = _prepare_permutation_blocks(
        score_frame,
        design,
        mission_col="mission",
        stratum_col="stratum",
        group_col="group",
        treatment_label="flight",
        control_label="control",
    )
    masks = [block.observed_treatment[None, :] for block in blocks]
    effects, variances = _mission_effect_batch(blocks, mission_labels, masks)
    estimate, tau2, standard_error, t_value = _batch_reml_mkh(
        effects[0], variances[0]
    )
    return estimate, tau2, standard_error, t_value


def _observed_meta(score_frame, design):
    blocks, mission_labels, axis_labels = _prepare_permutation_blocks(
        score_frame,
        design,
        mission_col="mission",
        stratum_col="stratum",
        group_col="group",
        treatment_label="flight",
        control_label="control",
    )
    masks = [block.observed_treatment[None, :] for block in blocks]
    effects, variances = _mission_effect_batch(blocks, mission_labels, masks)
    meta_rows = []
    mission_rows = []
    for axis_index, axis in enumerate(axis_labels):
        fit = random_effects_reml_mkh(
            effects[0, axis_index, :], variances[0, axis_index, :]
        )
        meta_rows.append(
            {
                "axis": axis,
                "estimate": fit.estimate,
                "tau2": fit.tau2,
                "standard_error_mkh": fit.standard_error,
                "t_mkh": fit.t,
                "df": fit.df,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh": fit.p,
                "i_squared": fit.i_squared,
            }
        )
        for mission_index, mission in enumerate(mission_labels):
            mission_rows.append(
                {
                    "axis": axis,
                    "mission": mission,
                    "estimate": effects[0, axis_index, mission_index],
                    "variance": variances[0, axis_index, mission_index],
                }
            )
    return pd.DataFrame(meta_rows), pd.DataFrame(mission_rows)


def run(args):
    config = yaml.safe_load(args.config.read_text())
    missions = load_primary_missions(config, REPO)
    threshold = float(config["eligibility"]["cpm_threshold"])
    common = set.intersection(
        *(cpm_eligible_genes(data, threshold) for data in missions.values())
    )
    common &= set.intersection(
        *(set(data.expression.index) for data in missions.values())
    )
    covariates, target, candidates = _reference_covariates(
        missions, common, args.gene_sets
    )
    standardized = robust_standardize(covariates[list(MATCHING_COLUMNS)])
    target_values = standardized.loc[target].to_numpy(dtype=float)
    candidate_values = standardized.loc[candidates].to_numpy(dtype=float)
    distance = euclidean_distance_matrix(target_values, candidate_values)
    gene_z, design = combined_gene_z(missions, sorted(set(target + candidates)))

    args.results.mkdir(parents=True, exist_ok=True)
    covariate_output = covariates.copy()
    covariate_output.insert(
        0,
        "matching_role",
        np.where(
            covariate_output.index.isin(target),
            "podocyte_target",
            "other_high_specificity_candidate",
        ),
    )
    for column in MATCHING_COLUMNS:
        covariate_output[f"standardized_{column}"] = standardized[column]
    covariate_output.reset_index().to_csv(
        args.results / "strict_podocyte_matching_covariates.tsv",
        sep="\t",
        index=False,
    )

    null_rows = []
    summary_rows = []
    balance_rows = []
    draw_seed_base = int(config["seed"]) + 410
    for pool_size in args.pool_sizes:
        draw_indices, attempts = draw_balanced_unique_panels(
            distance,
            target_values,
            candidate_values,
            pool_size=pool_size,
            n_draws=args.draws,
            seed=draw_seed_base + pool_size,
            balance_caliper=args.balance_caliper,
        )
        draws = np.asarray(candidates, dtype=object)[draw_indices]
        estimates, tau2, standard_errors, t_values = _panel_meta_statistics(
            gene_z, design, [target] + draws.tolist()
        )
        target_estimate = estimates[0]
        target_t = t_values[0]
        null_estimate = estimates[1:]
        null_t = t_values[1:]
        summary_rows.append(
            {
                "analysis": "same_tier_balanced_matched_panel_null",
                "pool_size_per_target": pool_size,
                "balance_caliper_robust_sd": args.balance_caliper,
                "n_target_genes": len(target),
                "n_candidate_genes": len(candidates),
                "n_matched_panels": len(draws),
                "draw_attempts": attempts,
                "acceptance_fraction": len(draws) / attempts,
                "target_meta_estimate": target_estimate,
                "target_meta_t_mkh": target_t,
                "matched_one_sided_p": (
                    1 + np.count_nonzero(null_estimate >= target_estimate)
                )
                / (len(draws) + 1),
                "matched_two_sided_p": (
                    1 + np.count_nonzero(
                        np.abs(null_estimate) >= abs(target_estimate)
                    )
                )
                / (len(draws) + 1),
                "matched_t_two_sided_p": (
                    1 + np.count_nonzero(np.abs(null_t) >= abs(target_t))
                )
                / (len(draws) + 1),
                "null_estimate_q95": np.quantile(null_estimate, 0.95),
                "null_estimate_q975": np.quantile(null_estimate, 0.975),
            }
        )
        panel_means = candidate_values[draw_indices].mean(axis=1)
        differences = panel_means - target_values.mean(axis=0)
        for column_index, column in enumerate(MATCHING_COLUMNS):
            balance_rows.append(
                {
                    "analysis": "same_tier_balanced_matched_panel_null",
                    "scheme": f"nearest_pool_{pool_size}",
                    "covariate": column,
                    "target_mean_robust_sd": target_values[:, column_index].mean(),
                    "mean_panel_difference_robust_sd": differences[
                        :, column_index
                    ].mean(),
                    "median_abs_panel_difference_robust_sd": np.median(
                        np.abs(differences[:, column_index])
                    ),
                    "q95_abs_panel_difference_robust_sd": np.quantile(
                        np.abs(differences[:, column_index]), 0.95
                    ),
                    "maximum_abs_panel_difference_robust_sd": np.max(
                        np.abs(differences[:, column_index])
                    ),
                }
            )
        for draw_index, (genes, estimate, tau, se, t_value) in enumerate(
            zip(draws, null_estimate, tau2[1:], standard_errors[1:], null_t)
        ):
            null_rows.append(
                {
                    "scheme": f"nearest_pool_{pool_size}",
                    "matched_set": f"matched_{draw_index:05d}",
                    "estimate": estimate,
                    "tau2": tau,
                    "standard_error_mkh": se,
                    "t_mkh": t_value,
                    "genes": "|".join(genes),
                }
            )

    pd.DataFrame(summary_rows).to_csv(
        args.results / "strict_podocyte_matched_panel_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(null_rows).to_csv(
        args.results / "strict_podocyte_matched_panel_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    candidate_groups = covariates.loc[candidates, "report_compartment"].to_numpy()
    candidate_nearest = cross_group_nearest_distances(
        candidate_values, candidate_groups
    )
    support_caliper = float(
        np.quantile(candidate_nearest, args.common_support_quantile)
    )
    target_nearest = distance.min(axis=1)
    schemes = {
        "full_common_observable": np.arange(len(target)),
        f"common_support_q{int(args.common_support_quantile * 1000):03d}": np.flatnonzero(
            target_nearest <= support_caliper
        ),
    }
    optimal_match_rows = []
    contrast_scores = {}
    descriptive_scores = {}
    for scheme, target_indices in schemes.items():
        subdistance = distance[target_indices]
        rows, matched_indices = optimal_unique_assignment(subdistance)
        if not np.array_equal(rows, np.arange(len(target_indices))):
            raise RuntimeError("unexpected target ordering from optimal assignment")
        target_genes = [target[index] for index in target_indices]
        matched_genes = [candidates[index] for index in matched_indices]
        target_score = gene_z[target_genes].mean(axis=1)
        matched_score = gene_z[matched_genes].mean(axis=1)
        contrast_scores[scheme] = target_score - matched_score
        descriptive_scores[f"{scheme}__target"] = target_score
        descriptive_scores[f"{scheme}__matched"] = matched_score
        descriptive_scores[f"{scheme}__difference"] = target_score - matched_score
        match_distances = subdistance[rows, matched_indices]
        matched_covariates = candidate_values[matched_indices]
        target_covariates = target_values[target_indices]
        mean_difference = matched_covariates.mean(axis=0) - target_covariates.mean(
            axis=0
        )
        for covariate_index, column in enumerate(MATCHING_COLUMNS):
            balance_rows.append(
                {
                    "analysis": "optimal_same_tier_contrast",
                    "scheme": scheme,
                    "covariate": column,
                    "target_mean_robust_sd": target_covariates[
                        :, covariate_index
                    ].mean(),
                    "mean_panel_difference_robust_sd": mean_difference[
                        covariate_index
                    ],
                    "median_abs_panel_difference_robust_sd": np.nan,
                    "q95_abs_panel_difference_robust_sd": np.nan,
                    "maximum_abs_panel_difference_robust_sd": np.nan,
                }
            )
        for target_position, (target_gene, matched_gene, match_distance) in enumerate(
            zip(target_genes, matched_genes, match_distances)
        ):
            optimal_match_rows.append(
                {
                    "scheme": scheme,
                    "target_gene": target_gene,
                    "matched_gene": matched_gene,
                    "matched_compartment": covariates.loc[
                        matched_gene, "report_compartment"
                    ],
                    "robust_euclidean_distance": match_distance,
                    "target_nearest_candidate_distance": target_nearest[
                        target_indices[target_position]
                    ],
                    "common_support_caliper": (
                        support_caliper
                        if scheme != "full_common_observable"
                        else np.nan
                    ),
                }
            )

    contrast_frame = pd.DataFrame(contrast_scores, index=gene_z.index)
    contrast_result = blocked_meta_permutation(
        contrast_frame,
        design,
        n_permutations=args.permutations,
        seed=int(config["seed"]) + 411,
        chunk_size=args.chunk_size,
    )
    descriptive_meta, descriptive_missions = _observed_meta(
        pd.DataFrame(descriptive_scores, index=gene_z.index), design
    )
    difference_lookup = contrast_result.observed_meta.reset_index().rename(
        columns={"axis": "scheme"}
    )
    difference_lookup = difference_lookup[
        [
            "scheme",
            "empirical_p_two_sided",
            "max_t_fwer",
        ]
    ]
    descriptive_meta["scheme"] = descriptive_meta["axis"].str.rsplit(
        "__", n=1
    ).str[0]
    descriptive_meta["score_type"] = descriptive_meta["axis"].str.rsplit(
        "__", n=1
    ).str[1]
    descriptive_meta = descriptive_meta.merge(
        difference_lookup, on="scheme", how="left"
    )
    descriptive_meta.loc[
        descriptive_meta["score_type"] != "difference",
        ["empirical_p_two_sided", "max_t_fwer"],
    ] = np.nan
    descriptive_meta["n_target_genes"] = descriptive_meta["scheme"].map(
        {scheme: len(indices) for scheme, indices in schemes.items()}
    )
    descriptive_meta["n_candidate_genes"] = len(candidates)
    descriptive_meta.to_csv(
        args.results / "strict_podocyte_optimal_contrast_meta.tsv",
        sep="\t",
        index=False,
    )
    descriptive_missions.to_csv(
        args.results / "strict_podocyte_optimal_contrast_mission_effects.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(optimal_match_rows).to_csv(
        args.results / "strict_podocyte_optimal_matches.tsv",
        sep="\t",
        index=False,
    )
    contrast_result.null_t.to_csv(
        args.results / "strict_podocyte_optimal_contrast_null_t.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.DataFrame(balance_rows).to_csv(
        args.results / "strict_podocyte_matching_balance.tsv",
        sep="\t",
        index=False,
    )

    manifest = {
        "analysis": "strict high-specificity-tier matched audit of podocyte RNA",
        "status": "post_hoc_adversarial_sensitivity",
        "config": str(args.config),
        "config_sha256": sha256(args.config),
        "gene_sets": str(args.gene_sets),
        "gene_sets_sha256": sha256(args.gene_sets),
        "n_common_observable_podocyte_targets": len(target),
        "n_common_observable_other_tier_candidates": len(candidates),
        "matching_covariates": list(MATCHING_COLUMNS),
        "candidate_rule": (
            "same frozen high-specificity atlas tier, non-podocyte compartment, "
            "observable in all five primary missions"
        ),
        "pool_sizes": list(args.pool_sizes),
        "n_draws_per_pool": args.draws,
        "aggregate_balance_caliper_robust_sd": args.balance_caliper,
        "common_support_definition": (
            "retain a podocyte target when its nearest non-podocyte candidate "
            "distance does not exceed the declared quantile of cross-compartment "
            "nearest-neighbour distances among candidate markers"
        ),
        "common_support_quantile": args.common_support_quantile,
        "common_support_caliper": support_caliper,
        "n_label_permutations": args.permutations,
        "draw_seed_base": draw_seed_base,
        "label_permutation_seed": int(config["seed"]) + 411,
        "interpretation": (
            "tests whether the podocyte high-specificity panel is unusual relative "
            "to equally selected, flight-label-blind markers; it does not establish "
            "podocyte cell counts, localization, injury, or function"
        ),
    }
    (args.results / "strict_podocyte_matching_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print(pd.DataFrame(summary_rows).to_string(index=False))
    print("\nOptimal same-tier contrasts")
    print(
        descriptive_meta[descriptive_meta["score_type"].eq("difference")][
            [
                "scheme",
                "n_target_genes",
                "estimate",
                "ci_low_mkh",
                "ci_high_mkh",
                "p_mkh",
                "empirical_p_two_sided",
                "max_t_fwer",
            ]
        ].to_string(index=False)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--gene-sets", type=Path, default=DEFAULT_GENE_SETS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--draws", type=int, default=10_000)
    parser.add_argument("--pool-sizes", type=int, nargs="+", default=[10, 25])
    parser.add_argument("--balance-caliper", type=float, default=0.25)
    parser.add_argument("--common-support-quantile", type=float, default=0.95)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--chunk-size", type=int, default=256)
    args = parser.parse_args()
    if not 0 < args.common_support_quantile <= 1:
        parser.error("--common-support-quantile must lie in (0, 1]")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
