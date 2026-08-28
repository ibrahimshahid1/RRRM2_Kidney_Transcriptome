#!/usr/bin/env python3
"""Matched-panel annotation null for the glomerular barrier transcript signal."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions  # noqa: E402
from src.clinical_axes.statistics import (  # noqa: E402
    _batch_reml_mkh,
    _mission_effect_batch,
    _prepare_permutation_blocks,
    genewise_z_scores,
)


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_RESULTS = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"
DEFAULT_ATLAS = (
    REPO / "data/processed/v13_dct_reference/mouse_kidney_atlas/atlas_compartment_expression.tsv.gz"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atlas_covariates(path: Path) -> pd.DataFrame:
    atlas = pd.read_csv(path, sep="\t")
    grouped = atlas.groupby("gene_symbol", sort=False)
    covariates = grouped.agg(
        atlas_mean_cpm=("mean_cpm", "mean"),
        atlas_max_cpm=("mean_cpm", "max"),
        atlas_mean_detection=("source_study_detection_fraction", "mean"),
    )
    breadth = atlas.assign(detected=atlas["mean_cpm"] >= 1.0).groupby(
        "gene_symbol"
    )["detected"].sum()
    covariates["atlas_breadth"] = breadth
    covariates["log1p_atlas_mean_cpm"] = np.log1p(covariates["atlas_mean_cpm"])
    covariates["log1p_atlas_max_cpm"] = np.log1p(covariates["atlas_max_cpm"])
    return covariates


def expression_covariates(missions, universe):
    mean_rows = []
    variance_rows = []
    for mission, data in missions.items():
        values = data.expression.loc[list(universe)]
        mean_rows.append(values.mean(axis=1).rename(mission))
        variance_rows.append(
            np.log(values.var(axis=1, ddof=1) + 1e-8).rename(mission)
        )
    means = pd.concat(mean_rows, axis=1)
    variances = pd.concat(variance_rows, axis=1)
    return pd.DataFrame(
        {
            "median_vst_mean": means.median(axis=1),
            "median_log_vst_variance": variances.median(axis=1),
        }
    )


def robust_standardize(table: pd.DataFrame) -> pd.DataFrame:
    median = table.median(axis=0)
    scale = (table - median).abs().median(axis=0) * 1.4826
    scale = scale.mask(scale <= 1e-12, table.std(axis=0, ddof=1))
    scale = scale.mask(scale <= 1e-12, 1.0)
    return (table - median) / scale


def nearest_pools(covariates, targets, pool_size):
    standardized = robust_standardize(covariates)
    candidates = [gene for gene in standardized.index if gene not in set(targets)]
    matrix = standardized.loc[candidates].to_numpy(dtype=float)
    pools = {}
    for gene in targets:
        vector = standardized.loc[gene].to_numpy(dtype=float)
        distance = np.sqrt(np.sum((matrix - vector) ** 2, axis=1))
        order = np.argsort(distance)[:pool_size]
        pools[gene] = np.asarray(candidates, dtype=object)[order]
    return pools


def draw_sets(pools, n_draws, seed):
    rng = np.random.default_rng(seed)
    targets = list(pools)
    draws = []
    for _ in range(n_draws):
        chosen = []
        for gene in targets:
            available = [candidate for candidate in pools[gene] if candidate not in chosen]
            if not available:
                raise RuntimeError(f"Matching pool exhausted for {gene}")
            chosen.append(available[rng.integers(0, len(available))])
        draws.append(chosen)
    return np.asarray(draws, dtype=object)


def combined_gene_z(missions, universe):
    frames = []
    design_frames = []
    for mission, data in missions.items():
        z = genewise_z_scores(data.expression.loc[list(universe)])
        z = z.T
        z.index = pd.Index(
            [f"{mission}::{sample}" for sample in z.index], name="sample"
        )
        design = pd.DataFrame(
            {
                "mission": mission,
                "stratum": data.metadata["block"].to_numpy(),
                "group": data.metadata["condition"].map(
                    {"FLT": "flight", "GC": "control"}
                ).to_numpy(),
            },
            index=z.index,
        )
        frames.append(z)
        design_frames.append(design)
    return pd.concat(frames), pd.concat(design_frames)


def run(args):
    config = yaml.safe_load(args.config.read_text())
    missions = load_primary_missions(config, REPO)
    threshold = float(config["eligibility"]["cpm_threshold"])
    target = list(
        config["primary_family"]["glomerular_barrier_identity_loss"]["subdomains"]
        ["barrier_core"]["genes"]
    )
    common = set.intersection(
        *(cpm_eligible_genes(data, threshold) for data in missions.values())
    )
    common &= set.intersection(
        *(set(data.expression.index) for data in missions.values())
    )
    atlas = atlas_covariates(args.atlas)
    universe = sorted(common.intersection(atlas.index))
    if not set(target).issubset(universe):
        raise RuntimeError(
            "Target genes absent from matched universe: "
            + ",".join(sorted(set(target) - set(universe)))
        )
    expression = expression_covariates(missions, universe)
    covariates = expression.join(
        atlas[
            [
                "log1p_atlas_mean_cpm",
                "log1p_atlas_max_cpm",
                "atlas_mean_detection",
                "atlas_breadth",
            ]
        ],
        how="inner",
    ).dropna()
    universe = sorted(covariates.index)
    pools = nearest_pools(covariates, target, args.pool_size)
    draws = draw_sets(pools, args.draws, int(config["seed"]) + 200)

    gene_z, design = combined_gene_z(missions, universe)
    gene_index = {gene: index for index, gene in enumerate(gene_z.columns)}
    draw_indices = np.asarray(
        [[gene_index[gene] for gene in draw] for draw in draws], dtype=int
    )
    target_indices = np.asarray([[gene_index[gene] for gene in target]], dtype=int)
    all_indices = np.vstack([target_indices, draw_indices])
    values = gene_z.to_numpy(dtype=float)
    # samples x panels; each row in all_indices defines an equal-size panel.
    score_array = values[:, all_indices].mean(axis=2)
    columns = ["target_barrier_core"] + [f"matched_{i:05d}" for i in range(args.draws)]
    score_frame = pd.DataFrame(score_array, index=gene_z.index, columns=columns)

    blocks, mission_labels, axes = _prepare_permutation_blocks(
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
    estimates, tau2, se, t_values = _batch_reml_mkh(
        effects[0], variances[0]
    )
    target_estimate = estimates[0]
    target_t = t_values[0]
    random_estimates = estimates[1:]
    random_t = t_values[1:]
    summary = pd.DataFrame(
        [
            {
                "target": "barrier_core_expression_increase",
                "n_target_genes": len(target),
                "n_universe_genes": len(universe),
                "n_matched_sets": args.draws,
                "nearest_pool_size_per_target": args.pool_size,
                "target_meta_estimate": target_estimate,
                "target_meta_t_mkh": target_t,
                "matched_one_sided_p": (
                    1 + np.count_nonzero(random_estimates >= target_estimate)
                )
                / (args.draws + 1),
                "matched_two_sided_p": (
                    1 + np.count_nonzero(np.abs(random_estimates) >= abs(target_estimate))
                )
                / (args.draws + 1),
                "matched_t_two_sided_p": (
                    1 + np.count_nonzero(np.abs(random_t) >= abs(target_t))
                )
                / (args.draws + 1),
                "random_estimate_q95": np.quantile(random_estimates, 0.95),
                "random_estimate_q975": np.quantile(random_estimates, 0.975),
            }
        ]
    )

    args.results.mkdir(parents=True, exist_ok=True)
    summary.to_csv(
        args.results / "barrier_matched_panel_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        {
            "matched_set": columns[1:],
            "estimate": random_estimates,
            "tau2": tau2[1:],
            "standard_error_mkh": se[1:],
            "t_mkh": random_t,
            "genes": ["|".join(draw) for draw in draws],
        }
    ).to_csv(
        args.results / "barrier_matched_panel_null.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pool_rows = []
    for gene, pool in pools.items():
        target_cov = covariates.loc[gene]
        for rank, candidate in enumerate(pool, start=1):
            row = {
                "target_gene": gene,
                "rank": rank,
                "candidate_gene": candidate,
            }
            for column in covariates:
                row[f"target_{column}"] = target_cov[column]
                row[f"candidate_{column}"] = covariates.loc[candidate, column]
            pool_rows.append(row)
    pd.DataFrame(pool_rows).to_csv(
        args.results / "barrier_matching_pools.tsv", sep="\t", index=False
    )
    manifest = {
        "analysis": "barrier-core observability and kidney-breadth matched-panel null",
        "status": "adversarial_secondary",
        "config": str(args.config),
        "config_sha256": sha256(args.config),
        "atlas": str(args.atlas),
        "atlas_sha256": sha256(args.atlas),
        "n_draws": int(args.draws),
        "pool_size": int(args.pool_size),
        "seed": int(config["seed"]) + 200,
        "matching_covariates": list(covariates.columns),
        "limitation": (
            "nearest pools do not perfectly reproduce the extreme atlas abundance "
            "of the barrier markers; this is a supportive annotation null"
        ),
    }
    (args.results / "barrier_matched_panel_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(summary.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--atlas", type=Path, default=DEFAULT_ATLAS)
    parser.add_argument("--draws", type=int, default=10_000)
    parser.add_argument("--pool-size", type=int, default=100)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
