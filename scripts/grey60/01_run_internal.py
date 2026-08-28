#!/usr/bin/env python3
"""Run the locked internal Grey60 adversarial tests."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.outliers_influence import variance_inflation_factor
import yaml

from src.grey60.adversarial import (
    age_specific_effects,
    familywise_p,
    hedges_g,
    max_t_permutation,
    mean_z_score,
    median_z_score,
    observed_saturated_t,
    pooled_iss_effect,
    stratified_bootstrap_iss_effect,
    weighted_mean_z_score,
    zscore_rows,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot JSON-serialize {type(value).__name__}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else REPO_ROOT / p


def load_expression(path: Path) -> pd.DataFrame:
    sep = "\t" if path.name.endswith((".tsv", ".tsv.gz")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df


def load_metadata(path: Path) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t")
    sample_col = next(
        c
        for c in ("Sample Name (raw_counts_colname)", "Sample Name", "sample")
        if c in meta.columns
    )
    meta = meta.set_index(sample_col, drop=False)
    meta["EnvGroup"] = meta["EnvGroup"].replace({"HGC": "GC", "VGC": "VIV"})
    meta["Age"] = meta["Age"].replace(
        {"Young": "YNG", "young": "YNG", "Yng": "YNG", "old": "OLD"}
    )
    meta["Arm"] = meta["Arm"].replace(
        {"ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T"}
    )
    return meta


def orient_to_reference(score: pd.Series, reference: pd.Series) -> tuple[pd.Series, float, float]:
    common = score.index.intersection(reference.index)
    corr = float(stats.pearsonr(score.loc[common], reference.loc[common]).statistic)
    direction = 1.0 if corr >= 0 else -1.0
    return score * direction, direction, corr


def score_effect_rows(scores: dict[str, pd.Series], meta: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for name, score in scores.items():
        for arm in ("ISS-T", "LAR"):
            effects = age_specific_effects(score, meta, arm=arm)
            rows.append(
                {
                    "score": name,
                    "arm": arm,
                    "age": "pooled",
                    "flight_minus_control": pooled_iss_effect(score, meta, arm=arm),
                    "n_flight": 10,
                    "n_control": 10,
                }
            )
            for age, effect in effects.items():
                rows.append(
                    {
                        "score": name,
                        "arm": arm,
                        "age": age,
                        "flight_minus_control": effect,
                        "n_flight": 5,
                        "n_control": 5,
                    }
                )
    return pd.DataFrame(rows)


def leave_one_mouse(
    score: pd.Series, meta: pd.DataFrame, arm: str = "ISS-T"
) -> pd.DataFrame:
    rows = []
    eligible = meta.index[meta["Arm"] == arm]
    for sample in eligible:
        keep = meta.index != sample
        effect = pooled_iss_effect(score.loc[keep], meta.loc[keep], arm=arm)
        rows.append({"omitted_sample": sample, "effect": effect})
    return pd.DataFrame(rows)


def matched_subset_null(
    target_genes: list[str],
    universe_genes: list[str],
    reference_expression: pd.DataFrame,
    per_gene_effect: pd.Series,
    *,
    n_draws: int,
    seed: int,
    exclude_target: bool,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Draw mean/variance-stratified gene sets from a declared universe."""
    genes = [
        g
        for g in universe_genes
        if g in reference_expression.index and g in per_gene_effect.index
    ]
    target = [g for g in target_genes if g in genes]
    if len(target) < 3:
        raise ValueError("Matched-set target has fewer than three eligible genes")
    ref = reference_expression.loc[genes].to_numpy(dtype=float)
    table = pd.DataFrame(
        {
            "gene": genes,
            "mean": np.nanmean(ref, axis=1),
            "log_variance": np.log(np.nanvar(ref, axis=1, ddof=1) + 1e-8),
        }
    ).set_index("gene")
    # Rank-based bins avoid duplicate-edge failures and lock approximately
    # equal-sized observability strata.
    table["mean_bin"] = pd.qcut(
        table["mean"].rank(method="first"), 5, labels=False
    )
    table["variance_bin"] = pd.qcut(
        table["log_variance"].rank(method="first"), 5, labels=False
    )
    table["bin"] = table["mean_bin"].astype(str) + "_" + table[
        "variance_bin"
    ].astype(str)
    target_counts = table.loc[target, "bin"].value_counts().to_dict()
    pools: dict[str, np.ndarray] = {}
    target_set = set(target)
    for bin_name, needed in target_counts.items():
        candidates = table.index[table["bin"] == bin_name].to_numpy(dtype=str)
        if exclude_target:
            candidates = np.array([g for g in candidates if g not in target_set])
        if len(candidates) < needed:
            candidates = table.index[table["bin"] == bin_name].to_numpy(dtype=str)
        if len(candidates) < needed:
            raise ValueError(f"Insufficient genes in matching bin {bin_name}")
        pools[bin_name] = candidates

    rng = np.random.default_rng(seed)
    out = np.empty(n_draws, dtype=float)
    for i in range(n_draws):
        sampled: list[str] = []
        for bin_name, needed in target_counts.items():
            sampled.extend(
                rng.choice(pools[bin_name], size=needed, replace=False).tolist()
            )
        out[i] = float(per_gene_effect.loc[sampled].mean())
    audit = table.copy()
    audit["is_target"] = audit.index.isin(target)
    return out, audit.reset_index()


def fit_score_model(
    score: pd.Series,
    meta: pd.DataFrame,
    nuisance: pd.DataFrame | None = None,
) -> tuple[sm.regression.linear_model.RegressionResultsWrapper, pd.DataFrame]:
    design = pd.DataFrame(index=meta.index)
    design["Flight"] = (meta["EnvGroup"] == "FLT").astype(float)
    design["AgeOld"] = (meta["Age"] == "OLD").astype(float)
    design["ArmLAR"] = (meta["Arm"] == "LAR").astype(float)
    design["Flight_ArmLAR"] = design["Flight"] * design["ArmLAR"]
    if nuisance is not None:
        design = pd.concat([design, nuisance.loc[design.index]], axis=1)
    design = sm.add_constant(design.astype(float), has_constant="add")
    fit = sm.OLS(score.loc[design.index].to_numpy(dtype=float), design).fit(
        cov_type="HC3"
    )
    return fit, design


def blocked_freedman_lane(
    y: pd.Series,
    full_design: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    tested_column: str,
    n_perm: int,
    seed: int,
) -> float:
    reduced = full_design.drop(columns=[tested_column])
    # The flight-by-arm column is undefined without Flight and is removed for
    # the reduced nuisance fit when testing the ISS-terminal flight coefficient.
    if tested_column == "Flight" and "Flight_ArmLAR" in reduced:
        reduced = reduced.drop(columns=["Flight_ArmLAR"])
    yr = y.loc[full_design.index].to_numpy(dtype=float)
    xr = reduced.to_numpy(dtype=float)
    xf = full_design.to_numpy(dtype=float)
    beta_r = np.linalg.pinv(xr) @ yr
    fitted = xr @ beta_r
    residual = yr - fitted
    pinv_f = np.linalg.pinv(xf)
    col_idx = list(full_design.columns).index(tested_column)
    observed = float((pinv_f @ yr)[col_idx])

    blocks = []
    for _, idx in metadata.loc[full_design.index].groupby(["Age", "Arm"]).groups.items():
        blocks.append(np.array([full_design.index.get_loc(x) for x in idx], dtype=int))
    rng = np.random.default_rng(seed)
    exceed = 0
    for _ in range(n_perm):
        rp = residual.copy()
        for idx in blocks:
            rp[idx] = residual[rng.permutation(idx)]
        beta = pinv_f @ (fitted + rp)
        exceed += int(abs(beta[col_idx]) >= abs(observed))
    return float((exceed + 1) / (n_perm + 1))


def composition_nuisance(
    meta: pd.DataFrame,
    clr_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    clr = pd.read_csv(clr_path, index_col=0).loc[meta.index]
    xclr = clr.to_numpy(dtype=float)
    xclr = (xclr - xclr.mean(axis=0)) / np.where(
        xclr.std(axis=0, ddof=1) < 1e-12, 1.0, xclr.std(axis=0, ddof=1)
    )
    u, s, _ = np.linalg.svd(xclr, full_matrices=False)
    pcs = pd.DataFrame(
        {"CompositionPC1": u[:, 0] * s[0], "CompositionPC2": u[:, 1] * s[1]},
        index=meta.index,
    )

    tech = pd.DataFrame(index=meta.index)
    batch = pd.get_dummies(meta["LibraryBatch"], prefix="Batch", drop_first=True)
    tech = pd.concat([tech, batch.astype(float)], axis=1)
    continuous = {
        "log_ReadDepth": np.log(pd.to_numeric(meta["ReadDepth"], errors="coerce")),
        "rRNA": pd.to_numeric(meta["rRNA"], errors="coerce"),
        "QA": pd.to_numeric(meta["QA"], errors="coerce"),
    }
    for name, values in continuous.items():
        values = values.fillna(values.median())
        sd = values.std(ddof=1)
        tech[name] = (values - values.mean()) / (sd if sd > 1e-12 else 1.0)
    nuisance = pd.concat([tech, pcs], axis=1)
    return nuisance, clr


def module_mean_scores(
    expression: pd.DataFrame,
    assignments: pd.DataFrame,
    samples: list[str],
) -> pd.DataFrame:
    scores = {}
    for color, sub in assignments[assignments["module_color"] != "grey"].groupby(
        "module_color"
    ):
        genes = [g for g in sub["gene"].astype(str) if g in expression.index]
        if len(genes) >= 3:
            scores[color] = mean_z_score(expression.loc[:, samples], genes)
    return pd.DataFrame(scores, index=samples)


def load_atlas_marker_genes(
    atlas_path: Path,
    id_map: pd.DataFrame,
    name: str,
) -> list[str]:
    """Load a flight-blind whole-kidney atlas comparator set as Ensembl IDs."""
    atlas = pd.read_csv(atlas_path, sep="\t")
    selected = atlas[
        (atlas["gene_set"] == name.casefold())
        & (atlas["status"] == "defined")
    ]
    symbol_to_ensembl = (
        id_map.dropna(subset=["mgi_symbol", "ensembl_gene_id"])
        .drop_duplicates("mgi_symbol")
        .set_index("mgi_symbol")["ensembl_gene_id"]
        .astype(str)
        .to_dict()
    )
    genes = [
        symbol_to_ensembl[symbol]
        for symbol in selected["gene_symbol"].astype(str)
        if symbol in symbol_to_ensembl
    ]
    if len(genes) < 50:
        raise ValueError(f"Atlas marker set {name} mapped only {len(genes)} genes")
    return list(dict.fromkeys(genes))


def run(args: argparse.Namespace) -> None:
    config_path = resolve(args.config)
    cfg = yaml.safe_load(config_path.read_text())
    inputs = {k: resolve(v) for k, v in cfg["inputs"].items()}
    outdir = resolve(args.outdir or cfg["output_dir"]) / "internal"
    outdir.mkdir(parents=True, exist_ok=True)
    for sub in ("selection", "influence", "composition", "compactness", "specificity"):
        (outdir / sub).mkdir(exist_ok=True)

    # Fail closed if a frozen input changed.
    hash_key_to_input = {
        "grey60_genes_sha256": "grey60_genes",
        "module_assignments_sha256": "module_assignments",
        "module_eigengenes_sha256": "module_eigengenes",
        "rrrm2_rtech_sha256": "rrrm2_rtech",
        "rrrm2_vst_sha256": "rrrm2_vst",
        "rrrm2_clr_sha256": "rrrm2_clr",
        "atlas_comparator_sets_sha256": "atlas_comparator_sets",
    }
    observed_hashes = {}
    for hash_key, input_key in hash_key_to_input.items():
        observed_hashes[hash_key] = sha256(inputs[input_key])
        expected = cfg["frozen_hashes"][hash_key]
        if observed_hashes[hash_key] != expected:
            raise RuntimeError(
                f"Frozen input hash mismatch for {input_key}: "
                f"{observed_hashes[hash_key]} != {expected}"
            )

    genes = [g.strip() for g in inputs["grey60_genes"].read_text().splitlines() if g.strip()]
    if len(genes) != 48 or len(set(genes)) != 48:
        raise RuntimeError("Frozen Grey60 membership is not exactly 48 unique genes")

    meta_all = load_metadata(inputs["rrrm2_meta"])
    rtech = load_expression(inputs["rrrm2_rtech"])
    vst = load_expression(inputs["rrrm2_vst"])
    samples_fg = [
        s for s in rtech.columns if s in meta_all.index and meta_all.loc[s, "EnvGroup"] in ("FLT", "GC")
    ]
    meta = meta_all.loc[samples_fg].copy()
    if len(meta) != 40:
        raise RuntimeError(f"Expected 40 FLT+GC samples; got {len(meta)}")

    assignments = pd.read_csv(inputs["module_assignments"])
    assignments["gene"] = assignments["gene"].astype(str)
    mes = pd.read_csv(inputs["module_eigengenes"], index_col=0).loc[samples_fg]
    module_map = (
        assignments.groupby("module_num", as_index=True)["module_color"].first().to_dict()
    )
    me17 = mes["ME17"].copy()
    if module_map.get(17) != "grey60":
        raise RuntimeError("Historical ME17 no longer maps to Grey60")
    non_grey_mes = pd.DataFrame(
        {
            me: mes[me]
            for me in mes.columns
            if module_map.get(int(me.replace("ME", ""))) != "grey"
        },
        index=samples_fg,
    )
    if non_grey_mes.shape[1] != 20:
        raise RuntimeError(f"Expected 20 non-grey eigengenes; got {non_grey_mes.shape[1]}")

    expr_fg = rtech.loc[:, samples_fg]
    kme = pd.read_csv(inputs["grey60_kme"]).set_index("gene")["kME"]
    primary, primary_dir, primary_corr = orient_to_reference(
        mean_z_score(expr_fg, genes), me17
    )
    median, median_dir, median_corr = orient_to_reference(
        median_z_score(expr_fg, genes), me17
    )
    weighted, weighted_dir, weighted_corr = orient_to_reference(
        weighted_mean_z_score(expr_fg, genes, kme), me17
    )
    scores = {
        "mean_z": primary,
        "median_z": median,
        "fixed_kme_weighted_mean_z": weighted,
        "saved_ME17": me17,
    }
    score_effect_rows(scores, meta).to_csv(
        outdir / "score_effects.tsv", sep="\t", index=False
    )

    # Gate A: selection-adjusted null over the exact historical inspection family.
    n_perm = int(cfg["resampling"]["maxT_permutations"])
    null_max = max_t_permutation(
        non_grey_mes,
        meta,
        n_perm=n_perm,
        seed=int(cfg["resampling"]["seed"]),
    )
    np.save(outdir / "selection" / "maxT_null.npy", null_max)
    historical_tests = observed_saturated_t(non_grey_mes, meta)
    historical_tests["maxT_fwer"] = historical_tests["t"].map(
        lambda x: familywise_p(null_max, x)
    )
    historical_tests.to_csv(
        outdir / "selection" / "historical_module_tests.tsv", sep="\t", index=False
    )
    primary_tests = observed_saturated_t(
        pd.DataFrame({"Grey60_mean_z": primary}, index=meta.index), meta
    )
    primary_tests["maxT_fwer"] = primary_tests["t"].map(
        lambda x: familywise_p(null_max, x)
    )
    primary_tests.to_csv(
        outdir / "selection" / "grey60_primary_tests.tsv", sep="\t", index=False
    )
    expected = primary_tests[primary_tests["term"].isin(["ISS-T_YNG", "ISS-T_OLD"])]
    expected_maxT_p = float(expected["maxT_fwer"].min())

    n_boot = int(cfg["resampling"]["animal_bootstrap"])
    boot = stratified_bootstrap_iss_effect(
        primary,
        meta,
        n_boot=n_boot,
        seed=int(cfg["resampling"]["seed"]) + 1,
    )
    pd.DataFrame({"pooled_iss_effect": boot}).to_csv(
        outdir / "selection" / "pooled_iss_bootstrap.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    ci_low, ci_high = np.quantile(boot, [0.025, 0.975])
    age_effect = age_specific_effects(primary, meta)

    # Gate C: animal and gene influence.
    loo_mouse = leave_one_mouse(primary, meta)
    loo_mouse.to_csv(
        outdir / "influence" / "leave_one_mouse.tsv", sep="\t", index=False
    )
    gene_z = pd.DataFrame(
        zscore_rows(expr_fg.loc[genes].to_numpy(dtype=float)) * primary_dir,
        index=genes,
        columns=samples_fg,
    )
    per_gene_effect = pd.Series(
        {
            gene: pooled_iss_effect(gene_z.loc[gene], meta)
            for gene in genes
        },
        name="pooled_iss_effect",
    )
    gene_influence_rows = []
    for gene in genes:
        score = gene_z.drop(index=gene).mean(axis=0)
        gene_influence_rows.append(
            {"omitted_gene": gene, "effect": pooled_iss_effect(score, meta)}
        )
    loo_gene = pd.DataFrame(gene_influence_rows)
    loo_gene.to_csv(
        outdir / "influence" / "leave_one_gene.tsv", sep="\t", index=False
    )
    id_map = pd.read_csv(inputs["id_map"], sep="\t", comment="#")
    symbol_map = dict(
        zip(id_map["ensembl_gene_id"].astype(str), id_map["mgi_symbol"].astype(str))
    )
    gene_contrasts = per_gene_effect.rename_axis("gene").reset_index()
    gene_contrasts["symbol"] = gene_contrasts["gene"].map(symbol_map)
    abs_sum = float(gene_contrasts["pooled_iss_effect"].abs().sum())
    gene_contrasts["absolute_contribution"] = (
        gene_contrasts["pooled_iss_effect"].abs() / abs_sum
    )
    gene_contrasts.to_csv(
        outdir / "influence" / "gene_contrasts.tsv", sep="\t", index=False
    )
    n_positive = int((per_gene_effect > 0).sum())
    sign_p = float(stats.binomtest(n_positive, len(per_gene_effect), 0.5, alternative="greater").pvalue)
    max_contribution = float(gene_contrasts["absolute_contribution"].max())

    ranked = [g for g in kme.sort_values(key=np.abs, ascending=False).index if g in genes]
    removal_rows = []
    for n_remove in (1, 5, 10):
        kept = [g for g in genes if g not in set(ranked[:n_remove])]
        score = gene_z.loc[kept].mean(axis=0)
        removal_rows.append(
            {
                "removal": f"top_{n_remove}_kME",
                "n_removed": n_remove,
                "effect": pooled_iss_effect(score, meta),
            }
        )
    marker_union: set[str] = set()
    for name in ("Endothelial", "Fibroblast"):
        marker_union |= set(
            load_atlas_marker_genes(
                inputs["atlas_comparator_sets"], id_map, name
            )
        )
    marker_overlap = sorted(set(genes) & marker_union)
    kept = [g for g in genes if g not in marker_union]
    marker_removed_effect = pooled_iss_effect(gene_z.loc[kept].mean(axis=0), meta)
    removal_rows.append(
        {
            "removal": "endothelial_fibroblast_marker_overlap",
            "n_removed": len(marker_overlap),
            "effect": marker_removed_effect,
        }
    )
    pd.DataFrame(removal_rows).to_csv(
        outdir / "influence" / "hub_marker_removals.tsv", sep="\t", index=False
    )
    (outdir / "influence" / "marker_overlap_genes.txt").write_text(
        "\n".join(marker_overlap) + ("\n" if marker_overlap else "")
    )

    # Gate D: raw VST, technical/composition covariates, and baseline falsification.
    vst_fg = vst.loc[:, samples_fg]
    raw_score, raw_dir, raw_corr = orient_to_reference(
        mean_z_score(vst_fg, genes), primary
    )
    nuisance_all, clr_all = composition_nuisance(meta_all, inputs["rrrm2_clr"])
    nuisance = nuisance_all.loc[meta.index]
    raw_unadjusted, design_unadjusted = fit_score_model(raw_score, meta)
    raw_adjusted, design_adjusted = fit_score_model(raw_score, meta, nuisance)
    fl_p = blocked_freedman_lane(
        raw_score,
        design_adjusted,
        meta,
        tested_column="Flight",
        n_perm=int(cfg["resampling"]["external_permutations"]),
        seed=int(cfg["resampling"]["seed"]) + 2,
    )
    vif_rows = []
    x_no_const = design_adjusted.drop(columns="const")
    for i, col in enumerate(x_no_const.columns):
        vif_rows.append(
            {
                "covariate": col,
                "vif": float(variance_inflation_factor(x_no_const.to_numpy(), i)),
            }
        )
    vif = pd.DataFrame(vif_rows)
    vif.to_csv(outdir / "composition" / "vif.tsv", sep="\t", index=False)
    omit_rows = []
    for col in nuisance.columns:
        fit_omit, _ = fit_score_model(raw_score, meta, nuisance.drop(columns=col))
        omit_rows.append(
            {
                "omitted_covariate": col,
                "flight_estimate": float(fit_omit.params["Flight"]),
            }
        )
    pd.DataFrame(omit_rows).to_csv(
        outdir / "composition" / "leave_one_covariate.tsv", sep="\t", index=False
    )
    composition_summary = pd.DataFrame(
        [
            {
                "model": "raw_vst_unadjusted",
                "flight_estimate": float(raw_unadjusted.params["Flight"]),
                "hc3_se": float(raw_unadjusted.bse["Flight"]),
                "hc3_p": float(raw_unadjusted.pvalues["Flight"]),
                "ci_low": float(raw_unadjusted.conf_int().loc["Flight", 0]),
                "ci_high": float(raw_unadjusted.conf_int().loc["Flight", 1]),
                "freedman_lane_p": np.nan,
            },
            {
                "model": "raw_vst_technical_composition",
                "flight_estimate": float(raw_adjusted.params["Flight"]),
                "hc3_se": float(raw_adjusted.bse["Flight"]),
                "hc3_p": float(raw_adjusted.pvalues["Flight"]),
                "ci_low": float(raw_adjusted.conf_int().loc["Flight", 0]),
                "ci_high": float(raw_adjusted.conf_int().loc["Flight", 1]),
                "freedman_lane_p": fl_p,
            },
        ]
    )
    composition_summary.to_csv(
        outdir / "composition" / "model_summary.tsv", sep="\t", index=False
    )
    adjusted_ratio = float(
        raw_adjusted.params["Flight"] / raw_unadjusted.params["Flight"]
    )

    # Baseline pseudo-contrast: VIV versus BSL inside the same four blocks.
    baseline_samples = [
        s
        for s in rtech.columns
        if s in meta_all.index and meta_all.loc[s, "EnvGroup"] in ("BSL", "VIV")
    ]
    baseline_meta = meta_all.loc[baseline_samples].copy()
    baseline_meta["EnvGroup"] = baseline_meta["EnvGroup"].replace(
        {"BSL": "GC", "VIV": "FLT"}
    )
    baseline_expr = rtech.loc[:, baseline_samples]
    baseline_score = mean_z_score(baseline_expr, genes) * primary_dir
    baseline_modules = module_mean_scores(
        baseline_expr, assignments, baseline_samples
    )
    baseline_null = max_t_permutation(
        baseline_modules,
        baseline_meta,
        n_perm=n_perm,
        seed=int(cfg["resampling"]["seed"]) + 3,
    )
    baseline_test = observed_saturated_t(
        pd.DataFrame({"Grey60": baseline_score}), baseline_meta
    )
    baseline_test["maxT_fwer"] = baseline_test["t"].map(
        lambda x: familywise_p(baseline_null, x)
    )
    baseline_test.to_csv(
        outdir / "composition" / "baseline_falsification_tests.tsv",
        sep="\t",
        index=False,
    )
    baseline_interaction = baseline_test.loc[
        baseline_test["term"] == "Flight:Arm"
    ].iloc[0]
    baseline_iss_idx = baseline_meta.index[baseline_meta["Arm"] == "ISS-T"]
    baseline_g = hedges_g(
        baseline_score.loc[
            baseline_iss_idx[
                baseline_meta.loc[baseline_iss_idx, "EnvGroup"] == "FLT"
            ]
        ],
        baseline_score.loc[
            baseline_iss_idx[
                baseline_meta.loc[baseline_iss_idx, "EnvGroup"] == "GC"
            ]
        ],
    )

    # Gate B compactness: existing true GC reference and matched subsets.
    gc_assign = pd.read_csv(inputs["gc_reference_assignments"])
    gc_assign["gene"] = gc_assign["gene"].astype(str)
    overlap_rows = []
    for color, sub in gc_assign.groupby("module_color_gc"):
        module_genes = set(sub["gene"])
        overlap = len(module_genes & set(genes))
        union = len(module_genes | set(genes))
        overlap_rows.append(
            {
                "gc_module": color,
                "module_size": len(module_genes),
                "overlap": overlap,
                "jaccard": overlap / union if union else 0.0,
            }
        )
    gc_overlap = pd.DataFrame(overlap_rows).sort_values(
        ["overlap", "jaccard"], ascending=False
    )
    gc_overlap.to_csv(
        outdir / "compactness" / "gc_reference_overlap.tsv",
        sep="\t",
        index=False,
    )
    best_gc = gc_overlap.iloc[0]
    gc_parent_genes = gc_assign.loc[
        gc_assign["module_color_gc"] == best_gc["gc_module"], "gene"
    ].tolist()
    target_gc = sorted(set(genes) & set(gc_parent_genes))
    remainder_gc = sorted(set(gc_parent_genes) - set(target_gc))
    all_gene_z = pd.DataFrame(
        zscore_rows(expr_fg.loc[[g for g in rtech.index if g in expr_fg.index]].to_numpy(dtype=float)),
        index=[g for g in rtech.index if g in expr_fg.index],
        columns=samples_fg,
    )
    per_gene_effect_all = pd.Series(
        {
            gene: pooled_iss_effect(all_gene_z.loc[gene] * primary_dir, meta)
            for gene in set(gc_parent_genes) | set(assignments["gene"])
            if gene in all_gene_z.index
        }
    )
    reference_samples = [
        s
        for s in vst.columns
        if s in meta_all.index and meta_all.loc[s, "EnvGroup"] in ("BSL", "VIV")
    ]
    gc_null, gc_match_audit = matched_subset_null(
        target_gc,
        gc_parent_genes,
        vst.loc[:, reference_samples],
        per_gene_effect_all,
        n_draws=int(cfg["resampling"]["matched_subsets"]),
        seed=int(cfg["resampling"]["seed"]) + 4,
        exclude_target=False,
    )
    np.save(outdir / "compactness" / "gc_green_matched_null.npy", gc_null)
    gc_match_audit.to_csv(
        outdir / "compactness" / "gc_green_matching_audit.tsv",
        sep="\t",
        index=False,
    )
    target_gc_effect = float(per_gene_effect_all.loc[target_gc].mean())
    remainder_gc_effect = float(per_gene_effect_all.loc[remainder_gc].mean())
    target_gc_percentile = float(np.mean(gc_null <= target_gc_effect))
    pd.DataFrame(
        [
            {
                "best_gc_module": best_gc["gc_module"],
                "best_gc_module_size": int(best_gc["module_size"]),
                "grey60_overlap": int(best_gc["overlap"]),
                "jaccard": float(best_gc["jaccard"]),
                "grey60_subset_effect": target_gc_effect,
                "gc_module_remainder_effect": remainder_gc_effect,
                "matched_null_mean": float(gc_null.mean()),
                "matched_null_q95": float(np.quantile(gc_null, 0.95)),
                "grey60_percentile": target_gc_percentile,
            }
        ]
    ).to_csv(
        outdir / "compactness" / "gc_green_competition.tsv",
        sep="\t",
        index=False,
    )

    # Gate F: independent marker adjustment and matched 48-gene universe.
    marker_scores = {}
    for name in ("Endothelial", "Fibroblast"):
        panel_genes = [
            g
            for g in load_atlas_marker_genes(
                inputs["atlas_comparator_sets"], id_map, name
            )
            if g not in set(genes)
        ]
        marker_scores[name] = mean_z_score(expr_fg, panel_genes)
    marker_df = pd.DataFrame(marker_scores, index=meta.index)
    generic_adjusted, _ = fit_score_model(primary, meta, marker_df)
    primary_unadjusted, _ = fit_score_model(primary, meta)
    generic_ratio = float(
        generic_adjusted.params["Flight"] / primary_unadjusted.params["Flight"]
    )
    whole_universe = [
        g for g in assignments["gene"].tolist() if g in per_gene_effect_all.index
    ]
    generic_null, generic_match_audit = matched_subset_null(
        genes,
        whole_universe,
        vst.loc[:, reference_samples],
        per_gene_effect_all,
        n_draws=int(cfg["resampling"]["matched_subsets"]),
        seed=int(cfg["resampling"]["seed"]) + 5,
        exclude_target=True,
    )
    np.save(outdir / "specificity" / "matched_48gene_null.npy", generic_null)
    generic_match_audit.to_csv(
        outdir / "specificity" / "matched_48gene_audit.tsv",
        sep="\t",
        index=False,
    )
    primary_effect = pooled_iss_effect(primary, meta)
    generic_percentile = float(np.mean(generic_null <= primary_effect))
    pd.DataFrame(
        [
            {
                "unadjusted_flight_estimate": float(primary_unadjusted.params["Flight"]),
                "marker_adjusted_flight_estimate": float(generic_adjusted.params["Flight"]),
                "adjusted_to_unadjusted_ratio": generic_ratio,
                "marker_adjusted_hc3_p": float(generic_adjusted.pvalues["Flight"]),
                "matched_null_mean": float(generic_null.mean()),
                "matched_null_q95": float(np.quantile(generic_null, 0.95)),
                "grey60_effect": primary_effect,
                "grey60_percentile": generic_percentile,
            }
        ]
    ).to_csv(
        outdir / "specificity" / "generic_state_competition.tsv",
        sep="\t",
        index=False,
    )

    thresholds = cfg["go_no_go"]
    gate_a = {
        "maxT": expected_maxT_p <= thresholds["gate_A"]["maxT_fwer_lte"],
        "bootstrap_ci": ci_low > 0,
        "age_directions": all(v > 0 for v in age_effect.values()),
    }
    gate_c = {
        "leave_one_mouse": bool((loo_mouse["effect"] > 0).all()),
        "leave_one_gene": bool((loo_gene["effect"] > 0).all()),
        "score_direction": all(
            pooled_iss_effect(s, meta) > 0 for s in scores.values()
        ),
        "top10_removal": bool(
            pd.DataFrame(removal_rows).loc[
                lambda x: x["removal"] == "top_10_kME", "effect"
            ].iloc[0]
            > 0
        ),
        "marker_removal": marker_removed_effect > 0,
        "gene_signs": (
            n_positive >= thresholds["gate_C"]["positive_gene_contrasts_gte"]
            and sign_p <= thresholds["gate_C"]["sign_test_p_lte"]
        ),
        "max_contribution": (
            max_contribution
            <= thresholds["gate_C"]["maximum_gene_absolute_contribution_lte"]
        ),
    }
    max_vif = float(vif["vif"].replace([np.inf, -np.inf], np.nan).max())
    omit_direction_stable = bool(
        (pd.DataFrame(omit_rows)["flight_estimate"] > 0).all()
    )
    gate_d = {
        "rtech_vst_direction": pooled_iss_effect(primary, meta)
        * float(raw_unadjusted.params["Flight"])
        > 0,
        "adjusted_effect": (
            raw_adjusted.params["Flight"] > 0
            and adjusted_ratio
            >= thresholds["gate_D"]["adjusted_to_unadjusted_ratio_gte"]
        ),
        "vif": max_vif < thresholds["gate_D"]["maximum_vif_lt"],
        "leave_one_covariate": omit_direction_stable,
        "baseline": (
            baseline_interaction["maxT_fwer"]
            > thresholds["gate_D"]["baseline_maxT_p_gt"]
            and abs(baseline_g.estimate)
            < thresholds["gate_D"]["baseline_abs_hedges_g_lt"]
        ),
    }
    gate_f = {
        "adjusted_effect": (
            generic_adjusted.params["Flight"] > 0
            and generic_ratio
            >= thresholds["gate_F"]["adjusted_to_unadjusted_ratio_gte"]
        ),
        "matched_percentile": (
            generic_percentile >= thresholds["gate_F"]["matched_set_percentile_gte"]
        ),
    }
    gate_b_partial = {
        "gc_subset_exceeds_remainder": target_gc_effect > remainder_gc_effect,
        "gc_matched_percentile": (
            target_gc_percentile
            >= thresholds["gate_B"]["matched_subset_percentile_gte"]
        ),
        "flight_blind_grid_pending": True,
    }

    gate_status = {
        "gate_A": {
            "pass": all(gate_a.values()),
            "components": gate_a,
            "expected_iss_strata_min_maxT_p": expected_maxT_p,
            "bootstrap_ci": [float(ci_low), float(ci_high)],
            "age_specific_effects": age_effect,
        },
        "gate_B_partial": {
            "pass_so_far": all(
                value for key, value in gate_b_partial.items() if key != "flight_blind_grid_pending"
            ),
            "components": gate_b_partial,
        },
        "gate_C": {
            "pass": all(gate_c.values()),
            "components": gate_c,
            "n_positive_genes": n_positive,
            "sign_test_p": sign_p,
            "maximum_absolute_contribution": max_contribution,
        },
        "gate_D": {
            "pass": all(gate_d.values()),
            "components": gate_d,
            "adjusted_ratio": adjusted_ratio,
            "maximum_vif": max_vif,
            "freedman_lane_p": fl_p,
            "baseline_interaction_maxT_p": float(baseline_interaction["maxT_fwer"]),
            "baseline_iss_hedges_g": baseline_g.estimate,
        },
        "gate_F": {
            "pass": all(gate_f.values()),
            "components": gate_f,
            "adjusted_ratio": generic_ratio,
            "matched_percentile": generic_percentile,
        },
    }
    (outdir / "internal_gate_status.json").write_text(
        json.dumps(gate_status, indent=2, default=json_default) + "\n"
    )
    manifest = {
        "analysis_id": cfg["analysis_id"],
        "config": str(config_path),
        "observed_hashes": observed_hashes,
        "n_grey60_genes": len(genes),
        "n_samples_flt_gc": len(meta),
        "score_orientation": {
            "mean_z": {"multiplier": primary_dir, "correlation_with_ME17": primary_corr},
            "median_z": {"multiplier": median_dir, "correlation_with_ME17": median_corr},
            "weighted": {"multiplier": weighted_dir, "correlation_with_ME17": weighted_corr},
            "raw_vst": {"multiplier": raw_dir, "correlation_with_Rtech_score": raw_corr},
        },
        "selection_family_tests": 220,
        "maxT_permutations": n_perm,
        "bootstrap_replicates": n_boot,
        "matched_subset_draws": int(cfg["resampling"]["matched_subsets"]),
        "claim_boundary": (
            "Bulk-kidney ISS-terminal-associated score shift; no causal, "
            "cell-localized, kidney-specific, or topology-rewiring inference."
        ),
    }
    (outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, default=json_default) + "\n"
    )
    print(json.dumps(gate_status, indent=2, default=json_default))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default="config/grey60_adversarial_reanalysis.yaml"
    )
    parser.add_argument("--outdir", default="")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
