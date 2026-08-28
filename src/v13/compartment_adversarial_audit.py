"""Post-hoc adversarial closure audit for OSD-462 kidney compartments.

This module deliberately separates reference-only set construction from
effect-aware post-processing.  ``prepare`` may be run before the exact label
permutation because it reads only the external kidney atlas and the frozen
KEGG structural-control source.  ``postprocess`` consumes the emitted exact
run and writes the artifact, contributor, observability, Grey60, clinical-axis,
and decision summaries.

The output object is parent-protein annotation enrichment in whole kidney.
Nothing here identifies a phosphosite's cell of origin or repairs the perfect
condition-to-reporter-block alias in OSD-462.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO / "config/v13_compartment_adversarial_audit.yaml"
DEFAULT_EXACT = REPO / "data/results/run_20260802_v13_compartment_exact"
DEFAULT_OUT = REPO / "data/results/run_20260802_v13_compartment_adversarial"

PRIMARY_SCORE = "median_negative_site_effect"
SCORE_NAMES = (
    "median_negative_site_effect",
    "mean_negative_site_effect",
    "one_sided_maxmean",
)
PRIMARY_TIER = "all_enriched"


CLINICAL_AXIS_PANELS: dict[str, tuple[str, ...]] = {
    "glomerular_barrier_podocyte": (
        "Nphs1", "Nphs2", "Podxl", "Wt1", "Synpo", "Magi2", "Actn4",
        "Cd2ap", "Plce1",
    ),
    "proximal_injury_stress": (
        "Havcr1", "Lcn2", "Fabp1", "B2m", "Cst3", "Spp1", "Lrp2", "Cubn",
    ),
    "TIMP2_IGFBP7_stress": ("Timp2", "Igfbp7"),
    "inflammatory_injury": ("Il18", "Ccl2", "Ccl14", "Lcn2", "Tnf", "Il6", "Nfkb1"),
    "fibrosis_remodeling": (
        "Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2", "Col1a1", "Col1a2",
        "Col3a1", "Fn1", "Postn", "Mmp7", "Dkk3", "Egf", "Col6a3",
        "Pdgfra", "Pdgfrb", "Acta2", "Dcn", "Lum",
    ),
    "endothelial_vascular": (
        "Pecam1", "Cdh5", "Kdr", "Flt1", "Emcn", "Vwf", "Klf2", "Nos3",
        "Tek", "Robo4",
    ),
    "TAL_transport": (
        "Slc12a1", "Kcnj1", "Umod", "Cldn16", "Cldn19", "Clcnka", "Clcnkb",
        "Bsnd", "Slc9a3",
    ),
    "DCT_ASDN_transport": (
        "Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Nedd4l", "Sgk1",
        "Nr3c2", "Hsd11b2", "Scnn1a", "Scnn1b", "Scnn1g", "Aqp2", "Trpv5",
        "Calb1", "Slc8a1",
    ),
}


def sha256(path: str | Path) -> str:
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_config(path: str | Path) -> dict[str, Any]:
    with Path(path).open() as handle:
        cfg = yaml.safe_load(handle)
    if not isinstance(cfg, dict):
        raise ValueError(f"Configuration is not a mapping: {path}")
    return cfg


def resolve(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else REPO / value


def bh_fdr(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(p.shape, np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    if not len(finite):
        return out
    order = finite[np.argsort(p[finite], kind="stable")]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out[order] = np.minimum(ranked, 1.0)
    return out


def competitive_statistic(values: np.ndarray, membership: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    membership = np.asarray(membership, dtype=bool)
    finite = np.isfinite(values)
    inside = finite & membership
    outside = finite & ~membership
    if not inside.any() or not outside.any():
        return float("nan")
    return float(values[inside].mean() - values[outside].mean())


def _case_map(symbols: Iterable[str]) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for symbol in map(str, symbols):
        mapping.setdefault(symbol.upper(), symbol)
    return mapping


def build_marker_tiers(
    atlas: pd.DataFrame,
    structural_terms: Mapping[str, Sequence[str]],
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    """Build all reference-only compartment tiers and the scaffold control."""

    required = {
        "gene_symbol", "compartment", "mean_cpm",
        "source_study_detection_fraction",
    }
    missing = required - set(atlas)
    if missing:
        raise ValueError(f"Atlas expression table lacks {sorted(missing)}")

    marker_cfg = cfg["marker_tiers"]
    mapping = cfg["compartment_mapping"]
    expression = atlas.pivot(
        index="gene_symbol", columns="compartment", values="mean_cpm"
    ).sort_index()
    detection = atlas.pivot(
        index="gene_symbol",
        columns="compartment",
        values="source_study_detection_fraction",
    ).reindex_like(expression)
    missing_compartments = set(mapping) - set(expression.columns)
    if missing_compartments:
        raise ValueError(
            f"Configured atlas compartments are absent: {sorted(missing_compartments)}"
        )

    atlas_case = _case_map(expression.index)
    source_terms = cfg["structural_scaffold_control"]["source_terms"]
    missing_terms = set(source_terms) - set(structural_terms)
    if missing_terms:
        raise ValueError(f"Structural control terms are absent: {sorted(missing_terms)}")
    scaffold: set[str] = set()
    for term in source_terms:
        for raw in structural_terms[term]:
            symbol = atlas_case.get(str(raw).upper())
            if symbol:
                scaffold.add(symbol)

    min_cpm = float(marker_cfg["target_mean_cpm_min"])
    min_detection = float(
        marker_cfg["target_source_study_detection_fraction_min"]
    )
    min_ratio = float(marker_cfg["all_enriched_log2_target_to_max_other_min"])
    high_ratio = float(
        marker_cfg["high_specificity_log2_target_to_max_other_min"]
    )
    high_detection = float(
        marker_cfg["high_specificity_target_detection_fraction_min"]
    )
    broad_min = int(marker_cfg["within_kidney_broad_non_target_count_min"])
    broad_cpm = float(
        marker_cfg["within_kidney_broad_non_target_mean_cpm_min"]
    )
    broad_detection = float(
        marker_cfg["within_kidney_broad_non_target_detection_fraction_min"]
    )

    rows: list[dict[str, Any]] = []
    for atlas_name, report_name in mapping.items():
        other = [name for name in expression.columns if name != atlas_name]
        target = expression[atlas_name]
        target_detection = detection[atlas_name]
        max_other = expression[other].max(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = np.log2(target / max_other)
        ratio = ratio.where(max_other > 0, np.inf)
        non_target_detected = (
            (expression[other] >= broad_cpm)
            & (detection[other] >= broad_detection)
        ).sum(axis=1)
        base = (
            (target >= min_cpm)
            & (target_detection >= min_detection)
            & (ratio >= min_ratio)
        )
        tier_masks = {
            "all_enriched": base,
            "within_kidney_not_broad": base & (non_target_detected < broad_min),
            "high_specificity": (
                base & (ratio >= high_ratio) & (target_detection >= high_detection)
            ),
            "broad_enriched": base & (non_target_detected >= broad_min),
            "scaffold_excluded": base & ~expression.index.isin(scaffold),
        }
        for tier, mask in tier_masks.items():
            for gene in expression.index[mask]:
                rows.append(
                    {
                        "gene_symbol": gene,
                        "gene_set": f"{report_name}__{tier}",
                        "atlas_compartment": atlas_name,
                        "report_compartment": report_name,
                        "tier": tier,
                        "category": "kidney_compartment_reference",
                        "definition_source": "mouse_kidney_atlas_flight_blind",
                        "status": "defined",
                        "final_for_testing": True,
                        "target_cpm": float(target.loc[gene]),
                        "max_other_cpm": float(max_other.loc[gene]),
                        "log2_target_to_max_other": float(ratio.loc[gene]),
                        "target_source_study_fraction": float(
                            target_detection.loc[gene]
                        ),
                        "n_non_target_compartments_detected": int(
                            non_target_detected.loc[gene]
                        ),
                        "structural_scaffold_member": gene in scaffold,
                    }
                )

    control_name = cfg["structural_scaffold_control"]["gene_set_name"]
    for gene in sorted(scaffold):
        rows.append(
            {
                "gene_symbol": gene,
                "gene_set": control_name,
                "atlas_compartment": "not_applicable",
                "report_compartment": "structural_scaffold_control",
                "tier": "all",
                "category": "protein_class_negative_control",
                "definition_source": "KEGG_2019_Mouse_frozen_union",
                "status": "defined",
                "final_for_testing": True,
                "target_cpm": np.nan,
                "max_other_cpm": np.nan,
                "log2_target_to_max_other": np.nan,
                "target_source_study_fraction": np.nan,
                "n_non_target_compartments_detected": np.nan,
                "structural_scaffold_member": True,
            }
        )

    table = pd.DataFrame(rows).drop_duplicates(["gene_set", "gene_symbol"])
    expected = list(map(str, cfg["set_test"]["primary_family"]))
    observed = set(table["gene_set"])
    missing_sets = set(expected) - observed
    unexpected = observed - set(expected)
    if missing_sets or unexpected:
        raise ValueError(
            "Reference-only set construction disagrees with frozen family: "
            f"missing={sorted(missing_sets)}, unexpected={sorted(unexpected)}"
        )
    order = {name: index for index, name in enumerate(expected)}
    table["_order"] = table["gene_set"].map(order)
    table = table.sort_values(["_order", "gene_symbol"], kind="stable").drop(
        columns="_order"
    )
    return table.reset_index(drop=True)


def prepare_reference(config_path: str | Path) -> Path:
    cfg = load_config(config_path)
    atlas_path = resolve(cfg["input"]["atlas_expression"])
    structural_path = resolve(cfg["input"]["structural_scaffold_reference"])
    atlas = pd.read_csv(atlas_path, sep="\t")
    structural = json.loads(structural_path.read_text())
    table = build_marker_tiers(atlas, structural, cfg)
    output = resolve(cfg["input"]["gene_sets"])
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)
    counts = (
        table.groupby(
            ["report_compartment", "tier", "gene_set"], dropna=False
        ).size().rename("n_reference_genes").reset_index()
    )
    counts.to_csv(output.with_name("frozen_compartment_tier_counts.tsv"), sep="\t", index=False)
    manifest = {
        "analysis_id": cfg["analysis"]["id"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "reference_only_before_new_exact_effect_inference",
        "config": str(Path(config_path)),
        "config_sha256": sha256(config_path),
        "inputs": {
            str(atlas_path.relative_to(REPO)): sha256(atlas_path),
            str(structural_path.relative_to(REPO)): sha256(structural_path),
        },
        "output": {
            str(output.relative_to(REPO)): sha256(output),
            str(output.with_name("frozen_compartment_tier_counts.tsv").relative_to(REPO)): sha256(
                output.with_name("frozen_compartment_tier_counts.tsv")
            ),
        },
        "n_sets": int(table["gene_set"].nunique()),
        "n_rows": int(len(table)),
        "terminology_boundary": cfg["marker_tiers"]["terminology_boundary"],
    }
    output.with_name("reference_freeze_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    return output


def load_primary_gene_covariates(cfg: Mapping[str, Any]) -> pd.DataFrame:
    source = resolve(cfg["input"]["source_exact_run"])
    calibration = pd.read_csv(
        source / "parent_gene_null_calibration.tsv", sep="\t"
    )
    genes = calibration[
        calibration["profile"].eq("primary")
        & calibration["gene_score"].eq(PRIMARY_SCORE)
        & calibration["fixed_null_universe_eligible"].astype(bool)
    ].copy()
    genes = genes.drop_duplicates("gene_symbol").reset_index(drop=True)

    sites = pd.read_csv(source / "site_universe_primary.tsv", sep="\t")
    sites["phosphopeptide_count"] = sites[
        ["n_pep_Samp1-5", "n_pep_Samp6-10"]
    ].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    peptide = (
        sites.groupby("gene_symbol")["phosphopeptide_count"]
        .median()
        .rename("median_phosphopeptide_count")
    )
    genes = genes.merge(peptide, left_on="gene_symbol", right_index=True, how="left")

    from src.v13.continuous_phospho_inference import PLEX1, PLEX2, parse_tmt_sheet

    protein = parse_tmt_sheet(
        resolve(cfg["input"]["protein_workbook"]),
        "protein_quant_2721",
        gene_col="Gene Symbol",
        peptide_cols={PLEX1: f"{PLEX1} Peptides", PLEX2: f"{PLEX2} Peptides"},
        id_col="Protein Id",
    )
    raw = protein.scaled.to_numpy(dtype=float)
    raw = np.where(raw > 0, raw, np.nan)
    with np.errstate(invalid="ignore"):
        row_abundance = np.nanmedian(np.log2(raw), axis=1)
    p1 = pd.to_numeric(protein.meta[f"n_pep_{PLEX1}"], errors="coerce")
    p2 = pd.to_numeric(protein.meta[f"n_pep_{PLEX2}"], errors="coerce")
    row_peptides = pd.concat([p1, p2], axis=1).mean(axis=1).fillna(0.0)
    row_frame = pd.DataFrame(
        {
            "gene_symbol": protein.meta["gene_symbol"].astype(str),
            "row_abundance": row_abundance,
            "row_peptides": row_peptides.to_numpy(dtype=float),
        }
    )
    parent_rows: list[dict[str, Any]] = []
    for symbol, part in row_frame.groupby("gene_symbol", sort=False):
        finite = np.isfinite(part["row_abundance"].to_numpy(dtype=float))
        weights = part["row_peptides"].to_numpy(dtype=float)
        weights = np.where(weights > 0, weights, 1.0)
        abundance = (
            float(np.average(part.loc[finite, "row_abundance"], weights=weights[finite]))
            if finite.any()
            else np.nan
        )
        parent_rows.append(
            {
                "gene_symbol": symbol,
                "parent_protein_log2_abundance": abundance,
                "parent_protein_peptide_count": float(
                    np.nansum(part["row_peptides"].to_numpy(dtype=float))
                ),
                "parent_protein_available": bool(finite.any()),
            }
        )
    parent = pd.DataFrame(parent_rows)
    genes = genes.merge(parent, on="gene_symbol", how="left")
    genes["parent_protein_available"] = genes["parent_protein_available"].fillna(False)
    genes["parent_protein_peptide_count"] = genes[
        "parent_protein_peptide_count"
    ].fillna(0.0)
    return genes


def _standardized_observability_matrix(frame: pd.DataFrame) -> tuple[np.ndarray, list[str]]:
    available = frame["parent_protein_available"].astype(bool).to_numpy()
    parent_abundance = pd.to_numeric(
        frame["parent_protein_log2_abundance"], errors="coerce"
    ).to_numpy(dtype=float)
    if np.isfinite(parent_abundance[available]).any():
        fill = float(np.nanmedian(parent_abundance[available]))
    else:
        fill = 0.0
    parent_abundance = np.where(np.isfinite(parent_abundance), parent_abundance, fill)
    columns = {
        "log1p_quantified_site_count": np.log1p(
            pd.to_numeric(frame["n_sites"], errors="coerce").fillna(0).to_numpy()
        ),
        "median_phosphopeptide_log2_signal": pd.to_numeric(
            frame["median_log2_signal"], errors="coerce"
        ).to_numpy(),
        "mean_missing_fraction": pd.to_numeric(
            frame["mean_missing_fraction"], errors="coerce"
        ).fillna(0).to_numpy(),
        "log1p_median_phosphopeptide_count": np.log1p(
            pd.to_numeric(
                frame["median_phosphopeptide_count"], errors="coerce"
            ).fillna(0).to_numpy()
        ),
        "parent_protein_available": available.astype(float),
        "parent_protein_log2_abundance": parent_abundance,
        "log1p_parent_protein_peptide_count": np.log1p(
            pd.to_numeric(
                frame["parent_protein_peptide_count"], errors="coerce"
            ).fillna(0).to_numpy()
        ),
    }
    names = list(columns)
    matrix = np.column_stack([columns[name] for name in names]).astype(float)
    for column in range(matrix.shape[1]):
        values = matrix[:, column]
        finite = np.isfinite(values)
        replacement = float(np.nanmedian(values[finite])) if finite.any() else 0.0
        values = np.where(finite, values, replacement)
        scale = float(np.std(values, ddof=0))
        matrix[:, column] = (values - float(np.mean(values))) / (scale if scale > 0 else 1.0)
    return matrix, names


def _propensity_strata(
    matrix: np.ndarray,
    membership: np.ndarray,
    parent_available: np.ndarray,
    *,
    n_bins: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.linear_model import LogisticRegression

    model = LogisticRegression(
        C=0.5,
        class_weight="balanced",
        max_iter=2000,
        random_state=seed,
        solver="lbfgs",
    )
    model.fit(matrix, membership.astype(int))
    propensity = model.predict_proba(matrix)[:, 1]
    strata = np.full(len(membership), -1, dtype=int)
    offset = 0
    for availability in (False, True):
        idx = np.flatnonzero(parent_available == availability)
        if not len(idx):
            continue
        order = idx[np.argsort(propensity[idx], kind="stable")]
        bins = min(n_bins, len(order))
        for bin_id, members in enumerate(np.array_split(order, bins)):
            strata[members] = offset + bin_id
        offset += bins
    if np.any(strata < 0):
        raise RuntimeError("Propensity strata did not cover the eligible universe")
    return propensity, strata


def _stratified_null_means(
    values: np.ndarray,
    strata: np.ndarray,
    target: np.ndarray,
    *,
    n_draws: int,
    rng: np.random.Generator,
    chunk_size: int = 250,
) -> np.ndarray:
    target_counts = pd.Series(strata[target]).value_counts().to_dict()
    null_sum = np.zeros(n_draws, dtype=float)
    for stratum, count_value in target_counts.items():
        count = int(count_value)
        pool = values[strata == int(stratum)]
        if count > len(pool):
            raise RuntimeError("Matched stratum requests more genes than available")
        if count == len(pool):
            null_sum += float(np.sum(pool))
            continue
        for start in range(0, n_draws, chunk_size):
            stop = min(n_draws, start + chunk_size)
            keys = rng.random((stop - start, len(pool)))
            take = np.argpartition(keys, count - 1, axis=1)[:, :count]
            null_sum[start:stop] += pool[take].sum(axis=1)
    return null_sum / int(target.sum())


def observability_matched_tests(
    genes: pd.DataFrame,
    membership: pd.DataFrame,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    matching = cfg["observability_matching"]
    set_names = list(map(str, matching["gene_sets"]))
    n_draws = int(matching["n_draws"])
    n_bins = int(matching["propensity_bins_per_parent_availability_class"])
    seed = int(matching["seed"])
    minimum = int(matching["minimum_observable_genes"])
    sets = {
        name: set(membership.loc[membership["gene_set"].eq(name), "gene_symbol"])
        for name in set_names
    }
    z = pd.to_numeric(genes["observed_gene_z"], errors="coerce").to_numpy(dtype=float)
    matrix, feature_names = _standardized_observability_matrix(genes)
    availability = genes["parent_protein_available"].astype(bool).to_numpy()
    rows: list[dict[str, Any]] = []
    for index, name in enumerate(set_names):
        target = genes["gene_symbol"].isin(sets[name]).to_numpy()
        row: dict[str, Any] = {
            "gene_set": name,
            "n_observable_genes": int(target.sum()),
            "n_draws": n_draws,
            "evaluation_status": "non_evaluable",
            "observed_mean_gene_z": np.nan,
            "matched_null_mean": np.nan,
            "matched_null_sd": np.nan,
            "matched_standardized_effect": np.nan,
            "matched_empirical_p_greater": np.nan,
            "max_abs_expected_matched_smd": np.nan,
            "matching_features": ";".join(feature_names),
        }
        if int(target.sum()) < minimum:
            rows.append(row)
            continue
        propensity, strata = _propensity_strata(
            matrix,
            target,
            availability,
            n_bins=n_bins,
            seed=seed + index,
        )
        rng = np.random.default_rng(seed + 1009 * (index + 1))
        null = _stratified_null_means(
            z, strata, target, n_draws=n_draws, rng=rng
        )
        observed = float(np.mean(z[target]))
        null_mean = float(np.mean(null))
        null_sd = float(np.std(null, ddof=0))
        expected = np.zeros(matrix.shape[1], dtype=float)
        target_counts = pd.Series(strata[target]).value_counts().to_dict()
        for stratum, count in target_counts.items():
            expected += int(count) * matrix[strata == int(stratum)].mean(axis=0)
        expected /= int(target.sum())
        target_mean = matrix[target].mean(axis=0)
        max_smd = float(np.max(np.abs(target_mean - expected)))
        p_value = float((1 + np.count_nonzero(null >= observed)) / (n_draws + 1))
        row.update(
            {
                "evaluation_status": "evaluated",
                "observed_mean_gene_z": observed,
                "matched_null_mean": null_mean,
                "matched_null_sd": null_sd,
                "matched_standardized_effect": (
                    (observed - null_mean) / null_sd if null_sd > 0 else np.nan
                ),
                "matched_empirical_p_greater": p_value,
                "max_abs_expected_matched_smd": max_smd,
                "propensity_min": float(propensity.min()),
                "propensity_max": float(propensity.max()),
                "n_matching_strata": int(pd.Series(strata).nunique()),
            }
        )
        rows.append(row)
    result = pd.DataFrame(rows)
    result["matched_bh_q"] = bh_fdr(result["matched_empirical_p_greater"])
    return result


def _set_metadata(membership: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "gene_set", "report_compartment", "tier", "category",
    ]
    return membership[columns].drop_duplicates("gene_set")


def compile_exact_family(
    exact_dir: Path,
    membership: pd.DataFrame,
    genes: pd.DataFrame,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    results = pd.read_csv(exact_dir / "set_level_permutation_inference.tsv", sep="\t")
    results = results[
        results["profile"].eq("primary")
        & results["exclusion"].eq("none")
        & results["set_test"].eq("competitive")
        & results["gene_score"].isin(SCORE_NAMES)
    ].copy()
    counts = membership.groupby("gene_set")["gene_symbol"].nunique().rename(
        "n_reference_genes"
    )
    observable = {
        name: int(
            genes["gene_symbol"].isin(
                set(part["gene_symbol"].astype(str))
            ).sum()
        )
        for name, part in membership.groupby("gene_set", sort=False)
    }
    family = pd.MultiIndex.from_product(
        [cfg["set_test"]["primary_family"], SCORE_NAMES],
        names=["gene_set", "gene_score"],
    ).to_frame(index=False)
    family = family.merge(results, on=["gene_set", "gene_score"], how="left")
    family = family.merge(counts, on="gene_set", how="left")
    family = family.merge(_set_metadata(membership), on="gene_set", how="left")
    family["n_observable_genes"] = family["n_observable_genes"].fillna(
        family["gene_set"].map(observable)
    )
    family["evaluation_status"] = np.where(
        family["observed_statistic"].notna(), "evaluated", "non_evaluable"
    )
    return family


def profile_direction_summary(
    cfg: Mapping[str, Any], membership: pd.DataFrame
) -> pd.DataFrame:
    source = resolve(cfg["input"]["source_exact_run"])
    calibration = pd.read_csv(source / "parent_gene_null_calibration.tsv", sep="\t")
    selected_sets = [
        f"{name}__{PRIMARY_TIER}"
        for name in cfg["compartment_mapping"].values()
    ] + [cfg["structural_scaffold_control"]["gene_set_name"]]
    sets = {
        name: set(membership.loc[membership["gene_set"].eq(name), "gene_symbol"])
        for name in selected_sets
    }
    rows: list[dict[str, Any]] = []
    for (profile, score), part in calibration.groupby(
        ["profile", "gene_score"], sort=False
    ):
        if score not in SCORE_NAMES:
            continue
        part = part[part["fixed_null_universe_eligible"].astype(bool)].copy()
        values = pd.to_numeric(part["observed_gene_z"], errors="coerce").to_numpy()
        for set_name in selected_sets:
            mask = part["gene_symbol"].isin(sets[set_name]).to_numpy()
            rows.append(
                {
                    "profile": profile,
                    "gene_score": score,
                    "gene_set": set_name,
                    "n_observable_genes": int(mask.sum()),
                    "observed_competitive_statistic": competitive_statistic(
                        values, mask
                    ),
                }
            )
    return pd.DataFrame(rows)


def leave_one_gene_out_tables(
    genes: pd.DataFrame, membership: pd.DataFrame, cfg: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    z = pd.to_numeric(genes["observed_gene_z"], errors="coerce").to_numpy(dtype=float)
    rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for set_name in cfg["set_test"]["primary_family"]:
        members = set(
            membership.loc[membership["gene_set"].eq(set_name), "gene_symbol"]
        )
        mask = genes["gene_symbol"].isin(members).to_numpy()
        idx = np.flatnonzero(mask)
        if len(idx) < 2:
            continue
        background_mean = float(np.mean(z[~mask]))
        observed = float(np.mean(z[mask]) - background_mean)
        loo_values: list[float] = []
        for gene_index in idx:
            keep = idx[idx != gene_index]
            statistic = float(np.mean(z[keep]) - background_mean)
            loo_values.append(statistic)
            rows.append(
                {
                    "gene_set": set_name,
                    "removed_gene": str(genes.iloc[gene_index]["gene_symbol"]),
                    "removed_gene_z": float(z[gene_index]),
                    "observed_full_statistic": observed,
                    "leave_one_gene_out_statistic": statistic,
                }
            )
        positive = np.sort(z[mask & (z > 0)])[::-1]
        total_positive = float(positive.sum()) if len(positive) else np.nan
        top1 = float(positive[:1].sum() / total_positive) if total_positive > 0 else np.nan
        top10 = float(positive[:10].sum() / total_positive) if total_positive > 0 else np.nan
        summary_rows.append(
            {
                "gene_set": set_name,
                "n_observable_genes": int(mask.sum()),
                "observed_full_statistic": observed,
                "minimum_leave_one_gene_out_statistic": float(np.min(loo_values)),
                "maximum_leave_one_gene_out_statistic": float(np.max(loo_values)),
                "all_leave_one_gene_out_positive": bool(np.min(loo_values) > 0),
                "n_positive_gene_z": int(len(positive)),
                "top_gene_fraction_of_positive_z": top1,
                "top_ten_fraction_of_positive_z": top10,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(summary_rows)


def leading_contributors(
    genes: pd.DataFrame, membership: pd.DataFrame, cfg: Mapping[str, Any]
) -> pd.DataFrame:
    metadata = membership[
        membership["tier"].eq(PRIMARY_TIER)
    ][
        [
            "gene_set", "gene_symbol", "report_compartment", "target_cpm",
            "max_other_cpm", "log2_target_to_max_other",
            "target_source_study_fraction", "n_non_target_compartments_detected",
            "structural_scaffold_member",
        ]
    ].copy()
    out = metadata.merge(genes, on="gene_symbol", how="inner")
    grey_path = resolve(cfg["input"]["grey60_gene_contrasts"])
    grey = set()
    if grey_path.exists():
        grey_table = pd.read_csv(grey_path, sep="\t")
        if "symbol" in grey_table:
            grey = set(grey_table["symbol"].dropna().astype(str))
    out["grey60_member"] = out["gene_symbol"].isin(grey)
    out["positive_suppression_z"] = out["observed_gene_z"].clip(lower=0)
    out = out.sort_values(
        ["report_compartment", "observed_gene_z"],
        ascending=[True, False],
        kind="stable",
    )
    out["contributor_rank"] = out.groupby("report_compartment").cumcount() + 1
    total = out.groupby("report_compartment")["positive_suppression_z"].transform("sum")
    out["fraction_of_compartment_positive_z"] = np.where(
        total > 0, out["positive_suppression_z"] / total, np.nan
    )
    return out.reset_index(drop=True)


def grey60_overlap(
    membership: pd.DataFrame,
    genes: pd.DataFrame,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    from scipy.stats import fisher_exact

    grey_table = pd.read_csv(resolve(cfg["input"]["grey60_gene_contrasts"]), sep="\t")
    grey = set(grey_table["symbol"].dropna().astype(str))
    atlas = pd.read_csv(resolve(cfg["input"]["atlas_expression"]), sep="\t")
    universe = set(atlas["gene_symbol"].astype(str))
    grey &= universe
    observable = set(genes["gene_symbol"].astype(str))
    rows: list[dict[str, Any]] = []
    primary = membership[membership["tier"].eq(PRIMARY_TIER)]
    for set_name, part in primary.groupby("gene_set", sort=False):
        members = set(part["gene_symbol"].astype(str)) & universe
        a = len(members & grey)
        b = len(members - grey)
        c = len(grey - members)
        d = len(universe - members - grey)
        odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "gene_set": set_name,
                "n_atlas_genes": len(members),
                "n_grey60_atlas_genes": len(grey),
                "n_overlap": a,
                "fisher_odds_ratio": float(odds),
                "fisher_p_greater": float(p_value),
                "n_overlap_phospho_observable": len(members & grey & observable),
                "observable_overlap_genes": ";".join(sorted(members & grey & observable)),
            }
        )
    return pd.DataFrame(rows)


def clinical_axis_observability(
    cfg: Mapping[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rna_files = {
        "UPX": REPO / "data/external/osdr/OSD-462/GLDS-462_rna_seq_differential_expression_UPX_GLbulkRNAseq.csv",
        "mRNA": REPO / "data/external/osdr/OSD-462/GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv",
        "totalRNA": REPO / "data/external/osdr/OSD-462/GLDS-462_rna_seq_differential_expression_totRNA_GLbulkRNAseq.csv",
    }
    rna: dict[str, pd.DataFrame] = {}
    for name, path in rna_files.items():
        table = pd.read_csv(path).dropna(subset=["SYMBOL"]).drop_duplicates("SYMBOL")
        rna[name] = table.set_index("SYMBOL")
    protein_path = REPO / "data/results/codex_review_osd462/osd462_anchor/protein_effects_gene_anyplex.tsv"
    protein = pd.read_csv(protein_path, sep="\t").drop_duplicates("gene_symbol").set_index("gene_symbol")
    calibration = pd.read_csv(
        resolve(cfg["input"]["source_exact_run"]) / "parent_gene_null_calibration.tsv",
        sep="\t",
    )
    phospho = calibration[
        calibration["profile"].eq("primary")
        & calibration["gene_score"].eq(PRIMARY_SCORE)
        & calibration["fixed_null_universe_eligible"].astype(bool)
    ].drop_duplicates("gene_symbol").set_index("gene_symbol")

    rows: list[dict[str, Any]] = []
    effect_column = "Log2fc_(Space Flight)v(Ground Control)"
    p_column = "P.value_(Space Flight)v(Ground Control)"
    q_column = "Adj.p.value_(Space Flight)v(Ground Control)"
    for axis, panel in CLINICAL_AXIS_PANELS.items():
        for gene in panel:
            row: dict[str, Any] = {"molecular_axis": axis, "gene_symbol": gene}
            for preparation, table in rna.items():
                row[f"rna_{preparation}_log2fc"] = (
                    float(table.loc[gene, effect_column])
                    if gene in table.index and pd.notna(table.loc[gene, effect_column])
                    else np.nan
                )
                row[f"rna_{preparation}_p"] = (
                    float(table.loc[gene, p_column])
                    if gene in table.index and p_column in table and pd.notna(table.loc[gene, p_column])
                    else np.nan
                )
                row[f"rna_{preparation}_q"] = (
                    float(table.loc[gene, q_column])
                    if gene in table.index and q_column in table and pd.notna(table.loc[gene, q_column])
                    else np.nan
                )
            row["protein_flight_effect"] = (
                float(protein.loc[gene, "flight_effect"])
                if gene in protein.index and pd.notna(protein.loc[gene, "flight_effect"])
                else np.nan
            )
            row["protein_abundance_log2"] = (
                float(protein.loc[gene, "abundance_log2"])
                if gene in protein.index and pd.notna(protein.loc[gene, "abundance_log2"])
                else np.nan
            )
            row["phospho_gene_z"] = (
                float(phospho.loc[gene, "observed_gene_z"])
                if gene in phospho.index and pd.notna(phospho.loc[gene, "observed_gene_z"])
                else np.nan
            )
            row["phospho_n_sites"] = (
                int(phospho.loc[gene, "n_sites"])
                if gene in phospho.index and pd.notna(phospho.loc[gene, "n_sites"])
                else 0
            )
            rows.append(row)
    detail = pd.DataFrame(rows)
    summaries: list[dict[str, Any]] = []
    for axis, part in detail.groupby("molecular_axis", sort=False):
        summary: dict[str, Any] = {
            "molecular_axis": axis,
            "n_panel_genes": len(part),
            "n_protein_observable": int(part["protein_flight_effect"].notna().sum()),
            "n_phospho_parent_observable": int(part["phospho_gene_z"].notna().sum()),
            "mean_protein_flight_effect": float(part["protein_flight_effect"].mean()),
            "mean_phospho_gene_z": float(part["phospho_gene_z"].mean()),
        }
        for preparation in rna:
            column = f"rna_{preparation}_log2fc"
            summary[f"n_rna_{preparation}_observable"] = int(part[column].notna().sum())
            summary[f"mean_rna_{preparation}_log2fc"] = float(part[column].mean())
            summary[f"n_rna_{preparation}_q_lt_0_05"] = int(
                (part[f"rna_{preparation}_q"] < 0.05).sum()
            )
        summaries.append(summary)
    return detail, pd.DataFrame(summaries)


def reporter_position_by_layer(cfg: Mapping[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = resolve(cfg["input"]["source_layer_block_diagnostic"])
    blocks = (("baseline", 0, 5), ("flight", 5, 10), ("ground", 10, 15))
    permutations = list(itertools.permutations(range(5)))
    rows: list[dict[str, Any]] = []
    for layer in ("phospho", "protein"):
        cache = np.load(source / f"_{layer}_cache.npz", allow_pickle=True)
        slopes: list[float] = []
        for plex in (1, 2):
            raw = cache[f"plex{plex}"]
            ok = np.all(np.isfinite(raw) & (raw > 0), axis=1)
            values = np.log2(raw[ok])
            values -= values.mean(axis=1, keepdims=True)
            for block_name, start, stop in blocks:
                block = values[:, start:stop]
                pos = np.arange(1, 6, dtype=float)
                pos -= pos.mean()
                denominator = float(np.sum(pos ** 2))
                observed = float(np.mean(np.sum(block * pos, axis=1) / denominator))
                null = np.asarray(
                    [
                        float(
                            np.mean(
                                np.sum(block[:, perm] * pos, axis=1) / denominator
                            )
                        )
                        for perm in permutations
                    ]
                )
                p_value = float(np.mean(np.abs(null) >= abs(observed)))
                slopes.append(observed)
                rows.append(
                    {
                        "layer": layer,
                        "plex": plex,
                        "block": block_name,
                        "n_complete_features": int(len(block)),
                        "slope_log2_per_channel_step": observed,
                        "exact_position_permutation_p_two_sided": p_value,
                        "n_exact_position_permutations": len(permutations),
                    }
                )
        pooled = float(np.mean(slopes))
        for row in rows:
            if row["layer"] == layer:
                row["layer_pooled_slope"] = pooled
                row["predicted_flight_minus_ground_from_linear_position"] = -5 * pooled
    marginal = pd.read_csv(source / "marginal_block_profile_by_layer.tsv", sep="\t")
    marginal["flight_is_block_minimum"] = (
        (marginal["flight_block_mean_log2"] < marginal["baseline_block_mean_log2"])
        & (marginal["flight_block_mean_log2"] < marginal["ground_block_mean_log2"])
    )
    marginal["baseline_flight_ground_order"] = marginal.apply(
        lambda row: " < ".join(
            sorted(
                ("baseline", "flight", "ground"),
                key=lambda name: row[f"{name}_block_mean_log2"],
            )
        ),
        axis=1,
    )
    return pd.DataFrame(rows), marginal


def decision_table(
    exact: pd.DataFrame,
    matched: pd.DataFrame,
    profiles: pd.DataFrame,
    loo_summary: pd.DataFrame,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    gate = cfg["post_hoc_candidate_gate"]
    alpha = float(gate["alpha"])
    rows: list[dict[str, Any]] = []
    for compartment in cfg["compartment_mapping"].values():
        base_name = f"{compartment}__all_enriched"
        high_name = f"{compartment}__high_specificity"
        base = exact[exact["gene_set"].eq(base_name)]
        primary = base[base["gene_score"].eq(PRIMARY_SCORE)]
        high = exact[
            exact["gene_set"].eq(high_name)
            & exact["gene_score"].eq(PRIMARY_SCORE)
        ]
        three = base.set_index("gene_score")["observed_statistic"]
        match_row = matched[matched["gene_set"].eq(base_name)]
        loo = loo_summary[loo_summary["gene_set"].eq(base_name)]
        profile = profiles[
            profiles["gene_set"].eq(base_name)
            & profiles["gene_score"].eq(PRIMARY_SCORE)
        ].set_index("profile")["observed_competitive_statistic"]
        primary_ok = bool(
            len(primary)
            and primary.iloc[0]["evaluation_status"] == "evaluated"
            and primary.iloc[0]["observed_statistic"] > 0
            and primary.iloc[0]["max_absT_fwer"] <= alpha
        )
        primary_evaluable = bool(
            len(primary) and primary.iloc[0]["evaluation_status"] == "evaluated"
        )
        three_ok = bool(
            primary_evaluable
            and all(
                name in three.index and pd.notna(three[name]) and three[name] > 0
                for name in SCORE_NAMES
            )
        )
        loo_ok = bool(
            primary_evaluable
            and len(loo)
            and loo.iloc[0]["all_leave_one_gene_out_positive"]
        )
        high_ok = bool(
            len(high)
            and high.iloc[0]["evaluation_status"] == "evaluated"
            and high.iloc[0]["observed_statistic"] > 0
        )
        matched_ok = bool(
            len(match_row)
            and match_row.iloc[0]["evaluation_status"] == "evaluated"
            and match_row.iloc[0]["matched_bh_q"] <= alpha
        )
        profile_ok = bool(
            primary_evaluable
            and "official_scaled_uncentered" in profile.index
            and "parent_protein_subtracted" in profile.index
            and profile["official_scaled_uncentered"] > 0
            and profile["parent_protein_subtracted"] > 0
        )
        concentration_ok = bool(
            primary_evaluable
            and len(loo)
            and loo.iloc[0]["top_gene_fraction_of_positive_z"]
            <= float(gate["maximum_top_gene_fraction_of_positive_z"])
            and loo.iloc[0]["top_ten_fraction_of_positive_z"]
            <= float(gate["maximum_top_ten_fraction_of_positive_z"])
        )
        numeric_pass = all(
            (primary_ok, three_ok, loo_ok, high_ok, matched_ok, profile_ok, concentration_ok)
        )
        rows.append(
            {
                "compartment": compartment,
                "n_primary_observable_genes": (
                    int(primary.iloc[0]["n_observable_genes"])
                    if len(primary)
                    else 0
                ),
                "primary_two_sided_family_gate": primary_ok,
                "three_gene_scores_direction_gate": three_ok,
                "leave_one_gene_out_gate": loo_ok,
                "high_specificity_gate": high_ok,
                "observability_matched_gate": matched_ok,
                "uncentered_parent_adjusted_direction_gate": profile_ok,
                "contributor_concentration_gate": concentration_ok,
                "adversarial_numeric_gate_pass": numeric_pass,
                "publication_claim_permitted": False,
                "publication_blocker": (
                    "condition is perfectly aliased with reporter-tag blocks; "
                    "whole-kidney parent proteins do not identify cell of origin"
                ),
                "decision": (
                    "NUMERIC_PATTERN_ONLY_DESIGN_BLOCKED"
                    if numeric_pass
                    else "FAIL_ADVERSARIAL_GATE"
                ),
            }
        )
    return pd.DataFrame(rows)


def plot_primary_family(exact: pd.DataFrame, output: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    primary = exact[
        exact["gene_score"].eq(PRIMARY_SCORE)
        & (
            exact["tier"].eq(PRIMARY_TIER)
            | exact["category"].eq("protein_class_negative_control")
        )
    ].copy()
    primary = primary[primary["evaluation_status"].eq("evaluated")]
    primary = primary.sort_values("observed_statistic")
    labels = []
    for value in primary["gene_set"]:
        if value == "broad_structural_scaffold_control__all":
            labels.append("structural/scaffold\ncontrol")
        else:
            labels.append(value.replace("__all_enriched", "").replace("_", " "))
    y = np.arange(len(primary))
    colors = np.where(primary["max_absT_fwer"] <= 0.05, "#9B2226", "#457B9D")
    fig, ax = plt.subplots(figsize=(8.0, max(4.5, 0.46 * len(primary) + 1.5)))
    ax.hlines(
        y,
        primary["null_ci_low"],
        primary["null_ci_high"],
        color="#A8ADB3",
        linewidth=2.0,
        label="central 95% exact-label null interval",
    )
    ax.scatter(primary["observed_statistic"], y, c=colors, s=48, zorder=3)
    ax.axvline(0, color="#333333", linewidth=1, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Competitive standardized parent-gene suppression statistic")
    ax.set_title("One frozen kidney-compartment family\n(two-sided max-|T| across 51 sets)")
    ax.legend(frameon=False, loc="lower right", fontsize=8)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    paths = [output / "primary_compartment_family.png", output / "primary_compartment_family.pdf"]
    fig.savefig(paths[0], dpi=220, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def write_decision_report(
    output: Path,
    decisions: pd.DataFrame,
    exact: pd.DataFrame,
    matched: pd.DataFrame,
    loo: pd.DataFrame,
    reporter: pd.DataFrame,
    cfg: Mapping[str, Any],
) -> Path:
    primary = exact[
        exact["gene_score"].eq(PRIMARY_SCORE)
        & exact["tier"].eq(PRIMARY_TIER)
    ][
        [
            "report_compartment", "n_observable_genes", "observed_statistic",
            "empirical_p_greater", "empirical_p_less", "max_absT_fwer",
        ]
    ].sort_values("observed_statistic", ascending=False)
    scaffold_name = cfg["structural_scaffold_control"]["gene_set_name"]
    scaffold = exact[
        exact["gene_set"].eq(scaffold_name)
        & exact["gene_score"].eq(PRIMARY_SCORE)
    ]
    significant_tiers = exact[
        exact["max_absT_fwer"].le(0.05)
    ][
        [
            "gene_set", "gene_score", "n_observable_genes",
            "observed_statistic", "empirical_p_greater", "empirical_p_less",
            "max_absT_fwer",
        ]
    ].sort_values(["max_absT_fwer", "gene_set"])
    lines = [
        "# V13 compartment adversarial audit: decision report",
        "",
        "This is a post-hoc closure audit. It does not repair the perfect "
        "condition-to-reporter-block alias and cannot localize parent proteins "
        "to a kidney cell type.",
        "",
        "## Frozen primary family",
        "",
        primary.to_markdown(index=False, floatfmt=".5g"),
        "",
        "## Adversarial gates",
        "",
        decisions.to_markdown(index=False),
        "",
        "No compartment passed every frozen adversarial gate. A tier-specific "
        "row may still pass the 51-set family; such a row is reported below and "
        "is not promoted when its all-enriched, normalization, observability, or "
        "contributor gates fail.",
        "",
        "## Tier/score rows passing their score-specific two-sided 51-set family",
        "",
        (
            significant_tiers.to_markdown(index=False, floatfmt=".5g")
            if len(significant_tiers)
            else "None."
        ),
        "",
        "## Broad structural/scaffold control",
        "",
        (
            scaffold[
                ["n_observable_genes", "observed_statistic", "max_absT_fwer"]
            ].to_markdown(index=False, floatfmt=".5g")
            if len(scaffold)
            else "Non-evaluable."
        ),
        "",
        "## Interpretation boundary",
        "",
        "- A positive row means flight-labelled channels are relatively lower "
        "  among parent proteins carrying an external annotation.",
        "- It does not mean that the annotated cells were injured, changed in "
        "  abundance, or were the cells in which the phosphosite was measured.",
        "- Tissue Lcn2, Havcr1, Cst3, Timp2/Igfbp7, collagen transcripts, or "
        "  podocyte annotations are not equivalent to urinary NGAL/KIM-1, GFR, "
        "  NephroCheck, fibrosis assays, albuminuria, or histology.",
        "- If every numeric gate fails, the frozen stop rule retires further "
        "  OSD-462 compartment mining.",
        "",
        "## Reporter-position audit",
        "",
        reporter.to_markdown(index=False, floatfmt=".5g"),
        "",
        "## Observability matched family",
        "",
        matched.to_markdown(index=False, floatfmt=".5g"),
        "",
        "## Contributor concentration",
        "",
        loo[loo["gene_set"].str.endswith(f"__{PRIMARY_TIER}")].to_markdown(
            index=False, floatfmt=".5g"
        ),
        "",
    ]
    path = output / "COMPARTMENT_ADVERSARIAL_DECISION_REPORT.md"
    path.write_text("\n".join(lines))
    return path


def postprocess(
    config_path: str | Path,
    exact_dir: str | Path,
    output_dir: str | Path,
) -> Path:
    cfg = load_config(config_path)
    membership_path = resolve(cfg["input"]["gene_sets"])
    membership = pd.read_csv(membership_path, sep="\t")
    exact_dir = resolve(exact_dir)
    output = resolve(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    genes = load_primary_gene_covariates(cfg)
    exact = compile_exact_family(exact_dir, membership, genes, cfg)
    matched = observability_matched_tests(genes, membership, cfg)
    profiles = profile_direction_summary(cfg, membership)
    loo_detail, loo_summary = leave_one_gene_out_tables(genes, membership, cfg)
    contributors = leading_contributors(genes, membership, cfg)
    grey = grey60_overlap(membership, genes, cfg)
    clinical_detail, clinical_summary = clinical_axis_observability(cfg)
    reporter, block_order = reporter_position_by_layer(cfg)
    decisions = decision_table(exact, matched, profiles, loo_summary, cfg)

    outputs: list[Path] = []
    tables = {
        "exact_family_all_scores.tsv": exact,
        "observability_matched_family.tsv": matched,
        "signal_profile_direction_summary.tsv": profiles,
        "leave_one_gene_out_all_sets.tsv": loo_detail,
        "leave_one_gene_out_summary.tsv": loo_summary,
        "leading_parent_gene_contributors.tsv": contributors,
        "grey60_compartment_overlap.tsv": grey,
        "clinical_axis_observability_gene.tsv": clinical_detail,
        "clinical_axis_observability_summary.tsv": clinical_summary,
        "reporter_position_by_layer.tsv": reporter,
        "reporter_block_order_by_layer.tsv": block_order,
        "compartment_decision_gates.tsv": decisions,
    }
    for name, table in tables.items():
        path = output / name
        table.to_csv(path, sep="\t", index=False)
        outputs.append(path)
    outputs.extend(plot_primary_family(exact, output))
    report = write_decision_report(
        output, decisions, exact, matched, loo_summary, reporter, cfg
    )
    outputs.append(report)

    input_paths = {
        "config": Path(config_path),
        "membership": membership_path,
        "exact_set_results": exact_dir / "set_level_permutation_inference.tsv",
        "exact_manifest": exact_dir / "manifest.json",
        "source_calibration": resolve(cfg["input"]["source_exact_run"])
        / "parent_gene_null_calibration.tsv",
        "source_site_universe": resolve(cfg["input"]["source_exact_run"])
        / "site_universe_primary.tsv",
        "protein_workbook": resolve(cfg["input"]["protein_workbook"]),
        "exact_engine_source": REPO / "src/v13/continuous_phospho_inference.py",
        "audit_source": Path(__file__),
    }
    manifest = {
        "analysis_id": cfg["analysis"]["id"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": cfg["analysis"]["status"],
        "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in input_paths.items()},
        "parameters": {
            "n_exact_label_assignments": 63504,
            "n_exact_family_sets": len(cfg["set_test"]["primary_family"]),
            "observability_draws": cfg["observability_matching"]["n_draws"],
            "observability_seed": cfg["observability_matching"]["seed"],
            "minimum_observable_genes": cfg["set_test"]["minimum_observable_genes"],
        },
        "outputs": {str(path.relative_to(output)): sha256(path) for path in outputs},
        "hard_boundaries": [
            "condition is perfectly aliased with reporter-tag blocks",
            "whole-kidney parent proteins do not identify cell of origin",
            "clinical urine, GFR, histology, electrolyte, and physiology endpoints are absent",
        ],
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("prepare", help="Freeze reference-only marker tiers.")
    post = sub.add_parser("postprocess", help="Compile the completed exact run.")
    post.add_argument("--exact-dir", type=Path, default=DEFAULT_EXACT)
    post.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "prepare":
        path = prepare_reference(args.config)
        print(path)
        return 0
    if args.command == "postprocess":
        path = postprocess(args.config, args.exact_dir, args.output_dir)
        print(path)
        return 0
    raise RuntimeError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
