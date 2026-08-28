#!/usr/bin/env python3
"""Run the frozen cross-mission clinically anchored renal tissue-axis study."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
import scipy
import yaml

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.analysis import (
    combined_score_design,
    descriptive_gene_effects,
    mission_effect_from_score,
    sample_manifest,
    score_mission_axes,
    technical_qc_audit,
)
from src.clinical_axes.data import load_moderator_missions, load_primary_missions
from src.clinical_axes.statistics import (
    blocked_meta_permutation,
    leave_one_mission_out,
    random_effects_reml_mkh,
)


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_OUT = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"


def resolve(path: Path | str) -> Path:
    value = Path(path)
    return value if value.is_absolute() else REPO / value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _holm(p_values: pd.Series) -> pd.Series:
    p = p_values.astype(float)
    order = np.argsort(p.to_numpy())
    adjusted_sorted = np.maximum.accumulate(
        (len(p) - np.arange(len(p))) * p.to_numpy()[order]
    )
    adjusted = np.empty(len(p), dtype=float)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)
    return pd.Series(adjusted, index=p.index, name="p_mkh_holm")


def _write_table(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False)


def _score_components(
    details: dict[str, dict[str, object]],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_rows = []
    subdomain_rows = []
    for mission, axes in details.items():
        for axis, result in axes.items():
            signed = result.signed_gene_z
            for gene in signed.index:
                for sample, value in signed.loc[gene].items():
                    gene_rows.append(
                        {
                            "mission": mission,
                            "axis": axis,
                            "animal": sample,
                            "gene": gene,
                            "signed_z": value,
                        }
                    )
            if result.subdomain_scores is not None:
                for sample, row in result.subdomain_scores.iterrows():
                    for subdomain, value in row.items():
                        subdomain_rows.append(
                            {
                                "mission": mission,
                                "axis": axis,
                                "animal": sample,
                                "subdomain": subdomain,
                                "score": value,
                            }
                        )
    return pd.DataFrame(gene_rows), pd.DataFrame(subdomain_rows)


def _annotate_meta(
    observed: pd.DataFrame,
    mission_effects: pd.DataFrame,
    family: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = observed.reset_index()
    max_weights = {}
    weight_rows = []
    for axis, sub in mission_effects.groupby("axis", sort=False):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        max_weights[axis] = float(fit.weights.max())
        for mission, weight in zip(sub["mission"], fit.weights):
            weight_rows.append(
                {"axis": axis, "mission": mission, "random_effect_weight": weight}
            )
    meta["maximum_weight"] = meta["axis"].map(max_weights)
    meta["role"] = meta["axis"].map(
        {axis: spec["role"] for axis, spec in family.items()}
    )
    holm = _holm(meta.set_index("axis")["p_mkh"])
    meta["p_mkh_holm"] = meta["axis"].map(holm)
    return meta, pd.DataFrame(weight_rows)


def _leave_one_mission(
    mission_effects: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for axis, sub in mission_effects.groupby("axis", sort=False):
        loo = leave_one_mission_out(
            sub["estimate"], sub["variance"], sub["mission"]
        )
        loo.insert(0, "axis", axis)
        rows.append(loo)
    return pd.concat(rows, ignore_index=True)


def _stratum_effects(
    missions,
    scores: pd.DataFrame,
    design: pd.DataFrame,
) -> pd.DataFrame:
    """Expose the age/duration cells that are fixed-effect pooled per mission."""
    rows = []
    for mission, data in missions.items():
        idx = design.index[design["mission"] == mission]
        for axis in scores.columns:
            local = pd.Series(
                scores.loc[idx, axis].to_numpy(),
                index=data.metadata.index,
                name=axis,
            )
            _, strata = mission_effect_from_score(local, data.metadata)
            for row in strata.to_dict(orient="records"):
                rows.append({"axis": axis, "mission": mission, **row})
    return pd.DataFrame(rows)


def _moderator_effects(
    moderators,
    family: dict[str, object],
    cpm_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    effect_rows = []
    score_rows = []
    coverage_rows = []
    for mission, data in moderators.items():
        scores, _, coverage = score_mission_axes(
            data, family, cpm_threshold=cpm_threshold
        )
        coverage_rows.append(coverage)
        for animal, row in scores.iterrows():
            for axis, value in row.items():
                score_rows.append(
                    {
                        "mission": mission,
                        "context": data.context,
                        "animal": animal,
                        "condition": data.metadata.loc[animal, "condition"],
                        "stratum": data.metadata.loc[animal, "block"],
                        "axis": axis,
                        "score": value,
                    }
                )
        for axis in scores:
            summary, strata = mission_effect_from_score(scores[axis], data.metadata)
            effect_rows.append(
                {
                    "mission": mission,
                    "context": data.context,
                    "axis": axis,
                    **summary,
                }
            )
    return (
        pd.DataFrame(effect_rows),
        pd.DataFrame(score_rows),
        pd.concat(coverage_rows, ignore_index=True),
    )


def run(args: argparse.Namespace) -> None:
    config_path = resolve(args.config)
    config = yaml.safe_load(config_path.read_text())
    outdir = resolve(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    family = config["primary_family"]
    threshold = float(config["eligibility"]["cpm_threshold"])

    primary = load_primary_missions(config, REPO)
    qc = technical_qc_audit(primary)
    _write_table(qc, outdir / "technical_qc.tsv")
    coverage_flags = qc[
        (qc["metric"] == config["technical_qc"]["primary_metric"])
        & qc["imbalance_flag"].fillna(False)
    ]["mission"].tolist()
    eligible_primary = {
        name: data for name, data in primary.items() if name not in coverage_flags
    }
    if len(eligible_primary) < 4:
        raise RuntimeError(
            f"Only {len(eligible_primary)} terminal missions passed coverage QC"
        )

    _write_table(sample_manifest(primary), outdir / "sample_manifest.tsv")
    scores, design, coverage, details = combined_score_design(
        eligible_primary,
        family,
        cpm_threshold=threshold,
        method="mean",
    )
    score_out = scores.join(design)
    score_out.index.name = "sample"
    _write_table(score_out.reset_index(), outdir / "animal_axis_scores.tsv")
    _write_table(coverage, outdir / "gene_coverage.tsv")
    gene_z, subdomains = _score_components(details)
    _write_table(gene_z, outdir / "signed_gene_z.tsv")
    _write_table(subdomains, outdir / "subdomain_scores.tsv")

    n_permutations = int(args.permutations or config["inference"]["permutations"])
    primary_result = blocked_meta_permutation(
        scores,
        design,
        n_permutations=n_permutations,
        seed=int(config["seed"]),
        chunk_size=int(args.chunk_size),
    )
    meta, weights = _annotate_meta(
        primary_result.observed_meta,
        primary_result.mission_effects,
        family,
    )
    _write_table(meta, outdir / "primary_meta_results.tsv")
    _write_table(primary_result.mission_effects, outdir / "primary_mission_effects.tsv")
    _write_table(
        _stratum_effects(eligible_primary, scores, design),
        outdir / "primary_stratum_effects.tsv",
    )
    _write_table(weights, outdir / "primary_mission_weights.tsv")
    _write_table(_leave_one_mission(primary_result.mission_effects), outdir / "leave_one_mission.tsv")
    primary_result.null_t.to_csv(
        outdir / "permutation_null_t.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    gene_effects = descriptive_gene_effects(eligible_primary, details)
    _write_table(gene_effects, outdir / "descriptive_signed_gene_effects.tsv")

    # Prespecified median-score sensitivity. It receives the same four-axis
    # family and complete blocked permutation rather than an unadjusted shortcut.
    median_scores, median_design, median_coverage, _ = combined_score_design(
        eligible_primary,
        family,
        cpm_threshold=threshold,
        method="median",
    )
    median_result = blocked_meta_permutation(
        median_scores,
        median_design,
        n_permutations=n_permutations,
        seed=int(config["seed"]) + 1,
        chunk_size=int(args.chunk_size),
    )
    median_meta, _ = _annotate_meta(
        median_result.observed_meta,
        median_result.mission_effects,
        family,
    )
    median_meta.insert(0, "sensitivity", "signed_median")
    _write_table(median_meta, outdir / "median_score_meta_results.tsv")

    moderators = load_moderator_missions(config, REPO)
    moderator_effects, moderator_scores, moderator_coverage = _moderator_effects(
        moderators, family, threshold
    )
    _write_table(moderator_effects, outdir / "recovery_moderator_effects.tsv")
    _write_table(moderator_scores, outdir / "recovery_moderator_scores.tsv")
    _write_table(moderator_coverage, outdir / "recovery_gene_coverage.tsv")

    input_paths = [config_path]
    for spec in config["primary_missions"].values():
        for key in ("vst", "counts", "runsheet", "metadata", "sample_table", "qc"):
            if key in spec:
                input_paths.append(resolve(spec[key]))
    input_paths.extend(
        [
            resolve(config["gene_mapping"]["path"]),
            resolve(config["gene_mapping"]["annotation_fallback"]),
        ]
    )
    hashes = {
        str(path.relative_to(REPO) if path.is_relative_to(REPO) else path): sha256(path)
        for path in dict.fromkeys(input_paths)
        if path.exists()
    }
    manifest = {
        "analysis_id": config["analysis_id"],
        "lock_date": str(config["lock_date"]),
        "status": config["status"],
        "config": str(config_path.relative_to(REPO)),
        "config_sha256": sha256(config_path),
        "primary_missions_loaded": list(primary),
        "primary_missions_excluded_by_coverage_qc": coverage_flags,
        "primary_missions_analyzed": list(eligible_primary),
        "primary_axes": list(family),
        "n_permutations": n_permutations,
        "seed_primary": int(config["seed"]),
        "seed_median": int(config["seed"]) + 1,
        "input_sha256": hashes,
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "interpretation_boundary": config["interpretation_boundary"],
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(meta.to_string(index=False))
    print(f"\nWrote {outdir}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--permutations", type=int, default=None)
    parser.add_argument("--chunk-size", type=int, default=2048)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
