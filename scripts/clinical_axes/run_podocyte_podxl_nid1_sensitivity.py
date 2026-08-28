#!/usr/bin/env python3
"""Force Podxl and Nid1 into the cross-mission podocyte RNA program.

This is a post-hoc, adversarial sensitivity analysis prompted by the PODXL/NID1
comparison with prior spaceflight literature.  The frozen high-specificity
podocyte set already contains Podxl; this script adds Nid1, verifies that both
genes pass the frozen CPM eligibility rule in every mission, and reruns the
animal-level blocked-label permutation.  It also reports Podxl and Nid1
separately so a stable set result cannot conceal discordant individual genes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.analysis import (  # noqa: E402
    combined_score_design,
    descriptive_gene_effects,
)
from src.clinical_axes.data import (  # noqa: E402
    cpm_eligible_genes,
    load_primary_missions,
)
from src.clinical_axes.statistics import (  # noqa: E402
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
    / "podxl_nid1_forced_sensitivity"
)
SOURCE_SET = "podocyte__high_specificity"
ORIGINAL_AXIS = "podocyte__high_specificity_original"
FORCED_AXIS = "podocyte__high_specificity_plus_Podxl_Nid1"
FORCED_GENES = ("Podxl", "Nid1")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_source_genes(path: Path) -> list[str]:
    table = pd.read_csv(path, sep="\t")
    selected = table[
        table["gene_set"].eq(SOURCE_SET)
        & table["final_for_testing"].astype(bool)
    ]
    genes = list(dict.fromkeys(selected["gene_symbol"].dropna().astype(str)))
    if len(genes) < 8:
        raise ValueError(f"{SOURCE_SET} contains only {len(genes)} genes")
    return genes


def _axis_spec(genes: list[str]) -> dict[str, object]:
    return {
        "role": "post_hoc_adversarial_sensitivity",
        "subdomains": {
            "atlas_markers": {
                "genes": {gene: 1 for gene in genes},
                "minimum_present": 8,
            }
        },
    }


def _forced_observability(missions, threshold: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for mission, data in missions.items():
        library = data.counts.sum(axis=0)
        cpm = data.counts.divide(library, axis=1) * 1_000_000.0
        eligible = cpm_eligible_genes(data, threshold)
        for gene in FORCED_GENES:
            in_expression = gene in data.expression.index
            in_counts = gene in data.counts.index
            row: dict[str, object] = {
                "mission": mission,
                "gene": gene,
                "in_expression": in_expression,
                "in_counts": in_counts,
                "cpm_eligible": gene in eligible,
                "cpm_threshold": threshold,
            }
            for condition in ("FLT", "GC"):
                samples = data.metadata.index[
                    data.metadata["condition"].eq(condition)
                ]
                if in_counts:
                    values = cpm.loc[gene, samples]
                    row[f"n_{condition.lower()}"] = len(samples)
                    row[f"n_{condition.lower()}_at_or_above_threshold"] = int(
                        values.ge(threshold).sum()
                    )
                    row[f"median_cpm_{condition.lower()}"] = float(values.median())
                else:
                    row[f"n_{condition.lower()}"] = len(samples)
                    row[f"n_{condition.lower()}_at_or_above_threshold"] = 0
                    row[f"median_cpm_{condition.lower()}"] = float("nan")
            rows.append(row)
    return pd.DataFrame(rows)


def _gene_meta(mission_gene: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for gene, sub in mission_gene.groupby("gene", sort=False):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        rows.append(
            {
                "gene": gene,
                "n_missions": int(sub["mission"].nunique()),
                "pooled_hedges_g": fit.estimate,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh_descriptive_unadjusted": fit.p,
                "tau2": fit.tau2,
                "i_squared": fit.i_squared,
                "prediction_low": fit.prediction_low,
                "prediction_high": fit.prediction_high,
                "n_positive_missions": int(sub["estimate"].gt(0).sum()),
                "n_negative_missions": int(sub["estimate"].lt(0).sum()),
            }
        )
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> None:
    config = yaml.safe_load(args.config.read_text())
    threshold = float(config["eligibility"]["cpm_threshold"])
    source_genes = _load_source_genes(args.gene_sets)
    forced_genes = list(dict.fromkeys([*source_genes, *FORCED_GENES]))

    missions = load_primary_missions(config, REPO)
    observability = _forced_observability(missions, threshold)
    failed = observability[
        ~(
            observability["in_expression"].astype(bool)
            & observability["in_counts"].astype(bool)
            & observability["cpm_eligible"].astype(bool)
        )
    ]
    if not failed.empty:
        detail = failed[["mission", "gene"]].to_dict(orient="records")
        raise ValueError(f"Forced genes are not observable in every mission: {detail}")

    family = {
        ORIGINAL_AXIS: _axis_spec(source_genes),
        FORCED_AXIS: _axis_spec(forced_genes),
    }
    scores, design, coverage, details = combined_score_design(
        missions,
        family,
        cpm_threshold=threshold,
    )
    forced_coverage = coverage[coverage["axis"].eq(FORCED_AXIS)]
    for row in forced_coverage.itertuples(index=False):
        used = set(str(row.genes_used).split("|"))
        missing_forced = sorted(set(FORCED_GENES) - used)
        if missing_forced:
            raise AssertionError(
                f"{row.mission}: forced genes absent from score: {missing_forced}"
            )

    permutation = blocked_meta_permutation(
        scores,
        design,
        n_permutations=args.permutations,
        seed=int(config["seed"]) + 503,
        chunk_size=args.chunk_size,
    )
    set_meta = permutation.observed_meta.reset_index()

    mission_gene = descriptive_gene_effects(missions, details)
    mission_gene = mission_gene[
        mission_gene["axis"].eq(FORCED_AXIS)
        & mission_gene["gene"].isin(FORCED_GENES)
    ].copy()
    gene_meta = _gene_meta(mission_gene)

    args.results.mkdir(parents=True, exist_ok=True)
    set_meta.to_csv(args.results / "set_meta_results.tsv", sep="\t", index=False)
    permutation.mission_effects.to_csv(
        args.results / "set_mission_effects.tsv", sep="\t", index=False
    )
    coverage.to_csv(args.results / "set_gene_coverage.tsv", sep="\t", index=False)
    observability.to_csv(
        args.results / "forced_gene_observability.tsv", sep="\t", index=False
    )
    mission_gene.to_csv(
        args.results / "forced_gene_mission_effects.tsv", sep="\t", index=False
    )
    gene_meta.to_csv(
        args.results / "forced_gene_meta_results.tsv", sep="\t", index=False
    )
    permutation.null_t.to_csv(
        args.results / "set_null_t.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    manifest = {
        "analysis": "Podxl/Nid1-forced podocyte program sensitivity",
        "status": "post_hoc_adversarial_sensitivity_not_new_confirmatory_family",
        "source_set": SOURCE_SET,
        "source_set_gene_count": len(source_genes),
        "source_set_contains_Podxl": "Podxl" in source_genes,
        "source_set_contains_Nid1": "Nid1" in source_genes,
        "forced_set_gene_count": len(forced_genes),
        "forced_genes": list(FORCED_GENES),
        "forced_genes_required_eligible_in_every_mission": True,
        "config": str(args.config),
        "config_sha256": _sha256(args.config),
        "gene_sets": str(args.gene_sets),
        "gene_sets_sha256": _sha256(args.gene_sets),
        "n_permutations": int(args.permutations),
        "seed": int(config["seed"]) + 503,
        "permutation_family": [ORIGINAL_AXIS, FORCED_AXIS],
        "multiplicity": "maximum absolute REML/mKH t over original and forced sets",
        "gene_level_status": (
            "descriptive unadjusted estimates included to expose concordance or "
            "discordance; not individual discoveries"
        ),
        "interpretation_boundary": (
            "bulk-kidney podocyte-associated transcript abundance; not urinary "
            "protein, podocyte localization, injury, or filtration function"
        ),
    }
    (args.results / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )

    print("Set-level results")
    print(
        set_meta[
            [
                "axis",
                "estimate",
                "ci_low_mkh",
                "ci_high_mkh",
                "p_mkh",
                "empirical_p_two_sided",
                "max_t_fwer",
                "i_squared",
            ]
        ].to_string(index=False)
    )
    print("\nForced-gene observability")
    print(observability.to_string(index=False))
    print("\nForced-gene mission effects")
    print(
        mission_gene[["mission", "gene", "estimate", "ci_low", "ci_high"]]
        .sort_values(["gene", "mission"])
        .to_string(index=False)
    )
    print("\nForced-gene pooled estimates")
    print(gene_meta.to_string(index=False))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--gene-sets", type=Path, default=DEFAULT_GENE_SETS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--chunk-size", type=int, default=512)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
