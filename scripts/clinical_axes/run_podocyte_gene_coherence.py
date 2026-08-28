#!/usr/bin/env python3
"""Describe gene-level coherence of the cross-mission podocyte RNA program.

This is a descriptive follow-up to the multiplicity-controlled compartment
family.  Gene-level intervals and p-values are retained for auditability but
must not be interpreted as a second confirmatory testing family.
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
from src.clinical_axes.data import load_primary_missions  # noqa: E402
from src.clinical_axes.statistics import random_effects_reml_mkh  # noqa: E402


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_GENE_SETS = (
    REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
)
DEFAULT_RESULTS = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"
SET_NAME = "podocyte__high_specificity"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(config_path: Path, gene_sets_path: Path, results: Path) -> None:
    config = yaml.safe_load(config_path.read_text())
    table = pd.read_csv(gene_sets_path, sep="\t")
    selected = table[
        table["gene_set"].eq(SET_NAME)
        & table["final_for_testing"].astype(bool)
    ].copy()
    genes = list(dict.fromkeys(selected["gene_symbol"].dropna().astype(str)))
    if len(genes) < 8:
        raise ValueError(f"{SET_NAME} contains only {len(genes)} genes")

    family = {
        SET_NAME: {
            "role": "descriptive_gene_coherence",
            "subdomains": {
                "atlas_markers": {
                    "genes": {gene: 1 for gene in genes},
                    "minimum_present": 8,
                }
            },
        }
    }
    missions = load_primary_missions(config, REPO)
    _, _, coverage, details = combined_score_design(
        missions,
        family,
        cpm_threshold=float(config["eligibility"]["cpm_threshold"]),
    )
    mission_gene = descriptive_gene_effects(missions, details)

    rows: list[dict[str, object]] = []
    for gene, sub in mission_gene.groupby("gene", sort=False):
        if sub["mission"].nunique() < 2:
            continue
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        rows.append(
            {
                "gene": gene,
                "n_missions": int(sub["mission"].nunique()),
                "pooled_hedges_g": fit.estimate,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh_descriptive_unadjusted": fit.p,
                "i_squared": fit.i_squared,
                "n_positive_missions": int(sub["estimate"].gt(0).sum()),
                "n_negative_missions": int(sub["estimate"].lt(0).sum()),
                "minimum_mission_effect": float(sub["estimate"].min()),
                "maximum_mission_effect": float(sub["estimate"].max()),
            }
        )
    gene_meta = pd.DataFrame(rows).sort_values(
        ["pooled_hedges_g", "gene"], ascending=[False, True]
    )

    summary = pd.DataFrame(
        [
            {
                "gene_set": SET_NAME,
                "n_defined": len(genes),
                "n_meta_analyzable": len(gene_meta),
                "n_observable_all_five_missions": int(
                    gene_meta["n_missions"].eq(len(missions)).sum()
                ),
                "median_pooled_hedges_g": gene_meta.loc[
                    gene_meta["n_missions"].eq(len(missions)), "pooled_hedges_g"
                ].median(),
                "n_positive_pooled_estimate": int(
                    gene_meta.loc[
                        gene_meta["n_missions"].eq(len(missions)),
                        "pooled_hedges_g",
                    ].gt(0).sum()
                ),
                "n_positive_at_least_four_missions": int(
                    gene_meta.loc[
                        gene_meta["n_missions"].eq(len(missions)),
                        "n_positive_missions",
                    ].ge(4).sum()
                ),
                "n_positive_all_five_missions": int(
                    gene_meta.loc[
                        gene_meta["n_missions"].eq(len(missions)),
                        "n_positive_missions",
                    ].eq(len(missions)).sum()
                ),
            }
        ]
    )

    results.mkdir(parents=True, exist_ok=True)
    mission_gene.to_csv(
        results / "podocyte_gene_coherence_mission_effects.tsv",
        sep="\t",
        index=False,
    )
    gene_meta.to_csv(
        results / "podocyte_gene_coherence_meta.tsv", sep="\t", index=False
    )
    summary.to_csv(
        results / "podocyte_gene_coherence_summary.tsv", sep="\t", index=False
    )
    coverage.to_csv(
        results / "podocyte_gene_coherence_coverage.tsv", sep="\t", index=False
    )
    manifest = {
        "analysis": "descriptive cross-mission podocyte gene coherence",
        "status": "post_hoc_descriptive_not_confirmatory",
        "gene_set": SET_NAME,
        "config": str(config_path),
        "config_sha256": _sha256(config_path),
        "gene_sets": str(gene_sets_path),
        "gene_sets_sha256": _sha256(gene_sets_path),
        "interpretation": (
            "Gene-level estimates describe the coherence and leading edge of the "
            "multiplicity-controlled set result. Individual p-values are unadjusted "
            "and are not a separate discovery family."
        ),
    }
    (results / "podocyte_gene_coherence_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(summary.to_string(index=False))
    print("\nLeading positive genes")
    print(
        gene_meta[gene_meta["n_missions"].eq(len(missions))]
        .head(25)
        .to_string(index=False)
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--gene-sets", type=Path, default=DEFAULT_GENE_SETS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    args = parser.parse_args()
    run(args.config, args.gene_sets, args.results)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
