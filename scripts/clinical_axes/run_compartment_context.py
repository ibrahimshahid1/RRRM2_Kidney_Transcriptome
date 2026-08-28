#!/usr/bin/env python3
"""Apply the frozen kidney-compartment atlas family to terminal RNA cohorts."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.analysis import combined_score_design  # noqa: E402
from src.clinical_axes.data import (  # noqa: E402
    cpm_eligible_genes,
    load_osd253_strain_sensitivity,
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
DEFAULT_RESULTS = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"

BARRIER_CORE_GENES = (
    "Nphs1",
    "Nphs2",
    "Synpo",
    "Ptpro",
    "Magi2",
    "Wt1",
)
BARRIER_EXPANDED_GENES = BARRIER_CORE_GENES + ("Podxl", "Cd2ap")
PODOCYTE_HIGH_SPECIFICITY = "podocyte__high_specificity"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_family(path: Path, minimum_defined: int = 8):
    table = pd.read_csv(path, sep="\t")
    table = table[table["final_for_testing"].astype(bool)].copy()
    family = {}
    audit = []
    for gene_set, sub in table.groupby("gene_set", sort=False):
        genes = list(dict.fromkeys(sub["gene_symbol"].dropna().astype(str)))
        evaluable = len(genes) >= minimum_defined
        audit.append(
            {
                "gene_set": gene_set,
                "report_compartment": sub["report_compartment"].iloc[0],
                "tier": sub["tier"].iloc[0],
                "n_defined": len(genes),
                "definition_evaluable": evaluable,
            }
        )
        if evaluable:
            family[gene_set] = {
                "role": "secondary_compartment_context",
                "subdomains": {
                    "atlas_markers": {
                        "genes": {gene: 1 for gene in genes},
                        "minimum_present": 8,
                    }
                },
            }
    return family, pd.DataFrame(audit)


def add_disjoint_podocyte_variants(family, definition_audit):
    """Add two conservative high-specificity podocyte sensitivity sets.

    The original atlas-defined family remains intact.  The two added sets use
    exactly the original high-specificity podocyte definition after removing
    (i) the six frozen barrier-core genes or (ii) those six plus the two
    expanded barrier markers.  Keeping the original family members and adding
    both sensitivities makes max-|T| inference conservative over the complete
    51-evaluable-set family rather than treating either sensitivity as a
    separately corrected test.
    """

    if PODOCYTE_HIGH_SPECIFICITY not in family:
        raise ValueError(
            f"Required source set {PODOCYTE_HIGH_SPECIFICITY!r} is not evaluable"
        )
    family = copy.deepcopy(family)
    audit_rows = []
    source_genes = family[PODOCYTE_HIGH_SPECIFICITY]["subdomains"][
        "atlas_markers"
    ]["genes"]
    variants = {
        "podocyte__high_specificity_disjoint_barrier_core": BARRIER_CORE_GENES,
        "podocyte__high_specificity_disjoint_barrier_expanded": (
            BARRIER_EXPANDED_GENES
        ),
    }
    for gene_set, requested_exclusions in variants.items():
        requested = set(requested_exclusions)
        retained_genes = {
            gene: sign for gene, sign in source_genes.items() if gene not in requested
        }
        present = sorted(requested.intersection(source_genes))
        absent = sorted(requested.difference(source_genes))
        spec = copy.deepcopy(family[PODOCYTE_HIGH_SPECIFICITY])
        spec["subdomains"]["atlas_markers"]["genes"] = retained_genes
        family[gene_set] = spec
        audit_rows.append(
            {
                "gene_set": gene_set,
                "report_compartment": "podocyte",
                "tier": gene_set.removeprefix("podocyte__"),
                "n_defined": len(retained_genes),
                "definition_evaluable": len(retained_genes) >= 8,
                "derivation": f"{PODOCYTE_HIGH_SPECIFICITY} minus frozen markers",
                "requested_excluded_genes": "|".join(requested_exclusions),
                "excluded_genes_present": "|".join(present),
                "excluded_genes_absent_from_source_set": "|".join(absent),
            }
        )
    return family, pd.concat(
        [definition_audit, pd.DataFrame(audit_rows)], ignore_index=True, sort=False
    )


def run(args):
    config = yaml.safe_load(args.config.read_text())
    threshold = float(config["eligibility"]["cpm_threshold"])
    missions = load_primary_missions(config, REPO)
    if args.osd253_strain is not None:
        missions["OSD-253"] = load_osd253_strain_sensitivity(
            config, REPO, strain=args.osd253_strain
        )
    family, definition_audit = load_family(args.gene_sets)
    if args.add_podocyte_disjoint_variants:
        family, definition_audit = add_disjoint_podocyte_variants(
            family, definition_audit
        )
    eligible_by_mission = {
        mission: cpm_eligible_genes(data, threshold)
        for mission, data in missions.items()
    }
    retained = {}
    cross_counts = {}
    for gene_set, spec in family.items():
        genes = set(spec["subdomains"]["atlas_markers"]["genes"])
        counts = {
            mission: len(genes.intersection(eligible))
            for mission, eligible in eligible_by_mission.items()
        }
        cross_counts[gene_set] = counts
        if min(counts.values()) >= 8:
            retained[gene_set] = spec
    family = retained
    definition_audit["cross_mission_evaluable"] = definition_audit["gene_set"].isin(
        family
    )
    for mission in missions:
        definition_audit[f"n_eligible_{mission}"] = definition_audit["gene_set"].map(
            lambda gene_set: cross_counts.get(gene_set, {}).get(mission, 0)
        )
    scores, design, coverage, _ = combined_score_design(
        missions, family, cpm_threshold=threshold
    )
    result = blocked_meta_permutation(
        scores,
        design,
        n_permutations=args.permutations,
        seed=int(config["seed"]) + 100,
        chunk_size=args.chunk_size,
    )
    meta = result.observed_meta.reset_index()
    lookup = definition_audit.set_index("gene_set")
    meta["report_compartment"] = meta["axis"].map(
        lookup["report_compartment"]
    )
    meta["tier"] = meta["axis"].map(lookup["tier"])
    weights = []
    for axis, sub in result.mission_effects.groupby("axis", sort=False):
        fit = random_effects_reml_mkh(sub["estimate"], sub["variance"])
        meta.loc[meta["axis"] == axis, "maximum_weight"] = fit.weights.max()
        for mission, weight in zip(sub["mission"], fit.weights):
            weights.append(
                {"axis": axis, "mission": mission, "random_effect_weight": weight}
            )
    args.results.mkdir(parents=True, exist_ok=True)
    meta.to_csv(
        args.results / "compartment_context_meta_results.tsv", sep="\t", index=False
    )
    result.mission_effects.to_csv(
        args.results / "compartment_context_mission_effects.tsv",
        sep="\t",
        index=False,
    )
    coverage.to_csv(
        args.results / "compartment_context_gene_coverage.tsv",
        sep="\t",
        index=False,
    )
    definition_audit.to_csv(
        args.results / "compartment_context_definition_audit.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(weights).to_csv(
        args.results / "compartment_context_weights.tsv", sep="\t", index=False
    )
    result.null_t.to_csv(
        args.results / "compartment_context_null_t.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    manifest = {
        "analysis": "frozen cross-mission kidney-compartment context family",
        "status": "secondary_compartment_context",
        "config": str(args.config),
        "config_sha256": sha256(args.config),
        "gene_sets": str(args.gene_sets),
        "gene_sets_sha256": sha256(args.gene_sets),
        "n_defined_sets": int(len(definition_audit)),
        "n_evaluable_sets": int(len(family)),
        "n_permutations": int(args.permutations),
        "seed": int(config["seed"]) + 100,
        "podocyte_disjoint_variants_added": bool(
            args.add_podocyte_disjoint_variants
        ),
        "osd253_strain": args.osd253_strain or "C57BL/6J_primary",
        "barrier_core_exclusions_requested": list(BARRIER_CORE_GENES),
        "barrier_expanded_exclusions_requested": list(BARRIER_EXPANDED_GENES),
        "multiplicity": "maximum absolute REML/mKH t over all evaluable frozen sets",
        "interpretation_boundary": (
            "bulk-kidney compartment-associated transcript abundance; not cell counts, "
            "cell localization, injury, or function"
        ),
    }
    (args.results / "compartment_context_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    cols = [
        "axis",
        "report_compartment",
        "tier",
        "estimate",
        "ci_low_mkh",
        "ci_high_mkh",
        "i_squared",
        "max_t_fwer",
    ]
    print(meta.sort_values("max_t_fwer")[cols].head(20).to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--gene-sets", type=Path, default=DEFAULT_GENE_SETS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--permutations", type=int, default=20_000)
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument(
        "--add-podocyte-disjoint-variants",
        action="store_true",
        help=(
            "Add high-specificity podocyte variants excluding the frozen "
            "six-gene barrier core and expanded eight-gene panel to the same "
            "max-|T| family."
        ),
    )
    parser.add_argument(
        "--osd253-strain",
        default=None,
        help=(
            "Sensitivity-only replacement strain for the OSD-253 original-"
            "control contrast (for example, C3H/HeJ)."
        ),
    )
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
