#!/usr/bin/env python3
"""Apply the frozen kidney-compartment atlas family to terminal RNA cohorts.

Default (``--method mean --gene-universe per_mission``) is the frozen primary
compartment-context family (seed = config seed + 100), written to ``--results``.

Family-wide scoring variants (podocyte sensitivity items 1 and family-level 2):

* ``--method median``: per-animal signed median of gene z-scores instead of the
  mean. The family (which sets are evaluable) is fixed before scoring, so it is
  identical to the mean family.
* ``--gene-universe common``: every set is restricted to genes that are
  CPM-eligible in every mission and present in every expression matrix; all of
  those genes are required in every mission and a set is kept only if at least
  8 remain.

Variants use their own seeds (config seed + 101 median, +102 common, +103
median+common), are written with the same filenames to
``{results}/compartment_context_{median|common|median_common}/`` (or to
``--results`` itself when it already names that folder), and refuse to
run unless the preregistration YAML is committed (``--allow-uncommitted-prereg``
is for tests only).
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Mapping

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
DEFAULT_PREREG = REPO / "config/clinical_axes_podocyte_sensitivities.yaml"

SCORING_METHODS = ("mean", "median")
GENE_UNIVERSES = ("per_mission", "common")
# Offsets added to the config seed; the frozen primary family keeps +100.
SEED_OFFSETS = {
    ("mean", "per_mission"): 100,
    ("median", "per_mission"): 101,
    ("mean", "common"): 102,
    ("median", "common"): 103,
}
MINIMUM_GENES = 8

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
    """Add barrier-depleted podocyte sensitivity sets into the same corrected max-|T| family, not as separate tests."""

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


def seed_offset(method: str = "mean", gene_universe: str = "per_mission") -> int:
    """Offset added to the config seed for a scoring variant (default 100)."""
    try:
        return SEED_OFFSETS[(method, gene_universe)]
    except KeyError as error:
        raise ValueError(
            f"unknown scoring variant method={method!r}, gene_universe={gene_universe!r}"
        ) from error


def is_default_variant(method: str = "mean", gene_universe: str = "per_mission") -> bool:
    return (method, gene_universe) == ("mean", "per_mission")


def variant_output_dir(
    results: Path, method: str = "mean", gene_universe: str = "per_mission"
) -> Path:
    """The primary family writes to ``results``; each variant to its own subfolder
    ``compartment_context_{median|common|median_common}``. If ``results`` already
    names that subfolder (registry stages pass it explicitly so no two stages share
    an output directory), it is used as is rather than nested again."""
    seed_offset(method, gene_universe)  # validates the pair
    results = Path(results)
    if is_default_variant(method, gene_universe):
        return results
    parts = []
    if method != "mean":
        parts.append(method)
    if gene_universe != "per_mission":
        parts.append(gene_universe)
    subdir = "compartment_context_" + "_".join(parts)
    return results if results.name == subdir else results / subdir


def common_gene_universe(missions: Mapping[str, object], threshold: float) -> set[str]:
    """Genes CPM-eligible in every mission and present in every expression matrix."""
    common = set.intersection(
        *(set(cpm_eligible_genes(data, threshold)) for data in missions.values())
    )
    common &= set.intersection(
        *(set(data.expression.index) for data in missions.values())
    )
    return common


def restrict_family_to_common(
    family, definition_audit: pd.DataFrame, common: set[str], minimum: int = MINIMUM_GENES
):
    """Restrict every set to ``set & common``, require all retained genes
    (minimum_present = number retained), keep a set only if >= ``minimum`` genes
    remain, and record genes dropped and the keep decision in the audit."""
    restricted = {}
    kept, dropped, n_left = {}, {}, {}
    for gene_set, spec in family.items():
        spec = copy.deepcopy(spec)
        n_before = 0
        n_after = 0
        for subdomain in spec["subdomains"].values():
            before = subdomain["genes"]
            after = {gene: sign for gene, sign in before.items() if gene in common}
            n_before += len(before)
            n_after += len(after)
            subdomain["genes"] = after
            subdomain["minimum_present"] = len(after)
        keep = n_after >= minimum and all(
            len(sub["genes"]) > 0 for sub in spec["subdomains"].values()
        )
        kept[gene_set] = keep
        dropped[gene_set] = n_before - n_after
        n_left[gene_set] = n_after
        if keep:
            restricted[gene_set] = spec
    audit = definition_audit.copy()
    audit["n_common_universe_genes"] = audit["gene_set"].map(n_left).astype("Int64")
    audit["n_dropped_common_universe"] = audit["gene_set"].map(dropped).astype("Int64")
    audit["common_universe_kept"] = audit["gene_set"].map(
        lambda gene_set: bool(kept.get(gene_set, False))
    )
    return restricted, audit


def prepare_family(
    missions: Mapping[str, object],
    gene_sets: Path,
    threshold: float,
    *,
    add_podocyte_disjoint: bool = False,
    gene_universe: str = "per_mission",
):
    """Build the evaluable family and its definition audit before any scoring.

    The scoring method plays no part here, so mean and median share one family.
    """
    if gene_universe not in GENE_UNIVERSES:
        raise ValueError(f"gene_universe must be one of {GENE_UNIVERSES}")
    family, definition_audit = load_family(gene_sets)
    if add_podocyte_disjoint:
        family, definition_audit = add_disjoint_podocyte_variants(
            family, definition_audit
        )
    if gene_universe == "common":
        family, definition_audit = restrict_family_to_common(
            family, definition_audit, common_gene_universe(missions, threshold)
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
        if min(counts.values()) >= MINIMUM_GENES:
            retained[gene_set] = spec
    family = retained
    definition_audit["cross_mission_evaluable"] = definition_audit["gene_set"].isin(
        family
    )
    for mission in missions:
        definition_audit[f"n_eligible_{mission}"] = definition_audit["gene_set"].map(
            lambda gene_set: cross_counts.get(gene_set, {}).get(mission, 0)
        )
    return family, definition_audit


def _git_head() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True, check=False
        )
    except OSError:
        return None
    return result.stdout.strip() or None


def prereg_guard(path: Path, allow_uncommitted: bool = False) -> dict[str, object]:
    """Refuse a variant run unless the preregistration YAML exists and has no
    uncommitted changes (``git status --porcelain <path>`` empty). Returns the
    manifest record (YAML sha256, git HEAD, guard outcome)."""
    path = Path(path)
    problem = None
    if not path.exists():
        problem = "preregistration YAML is missing"
    else:
        try:
            status = subprocess.run(
                ["git", "status", "--porcelain", "--", str(path)],
                cwd=REPO,
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as error:
            problem = f"git status failed: {error}"
        else:
            if status.returncode != 0:
                problem = "git status failed: " + status.stderr.strip()
            elif status.stdout.strip():
                problem = "uncommitted changes: " + status.stdout.strip()
            else:
                tracked = subprocess.run(
                    ["git", "ls-files", "--error-unmatch", "--", str(path)],
                    cwd=REPO,
                    capture_output=True,
                    text=True,
                    check=False,
                )
                if tracked.returncode != 0:
                    problem = "preregistration YAML is not tracked by git"
    if problem is not None and not allow_uncommitted:
        raise SystemExit(
            f"refusing to run a preregistered variant: {problem} ({path}); commit the "
            "preregistration first (--allow-uncommitted-prereg is for tests only)"
        )
    return {
        "prereg_path": str(path),
        "prereg_sha256": sha256(path) if path.exists() else None,
        "git_head": _git_head(),
        "prereg_guard": "clean" if problem is None else f"bypassed: {problem}",
    }


def apply_gene_map_override(config: dict, gene_map: Path | None) -> dict:
    """Point config["gene_mapping"]["path"] at an alternative ID map, in memory only
    (the config file is never edited). Returns the config for chaining."""
    if gene_map is not None:
        resolved = Path(gene_map).expanduser().resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"--gene-map {resolved} does not exist")
        config["gene_mapping"]["path"] = str(resolved)
    return config


def run(args):
    method = getattr(args, "method", "mean")
    gene_universe = getattr(args, "gene_universe", "per_mission")
    offset = seed_offset(method, gene_universe)
    prereg_record = None
    if not is_default_variant(method, gene_universe):
        prereg_record = prereg_guard(
            getattr(args, "prereg", DEFAULT_PREREG),
            bool(getattr(args, "allow_uncommitted_prereg", False)),
        )
    config = yaml.safe_load(args.config.read_text())
    config = apply_gene_map_override(config, getattr(args, "gene_map", None))
    gene_map_path = REPO / str(config["gene_mapping"]["path"])
    threshold = float(config["eligibility"]["cpm_threshold"])
    seed = int(config["seed"]) + offset
    missions = load_primary_missions(config, REPO)
    if args.osd253_strain is not None:
        missions["OSD-253"] = load_osd253_strain_sensitivity(
            config, REPO, strain=args.osd253_strain
        )
    family, definition_audit = prepare_family(
        missions,
        args.gene_sets,
        threshold,
        add_podocyte_disjoint=bool(args.add_podocyte_disjoint_variants),
        gene_universe=gene_universe,
    )
    scores, design, coverage, _ = combined_score_design(
        missions, family, cpm_threshold=threshold, method=method
    )
    result = blocked_meta_permutation(
        scores,
        design,
        n_permutations=args.permutations,
        seed=seed,
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
    out = variant_output_dir(args.results, method, gene_universe)
    out.mkdir(parents=True, exist_ok=True)
    meta.to_csv(
        out / "compartment_context_meta_results.tsv", sep="\t", index=False
    )
    result.mission_effects.to_csv(
        out / "compartment_context_mission_effects.tsv",
        sep="\t",
        index=False,
    )
    coverage.to_csv(
        out / "compartment_context_gene_coverage.tsv",
        sep="\t",
        index=False,
    )
    definition_audit.to_csv(
        out / "compartment_context_definition_audit.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(weights).to_csv(
        out / "compartment_context_weights.tsv", sep="\t", index=False
    )
    result.null_t.to_csv(
        out / "compartment_context_null_t.tsv.gz",
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
        "seed": seed,
        "seed_offset": offset,
        "scoring_method": method,
        "gene_universe": gene_universe,
        "gene_universe_definition": (
            "per_mission: genes CPM-eligible and present in each mission separately"
            if gene_universe == "per_mission"
            else "common: set restricted to genes CPM-eligible in every mission and "
            "present in every expression matrix; all retained genes required; set "
            f"kept only if >= {MINIMUM_GENES} genes remain"
        ),
        "gene_map": str(gene_map_path),
        "gene_map_sha256": sha256(gene_map_path) if gene_map_path.exists() else None,
        "gene_map_overridden": getattr(args, "gene_map", None) is not None,
        "output_dir": str(out),
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
    if prereg_record is not None:
        manifest["preregistration"] = prereg_record
    (out / "compartment_context_manifest.json").write_text(
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
    parser.add_argument(
        "--method",
        choices=SCORING_METHODS,
        default="mean",
        help="Per-animal aggregation of signed gene z-scores (median: item 1 sensitivity).",
    )
    parser.add_argument(
        "--gene-universe",
        choices=GENE_UNIVERSES,
        default="per_mission",
        help=(
            "per_mission (primary) or common: restrict every set to genes eligible "
            "and observed in all missions (family-level item 2 sensitivity)."
        ),
    )
    parser.add_argument(
        "--gene-map",
        type=Path,
        default=None,
        help=(
            "Alternative Ensembl->symbol ID map; replaces config gene_mapping.path in "
            "memory only (the config file is not edited)."
        ),
    )
    parser.add_argument(
        "--prereg",
        type=Path,
        default=DEFAULT_PREREG,
        help="Preregistration YAML that must be committed before a variant run.",
    )
    parser.add_argument(
        "--allow-uncommitted-prereg",
        action="store_true",
        help="Tests only: run a variant although the preregistration is missing or dirty.",
    )
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
