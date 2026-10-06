#!/usr/bin/env python3
"""Item 8: leave-one-atlas-source-study-out reconstruction of the compartment marker tiers.

The 51 frozen marker sets (the 49-set evaluable family that contains
``podocyte__high_specificity`` plus the structural scaffold control) were built from
the Mouse Kidney Atlas pooled over its eight ``Origin`` source studies.  This script
asks whether the podocyte high-specificity definition, and the flight effect measured
with it, depends on any one source study.

For each omitted study it

1. rebuilds the compartment expression table from the source-study pseudobulks
   (``compartment_expression_from_pseudobulk``; exactly the arithmetic of
   ``03_atlas_pseudobulk.py``) and the marker tiers with the frozen thresholds
   (``build_marker_tiers``; thresholds are never tuned);
2. applies the frozen atlas QC rule (>= 3 source studies per required compartment);
   a failing rebuild is computed but excluded from the gate; a rebuild without a
   usable podocyte set (< 8 defined or < 8 observable genes) is NOT_EVALUABLE;
3. (scoring only) scores the rebuilt family in the five primary missions through the
   same code path as ``run_compartment_context.py`` (blocked max-|T| permutation,
   seed = parent + 800 + j), and K0 model P1 (rebuilt podocyte set adjusted for the
   structural scaffold minus the rebuilt podocyte genes; HC3 per mission, REML/mKH).

Steps 1-2 are flight-label blind and run with ``--skip-scoring`` at any time.  Step 3
reads mission data and, like every real-data script of the sensitivity programme,
refuses to run unless the preregistration YAML exists and has no uncommitted changes
(``--allow-uncommitted-prereg`` is for tests only).

Dropping a study lowers the number of studies contributing to each compartment, so
the frozen fractional rules become effectively stricter (the high-specificity
detection fraction >= 0.75 is 3-of-4 studies for the full podocyte compartment but
3-of-3 once one of the four podocyte-contributing studies is omitted, and those
rebuilds sit at the frozen QC minimum of 3 studies).  Thresholds are left unchanged;
results are therefore reported separately for podocyte-contributing and
non-contributing omitted studies.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import datetime
import json
from pathlib import Path
import re
import sys
from typing import Any, Mapping

import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.clinical_axes.run_compartment_context import load_family  # noqa: E402
from scripts.clinical_axes.run_podocyte_scaffold_specificity import (  # noqa: E402
    _frozen_genes,
    _regression_rows,
)
from src.clinical_axes import atlas_loso as al  # noqa: E402
from src.clinical_axes.analysis import combined_score_design  # noqa: E402
from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions  # noqa: E402
from src.clinical_axes.statistics import (  # noqa: E402
    blocked_meta_permutation,
    random_effects_reml_mkh,
    score_signed_axis,
)
from src.v13.compartment_adversarial_audit import load_config  # noqa: E402

ANALYSIS_ID = "clinical_axes_atlas_loso_markers_v1"
DEFAULT_AUDIT_CONFIG = REPO / "config/v13_compartment_adversarial_audit.yaml"
DEFAULT_FREEZE_CONFIG = REPO / "config/dct_subtype_reference_freeze_v1.yaml"
DEFAULT_ANALYSIS_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_ATLAS_DIR = REPO / "data/processed/v13_dct_reference/mouse_kidney_atlas"
DEFAULT_FROZEN_TIERS = (
    REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
)
DEFAULT_PREREG = REPO / "config/clinical_axes_podocyte_sensitivities.yaml"
DEFAULT_RESULTS = REPO / "data/results/atlas_loso_markers"

REFERENCE_LABEL = "FULL_ATLAS_REFERENCE"
MINIMUM_GENES = al.MINIMUM_SET_GENES

SUMMARY_COLUMNS = [
    "omitted_origin",
    "omit_index",
    "seed",
    "podocyte_contributing",
    "status",
    "status_reason",
    "qc_pass",
    "qc_failed_compartments",
    "n_source_studies_podocyte",
    "effective_podocyte_hs_rule",
    "n_sets_defined",
    "n_sets_ge8_genes",
    "n_sets_missing_from_rebuild",
    "n_full_podocyte_hs",
    "n_podocyte_hs",
    "n_retained_from_full",
    "jaccard_vs_full",
    "n_observable_all_missions",
    "estimate",
    "ci_low_mkh",
    "ci_high_mkh",
    "i_squared",
    "max_t_fwer",
    "rank_by_abs_t",
    "n_family_sets",
    "reference_estimate",
    "retention_ratio",
    "direction_retained",
    "k0_p1_coef_score_units",
    "k0_p1_ci_low",
    "k0_p1_ci_high",
    "k0_p1_p_mkh",
    "k0_p1_positive",
    "k0_status",
    "gate_included",
    "gate_rule",
    "gate_pass",
]
BOOLEAN_COLUMNS = [
    "qc_pass",
    "direction_retained",
    "k0_p1_positive",
    "gate_included",
    "gate_pass",
]


# --- scoring pieces (flight data; only reached without --skip-scoring) --------------


def observable_by_mission(missions: Mapping[str, Any], threshold: float):
    """CPM-eligible genes per mission and the subset that is also in the expression matrix."""

    eligible = {
        mission: set(cpm_eligible_genes(data, threshold))
        for mission, data in missions.items()
    }
    observable = {
        mission: eligible[mission] & set(data.expression.index.astype(str))
        for mission, data in missions.items()
    }
    return eligible, observable


def evaluate_rebuild(
    tier_path: Path,
    missions: Mapping[str, Any],
    *,
    threshold: float,
    n_permutations: int,
    seed: int,
    chunk_size: int,
    eligible: Mapping[str, set[str]] | None = None,
    observable: Mapping[str, set[str]] | None = None,
) -> dict[str, Any]:
    """Score one tier table exactly as ``run_compartment_context`` does and extract the podocyte row.

    Evaluability follows the frozen rule (a set needs >= 8 observable genes in every
    mission); the observable genes are CPM-eligible AND present in the expression
    matrix, which can only differ from the compartment-context rule where that script
    would otherwise abort.
    """

    family, audit = load_family(Path(tier_path))
    if al.PODOCYTE_SET not in family:
        n_defined = int(
            audit.loc[audit["gene_set"] == al.PODOCYTE_SET, "n_defined"].sum()
        )
        return {
            "status": al.STATUS_NOT_EVALUABLE,
            "reason": (
                f"{al.PODOCYTE_SET} defined with {n_defined} genes (< {MINIMUM_GENES})"
            ),
        }
    if eligible is None or observable is None:
        eligible, observable = observable_by_mission(missions, threshold)
    retained = {}
    for gene_set, spec in family.items():
        genes = set(spec["subdomains"]["atlas_markers"]["genes"])
        if min(len(genes & observable[mission]) for mission in missions) >= MINIMUM_GENES:
            retained[gene_set] = spec
    if al.PODOCYTE_SET not in retained:
        return {
            "status": al.STATUS_NOT_EVALUABLE,
            "reason": (
                f"{al.PODOCYTE_SET} has < {MINIMUM_GENES} observable genes in at "
                "least one primary mission"
            ),
        }
    scores, design, _, _ = combined_score_design(
        missions, retained, cpm_threshold=threshold
    )
    result = blocked_meta_permutation(
        scores,
        design,
        n_permutations=int(n_permutations),
        seed=int(seed),
        chunk_size=int(chunk_size),
    )
    meta = result.observed_meta.reset_index()
    meta["abs_t_mkh"] = meta["t_mkh"].abs()
    meta["rank_by_abs_t"] = (
        meta["abs_t_mkh"].rank(ascending=False, method="min").astype(int)
    )
    lookup = audit.set_index("gene_set")
    meta["report_compartment"] = meta["axis"].map(lookup["report_compartment"])
    meta["tier"] = meta["axis"].map(lookup["tier"])
    meta["n_defined_genes"] = meta["axis"].map(lookup["n_defined"])
    meta["n_permutations"] = int(n_permutations)
    meta["seed"] = int(seed)
    row = meta.loc[meta["axis"] == al.PODOCYTE_SET].iloc[0]
    podocyte_genes = list(
        retained[al.PODOCYTE_SET]["subdomains"]["atlas_markers"]["genes"]
    )
    common = set.intersection(*(observable[mission] for mission in missions))
    return {
        "status": al.STATUS_EVALUABLE,
        "reason": "",
        "meta": meta,
        "mission_effects": result.mission_effects,
        "podocyte_row": row,
        "n_family_sets": int(len(retained)),
        "n_observable_all_missions": int(len(set(podocyte_genes) & common)),
    }


def k0_p1(
    missions: Mapping[str, Any],
    podocyte_genes: list[str],
    structural_genes: list[str],
    eligible: Mapping[str, set[str]],
) -> dict[str, Any]:
    """K0 model P1 for a (rebuilt) podocyte set: HC3 per mission, REML/mKH across missions.

    The structural comparator is the scaffold control minus the podocyte genes.  Same
    scoring and regression helpers as ``run_podocyte_scaffold_specificity.py``; the
    flight coefficient is in score units (``coef_score_units``), not a standardized g.
    """

    pod_set = set(podocyte_genes)
    structural_disjoint = [gene for gene in structural_genes if gene not in pod_set]
    estimates: list[float] = []
    variances: list[float] = []
    for mission, data in missions.items():
        index = set(data.expression.index.astype(str))
        pod_used = [g for g in podocyte_genes if g in eligible[mission] and g in index]
        struc_used = [
            g for g in structural_disjoint if g in eligible[mission] and g in index
        ]
        if len(pod_used) < MINIMUM_GENES or len(struc_used) < MINIMUM_GENES:
            return {
                "k0_status": (
                    f"NOT_EVALUABLE: {mission} podocyte={len(pod_used)}, "
                    f"structural={len(struc_used)} (< {MINIMUM_GENES})"
                )
            }
        podocyte = score_signed_axis(data.expression, {g: 1 for g in pod_used}).scores
        structural = score_signed_axis(
            data.expression, {g: 1 for g in struc_used}
        ).scores
        hc3 = next(
            row
            for row in _regression_rows(
                mission, data.metadata, podocyte, structural, "P1_podocyte_adjusted"
            )
            if row["variance_type"] == "HC3"
        )
        estimates.append(hc3["estimate"])
        variances.append(hc3["variance"])
    fit = random_effects_reml_mkh(np.array(estimates), np.array(variances))
    return {
        "k0_status": "ok",
        "k0_p1_coef_score_units": float(fit.estimate),
        "k0_p1_ci_low": float(fit.ci_low),
        "k0_p1_ci_high": float(fit.ci_high),
        "k0_p1_p_mkh": float(fit.p),
        "k0_p1_positive": bool(fit.estimate > 0),
    }


# --- helpers -----------------------------------------------------------------------


def _safe_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", str(name))


def apply_gene_map_override(config: dict, gene_map: Path | None) -> dict:
    """Point ``config["gene_mapping"]["path"]`` at an alternative ID map, in memory only.

    Never edits a config file (same semantics as the other clinical-axes entry scripts).
    """

    if gene_map is not None:
        resolved = Path(gene_map).expanduser().resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"--gene-map {resolved} does not exist")
        config["gene_mapping"]["path"] = str(resolved)
    return config


def _join_reasons(*reasons: str) -> str:
    return "; ".join(reason for reason in reasons if reason)


def _file_record(path: Path) -> dict[str, Any]:
    path = Path(path)
    try:
        shown = str(path.resolve().relative_to(REPO))
    except ValueError:
        shown = str(path)
    if not path.is_file():
        return {"path": shown, "sha256": None}
    return {"path": shown, "sha256": al.sha256_file(path)}


def reference_crosscheck(
    meta_path: Path, estimate: float, max_t_fwer: float, *, tolerance: float = 1e-6
) -> dict[str, Any]:
    """Compare the in-process full-atlas podocyte row with a ``compartment_context`` run.

    The in-process reference (frozen tiers, same missions and gene map, seed parent+100)
    is what the retention ratios use; a reference-run ``compartment_context_meta_results.tsv``
    only cross-checks it.  A mismatch is recorded, not fatal (it signals a different
    gene map or tier file).
    """

    meta_path = Path(meta_path)
    table = pd.read_csv(meta_path, sep="\t")
    row = table.loc[table["axis"] == al.PODOCYTE_SET]
    if row.empty:
        return {**_file_record(meta_path), "estimate_matches": False, "note": "podocyte row absent"}
    row = row.iloc[0]
    diff = abs(float(row["estimate"]) - estimate)
    out: dict[str, Any] = {
        **_file_record(meta_path),
        "estimate_in_file": float(row["estimate"]),
        "estimate_in_process": estimate,
        "abs_diff_estimate": diff,
        "estimate_matches": bool(diff <= tolerance),
    }
    if "max_t_fwer" in row.index and pd.notna(row["max_t_fwer"]):
        out["max_t_fwer_in_file"] = float(row["max_t_fwer"])
        out["max_t_fwer_in_process"] = max_t_fwer
        out["abs_diff_max_t_fwer"] = abs(float(row["max_t_fwer"]) - max_t_fwer)
    return out


def baseline_contract(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
    structural_terms: Mapping[str, Any],
    audit_cfg: Mapping[str, Any],
    qc_cfg: Mapping[str, Any],
    frozen_tiers: pd.DataFrame,
    expression_path: Path | None,
    *,
    min_cpm: float,
) -> dict[str, Any]:
    """The no-exclusion rebuild must reproduce the frozen expression table and tiers."""

    rebuild = al.rebuild_tiers(
        counts,
        sample_meta,
        structural_terms,
        audit_cfg,
        qc_cfg,
        omitted=None,
        min_cpm=min_cpm,
        keep_expression=True,
    )
    contract: dict[str, Any] = {"rebuild_status": rebuild.status}
    if rebuild.tiers is None:
        contract["tiers_identical_to_frozen"] = False
        contract["detail"] = rebuild.reason
        return contract
    contract["tiers_comparison"] = al.tiers_equal_setwise(rebuild.tiers, frozen_tiers)
    contract["tiers_identical_to_frozen"] = bool(
        contract["tiers_comparison"]["identical"]
    )
    contract["expression_max_abs_diff"] = None
    if expression_path is not None and Path(expression_path).is_file():
        stored = pd.read_csv(expression_path, sep="\t", float_precision="round_trip")
        rebuilt = rebuild.expression
        same_layout = (
            len(stored) == len(rebuilt)
            and (stored["gene_symbol"].to_numpy() == rebuilt["gene_symbol"].to_numpy()).all()
            and (stored["compartment"].to_numpy() == rebuilt["compartment"].to_numpy()).all()
        )
        contract["expression_same_layout"] = bool(same_layout)
        if same_layout:
            contract["expression_max_abs_diff"] = float(
                max(
                    np.abs(stored[col].to_numpy() - rebuilt[col].to_numpy()).max()
                    for col in (
                        "mean_cpm",
                        "median_cpm",
                        "source_study_detection_fraction",
                        "n_source_studies",
                    )
                )
            )
    return contract


def _format_table(table: pd.DataFrame) -> str:
    return table.to_string(index=False, float_format=lambda v: f"{v:.3f}")


# --- main run ------------------------------------------------------------------------


def run(args: argparse.Namespace, *, missions: Mapping[str, Any] | None = None) -> dict[str, Any]:
    results = Path(args.results)
    skip = bool(args.skip_scoring)

    audit_cfg = load_config(args.audit_config)
    freeze_cfg = yaml.safe_load(Path(args.freeze_config).read_text())
    builder = freeze_cfg["reference_builder"]
    qc_cfg = builder["atlas_qc"]
    min_cpm = float(builder["broad_expression"]["min_cpm"])
    marker_cfg = audit_cfg["marker_tiers"]
    hs_fraction = float(marker_cfg["high_specificity_target_detection_fraction_min"])
    analysis_cfg = yaml.safe_load(Path(args.analysis_config).read_text())
    parent_seed = int(analysis_cfg["seed"])
    threshold = float(analysis_cfg["eligibility"]["cpm_threshold"])

    # Preregistration guard: only the scoring path reads flight data.
    if skip:
        prereg_record: dict[str, Any] = {
            "applicable": False,
            "reason": "--skip-scoring: tiers only, no mission data read",
        }
        retention_min = al.RETENTION_RATIO_MIN
    else:
        prereg_path = Path(args.prereg)
        if not prereg_path.is_absolute():
            prereg_path = REPO / prereg_path
        prereg_record = al.check_preregistration(
            prereg_path, repo=REPO, allow_uncommitted=args.allow_uncommitted_prereg
        )
        prereg_record["applicable"] = True
        prereg_yaml = yaml.safe_load(prereg_path.read_text()) or {}
        retention_min = float(prereg_yaml.get("retention_ratio_min", al.RETENTION_RATIO_MIN))
        prereg_record["analysis_id"] = prereg_yaml.get("analysis_id")
        prereg_record["retention_ratio_min"] = retention_min

    results.mkdir(parents=True, exist_ok=True)
    counts, sample_meta = al.load_pseudobulk(args.pseudobulk_counts, args.pseudobulk_meta)
    structural_path = Path(audit_cfg["input"]["structural_scaffold_reference"])
    if not structural_path.is_absolute():
        structural_path = REPO / structural_path
    structural_terms = json.loads(structural_path.read_text())
    frozen_tiers = pd.read_csv(args.frozen_tiers, sep="\t")
    full_sets = al.set_membership(frozen_tiers)

    contract = baseline_contract(
        counts,
        sample_meta,
        structural_terms,
        audit_cfg,
        qc_cfg,
        frozen_tiers,
        args.expression,
        min_cpm=min_cpm,
    )
    if not contract["tiers_identical_to_frozen"] and not args.allow_baseline_mismatch:
        raise RuntimeError(
            "the no-exclusion rebuild does not reproduce the frozen tiers set-for-set "
            f"({contract}); item 8 would not be a leave-out of the frozen definition. "
            "Pass --allow-baseline-mismatch to override."
        )

    rebuilds = al.leave_one_study_out_tiers(
        counts,
        sample_meta,
        structural_terms,
        audit_cfg,
        qc_cfg,
        origins=args.origins,
        min_cpm=min_cpm,
    )

    summary_rows: list[dict[str, Any]] = []
    set_rows: list[pd.DataFrame] = []
    comparison_rows: list[pd.DataFrame] = []
    tier_paths: dict[str, Path] = {}
    for rebuild in rebuilds:
        origin = str(rebuild.omitted_origin)
        row = al.tier_summary_row(rebuild, full_sets, hs_detection_fraction=hs_fraction)
        row["seed"] = al.loso_seed(parent_seed, int(rebuild.index))
        if rebuild.tiers is not None:
            path = results / f"atlas_loso_tiers_{_safe_name(origin)}.tsv"
            rebuild.tiers.to_csv(path, sep="\t", index=False)
            tier_paths[origin] = path
            rebuilt_sets = al.set_membership(rebuild.tiers)
            comparison_rows.append(al.compare_sets(full_sets, rebuilt_sets, origin))
            members = rebuild.tiers[
                ["gene_set", "report_compartment", "tier", "gene_symbol"]
            ].copy()
            members.insert(0, "omitted_origin", origin)
            members["in_full_set"] = [
                gene in full_sets.get(name, frozenset())
                for name, gene in zip(members["gene_set"], members["gene_symbol"])
            ]
            set_rows.append(members)
        summary_rows.append(row)

    family_meta: list[pd.DataFrame] = []
    reference: dict[str, Any] = {}
    gene_map_record: dict[str, Any] = {"used": False}
    if skip:
        for row in summary_rows:
            row["status"] = (
                al.STATUS_NOT_EVALUABLE
                if not row["tiers_evaluable"]
                else al.STATUS_NOT_SCORED
            )
            row["status_reason"] = row["tier_status_reason"]
            if row["status"] == al.STATUS_NOT_SCORED:
                row["status_reason"] = _join_reasons(
                    row["tier_status_reason"], "--skip-scoring: tiers only"
                )
    else:
        config = apply_gene_map_override(deepcopy(analysis_cfg), args.gene_map)
        if missions is None:
            missions = load_primary_missions(config, REPO)
        map_path = Path(config["gene_mapping"]["path"])
        if not map_path.is_absolute():
            map_path = REPO / map_path
        gene_map_record = {
            "used": True,
            "override": args.gene_map is not None,
            **_file_record(map_path),
        }
        eligible, observable = observable_by_mission(missions, threshold)
        ref = evaluate_rebuild(
            args.frozen_tiers,
            missions,
            threshold=threshold,
            n_permutations=args.permutations,
            seed=parent_seed + al.SEED_OFFSET_REFERENCE,
            chunk_size=args.chunk_size,
            eligible=eligible,
            observable=observable,
        )
        if ref["status"] != al.STATUS_EVALUABLE:
            raise RuntimeError(f"full-atlas reference is not evaluable: {ref['reason']}")
        ref_row = ref["podocyte_row"]
        reference_estimate = float(ref_row["estimate"])
        ref_pod = _frozen_genes(frozen_tiers, al.PODOCYTE_SET)
        ref_struct = _frozen_genes(frozen_tiers, al.STRUCTURAL_SET)
        ref_k0 = k0_p1(missions, ref_pod, ref_struct, eligible)
        reference = {
            "seed": parent_seed + al.SEED_OFFSET_REFERENCE,
            "estimate": reference_estimate,
            "ci_low_mkh": float(ref_row["ci_low_mkh"]),
            "ci_high_mkh": float(ref_row["ci_high_mkh"]),
            "i_squared": float(ref_row["i_squared"]),
            "max_t_fwer": float(ref_row["max_t_fwer"]),
            "rank_by_abs_t": int(ref_row["rank_by_abs_t"]),
            "n_family_sets": ref["n_family_sets"],
            "n_observable_all_missions": ref["n_observable_all_missions"],
            "k0_p1": ref_k0,
        }
        if args.reference_meta is not None:
            reference["crosscheck"] = reference_crosscheck(
                args.reference_meta, reference_estimate, float(ref_row["max_t_fwer"])
            )
            if not reference["crosscheck"]["estimate_matches"]:
                print(
                    "WARNING: in-process full-atlas podocyte estimate differs from "
                    f"{args.reference_meta} (check the gene map / tiers)",
                    file=sys.stderr,
                )
        ref_meta = ref["meta"].copy()
        ref_meta.insert(0, "omitted_origin", REFERENCE_LABEL)
        family_meta.append(ref_meta)

        for row, rebuild in zip(summary_rows, rebuilds):
            origin = str(rebuild.omitted_origin)
            row["reference_estimate"] = reference_estimate
            if not row["tiers_evaluable"]:
                row["status"] = al.STATUS_NOT_EVALUABLE
                row["status_reason"] = row["tier_status_reason"]
                continue
            scored = evaluate_rebuild(
                tier_paths[origin],
                missions,
                threshold=threshold,
                n_permutations=args.permutations,
                seed=int(row["seed"]),
                chunk_size=args.chunk_size,
                eligible=eligible,
                observable=observable,
            )
            row["status"] = scored["status"]
            row["status_reason"] = _join_reasons(
                row["tier_status_reason"], scored["reason"]
            )
            if scored["status"] != al.STATUS_EVALUABLE:
                continue
            prow = scored["podocyte_row"]
            row.update(
                {
                    "n_observable_all_missions": scored["n_observable_all_missions"],
                    "estimate": float(prow["estimate"]),
                    "ci_low_mkh": float(prow["ci_low_mkh"]),
                    "ci_high_mkh": float(prow["ci_high_mkh"]),
                    "i_squared": float(prow["i_squared"]),
                    "max_t_fwer": float(prow["max_t_fwer"]),
                    "rank_by_abs_t": int(prow["rank_by_abs_t"]),
                    "n_family_sets": scored["n_family_sets"],
                }
            )
            row.update(
                al.retention_gate(
                    row["estimate"], reference_estimate, retention_min=retention_min
                )
            )
            rebuilt_pod = _frozen_genes(rebuild.tiers, al.PODOCYTE_SET)
            rebuilt_struct = _frozen_genes(rebuild.tiers, al.STRUCTURAL_SET)
            row.update(k0_p1(missions, rebuilt_pod, rebuilt_struct, eligible))
            meta = scored["meta"].copy()
            meta.insert(0, "omitted_origin", origin)
            family_meta.append(meta)

    gate_rule = al.GATE_RULE_TEXT.format(minimum=retention_min)
    for row in summary_rows:
        included = bool(
            row["qc_pass"]
            and row.get("status") == al.STATUS_EVALUABLE
            and row.get("n_podocyte_hs", 0) >= MINIMUM_GENES
        )
        row["gate_included"] = included
        row["gate_rule"] = gate_rule
        if not included:
            row["gate_pass"] = pd.NA

    summary = pd.DataFrame(summary_rows).reindex(columns=SUMMARY_COLUMNS)
    for column in BOOLEAN_COLUMNS:
        summary[column] = summary[column].astype("boolean")

    summary.to_csv(results / "atlas_loso_summary.tsv", sep="\t", index=False)
    if set_rows:
        pd.concat(set_rows, ignore_index=True).to_csv(
            results / "atlas_loso_sets.tsv.gz", sep="\t", index=False, compression="gzip"
        )
    if comparison_rows:
        pd.concat(comparison_rows, ignore_index=True).to_csv(
            results / "atlas_loso_set_comparison.tsv", sep="\t", index=False
        )
    if family_meta:
        pd.concat(family_meta, ignore_index=True).to_csv(
            results / "atlas_loso_family_meta.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

    manifest = build_manifest(
        args=args,
        skip=skip,
        summary=summary,
        contract=contract,
        reference=reference,
        prereg_record=prereg_record,
        gene_map_record=gene_map_record,
        parent_seed=parent_seed,
        min_cpm=min_cpm,
        gate_rule=gate_rule,
        retention_min=retention_min,
        marker_cfg=marker_cfg,
        qc_cfg=qc_cfg,
        sample_meta=sample_meta,
    )
    (results / "atlas_loso_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=_json_default) + "\n"
    )
    show = summary[
        [
            c
            for c in (
                "omitted_origin",
                "podocyte_contributing",
                "status",
                "qc_pass",
                "n_source_studies_podocyte",
                "n_podocyte_hs",
                "n_retained_from_full",
                "jaccard_vs_full",
                "estimate",
                "retention_ratio",
                "max_t_fwer",
                "rank_by_abs_t",
                "k0_p1_coef_score_units",
                "gate_pass",
            )
            if c in summary.columns
        ]
    ]
    print("=== item 8: leave-one-atlas-source-study-out ===")
    print(_format_table(show))
    return {"summary": summary, "manifest": manifest}


def _json_default(value: Any):
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if np.isnan(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if value is pd.NA:
        return None
    if isinstance(value, datetime.date):  # YAML parses unquoted dates (lock_date)
        return value.isoformat()
    raise TypeError(f"not JSON serializable: {type(value)}")


def _group_gate(summary: pd.DataFrame) -> dict[str, Any]:
    included = summary[summary["gate_included"].fillna(False).astype(bool)]
    if included.empty:
        return {"n_included": 0, "n_pass": 0, "all_pass": None}
    passed = included["gate_pass"].fillna(False).astype(bool)
    return {
        "n_included": int(len(included)),
        "n_pass": int(passed.sum()),
        "all_pass": bool(passed.all()),
        "failed_origins": list(included.loc[~passed, "omitted_origin"]),
        "min_retention_ratio": float(included["retention_ratio"].min()),
        "max_retention_ratio": float(included["retention_ratio"].max()),
    }


def build_manifest(
    *,
    args: argparse.Namespace,
    skip: bool,
    summary: pd.DataFrame,
    contract: Mapping[str, Any],
    reference: Mapping[str, Any],
    prereg_record: Mapping[str, Any],
    gene_map_record: Mapping[str, Any],
    parent_seed: int,
    min_cpm: float,
    gate_rule: str,
    retention_min: float,
    marker_cfg: Mapping[str, Any],
    qc_cfg: Mapping[str, Any],
    sample_meta: pd.DataFrame,
) -> dict[str, Any]:
    n_by_status = summary["status"].value_counts().to_dict()
    manifest: dict[str, Any] = {
        "analysis_id": ANALYSIS_ID,
        "item": "8_atlas_loso",
        "git_head": al.git_head(REPO),
        "skip_scoring": skip,
        "blinding": (
            "tiers only; no flight/ground data read" if skip else
            "scored after the preregistration was committed (see preregistration)"
        ),
        "preregistration": dict(prereg_record),
        "inputs": {
            "pseudobulk_counts": _file_record(Path(args.pseudobulk_counts)),
            "pseudobulk_sample_metadata": _file_record(Path(args.pseudobulk_meta)),
            "frozen_tiers": _file_record(Path(args.frozen_tiers)),
            "atlas_compartment_expression": (
                None if args.expression is None else _file_record(Path(args.expression))
            ),
            "audit_config": _file_record(Path(args.audit_config)),
            "freeze_config": _file_record(Path(args.freeze_config)),
            "analysis_config": _file_record(Path(args.analysis_config)),
        },
        "gene_map": dict(gene_map_record),
        "baseline_contract": dict(contract),
        "thresholds": {
            "status": "unchanged_frozen",
            "marker_tiers": dict(marker_cfg),
            "atlas_qc": dict(qc_cfg),
            "min_cpm_detection": min_cpm,
        },
        "effective_stringency_disclosure": (
            "Dropping a source study lowers the number of studies contributing to each "
            "compartment, so fractional frozen rules become effectively stricter: the "
            "high-specificity detection fraction >= 0.75 is 3-of-4 for the full podocyte "
            "compartment (4 contributing studies) and 3-of-3 once a podocyte-contributing "
            "study is omitted, and those rebuilds sit at the QC minimum of 3 studies. "
            "Omitting a non-contributing study changes only the max-other denominators "
            "(and other compartments' study counts). Thresholds are not retuned; results "
            "are reported separately by podocyte_contributing."
        ),
        "seeds": {
            "parent_seed": parent_seed,
            "reference_offset": al.SEED_OFFSET_REFERENCE,
            "loso_offset_base": al.SEED_OFFSET_LOSO,
            "loso_seed_by_origin": {
                str(o): al.loso_seed(parent_seed, i)
                for i, o in enumerate(sorted(sample_meta["source_study_id"].astype(str).unique()))
            },
            "n_permutations": None if skip else int(args.permutations),
            "chunk_size": None if skip else int(args.chunk_size),
        },
        "gate": {
            "rule": gate_rule,
            "retention_ratio_min": retention_min,
            "applies_to": "QC-passing rebuilds with >= 8 podocyte high-specificity genes (EVALUABLE)",
            "not_evaluable_is_not_a_failure": True,
            "n_by_status": {str(k): int(v) for k, v in n_by_status.items()},
        },
        "reference_full_atlas": dict(reference),
    }
    if not skip:
        by_group = {}
        for label in ("yes", "no"):
            sub = summary[summary["podocyte_contributing"] == label]
            by_group[f"podocyte_contributing_{label}"] = _group_gate(sub)
        manifest["gate"]["overall"] = _group_gate(summary)
        manifest["gate"]["by_podocyte_contributing"] = by_group
        k0_nonpositive = summary.loc[
            summary["k0_p1_positive"].notna() & (~summary["k0_p1_positive"].fillna(True).astype(bool)),
            "omitted_origin",
        ]
        manifest["k0_rule"] = {
            "statement": (
                "No rebuild can upgrade K0's verdict. If P1 <= 0 for an omitted study, report "
                "'podocyte-leaning tail not robust to omitting <study>'."
            ),
            "omitted_origins_with_nonpositive_p1": list(k0_nonpositive),
        }
    manifest["outputs"] = sorted(
        [
            "atlas_loso_summary.tsv",
            "atlas_loso_sets.tsv.gz",
            "atlas_loso_set_comparison.tsv",
            "atlas_loso_manifest.json",
            "atlas_loso_tiers_<origin>.tsv",
        ]
        + ([] if skip else ["atlas_loso_family_meta.tsv.gz"])
    )
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--audit-config", type=Path, default=DEFAULT_AUDIT_CONFIG)
    parser.add_argument("--freeze-config", type=Path, default=DEFAULT_FREEZE_CONFIG)
    parser.add_argument("--analysis-config", type=Path, default=DEFAULT_ANALYSIS_CONFIG)
    parser.add_argument(
        "--pseudobulk-counts",
        type=Path,
        default=DEFAULT_ATLAS_DIR / "atlas_pseudobulk_counts.tsv.gz",
    )
    parser.add_argument(
        "--pseudobulk-meta",
        type=Path,
        default=DEFAULT_ATLAS_DIR / "atlas_pseudobulk_sample_metadata.tsv",
    )
    parser.add_argument(
        "--expression",
        type=Path,
        default=DEFAULT_ATLAS_DIR / "atlas_compartment_expression.tsv.gz",
        help="Stored compartment expression table used only for the baseline contract check.",
    )
    parser.add_argument("--frozen-tiers", type=Path, default=DEFAULT_FROZEN_TIERS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--prereg", type=Path, default=DEFAULT_PREREG)
    parser.add_argument(
        "--allow-uncommitted-prereg",
        action="store_true",
        help="Tests only: run although the preregistration YAML has uncommitted changes.",
    )
    parser.add_argument(
        "--skip-scoring",
        action="store_true",
        help=(
            "Only rebuild the tiers and report flight-blind set statistics; reads no "
            "mission data and needs no committed preregistration."
        ),
    )
    parser.add_argument(
        "--gene-map",
        type=Path,
        default=None,
        help=(
            "Override gene_mapping.path (Ensembl -> symbol map) in memory for the mission "
            "loaders; never edits a config file."
        ),
    )
    parser.add_argument(
        "--reference-meta",
        type=Path,
        default=None,
        help=(
            "Optional compartment_context_meta_results.tsv of the reference run; the "
            "full-atlas podocyte row is cross-checked against it in the manifest."
        ),
    )
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument(
        "--origins",
        nargs="+",
        default=None,
        help="Restrict the leave-out loop to these source studies (default: all).",
    )
    parser.add_argument(
        "--allow-baseline-mismatch",
        action="store_true",
        help="Proceed although the no-exclusion rebuild differs from the frozen tiers.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        run(args)
    except al.PreregistrationError as exc:
        print(f"REFUSED: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
