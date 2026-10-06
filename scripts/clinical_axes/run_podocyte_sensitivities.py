#!/usr/bin/env python3
"""G1 sensitivities of the frozen 157-gene podocyte high-specificity set (items 1-7).

Targets: the podocyte set (gated), the structural scaffold control and its podocyte-disjoint
part (descriptive comparators), and the K0 models P0, P1, S0, S1, D (observed HC3 coefficients
pooled by REML/mKH, in score units, never called g). Every gate rule, threshold, label and seed
offset is read from the preregistration YAML at run time, and the script refuses to run unless
that YAML is committed (``--allow-uncommitted-prereg`` exists for tests only).

Items: 1 signed median, 2 common-gene intersection, 3 leave-one-mission, 4 gene influence
(leave-one-gene, exact signed contributions, leave-top-k), 5 OSD-462 preparation, 6 OSD-253
white-light rerun control, 7 OSD-163 mapping-rate adjustment. Outputs go to
``{results}/podocyte_sensitivities/``. Sensitivities are not new discovery families.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import re
import subprocess
import sys
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.clinical_axes.run_compartment_context import load_family  # noqa: E402
from scripts.clinical_axes.run_podocyte_scaffold_specificity import (  # noqa: E402
    PODOCYTE_SET,
    STRUCTURAL_SET,
    _frozen_genes,
)
from scripts.clinical_axes.run_sensitivities import (  # noqa: E402
    leave_one_gene,
    secondary_qc_covariate_sensitivity,
)
from src.clinical_axes.analysis import combined_score_design  # noqa: E402
from src.clinical_axes.data import (  # noqa: E402
    MissionData,
    cpm_eligible_genes,
    load_osd253_control_sensitivity,
    load_osd462_preparation,
    load_primary_missions,
)
from src.clinical_axes.sensitivity import (  # noqa: E402
    DEFAULT_TOP_KS,
    K0_MODEL_CODES,
    K0_MODEL_SPECS,
    ITEM4_GATE_KEYS,
    ITEM4_SECONDARY_KEYS,
    GateRule,
    PreregistrationError,
    ancova_standardized_effect,
    annotate_leave_one_gene,
    ci_excludes_zero,
    classify_overall,
    combine_passes,
    common_observable_family,
    common_observable_genes,
    conservative_verdict,
    direction_retained,
    estimate_verdict,
    evaluate_gate,
    evaluate_gene_influence,
    family_genes,
    genes_used_by_mission,
    k0_observed,
    leave_top_k_curve,
    parse_gate_rule,
    residualized_mission_effect,
    retention_ratio,
    signed_contributions,
    single_set_family,
    standardized_qc_covariate,
    variant_meta,
    verdict_label,
)
from src.clinical_axes.statistics import (  # noqa: E402
    blocked_meta_permutation,
    leave_one_mission_out,
    random_effects_reml_mkh,
)


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
DEFAULT_PREREG = REPO / "config/clinical_axes_podocyte_sensitivities.yaml"
DEFAULT_RESULTS = REPO / "data/results/clinical-axes/podocyte_sensitivities_run"
OUTPUT_SUBDIR = "podocyte_sensitivities"
STRUCTURAL_DISJOINT = "structural_disjoint"
QC_MISSION = "OSD-163"
QC_METRIC = "uniquely_mapped_percent"
CONSISTENCY_TOLERANCE = 1e-6

ITEM_SLUGS = {
    1: "signed_median",
    2: "common_intersection",
    3: "leave_one_mission",
    4: "gene_influence",
    5: "osd462_preparation",
    6: "osd253_rerun_control",
    7: "osd163_mapping_rate",
}
SUPPORTED_OSD462_PREPARATIONS = ("mRNA", "UPX")
SUPPORTED_OSD253_VARIANTS = ("rerun_white",)
SUPPORTED_ITEM7_PRIMARY = ("label_blind_residualization",)
SUPPORTED_ITEM7_SECONDARY = ("ancova_descriptive",)
DEFAULT_K0_NOTE = "podocyte-leaning tail not robust to <variant>"

SUMMARY_COLUMNS = [
    "item",
    "variant",
    "target",
    "scoring_method",
    "estimate",
    "ci_low_mkh",
    "ci_high_mkh",
    "p_mkh",
    "i_squared",
    "tau2",
    "k_missions",
    "n_positive_missions",
    "reference_estimate",
    "retention_ratio",
    "direction_retained",
    "ci_excludes_zero",
    "gate_rule",
    "gate_pass",
    "notes",
]
MISSION_COLUMNS = ["item", "variant", "target", "mission", "estimate", "variance", "n_genes_used"]
LOMO_COLUMNS = [
    "target",
    "omitted_mission",
    "estimate",
    "ci_low_mkh",
    "ci_high_mkh",
    "p_mkh",
    "i_squared",
    "direction_retained",
    "retention_ratio",
    "lomo_max_t_fwer",
]
K0_COLUMNS = [
    "item",
    "variant",
    "model",
    "scoring_method",
    "coef_score_units",
    "standard_error_mkh",
    "ci_low_mkh",
    "ci_high_mkh",
    "p_mkh",
    "i_squared",
    "tau2",
    "k_missions",
    "n_positive_missions",
    "reference_coef_score_units",
    "direction_retained",
    "ci_excludes_zero",
    "notes",
]


# --------------------------------------------------------------------------------------
# Preregistration guard and parsing
# --------------------------------------------------------------------------------------


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args], cwd=repo, capture_output=True, text=True, check=False
    )


def preregistration_guard(
    path: Path, *, repo: Path = REPO, allow_uncommitted: bool = False
) -> dict[str, Any]:
    """Refuse unless the preregistration exists, is tracked and has no uncommitted changes.

    Read-only git calls only. ``allow_uncommitted`` (tests) records the dirty state instead of
    refusing; a missing file is always refused because its rules cannot be read.
    """

    path = Path(path)
    resolved = path if path.is_absolute() else Path(repo) / path
    if not resolved.is_file():
        raise PreregistrationError(f"preregistration {path} does not exist; refusing to run")
    record: dict[str, Any] = {
        "path": str(path),
        "sha256": sha256(resolved),
        "allow_uncommitted_prereg": bool(allow_uncommitted),
    }
    try:
        head = _git(repo, "rev-parse", "HEAD")
        status = _git(repo, "status", "--porcelain", "--", str(resolved))
        tracked = _git(repo, "ls-files", "--error-unmatch", "--", str(resolved))
        last = _git(repo, "log", "-1", "--format=%H", "--", str(resolved))
    except FileNotFoundError as error:  # git not installed
        record.update({"git_head": None, "committed": False, "git_error": str(error)})
    else:
        record["git_head"] = head.stdout.strip() if head.returncode == 0 else None
        record["git_status_porcelain"] = (
            status.stdout.strip()
            if status.returncode == 0
            else f"git status failed: {status.stderr.strip()}"
        )
        record["tracked"] = tracked.returncode == 0
        record["prereg_last_commit"] = (
            last.stdout.strip() or None if last.returncode == 0 else None
        )
        record["committed"] = bool(
            head.returncode == 0
            and status.returncode == 0
            and tracked.returncode == 0
            and not status.stdout.strip()
        )
    if not record["committed"] and not allow_uncommitted:
        detail = record.get("git_status_porcelain") or record.get("git_error") or "untracked"
        raise PreregistrationError(
            f"preregistration {path} is not committed ({detail}); refusing to run. "
            "Commit it first (tests may pass --allow-uncommitted-prereg)."
        )
    return record


@dataclass(frozen=True)
class PodocytePrereg:
    """The parts of the preregistration this driver evaluates, validated."""

    analysis_id: str
    primary_target: str
    comparators: tuple[str, ...]
    k0_models: tuple[str, ...]
    retention_ratio_min: float
    item_keys: dict[int, str]
    rules: dict[int, GateRule]
    gene_influence_gates: dict[str, Any]
    gene_influence_secondary: dict[str, Any]
    osd462_variants: tuple[str, ...]
    osd253_variant: str
    item7_primary: str
    item7_secondary: str | None
    required_items: tuple[int, ...]
    robust_label: str
    fragile_label: str
    seed_offsets: dict[str, Any]
    seed_base_from_parent: bool
    k0_rule_text: str
    k0_note_template: str
    k0_rule_model: str
    raw: dict[str, Any] = field(repr=False)


def _need(mapping: Mapping[str, Any], key: str, where: str) -> Any:
    if not isinstance(mapping, Mapping) or key not in mapping:
        raise PreregistrationError(f"preregistration lacks {where}.{key}")
    return mapping[key]


def _parse_overall(overall: Mapping[str, Any]) -> tuple[str, str, tuple[int, ...]]:
    robust = fragile = None
    for label, text in overall.items():
        normalized = " ".join(str(text).split())
        match = re.fullmatch(
            r"items\s+(\d+)\s*-\s*(\d+)\s+all\s+pass", normalized, re.IGNORECASE
        )
        if match:
            robust = (str(label), int(match.group(1)), int(match.group(2)))
        match = re.match(
            r"any\s+of\s+items\s+(\d+)\s*-\s*(\d+)\s+fails", normalized, re.IGNORECASE
        )
        if match:
            fragile = (str(label), int(match.group(1)), int(match.group(2)))
    if robust is None or fragile is None:
        raise PreregistrationError(
            "overall must define 'items A-B all pass' and 'any of items A-B fails'"
        )
    if robust[1:] != fragile[1:]:
        raise PreregistrationError("robust and fragile item ranges disagree")
    return robust[0], fragile[0], tuple(range(robust[1], robust[2] + 1))


def parse_preregistration(raw: Mapping[str, Any]) -> PodocytePrereg:
    """Validate and parse the G1 preregistration; anything unrecognised is refused."""

    if not isinstance(raw, Mapping):
        raise PreregistrationError("preregistration must be a YAML mapping")
    targets = _need(raw, "targets", "")
    primary = str(_need(targets, "primary", "targets"))
    comparators = tuple(str(name) for name in targets.get("descriptive_comparators") or [])
    codes = targets.get("k0_models") or list(K0_MODEL_CODES)
    try:
        k0_models = tuple(K0_MODEL_CODES[str(code)] for code in codes)
    except KeyError as error:
        raise PreregistrationError(f"unknown K0 model code {error}") from None
    retention_min = float(_need(raw, "retention_ratio_min", ""))
    default_text = str(_need(raw, "gate_rule_default", ""))
    aliases = {"gate_rule_default": default_text}
    default_rule = parse_gate_rule(default_text)
    items = _need(raw, "items", "")
    item_keys: dict[int, str] = {}
    for key in items:
        match = re.match(r"(\d+)_", str(key))
        if match:
            number = int(match.group(1))
            if number in item_keys:
                raise PreregistrationError(f"duplicate item number {number}")
            item_keys[number] = str(key)
    robust_label, fragile_label, required = _parse_overall(_need(raw, "overall", ""))
    for number in required:
        if number not in ITEM_SLUGS:
            raise PreregistrationError(f"item {number} is required but not run by this driver")
        if number not in item_keys:
            raise PreregistrationError(f"required item {number} is missing from items")

    rules: dict[int, GateRule] = {}
    for number in (1, 2, 3, 5, 6, 7):
        spec = items[item_keys[number]]
        rules[number] = parse_gate_rule(
            _need(spec, "gate", f"items.{item_keys[number]}"), aliases=aliases
        )
    for number, rule in {0: default_rule, **rules}.items():
        for clause in rule.clauses:
            if clause.kind == "ratio" and abs(clause.threshold - retention_min) > 1e-12:
                raise PreregistrationError(
                    f"item {number} ratio threshold {clause.threshold} differs from "
                    f"retention_ratio_min {retention_min}"
                )
    if not any(clause.kind == "all" for clause in rules[3].clauses):
        raise PreregistrationError("item 3 gate must be an 'all N LOMO estimates' rule")

    spec4 = items[item_keys[4]]
    gates4 = dict(_need(spec4, "gates", f"items.{item_keys[4]}"))
    secondary4 = dict(spec4.get("secondary_bounds_v13_precedent") or {})
    if set(gates4) - set(ITEM4_GATE_KEYS) or not gates4:
        raise PreregistrationError(f"item 4 gates must be among {ITEM4_GATE_KEYS}")
    if set(secondary4) - set(ITEM4_SECONDARY_KEYS):
        raise PreregistrationError(f"item 4 secondary bounds must be among {ITEM4_SECONDARY_KEYS}")

    spec5 = items[item_keys[5]]
    variants5 = tuple(str(v) for v in _need(spec5, "variants", f"items.{item_keys[5]}"))
    if not variants5 or set(variants5) - set(SUPPORTED_OSD462_PREPARATIONS):
        raise PreregistrationError(f"item 5 variants must be among {SUPPORTED_OSD462_PREPARATIONS}")
    variant6 = str(_need(items[item_keys[6]], "variant", f"items.{item_keys[6]}"))
    if variant6 not in SUPPORTED_OSD253_VARIANTS:
        raise PreregistrationError(f"item 6 variant must be one of {SUPPORTED_OSD253_VARIANTS}")
    spec7 = items[item_keys[7]]
    primary7 = str(_need(spec7, "primary_method", f"items.{item_keys[7]}"))
    if primary7 not in SUPPORTED_ITEM7_PRIMARY:
        raise PreregistrationError(f"item 7 primary_method must be one of {SUPPORTED_ITEM7_PRIMARY}")
    secondary7 = spec7.get("secondary")
    if secondary7 is not None and str(secondary7) not in SUPPORTED_ITEM7_SECONDARY:
        raise PreregistrationError(f"item 7 secondary must be one of {SUPPORTED_ITEM7_SECONDARY}")

    k0_rule = " ".join(str(raw.get("k0_rule", "")).split())
    template = re.search(r'"([^"]*<variant>[^"]*)"', k0_rule)
    model = re.search(r"If\s+(P0|P1|S0|S1|D)\s*<=\s*0", k0_rule)
    return PodocytePrereg(
        analysis_id=str(raw.get("analysis_id", "")),
        primary_target=primary,
        comparators=comparators,
        k0_models=k0_models,
        retention_ratio_min=retention_min,
        item_keys=item_keys,
        rules=rules,
        gene_influence_gates=gates4,
        gene_influence_secondary=secondary4,
        osd462_variants=variants5,
        osd253_variant=variant6,
        item7_primary=primary7,
        item7_secondary=None if secondary7 is None else str(secondary7),
        required_items=required,
        robust_label=robust_label,
        fragile_label=fragile_label,
        seed_offsets=dict(raw.get("seed_offsets") or {}),
        seed_base_from_parent=bool(raw.get("seed_base_from_parent", True)),
        k0_rule_text=k0_rule,
        k0_note_template=template.group(1) if template else DEFAULT_K0_NOTE,
        k0_rule_model=K0_MODEL_CODES[model.group(1)] if model else K0_MODEL_CODES["P1"],
        raw=dict(raw),
    )


def load_preregistration(
    path: Path, *, repo: Path = REPO, allow_uncommitted: bool = False
) -> tuple[PodocytePrereg, dict[str, Any]]:
    """Guard first (refuse a missing or dirty file), then parse."""

    record = preregistration_guard(path, repo=repo, allow_uncommitted=allow_uncommitted)
    resolved = Path(path) if Path(path).is_absolute() else Path(repo) / path
    raw = yaml.safe_load(resolved.read_text())
    prereg = parse_preregistration(raw)
    record["analysis_id"] = prereg.analysis_id
    record["lock_date"] = (raw or {}).get("lock_date")
    return prereg, record


def config_with_gene_map(config: Mapping[str, Any], gene_map: Path | None) -> dict[str, Any]:
    """Copy of the frozen config with ``gene_mapping.path`` replaced in memory only."""

    out = deepcopy(dict(config))
    if gene_map is not None:
        resolved = Path(gene_map).expanduser().resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"gene-ID map {resolved} does not exist")
        out["gene_mapping"] = dict(out["gene_mapping"])
        out["gene_mapping"]["path"] = str(resolved)
    return out


def gene_map_record(config: Mapping[str, Any], label: str) -> dict[str, Any]:
    path = Path(str(config["gene_mapping"]["path"]))
    resolved = path if path.is_absolute() else REPO / path
    return {
        "label": label,
        "path": str(path),
        "sha256": sha256(resolved) if resolved.is_file() else None,
    }


# --------------------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------------------


@dataclass
class PipelineInputs:
    missions: dict[str, MissionData]
    alternates: dict[str, tuple[str, Callable[[], MissionData]]]
    tiers_path: Path
    threshold: float
    seed_base: int


@dataclass
class PipelineOptions:
    lomo_permutations: int = 0
    chunk_size: int = 256
    k0_leave_one_gene: bool = False
    structural_gene_influence: bool = True
    median_results: Path | None = None
    reference_run: Path | None = None
    top_ks: tuple[int | str, ...] = DEFAULT_TOP_KS


@dataclass
class PipelineResult:
    tables: dict[str, pd.DataFrame]
    gates: dict[str, Any]
    info: dict[str, Any]


def build_inputs(config: Mapping[str, Any], tiers_path: Path, prereg: PodocytePrereg) -> PipelineInputs:
    """Real-data inputs: five primary missions plus lazy loaders for items 5 and 6."""

    missions = load_primary_missions(config, REPO)
    alternates: dict[str, tuple[str, Callable[[], MissionData]]] = {}
    for preparation in prereg.osd462_variants:
        alternates[f"OSD462_{preparation}"] = (
            "OSD-462",
            lambda preparation=preparation: load_osd462_preparation(config, REPO, preparation),
        )
    alternates[f"OSD253_{prereg.osd253_variant}"] = (
        "OSD-253",
        lambda: load_osd253_control_sensitivity(config, REPO),
    )
    return PipelineInputs(
        missions=missions,
        alternates=alternates,
        tiers_path=Path(tiers_path),
        threshold=float(config["eligibility"]["cpm_threshold"]),
        seed_base=int(config["seed"]),
    )


def resolve_targets(
    tiers_path: Path, prereg: PodocytePrereg
) -> tuple[dict[str, list[str]], list[str], list[str], dict[str, Any]]:
    """Target gene lists (frozen definitions) and the full frozen set family."""

    family_all, _ = load_family(Path(tiers_path))
    tiers = pd.read_csv(tiers_path, sep="\t")
    if prereg.primary_target != PODOCYTE_SET:
        raise PreregistrationError(
            f"primary target {prereg.primary_target!r} is not K0's {PODOCYTE_SET!r}"
        )
    podocyte = _frozen_genes(tiers, PODOCYTE_SET)
    structural = _frozen_genes(tiers, STRUCTURAL_SET)
    disjoint = [gene for gene in structural if gene not in set(podocyte)]
    targets: dict[str, list[str]] = {PODOCYTE_SET: podocyte}
    known = set(tiers["gene_set"].astype(str))
    for name in prereg.comparators:
        if name == STRUCTURAL_DISJOINT:
            targets[name] = disjoint
        elif name in known:
            targets[name] = _frozen_genes(tiers, name)
        else:
            raise PreregistrationError(f"unknown comparator target {name!r}")
    for name, genes in targets.items():
        if name in family_all and family_genes(family_all, name) != genes:
            raise RuntimeError(f"{name}: load_family and _frozen_genes disagree")
    return targets, podocyte, disjoint, family_all


def evaluable_family(
    family_all: Mapping[str, Any], missions: Mapping[str, MissionData], threshold: float
) -> dict[str, Any]:
    """The compartment-context evaluability rule: >= 8 CPM-eligible genes in every mission."""

    eligible = {m: set(cpm_eligible_genes(d, threshold)) for m, d in missions.items()}
    return {
        name: spec
        for name, spec in family_all.items()
        if min(
            len(set(spec["subdomains"]["atlas_markers"]["genes"]) & genes)
            for genes in eligible.values()
        )
        >= 8
    }


class _Collector:
    """Accumulates rows for the summary, mission-effect and K0 tables."""

    def __init__(self, prereg: PodocytePrereg, mission_order: Sequence[str]):
        self.prereg = prereg
        self.order = list(mission_order)
        self.summary: list[dict[str, Any]] = []
        self.missions: list[dict[str, Any]] = []
        self.k0: list[dict[str, Any]] = []
        self.k0_missions: list[pd.DataFrame] = []
        self.coverage: list[pd.DataFrame] = []
        self.reference: dict[str, float] = {}
        self.k0_reference: dict[str, float] = {}
        self.k0_notes: list[str] = []

    def meta_row(
        self,
        item: str,
        variant: str,
        target: str,
        method: str,
        row: Mapping[str, Any] | None,
        *,
        gate: dict[str, Any] | None = None,
        notes: str = "",
    ) -> dict[str, Any]:
        row = {} if row is None else dict(row)
        estimate = float(row.get("estimate", np.nan))
        low = float(row.get("ci_low_mkh", np.nan))
        high = float(row.get("ci_high_mkh", np.nan))
        reference = self.reference.get(target, np.nan)
        out = {
            "item": item,
            "variant": variant,
            "target": target,
            "scoring_method": method,
            "estimate": estimate,
            "ci_low_mkh": low,
            "ci_high_mkh": high,
            "p_mkh": float(row.get("p_mkh", np.nan)),
            "i_squared": float(row.get("i_squared", np.nan)),
            "tau2": float(row.get("tau2", np.nan)),
            "k_missions": row.get("k_missions", np.nan),
            "n_positive_missions": row.get("n_positive_missions", np.nan),
            "reference_estimate": reference,
            "retention_ratio": retention_ratio(estimate, reference),
            "direction_retained": direction_retained(estimate, reference),
            "ci_excludes_zero": ci_excludes_zero(low, high),
            "gate_rule": "" if gate is None else gate["rule"],
            "gate_pass": None if gate is None else gate["passed"],
            "notes": notes,
        }
        self.summary.append(out)
        return out

    def mission_effects(
        self,
        item: str,
        variant: str,
        target: str,
        effects: pd.DataFrame,
        coverage: pd.DataFrame | None,
    ) -> None:
        used = {}
        if coverage is not None and len(coverage):
            sub = coverage[coverage["axis"] == target]
            used = sub.groupby("mission", sort=False)["n_used"].sum().to_dict()
        for _, row in effects.iterrows():
            self.missions.append(
                {
                    "item": item,
                    "variant": variant,
                    "target": target,
                    "mission": row["mission"],
                    "estimate": float(row["estimate"]),
                    "variance": float(row["variance"]),
                    "n_genes_used": used.get(row["mission"], np.nan),
                }
            )
        if coverage is not None and len(coverage):
            self.coverage.append(
                coverage[coverage["axis"] == target].assign(item=item, variant=variant)
            )

    def k0_rows(self, item: str, variant: str, result, *, notes: str = "") -> None:
        meta = result.meta.set_index("model")
        for model in self.prereg.k0_models:
            row = meta.loc[model]
            coef = float(row["coef_score_units"])
            reference = self.k0_reference.get(model, np.nan)
            note = notes
            if model == self.prereg.k0_rule_model and coef <= 0:
                message = self.prereg.k0_note_template.replace("<variant>", f"{item}/{variant}")
                note = "; ".join(part for part in (note, message) if part)
                self.k0_notes.append(message)
            self.k0.append(
                {
                    "item": item,
                    "variant": variant,
                    "model": model,
                    "scoring_method": row["scoring_method"],
                    "coef_score_units": coef,
                    "standard_error_mkh": float(row["standard_error_mkh"]),
                    "ci_low_mkh": float(row["ci_low_mkh"]),
                    "ci_high_mkh": float(row["ci_high_mkh"]),
                    "p_mkh": float(row["p_mkh"]),
                    "i_squared": float(row["i_squared"]),
                    "tau2": float(row["tau2"]),
                    "k_missions": int(row["k_missions"]),
                    "n_positive_missions": int(row["n_positive_missions"]),
                    "reference_coef_score_units": reference,
                    "direction_retained": direction_retained(coef, reference),
                    "ci_excludes_zero": ci_excludes_zero(
                        float(row["ci_low_mkh"]), float(row["ci_high_mkh"])
                    ),
                    "notes": note,
                }
            )
        self.k0_missions.append(result.mission_effects.assign(item=item, variant=variant))

    def k0_failure(self, item: str, variant: str, reason: str) -> None:
        for model in self.prereg.k0_models:
            self.k0.append(
                {
                    "item": item,
                    "variant": variant,
                    "model": model,
                    "coef_score_units": np.nan,
                    "reference_coef_score_units": self.k0_reference.get(model, np.nan),
                    "notes": f"not evaluable: {reason}",
                }
            )


_VARIANT_ERRORS = (ValueError, RuntimeError, KeyError, FileNotFoundError)


def _combine_item(results: Sequence[dict[str, Any]]) -> dict[str, Any]:
    passed = combine_passes([result["passed"] for result in results])
    return {"passed": passed, "verdict": verdict_label(passed)}


def run_pipeline(
    inputs: PipelineInputs, prereg: PodocytePrereg, options: PipelineOptions
) -> PipelineResult:
    """Steps 1-9 on in-memory missions; no file is read except the tiers table and,
    optionally, the median and reference-run tables for consistency checks."""

    missions = dict(inputs.missions)
    order = list(missions)
    threshold = inputs.threshold
    targets, podocyte_genes, structural_disjoint, family_all = resolve_targets(
        inputs.tiers_path, prereg
    )
    primary = prereg.primary_target
    families = {name: single_set_family(genes, name) for name, genes in targets.items()}
    collect = _Collector(prereg, order)
    info: dict[str, Any] = {
        "missions": order,
        "targets": {name: len(genes) for name, genes in targets.items()},
        "structural_overlap_removed": sorted(
            set(_frozen_genes(pd.read_csv(inputs.tiers_path, sep="\t"), STRUCTURAL_SET))
            & set(podocyte_genes)
        ),
        "checks": {},
    }
    items: dict[int, dict[str, Any]] = {}
    item_key = prereg.item_keys

    # ---- Step 1: reference (mean scoring, five primary missions) --------------------------
    reference_meta: dict[str, pd.DataFrame] = {}
    reference_effects: dict[str, pd.DataFrame] = {}
    reference_details: dict[str, Any] = {}
    reference_coverage: dict[str, pd.DataFrame] = {}
    used_by_target: dict[str, dict[str, list[str]]] = {}
    for name, family in families.items():
        result = variant_meta(missions, family, threshold=threshold)
        _, _, _, details = combined_score_design(missions, family, cpm_threshold=threshold)
        reference_meta[name] = result.meta
        reference_effects[name] = result.mission_effects
        reference_details[name] = details
        reference_coverage[name] = result.coverage
        used_by_target[name] = genes_used_by_mission(result.coverage, name)
        collect.reference[name] = float(result.meta.iloc[0]["estimate"])
        collect.meta_row("0_reference", "primary_mean", name, "mean", result.meta.iloc[0])
        collect.mission_effects("0_reference", "primary_mean", name, result.mission_effects, result.coverage)
    k0_reference = k0_observed(missions, podocyte_genes, structural_disjoint, threshold=threshold)
    collect.k0_reference = (
        k0_reference.meta.set_index("model")["coef_score_units"].astype(float).to_dict()
    )
    collect.k0_rows("0_reference", "primary_mean", k0_reference)
    info["checks"]["reference_run"] = _check_reference_run(options.reference_run, reference_meta)

    # ---- Step 2 / item 2: common-gene intersection ----------------------------------------
    common = common_observable_genes(missions, threshold)
    info["n_common_observable_genes"] = {
        name: len(set(genes) & common) for name, genes in targets.items()
    }
    for name, family in families.items():
        gate = None
        try:
            common_family = common_observable_family(missions, family, threshold)
            result = variant_meta(missions, common_family, threshold=threshold)
        except ValueError as error:
            row, notes, result = None, f"not evaluable: {error}", None
        else:
            row, notes = result.meta.iloc[0], f"n_common_genes={len(family_genes(common_family, name))}"
        if name == primary:
            gate = evaluate_gate(
                prereg.rules[2],
                estimate=np.nan if row is None else float(row["estimate"]),
                reference=collect.reference[name],
            )
            items[2] = {**gate, "statistic": None if row is None else float(row["estimate"])}
        collect.meta_row(item_key[2], "common", name, "mean", row, gate=gate, notes=notes)
        if result is not None:
            collect.mission_effects(item_key[2], "common", name, result.mission_effects, result.coverage)
    try:
        k0_common = k0_observed(
            missions,
            [g for g in podocyte_genes if g in common],
            [g for g in structural_disjoint if g in common],
            threshold=threshold,
        )
    except _VARIANT_ERRORS as error:
        collect.k0_failure(item_key[2], "common", str(error))
    else:
        collect.k0_rows(item_key[2], "common", k0_common)

    # ---- Step 3 / item 3: leave one mission out -------------------------------------------
    lomo_fwer, lomo_info = _lomo_family_permutations(
        missions, family_all, targets, threshold, inputs.seed_base, prereg, options
    )
    info["lomo_family_permutation"] = lomo_info
    lomo_rows = []
    rule3 = prereg.rules[3]
    all_clause = next(clause for clause in rule3.clauses if clause.kind == "all")
    for name in families:
        effects = reference_effects[name].set_index("mission").loc[order]
        table = leave_one_mission_out(
            effects["estimate"].to_numpy(), effects["variance"].to_numpy(), order
        )
        for _, row in table.iterrows():
            omitted = row["omitted_mission"]
            fwer = lomo_fwer.get((name, omitted), np.nan)
            reference = collect.reference[name]
            n_positive = int((effects.drop(index=omitted)["estimate"] > 0).sum())
            row = pd.Series({**row.to_dict(), "n_positive_missions": n_positive})
            lomo_rows.append(
                {
                    "target": name,
                    "omitted_mission": omitted,
                    "estimate": float(row["estimate"]),
                    "ci_low_mkh": float(row["ci_low_mkh"]),
                    "ci_high_mkh": float(row["ci_high_mkh"]),
                    "p_mkh": float(row["p_mkh"]),
                    "i_squared": float(row["i_squared"]),
                    "tau2": float(row["tau2"]),
                    "k_missions": int(row["k_missions"]),
                    "n_positive_missions": n_positive,
                    "direction_retained": direction_retained(row["estimate"], reference),
                    "retention_ratio": retention_ratio(row["estimate"], reference),
                    "lomo_max_t_fwer": fwer,
                }
            )
            row_gate = None
            if name == primary:
                passed = bool(
                    _compare(all_clause.op, float(row["estimate"]), all_clause.threshold)
                )
                row_gate = {"rule": rule3.text, "passed": passed}
            collect.meta_row(
                item_key[3],
                f"omit_{omitted}",
                name,
                "mean",
                row,
                gate=row_gate,
                notes="" if np.isnan(fwer) else f"lomo_max_t_fwer={fwer:.6g} (descriptive)",
            )
        if name == primary:
            gate = evaluate_gate(
                rule3,
                lomo_estimates=table["estimate"].tolist(),
                reference=collect.reference[name],
            )
            items[3] = {**gate, "statistic": float(table["estimate"].min())}
    leave_one_mission_table = pd.DataFrame(lomo_rows)
    for key, value in lomo_info.get("observed_estimates", {}).items():
        name, omitted = key.split("|", 1)
        ours = leave_one_mission_table[
            (leave_one_mission_table["target"] == name)
            & (leave_one_mission_table["omitted_mission"] == omitted)
        ]["estimate"].iloc[0]
        if abs(float(ours) - value) > 1e-8:
            raise RuntimeError(f"{name}/{omitted}: LOMO family estimate disagrees ({ours} vs {value})")

    # ---- Step 4 / item 4: gene influence ---------------------------------------------------
    influence_targets = [primary] + (
        [name for name in families if name != primary]
        if options.structural_gene_influence
        else []
    )
    logo_frames = []
    contribution_frames, concentration_frames, strata_frames, top_k_frames = [], [], [], []
    top_k_summary: dict[str, dict[str, Any]] = {}
    for name in influence_targets:
        logo_frames.append(
            leave_one_gene(missions, families[name], threshold, reference_meta[name])
        )
        contributions = signed_contributions(
            missions, families[name], reference_details[name]
        )
        contribution_frames.append(contributions.genes)
        concentration_frames.append(contributions.concentration)
        strata_frames.append(contributions.strata)
        curve = leave_top_k_curve(
            missions,
            targets[name],
            contributions.genes.set_index("gene")["signed_contribution"],
            threshold=threshold,
            ks=options.top_ks,
            name=name,
        )
        top_k_frames.append(curve.curve)
        top_k_summary[name] = {"k50": curve.k50, "k_ci0": curve.k_ci0}
    logo = annotate_leave_one_gene(
        pd.concat(logo_frames, ignore_index=True), collect.reference, used_by_target
    )
    concentration = pd.concat(concentration_frames, ignore_index=True)
    concentration["k50"] = concentration["target"].map(
        {name: value["k50"] for name, value in top_k_summary.items()}
    )
    concentration["k_ci0"] = concentration["target"].map(
        {name: value["k_ci0"] for name, value in top_k_summary.items()}
    )
    counted = logo[(logo["target"] == primary) & logo["counted_in_gate"]]
    primary_concentration = concentration[concentration["target"] == primary].iloc[0].to_dict()
    item4 = evaluate_gene_influence(
        prereg.gene_influence_gates,
        prereg.gene_influence_secondary,
        logo=counted,
        concentration=primary_concentration,
    )
    item4["descriptive"] = {
        "k50": top_k_summary[primary]["k50"],
        "k_ci0": top_k_summary[primary]["k_ci0"],
        "n_eff_positive": primary_concentration.get("n_eff_positive"),
        "max_absolute_share": primary_concentration.get("max_absolute_share"),
        "n_genes_counted": int(len(counted)),
        "n_genes_never_eligible": int(
            ((logo["target"] == primary) & ~logo["gene_used_in_any_mission"]).sum()
        ),
    }
    item4["statistic"] = float(counted["retention_ratio"].min()) if len(counted) else np.nan
    items[4] = item4
    if len(counted):
        worst = counted.loc[counted["retention_ratio"].idxmin()]
        collect.meta_row(
            item_key[4],
            "leave_one_gene_min_retention",
            primary,
            "mean",
            {**worst.to_dict(), "k_missions": len(order)},
            gate={"rule": "item-4 gates (see podocyte_sensitivity_gates.json)", "passed": item4["passed"]},
            notes=f"omitted_gene={worst['omitted_gene']}; "
            f"max_signed_share={primary_concentration.get('max_signed_share'):.6g} "
            f"({primary_concentration.get('max_signed_share_gene')})",
        )
    k0_logo = _k0_leave_one_gene(
        missions, podocyte_genes, structural_disjoint, threshold, used_by_target[primary], collect
    ) if options.k0_leave_one_gene else pd.DataFrame()

    # ---- Steps 5-6 / items 5-6: alternate OSD-462 preparations and OSD-253 control ----------
    alternate_plan = [(5, f"OSD462_{prep}") for prep in prereg.osd462_variants] + [
        (6, f"OSD253_{prereg.osd253_variant}")
    ]
    alternate_gates: dict[int, list[dict[str, Any]]] = {5: [], 6: []}
    for number, variant in alternate_plan:
        replaced, loader = inputs.alternates.get(variant, (None, None))
        alternate = None
        error_text = ""
        if loader is None:
            error_text = f"no loader for {variant}"
        else:
            try:
                alternate = loader()
            except _VARIANT_ERRORS as error:
                error_text = f"{type(error).__name__}: {error}"
        missions_variant = dict(missions)
        if alternate is not None:
            missions_variant[replaced] = alternate
        for name, family in families.items():
            row, result, notes = None, None, ""
            if alternate is None:
                notes = f"not evaluable: {error_text}"
            else:
                try:
                    result = variant_meta(missions_variant, family, threshold=threshold)
                except _VARIANT_ERRORS as error:
                    notes = f"not evaluable: {error}"
                else:
                    row = result.meta.iloc[0]
                    effect = result.mission_effects.set_index("mission").loc[replaced]
                    notes = f"{replaced} mission effect {float(effect['estimate']):.6g}"
            gate = None
            if name == primary:
                gate = evaluate_gate(
                    prereg.rules[number],
                    estimate=np.nan if row is None else float(row["estimate"]),
                    reference=collect.reference[name],
                )
                alternate_gates[number].append({"variant": variant, **gate})
            collect.meta_row(item_key[number], variant, name, "mean", row, gate=gate, notes=notes)
            if result is not None:
                collect.mission_effects(
                    item_key[number], variant, name, result.mission_effects, result.coverage
                )
        if alternate is None:
            collect.k0_failure(item_key[number], variant, error_text)
        else:
            try:
                k0_variant = k0_observed(
                    missions_variant, podocyte_genes, structural_disjoint, threshold=threshold
                )
            except _VARIANT_ERRORS as error:
                collect.k0_failure(item_key[number], variant, str(error))
            else:
                collect.k0_rows(item_key[number], variant, k0_variant)
    for number in (5, 6):
        combined = _combine_item(alternate_gates[number])
        ratios = [
            clause["value"]
            for gate in alternate_gates[number]
            for clause in gate["clauses"]
            if clause["kind"] == "ratio" and clause["value"] is not None
        ]
        items[number] = {
            **combined,
            "rule": prereg.rules[number].text,
            "variants": alternate_gates[number],
            "statistic": float(min(ratios)) if ratios else np.nan,
        }

    # ---- Step 7 / item 7: OSD-163 mapping rate ---------------------------------------------
    items[7] = _item7(
        missions,
        families,
        reference_effects,
        reference_details,
        reference_coverage,
        prereg,
        collect,
        threshold,
    )
    try:
        covariate = standardized_qc_covariate(missions[QC_MISSION], QC_METRIC)
        k0_qc = k0_observed(
            missions,
            podocyte_genes,
            structural_disjoint,
            threshold=threshold,
            extra_covariates={QC_MISSION: pd.DataFrame({f"z_{QC_METRIC}": covariate})},
        )
    except _VARIANT_ERRORS as error:
        collect.k0_failure(item_key[7], "osd163_mapping_rate_covariate", str(error))
    else:
        collect.k0_rows(item_key[7], "osd163_mapping_rate_covariate", k0_qc)

    # ---- Step 8 / item 1: signed median ----------------------------------------------------
    median_file = _median_reference(options.median_results)
    info["checks"]["median_family_file"] = median_file["status"]
    for name, family in families.items():
        gate, notes = None, ""
        try:
            result = variant_meta(missions, family, threshold=threshold, method="median")
        except _VARIANT_ERRORS as error:
            row, result, notes = None, None, f"not evaluable: {error}"
        else:
            row = result.meta.iloc[0]
            notes = _median_file_note(median_file, name, float(row["estimate"]), info)
        if name == primary:
            gate = evaluate_gate(
                prereg.rules[1],
                estimate=np.nan if row is None else float(row["estimate"]),
                reference=collect.reference[name],
            )
            items[1] = {**gate, "statistic": None if row is None else float(row["estimate"])}
        collect.meta_row(item_key[1], "median", name, "median", row, gate=gate, notes=notes)
        if result is not None:
            collect.mission_effects(item_key[1], "median", name, result.mission_effects, result.coverage)
    try:
        k0_median = k0_observed(
            missions, podocyte_genes, structural_disjoint, threshold=threshold, method="median"
        )
    except _VARIANT_ERRORS as error:
        collect.k0_failure(item_key[1], "median", str(error))
    else:
        collect.k0_rows(item_key[1], "median", k0_median)

    # ---- Step 9: gates ----------------------------------------------------------------------
    verdicts = {number: items[number]["verdict"] for number in items}
    overall = classify_overall(
        verdicts,
        prereg.required_items,
        robust_label=prereg.robust_label,
        fragile_label=prereg.fragile_label,
    )
    overall["failed_item_keys"] = [item_key[n] for n in overall["failed_items"]]
    overall["not_evaluable_item_keys"] = [item_key[n] for n in overall["not_evaluable_items"]]
    reference_p1 = k0_reference.meta.set_index("model").loc[prereg.k0_rule_model]
    gates = {
        "analysis_id": prereg.analysis_id,
        "primary_target": primary,
        "reference_estimate": collect.reference[primary],
        "retention_ratio_min": prereg.retention_ratio_min,
        "items": {
            item_key[number]: {"item": number, **_jsonable(items[number])}
            for number in sorted(items)
        },
        "overall": overall,
        "k0": {
            "rule": prereg.k0_rule_text,
            "model": prereg.k0_rule_model,
            "reference_coef_score_units": float(reference_p1["coef_score_units"]),
            "reference_ci_low_mkh": float(reference_p1["ci_low_mkh"]),
            "reference_ci_high_mkh": float(reference_p1["ci_high_mkh"]),
            "reference_ci_crosses_zero": bool(
                reference_p1["ci_low_mkh"] <= 0.0 <= reference_p1["ci_high_mkh"]
            ),
            "verdict_upgrade_permitted": False,
            "nonpositive_variant_notes": list(dict.fromkeys(collect.k0_notes)),
        },
        "descriptive_comparators_not_gated": [name for name in families if name != primary],
        "items_not_run_by_this_driver": sorted(
            key for number, key in prereg.item_keys.items() if number not in ITEM_SLUGS
        ),
        "multiplicity": prereg.raw.get("multiplicity", ""),
    }
    tables = {
        "summary": _ordered(_by_item(pd.DataFrame(collect.summary)), SUMMARY_COLUMNS),
        "mission_effects": _ordered(pd.DataFrame(collect.missions), MISSION_COLUMNS),
        "leave_one_mission": _ordered(leave_one_mission_table, LOMO_COLUMNS),
        "leave_one_gene": logo,
        "gene_contributions": pd.concat(contribution_frames, ignore_index=True),
        "contributor_concentration": concentration,
        "contribution_strata_audit": pd.concat(strata_frames, ignore_index=True),
        "leave_top_k": pd.concat(top_k_frames, ignore_index=True),
        "k0": _ordered(_by_item(pd.DataFrame(collect.k0)), K0_COLUMNS),
        "k0_mission_effects": pd.concat(collect.k0_missions, ignore_index=True),
        "gene_coverage": pd.concat(collect.coverage, ignore_index=True)
        if collect.coverage
        else pd.DataFrame(),
    }
    if len(k0_logo):
        tables["k0_leave_one_gene"] = k0_logo
    return PipelineResult(tables=tables, gates=gates, info=info)


def _by_item(frame: pd.DataFrame) -> pd.DataFrame:
    """Stable sort by the leading item number (0_reference first)."""

    number = frame["item"].astype(str).str.extract(r"^(\d+)_", expand=False).astype(int)
    return frame.iloc[np.argsort(number.to_numpy(), kind="stable")].reset_index(drop=True)


def _compare(op: str, value: float, threshold: float) -> bool:
    return {">": value > threshold, ">=": value >= threshold, "<": value < threshold, "<=": value <= threshold}[op]


def _ordered(frame: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    frame = frame.copy()
    for column in columns:
        if column not in frame.columns:
            frame[column] = np.nan
    extra = [column for column in frame.columns if column not in columns]
    return frame[list(columns) + extra]


def _check_reference_run(
    reference_run: Path | None, reference_meta: Mapping[str, pd.DataFrame]
) -> dict[str, Any]:
    """Refuse to mix baselines: an existing reference run must reproduce step 1 exactly."""

    if reference_run is None:
        return {"status": "not_requested"}
    path = Path(reference_run) / "compartment_context_meta_results.tsv"
    if not path.is_file():
        return {"status": "absent", "path": str(path)}
    table = pd.read_csv(path, sep="\t").set_index("axis")
    compared = {}
    for name, meta in reference_meta.items():
        if name not in table.index:
            continue
        ours = meta.iloc[0]
        delta = max(
            abs(float(ours[column]) - float(table.loc[name, column]))
            for column in ("estimate", "ci_low_mkh", "ci_high_mkh")
        )
        compared[name] = delta
        if delta > CONSISTENCY_TOLERANCE:
            raise RuntimeError(
                f"{name}: step-1 reference differs from {path} by {delta:.3g}; "
                "inputs (gene-ID map, tiers or data) differ from the reference run"
            )
    return {"status": "consistent", "path": str(path), "max_abs_delta": compared}


def _lomo_family_permutations(
    missions: Mapping[str, MissionData],
    family_all: Mapping[str, Any],
    targets: Mapping[str, Sequence[str]],
    threshold: float,
    seed_base: int,
    prereg: PodocytePrereg,
    options: PipelineOptions,
) -> tuple[dict[tuple[str, str], float], dict[str, Any]]:
    """Descriptive LOMO family-wise p over the fixed evaluable family on each 4-mission subset."""

    if options.lomo_permutations <= 0:
        return {}, {"status": "not_requested"}
    offsets = list(prereg.seed_offsets.get("lomo_family") or [])
    order = list(missions)
    if len(offsets) != len(order):
        raise PreregistrationError(
            f"seed_offsets.lomo_family needs {len(order)} offsets, found {len(offsets)}"
        )
    family = evaluable_family(family_all, missions, threshold)
    scores, design, _, _ = combined_score_design(missions, family, cpm_threshold=threshold)
    out: dict[tuple[str, str], float] = {}
    seeds = {}
    estimates = {}
    for omitted, offset in zip(order, offsets):
        keep = design["mission"] != omitted
        seed = int(seed_base) + int(offset)
        seeds[omitted] = seed
        result = blocked_meta_permutation(
            scores.loc[keep],
            design.loc[keep],
            n_permutations=int(options.lomo_permutations),
            seed=seed,
            chunk_size=int(options.chunk_size),
        )
        for name in targets:
            if name in result.max_t_fwer.index:
                out[(name, omitted)] = float(result.max_t_fwer[name])
                estimates[f"{name}|{omitted}"] = float(result.observed_meta.loc[name, "estimate"])
    return out, {
        "status": "computed",
        "n_permutations": int(options.lomo_permutations),
        "n_family_sets": len(family),
        "seeds_by_omitted_mission": seeds,
        "observed_estimates": estimates,
        "interpretation": "descriptive; not a new discovery family",
    }


def _k0_leave_one_gene(
    missions, podocyte_genes, structural_disjoint, threshold, used_by_mission, collect
) -> pd.DataFrame:
    """Descriptive: P1 after omitting each podocyte gene used in any mission."""

    used = {gene for genes in used_by_mission.values() for gene in genes}
    model = collect.prereg.k0_rule_model
    spec = [entry for entry in K0_MODEL_SPECS if entry[0] == model]
    reference = collect.k0_reference.get(model, np.nan)
    rows = []
    for gene in podocyte_genes:
        if gene not in used:
            continue
        try:
            result = k0_observed(
                missions,
                [g for g in podocyte_genes if g != gene],
                structural_disjoint,
                threshold=threshold,
                models=spec,
            )
        except _VARIANT_ERRORS as error:
            rows.append({"omitted_gene": gene, "model": model, "evaluable": False, "reason": str(error)})
            continue
        row = result.meta.iloc[0]
        rows.append(
            {
                "omitted_gene": gene,
                "model": model,
                "evaluable": True,
                "coef_score_units": float(row["coef_score_units"]),
                "ci_low_mkh": float(row["ci_low_mkh"]),
                "ci_high_mkh": float(row["ci_high_mkh"]),
                "reference_coef_score_units": reference,
                "direction_retained": direction_retained(row["coef_score_units"], reference),
                "retention_ratio": retention_ratio(row["coef_score_units"], reference),
                "status": "descriptive",
            }
        )
    return pd.DataFrame(rows)


def _item7(
    missions,
    families,
    reference_effects,
    reference_details,
    reference_coverage,
    prereg,
    collect,
    threshold,
):
    """OSD-163 mapping rate: label-blind residualization (gated) and ANCOVA (descriptive)."""

    primary = prereg.primary_target
    key = prereg.item_keys[7]
    gate_result = None
    failure = ""
    try:
        covariate = standardized_qc_covariate(missions[QC_MISSION], QC_METRIC)
    except _VARIANT_ERRORS as error:
        covariate = None
        failure = str(error)
    metadata = missions[QC_MISSION].metadata
    for name, family in families.items():
        gate = None
        row, notes, effects = None, "", None
        if covariate is not None:
            try:
                contract = secondary_qc_covariate_sensitivity(missions, family, threshold).iloc[0]
                score = reference_details[name][QC_MISSION][name].scores
                adjusted = residualized_mission_effect(score, metadata, covariate)
                if abs(adjusted["estimate"] - float(contract["osd163_adjusted_estimate"])) > 1e-10:
                    raise RuntimeError("residualized OSD-163 effect disagrees with the frozen contract")
                effects = reference_effects[name].copy()
                mask = effects["mission"] == QC_MISSION
                effects.loc[mask, "estimate"] = adjusted["estimate"]
                effects.loc[mask, "variance"] = adjusted["variance"]
                check = random_effects_reml_mkh(effects["estimate"], effects["variance"])
                if abs(check.estimate - float(contract["meta_estimate"])) > 1e-10:
                    raise RuntimeError("item-7 meta disagrees with the frozen contract")
                row = {
                    "estimate": float(contract["meta_estimate"]),
                    "ci_low_mkh": float(contract["meta_ci_low_mkh"]),
                    "ci_high_mkh": float(contract["meta_ci_high_mkh"]),
                    "p_mkh": float(contract["meta_p_mkh"]),
                    "tau2": float(contract["meta_tau2"]),
                    "i_squared": float(contract["meta_i_squared"]),
                    "k_missions": int(check.k),
                    "n_positive_missions": int((effects["estimate"] > 0).sum()),
                }
                notes = (
                    f"OSD-163 raw {float(contract['osd163_raw_estimate']):.6g} -> adjusted "
                    f"{float(contract['osd163_adjusted_estimate']):.6g}; outcome-QC r "
                    f"{float(contract['outcome_qc_pearson_r']):.3g}"
                )
            except _VARIANT_ERRORS as error:
                row, notes, effects = None, f"not evaluable: {error}", None
        else:
            notes = f"not evaluable: {failure}"
        if name == primary:
            gate = evaluate_gate(
                prereg.rules[7],
                estimate=np.nan if row is None else row["estimate"],
                reference=collect.reference[name],
            )
            gate_result = gate
        collect.meta_row(key, prereg.item7_primary, name, "mean", row, gate=gate, notes=notes)
        if effects is not None:
            collect.mission_effects(
                key, prereg.item7_primary, name, effects, reference_coverage[name]
            )

        if prereg.item7_secondary is not None and covariate is not None:
            try:
                score = reference_details[name][QC_MISSION][name].scores
                ancova = ancova_standardized_effect(score, metadata, covariate)
                effects2 = reference_effects[name].copy()
                mask = effects2["mission"] == QC_MISSION
                effects2.loc[mask, "estimate"] = ancova["estimate"]
                effects2.loc[mask, "variance"] = ancova["variance"]
                fit = random_effects_reml_mkh(effects2["estimate"], effects2["variance"])
            except _VARIANT_ERRORS as error:
                collect.meta_row(key, prereg.item7_secondary, name, "mean", None, notes=f"not evaluable: {error}")
            else:
                collect.meta_row(
                    key,
                    prereg.item7_secondary,
                    name,
                    "mean",
                    {
                        "estimate": fit.estimate,
                        "ci_low_mkh": fit.ci_low,
                        "ci_high_mkh": fit.ci_high,
                        "p_mkh": fit.p,
                        "tau2": fit.tau2,
                        "i_squared": fit.i_squared,
                        "k_missions": fit.k,
                        "n_positive_missions": int((effects2["estimate"] > 0).sum()),
                    },
                    notes=(
                        f"descriptive; OSD-163 ANCOVA g {ancova['estimate']:.6g} "
                        f"(df {ancova['df']:.0f}, J-corrected, unadjusted pooled within-group SD)"
                    ),
                )
                collect.mission_effects(
                    key, prereg.item7_secondary, name, effects2, reference_coverage[name]
                )
    statistic = None
    if gate_result is not None:
        ratios = [c["value"] for c in gate_result["clauses"] if c["kind"] == "ratio"]
        statistic = ratios[0] if ratios else None
    return {**gate_result, "statistic": statistic, "method": prereg.item7_primary}


def _median_reference(directory: Path | None) -> dict[str, Any]:
    if directory is None:
        return {"status": "not_requested", "table": None}
    path = Path(directory) / "compartment_context_meta_results.tsv"
    if not path.is_file():
        return {"status": f"absent ({path}); median computed here", "table": None}
    return {"status": f"read {path}", "table": pd.read_csv(path, sep="\t"), "path": str(path)}


def _median_file_note(median_file, name, estimate, info) -> str:
    table = median_file.get("table")
    if table is None or name not in set(table["axis"]):
        return "median computed here (family file absent or lacks this target)"
    family = table.set_index("axis")
    delta = abs(float(family.loc[name, "estimate"]) - estimate)
    info["checks"].setdefault("median_family_delta", {})[name] = delta
    if delta > CONSISTENCY_TOLERANCE:
        return (
            f"median family file disagrees by {delta:.3g} (different inputs); "
            "its descriptive columns were not used"
        )
    parts = ["median family file agrees"]
    if "max_t_fwer" in family.columns:
        parts.append(f"max_t_fwer={float(family.loc[name, 'max_t_fwer']):.6g}")
    if "t_mkh" in family.columns:
        rank = family["t_mkh"].abs().rank(ascending=False, method="min")
        parts.append(f"rank_abs_t={int(rank[name])}/{len(family)}")
    return "; ".join(parts) + " (descriptive)"


# --------------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------------

OUTPUT_FILES = {
    "summary": "podocyte_sensitivity_summary.tsv",
    "mission_effects": "podocyte_sensitivity_mission_effects.tsv",
    "leave_one_mission": "podocyte_leave_one_mission.tsv",
    "leave_one_gene": "podocyte_leave_one_gene.tsv",
    "gene_contributions": "podocyte_gene_contributions.tsv",
    "contributor_concentration": "podocyte_contributor_concentration.tsv",
    "contribution_strata_audit": "podocyte_contribution_strata_audit.tsv",
    "leave_top_k": "podocyte_leave_top_k.tsv",
    "k0": "podocyte_k0_sensitivities.tsv",
    "k0_mission_effects": "podocyte_k0_mission_effects.tsv",
    "k0_leave_one_gene": "podocyte_k0_leave_one_gene.tsv",
    "gene_coverage": "podocyte_sensitivity_gene_coverage.tsv",
}
GATES_FILE = "podocyte_sensitivity_gates.json"
MANIFEST_FILE = "podocyte_sensitivity_manifest.json"


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, Path):
        return str(value)
    if value is pd.NA:
        return None
    return value


def write_outputs(result: PipelineResult, outdir: Path, manifest: Mapping[str, Any]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    written = []
    for key, filename in OUTPUT_FILES.items():
        table = result.tables.get(key)
        if table is None:
            continue
        table.to_csv(outdir / filename, sep="\t", index=False)
        written.append(filename)
    (outdir / GATES_FILE).write_text(json.dumps(_jsonable(result.gates), indent=2) + "\n")
    full = {**manifest, "outputs": written + [GATES_FILE], "pipeline": result.info}
    (outdir / MANIFEST_FILE).write_text(json.dumps(_jsonable(full), indent=2) + "\n")


def headline_rows(result: PipelineResult, map_info: Mapping[str, Any]) -> list[dict[str, Any]]:
    """APPENDIX A item 3 headline outputs of one run."""

    gates = result.gates
    primary = gates["primary_target"]
    summary = result.tables["summary"]
    reference = summary[(summary["item"] == "0_reference") & (summary["target"] == primary)].iloc[0]
    k0 = result.tables["k0"]
    model = gates["k0"]["model"]
    p1 = k0[(k0["item"] == "0_reference") & (k0["model"] == model)].iloc[0]
    base = {"map": map_info["label"], "map_path": map_info["path"], "map_sha256": map_info["sha256"]}
    rows = [
        {
            **base,
            "output": f"{primary}__primary_estimate",
            "value": float(reference["estimate"]),
            "ci_low": float(reference["ci_low_mkh"]),
            "ci_high": float(reference["ci_high_mkh"]),
            "verdict": estimate_verdict(
                reference["estimate"], reference["ci_low_mkh"], reference["ci_high_mkh"]
            ),
        },
        {
            **base,
            "output": f"k0_{model}__coef_score_units",
            "value": float(p1["coef_score_units"]),
            "ci_low": float(p1["ci_low_mkh"]),
            "ci_high": float(p1["ci_high_mkh"]),
            "verdict": estimate_verdict(p1["coef_score_units"], p1["ci_low_mkh"], p1["ci_high_mkh"]),
        },
    ]
    for key, item in gates["items"].items():
        statistic = item.get("statistic")
        rows.append(
            {
                **base,
                "output": f"gate__{key}",
                "value": np.nan if statistic is None else statistic,
                "ci_low": np.nan,
                "ci_high": np.nan,
                "verdict": item["verdict"],
            }
        )
    rows.append(
        {
            **base,
            "output": "overall_classification",
            "value": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "verdict": gates["overall"]["classification"],
        }
    )
    return rows


def idmap_sensitivity_table(
    rows: Sequence[Mapping[str, Any]], *, robust_label: str, fragile_label: str
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Flag verdict flips across gene-ID maps; the headline takes the more conservative verdict."""

    table = pd.DataFrame(rows)
    order = {fragile_label: 0, "INCOMPLETE": 1, robust_label: 2}
    flips, conservative = [], {}
    for output, sub in table.groupby("output", sort=False):
        verdicts = sub["verdict"].astype(str).tolist()
        flipped = len(set(verdicts)) > 1
        if flipped:
            flips.append(output)
        conservative[output] = conservative_verdict(
            verdicts, order=order if output == "overall_classification" else None
        )
    table["verdict_flip"] = table["output"].isin(flips)
    table["conservative_verdict"] = table["output"].map(conservative)
    summary = {
        "flipped_outputs": flips,
        "robust_to_gene_id_reconstruction": not flips,
        "statement": (
            "not robust to gene-ID reconstruction: " + ", ".join(flips)
            if flips
            else "no verdict changes between the winner and runner-up gene-ID maps"
        ),
        "conservative_overall": conservative.get("overall_classification"),
        "rule": "any flip => not robust to gene-ID reconstruction; headline takes the more conservative verdict",
    }
    return table, summary


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------


def run(args: argparse.Namespace) -> dict[str, PipelineResult]:
    prereg, guard = load_preregistration(
        args.prereg, repo=REPO, allow_uncommitted=args.allow_uncommitted_prereg
    )
    if args.validate_prereg_only:
        _print_prereg(prereg, guard)
        return {}
    config = yaml.safe_load(Path(args.config).read_text())
    if not prereg.seed_base_from_parent:
        raise PreregistrationError("only seed_base_from_parent: true is supported")
    if args.gene_maps:
        maps = [("winner", Path(args.gene_maps[0])), ("runner_up", Path(args.gene_maps[1]))]
    else:
        maps = [("primary", None if args.gene_map is None else Path(args.gene_map))]
    root = Path(args.results) / OUTPUT_SUBDIR
    head = _git(REPO, "rev-parse", "HEAD")
    results: dict[str, PipelineResult] = {}
    rows: list[dict[str, Any]] = []
    for label, gene_map in maps:
        map_config = config_with_gene_map(config, gene_map)
        map_info = gene_map_record(map_config, label)
        inputs = build_inputs(map_config, Path(args.tiers), prereg)
        baseline_run = label in ("primary", "winner")
        options = PipelineOptions(
            lomo_permutations=int(args.lomo_permutations),
            chunk_size=int(args.chunk_size),
            k0_leave_one_gene=bool(args.k0_leave_one_gene),
            structural_gene_influence=not args.skip_structural_gene_influence,
            median_results=(
                Path(args.median_results)
                if args.median_results
                else Path(args.results) / "compartment_context_median"
            )
            if baseline_run
            else None,
            reference_run=(Path(args.reference_run) if args.reference_run else Path(args.results))
            if baseline_run
            else None,
        )
        result = run_pipeline(inputs, prereg, options)
        outdir = root if baseline_run else root / f"idmap_{label}"
        manifest = {
            "analysis": "G1 sensitivities of the frozen podocyte high-specificity set (items 1-7)",
            "status": prereg.raw.get("status", "post_hoc_sensitivities_preregistered_before_computation"),
            "preregistration": guard,
            "git_head": head.stdout.strip() if head.returncode == 0 else None,
            "config": str(args.config),
            "config_sha256": sha256(Path(args.config)),
            "gene_map": map_info,
            "tiers": str(args.tiers),
            "tiers_sha256": sha256(Path(args.tiers)),
            "cpm_threshold": inputs.threshold,
            "seed_base": inputs.seed_base,
            "options": {
                "lomo_permutations": options.lomo_permutations,
                "chunk_size": options.chunk_size,
                "k0_leave_one_gene": options.k0_leave_one_gene,
                "structural_gene_influence": options.structural_gene_influence,
                "top_ks": list(options.top_ks),
                "median_results": options.median_results,
                "reference_run": options.reference_run,
            },
            "conventions": {
                "k0": "coef_score_units are raw OLS score-unit coefficients (HC3 per mission, REML/mKH), not Hedges g",
                "retention_ratio": "signed estimate / reference; NA when the reference is zero",
                "leave_top_k_percent": "ceil(fraction x genes with a signed contribution)",
                "item7_ancova": "J-corrected, unadjusted pooled within-group SD, Hedges variance with df = n - rank(design)",
                "multiplicity": prereg.raw.get("multiplicity", ""),
            },
        }
        write_outputs(result, outdir, manifest)
        results[label] = result
        rows.extend(headline_rows(result, map_info))
        _print_verdicts(label, outdir, result)
    if len(maps) == 2:
        table, summary = idmap_sensitivity_table(
            rows, robust_label=prereg.robust_label, fragile_label=prereg.fragile_label
        )
        table.to_csv(root / "idmap_sensitivity.tsv", sep="\t", index=False)
        (root / "idmap_sensitivity_summary.json").write_text(
            json.dumps(_jsonable(summary), indent=2) + "\n"
        )
        print(summary["statement"])
    return results


def _print_prereg(prereg: PodocytePrereg, guard: Mapping[str, Any]) -> None:
    """No data are read: show how each preregistered rule was parsed."""

    print(f"preregistration {guard['path']} sha256 {guard['sha256']} committed={guard['committed']}")
    print(f"  required items {list(prereg.required_items)} -> {prereg.robust_label} / {prereg.fragile_label}")
    for number, rule in sorted(prereg.rules.items()):
        clauses = "; ".join(
            f"{c.kind}({c.op or ''}{'' if c.threshold is None else c.threshold}"
            f"{'' if c.count is None else ', n=' + str(c.count)})"
            for c in rule.clauses
        )
        print(f"  {prereg.item_keys[number]:28s} {rule.text}  ->  {clauses}")
    print(f"  {prereg.item_keys[4]:28s} gates {prereg.gene_influence_gates}; "
          f"secondary (reported only) {prereg.gene_influence_secondary}")
    print(f"  K0 rule model {prereg.k0_rule_model}; note '{prereg.k0_note_template}'")
    print(f"  LOMO seed offsets {prereg.seed_offsets.get('lomo_family')}")


def _print_verdicts(label: str, outdir: Path, result: PipelineResult) -> None:
    gates = result.gates
    print(f"=== podocyte sensitivities [{label}] -> {outdir}")
    for key, item in gates["items"].items():
        print(f"  {key:28s} {item['verdict']}")
    print(f"  overall: {gates['overall']['classification']}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS)
    parser.add_argument(
        "--results",
        type=Path,
        default=DEFAULT_RESULTS,
        help="run directory; outputs go to <results>/podocyte_sensitivities/",
    )
    parser.add_argument("--prereg", type=Path, default=DEFAULT_PREREG)
    parser.add_argument(
        "--allow-uncommitted-prereg",
        action="store_true",
        help="tests only: run although the preregistration YAML is not committed",
    )
    parser.add_argument(
        "--validate-prereg-only",
        action="store_true",
        help="parse the preregistration and print the rules; reads no data",
    )
    maps = parser.add_mutually_exclusive_group()
    maps.add_argument(
        "--gene-map",
        type=Path,
        default=None,
        help="replace config gene_mapping.path in memory (the config file is never edited)",
    )
    maps.add_argument(
        "--gene-maps",
        type=Path,
        nargs=2,
        metavar=("WINNER", "RUNNER_UP"),
        default=None,
        help="run under both maps and write idmap_sensitivity.tsv (APPENDIX A item 3)",
    )
    parser.add_argument(
        "--lomo-permutations",
        type=int,
        default=0,
        help="optional descriptive LOMO family-wise p over the fixed evaluable family",
    )
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument(
        "--k0-leave-one-gene",
        action="store_true",
        help="descriptive: K0 P1 after omitting each podocyte gene",
    )
    parser.add_argument(
        "--skip-structural-gene-influence",
        action="store_true",
        help="skip the descriptive item-4 diagnostics for the structural comparators",
    )
    parser.add_argument(
        "--median-results",
        type=Path,
        default=None,
        help="WP-S1 median family directory (default <results>/compartment_context_median)",
    )
    parser.add_argument(
        "--reference-run",
        type=Path,
        default=None,
        help=(
            "directory with compartment_context_meta_results.tsv that step 1 must reproduce "
            "(default <results>; skipped when absent)"
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        run(args)
    except PreregistrationError as error:
        print(f"refusing to run: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
