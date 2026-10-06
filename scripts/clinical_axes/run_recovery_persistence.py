#!/usr/bin/env python3
"""G2: does the terminal-flight podocyte / structural / compartment shift persist
through live-animal return (OSD-771 LAR, 32 d flight + 24 d on Earth) and early
post-landing readaptation (OSD-513)?

Preregistered in ``config/clinical_axes_recovery_persistence.yaml`` (plan §4.2). The
script refuses to touch expression data unless that YAML exists, is tracked by git
and has no uncommitted changes (``--allow-uncommitted-prereg`` is for tests only);
its sha256, the git HEAD and the YAML's last commit are written to the manifest.

Layers (all statistics in ``src/clinical_axes/persistence.py``):
  1. exact single-cohort tests: ISS-T and LAR (63,504 within-age allocations each) and
     OSD-513 (48,620), own eligibility and own z-scores, plus per-arm P1 (HC3, coef in
     score units) and D = podocyte - structural_disjoint;
  2. joint OSD-771 model (primary): intersection eligibility, z over the 40 animals,
     score ~ 0 + arm x age + F_T + F_L, HC3, interaction (kappa 1) and margin
     (kappa 0.5) with Satterthwaite df, Freedman--Lane within arm x age blocks
     (one permutation shared by every set), Fieller / bootstrap / flat-prior ratio
     summaries and the preregistered classification;
  3. structural-adjusted persistence for the podocyte set (arm-specific slopes);
  4. pattern concordance with exact label nulls; plus descriptive duration context.

Modes: full run (default); ``--power-only`` (simulation only, loads NO expression
data); ``--idmap-sensitivity MAP [MAP ...]`` (headline T1-T3 under each gene-ID map,
APPENDIX A item 3). ``--gene-map`` replaces ``gene_mapping.path`` in memory only.

Interpretation boundary: recovery is aliased with shorter flight (32 vs 53-56 d),
on-ISS vs on-Earth euthanasia and handling, and re-entry/landing; the podocyte set was
selected partly on OSD-771, so theta_ISS-T is winner's-curse inflated and rho is
biased low. OSD-513 is ~1 d post-landing with a 3 d ground-control offset.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
from typing import Callable, Mapping, Sequence

import numpy as np
import pandas as pd
import scipy
from scipy import stats
import yaml

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes import persistence as P  # noqa: E402
from src.clinical_axes.analysis import (  # noqa: E402
    mission_effect_from_score,
    sample_manifest,
    technical_qc_audit,
)
from src.clinical_axes.data import MissionData, cpm_eligible_genes  # noqa: E402
from src.clinical_axes.statistics import random_effects_reml_mkh  # noqa: E402

DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_PREREG = REPO / "config/clinical_axes_recovery_persistence.yaml"
DEFAULT_TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
DEFAULT_RESULTS = REPO / "data/results/clinical-axes/recovery_persistence_manual"
OUTPUT_SUBDIR = "recovery_persistence"

DEFAULT_SEED_OFFSETS = {
    "freedman_lane": 700,
    "bootstrap": 701,
    "power_simulation": 702,
    "posterior": 703,
}
DEFAULT_N_PERM = 100_000
DEFAULT_N_BOOT = 10_000
DEFAULT_POSTERIOR_DRAWS = 1_000_000
POWER_THETAS = (0.689, 1.174)
POWER_RHOS = (1.0, 0.5, 0.0, -0.5)

T1, T2, T3 = P.PODOCYTE_SET, P.STRUCTURAL_DISJOINT, P.D_CONTRAST
L1_COHORTS = ("OSD-771", "OSD-771-LAR", "OSD-513")


# --------------------------------------------------------------------------- provenance


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class PreregistrationError(SystemExit):
    """The preregistration YAML is missing, untracked or has uncommitted changes."""


def _git(runner: Callable, args: Sequence[str], cwd: Path) -> subprocess.CompletedProcess | None:
    try:
        return runner(["git", *args], cwd=str(cwd), capture_output=True, text=True, check=False)
    except OSError:
        return None


def prereg_guard(
    path: Path,
    allow_uncommitted: bool = False,
    *,
    runner: Callable = subprocess.run,
) -> dict[str, object]:
    """Refuse unless the YAML exists, is tracked, and ``git status --porcelain <path>`` is empty.

    git runs in the YAML's own directory, so the check applies to whichever repository
    holds the file. Returns the manifest record (sha256, git HEAD, last commit, outcome).
    """
    path = Path(path).resolve()
    problem = None
    cwd = path.parent if path.parent.exists() else REPO
    if not path.exists():
        problem = "preregistration YAML is missing"
    else:
        status = _git(runner, ["status", "--porcelain", "--", str(path)], cwd)
        if status is None or status.returncode != 0:
            problem = "git status failed: " + ("" if status is None else status.stderr.strip())
        elif status.stdout.strip():
            problem = "uncommitted changes: " + status.stdout.strip()
        else:
            tracked = _git(runner, ["ls-files", "--error-unmatch", "--", str(path)], cwd)
            if tracked is None or tracked.returncode != 0:
                problem = "preregistration YAML is not tracked by git"
    if problem is not None and not allow_uncommitted:
        raise PreregistrationError(
            f"refusing to run G2 recovery/persistence: {problem} ({path}). Commit the "
            "preregistration first; --allow-uncommitted-prereg is for tests only."
        )
    head = _git(runner, ["rev-parse", "HEAD"], cwd)
    last = _git(runner, ["log", "-1", "--format=%H", "--", str(path)], cwd) if path.exists() else None
    return {
        "prereg_path": str(path),
        "prereg_sha256": sha256(path) if path.exists() else None,
        "git_head": head.stdout.strip() if head is not None and head.returncode == 0 else None,
        "prereg_last_commit": (
            last.stdout.strip() if last is not None and last.returncode == 0 and last.stdout.strip() else None
        ),
        "prereg_guard": "clean" if problem is None else f"bypassed: {problem}",
    }


def load_prereg(path: Path) -> dict[str, object]:
    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError(f"{path}: preregistration must be a YAML mapping")
    return data


def apply_gene_map_override(config: Mapping[str, object], gene_map: Path | str | None) -> dict:
    """Return a deep copy whose ``gene_mapping.path`` points at ``gene_map`` (config files are never edited)."""
    updated = copy.deepcopy(dict(config))
    if gene_map is not None:
        resolved = Path(gene_map).expanduser()
        resolved = resolved if resolved.is_absolute() else (REPO / resolved)
        resolved = resolved.resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"--gene-map {resolved} does not exist")
        updated["gene_mapping"] = dict(updated["gene_mapping"])
        updated["gene_mapping"]["path"] = str(resolved)
    return updated


# --------------------------------------------------------------------------- settings


@dataclass(frozen=True)
class Settings:
    """Run settings, read from the preregistration (CLI overrides are recorded as deviations)."""

    seed: int
    offsets: Mapping[str, int]
    n_perm: int
    n_boot: int
    posterior_draws: int
    min_genes: int
    cpm_threshold: float
    chunk: int = 2048
    primary_target: str = T1
    key_secondary: tuple[str, ...] = (T2, T3)
    deviations: tuple[str, ...] = field(default_factory=tuple)
    margin_kappas: tuple[float, ...] = (0.5,)
    primary_margin_kappa: float = 0.5
    sweep_sets: tuple[str, ...] = (T1, T2, T3)
    sweep_all_sets: bool = False

    def seed_for(self, name: str) -> int:
        return int(self.seed) + int(self.offsets[name])


def _first_int(pattern: str, text: object) -> int | None:
    if text is None:
        return None
    match = re.search(pattern, str(text), flags=re.IGNORECASE)
    if match is None:
        return None
    return int(match.group(1).replace(",", "").replace("_", ""))


def prereg_settings(
    prereg: Mapping[str, object] | None,
    config: Mapping[str, object],
    *,
    n_perm: int | None = None,
    n_boot: int | None = None,
    posterior_draws: int | None = None,
    chunk: int = 2048,
    sweep_all_sets: bool = False,
) -> Settings:
    prereg = prereg or {}
    offsets = dict(DEFAULT_SEED_OFFSETS)
    offsets.update({str(k): int(v) for k, v in (prereg.get("seed_offsets") or {}).items()})
    inference = prereg.get("inference") or {}
    pre_perm = _first_int(r"(\d{1,3}(?:,\d{3})+|\d{4,})", inference.get("permutation")) or DEFAULT_N_PERM
    pre_boot = _first_int(r"bootstrap\s+(\d{1,3}(?:,\d{3})+|\d+)", inference.get("ratio")) or DEFAULT_N_BOOT
    min_genes = _first_int(r">=\s*(\d+)", prereg.get("eligibility")) or 8
    deviations = []
    if n_perm is not None and int(n_perm) != pre_perm:
        deviations.append(f"Freedman-Lane permutations {n_perm} instead of preregistered {pre_perm}")
    if n_boot is not None and int(n_boot) != pre_boot:
        deviations.append(f"bootstrap resamples {n_boot} instead of preregistered {pre_boot}")
    if posterior_draws is not None and int(posterior_draws) != DEFAULT_POSTERIOR_DRAWS:
        deviations.append(f"posterior draws {posterior_draws} instead of {DEFAULT_POSTERIOR_DRAWS}")
    targets = prereg.get("targets") or {}
    primary = str(targets.get("primary", T1))
    key_secondary = tuple(str(t) for t in targets.get("key_secondary", (T2, T3)))
    if primary != T1:
        raise ValueError(f"preregistered primary target {primary!r} is not {T1!r}")
    if set(key_secondary) != {T2, T3}:
        raise ValueError(f"preregistered key secondary targets {key_secondary} are not {(T2, T3)}")
    contrasts = prereg.get("contrasts") or P.PLAN_4_2_CONTRASTS
    primary_kappa = P._kappa(contrasts["margin"], "margin")
    kappas = P.parse_margin_sweep(prereg, primary_kappa)
    return Settings(
        seed=int(config["seed"]),
        offsets=offsets,
        n_perm=int(n_perm if n_perm is not None else pre_perm),
        n_boot=int(n_boot if n_boot is not None else pre_boot),
        posterior_draws=int(posterior_draws if posterior_draws is not None else DEFAULT_POSTERIOR_DRAWS),
        min_genes=int(min_genes),
        cpm_threshold=float(config["eligibility"]["cpm_threshold"]),
        chunk=int(chunk),
        primary_target=primary,
        key_secondary=key_secondary,
        deviations=tuple(deviations),
        margin_kappas=kappas,
        primary_margin_kappa=float(primary_kappa),
        sweep_sets=P.margin_sweep_targets(prereg),
        sweep_all_sets=bool(sweep_all_sets),
    )


def _resolve(path: Path | str) -> Path:
    value = Path(path).expanduser()
    return (value if value.is_absolute() else REPO / value).resolve()


def gene_map_record(map_path: Path | str, prereg: Mapping[str, object] | None) -> dict[str, object]:
    """Role of the analysis gene-ID map relative to the preregistered ``gene_map`` block and hash check."""
    resolved = _resolve(map_path)
    actual = sha256(resolved) if resolved.exists() else None
    block = (prereg or {}).get("gene_map") or {}
    record: dict[str, object] = {"path": str(resolved), "sha256": actual, "role": "not_declared_in_prereg",
                                 "expected_sha256": None, "sha256_matches": None}
    for role in ("baseline", "runner_up"):
        declared = block.get(role)
        if declared and _resolve(declared) == resolved:
            expected = block.get(f"{role}_sha256")
            record.update(role=role, expected_sha256=expected,
                          sha256_matches=None if expected is None else actual == expected)
            return record
    if block:
        record["role"] = "not_a_preregistered_map"
    return record


def check_gene_map(record: Mapping[str, object], allow: bool) -> list[str]:
    """Refuse a declared map whose bytes differ from the preregistered hash; flag undeclared maps."""
    deviations = []
    if record["sha256_matches"] is False:
        message = (f"gene map {record['path']} ({record['role']}) sha256 {record['sha256']} differs from the "
                   f"preregistered {record['expected_sha256']}")
        if not allow:
            raise PreregistrationError("refusing to run G2 recovery/persistence: " + message)
        deviations.append(message)
    if record["role"] == "not_a_preregistered_map":
        deviations.append(f"analysis gene map {record['path']} is neither the preregistered baseline nor runner-up")
    return deviations


# --------------------------------------------------------------------------- set definitions


@dataclass
class SetDefinitions:
    """The 49-set compartment family, structural_disjoint and the 4 frozen axes."""

    compartment: dict[str, Mapping[str, object]]
    structural_disjoint: Mapping[str, object]
    axes: dict[str, Mapping[str, object]]
    definition_audit: pd.DataFrame
    primary_counts: pd.DataFrame
    structural_overlap: list[str]

    def layer_family(self, sets: Sequence[str] | None = None) -> dict[str, Mapping[str, object]]:
        family: dict[str, Mapping[str, object]] = {}
        if T1 in self.compartment:
            family[T1] = self.compartment[T1]
        family[T2] = self.structural_disjoint
        for name, spec in self.compartment.items():
            family.setdefault(name, spec)
        for name, spec in self.axes.items():
            family.setdefault(name, spec)
        if sets is not None:
            wanted = set(sets) | ({T1, T2} if T3 in sets else set())
            family = {name: spec for name, spec in family.items() if name in wanted}
        return family

    def scope(self, name: str) -> str:
        if name == T1:
            return "primary"
        if name in (T2, T3):
            return "key_secondary"
        if name in self.compartment:
            return "exploratory"
        if name in self.axes:
            return "context"
        return "other"


def build_set_definitions(
    tiers_path: Path,
    primary: Mapping[str, MissionData],
    axes_family: Mapping[str, Mapping[str, object]],
    threshold: float,
    min_genes: int = 8,
) -> SetDefinitions:
    # Lazy import: run_compartment_context is a stable contract but may be mid-edit.
    from scripts.clinical_axes.run_compartment_context import load_family

    family, audit = load_family(Path(tiers_path))
    compartment, counts = P.primary_evaluable_family(family, primary, threshold, min_genes)
    if T1 not in compartment:
        raise RuntimeError(f"{T1} is not evaluable in the primary missions")
    disjoint, overlap = P.structural_disjoint_spec(family)
    return SetDefinitions(
        compartment=compartment,
        structural_disjoint=disjoint,
        axes=dict(axes_family),
        definition_audit=audit,
        primary_counts=counts,
        structural_overlap=overlap,
    )


def _ordered_columns(columns: Sequence[str], defs: SetDefinitions) -> list[str]:
    first = [c for c in (T1, T2, T3) if c in columns]
    rest = [c for c in columns if c not in first and c in defs.compartment]
    axes = [c for c in columns if c not in first and c not in rest]
    return first + rest + axes


# --------------------------------------------------------------------------- terminal reference


def pooled_terminal(
    primary: Mapping[str, MissionData],
    family: Mapping[str, Mapping[str, object]],
    threshold: float,
    min_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Five-mission REML/mKH pooled g per set (observed only; already-known terminal data)."""
    rows = []
    for mission, data in primary.items():
        scores, _ = P.score_arm(data, family, threshold, min_genes=min_genes)
        observed = P.observed_arm_effects(scores, data.metadata)
        for name, row in observed.iterrows():
            rows.append({"mission": mission, "set": name, "estimate": row["estimate_g"], "variance": row["variance"]})
    effects = pd.DataFrame(rows)
    pooled = []
    for name, sub in effects.groupby("set", sort=False):
        if len(sub) < 2:
            continue
        fit = random_effects_reml_mkh(sub["estimate"].to_numpy(), sub["variance"].to_numpy())
        pooled.append(
            {
                "set": name,
                "pooled_terminal_estimate": fit.estimate,
                "pooled_terminal_ci_low_mkh": fit.ci_low,
                "pooled_terminal_ci_high_mkh": fit.ci_high,
                "pooled_terminal_k": fit.k,
                "pooled_terminal_n_positive": int((sub["estimate"] > 0).sum()),
            }
        )
    return pd.DataFrame(pooled).set_index("set"), effects


def orientation_for(columns: Sequence[str], terminal: pd.DataFrame) -> pd.Series:
    signs = {}
    for name in columns:
        value = terminal["pooled_terminal_estimate"].get(name, np.nan)
        signs[name] = -1.0 if np.isfinite(value) and value < 0 else 1.0
    return pd.Series(signs, name="orientation")


# --------------------------------------------------------------------------- layer 1


def arm_p1(
    scores: pd.DataFrame, meta: pd.DataFrame, n_perm: int, seed: object, chunk: int
) -> dict[str, object]:
    """Per-arm P1: podocyte ~ 0 + block + flight + structural_disjoint (HC3, coef in score units)."""
    blocks = pd.get_dummies(meta["block"], dtype=float).to_numpy()
    flight = meta["condition"].eq("FLT").to_numpy(dtype=float)
    structural = scores.loc[meta.index, T2].to_numpy(dtype=float)
    y = scores.loc[meta.index, T1].to_numpy(dtype=float)
    x = np.column_stack([blocks, flight, structural])
    fit = P.ols_hc3(x, y)
    index = blocks.shape[1]
    coef = float(fit.beta[index, 0])
    se = float(np.sqrt(fit.cov_hc3[0, index, index]))
    df = fit.df_resid
    q = float(P._t_quantile(0.975, df))
    fl = P.freedman_lane(
        y, np.column_stack([blocks, structural]), flight, meta["block"].to_numpy(), n_perm, seed, chunk, keep_null=False
    )
    return {
        "coef_score_units": coef,
        "coef_se_hc3": se,
        "coef_df": df,
        "coef_ci_low_hc3": coef - q * se,
        "coef_ci_high_hc3": coef + q * se,
        "p_hc3": float(2.0 * stats.t.sf(abs(coef / se), df)),
        "p_fl": float(fl.p_two_sided[0]),
        "fl_t_identity_matches_hc3": bool(abs(fl.t_obs[0] - coef / se) < 1e-9),
    }


def layer1(
    cohorts: Mapping[str, object],
    defs: SetDefinitions,
    settings: Settings,
    flagged: Sequence[str],
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame], pd.DataFrame, dict[str, P.AllocationSet]]:
    family = defs.layer_family()
    data_by_cohort = {"OSD-771": cohorts["iss"], "OSD-771-LAR": cohorts["lar"], "OSD-513": cohorts["osd513"]}
    tables, scores_by_cohort, coverage, allocations = [], {}, [], {}
    for index, (cohort, data) in enumerate(data_by_cohort.items(), start=1):
        scores, cov = P.score_arm(data, family, settings.cpm_threshold, min_genes=settings.min_genes)
        scores = scores[_ordered_columns(list(scores.columns), defs)]
        scores_by_cohort[cohort] = scores
        coverage.append(cov)
        meta = data.metadata
        n_flight = {str(b): int(s["condition"].eq("FLT").sum()) for b, s in meta.groupby("block", sort=False)}
        alloc = P.enumerate_block_allocations(meta["block"], n_flight, observed=meta["condition"].eq("FLT").to_numpy())
        allocations[cohort] = alloc
        families = {
            "compartment_49": [s for s in scores.columns if s in defs.compartment],
            "key_secondary": [s for s in (T1, T2, T3) if s in scores.columns],
            "clinical_axes_4": [s for s in scores.columns if s in defs.axes],
        }
        table = P.arm_exact_test(scores, meta, families, allocations=alloc)
        table.insert(0, "context", data.context)
        table.insert(0, "cohort", cohort)
        table.insert(2, "scope", [defs.scope(s) for s in table["set"]])
        table.insert(3, "statistic", "hedges_g_fixed_effect_over_strata")
        table["technical_flag_primary_metric"] = cohort in flagged
        if T1 in scores and T2 in scores:
            p1 = arm_p1(scores, meta, settings.n_perm, np.random.SeedSequence([settings.seed_for("freedman_lane"), index]), settings.chunk)
            p1_row = {
                "cohort": cohort,
                "context": data.context,
                "set": P.P1_ADJUSTED,
                "scope": "secondary",
                "statistic": "hc3_coef_score_units",
                "n_flight": int(meta["condition"].eq("FLT").sum()),
                "n_control": int(meta["condition"].eq("GC").sum()),
                "n_strata": meta["block"].nunique(),
                "technical_flag_primary_metric": cohort in flagged,
                **p1,
            }
            table = pd.concat([table, pd.DataFrame([p1_row])], ignore_index=True, sort=False)
        tables.append(table)
    return pd.concat(tables, ignore_index=True, sort=False), scores_by_cohort, pd.concat(coverage, ignore_index=True), allocations


# --------------------------------------------------------------------------- layer 2


def _classifier_inputs(table: pd.DataFrame, technical_flag: bool, p_persist: str, p_decay: str) -> dict[str, np.ndarray]:
    return {
        "technical_flag": np.full(len(table), float(technical_flag)),
        "theta_isst": table["theta_isst"].to_numpy(float),
        "se_isst": table["se_isst"].to_numpy(float),
        "df_isst": table["df_isst"].to_numpy(float),
        "theta_lar": table["theta_lar"].to_numpy(float),
        "se_lar": table["se_lar"].to_numpy(float),
        "df_lar": table["df_lar"].to_numpy(float),
        "margin": table["margin"].to_numpy(float),
        "margin_se": table["margin_se"].to_numpy(float),
        "margin_df": table["margin_df"].to_numpy(float),
        "p_persist": table[p_persist].to_numpy(float),
        "p_decay": table[p_decay].to_numpy(float),
        "pooled_terminal_estimate": table["pooled_terminal_estimate_oriented"].to_numpy(float),
    }


def _classify_table(table: pd.DataFrame, rules: P.PersistenceRules, technical_flag: bool) -> pd.DataFrame:
    inputs = _classifier_inputs(table, technical_flag, "p_persist_rule", "p_decay_rule")
    decision = P.classify_persistence_arrays(inputs, rules)
    labels, reasons = [], []
    for i in range(len(table)):
        row = {key: value[i] for key, value in inputs.items()}
        label, reason = P.classify_persistence(row, rules)
        labels.append(label)
        reasons.append(reason)
    table["isst_effect_established"] = ~(decision["includes_zero"] | decision["sign_mismatch"])
    table["classification"] = labels
    table["classification_base"] = [label.split(":")[0] for label in labels]
    table["classification_sublabel"] = [label.split(":")[1] if ":" in label else "" for label in labels]
    table["classification_reason"] = reasons
    return table


def layer2(
    iss: MissionData,
    lar: MissionData,
    defs: SetDefinitions,
    terminal: pd.DataFrame,
    settings: Settings,
    rules: P.PersistenceRules,
    technical_flag: bool,
    *,
    sets: Sequence[str] | None = None,
) -> dict[str, object]:
    expr40, _, meta40 = P.joint_frame(iss, lar)
    eligible = P.intersection_eligible(iss, lar, settings.cpm_threshold)
    family = defs.layer_family(sets)
    joint, coverage = P.joint_scores(expr40, family, eligible, min_genes=settings.min_genes)
    if T1 in joint and T2 in joint:
        joint[T3] = joint[T1] - joint[T2]
    if sets is not None:
        joint = joint[[c for c in joint.columns if c in set(sets)]]
    joint = joint[_ordered_columns(list(joint.columns), defs)]
    orientation = orientation_for(joint.columns, terminal)
    oriented = joint * orientation
    families = {
        "key_secondary": [s for s in (T1, T2, T3) if s in joint.columns],
        "exploratory": [s for s in joint.columns if s in defs.compartment],
    }
    contrasts = {"interaction": rules.interaction_kappa, "margin": rules.margin_kappa}
    table, fl = P.joint_hc3_contrasts(
        oriented,
        meta40,
        contrasts=contrasts,
        n_perm=settings.n_perm,
        seed=settings.seed_for("freedman_lane"),
        chunk=settings.chunk,
        families=families,
    )
    table = table.reset_index()
    table.insert(1, "family", [defs.scope(s) for s in table["set"]])
    table.insert(2, "orientation", table["set"].map(orientation).to_numpy())
    pooled = terminal.reindex(table["set"])
    table["pooled_terminal_estimate"] = pooled["pooled_terminal_estimate"].to_numpy()
    table["pooled_terminal_estimate_oriented"] = table["pooled_terminal_estimate"] * table["orientation"]
    n_used = coverage.set_index("set")["n_used_joint"]
    table["n_genes_joint"] = [n_used.get(s, np.nan) if s != T3 else n_used.get(T1, 0) + n_used.get(T2, 0) for s in table["set"]]
    table = table.rename(columns={"p_margin_fl_upper": "p_persist_fl", "p_margin_fl_lower": "p_decay_fl"})
    key = table["set"].isin([T1, T2, T3])
    expl = table["set"].isin(families["exploratory"])
    table["interaction_maxT_fwer_fl"] = np.where(
        key, table["interaction_maxT_fwer_fl__key_secondary"], table["interaction_maxT_fwer_fl__exploratory"]
    )
    table["interaction_bh_q"] = np.nan
    table.loc[expl, "interaction_bh_q"] = P.benjamini_hochberg(table.loc[expl, "p_interaction_fl"].to_numpy())
    holm_rows = table["set"].isin([T2, T3])
    for column in ("p_persist_fl", "p_decay_fl"):
        table[f"{column}_holm_key_secondary"] = np.nan
        table.loc[holm_rows, f"{column}_holm_key_secondary"] = P.holm(table.loc[holm_rows, column].to_numpy())
    table["p_persist_rule"] = np.where(holm_rows, table["p_persist_fl_holm_key_secondary"], table["p_persist_fl"])
    table["p_decay_rule"] = np.where(holm_rows, table["p_decay_fl_holm_key_secondary"], table["p_decay_fl"])
    with np.errstate(divide="ignore", invalid="ignore"):
        table["rho_hat"] = table["theta_lar"] / table["theta_isst"]
    ratio_df = np.minimum(table["df_isst"].to_numpy(float), table["df_lar"].to_numpy(float))
    for level in (90, 95):
        bounds = [
            P.fieller_ratio(tT, vT, tL, vL, df, level / 100.0)
            for tT, vT, tL, vL, df in zip(
                table["theta_isst"], table["se_isst"] ** 2, table["theta_lar"], table["se_lar"] ** 2, ratio_df
            )
        ]
        table[f"rho_fieller{level}_low"] = [b[0] for b in bounds]
        table[f"rho_fieller{level}_high"] = [b[1] for b in bounds]
        table[f"rho_fieller{level}_bounded"] = [b[2] for b in bounds]
    # Fieller 90% lower bound = largest margin compatible with persistence at one-sided 5%
    # (inversion of the HC3 margin test with the per-arm df); NaN when the set is unbounded.
    table["rho_lower90_fieller"] = table["rho_fieller90_low"]
    boot = P.stratified_bootstrap_ratio(oriented, meta40, B=settings.n_boot, seed=settings.seed_for("bootstrap"))
    for column in ("rho_boot95_low", "rho_boot95_high", "boot_frac_thetaT_le0", "bootstrap_unreliable"):
        table[column] = table["set"].map(boot[column]).to_numpy()
    post = P.posterior_summaries(
        table["theta_isst"].to_numpy(float),
        table["se_isst"].to_numpy(float) ** 2,
        table["theta_lar"].to_numpy(float),
        table["se_lar"].to_numpy(float) ** 2,
        draws=settings.posterior_draws,
        seed=settings.seed_for("posterior"),
        kappa=rules.margin_kappa,
    )
    for column in ("post_p_rho_ge_0_5", "post_p_thetaL_le_0", "post_p_rho_ge_1"):
        table[column] = post[column].to_numpy()
    kappa = rules.margin_kappa
    table["margin_g"] = table["g_lar_joint"] - kappa * table["g_isst_joint"]
    table["margin_g_var"] = table["g_lar_joint_var"] + kappa**2 * table["g_isst_joint_var"]
    z90 = stats.norm.ppf(0.95)
    table["margin_g_ci90_low"] = table["margin_g"] - z90 * np.sqrt(table["margin_g_var"])
    table["margin_g_ci90_high"] = table["margin_g"] + z90 * np.sqrt(table["margin_g_var"])
    table["technical_flag"] = bool(technical_flag)
    table = _classify_table(table, rules, technical_flag)
    table["classification_g_scale"] = [
        P.classify_g_scale(m, v, label, rules)
        for m, v, label in zip(table["margin_g"], table["margin_g_var"], table["classification"])
    ]
    table["g_scale_label_agrees"] = table["classification_g_scale"] == table["classification_base"]
    table["classification_scope"] = table["family"].map(
        {
            "primary": "confirmatory_primary",
            "key_secondary": "key_secondary_holm",
            "exploratory": "exploratory_descriptive",
            "context": "context_descriptive",
        }
    )
    return {
        "table": table,
        "fl": fl,
        "joint": joint,
        "oriented": oriented,
        "orientation": orientation,
        "coverage": coverage,
        "meta40": meta40,
        "eligible": eligible,
    }


CONTRAST_COLUMNS = [
    "set", "family", "classification_scope", "orientation", "n_genes_joint",
    "pooled_terminal_estimate", "theta_isst", "se_isst", "df_isst", "theta_isst_ci95_low", "theta_isst_ci95_high",
    "theta_lar", "se_lar", "df_lar", "theta_lar_ci95_low", "theta_lar_ci95_high", "cov_isst_lar_hc3",
    "g_isst_joint", "g_lar_joint",
    "interaction", "interaction_se", "interaction_df", "interaction_ci95_low", "interaction_ci95_high",
    "p_interaction_hc3", "p_interaction_fl", "interaction_maxT_fwer_fl",
    "interaction_maxT_fwer_fl__key_secondary", "interaction_maxT_fwer_fl__exploratory", "interaction_bh_q",
    "margin", "margin_se", "margin_df", "margin_ci90_low", "margin_ci90_high",
    "p_margin_hc3_upper", "p_margin_hc3_lower", "p_persist_fl", "p_decay_fl",
    "p_persist_fl_holm_key_secondary", "p_decay_fl_holm_key_secondary", "p_persist_rule", "p_decay_rule",
    "rho_hat", "rho_lower90_fieller", "rho_fieller90_low", "rho_fieller90_high", "rho_fieller90_bounded",
    "rho_fieller95_low", "rho_fieller95_high", "rho_fieller95_bounded",
    "rho_boot95_low", "rho_boot95_high", "boot_frac_thetaT_le0", "bootstrap_unreliable",
    "post_p_rho_ge_0_5", "post_p_thetaL_le_0", "post_p_rho_ge_1",
    "margin_g", "margin_g_ci90_low", "margin_g_ci90_high", "classification_g_scale", "g_scale_label_agrees",
    "technical_flag", "isst_effect_established", "classification", "classification_base",
    "classification_sublabel", "classification_reason",
]


# --------------------------------------------------------------------------- margin sweep


SWEEP_COLUMNS = [
    "set", "family", "kappa", "is_primary", "margin", "se", "df", "ci90_low", "ci90_high",
    "p_persist_fl", "p_decay_fl", "label_descriptive",
]


def margin_sweep_table(
    l2: Mapping[str, object],
    defs: SetDefinitions,
    rules: P.PersistenceRules,
    settings: Settings,
    technical_flag: bool,
) -> pd.DataFrame:
    """Preregistered margin sweep (secondary, descriptive): M_k = theta_LAR - k theta_ISS-T.

    Only k = the primary margin kappa sets the classification; the per-k labels apply the
    same rule (technical / terminal-effect gates, FL one-sided p AND HC3 90% CI, Holm over
    {structural_disjoint, D} per direction within each k) and are descriptive. k = 0 is
    "any retention" (theta_LAR > 0), k = 1 "no decay" (the interaction test).
    """
    oriented = l2["oriented"]
    sets = list(oriented.columns) if settings.sweep_all_sets else [s for s in settings.sweep_sets if s in oriented]
    sweep = P.margin_sweep(
        oriented[sets],
        l2["meta40"],
        settings.margin_kappas,
        primary_kappa=settings.primary_margin_kappa,
        n_perm=settings.n_perm,
        seed=settings.seed_for("freedman_lane"),
        chunk=settings.chunk,
    )
    sweep.insert(1, "family", [defs.scope(s) for s in sweep["set"]])
    terminal = l2["table"].set_index("set")["pooled_terminal_estimate_oriented"]
    sweep["pooled_terminal_estimate_oriented"] = sweep["set"].map(terminal).to_numpy()
    sweep["p_persist_rule"] = sweep["p_persist_fl"]
    sweep["p_decay_rule"] = sweep["p_decay_fl"]
    for kappa in settings.margin_kappas:
        rows = (sweep["kappa"] == kappa) & sweep["set"].isin([T2, T3])
        for column in ("persist", "decay"):
            sweep.loc[rows, f"p_{column}_rule"] = P.holm(sweep.loc[rows, f"p_{column}_fl"].to_numpy())
    labels = sweep.rename(columns={"se": "margin_se", "df": "margin_df"})
    labels["technical_flag"] = bool(technical_flag)
    labels = _classify_table(labels, rules, technical_flag)
    sweep["label_descriptive"] = labels["classification"].to_numpy()
    sweep["label_reason"] = labels["classification_reason"].to_numpy()
    # Prereg multiplicity.margin_sweep = none: no correction across kappas. The per-target
    # rule (Holm over the two key-secondary targets) is kept so the primary kappa
    # reproduces the classification; the fully unadjusted label is shown alongside.
    raw = labels.copy()
    raw["p_persist_rule"], raw["p_decay_rule"] = raw["p_persist_fl"], raw["p_decay_fl"]
    sweep["label_unadjusted"] = _classify_table(raw, rules, technical_flag)["classification"].to_numpy()
    sweep["scope"] = np.where(sweep["is_primary"], "primary_margin_sets_classification", "secondary_descriptive")
    return sweep[SWEEP_COLUMNS + [c for c in sweep.columns if c not in SWEEP_COLUMNS]]


# --------------------------------------------------------------------------- layer 3


def layer3(
    l2: Mapping[str, object],
    settings: Settings,
    rules: P.PersistenceRules,
    technical_flag: bool,
    terminal: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, P.FreedmanLaneResult]]:
    """P1-type adjusted persistence for the podocyte set (arm-specific structural slopes)."""
    oriented = l2["oriented"]
    meta40 = l2["meta40"]
    if T1 not in oriented or T2 not in oriented:
        return pd.DataFrame(), {}
    iss = meta40["arm"].eq("ISS-T").to_numpy(dtype=float)
    structural = l2["joint"][T2].loc[meta40.index].to_numpy(dtype=float)
    extra = np.column_stack([structural * iss, structural * (1.0 - iss)])
    table, fl = P.joint_hc3_contrasts(
        oriented[[T1]].rename(columns={T1: P.P1_ADJUSTED}),
        meta40,
        contrasts={"interaction": rules.interaction_kappa, "margin": rules.margin_kappa},
        n_perm=settings.n_perm,
        seed=settings.seed_for("freedman_lane"),
        chunk=settings.chunk,
        extra_nuisance=extra,
    )
    table = table.reset_index().rename(columns={"p_margin_fl_upper": "p_persist_fl", "p_margin_fl_lower": "p_decay_fl"})
    table.insert(1, "adjuster", f"{T2} (arm-specific slopes)")
    table.insert(2, "orientation", float(l2["orientation"][T1]))
    table["pooled_terminal_estimate_oriented"] = abs(float(terminal["pooled_terminal_estimate"].get(T1, np.nan)))
    table["p_persist_rule"] = table["p_persist_fl"]
    table["p_decay_rule"] = table["p_decay_fl"]
    with np.errstate(divide="ignore", invalid="ignore"):
        table["rho_hat"] = table["theta_lar"] / table["theta_isst"]
    table["technical_flag"] = bool(technical_flag)
    table = _classify_table(table, rules, technical_flag)
    table["classification_scope"] = "secondary"
    drop = [c for c in table.columns if c.startswith("g_") or c.startswith("resid_sd")]
    return table.drop(columns=drop), fl


# --------------------------------------------------------------------------- layer 4 / context


def layer4(
    arm_effects: pd.DataFrame,
    arm_scores: Mapping[str, pd.DataFrame],
    cohorts: Mapping[str, object],
    defs: SetDefinitions,
    terminal: pd.DataFrame,
    allocations: Mapping[str, P.AllocationSet],
) -> pd.DataFrame:
    rows = []
    iss_g = arm_effects[(arm_effects["cohort"] == "OSD-771") & arm_effects["set"].isin(defs.compartment)].set_index("set")["estimate_g"]
    lar_scores = arm_scores["OSD-771-LAR"]
    sets_a = [s for s in lar_scores.columns if s in defs.compartment and s in iss_g.index]
    result = P.pattern_concordance_exact(
        iss_g.loc[sets_a], lar_scores[sets_a], cohorts["lar"].metadata, allocations=allocations["OSD-771-LAR"]
    )
    rows.append({"comparison": "ISS-T_per_arm_g_vs_LAR_per_arm_g", "null": "exact LAR within-age label enumeration", **result})
    osd_scores = arm_scores["OSD-513"]
    ref = terminal["pooled_terminal_estimate"]
    sets_b = [s for s in osd_scores.columns if s in defs.compartment and s in ref.index]
    result = P.pattern_concordance_exact(
        ref.loc[sets_b], osd_scores[sets_b], cohorts["osd513"].metadata, allocations=allocations["OSD-513"]
    )
    rows.append({"comparison": "OSD-513_g_vs_pooled_5_mission_terminal", "null": "exact OSD-513 label enumeration", **result})
    table = pd.DataFrame(rows)
    table["sets"] = table["sets"].map(lambda values: "|".join(values))
    table["scope"] = "descriptive_pattern_concordance"
    return table


def duration_context(
    cohorts: Mapping[str, object],
    defs: SetDefinitions,
    settings: Settings,
    arm_scores: Mapping[str, pd.DataFrame],
    l2: Mapping[str, object] | None,
) -> pd.DataFrame:
    rows = []
    osd253 = cohorts["primary"].get("OSD-253")
    targets = (T1, T2, T3)
    if osd253 is not None:
        scores, _ = P.score_arm(osd253, defs.layer_family([T1, T2]), settings.cpm_threshold, min_genes=settings.min_genes)
        for name in [t for t in targets if t in scores]:
            _, strata = mission_effect_from_score(scores[name], osd253.metadata)
            for _, stratum in strata.iterrows():
                rows.append(
                    {"context": "OSD-253_duration", "set": name, "stratum": stratum["stratum"], "quantity": "hedges_g",
                     "estimate": stratum["estimate"], "variance": stratum["variance"], "ci_low": stratum["ci_low"],
                     "ci_high": stratum["ci_high"], "n_flight": stratum["n_treatment"], "n_control": stratum["n_control"]}
                )
            by = strata.set_index("stratum")
            if {"day25", "day75"} <= set(by.index):
                diff = by.loc["day75", "estimate"] - by.loc["day25", "estimate"]
                var = by.loc["day75", "variance"] + by.loc["day25", "variance"]
                rows.append(
                    {"context": "OSD-253_duration", "set": name, "stratum": "day75_minus_day25", "quantity": "hedges_g_difference",
                     "estimate": diff, "variance": var, "ci_low": diff - 1.959964 * np.sqrt(var), "ci_high": diff + 1.959964 * np.sqrt(var)}
                )
    for cohort, data in (("OSD-771", cohorts["iss"]), ("OSD-771-LAR", cohorts["lar"])):
        scores = arm_scores.get(cohort)
        if scores is None:
            continue
        for name in [t for t in targets if t in scores]:
            _, strata = mission_effect_from_score(scores[name], data.metadata)
            for _, stratum in strata.iterrows():
                rows.append(
                    {"context": f"{cohort}_age_stratum", "set": name, "stratum": stratum["stratum"], "quantity": "hedges_g_own_arm_z",
                     "estimate": stratum["estimate"], "variance": stratum["variance"], "ci_low": stratum["ci_low"],
                     "ci_high": stratum["ci_high"], "n_flight": stratum["n_treatment"], "n_control": stratum["n_control"]}
                )
    if l2 is not None:
        oriented, meta40 = l2["oriented"], l2["meta40"]
        for name in [t for t in targets if t in oriented]:
            for age in P.AGES:
                cell = {}
                for arm in P.ARMS:
                    mask = meta40["arm"].eq(arm) & meta40["age"].eq(age)
                    flight = oriented.loc[meta40.index[mask & meta40["condition"].eq("FLT")], name]
                    ground = oriented.loc[meta40.index[mask & meta40["condition"].eq("GC")], name]
                    diff = float(flight.mean() - ground.mean())
                    var = float(flight.var(ddof=1) / len(flight) + ground.var(ddof=1) / len(ground))
                    cell[arm] = (diff, var)
                    rows.append(
                        {"context": "OSD-771_joint_age_stratum", "set": name, "stratum": f"{arm}_{age}",
                         "quantity": "raw_mean_difference_oriented_joint_z", "estimate": diff, "variance": var,
                         "ci_low": diff - 1.959964 * np.sqrt(var), "ci_high": diff + 1.959964 * np.sqrt(var),
                         "n_flight": len(flight), "n_control": len(ground)}
                    )
                est = cell["LAR"][0] - cell["ISS-T"][0]
                var = cell["LAR"][1] + cell["ISS-T"][1]
                rows.append(
                    {"context": "OSD-771_joint_age_stratum", "set": name, "stratum": f"interaction_{age}",
                     "quantity": "LAR_minus_ISS-T_raw_oriented", "estimate": est, "variance": var,
                     "ci_low": est - 1.959964 * np.sqrt(var), "ci_high": est + 1.959964 * np.sqrt(var)}
                )
    table = pd.DataFrame(rows)
    if not table.empty:
        table["scope"] = "descriptive"
    return table


def gene_coverage(cohorts: Mapping[str, object], defs: SetDefinitions, settings: Settings, eligible_both: set[str]) -> pd.DataFrame:
    family = defs.layer_family()
    threshold = settings.cpm_threshold
    elig = {
        "isst": cpm_eligible_genes(cohorts["iss"], threshold),
        "lar": cpm_eligible_genes(cohorts["lar"], threshold),
        "osd513": cpm_eligible_genes(cohorts["osd513"], threshold),
    }
    index = cohorts["iss"].expression.index
    index_513 = cohorts["osd513"].expression.index
    primary_flag = defs.primary_counts.set_index("set")["primary_evaluable"]
    rows = []
    for name, spec in family.items():
        genes = P.set_genes(spec)
        evaluable = {}
        for label, eligible, idx in (("isst", elig["isst"], index), ("lar", elig["lar"], index),
                                     ("osd513", elig["osd513"], index_513), ("layer2", eligible_both, index)):
            _, _, ok, _, _ = P._set_directions(spec, eligible, idx, settings.min_genes)
            evaluable[label] = ok
        rows.append(
            {
                "set": name,
                "scope": defs.scope(name),
                "n_defined": len(genes),
                "n_elig_isst": sum(g in elig["isst"] and g in index for g in genes),
                "n_elig_lar": sum(g in elig["lar"] and g in index for g in genes),
                "n_elig_both": sum(g in eligible_both and g in index for g in genes),
                "n_elig_osd513": sum(g in elig["osd513"] and g in index_513 for g in genes),
                "primary_evaluable": bool(primary_flag.get(name, name == T2 or name in defs.axes)),
                "evaluable_isst": evaluable["isst"],
                "evaluable_lar": evaluable["lar"],
                "evaluable_osd513": evaluable["osd513"],
                "evaluable_layer2": evaluable["layer2"],
            }
        )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- orchestration


def technical_flags(qc: pd.DataFrame, metric: str) -> list[str]:
    if qc.empty or "imbalance_flag" not in qc:
        return []
    mask = qc["metric"].eq(metric) & qc["imbalance_flag"].astype("boolean").fillna(False).astype(bool)
    return qc.loc[mask, "mission"].astype(str).tolist()


def analyze(
    cohorts: Mapping[str, object],
    defs: SetDefinitions,
    config: Mapping[str, object],
    rules: P.PersistenceRules,
    settings: Settings,
    *,
    log: Callable[[str], None] = print,
) -> dict[str, object]:
    """Run every layer on loaded cohorts (no file IO)."""
    metric = str(config["technical_qc"]["primary_metric"])
    l1_missions = {"OSD-771": cohorts["iss"], "OSD-771-LAR": cohorts["lar"], "OSD-513": cohorts["osd513"]}
    qc = technical_qc_audit(l1_missions)
    flagged = technical_flags(qc, metric)
    lar_flag = "OSD-771-LAR" in flagged
    log(f"technical QC ({metric}) flagged: {flagged or 'none'}")
    log("pooled terminal reference (5 primary missions)")
    # score_arm appends D = podocyte - structural_disjoint per mission, so D gets its own reference.
    terminal, terminal_effects = pooled_terminal(
        cohorts["primary"], defs.layer_family(), settings.cpm_threshold, settings.min_genes
    )
    log("layer 1: exact single-cohort enumeration")
    arm_effects, arm_scores, arm_coverage, allocations = layer1(cohorts, defs, settings, flagged)
    log("layer 2: joint OSD-771 HC3 + Freedman-Lane")
    l2 = layer2(cohorts["iss"], cohorts["lar"], defs, terminal, settings, rules, lar_flag)
    log(f"margin sweep: kappas {list(settings.margin_kappas)}")
    sweep = margin_sweep_table(l2, defs, rules, settings, lar_flag)
    log("layer 3: structural-adjusted persistence")
    adjusted, fl_adjusted = layer3(l2, settings, rules, lar_flag, terminal)
    log("layer 4: pattern concordance")
    concordance = layer4(arm_effects, arm_scores, cohorts, defs, terminal, allocations)
    context = duration_context(cohorts, defs, settings, arm_scores, l2)
    coverage = gene_coverage(cohorts, defs, settings, l2["eligible"])

    contrasts = l2["table"]
    contrast_out = contrasts[[c for c in CONTRAST_COLUMNS if c in contrasts.columns]
                             + [c for c in contrasts.columns if c not in CONTRAST_COLUMNS]]
    null_columns = {}
    columns = list(l2["oriented"].columns)
    for index, name in enumerate(columns):
        null_columns[f"interaction::{name}"] = l2["fl"]["interaction"].t_null[:, index]
        if name in (T1, T2, T3):
            null_columns[f"margin::{name}"] = l2["fl"]["margin"].t_null[:, index]
    for key, result in fl_adjusted.items():
        null_columns[f"adjusted_{key}::{P.P1_ADJUSTED}"] = result.t_null[:, 0]
    score_rows = []
    for cohort, scores in arm_scores.items():
        data = l1_missions[cohort]
        long = scores.reset_index(names="animal").melt(id_vars="animal", var_name="set", value_name="score")
        long["layer"] = "single_cohort_own_z"
        long["cohort"] = cohort
        long["condition"] = long["animal"].map(data.metadata["condition"])
        long["block"] = long["animal"].map(data.metadata["block"])
        long["orientation"] = np.nan
        score_rows.append(long)
    joint = l2["joint"].reset_index(names="animal").melt(id_vars="animal", var_name="set", value_name="score")
    joint["layer"] = "joint_osd771_z40"
    joint["cohort"] = "OSD-771_joint"
    joint["condition"] = joint["animal"].map(l2["meta40"]["condition"])
    joint["block"] = joint["animal"].map(l2["meta40"]["block"])
    joint["orientation"] = joint["set"].map(l2["orientation"])
    score_rows.append(joint)
    return {
        "technical_qc": qc,
        "flagged": flagged,
        "lar_technical_flag": lar_flag,
        "terminal": terminal.reset_index(),
        "terminal_effects": terminal_effects,
        "arm_effects": arm_effects,
        "arm_coverage": arm_coverage,
        "contrasts": contrast_out,
        "margin_sweep": sweep,
        "adjusted": adjusted,
        "concordance": concordance,
        "duration_context": context,
        "gene_coverage": coverage,
        "scores": pd.concat(score_rows, ignore_index=True)[
            ["layer", "cohort", "animal", "condition", "block", "set", "score", "orientation"]
        ],
        "null_t": pd.DataFrame(null_columns),
        "cohort_manifest": sample_manifest(l1_missions),
        "families": {
            "compartment_49": list(defs.compartment),
            "key_secondary": [s for s in (T1, T2, T3) if s in columns],
            "exploratory_layer2": [s for s in columns if s in defs.compartment],
            "context_axes": list(defs.axes),
        },
        "allocation_counts": {cohort: int(len(a.masks)) for cohort, a in allocations.items()},
    }


def headline(
    cohorts: Mapping[str, object],
    defs: SetDefinitions,
    config: Mapping[str, object],
    rules: P.PersistenceRules,
    settings: Settings,
) -> pd.DataFrame:
    """T1-T3 Layer-2 rows only (identical to the full run's rows; used for map sensitivity)."""
    metric = str(config["technical_qc"]["primary_metric"])
    qc = technical_qc_audit({"OSD-771-LAR": cohorts["lar"]})
    lar_flag = "OSD-771-LAR" in technical_flags(qc, metric)
    terminal, _ = pooled_terminal(
        cohorts["primary"], defs.layer_family([T1, T2]), settings.cpm_threshold, settings.min_genes
    )
    l2 = layer2(cohorts["iss"], cohorts["lar"], defs, terminal, settings, rules, lar_flag, sets=[T1, T2, T3])
    return l2["table"]


HEADLINE_OUTPUTS = (
    "classification", "theta_isst", "theta_lar", "interaction", "p_interaction_fl", "margin",
    "margin_ci90_low", "margin_ci90_high", "p_persist_rule", "p_decay_rule", "n_genes_joint",
)


def idmap_sensitivity_table(per_map: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    """Long table (output, map, value, verdict) plus a robustness/conservative-verdict row per target."""
    rows = []
    for map_name, table in per_map.items():
        indexed = table.set_index("set")
        for target in [t for t in (T1, T2, T3) if t in indexed.index]:
            verdict = str(indexed.loc[target, "classification"])
            for output in HEADLINE_OUTPUTS:
                rows.append({"output": f"{target}::{output}", "map": map_name,
                             "value": indexed.loc[target, output], "verdict": verdict})
    for target in (T1, T2, T3):
        labels = [str(t.set_index("set").loc[target, "classification"]) for t in per_map.values()
                  if target in set(t["set"])]
        if not labels:
            continue
        verdict, robust = P.conservative_verdict(labels)
        rows.append({"output": f"{target}::robust_to_gene_id_reconstruction", "map": "ALL",
                     "value": robust, "verdict": verdict if robust else f"{verdict} (not robust to gene-ID reconstruction)"})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- IO


def _write(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False)


def write_outputs(outdir: Path, result: Mapping[str, object], manifest: Mapping[str, object]) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    _write(result["cohort_manifest"], outdir / "persistence_cohort_manifest.tsv")
    _write(result["technical_qc"], outdir / "persistence_technical_qc.tsv")
    _write(result["gene_coverage"], outdir / "persistence_gene_coverage.tsv")
    _write(result["arm_effects"], outdir / "persistence_arm_effects.tsv")
    _write(result["contrasts"], outdir / "persistence_contrasts.tsv")
    _write(result["margin_sweep"], outdir / "persistence_margin_sweep.tsv")
    _write(result["adjusted"], outdir / "persistence_adjusted.tsv")
    _write(result["concordance"], outdir / "persistence_pattern_concordance.tsv")
    _write(result["duration_context"], outdir / "persistence_duration_context.tsv")
    _write(result["scores"], outdir / "persistence_scores.tsv")
    _write(result["terminal"], outdir / "persistence_terminal_reference.tsv")
    result["null_t"].to_csv(
        outdir / "persistence_null_t.tsv.gz", sep="\t", index=False, compression="gzip", float_format="%.6g"
    )
    if "power" in result:
        _write(result["power"], outdir / "persistence_power.tsv")
    (outdir / "persistence_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")


def _software() -> dict[str, str]:
    return {
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
    }


def _rules_for(prereg: Mapping[str, object] | None, prereg_path: Path) -> P.PersistenceRules:
    if prereg is None:
        return P.builtin_plan_rules()
    return P.parse_classification_rules(prereg, source=str(prereg_path))


def power_table(settings: Settings, rules: P.PersistenceRules, n_sim: int, n_perm: int) -> pd.DataFrame:
    scenarios = [(theta, rho) for theta in POWER_THETAS for rho in POWER_RHOS]
    return P.simulate_power(
        scenarios, n_per_cell=5, n_sim=n_sim, n_perm=n_perm, seed=settings.seed_for("power_simulation"), rules=rules
    )


def run_power_only(args: argparse.Namespace) -> pd.DataFrame:
    """Power simulation only; loads NO expression data (safe before the preregistration lock)."""
    config = yaml.safe_load(Path(args.config).read_text())
    prereg_path = Path(args.prereg)
    prereg = load_prereg(prereg_path) if prereg_path.exists() else None
    rules = _rules_for(prereg, prereg_path)
    settings = prereg_settings(prereg, config, chunk=args.chunk_size)
    table = power_table(settings, rules, args.power_sims, args.power_permutations)
    outdir = Path(args.results) / OUTPUT_SUBDIR
    outdir.mkdir(parents=True, exist_ok=True)
    _write(table, outdir / "persistence_power.tsv")
    manifest = {
        "analysis": "G2 recovery/persistence power simulation (no expression data loaded)",
        "mode": "power_only",
        "prereg_path": str(prereg_path),
        "prereg_present": prereg is not None,
        "prereg_sha256": sha256(prereg_path) if prereg is not None else None,
        "rules": rules.as_dict(),
        "seed": settings.seed_for("power_simulation"),
        "n_sim": int(args.power_sims),
        "n_perm": int(args.power_permutations),
        "scenarios": {"theta_isst": list(POWER_THETAS), "rho": list(POWER_RHOS)},
        "noise": "N(0,1) per animal plus fixed arm x age block means; n_per_cell 5",
        "software": _software(),
    }
    (outdir / "persistence_power_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return table


def _input_hashes(config: Mapping[str, object], tiers: Path) -> dict[str, str]:
    paths = [Path(tiers)]
    for group in ("primary_missions", "moderator_cohorts"):
        for spec in (config.get(group) or {}).values():
            for key in ("vst", "counts", "runsheet", "metadata", "qc"):
                if key in spec:
                    paths.append(REPO / str(spec[key]))
    paths.append(REPO / str(config["gene_mapping"]["path"]))
    paths.append(REPO / str(config["gene_mapping"]["annotation_fallback"]))
    out = {}
    for path in dict.fromkeys(paths):
        if path.exists():
            key = str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path)
            out[key] = sha256(path)
    return out


def run_idmap_sensitivity(
    args: argparse.Namespace,
    *,
    loader: Callable[[Mapping[str, object], Path], dict] = P.load_persistence_cohorts,
) -> pd.DataFrame:
    guard = prereg_guard(Path(args.prereg), args.allow_uncommitted_prereg)
    config = yaml.safe_load(Path(args.config).read_text())
    prereg = load_prereg(Path(args.prereg))
    rules = _rules_for(prereg, Path(args.prereg))
    settings = prereg_settings(prereg, config, n_perm=args.permutations, n_boot=args.bootstrap,
                               posterior_draws=args.posterior_draws, chunk=args.chunk_size)
    maps = list(args.idmap_sensitivity or [])
    if not maps:
        block = prereg.get("gene_map") or {}
        if not block.get("baseline") or not block.get("runner_up"):
            raise ValueError("--idmap-sensitivity without maps needs gene_map.baseline and runner_up in the prereg")
        maps = [Path(block["baseline"]), Path(block["runner_up"])]
    deviations = []
    for map_path in maps:
        deviations += check_gene_map(gene_map_record(map_path, prereg), args.allow_uncommitted_prereg)
    args.idmap_sensitivity = maps
    per_map = {}
    for map_path in maps:
        config_map = apply_gene_map_override(config, map_path)
        cohorts = loader(config_map, REPO)
        defs = build_set_definitions(args.tiers, cohorts["primary"], config_map["primary_family"],
                                     settings.cpm_threshold, settings.min_genes)
        per_map[str(map_path)] = headline(cohorts, defs, config_map, rules, settings)
    table = idmap_sensitivity_table(per_map)
    outdir = Path(args.results) / OUTPUT_SUBDIR
    outdir.mkdir(parents=True, exist_ok=True)
    _write(table, outdir / "idmap_sensitivity.tsv")
    manifest = {**guard, "mode": "idmap_sensitivity", "maps": {str(m): sha256(Path(apply_gene_map_override(config, m)["gene_mapping"]["path"])) for m in args.idmap_sensitivity},
                "map_records": [gene_map_record(m, prereg) for m in maps],
                "rules": rules.as_dict(), "n_perm": settings.n_perm, "deviations": list(settings.deviations) + deviations,
                "software": _software()}
    (outdir / "idmap_sensitivity_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return table


PREREG_KEYS_USED = (
    "analysis_id", "status", "lock_date", "seed_offsets", "gene_map", "targets", "eligibility", "contrasts",
    "inference", "classification", "margin_sweep",
)
PREREG_KEYS_RECORDED = ("multiplicity", "reporting_language")


def run_validate_prereg(args: argparse.Namespace) -> dict[str, object]:
    """Parse and check the preregistration exactly as a real run would; loads NO expression data."""
    guard = prereg_guard(Path(args.prereg), args.allow_uncommitted_prereg)
    prereg = load_prereg(Path(args.prereg))
    rules = _rules_for(prereg, Path(args.prereg))
    config = apply_gene_map_override(yaml.safe_load(Path(args.config).read_text()), args.gene_map)
    settings = prereg_settings(prereg, config, n_perm=args.permutations, n_boot=args.bootstrap,
                               posterior_draws=args.posterior_draws, chunk=args.chunk_size,
                               sweep_all_sets=args.sweep_all_sets)
    analysis_map = gene_map_record(config["gene_mapping"]["path"], prereg)
    block = prereg.get("gene_map") or {}
    declared = {role: gene_map_record(block[role], prereg) for role in ("baseline", "runner_up") if block.get(role)}
    deviations = list(settings.deviations)
    for record in (analysis_map, *declared.values()):
        deviations += check_gene_map(record, args.allow_uncommitted_prereg)
    summary = {
        **guard,
        "mode": "validate_prereg_only (no expression data loaded)",
        "analysis_id": prereg.get("analysis_id"),
        "status": prereg.get("status"),
        "lock_date": prereg.get("lock_date"),
        "rules": rules.as_dict(),
        "seeds": {name: settings.seed_for(name) for name in settings.offsets},
        "n_permutations_freedman_lane": settings.n_perm,
        "n_bootstrap": settings.n_boot,
        "posterior_draws": settings.posterior_draws,
        "min_genes": settings.min_genes,
        "cpm_threshold": settings.cpm_threshold,
        "targets": {"primary": settings.primary_target, "key_secondary": list(settings.key_secondary)},
        "margin_sweep": {"kappas": list(settings.margin_kappas), "primary_kappa": settings.primary_margin_kappa,
                         "applies_to": list(settings.sweep_sets)},
        "analysis_gene_map": analysis_map,
        "declared_gene_maps": declared,
        "deviations": deviations,
        "prereg_keys_used": [k for k in PREREG_KEYS_USED if k in prereg],
        "prereg_keys_recorded_only": [k for k in PREREG_KEYS_RECORDED if k in prereg],
        "prereg_keys_not_used": [k for k in prereg if k not in PREREG_KEYS_USED + PREREG_KEYS_RECORDED],
    }
    print(json.dumps(summary, indent=2, default=str))
    return summary


def run(args: argparse.Namespace, *, loader: Callable[[Mapping[str, object], Path], dict] = P.load_persistence_cohorts) -> dict[str, object]:
    guard = prereg_guard(Path(args.prereg), args.allow_uncommitted_prereg)
    prereg = load_prereg(Path(args.prereg))
    rules = _rules_for(prereg, Path(args.prereg))
    config = apply_gene_map_override(yaml.safe_load(Path(args.config).read_text()), args.gene_map)
    settings = prereg_settings(prereg, config, n_perm=args.permutations, n_boot=args.bootstrap,
                               posterior_draws=args.posterior_draws, chunk=args.chunk_size,
                               sweep_all_sets=args.sweep_all_sets)
    map_record = gene_map_record(config["gene_mapping"]["path"], prereg)
    map_deviations = check_gene_map(map_record, args.allow_uncommitted_prereg)
    cohorts = loader(config, REPO)
    defs = build_set_definitions(args.tiers, cohorts["primary"], config["primary_family"],
                                 settings.cpm_threshold, settings.min_genes)
    result = analyze(cohorts, defs, config, rules, settings)
    if args.with_power:
        result["power"] = power_table(settings, rules, args.power_sims, args.power_permutations)
    manifest = {
        "analysis_id": prereg.get("analysis_id"),
        "analysis": "G2 recovery/persistence of the terminal-flight shift",
        "mode": "full",
        **guard,
        "prereg_lock_date": prereg.get("lock_date"),
        "prereg_status": prereg.get("status"),
        "config": str(args.config),
        "config_sha256": sha256(Path(args.config)),
        "gene_map": map_record,
        "gene_map_override": args.gene_map is not None,
        "tiers": str(args.tiers),
        "rules": rules.as_dict(),
        "seeds": {name: settings.seed_for(name) for name in settings.offsets},
        "n_permutations_freedman_lane": settings.n_perm,
        "margin_sweep": {
            "kappas": list(settings.margin_kappas),
            "primary_kappa": settings.primary_margin_kappa,
            "role": "secondary_descriptive; only the primary kappa sets the classification",
            "sets": "all layer-2 sets" if settings.sweep_all_sets else list(settings.sweep_sets),
        },
        "n_bootstrap": settings.n_boot,
        "posterior_draws": settings.posterior_draws,
        "min_genes": settings.min_genes,
        "cpm_threshold": settings.cpm_threshold,
        "permutation_stream": "one within-arm x age-block permutation stream (seed+700, chunk) shared by every set, both contrasts and layer 3",
        "orientation": "layer-2/3 effects multiplied by the sign of the pooled 5-mission terminal estimate",
        "lar_technical_flag": result["lar_technical_flag"],
        "technical_flags": result["flagged"],
        "families": result["families"],
        "allocation_counts": result["allocation_counts"],
        "multiplicity": {
            "primary": "single target, no correction",
            "key_secondary": "Holm over {structural_disjoint, D} per margin direction; interaction max-T over {T1,T2,T3}",
            "exploratory": "interaction max-T + BH over evaluable 49-set members; labels descriptive",
        },
        "deviations": list(settings.deviations) + map_deviations,
        "prereg_multiplicity": prereg.get("multiplicity"),
        "reporting_language": prereg.get("reporting_language"),
        "not_done": ["REML/mKH across recovery cohorts (k=2, df=1)", "confirmatory 7-cohort meta-regression"],
        "input_sha256": _input_hashes(config, Path(args.tiers)),
        "software": _software(),
        "interpretation_boundary": (
            "recovery aliased with 32 vs 53-56 d flight, on-ISS vs on-Earth euthanasia/handling, re-entry and "
            "landing; theta_ISS-T for the selected podocyte set is winner's-curse inflated, biasing rho low"
        ),
    }
    outdir = Path(args.results) / OUTPUT_SUBDIR
    write_outputs(outdir, result, manifest)
    show = result["contrasts"].set_index("set").loc[[s for s in (T1, T2, T3) if s in set(result["contrasts"]["set"])]]
    print(show[["theta_isst", "theta_lar", "margin", "margin_ci90_low", "p_persist_rule", "p_decay_rule", "classification"]].to_string())
    print(f"\nWrote {outdir}")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--prereg", type=Path, default=DEFAULT_PREREG)
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS, help=f"run directory; outputs go to <results>/{OUTPUT_SUBDIR}/")
    parser.add_argument("--gene-map", type=Path, default=None, help="alternative id_map.tsv (replaces gene_mapping.path in memory)")
    parser.add_argument("--allow-uncommitted-prereg", action="store_true", help="tests only")
    parser.add_argument("--power-only", action="store_true", help="power simulation only; no expression data is loaded")
    parser.add_argument("--with-power", action="store_true", help="also run the power simulation in a full run")
    parser.add_argument("--idmap-sensitivity", nargs="*", type=Path, default=None, metavar="MAP",
                        help="headline T1-T3 under each gene-ID map (APPENDIX A item 3); "
                             "no MAP = the prereg's gene_map baseline and runner_up")
    parser.add_argument("--validate-prereg-only", action="store_true",
                        help="parse and check the preregistration as a real run would; loads no expression data")
    parser.add_argument("--permutations", type=int, default=None, help="Freedman-Lane B (default: preregistered)")
    parser.add_argument("--bootstrap", type=int, default=None, help="bootstrap B (default: preregistered)")
    parser.add_argument("--posterior-draws", type=int, default=None)
    parser.add_argument("--chunk-size", type=int, default=2048)
    parser.add_argument("--sweep-all-sets", action="store_true", help="margin sweep over every layer-2 set (default T1-T3)")
    parser.add_argument("--power-sims", type=int, default=4000)
    parser.add_argument("--power-permutations", type=int, default=999)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.power_only:
        table = run_power_only(args)
        print(table.to_string(index=False))
        return 0
    if args.validate_prereg_only:
        run_validate_prereg(args)
        return 0
    if args.idmap_sensitivity is not None:
        print(run_idmap_sensitivity(args).to_string(index=False))
        return 0
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
