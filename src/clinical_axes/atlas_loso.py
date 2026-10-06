"""Leave-one-atlas-source-study-out reconstruction of the frozen compartment marker tiers.

Item 8 of the podocyte sensitivity programme asks whether the 157-gene podocyte
high-specificity definition (and the 49-set family it sits in) depends on any
single Mouse Kidney Atlas source study.  The atlas is flight-label blind, so
rebuilding the marker tiers from the source-study pseudobulks never touches
mission data; this module contains only that flight-blind part:

* :func:`compartment_expression_from_pseudobulk` reproduces, from the
  ``atlas_pseudobulk_counts`` / ``atlas_pseudobulk_sample_metadata`` outputs of
  ``scripts/subtype_reference/03_atlas_pseudobulk.py``, the compartment
  expression table ``atlas_compartment_expression.tsv.gz`` (same arithmetic,
  same schema), optionally with whole source studies removed;
* :func:`rebuild_tiers` / :func:`leave_one_study_out_tiers` feed that table to the
  frozen ``build_marker_tiers`` with the frozen thresholds (never tuned) and apply
  the frozen atlas QC rule (>= 3 source studies for every required compartment);
* set-comparison, effective-stringency, gate and preregistration-guard helpers.

Scoring the rebuilt sets against flight data lives in
``scripts/clinical_axes/run_atlas_loso_markers.py`` and may run only after the
preregistration is committed.
"""

from __future__ import annotations

import ast
import copy
from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path
import re
import subprocess
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from src.v13.compartment_adversarial_audit import build_marker_tiers

EXPRESSION_COLUMNS = (
    "gene_symbol",
    "compartment",
    "mean_cpm",
    "median_cpm",
    "source_study_detection_fraction",
    "n_source_studies",
    "reference_label",
)
PODOCYTE_COMPARTMENT = "podocyte"
PODOCYTE_SET = "podocyte__high_specificity"
STRUCTURAL_SET = "broad_structural_scaffold_control__all"
MINIMUM_SET_GENES = 8

# Seeds are the parent seed (config) plus an offset (plan section 0).
SEED_OFFSET_REFERENCE = 100  # primary compartment-context seed; reproduces the primary family
SEED_OFFSET_LOSO = 800  # omitted study j (alphabetical by source_study_id) uses +800+j

RETENTION_RATIO_MIN = 0.50
GATE_RULE_TEXT = (
    "sign(estimate) == sign(reference) and estimate / reference >= {minimum:.2f}"
)

STATUS_NOT_EVALUABLE = "NOT_EVALUABLE"
STATUS_EVALUABLE = "EVALUABLE"
STATUS_NOT_SCORED = "NOT_SCORED"


# --- preregistration guard -----------------------------------------------------


class PreregistrationError(RuntimeError):
    """The preregistration YAML is missing or has uncommitted changes."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head(repo: Path) -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=repo,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return None


def _git_lines(repo: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args], cwd=repo, capture_output=True, text=True
    )


def check_preregistration(
    prereg: Path, *, repo: Path, allow_uncommitted: bool = False
) -> dict[str, Any]:
    """Refuse unless the preregistration YAML exists and is committed without local changes.

    "Uncommitted" means a non-empty ``git status --porcelain -- <path>`` (modified,
    staged or untracked), a path that git does not track (for example one that is
    ignored), or a path outside the repository (commit state cannot be verified).
    ``allow_uncommitted`` (tests only) waives that refusal; a missing file is always
    refused.  Returns the record that goes into the run manifest (path, sha256, git
    HEAD, whether the file was uncommitted).
    """

    repo = Path(repo)
    path = Path(prereg)
    if not path.is_absolute():
        path = repo / path
    if not path.is_file():
        raise PreregistrationError(f"preregistration YAML is missing: {path}")
    resolved_repo = repo.resolve()
    inside = True
    try:
        shown = str(path.resolve().relative_to(resolved_repo))
    except ValueError:
        inside = False
        shown = str(path)
    if inside:
        status = _git_lines(repo, "status", "--porcelain", "--", str(path))
        if status.returncode != 0:
            raise PreregistrationError(
                f"could not determine git status of {path}: {status.stderr.strip()}"
            )
        tracked = _git_lines(repo, "ls-files", "--error-unmatch", "--", str(path))
        dirty = bool(status.stdout.strip()) or tracked.returncode != 0
        detail = status.stdout.strip() or ("not tracked by git" if dirty else "")
    else:
        dirty = True
        detail = "outside the repository; commit state cannot be verified"
    if dirty and not allow_uncommitted:
        raise PreregistrationError(
            f"preregistration YAML is not committed cleanly ({detail}); commit it before "
            "any real-data run (or pass --allow-uncommitted-prereg in tests)"
        )
    return {
        "path": shown,
        "sha256": sha256_file(path),
        "git_head": git_head(repo),
        "uncommitted_changes": dirty,
        "allow_uncommitted_prereg": bool(allow_uncommitted),
    }


# --- pseudobulk -> compartment expression --------------------------------------


def _as_bool(values: pd.Series) -> pd.Series:
    if values.dtype == bool:
        return values
    return values.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def load_pseudobulk(
    counts_path: Path, sample_meta_path: Path
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the ``03_atlas_pseudobulk.py`` count matrix (genes x samples) and sample table."""

    counts = pd.read_csv(counts_path, sep="\t", index_col=0)
    sample_meta = pd.read_csv(sample_meta_path, sep="\t")
    return counts, sample_meta


def _kept_samples(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
    exclude_sources: Iterable[str] = (),
) -> pd.DataFrame:
    """Included pseudobulks (>= min cells) minus the excluded source studies, in table order."""

    required = {"sample_id", "source_study_id", "compartment"}
    missing = required - set(sample_meta.columns)
    if missing:
        raise ValueError(f"sample metadata lacks columns {sorted(missing)}")
    meta = sample_meta.copy()
    if "included" in meta.columns:
        meta = meta[_as_bool(meta["included"])]
    exclude = {str(value) for value in exclude_sources}
    unknown = exclude - set(sample_meta["source_study_id"].astype(str))
    if unknown:
        raise ValueError(f"cannot exclude unknown source studies: {sorted(unknown)}")
    meta = meta[~meta["source_study_id"].astype(str).isin(exclude)]
    absent = [sample for sample in meta["sample_id"] if sample not in counts.columns]
    if absent:
        raise ValueError(f"pseudobulk counts lack samples: {absent[:5]}")
    return meta.reset_index(drop=True)


def compartment_expression_from_pseudobulk(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
    *,
    min_cpm: float = 1.0,
    exclude_sources: Iterable[str] = (),
    reference_label: str = "mouse_kidney_atlas",
) -> pd.DataFrame:
    """Compartment expression table, identical in arithmetic to ``03_atlas_pseudobulk.py``.

    ``counts`` is genes x ``Origin::compartment`` pseudobulks (integer counts over the
    full gene universe, so library size is the column sum); ``sample_meta`` is the
    pseudobulk sample table (``sample_id``, ``source_study_id``, ``compartment`` and,
    optionally, ``included``).  Included pseudobulks of ``exclude_sources`` are
    dropped.  CPM is computed per pseudobulk column; per compartment the output holds
    ``mean_cpm`` and ``median_cpm`` over its pseudobulks, ``source_study_detection_fraction``
    = mean(CPM >= ``min_cpm``) and ``n_source_studies`` (its number of pseudobulks, one
    per source study).  Compartments are emitted alphabetically, genes in input order.
    """

    kept = _kept_samples(counts, sample_meta, exclude_sources)
    if kept.empty:
        raise ValueError("no pseudobulks remain after exclusion")
    dense = np.ascontiguousarray(
        counts.loc[:, kept["sample_id"]].to_numpy(dtype=np.float64).T
    )
    library = dense.sum(axis=1)
    cpm = np.divide(
        dense * 1e6,
        library[:, None],
        out=np.zeros_like(dense),
        where=library[:, None] > 0,
    )
    genes = pd.Index(counts.index.astype(str), name="gene_symbol")
    frames: list[pd.DataFrame] = []
    for compartment, indexes in kept.groupby("compartment", observed=True).groups.items():
        idx = np.asarray(list(indexes), dtype=int)
        frames.append(
            pd.DataFrame(
                {
                    "gene_symbol": genes,
                    "compartment": compartment,
                    "mean_cpm": cpm[idx, :].mean(axis=0),
                    "median_cpm": np.median(cpm[idx, :], axis=0),
                    "source_study_detection_fraction": (
                        cpm[idx, :] >= min_cpm
                    ).mean(axis=0),
                    "n_source_studies": len(idx),
                    "reference_label": reference_label,
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def studies_by_compartment(
    sample_meta: pd.DataFrame,
    exclude_sources: Iterable[str] = (),
) -> dict[str, int]:
    """Number of distinct source studies with an included pseudobulk, per compartment."""

    meta = sample_meta.copy()
    if "included" in meta.columns:
        meta = meta[_as_bool(meta["included"])]
    exclude = {str(value) for value in exclude_sources}
    meta = meta[~meta["source_study_id"].astype(str).isin(exclude)]
    return {
        str(compartment): int(sub["source_study_id"].nunique())
        for compartment, sub in meta.groupby("compartment", observed=True)
    }


def atlas_qc(
    sample_meta: pd.DataFrame,
    qc_cfg: Mapping[str, Any],
    exclude_sources: Iterable[str] = (),
) -> dict[str, Any]:
    """Frozen coverage QC: every required compartment needs >= N source studies.

    ``qc_cfg`` is ``reference_builder.atlas_qc`` of the frozen DCT reference config.
    The cell-level mapping-fraction rule is a property of the whole atlas and is not
    affected by dropping a study's pseudobulks.
    """

    studies = studies_by_compartment(sample_meta, exclude_sources)
    minimum = int(qc_cfg["minimum_source_studies_per_compartment"])
    required = list(map(str, qc_cfg["required_compartments"]))
    failed = {
        compartment: int(studies.get(compartment, 0))
        for compartment in required
        if int(studies.get(compartment, 0)) < minimum
    }
    return {
        "studies_by_compartment": studies,
        "minimum_source_studies": minimum,
        "required_compartments": required,
        "failed_compartments": failed,
        "qc_pass": not failed,
    }


# --- marker tiers --------------------------------------------------------------

_FAMILY_DISAGREEMENT = re.compile(r"missing=(\[[^\]]*\]), unexpected=(\[[^\]]*\])")


def _parse_family_disagreement(message: str) -> tuple[list[str], list[str]] | None:
    found = _FAMILY_DISAGREEMENT.search(message)
    if found is None:
        return None
    return ast.literal_eval(found.group(1)), ast.literal_eval(found.group(2))


def build_tiers_tolerant(
    atlas: pd.DataFrame,
    structural_terms: Mapping[str, Sequence[str]],
    audit_cfg: Mapping[str, Any],
) -> tuple[pd.DataFrame | None, dict[str, Any]]:
    """``build_marker_tiers`` that tolerates sets that become empty in a rebuild.

    The frozen builder insists that the produced sets equal the 51 frozen names.  In
    a leave-one-study-out rebuild a small set (for example the 7-gene
    DCT2/CNT high-specificity set) can legitimately become empty; those sets are then
    reported as missing instead of aborting the rebuild.  If a mapped atlas
    compartment is absent altogether the rebuild is NOT_EVALUABLE and no tiers are
    returned.  Thresholds are never altered.
    """

    mapping = audit_cfg["compartment_mapping"]
    absent = sorted(set(mapping) - set(atlas["compartment"].astype(str)))
    if absent:
        return None, {
            "status": STATUS_NOT_EVALUABLE,
            "reason": f"atlas compartments absent after exclusion: {absent}",
            "missing_sets": [],
        }
    try:
        tiers = build_marker_tiers(atlas, structural_terms, audit_cfg)
        return tiers, {"status": "OK", "reason": "", "missing_sets": []}
    except ValueError as exc:
        parsed = _parse_family_disagreement(str(exc))
        if parsed is None:
            raise
        missing, unexpected = parsed
        if unexpected:
            raise
    relaxed = copy.deepcopy(dict(audit_cfg))
    relaxed["set_test"] = dict(relaxed["set_test"])
    relaxed["set_test"]["primary_family"] = [
        name for name in audit_cfg["set_test"]["primary_family"] if name not in set(missing)
    ]
    tiers = build_marker_tiers(atlas, structural_terms, relaxed)
    return tiers, {
        "status": "OK",
        "reason": f"{len(missing)} frozen set(s) empty in this rebuild",
        "missing_sets": sorted(missing),
    }


def set_membership(tiers: pd.DataFrame) -> dict[str, frozenset[str]]:
    """gene_set -> frozenset of genes (rows flagged ``final_for_testing`` when present)."""

    table = tiers
    if "final_for_testing" in table.columns:
        table = table[table["final_for_testing"].astype(bool)]
    return {
        str(name): frozenset(sub["gene_symbol"].dropna().astype(str))
        for name, sub in table.groupby("gene_set", sort=False)
    }


def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    a, b = set(left), set(right)
    union = a | b
    if not union:
        return float("nan")
    return len(a & b) / len(union)


def compare_sets(
    full_sets: Mapping[str, frozenset[str]],
    rebuilt_sets: Mapping[str, frozenset[str]],
    omitted_origin: str,
) -> pd.DataFrame:
    """Per-set size / overlap / Jaccard of a rebuild against the full-atlas sets."""

    rows = []
    for name in full_sets:
        full = full_sets[name]
        rebuilt = rebuilt_sets.get(name, frozenset())
        shared = full & rebuilt
        rows.append(
            {
                "omitted_origin": omitted_origin,
                "gene_set": name,
                "n_full": len(full),
                "n_rebuilt": len(rebuilt),
                "n_shared": len(shared),
                "n_gained": len(rebuilt - full),
                "n_lost": len(full - rebuilt),
                "jaccard_vs_full": jaccard(full, rebuilt),
            }
        )
    return pd.DataFrame(rows)


def tiers_equal_setwise(
    left: pd.DataFrame, right: pd.DataFrame
) -> dict[str, Any]:
    """Set-for-set equality of two tier tables (names and members, not row order)."""

    a, b = set_membership(left), set_membership(right)
    differing = sorted(name for name in set(a) & set(b) if a[name] != b[name])
    return {
        "identical": set(a) == set(b) and not differing,
        "only_in_left": sorted(set(a) - set(b)),
        "only_in_right": sorted(set(b) - set(a)),
        "differing_sets": differing,
    }


def minimum_studies_for_fraction(n_studies: int, fraction: float) -> int:
    """Smallest number of studies that satisfies ``detection_fraction >= fraction``."""

    if n_studies < 1:
        return 0
    return int(math.ceil(fraction * n_studies - 1e-9))


def effective_detection_rule(n_studies: int, fraction: float) -> str:
    """Human-readable effective rule, e.g. ``3-of-4`` (0.75 with four studies)."""

    return f"{minimum_studies_for_fraction(n_studies, fraction)}-of-{n_studies}"


def retention_gate(
    estimate: float | None,
    reference: float | None,
    *,
    retention_min: float = RETENTION_RATIO_MIN,
) -> dict[str, Any]:
    """Default gate: same sign as the reference AND estimate / reference >= ``retention_min``.

    A zero or non-finite reference (or estimate) gives NA values rather than a verdict.
    """

    na = {"retention_ratio": float("nan"), "direction_retained": pd.NA, "gate_pass": pd.NA}
    if estimate is None or reference is None:
        return na
    if not (math.isfinite(estimate) and math.isfinite(reference)) or reference == 0:
        return na
    ratio = estimate / reference
    direction = bool(estimate != 0 and math.copysign(1.0, estimate) == math.copysign(1.0, reference))
    return {
        "retention_ratio": float(ratio),
        "direction_retained": direction,
        "gate_pass": bool(direction and ratio >= retention_min),
    }


def loso_seed(parent_seed: int, index: int) -> int:
    """Permutation seed for the ``index``-th omitted study (alphabetical): parent + 800 + index."""

    return int(parent_seed) + SEED_OFFSET_LOSO + int(index)


# --- leave-one-study-out tier rebuild ------------------------------------------


@dataclass(frozen=True)
class TierRebuild:
    """One atlas rebuild (``omitted_origin`` is ``None`` for the no-exclusion baseline)."""

    omitted_origin: str | None
    index: int | None
    podocyte_contributing: bool | None
    qc: dict[str, Any]
    status: str
    reason: str
    tiers: pd.DataFrame | None
    missing_sets: tuple[str, ...] = ()
    expression: pd.DataFrame | None = field(default=None, repr=False)


def podocyte_contributes(
    sample_meta: pd.DataFrame, origin: str, compartment: str = PODOCYTE_COMPARTMENT
) -> bool:
    """True when ``origin`` has an included pseudobulk for the compartment."""

    meta = sample_meta.copy()
    if "included" in meta.columns:
        meta = meta[_as_bool(meta["included"])]
    return bool(
        (
            (meta["source_study_id"].astype(str) == str(origin))
            & (meta["compartment"].astype(str) == compartment)
        ).any()
    )


def rebuild_tiers(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
    structural_terms: Mapping[str, Sequence[str]],
    audit_cfg: Mapping[str, Any],
    qc_cfg: Mapping[str, Any],
    *,
    omitted: str | None = None,
    index: int | None = None,
    min_cpm: float = 1.0,
    reference_label: str = "mouse_kidney_atlas",
    keep_expression: bool = False,
) -> TierRebuild:
    """Rebuild the marker tiers with ``omitted`` source study removed (``None`` = baseline)."""

    exclude = () if omitted is None else (omitted,)
    qc = atlas_qc(sample_meta, qc_cfg, exclude)
    contributing = (
        None if omitted is None else podocyte_contributes(sample_meta, omitted)
    )
    try:
        atlas = compartment_expression_from_pseudobulk(
            counts,
            sample_meta,
            min_cpm=min_cpm,
            exclude_sources=exclude,
            reference_label=reference_label,
        )
    except ValueError as exc:
        return TierRebuild(omitted, index, contributing, qc, STATUS_NOT_EVALUABLE, str(exc), None)
    tiers, info = build_tiers_tolerant(atlas, structural_terms, audit_cfg)
    return TierRebuild(
        omitted_origin=omitted,
        index=index,
        podocyte_contributing=contributing,
        qc=qc,
        status=info["status"],
        reason=info["reason"],
        tiers=tiers,
        missing_sets=tuple(info["missing_sets"]),
        expression=atlas if keep_expression else None,
    )


def leave_one_study_out_tiers(
    counts: pd.DataFrame,
    sample_meta: pd.DataFrame,
    structural_terms: Mapping[str, Sequence[str]],
    audit_cfg: Mapping[str, Any],
    qc_cfg: Mapping[str, Any],
    *,
    origins: Sequence[str] | None = None,
    min_cpm: float = 1.0,
    reference_label: str = "mouse_kidney_atlas",
) -> list[TierRebuild]:
    """One :class:`TierRebuild` per source study (alphabetical; index = position in that order)."""

    all_origins = sorted(sample_meta["source_study_id"].astype(str).unique())
    chosen = all_origins if origins is None else [str(o) for o in origins]
    unknown = sorted(set(chosen) - set(all_origins))
    if unknown:
        raise ValueError(f"unknown source studies: {unknown}")
    return [
        rebuild_tiers(
            counts,
            sample_meta,
            structural_terms,
            audit_cfg,
            qc_cfg,
            omitted=origin,
            index=all_origins.index(origin),
            min_cpm=min_cpm,
            reference_label=reference_label,
        )
        for origin in chosen
    ]


def tier_summary_row(
    rebuild: TierRebuild,
    full_sets: Mapping[str, frozenset[str]],
    *,
    hs_detection_fraction: float,
    minimum_podocyte_genes: int = MINIMUM_SET_GENES,
    podocyte_compartment: str = PODOCYTE_COMPARTMENT,
    podocyte_set: str = PODOCYTE_SET,
) -> dict[str, Any]:
    """Flight-blind summary of a rebuild: QC, podocyte set size, overlap with the full set."""

    sets = set_membership(rebuild.tiers) if rebuild.tiers is not None else {}
    full_podocyte = full_sets.get(podocyte_set, frozenset())
    podocyte = sets.get(podocyte_set, frozenset())
    n_podocyte_studies = int(
        rebuild.qc["studies_by_compartment"].get(podocyte_compartment, 0)
    )
    evaluable_tiers = (
        rebuild.status != STATUS_NOT_EVALUABLE
        and len(podocyte) >= minimum_podocyte_genes
    )
    if rebuild.status == STATUS_NOT_EVALUABLE:
        reason = rebuild.reason
    elif len(podocyte) < minimum_podocyte_genes:
        reason = (
            f"podocyte high-specificity set has {len(podocyte)} genes "
            f"(< {minimum_podocyte_genes})"
        )
    else:
        reason = rebuild.reason
    contributing = rebuild.podocyte_contributing
    return {
        "omitted_origin": rebuild.omitted_origin,
        "omit_index": rebuild.index,
        "podocyte_contributing": (
            pd.NA if contributing is None else ("yes" if contributing else "no")
        ),
        "qc_pass": bool(rebuild.qc["qc_pass"]),
        "qc_failed_compartments": "|".join(
            f"{name}:{n}" for name, n in rebuild.qc["failed_compartments"].items()
        ),
        "n_source_studies_podocyte": n_podocyte_studies,
        "effective_podocyte_hs_rule": effective_detection_rule(
            n_podocyte_studies, hs_detection_fraction
        ),
        "n_sets_defined": len(sets),
        "n_sets_ge8_genes": int(sum(len(genes) >= MINIMUM_SET_GENES for genes in sets.values())),
        "n_sets_missing_from_rebuild": len(rebuild.missing_sets),
        "n_full_podocyte_hs": len(full_podocyte),
        "n_podocyte_hs": len(podocyte),
        "n_retained_from_full": len(podocyte & full_podocyte),
        "jaccard_vs_full": jaccard(podocyte, full_podocyte) if sets else float("nan"),
        "tiers_evaluable": bool(evaluable_tiers),
        "tier_status_reason": reason,
    }
