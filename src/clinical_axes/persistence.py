"""Recovery/persistence (G2): does the terminal-flight shift persist after live return?

Preregistered in ``config/clinical_axes_recovery_persistence.yaml`` (plan §3 WP-P, §4.2).
Pure functions only; the CLI driver is ``scripts/clinical_axes/run_recovery_persistence.py``.

Layers
------
1. Exact single-arm / single-cohort tests (ISS-T, LAR, OSD-513): each cohort keeps its
   own CPM eligibility and its own gene z-scores (frozen rule); per-stratum Hedges g is
   fixed-effect combined; p-values come from complete enumeration of the within-block
   label allocations (63,504 per OSD-771 arm, 48,620 for OSD-513).
2. Joint OSD-771 model (40 animals, raw score units): genes eligible in BOTH arms,
   z-scored once over all 40 animals, ``score ~ 0 + block(arm x age) + F_T + F_L``,
   OLS with HC3; interaction ``I = theta_L - theta_T`` and margin
   ``M = theta_L - 0.5 theta_T``; Freedman--Lane residual permutation within the four
   arm x age blocks with one permutation shared by every set (valid max-T).
3. Structural-adjusted persistence (P1-type): arm-specific structural slopes added to
   the nuisance design.
4. Pattern concordance across sets with exact label-enumeration nulls.

All Layer-2/3 inputs are oriented by the sign of the pooled five-mission terminal
estimate, so ``margin > 0`` always means "more than half of the terminal shift remains".

Inference conventions (shared with the rest of ``src.clinical_axes``): plus-one
Monte-Carlo p-values, complete enumeration where it is cheap, HC3 without a
finite-sample df correction (matching statsmodels), and per-arm residual df for
t quantiles with Satterthwaite df for contrasts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import itertools
import math
import re
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .analysis import score_mission_axes
from .data import MissionData, cpm_eligible_genes, load_moderator_missions, load_primary_missions
from .statistics import _block_hedges_g_batch, score_signed_axis


PODOCYTE_SET = "podocyte__high_specificity"
STRUCTURAL_FULL_SET = "broad_structural_scaffold_control__all"
STRUCTURAL_DISJOINT = "structural_disjoint"
D_CONTRAST = "D_podocyte_minus_structural_disjoint"
P1_ADJUSTED = "P1_type_adjusted_podocyte"
KEY_TARGETS = (PODOCYTE_SET, STRUCTURAL_DISJOINT, D_CONTRAST)

ARMS = ("ISS-T", "LAR")
AGES = ("YNG", "OLD")
BLOCK_COLUMNS = tuple(f"{arm}_{age}" for arm in ARMS for age in AGES)

NOT_EVALUABLE_TECHNICAL = "NOT_EVALUABLE_TECHNICAL"
NOT_EVALUABLE_NO_TERMINAL_EFFECT = "NOT_EVALUABLE_NO_TERMINAL_EFFECT"
PERSISTS = "PERSISTS"
DECAYS = "DECAYS"
INCONCLUSIVE = "INCONCLUSIVE"
REVERSED = "REVERSED"
LABELS = (
    NOT_EVALUABLE_TECHNICAL,
    NOT_EVALUABLE_NO_TERMINAL_EFFECT,
    PERSISTS,
    DECAYS,
    INCONCLUSIVE,
)
DECAYS_REVERSED = f"{DECAYS}:{REVERSED}"

# Verbatim plan §4.2 draft. Used only when no preregistration YAML exists yet
# (``--power-only`` before the lock); every real-data run reads the committed YAML.
PLAN_4_2_CLASSIFICATION: dict[str, str] = {
    NOT_EVALUABLE_TECHNICAL: "LAR primary coverage metric flagged by the frozen QC rule",
    NOT_EVALUABLE_NO_TERMINAL_EFFECT: (
        "HC3 95% CI of theta_ISS-T includes 0, or sign differs from pooled terminal"
    ),
    PERSISTS: "FL one-sided p(margin>0) <= 0.05 AND HC3 90% CI lower(margin) > 0",
    DECAYS: (
        "FL one-sided p(margin<0) <= 0.05 AND HC3 90% CI upper(margin) < 0  "
        "(sublabel REVERSED if theta_LAR 95% CI upper < 0)"
    ),
    INCONCLUSIVE: "otherwise, including FL/HC3 disagreement",
}
PLAN_4_2_CONTRASTS: dict[str, str] = {
    "interaction": "theta_LAR - theta_ISS-T",
    "margin": "theta_LAR - 0.5*theta_ISS-T",
}

_TIE_RTOL = 1e-10


class CohortContractError(ValueError):
    """A persistence cohort does not have the frozen design."""


# --------------------------------------------------------------------------- cohorts


def _check_arm(data: MissionData, arm: str) -> None:
    meta = data.metadata
    counts = meta["condition"].value_counts()
    if int(counts.get("FLT", 0)) != 10 or int(counts.get("GC", 0)) != 10:
        raise CohortContractError(f"{arm}: expected 10 FLT / 10 GC, found {dict(counts)}")
    if set(meta["block"]) != set(AGES):
        raise CohortContractError(f"{arm}: expected blocks {AGES}, found {sorted(set(meta['block']))}")
    for block, sub in meta.groupby("block", sort=False):
        cell = sub["condition"].value_counts()
        if int(cell.get("FLT", 0)) != 5 or int(cell.get("GC", 0)) != 5:
            raise CohortContractError(f"{arm}/{block}: expected 5 FLT / 5 GC, found {dict(cell)}")


def check_persistence_cohorts(iss: MissionData, lar: MissionData, osd513: MissionData) -> None:
    """Raise unless ISS-T/LAR are 10/10 with YNG/OLD 5/5 blocks and OSD-513 is 9/9 in block 'all'."""
    _check_arm(iss, "ISS-T")
    _check_arm(lar, "LAR")
    meta = osd513.metadata
    counts = meta["condition"].value_counts()
    if int(counts.get("FLT", 0)) != 9 or int(counts.get("GC", 0)) != 9:
        raise CohortContractError(f"OSD-513: expected 9 FLT / 9 GC, found {dict(counts)}")
    if set(meta["block"]) != {"all"}:
        raise CohortContractError(f"OSD-513: expected the single block 'all', found {sorted(set(meta['block']))}")
    if not iss.expression.index.equals(lar.expression.index):
        raise CohortContractError("ISS-T and LAR must share one expression (VST) gene index")
    if not iss.counts.index.equals(lar.counts.index):
        raise CohortContractError("ISS-T and LAR must share one count-matrix gene index")
    overlap = set(iss.metadata.index) & set(lar.metadata.index)
    if overlap:
        raise CohortContractError(f"ISS-T and LAR animals overlap: {sorted(overlap)[:5]}")


def load_persistence_cohorts(config: Mapping[str, object], root: Path) -> dict[str, object]:
    """Load ISS-T (primary OSD-771), LAR, OSD-513 and the five primary missions.

    Loading is label-blind; no effect is computed here. The returned dict has keys
    ``iss``, ``lar``, ``osd513`` (MissionData) and ``primary`` (dict of MissionData).
    """
    primary = load_primary_missions(config, root)
    moderators = load_moderator_missions(config, root)
    iss = primary["OSD-771"]
    lar = moderators["OSD-771-LAR"]
    osd513 = moderators["OSD-513"]
    check_persistence_cohorts(iss, lar, osd513)
    return {"iss": iss, "lar": lar, "osd513": osd513, "primary": primary}


def joint_frame(
    iss: MissionData, lar: MissionData
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Concatenate the two OSD-771 arms into one 40-animal analysis set.

    ``MissionData.validate`` requires exactly FLT/GC, so the 2 x 2 is assembled here
    instead of in a loader. meta40 columns: condition, arm, age, block (= arm_age),
    source_sample.
    """
    if not iss.expression.index.equals(lar.expression.index):
        raise CohortContractError("ISS-T and LAR must share one expression gene index")
    if not iss.counts.index.equals(lar.counts.index):
        raise CohortContractError("ISS-T and LAR must share one count gene index")
    if set(iss.metadata.index) & set(lar.metadata.index):
        raise CohortContractError("ISS-T and LAR animal IDs must be disjoint")
    expr40 = pd.concat([iss.expression, lar.expression], axis=1)
    counts40 = pd.concat([iss.counts, lar.counts], axis=1)
    frames = []
    for arm, data in zip(ARMS, (iss, lar)):
        meta = pd.DataFrame(
            {
                "condition": data.metadata["condition"].astype(str),
                "arm": arm,
                "age": data.metadata["block"].astype(str),
                "source_sample": data.metadata["source_sample"].astype(str),
            },
            index=data.metadata.index,
        )
        meta["block"] = arm + "_" + meta["age"]
        frames.append(meta)
    meta40 = pd.concat(frames, axis=0)
    meta40.index.name = "animal"
    if list(expr40.columns) != list(meta40.index):
        raise CohortContractError("joint expression/metadata order mismatch")
    return expr40, counts40, meta40


def intersection_eligible(iss: MissionData, lar: MissionData, threshold: float) -> set[str]:
    """Genes CPM-eligible in BOTH arms, each under its own FLT/GC half-arm rule."""
    return cpm_eligible_genes(iss, threshold) & cpm_eligible_genes(lar, threshold)


# --------------------------------------------------------------------------- set definitions


def single_set_spec(genes: Iterable[str], *, role: str, minimum: int = 8) -> dict[str, object]:
    """Family entry in the ``load_family`` shape (one equal-weight subdomain)."""
    return {
        "role": role,
        "subdomains": {
            "atlas_markers": {
                "genes": {str(gene): 1 for gene in dict.fromkeys(genes)},
                "minimum_present": int(minimum),
            }
        },
    }


def set_genes(spec: Mapping[str, object]) -> list[str]:
    genes: list[str] = []
    for subdomain in spec["subdomains"].values():
        genes.extend(str(gene) for gene in subdomain["genes"])
    return list(dict.fromkeys(genes))


def structural_disjoint_spec(family: Mapping[str, Mapping[str, object]]) -> tuple[dict[str, object], list[str]]:
    """Structural scaffold control minus every podocyte high-specificity gene (K0 comparator)."""
    if PODOCYTE_SET not in family or STRUCTURAL_FULL_SET not in family:
        raise ValueError(f"family must contain {PODOCYTE_SET!r} and {STRUCTURAL_FULL_SET!r}")
    podocyte = set(set_genes(family[PODOCYTE_SET]))
    structural = set_genes(family[STRUCTURAL_FULL_SET])
    overlap = sorted(podocyte.intersection(structural))
    disjoint = [gene for gene in structural if gene not in podocyte]
    return single_set_spec(disjoint, role="key_secondary_structural_comparator"), overlap


def primary_evaluable_family(
    family: Mapping[str, Mapping[str, object]],
    missions: Mapping[str, MissionData],
    threshold: float,
    minimum: int = 8,
) -> tuple[dict[str, Mapping[str, object]], pd.DataFrame]:
    """The compartment-context evaluability rule: >= ``minimum`` CPM-eligible genes in every mission."""
    eligible = {name: cpm_eligible_genes(data, threshold) for name, data in missions.items()}
    retained: dict[str, Mapping[str, object]] = {}
    rows = []
    for name, spec in family.items():
        genes = set(set_genes(spec))
        counts = {mission: len(genes & elig) for mission, elig in eligible.items()}
        keep = min(counts.values()) >= minimum
        if keep:
            retained[name] = spec
        rows.append(
            {
                "set": name,
                "primary_evaluable": keep,
                **{f"n_eligible_{mission}": value for mission, value in counts.items()},
            }
        )
    return retained, pd.DataFrame(rows)


def _set_directions(
    spec: Mapping[str, object],
    eligible: set[str],
    expression_index: pd.Index,
    min_genes: int,
) -> tuple[dict[str, float], dict[str, list[str]], bool, int, int]:
    directions: dict[str, float] = {}
    subdomains: dict[str, list[str]] = {}
    ok = True
    n_defined = 0
    for name, subdomain in spec["subdomains"].items():
        requested = [str(gene) for gene in subdomain["genes"]]
        n_defined += len(requested)
        used = [gene for gene in requested if gene in eligible and gene in expression_index]
        required = int(subdomain.get("minimum_present", min_genes))
        ok &= len(used) >= required
        subdomains[str(name)] = used
        directions.update({gene: float(subdomain["genes"][gene]) for gene in used})
    return directions, subdomains, bool(ok), n_defined, len(directions)


def joint_scores(
    expr40: pd.DataFrame,
    family: Mapping[str, Mapping[str, object]],
    eligible: set[str],
    min_genes: int = 8,
    method: str = "mean",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score every evaluable set on the 40-animal joint set.

    Genes are z-scored over all 40 animals by ``score_signed_axis`` (frozen rule:
    z over the complete analysis set, never within blocks). A subdomain needs its
    ``minimum_present`` genes (``min_genes`` when absent); the call mirrors
    ``score_mission_axes`` (equal-weight subdomains, ``require_all_genes``).
    Returns (scores 40 x S, coverage).
    """
    scores: dict[str, pd.Series] = {}
    rows = []
    for name, spec in family.items():
        directions, subdomains, ok, n_defined, n_used = _set_directions(
            spec, eligible, expr40.index, min_genes
        )
        rows.append(
            {
                "set": name,
                "n_defined": n_defined,
                "n_used_joint": n_used,
                "evaluable_layer2": ok,
            }
        )
        if not ok:
            continue
        result = score_signed_axis(
            expr40,
            directions,
            method=method,
            subdomains=subdomains,
            min_genes_per_subdomain=1,
            require_all_genes=True,
        )
        if result.scores.isna().any():
            raise ValueError(f"joint score for {name} is not finite for every animal")
        scores[name] = result.scores
    return pd.DataFrame(scores, index=expr40.columns), pd.DataFrame(rows)


def arm_evaluable_family(
    data: MissionData,
    family: Mapping[str, Mapping[str, object]],
    threshold: float,
    min_genes: int = 8,
) -> tuple[dict[str, Mapping[str, object]], pd.DataFrame]:
    """Subset of ``family`` that meets every subdomain minimum in this cohort (own eligibility)."""
    eligible = cpm_eligible_genes(data, threshold)
    retained: dict[str, Mapping[str, object]] = {}
    rows = []
    for name, spec in family.items():
        _, _, ok, n_defined, n_used = _set_directions(spec, eligible, data.expression.index, min_genes)
        rows.append({"set": name, "n_defined": n_defined, "n_used": n_used, "evaluable": ok})
        if ok:
            retained[name] = spec
    return retained, pd.DataFrame(rows)


def score_arm(
    data: MissionData,
    family: Mapping[str, Mapping[str, object]],
    threshold: float,
    *,
    method: str = "mean",
    min_genes: int = 8,
    add_difference: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Frozen single-cohort scoring (``score_mission_axes``) on the cohort-evaluable subset.

    When both the podocyte set and ``structural_disjoint`` are scored, the column
    ``D_podocyte_minus_structural_disjoint`` (difference of the two scores) is added.
    """
    retained, coverage = arm_evaluable_family(data, family, threshold, min_genes)
    if not retained:
        raise ValueError(f"{data.mission}: no set is evaluable")
    scores, _, _ = score_mission_axes(data, retained, cpm_threshold=threshold, method=method)
    if add_difference and PODOCYTE_SET in scores and STRUCTURAL_DISJOINT in scores:
        scores[D_CONTRAST] = scores[PODOCYTE_SET] - scores[STRUCTURAL_DISJOINT]
    coverage.insert(0, "cohort", data.mission)
    return scores, coverage


# --------------------------------------------------------------------------- joint OLS / HC3


def design_matrix(meta40: pd.DataFrame) -> pd.DataFrame:
    """X (40 x 6): four arm x age block indicators (no intercept), F_T and F_L."""
    blocks = set(meta40["block"].astype(str))
    if blocks != set(BLOCK_COLUMNS):
        raise ValueError(f"meta40 blocks must be {BLOCK_COLUMNS}; found {sorted(blocks)}")
    flight = meta40["condition"].astype(str).eq("FLT")
    design = pd.DataFrame(index=meta40.index)
    for block in BLOCK_COLUMNS:
        design[block] = meta40["block"].astype(str).eq(block).astype(float)
    design["F_T"] = (flight & meta40["arm"].eq("ISS-T")).astype(float)
    design["F_L"] = (flight & meta40["arm"].eq("LAR")).astype(float)
    return design


@dataclass(frozen=True)
class OLSHC3Result:
    """OLS fit of several outcomes on one design with HC3 sandwich covariances."""

    beta: np.ndarray  # p x S
    cov_hc3: np.ndarray  # S x p x p
    leverage: np.ndarray  # n
    residuals: np.ndarray  # n x S
    xtx_inv: np.ndarray  # p x p
    df_resid: int


def _as_2d(values: np.ndarray | pd.DataFrame | pd.Series) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        array = array[:, None]
    if array.ndim != 2:
        raise ValueError("expected a vector or an n x S matrix")
    return array


def ols_hc3(X: np.ndarray | pd.DataFrame, Y: np.ndarray | pd.DataFrame) -> OLSHC3Result:
    """beta = (X'X)^-1 X'Y; HC3 = (X'X)^-1 X' diag(r^2/(1-h)^2) X (X'X)^-1 per outcome."""
    x = np.asarray(X, dtype=float)
    y = _as_2d(Y)
    n, p = x.shape
    if y.shape[0] != n:
        raise ValueError("X and Y must have the same number of rows")
    if np.linalg.matrix_rank(x) != p:
        raise ValueError("design matrix is rank deficient")
    if not np.all(np.isfinite(y)):
        raise ValueError("outcomes must be finite")
    xtx_inv = np.linalg.inv(x.T @ x)
    a = xtx_inv @ x.T
    beta = a @ y
    residuals = y - x @ beta
    leverage = np.einsum("ij,ji->i", x, a)
    if np.any(leverage >= 1.0 - 1e-10):
        raise ValueError("HC3 is undefined for observations with leverage 1")
    scale = residuals**2 / (1.0 - leverage)[:, None] ** 2  # n x S
    weighted = a[None, :, :] * scale.T[:, None, :]  # S x p x n
    cov = weighted @ a.T  # S x p x p
    return OLSHC3Result(
        beta=beta,
        cov_hc3=cov,
        leverage=leverage,
        residuals=residuals,
        xtx_inv=xtx_inv,
        df_resid=int(n - p),
    )


def arm_residual_df(X: np.ndarray, arm_mask: np.ndarray) -> int:
    """Residual df within one arm of a block-diagonal design (n_arm minus supported columns)."""
    sub = np.asarray(X, dtype=float)[np.asarray(arm_mask, dtype=bool)]
    used = np.any(np.abs(sub) > 0, axis=0)
    return int(sub.shape[0] - np.linalg.matrix_rank(sub[:, used]))


def satterthwaite_df(var_t: np.ndarray, df_t: float, var_l: np.ndarray, df_l: float, kappa: float) -> np.ndarray:
    """df for theta_L - kappa theta_T: (vL + k^2 vT)^2 / (vL^2/df_L + k^4 vT^2/df_T)."""
    var_t = np.asarray(var_t, dtype=float)
    var_l = np.asarray(var_l, dtype=float)
    numerator = (var_l + kappa**2 * var_t) ** 2
    denominator = var_l**2 / df_l + kappa**4 * var_t**2 / df_t
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(denominator > 0, numerator / denominator, np.nan)


def _t_quantile(prob: float, df: np.ndarray | float) -> np.ndarray:
    """Student-t quantile; df = inf gives the normal quantile, df <= 0 or NaN gives NaN."""
    df = np.asarray(df, dtype=float)
    with np.errstate(invalid="ignore"):
        valid = df > 0
    finite = valid & np.isfinite(df)
    out = np.asarray(stats.t.ppf(prob, np.where(finite, df, 1.0)), dtype=float)
    out = np.where(valid & ~np.isfinite(df), float(stats.norm.ppf(prob)), out)
    return np.where(valid, out, np.nan)


# --------------------------------------------------------------------------- Freedman--Lane


@dataclass(frozen=True)
class FreedmanLaneResult:
    """Freedman--Lane permutation of a studentized (HC3) coefficient, one per outcome."""

    estimate: np.ndarray
    standard_error: np.ndarray
    t_obs: np.ndarray
    t_null: np.ndarray | None
    p_upper: np.ndarray
    p_lower: np.ndarray
    p_two_sided: np.ndarray
    p_max_t: np.ndarray | None
    family: tuple[int, ...] | None
    n_permutations: int
    seed: object


def _orthonormal_basis(matrix: np.ndarray, tol: float = 1e-10) -> np.ndarray:
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    if s.size == 0:
        return np.zeros((matrix.shape[0], 0))
    rank = int(np.sum(s > tol * max(1.0, float(s[0]))))
    return u[:, :rank]


def block_positions(block_ids: Sequence[object]) -> list[np.ndarray]:
    labels = pd.Series(np.asarray(block_ids).astype(str))
    return [np.flatnonzero(labels.to_numpy() == block) for block in pd.unique(labels)]


def within_block_permutations(
    positions: Sequence[np.ndarray], n: int, batch: int, rng: np.random.Generator
) -> np.ndarray:
    """(batch x n) index arrays; row b sends position i to a position in i's block."""
    out = np.empty((batch, n), dtype=np.intp)
    covered = np.zeros(n, dtype=bool)
    for pos in positions:
        keys = rng.random((batch, len(pos)))
        out[:, pos] = pos[np.argsort(keys, axis=1)]
        covered[pos] = True
    if not covered.all():
        raise ValueError("block positions do not cover every observation")
    return out


def _rng(seed: object) -> np.random.Generator:
    if isinstance(seed, np.random.SeedSequence):
        return np.random.default_rng(seed)
    return np.random.default_rng(np.random.SeedSequence(seed))


def freedman_lane(
    Y: np.ndarray | pd.DataFrame,
    Z: np.ndarray,
    x: np.ndarray,
    block_ids: Sequence[object],
    n_perm: int,
    seed: object,
    chunk: int = 2048,
    *,
    family: Sequence[int] | None = None,
    keep_null: bool = True,
    memory_budget: int = 4_000_000,
) -> FreedmanLaneResult:
    """Freedman--Lane test of the coefficient on ``x`` in ``Y ~ [Z, x]`` with HC3 studentization.

    e = M_Z Y; a = M_Z x / (x' M_Z x); for each within-block permutation pi (the same pi
    for every outcome column): delta* = a' e[pi], r* = M_full e[pi],
    var* = sum_i a_i^2 r*_i^2 / (1 - h_i)^2, t* = delta*/sqrt(var*). With pi = identity
    this is exactly the HC3 t of the full model. The permutation stream depends only on
    (seed, block structure, chunk), so calls sharing those share their permutations.
    """
    if n_perm < 1 or chunk < 1:
        raise ValueError("n_perm and chunk must be positive")
    y = _as_2d(Y)
    n, n_out = y.shape
    z = np.asarray(Z, dtype=float).reshape(n, -1)
    xv = np.asarray(x, dtype=float).reshape(n)
    full = np.column_stack([z, xv])
    basis_z = _orthonormal_basis(z)
    basis_full = _orthonormal_basis(full)
    if basis_full.shape[1] != basis_z.shape[1] + 1:
        raise ValueError("tested column lies in the span of the nuisance design")
    e = y - basis_z @ (basis_z.T @ y)
    x_res = xv - basis_z @ (basis_z.T @ xv)
    a = x_res / float(x_res @ x_res)
    leverage = np.sum(basis_full**2, axis=1)
    if np.any(leverage >= 1.0 - 1e-10):
        raise ValueError("HC3 is undefined for observations with leverage 1")
    w = a**2 / (1.0 - leverage) ** 2

    def statistic(err: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        # err: (c, n, S); returns delta*, HC3 var*, t* each (c, S)
        delta = np.tensordot(err, a, axes=([1], [0]))
        resid = err - basis_full @ (basis_full.T @ err)
        var = np.tensordot(resid * resid, w, axes=([1], [0]))
        with np.errstate(divide="ignore", invalid="ignore"):
            t = delta / np.sqrt(var)
        return delta, var, t

    delta_2d, var_2d, t_2d = statistic(e[None, :, :])
    delta_obs, se_obs, t_obs = delta_2d[0], np.sqrt(var_2d[0]), t_2d[0]
    if not np.all(np.isfinite(t_obs)):
        raise ValueError("observed Freedman-Lane statistic is not finite (zero residual variance?)")
    tol = _TIE_RTOL * np.maximum(1.0, np.abs(t_obs))

    fam = None if family is None else tuple(int(i) for i in family)
    positions = block_positions(block_ids)
    rng = _rng(seed)
    upper = np.zeros(n_out, dtype=np.int64)
    lower = np.zeros(n_out, dtype=np.int64)
    two = np.zeros(n_out, dtype=np.int64)
    max_abs = np.empty(n_perm, dtype=float) if fam else None
    null = np.empty((n_perm, n_out), dtype=float) if keep_null else None
    sub = max(1, min(chunk, int(memory_budget // max(1, n * n_out))))
    for start in range(0, n_perm, chunk):
        stop = min(start + chunk, n_perm)
        perms = within_block_permutations(positions, n, stop - start, rng)
        for s0 in range(0, stop - start, sub):
            s1 = min(s0 + sub, stop - start)
            _, _, t = statistic(e[perms[s0:s1]])
            t = np.where(np.isfinite(t), t, 0.0)
            upper += np.count_nonzero(t >= t_obs - tol, axis=0)
            lower += np.count_nonzero(t <= t_obs + tol, axis=0)
            two += np.count_nonzero(np.abs(t) >= np.abs(t_obs) - tol, axis=0)
            if fam:
                max_abs[start + s0 : start + s1] = np.max(np.abs(t[:, fam]), axis=1)
            if null is not None:
                null[start + s0 : start + s1] = t
    denom = n_perm + 1.0
    p_max = None
    if fam:
        p_max = np.full(n_out, np.nan)
        for index in fam:
            p_max[index] = (
                1.0 + np.count_nonzero(max_abs >= abs(t_obs[index]) - tol[index])
            ) / denom
    return FreedmanLaneResult(
        estimate=delta_obs,
        standard_error=se_obs,
        t_obs=t_obs,
        t_null=null,
        p_upper=(1.0 + upper) / denom,
        p_lower=(1.0 + lower) / denom,
        p_two_sided=(1.0 + two) / denom,
        p_max_t=p_max,
        family=fam,
        n_permutations=int(n_perm),
        seed=seed,
    )


def freedman_lane_contrast(
    Y: np.ndarray | pd.DataFrame,
    X_blocks: np.ndarray,
    f_t: np.ndarray,
    f_l: np.ndarray,
    kappa: float,
    block_ids: Sequence[object],
    n_perm: int,
    seed: object,
    chunk: int = 2048,
    *,
    family: Sequence[int] | None = None,
    keep_null: bool = True,
) -> FreedmanLaneResult:
    """Test delta = theta_L - kappa theta_T by Freedman--Lane.

    Reduced design Z = [X_blocks, F_T + kappa F_L]; tested column F_L. Because
    [Z, F_L] spans the same space as [X_blocks, F_T, F_L], the F_L coefficient of the
    reparameterized fit is delta = theta_L - kappa theta_T. ``X_blocks`` may carry extra
    nuisance columns (e.g. arm-specific structural slopes for Layer 3).
    """
    f_t = np.asarray(f_t, dtype=float)
    f_l = np.asarray(f_l, dtype=float)
    z = np.column_stack([np.asarray(X_blocks, dtype=float), f_t + kappa * f_l])
    return freedman_lane(
        Y, z, f_l, block_ids, n_perm, seed, chunk, family=family, keep_null=keep_null
    )


# --------------------------------------------------------------------------- exact single-arm tests


@dataclass(frozen=True)
class AllocationSet:
    """Within-block flight allocations (rows) over the cohort's animals (columns)."""

    masks: np.ndarray  # N x n boolean
    exact: bool
    n_total: int  # number of distinct allocations in the reference set
    observed_index: int | None


def enumerate_block_allocations(
    block_labels: Sequence[object],
    n_flight_per_block: Mapping[object, int],
    *,
    observed: Sequence[bool] | None = None,
    max_exact: int = 2_000_000,
    n_monte_carlo: int = 100_000,
    seed: object = 0,
) -> AllocationSet:
    """Cartesian product of per-block ``itertools.combinations`` as an N x n boolean matrix.

    Exact when the product has at most ``max_exact`` members. Otherwise the observed
    allocation is row 0 and ``n_monte_carlo`` uniform within-block allocations follow,
    so ``#(stat* >= stat_obs) / N`` is the plus-one Monte-Carlo p-value.
    """
    labels = np.asarray(block_labels).astype(str)
    n = len(labels)
    blocks = list(pd.unique(pd.Series(labels)))
    positions = [np.flatnonzero(labels == block) for block in blocks]
    lookup = {str(key): int(value) for key, value in n_flight_per_block.items()}
    counts = []
    for block, pos in zip(blocks, positions):
        if block not in lookup:
            raise ValueError(f"no flight count for block {block!r}")
        k = lookup[block]
        if not 0 <= k <= len(pos):
            raise ValueError(f"block {block!r}: cannot place {k} flights among {len(pos)} animals")
        counts.append(k)
    sizes = [math.comb(len(pos), k) for pos, k in zip(positions, counts)]
    total = int(np.prod(sizes, dtype=object))
    obs = None if observed is None else np.asarray(observed, dtype=bool)
    if obs is not None:
        if obs.shape != (n,):
            raise ValueError("observed allocation has the wrong length")
        for pos, k in zip(positions, counts):
            if int(obs[pos].sum()) != k:
                raise ValueError("observed allocation does not match the block flight counts")
    if total <= max_exact:
        per_block = []
        for pos, k in zip(positions, counts):
            combos = np.array(list(itertools.combinations(range(len(pos)), k)), dtype=np.intp)
            combos = combos.reshape(len(combos), k)
            mask = np.zeros((len(combos), len(pos)), dtype=bool)
            if k:
                mask[np.arange(len(combos))[:, None], combos] = True
            per_block.append(mask)
        grid = np.indices(sizes).reshape(len(sizes), -1)
        masks = np.zeros((total, n), dtype=bool)
        for index, (pos, mask) in enumerate(zip(positions, per_block)):
            masks[:, pos] = mask[grid[index]]
        observed_index = None
        if obs is not None:
            hits = np.flatnonzero(np.all(masks == obs[None, :], axis=1))
            if len(hits) != 1:
                raise RuntimeError("observed allocation not found exactly once in the enumeration")
            observed_index = int(hits[0])
        return AllocationSet(masks=masks, exact=True, n_total=total, observed_index=observed_index)
    if obs is None:
        raise ValueError("Monte-Carlo allocation sets need the observed allocation")
    rng = _rng(seed)
    masks = np.zeros((n_monte_carlo + 1, n), dtype=bool)
    masks[0] = obs
    for pos, k in zip(positions, counts):
        keys = rng.random((n_monte_carlo, len(pos)))
        chosen = np.argsort(keys, axis=1)[:, :k]
        block_mask = np.zeros((n_monte_carlo, len(pos)), dtype=bool)
        block_mask[np.arange(n_monte_carlo)[:, None], chosen] = True
        masks[1:, pos] = block_mask
    return AllocationSet(masks=masks, exact=False, n_total=total, observed_index=0)


def _fixed_effect_batch(
    block_values: Sequence[np.ndarray], block_masks: Sequence[np.ndarray]
) -> tuple[np.ndarray, np.ndarray]:
    """Per-stratum Hedges g (``_block_hedges_g_batch``) combined by inverse-variance weights."""
    sum_w = None
    sum_wg = None
    for values, mask in zip(block_values, block_masks):
        g, v = _block_hedges_g_batch(values, mask)
        w = 1.0 / v
        sum_w = w if sum_w is None else sum_w + w
        sum_wg = w * g if sum_wg is None else sum_wg + w * g
    return sum_wg / sum_w, 1.0 / sum_w


def _arm_blocks(scores: pd.DataFrame, meta: pd.DataFrame) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray], np.ndarray]:
    if not set(scores.index) == set(meta.index):
        raise ValueError("scores and metadata must cover the same animals")
    values = scores.loc[meta.index].to_numpy(dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("scores must be finite")
    flight = meta["condition"].astype(str).eq("FLT").to_numpy()
    positions = block_positions(meta["block"])
    block_values = [values[pos] for pos in positions]
    return values, positions, block_values, flight


def _family_columns(
    family_sets: Mapping[str, Sequence[str]] | Sequence[str] | None, columns: Sequence[str]
) -> dict[str, list[int]]:
    if family_sets is None:
        return {}
    if not isinstance(family_sets, Mapping):
        family_sets = {"family": list(family_sets)}
    index = {name: i for i, name in enumerate(columns)}
    return {
        str(name): [index[s] for s in members if s in index]
        for name, members in family_sets.items()
        if any(s in index for s in members)
    }


def observed_arm_effects(scores_arm: pd.DataFrame, meta_arm: pd.DataFrame) -> pd.DataFrame:
    """Observed fixed-effect Hedges g per set (no enumeration)."""
    _, positions, block_values, flight = _arm_blocks(scores_arm, meta_arm)
    g, v = _fixed_effect_batch(block_values, [flight[pos][None, :] for pos in positions])
    return pd.DataFrame(
        {"estimate_g": g[0], "variance": v[0]}, index=pd.Index(scores_arm.columns, name="set")
    )


def arm_exact_test(
    scores_arm: pd.DataFrame,
    meta_arm: pd.DataFrame,
    family_sets: Mapping[str, Sequence[str]] | Sequence[str] | None = None,
    *,
    allocations: AllocationSet | None = None,
    chunk: int = 16_384,
) -> pd.DataFrame:
    """Exact within-block label enumeration for a single cohort.

    For each allocation: per-stratum Hedges g (``_block_hedges_g_batch``), inverse-variance
    fixed effect, z = g/sqrt(v). p_exact = #(|z*| >= |z_obs|)/N with the observed allocation
    included; max-T family-wise p = #(max_{s in F}|z*_s| >= |z_obs,s|)/N per named family.
    """
    _, positions, block_values, flight = _arm_blocks(scores_arm, meta_arm)
    columns = [str(c) for c in scores_arm.columns]
    if allocations is None:
        n_flight = {
            str(block): int(sub["condition"].eq("FLT").sum())
            for block, sub in meta_arm.groupby("block", sort=False)
        }
        allocations = enumerate_block_allocations(meta_arm["block"], n_flight, observed=flight)
    g_obs, v_obs = _fixed_effect_batch(block_values, [flight[pos][None, :] for pos in positions])
    g_obs, v_obs = g_obs[0], v_obs[0]
    z_obs = g_obs / np.sqrt(v_obs)
    abs_obs = np.abs(z_obs)
    tol = _TIE_RTOL * np.maximum(1.0, abs_obs)
    families = _family_columns(family_sets, columns)
    exceed = np.zeros(len(columns), dtype=np.int64)
    fam_exceed = {name: np.zeros(len(cols), dtype=np.int64) for name, cols in families.items()}
    masks = allocations.masks
    for start in range(0, len(masks), chunk):
        part = masks[start : start + chunk]
        g, v = _fixed_effect_batch(block_values, [part[:, pos] for pos in positions])
        z = np.abs(g / np.sqrt(v))
        exceed += np.count_nonzero(z >= abs_obs - tol, axis=0)
        for name, cols in families.items():
            fam_max = np.max(z[:, cols], axis=1)
            fam_exceed[name] += np.count_nonzero(
                fam_max[:, None] >= abs_obs[cols][None, :] - tol[cols][None, :], axis=0
            )
    n_alloc = len(masks)
    se = np.sqrt(v_obs)
    table = pd.DataFrame(
        {
            "set": columns,
            "estimate_g": g_obs,
            "variance": v_obs,
            "standard_error": se,
            "ci_low": g_obs - stats.norm.ppf(0.975) * se,
            "ci_high": g_obs + stats.norm.ppf(0.975) * se,
            "z": z_obs,
            "p_exact": exceed / n_alloc,
            "max_t_fwer_exact": np.nan,
            "max_t_family": "",
            "n_allocations": n_alloc,
            "allocation_exact": allocations.exact,
            "n_flight": int(flight.sum()),
            "n_control": int((~flight).sum()),
            "n_strata": len(positions),
        }
    )
    for name, cols in families.items():
        column = f"max_t_fwer_exact__{name}"
        table[column] = np.nan
        table.loc[cols, column] = fam_exceed[name] / n_alloc
        unset = table.index.isin(cols) & table["max_t_family"].eq("")
        table.loc[unset, "max_t_fwer_exact"] = table.loc[unset, column]
        table.loc[unset, "max_t_family"] = name
    return table


def pattern_concordance_exact(
    ref_effects: pd.Series | Mapping[str, float],
    scores_arm: pd.DataFrame,
    meta_arm: pd.DataFrame,
    *,
    allocations: AllocationSet | None = None,
    chunk: int = 16_384,
) -> dict[str, object]:
    """Pearson r across sets between reference effects and the arm's per-set g.

    Null: every within-block allocation of the arm's labels (reference held fixed).
    One-sided p = #(r* >= r_obs)/N with the observed allocation included.
    """
    ref = pd.Series(ref_effects, dtype=float)
    common = [c for c in scores_arm.columns if c in ref.index and np.isfinite(ref[c])]
    if len(common) < 3:
        raise ValueError("pattern concordance needs at least three shared sets")
    scores = scores_arm[common]
    _, positions, block_values, flight = _arm_blocks(scores, meta_arm)
    if allocations is None:
        n_flight = {
            str(block): int(sub["condition"].eq("FLT").sum())
            for block, sub in meta_arm.groupby("block", sort=False)
        }
        allocations = enumerate_block_allocations(meta_arm["block"], n_flight, observed=flight)
    ref_c = ref.loc[common].to_numpy(dtype=float)
    ref_c = ref_c - ref_c.mean()
    ref_norm = float(np.sqrt(ref_c @ ref_c))
    if ref_norm <= 0:
        raise ValueError("reference effects are constant")

    def pearson(g: np.ndarray) -> np.ndarray:
        centered = g - g.mean(axis=1, keepdims=True)
        norm = np.sqrt(np.sum(centered**2, axis=1))
        with np.errstate(divide="ignore", invalid="ignore"):
            return (centered @ ref_c) / (norm * ref_norm)

    g_obs, _ = _fixed_effect_batch(block_values, [flight[pos][None, :] for pos in positions])
    r_obs = float(pearson(g_obs)[0])
    tol = _TIE_RTOL * max(1.0, abs(r_obs))
    exceed = 0
    for start in range(0, len(allocations.masks), chunk):
        part = allocations.masks[start : start + chunk]
        g, _ = _fixed_effect_batch(block_values, [part[:, pos] for pos in positions])
        r = pearson(g)
        exceed += int(np.count_nonzero(np.nan_to_num(r, nan=-np.inf) >= r_obs - tol))
    spearman = float(stats.spearmanr(g_obs[0], ref.loc[common].to_numpy(dtype=float)).statistic)
    return {
        "n_sets": len(common),
        "pearson_r": r_obs,
        "p_one_sided_exact": exceed / len(allocations.masks),
        "n_allocations": len(allocations.masks),
        "allocation_exact": allocations.exact,
        "spearman_r_descriptive": spearman,
        "sets": common,
    }


# --------------------------------------------------------------------------- ratio intervals


def fieller_ratio(
    tT: float, vT: float, tL: float, vL: float, df: float, level: float, cov: float = 0.0
) -> tuple[float, float, bool]:
    """Fieller set {rho: (tL - rho tT)^2 <= q^2 (vL + rho^2 vT - 2 rho cov)}, q = t_{df,1-(1-level)/2}.

    a = tT^2 - q^2 vT, b = -2 (tT tL - q^2 cov), c = tL^2 - q^2 vL. If a <= 0 the set is
    not a bounded interval (theta_T not significantly non-zero) and (nan, nan, False)
    is returned.

    Inversion: rho0 lies outside the set exactly when the t test of
    theta_L - rho0 theta_T = 0 rejects at two-sided 1 - level. So the 90% lower bound is
    the largest margin compatible with persistence at one-sided 5% (reported as
    ``rho_lower90_fieller``): every margin k below it is rejected in the persist direction
    by the HC3 test with the same df.
    """
    if not 0 < level < 1:
        raise ValueError("level must lie in (0, 1)")
    q = float(_t_quantile(1.0 - (1.0 - level) / 2.0, df))
    a = tT**2 - q**2 * vT
    b = -2.0 * (tT * tL - q**2 * cov)
    c = tL**2 - q**2 * vL
    if not np.isfinite(a) or a <= 0:
        return float("nan"), float("nan"), False
    disc = b**2 - 4.0 * a * c
    if disc < 0:  # pragma: no cover - impossible with cov = 0 and a > 0
        return float("nan"), float("nan"), False
    root = math.sqrt(disc)
    low, high = sorted(((-b - root) / (2.0 * a), (-b + root) / (2.0 * a)))
    return float(low), float(high), True


def _cell_positions(meta40: pd.DataFrame) -> list[np.ndarray]:
    keys = meta40["arm"].astype(str) + "|" + meta40["age"].astype(str) + "|" + meta40["condition"].astype(str)
    return block_positions(keys)


def stratified_bootstrap_ratio(
    scores: pd.DataFrame,
    meta40: pd.DataFrame,
    B: int = 10_000,
    seed: object = 0,
    *,
    chunk: int = 1_000,
    kappa_reference: float = 0.0,
) -> pd.DataFrame:
    """Resample animals within the 8 arm x age x condition cells, refit OLS, rho* = theta_L*/theta_T*.

    Scores should be terminal-oriented. Returns the percentile 95% interval for rho,
    the fraction of resamples with theta_T* <= 0 and ``bootstrap_unreliable`` when that
    fraction exceeds 0.025 (rho is then dominated by near-zero denominators).
    """
    design = design_matrix(meta40)
    x = design.to_numpy(dtype=float)
    a = np.linalg.inv(x.T @ x) @ x.T
    a_t = a[list(design.columns).index("F_T")]
    a_l = a[list(design.columns).index("F_L")]
    values = scores.loc[meta40.index].to_numpy(dtype=float)
    cells = _cell_positions(meta40)
    rng = _rng(seed)
    n = len(meta40)
    theta_t = np.empty((B, values.shape[1]))
    theta_l = np.empty((B, values.shape[1]))
    for start in range(0, B, chunk):
        stop = min(start + chunk, B)
        idx = np.empty((stop - start, n), dtype=np.intp)
        for pos in cells:
            idx[:, pos] = pos[rng.integers(0, len(pos), size=(stop - start, len(pos)))]
        sample = values[idx]  # c x n x S
        theta_t[start:stop] = np.tensordot(sample, a_t, axes=([1], [0]))
        theta_l[start:stop] = np.tensordot(sample, a_l, axes=([1], [0]))
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = theta_l / theta_t
    frac = np.mean(theta_t <= 0.0, axis=0)
    rho = np.where(np.isnan(rho), np.inf, rho)
    low, high = np.quantile(rho, [0.025, 0.975], axis=0)
    return pd.DataFrame(
        {
            "rho_boot95_low": low,
            "rho_boot95_high": high,
            "boot_frac_thetaT_le0": frac,
            "bootstrap_unreliable": frac > 0.025,
            "n_bootstrap": B,
        },
        index=pd.Index(scores.columns, name="set"),
    )


def posterior_summaries(
    tT: np.ndarray | float,
    vT: np.ndarray | float,
    tL: np.ndarray | float,
    vL: np.ndarray | float,
    draws: int = 1_000_000,
    seed: object = 0,
    *,
    kappa: float = 0.5,
    chunk: int = 50_000,
) -> pd.DataFrame:
    """Flat priors with normal likelihoods (independent arms); descriptive only.

    Returns P(rho >= kappa and theta_T > 0), P(theta_L <= 0) and P(rho >= 1 and theta_T > 0)
    from common random numbers shared across sets.
    """
    tT, vT, tL, vL = (np.atleast_1d(np.asarray(v, dtype=float)) for v in (tT, vT, tL, vL))
    rng = _rng(seed)
    persist = np.zeros(len(tT))
    lar_le0 = np.zeros(len(tT))
    full = np.zeros(len(tT))
    sd_t, sd_l = np.sqrt(vT), np.sqrt(vL)
    for start in range(0, draws, chunk):
        size = min(chunk, draws - start)
        z1 = rng.standard_normal(size)[:, None]
        z2 = rng.standard_normal(size)[:, None]
        th_t = tT[None, :] + sd_t[None, :] * z1
        th_l = tL[None, :] + sd_l[None, :] * z2
        positive = th_t > 0
        persist += np.count_nonzero(positive & (th_l >= kappa * th_t), axis=0)
        lar_le0 += np.count_nonzero(th_l <= 0, axis=0)
        full += np.count_nonzero(positive & (th_l >= th_t), axis=0)
    return pd.DataFrame(
        {
            "post_p_rho_ge_0_5": persist / draws,
            "post_p_thetaL_le_0": lar_le0 / draws,
            "post_p_rho_ge_1": full / draws,
            "posterior_draws": draws,
        }
    )


# --------------------------------------------------------------------------- joint contrasts


def joint_hc3_contrasts(
    scores: pd.DataFrame,
    meta40: pd.DataFrame,
    *,
    contrasts: Mapping[str, float],
    n_perm: int,
    seed: object,
    chunk: int = 2048,
    extra_nuisance: np.ndarray | None = None,
    families: Mapping[str, Sequence[str]] | None = None,
    keep_null: bool = True,
) -> tuple[pd.DataFrame, dict[str, FreedmanLaneResult]]:
    """Layer-2/3 core: OLS + HC3 per set, contrasts theta_L - k theta_T with Satterthwaite df,
    and one Freedman--Lane run per contrast (same permutations for every set and contrast).

    ``scores`` must already be terminal-oriented. ``families`` maps a family name to the
    member sets for max-|T| family-wise p-values of each contrast.
    """
    design = design_matrix(meta40)
    blocks = design[list(BLOCK_COLUMNS)].to_numpy(dtype=float)
    f_t = design["F_T"].to_numpy(dtype=float)
    f_l = design["F_L"].to_numpy(dtype=float)
    x = np.column_stack([blocks, f_t, f_l])
    if extra_nuisance is not None:
        x = np.column_stack([x, np.asarray(extra_nuisance, dtype=float)])
    y = scores.loc[meta40.index].to_numpy(dtype=float)
    fit = ols_hc3(x, y)
    it, il = 4, 5
    theta_t, theta_l = fit.beta[it], fit.beta[il]
    var_t, var_l = fit.cov_hc3[:, it, it], fit.cov_hc3[:, il, il]
    cov_tl = fit.cov_hc3[:, it, il]
    if np.any(np.abs(cov_tl) > 1e-12):
        raise AssertionError("cov(theta_T, theta_L) must vanish for an arm-block-diagonal design")
    iss_mask = meta40["arm"].eq("ISS-T").to_numpy()
    df_t = arm_residual_df(x, iss_mask)
    df_l = arm_residual_df(x, ~iss_mask)
    rss_t = np.sum(fit.residuals[iss_mask] ** 2, axis=0)
    rss_l = np.sum(fit.residuals[~iss_mask] ** 2, axis=0)
    sd_t, sd_l = np.sqrt(rss_t / df_t), np.sqrt(rss_l / df_l)
    q95_t = _t_quantile(0.975, df_t)
    q95_l = _t_quantile(0.975, df_l)
    se_t, se_l = np.sqrt(var_t), np.sqrt(var_l)
    table = pd.DataFrame(index=pd.Index(scores.columns, name="set"))
    table["theta_isst"] = theta_t
    table["se_isst"] = se_t
    table["df_isst"] = df_t
    table["theta_isst_ci95_low"] = theta_t - q95_t * se_t
    table["theta_isst_ci95_high"] = theta_t + q95_t * se_t
    table["p_theta_isst_hc3"] = 2.0 * stats.t.sf(np.abs(theta_t / se_t), df_t)
    table["theta_lar"] = theta_l
    table["se_lar"] = se_l
    table["df_lar"] = df_l
    table["theta_lar_ci95_low"] = theta_l - q95_l * se_l
    table["theta_lar_ci95_high"] = theta_l + q95_l * se_l
    table["cov_isst_lar_hc3"] = cov_tl
    table["resid_sd_isst"] = sd_t
    table["resid_sd_lar"] = sd_l
    n_t_f = int(np.sum(iss_mask & (f_t > 0)))
    n_t_c = int(np.sum(iss_mask & (f_t == 0)))
    n_l_f = int(np.sum(~iss_mask & (f_l > 0)))
    n_l_c = int(np.sum(~iss_mask & (f_l == 0)))
    g_t, g_l = theta_t / sd_t, theta_l / sd_l
    table["g_isst_joint"] = g_t
    table["g_lar_joint"] = g_l
    table["g_isst_joint_var"] = (n_t_f + n_t_c) / (n_t_f * n_t_c) + g_t**2 / (2.0 * df_t)
    table["g_lar_joint_var"] = (n_l_f + n_l_c) / (n_l_f * n_l_c) + g_l**2 / (2.0 * df_l)
    nuisance = blocks if extra_nuisance is None else np.column_stack([blocks, extra_nuisance])
    block_ids = meta40["block"].to_numpy()
    fl_results: dict[str, FreedmanLaneResult] = {}
    columns = list(scores.columns)
    for name, kappa in contrasts.items():
        est = theta_l - kappa * theta_t
        var = var_l + kappa**2 * var_t
        se = np.sqrt(var)
        df = satterthwaite_df(var_t, df_t, var_l, df_l, kappa)
        t = est / se
        table[name] = est
        table[f"{name}_se"] = se
        table[f"{name}_df"] = df
        table[f"{name}_t"] = t
        table[f"p_{name}_hc3"] = 2.0 * stats.t.sf(np.abs(t), df)
        table[f"p_{name}_hc3_upper"] = stats.t.sf(t, df)
        table[f"p_{name}_hc3_lower"] = stats.t.cdf(t, df)
        for level in (95, 90):
            q = _t_quantile(1.0 - (1.0 - level / 100.0) / 2.0, df)
            table[f"{name}_ci{level}_low"] = est - q * se
            table[f"{name}_ci{level}_high"] = est + q * se
        fl = freedman_lane_contrast(
            y,
            nuisance,
            f_t,
            f_l,
            kappa,
            block_ids,
            n_perm,
            seed,
            chunk,
            keep_null=bool(keep_null or families),
        )
        fl_results[name] = fl
        table[f"{name}_t_fl_identity"] = fl.t_obs
        table[f"p_{name}_fl"] = fl.p_two_sided
        table[f"p_{name}_fl_upper"] = fl.p_upper
        table[f"p_{name}_fl_lower"] = fl.p_lower
        for family_name, members in (families or {}).items():
            column = f"{name}_maxT_fwer_fl__{family_name}"
            cols = [columns.index(s) for s in members if s in columns]
            values = np.full(len(columns), np.nan)
            if cols:
                values[cols] = max_t_pvalues(fl.t_obs, fl.t_null, cols)
            table[column] = values
    return table, fl_results


def parse_margin_sweep(prereg: Mapping[str, object] | None, primary_kappa: float) -> tuple[float, ...]:
    """Kappas of the preregistered margin sweep (``margin_sweep: {kappas: [...]}``).

    Without a ``margin_sweep`` block the sweep is the primary margin alone. The list must
    contain the primary kappa (only that kappa sets the classification).
    """
    block = (prereg or {}).get("margin_sweep")
    if block is None:
        return (float(primary_kappa),)
    values = block.get("kappas") if isinstance(block, Mapping) else block
    if not isinstance(values, Sequence) or isinstance(values, str) or not values:
        raise RuleParseError("margin_sweep.kappas must be a non-empty list of numbers")
    kappas = tuple(float(value) for value in values)
    if len(set(kappas)) != len(kappas) or not all(np.isfinite(kappas)):
        raise RuleParseError(f"margin_sweep.kappas must be distinct finite numbers: {values}")
    if not any(abs(k - primary_kappa) < 1e-12 for k in kappas):
        raise RuleParseError(f"margin_sweep.kappas {kappas} must include the primary margin kappa {primary_kappa}")
    if isinstance(block, Mapping) and block.get("primary_kappa") is not None:
        declared = float(block["primary_kappa"])
        if abs(declared - primary_kappa) > 1e-12:
            raise RuleParseError(
                f"margin_sweep.primary_kappa {declared} disagrees with contrasts.margin kappa {primary_kappa}"
            )
    return tuple(sorted(kappas))


def margin_sweep_targets(prereg: Mapping[str, object] | None) -> tuple[str, ...]:
    """Sets the sweep applies to (``margin_sweep.applies_to``; default T1-T3)."""
    block = (prereg or {}).get("margin_sweep")
    targets = block.get("applies_to") if isinstance(block, Mapping) else None
    if targets is None:
        return KEY_TARGETS
    targets = tuple(str(t) for t in targets)
    unknown = [t for t in targets if t not in KEY_TARGETS]
    if unknown:
        raise RuleParseError(f"margin_sweep.applies_to names unknown targets {unknown}; expected a subset of {KEY_TARGETS}")
    return targets


def margin_sweep(
    scores: pd.DataFrame,
    meta40: pd.DataFrame,
    kappas: Sequence[float],
    *,
    primary_kappa: float,
    n_perm: int,
    seed: object,
    chunk: int = 2048,
) -> pd.DataFrame:
    """Margin M_k = theta_L - k theta_T for every k (terminal-oriented scores), long format.

    Per k: HC3 SE sqrt(vL + k^2 vT), Satterthwaite df, 90% CI and Freedman--Lane one-sided
    p-values for M_k > 0 (persist) and M_k < 0 (decay). Every k reuses the same permutation
    stream (same seed, blocks and chunk), so k = primary reproduces the primary margin row
    exactly. k = 0 tests "any retention" (theta_L > 0); k = 1 is "no decay" and equals the
    interaction test. The Fieller 90% lower bound of rho inverts the HC3 half of this test
    (with the per-arm df instead of Satterthwaite): every k below it is rejected in the
    persist direction at one-sided 5%.
    """
    contrasts = {f"k{index}": float(k) for index, k in enumerate(kappas)}
    table, _ = joint_hc3_contrasts(
        scores, meta40, contrasts=contrasts, n_perm=n_perm, seed=seed, chunk=chunk, keep_null=False
    )
    rows = []
    for name, kappa in contrasts.items():
        part = pd.DataFrame(
            {
                "set": table.index,
                "kappa": kappa,
                "is_primary": abs(kappa - primary_kappa) < 1e-12,
                "margin": table[name].to_numpy(),
                "se": table[f"{name}_se"].to_numpy(),
                "df": table[f"{name}_df"].to_numpy(),
                "ci90_low": table[f"{name}_ci90_low"].to_numpy(),
                "ci90_high": table[f"{name}_ci90_high"].to_numpy(),
                "t": table[f"{name}_t"].to_numpy(),
                "p_persist_fl": table[f"p_{name}_fl_upper"].to_numpy(),
                "p_decay_fl": table[f"p_{name}_fl_lower"].to_numpy(),
                "theta_isst": table["theta_isst"].to_numpy(),
                "se_isst": table["se_isst"].to_numpy(),
                "df_isst": table["df_isst"].to_numpy(),
                "theta_lar": table["theta_lar"].to_numpy(),
                "se_lar": table["se_lar"].to_numpy(),
                "df_lar": table["df_lar"].to_numpy(),
            }
        )
        rows.append(part)
    return pd.concat(rows, ignore_index=True)


def max_t_pvalues(t_obs: np.ndarray, t_null: np.ndarray, columns: Sequence[int]) -> np.ndarray:
    """Plus-one max-|T| family-wise p over ``columns`` of a shared-permutation null."""
    cols = list(columns)
    observed = np.abs(np.asarray(t_obs, dtype=float)[cols])
    null = np.where(np.isfinite(t_null[:, cols]), np.abs(t_null[:, cols]), 0.0)
    maximum = np.max(null, axis=1)
    tol = _TIE_RTOL * np.maximum(1.0, observed)
    return np.array(
        [(1.0 + np.count_nonzero(maximum >= obs - eps)) / (len(maximum) + 1.0) for obs, eps in zip(observed, tol)]
    )


def holm(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    if finite.size == 0:
        return out
    order = finite[np.argsort(p[finite])]
    m = len(order)
    adjusted = np.maximum.accumulate((m - np.arange(m)) * p[order])
    out[order] = np.minimum(adjusted, 1.0)
    return out


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    out = np.full(len(p), np.nan)
    finite = np.flatnonzero(np.isfinite(p))
    if finite.size == 0:
        return out
    order = finite[np.argsort(p[finite])]
    m = len(order)
    scaled = p[order] * m / np.arange(1, m + 1)
    out[order] = np.minimum(np.minimum.accumulate(scaled[::-1])[::-1], 1.0)
    return out


# --------------------------------------------------------------------------- classification


@dataclass(frozen=True)
class PersistenceRules:
    """Machine form of the preregistered classification (read from the YAML at runtime)."""

    terminal_ci_level: float
    require_terminal_sign_match: bool
    persist_alpha: float
    persist_ci_level: float
    decay_alpha: float
    decay_ci_level: float
    reversed_ci_level: float | None
    margin_kappa: float
    interaction_kappa: float
    label_order: tuple[str, ...]
    source: str
    text: Mapping[str, str] = field(default_factory=dict)

    def as_dict(self) -> dict[str, object]:
        return {
            "terminal_ci_level": self.terminal_ci_level,
            "require_terminal_sign_match": self.require_terminal_sign_match,
            "persist_alpha": self.persist_alpha,
            "persist_ci_level": self.persist_ci_level,
            "decay_alpha": self.decay_alpha,
            "decay_ci_level": self.decay_ci_level,
            "reversed_ci_level": self.reversed_ci_level,
            "margin_kappa": self.margin_kappa,
            "interaction_kappa": self.interaction_kappa,
            "label_order": list(self.label_order),
            "source": self.source,
            "text": dict(self.text),
        }


class RuleParseError(ValueError):
    """The preregistered classification text does not match the supported grammar."""


_NUM = r"([0-9]*\.?[0-9]+)"


def _level(value: str) -> float:
    number = float(value)
    return number / 100.0 if number > 1 else number


def _search(pattern: str, text: str, label: str) -> re.Match[str]:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if match is None:
        raise RuleParseError(f"{label}: cannot parse rule text {text!r}")
    return match


def _kappa(text: object, label: str) -> float:
    if isinstance(text, Mapping) and "kappa" in text:
        return float(text["kappa"])
    string = str(text).replace(" ", "")
    if re.fullmatch(r"theta_LAR-theta_ISS-T", string, flags=re.IGNORECASE):
        return 1.0
    match = re.fullmatch(r"theta_LAR-" + _NUM + r"\*theta_ISS-T", string, flags=re.IGNORECASE)
    if match is None:
        raise RuleParseError(f"contrast {label}: cannot parse {text!r}")
    return float(match.group(1))


def parse_classification_rules(prereg: Mapping[str, object], *, source: str = "prereg") -> PersistenceRules:
    """Read the §4.2 ``classification`` and ``contrasts`` blocks into ``PersistenceRules``.

    Each label may be the §4.2 sentence (parsed with a strict grammar; anything that
    does not parse raises) or a mapping with explicit keys
    (``hc3_ci_level``, ``require_sign_match_pooled_terminal``, ``fl_one_sided_alpha``,
    ``reversed_theta_lar_ci_level``). The label set must be exactly the five labels.
    """
    block = prereg.get("classification")
    if not isinstance(block, Mapping):
        raise RuleParseError("preregistration has no 'classification' mapping")
    labels = tuple(str(label) for label in block)
    if set(labels) != set(LABELS):
        raise RuleParseError(f"classification labels must be exactly {LABELS}; found {labels}")
    # Mapping order carries no meaning in YAML; precedence is the fixed §4.2 order (LABELS).
    labels = LABELS
    text = {label: block[label] for label in labels}

    tech = text[NOT_EVALUABLE_TECHNICAL]
    if not isinstance(tech, Mapping) and not re.search(r"flag", str(tech), flags=re.IGNORECASE):
        raise RuleParseError(f"{NOT_EVALUABLE_TECHNICAL}: expected a QC-flag rule, got {tech!r}")

    terminal = text[NOT_EVALUABLE_NO_TERMINAL_EFFECT]
    if isinstance(terminal, Mapping):
        terminal_level = _level(str(terminal["hc3_ci_level"]))
        sign_match = bool(terminal.get("require_sign_match_pooled_terminal", True))
    else:
        match = _search(r"HC3\s+" + _NUM + r"\s*%\s*CI\s+of\s+theta_ISS-T\s+includes\s+0", str(terminal), NOT_EVALUABLE_NO_TERMINAL_EFFECT)
        terminal_level = _level(match.group(1))
        sign_match = re.search(r"sign\s+differs\s+from\s+pooled\s+terminal", str(terminal), flags=re.IGNORECASE) is not None

    def margin_rule(label: str, direction: str, bound: str, comparator: str) -> tuple[float, float]:
        rule = text[label]
        if isinstance(rule, Mapping):
            return float(rule["fl_one_sided_alpha"]), _level(str(rule["hc3_ci_level"]))
        cmp_dir = re.escape(direction)
        match = _search(
            r"FL\s+one-sided\s+p\(\s*margin\s*" + cmp_dir + r"\s*0\s*\)\s*<=\s*" + _NUM
            + r"\s+AND\s+HC3\s+" + _NUM + r"\s*%\s*CI\s+" + bound + r"\(\s*margin\s*\)\s*"
            + re.escape(comparator) + r"\s*0",
            str(rule),
            label,
        )
        return float(match.group(1)), _level(match.group(2))

    persist_alpha, persist_level = margin_rule(PERSISTS, ">", "lower", ">")
    decay_alpha, decay_level = margin_rule(DECAYS, "<", "upper", "<")
    decay_text = text[DECAYS]
    if isinstance(decay_text, Mapping):
        reversed_level = decay_text.get("reversed_theta_lar_ci_level")
        reversed_level = None if reversed_level is None else _level(str(reversed_level))
    else:
        match = re.search(
            r"REVERSED\s+if\s+theta_LAR\s+" + _NUM + r"\s*%\s*CI\s+upper\s*<\s*0",
            str(decay_text),
            flags=re.IGNORECASE,
        )
        reversed_level = None if match is None else _level(match.group(1))
    inconclusive = text[INCONCLUSIVE]
    if not isinstance(inconclusive, Mapping) and "otherwise" not in str(inconclusive).lower():
        raise RuleParseError(f"{INCONCLUSIVE}: expected the 'otherwise' catch-all, got {inconclusive!r}")

    contrasts = prereg.get("contrasts") or PLAN_4_2_CONTRASTS
    if not isinstance(contrasts, Mapping) or "margin" not in contrasts or "interaction" not in contrasts:
        raise RuleParseError("preregistration 'contrasts' must define interaction and margin")
    return PersistenceRules(
        terminal_ci_level=terminal_level,
        require_terminal_sign_match=sign_match,
        persist_alpha=persist_alpha,
        persist_ci_level=persist_level,
        decay_alpha=decay_alpha,
        decay_ci_level=decay_level,
        reversed_ci_level=reversed_level,
        margin_kappa=_kappa(contrasts["margin"], "margin"),
        interaction_kappa=_kappa(contrasts["interaction"], "interaction"),
        label_order=labels,
        source=source,
        text={label: str(value) for label, value in text.items()},
    )


def builtin_plan_rules() -> PersistenceRules:
    """The §4.2 draft rules; only for ``--power-only`` before the YAML is committed."""
    return parse_classification_rules(
        {"classification": PLAN_4_2_CLASSIFICATION, "contrasts": PLAN_4_2_CONTRASTS},
        source="builtin_plan_section_4.2_draft",
    )


_CLASSIFIER_FIELDS = (
    "technical_flag",
    "theta_isst",
    "se_isst",
    "df_isst",
    "theta_lar",
    "se_lar",
    "df_lar",
    "margin",
    "margin_se",
    "margin_df",
    "p_persist",
    "p_decay",
    "pooled_terminal_estimate",
)


def classify_persistence_arrays(values: Mapping[str, object], rules: PersistenceRules) -> dict[str, np.ndarray]:
    """Vectorized decision rule (single source of truth for the script and the power simulation).

    Inputs are in terminal orientation. Precedence: NOT_EVALUABLE_TECHNICAL,
    NOT_EVALUABLE_NO_TERMINAL_EFFECT, PERSISTS, DECAYS, INCONCLUSIVE. PERSISTS/DECAYS need
    the Freedman--Lane one-sided p AND the HC3 margin interval; exactly one of the two is
    an FL/HC3 disagreement and yields INCONCLUSIVE.
    """
    missing = [name for name in _CLASSIFIER_FIELDS if name not in values]
    if missing:
        raise KeyError(f"classifier inputs missing: {missing}")
    arr = {name: np.atleast_1d(np.asarray(values[name], dtype=float)) for name in _CLASSIFIER_FIELDS}
    n = max(len(v) for v in arr.values())
    arr = {name: np.broadcast_to(v, (n,)) for name, v in arr.items()}
    technical = arr["technical_flag"] > 0
    q_term = _t_quantile(1.0 - (1.0 - rules.terminal_ci_level) / 2.0, arr["df_isst"])
    t_low = arr["theta_isst"] - q_term * arr["se_isst"]
    t_high = arr["theta_isst"] + q_term * arr["se_isst"]
    terminal_missing = ~(np.isfinite(t_low) & np.isfinite(t_high))
    includes_zero = terminal_missing | ((t_low <= 0.0) & (t_high >= 0.0))
    sign_mismatch = np.zeros(n, dtype=bool)
    if rules.require_terminal_sign_match:
        sign_mismatch = np.sign(arr["theta_isst"]) != np.sign(arr["pooled_terminal_estimate"])
    no_terminal = includes_zero | sign_mismatch
    q_p = _t_quantile(1.0 - (1.0 - rules.persist_ci_level) / 2.0, arr["margin_df"])
    q_d = _t_quantile(1.0 - (1.0 - rules.decay_ci_level) / 2.0, arr["margin_df"])
    margin_low = arr["margin"] - q_p * arr["margin_se"]
    margin_high = arr["margin"] + q_d * arr["margin_se"]
    with np.errstate(invalid="ignore"):
        fl_persist = arr["p_persist"] <= rules.persist_alpha
        ci_persist = margin_low > 0.0
        fl_decay = arr["p_decay"] <= rules.decay_alpha
        ci_decay = margin_high < 0.0
    margin_missing = ~(np.isfinite(margin_low) & np.isfinite(margin_high) & np.isfinite(arr["p_persist"]) & np.isfinite(arr["p_decay"]))
    persists = fl_persist & ci_persist
    decays = fl_decay & ci_decay
    disagree = (fl_persist ^ ci_persist) | (fl_decay ^ ci_decay)
    lar_high = np.full(n, np.nan)
    reversed_flag = np.zeros(n, dtype=bool)
    if rules.reversed_ci_level is not None:
        q_r = _t_quantile(1.0 - (1.0 - rules.reversed_ci_level) / 2.0, arr["df_lar"])
        lar_high = arr["theta_lar"] + q_r * arr["se_lar"]
        with np.errstate(invalid="ignore"):
            reversed_flag = decays & (lar_high < 0.0)
    code = np.full(n, 4, dtype=int)
    code[decays] = 3
    code[persists] = 2
    code[no_terminal] = 1
    code[technical] = 0
    labels = np.array([LABELS[c] for c in code], dtype=object)
    labels[(code == 3) & reversed_flag] = DECAYS_REVERSED
    return {
        "label": labels,
        "code": code,
        "technical": technical,
        "terminal_low": t_low,
        "terminal_high": t_high,
        "terminal_missing": terminal_missing,
        "includes_zero": includes_zero,
        "sign_mismatch": sign_mismatch,
        "margin_low": margin_low,
        "margin_high": margin_high,
        "fl_persist": fl_persist,
        "ci_persist": ci_persist,
        "fl_decay": fl_decay,
        "ci_decay": ci_decay,
        "disagree": disagree & (code == 4),
        "margin_missing": margin_missing,
        "lar_high": lar_high,
        "reversed": reversed_flag & (code == 3),
    }


def classify_persistence(row: Mapping[str, object], rules: PersistenceRules) -> tuple[str, str]:
    """(label, reason) for one target under the preregistered rules (inputs terminal-oriented)."""
    d = {k: v[0] for k, v in classify_persistence_arrays(row, rules).items()}
    pct = lambda level: f"{100 * level:g}%"  # noqa: E731
    if d["code"] == 0:
        return NOT_EVALUABLE_TECHNICAL, "LAR primary coverage metric flagged by the frozen technical-QC rule"
    if d["code"] == 1:
        if d["terminal_missing"]:
            return NOT_EVALUABLE_NO_TERMINAL_EFFECT, "theta_ISS-T interval not estimable"
        if d["includes_zero"]:
            return NOT_EVALUABLE_NO_TERMINAL_EFFECT, (
                f"HC3 {pct(rules.terminal_ci_level)} CI of theta_ISS-T "
                f"[{d['terminal_low']:.4g}, {d['terminal_high']:.4g}] includes 0"
            )
        return NOT_EVALUABLE_NO_TERMINAL_EFFECT, (
            f"theta_ISS-T ({float(row['theta_isst']):.4g}) sign differs from the pooled terminal "
            f"estimate ({float(row['pooled_terminal_estimate']):.4g})"
        )
    p_persist, p_decay = float(row["p_persist"]), float(row["p_decay"])
    if d["code"] == 2:
        return PERSISTS, (
            f"FL p(margin>0)={p_persist:.4g}<={rules.persist_alpha:g} and HC3 "
            f"{pct(rules.persist_ci_level)} CI lower(margin)={d['margin_low']:.4g}>0"
        )
    if d["code"] == 3:
        reason = (
            f"FL p(margin<0)={p_decay:.4g}<={rules.decay_alpha:g} and HC3 "
            f"{pct(rules.decay_ci_level)} CI upper(margin)={d['margin_high']:.4g}<0"
        )
        if d["reversed"]:
            reason += (
                f"; REVERSED: theta_LAR {pct(rules.reversed_ci_level)} CI upper="
                f"{d['lar_high']:.4g}<0"
            )
            return DECAYS_REVERSED, reason
        return DECAYS, reason
    if d["margin_missing"]:
        return INCONCLUSIVE, "margin test inputs not available"
    if d["disagree"]:
        parts = []
        if d["fl_persist"] != d["ci_persist"]:
            parts.append(
                f"persist direction: FL p={p_persist:.4g} ({'<=' if d['fl_persist'] else '>'}"
                f"{rules.persist_alpha:g}) vs HC3 lower={d['margin_low']:.4g}"
            )
        if d["fl_decay"] != d["ci_decay"]:
            parts.append(
                f"decay direction: FL p={p_decay:.4g} ({'<=' if d['fl_decay'] else '>'}"
                f"{rules.decay_alpha:g}) vs HC3 upper={d['margin_high']:.4g}"
            )
        return INCONCLUSIVE, "FL/HC3 disagreement (" + "; ".join(parts) + ")"
    return INCONCLUSIVE, (
        f"neither margin test met: FL p(margin>0)={p_persist:.4g}, p(margin<0)={p_decay:.4g}; "
        f"HC3 margin CI [{d['margin_low']:.4g}, {d['margin_high']:.4g}] includes 0"
    )


def classify_g_scale(margin_g: float, margin_g_var: float, primary_label: str, rules: PersistenceRules) -> str:
    """Descriptive g-scale companion: normal CI of g_L - k g_T at the persist/decay levels."""
    base = primary_label.split(":")[0]
    if base in (NOT_EVALUABLE_TECHNICAL, NOT_EVALUABLE_NO_TERMINAL_EFFECT):
        return base
    if not (np.isfinite(margin_g) and np.isfinite(margin_g_var)):
        return INCONCLUSIVE
    se = math.sqrt(margin_g_var)
    low = margin_g - stats.norm.ppf(1.0 - (1.0 - rules.persist_ci_level) / 2.0) * se
    high = margin_g + stats.norm.ppf(1.0 - (1.0 - rules.decay_ci_level) / 2.0) * se
    if low > 0:
        return PERSISTS
    if high < 0:
        return DECAYS
    return INCONCLUSIVE


_CONSERVATIVE_RANK = {
    NOT_EVALUABLE_TECHNICAL: 0,
    NOT_EVALUABLE_NO_TERMINAL_EFFECT: 1,
    INCONCLUSIVE: 2,
    PERSISTS: 3,
    DECAYS: 3,
}


def conservative_verdict(labels: Sequence[str]) -> tuple[str, bool]:
    """Headline label across gene-ID maps (APPENDIX A item 3) and whether it is robust.

    Agreement keeps the label; any flip is 'not robust' and takes the most conservative
    label (NOT_EVALUABLE < INCONCLUSIVE < PERSISTS/DECAYS; PERSISTS vs DECAYS -> INCONCLUSIVE).
    """
    labels = [str(label) for label in labels]
    if not labels:
        raise ValueError("no labels")
    if len(set(labels)) == 1:
        return labels[0], True
    bases = [label.split(":")[0] for label in labels]
    rank = min(_CONSERVATIVE_RANK[base] for base in bases)
    if rank == 3:  # PERSISTS vs DECAYS, or DECAYS vs DECAYS:REVERSED
        if set(bases) == {DECAYS}:
            return DECAYS, False
        return INCONCLUSIVE, False
    candidates = [base for base in bases if _CONSERVATIVE_RANK[base] == rank]
    return candidates[0], False


# --------------------------------------------------------------------------- power


def synthetic_meta40(n_per_cell: int = 5) -> pd.DataFrame:
    rows = []
    for arm in ARMS:
        for age in AGES:
            for condition in ("FLT", "GC"):
                for rep in range(n_per_cell):
                    rows.append(
                        {
                            "animal": f"{arm}_{age}_{condition}_{rep}",
                            "condition": condition,
                            "arm": arm,
                            "age": age,
                            "block": f"{arm}_{age}",
                            "source_sample": f"{arm}_{age}_{condition}_{rep}",
                        }
                    )
    return pd.DataFrame(rows).set_index("animal")


def analytic_power(theta_isst: float, rho: float, n_per_cell: int = 5, alpha: float = 0.05, kappa: float = 0.5) -> dict[str, float]:
    """Normal-approximation powers in score-SD units (the §4.2 analytic numbers)."""
    se_theta = math.sqrt(2.0 / n_per_cell / len(AGES))
    se_int = math.sqrt(2.0) * se_theta
    se_margin = math.sqrt(1.0 + kappa**2) * se_theta
    theta_lar = rho * theta_isst
    interaction = theta_lar - theta_isst
    margin = theta_lar - kappa * theta_isst
    z2 = stats.norm.ppf(1.0 - alpha / 2.0)
    z1 = stats.norm.ppf(1.0 - alpha)
    return {
        "analytic_se_theta": se_theta,
        "analytic_se_interaction": se_int,
        "analytic_se_margin": se_margin,
        "analytic_power_interaction": float(
            stats.norm.cdf(abs(interaction) / se_int - z2) + stats.norm.cdf(-abs(interaction) / se_int - z2)
        ),
        "analytic_power_margin_persist": float(stats.norm.cdf(margin / se_margin - z1)),
        "analytic_power_margin_decay": float(stats.norm.cdf(-margin / se_margin - z1)),
    }


def simulate_power(
    scenarios: Iterable[tuple[float, float]],
    n_per_cell: int = 5,
    n_sim: int = 4000,
    n_perm: int = 999,
    seed: object = 20260811 + 702,
    *,
    rules: PersistenceRules | None = None,
    alpha_interaction: float = 0.05,
    chunk: int = 1024,
    block_means: Sequence[float] = (0.3, -0.2, 0.1, 0.4),
) -> pd.DataFrame:
    """Apply the real HC3 + Freedman--Lane decision rule to N(0, 1) noise plus block means.

    Scenarios are (theta_ISS-T, rho) with theta_LAR = rho theta_ISS-T, in score-SD units.
    Simulated datasets are the columns of one outcome matrix, so they share the
    permutation set within a scenario (each p-value is still a valid FL p-value; this
    only correlates Monte-Carlo errors). The pooled terminal sign is taken as positive.
    """
    rules = rules or builtin_plan_rules()
    meta = synthetic_meta40(n_per_cell)
    design = design_matrix(meta)
    blocks = design[list(BLOCK_COLUMNS)].to_numpy(dtype=float)
    f_t = design["F_T"].to_numpy()
    f_l = design["F_L"].to_numpy()
    x = design.to_numpy(dtype=float)
    iss_mask = meta["arm"].eq("ISS-T").to_numpy()
    df_t, df_l = arm_residual_df(x, iss_mask), arm_residual_df(x, ~iss_mask)
    scenarios = [(float(t), float(r)) for t, r in scenarios]
    children = np.random.SeedSequence(seed).spawn(len(scenarios))
    rows = []
    for index, ((theta_t, rho), child) in enumerate(zip(scenarios, children)):
        data_seed, perm_seed = child.spawn(2)
        rng = np.random.default_rng(data_seed)
        theta_l = rho * theta_t
        mean = blocks @ np.asarray(block_means, dtype=float) + theta_t * f_t + theta_l * f_l
        y = mean[:, None] + rng.standard_normal((len(meta), n_sim))
        fit = ols_hc3(x, y)
        th_t, th_l = fit.beta[4], fit.beta[5]
        v_t, v_l = fit.cov_hc3[:, 4, 4], fit.cov_hc3[:, 5, 5]
        int_df = satterthwaite_df(v_t, df_t, v_l, df_l, rules.interaction_kappa)
        int_t = (th_l - rules.interaction_kappa * th_t) / np.sqrt(v_l + rules.interaction_kappa**2 * v_t)
        p_int_hc3 = 2.0 * stats.t.sf(np.abs(int_t), int_df)
        margin = th_l - rules.margin_kappa * th_t
        margin_se = np.sqrt(v_l + rules.margin_kappa**2 * v_t)
        margin_df = satterthwaite_df(v_t, df_t, v_l, df_l, rules.margin_kappa)
        fl_int = freedman_lane_contrast(
            y, blocks, f_t, f_l, rules.interaction_kappa, meta["block"].to_numpy(), n_perm, perm_seed, chunk, keep_null=False
        )
        fl_margin = freedman_lane_contrast(
            y, blocks, f_t, f_l, rules.margin_kappa, meta["block"].to_numpy(), n_perm, perm_seed, chunk, keep_null=False
        )
        decision = classify_persistence_arrays(
            {
                "technical_flag": 0.0,
                "theta_isst": th_t,
                "se_isst": np.sqrt(v_t),
                "df_isst": float(df_t),
                "theta_lar": th_l,
                "se_lar": np.sqrt(v_l),
                "df_lar": float(df_l),
                "margin": margin,
                "margin_se": margin_se,
                "margin_df": margin_df,
                "p_persist": fl_margin.p_upper,
                "p_decay": fl_margin.p_lower,
                "pooled_terminal_estimate": 1.0,
            },
            rules,
        )
        labels = pd.Series(decision["label"])
        established = decision["code"] >= 2
        base = labels.str.split(":").str[0]
        row = {
            "scenario": index,
            "theta_isst": theta_t,
            "rho": rho,
            "theta_lar": theta_l,
            "n_per_cell": n_per_cell,
            "n_sim": n_sim,
            "n_perm": n_perm,
            "power_interaction_fl": float(np.mean(fl_int.p_two_sided <= alpha_interaction)),
            "power_interaction_hc3": float(np.mean(p_int_hc3 <= alpha_interaction)),
            "power_margin_persist_fl": float(np.mean(fl_margin.p_upper <= rules.persist_alpha)),
            "power_margin_persist_fl_and_hc3": float(np.mean(decision["fl_persist"] & decision["ci_persist"])),
            "power_margin_decay_fl": float(np.mean(fl_margin.p_lower <= rules.decay_alpha)),
            "power_margin_decay_fl_and_hc3": float(np.mean(decision["fl_decay"] & decision["ci_decay"])),
            "p_isst_effect_established": float(np.mean(established)),
        }
        for label in LABELS:
            row[f"p_label_{label}"] = float(np.mean(base == label))
        row[f"p_label_{DECAYS_REVERSED}"] = float(np.mean(labels == DECAYS_REVERSED))
        row["p_label_FL_HC3_disagreement"] = float(np.mean(decision["disagree"]))
        with np.errstate(invalid="ignore", divide="ignore"):
            row["p_persists_given_established"] = float(np.mean(base[established] == PERSISTS)) if established.any() else np.nan
            row["p_decays_given_established"] = float(np.mean(base[established] == DECAYS)) if established.any() else np.nan
        row["modal_label"] = str(base.value_counts().idxmax())
        row["mc_se_max"] = math.sqrt(0.25 / n_sim)
        row.update(analytic_power(theta_t, rho, n_per_cell, alpha_interaction, rules.margin_kappa))
        rows.append(row)
    table = pd.DataFrame(rows)
    table.attrs["seed"] = str(seed)
    return table
