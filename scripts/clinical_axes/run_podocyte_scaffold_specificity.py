#!/usr/bin/env python3
"""K0: is the cross-mission podocyte flight signal separable from the broad
structural-scaffold program, or is it the upper tail of a shared drift?

For every terminal mission we score two label-blind programs per animal
(within-mission signed gene z-scores, equal weight):

  * podocyte      = podocyte__high_specificity            (157 frozen genes)
  * structural    = broad_structural_scaffold_control__all with the genes it
                    shares with the podocyte set removed  (disjoint comparator)

We then fit, within each mission, HC3-robust regressions blocked by the mission's
exchangeability stratum and pool the flight coefficient across missions with
REML / modified Hartung-Knapp:

  P0  podocyte   ~ flight + block                 (unadjusted)
  P1  podocyte   ~ flight + block + structural     (adjusted)          <- decisive
  S0  structural ~ flight + block                 (unadjusted)
  S1  structural ~ flight + block + podocyte       (reciprocal adjusted)
  D   (podocyte - structural) ~ flight + block     (direct contrast)

The pooled flight coefficients are raw score-unit OLS coefficients
(``coef_score_units``); they are NOT Hedges g. The ``estimate`` columns of the
legacy tables carry the same numbers and are kept for backward compatibility.

A blocked permutation (within mission x stratum) calibrates a max-|T| family-wise
test over {P1, S1, D}. Three schemes are emitted side by side:

  legacy_independent_streams  the pre-correction code path, kept for audit. One RNG
                              per (model, mission): each family model is shuffled
                              with its own label permutations and the "max" is taken
                              over independent draws, which is not a joint null
                              (by Sidak's inequality it is conservative).
  joint_label                 one RNG per mission; the same within-block label
                              permutation is applied to every family model
                              (permute-X: adjuster held fixed).
  freedman_lane_joint         declared primary going forward. Freedman & Lane (1983):
                              residualize the outcome on the reduced design
                              Z = (const, block, adjuster), permute the residuals
                              within blocks with the same pi per mission across
                              P1, S1 and D, and refit the full model. Unlike
                              permute-X it respects that flight also shifts the
                              adjuster (Anderson & Legendre 1999; Winkler et al. 2014).

joint_label and freedman_lane_joint share the same within-block permutations
(common random numbers). Their per-mission streams are SeedSequence(seed).spawn(k),
which are also the legacy P1 streams, so the P1 marginal null of joint_label is
bit-identical to legacy's and only the coupling across models changes.

Leave-one-mission repeats the HC3 synthesis for P1, S1 and D. The K0 verdict rests
on the HC3 interval for P1, which no permutation scheme changes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions  # noqa: E402
from src.clinical_axes.statistics import (  # noqa: E402
    _batch_reml_mkh,
    leave_one_mission_out,
    random_effects_reml_mkh,
    score_signed_axis,
)

DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_TIERS = REPO / "data/processed/v13_compartment_audit/frozen_compartment_tiers.tsv"
DEFAULT_RESULTS = REPO / "data/results/podocyte_scaffold_specificity"

PODOCYTE_SET = "podocyte__high_specificity"
STRUCTURAL_SET = "broad_structural_scaffold_control__all"

# The three tests that form the max-|T| family. Each is (label, outcome, adjuster).
# outcome/adjuster name a scored program, or "diff" for the podocyte-minus-structural
# contrast; adjuster None means unadjusted.
MODELS = [
    ("P0_podocyte_unadjusted", "podocyte", None),
    ("P1_podocyte_adjusted", "podocyte", "structural"),
    ("S0_structural_unadjusted", "structural", None),
    ("S1_structural_adjusted", "structural", "podocyte"),
    ("D_direct_difference", "diff", None),
]
FAMILY = ["P1_podocyte_adjusted", "S1_structural_adjusted", "D_direct_difference"]
# HC3 leave-one-mission syntheses (D added by the K0 permutation correction).
LOMO_MODELS = ("P1_podocyte_adjusted", "S1_structural_adjusted", "D_direct_difference")

LEGACY_SCHEME = "legacy_independent_streams"
JOINT_LABEL_SCHEME = "joint_label"
FREEDMAN_LANE_SCHEME = "freedman_lane_joint"
PERMUTATION_SCHEMES = (LEGACY_SCHEME, JOINT_LABEL_SCHEME, FREEDMAN_LANE_SCHEME)
PRIMARY_PERMUTATION_SCHEME = FREEDMAN_LANE_SCHEME
SCHEME_DESCRIPTIONS = {
    LEGACY_SCHEME: (
        "pre-correction path kept for audit: within-block flight-label permutation "
        "(adjuster fixed) with one independent RNG per (model, mission); max-|T| is "
        "taken over independent draws, not a joint null"
    ),
    JOINT_LABEL_SCHEME: (
        "within-block flight-label permutation (adjuster fixed), one RNG per mission, "
        "the same permutation applied to every family model"
    ),
    FREEDMAN_LANE_SCHEME: (
        "Freedman-Lane: residuals of outcome ~ const + block + adjuster permuted within "
        "blocks, the same permutation per mission across family models; declared primary"
    ),
}
COEFFICIENT_UNITS = "coef_score_units"


def _frozen_genes(tiers: pd.DataFrame, gene_set: str) -> list[str]:
    mask = (tiers["gene_set"] == gene_set) & tiers["final_for_testing"].astype(bool)
    return list(dict.fromkeys(tiers.loc[mask, "gene_symbol"].astype(str)))


def _design(metadata: pd.DataFrame, adjuster: pd.Series | None) -> pd.DataFrame:
    design = pd.DataFrame(
        {"flight": (metadata["condition"] == "FLT").astype(float)},
        index=metadata.index,
    )
    block = pd.get_dummies(metadata["block"], drop_first=True, dtype=float)
    design = design.join(block)
    if adjuster is not None:
        design["adjuster"] = adjuster.loc[metadata.index].astype(float)
    return sm.add_constant(design, has_constant="add")


def _regression_rows(mission, metadata, outcome, adjuster, model_label):
    design = _design(metadata, adjuster)
    fit = sm.OLS(outcome, design).fit()
    robust = fit.get_robustcov_results(cov_type="HC3")
    idx = list(design.columns).index("flight")
    rows = []
    for variance_type, f in (("model_based", fit), ("HC3", robust)):
        est = float(f.params[idx] if variance_type == "HC3" else f.params["flight"])
        se = float(f.bse[idx] if variance_type == "HC3" else f.bse["flight"])
        rows.append(
            {
                "mission": mission,
                "model": model_label,
                "variance_type": variance_type,
                "estimate": est,
                "standard_error": se,
                "variance": se**2,
                "p": float(f.pvalues[idx] if variance_type == "HC3" else f.pvalues["flight"]),
                "n": int(len(outcome)),
            }
        )
    return rows


# --- closed-form blocked permutation (Frisch-Waugh-Lovell) ---------------------

class _MissionPerm:
    """Precomputed per-mission quantities for a fast blocked permutation of one
    adjusted model (nuisance columns Z = const + block dummies [+ adjuster]).

    Two permutation statistics are available:

    * ``coeff_var(flights)``: permute-X. The flight indicator is permuted within
      blocks while y and Z stay fixed.
    * ``freedman_lane_coeff_var(perm)``: Freedman-Lane. The reduced-model residuals
      e = M_Z y are permuted within blocks, y* = Z g_hat + e[pi], and the full model
      is refitted. For pi = identity this is the observed model-based OLS fit.
    """

    def __init__(self, metadata: pd.DataFrame, outcome: pd.Series, adjuster: pd.Series | None):
        design = _design(metadata, adjuster)
        # Nuisance Z = everything except the flight column.
        z = design.drop(columns=["flight"]).to_numpy(dtype=float)
        y = outcome.to_numpy(dtype=float)
        n = len(y)
        # Residual maker M_Z = I - Z (Z'Z)^+ Z'
        zpz_pinv = np.linalg.pinv(z.T @ z)
        hat = z @ zpz_pinv @ z.T
        m_z = np.eye(n) - hat
        self.m_z = m_z
        self.y_res = m_z @ y  # residualized outcome (fixed across permutations)
        self.yy = float(self.y_res @ self.y_res)
        self.p_full = int(np.linalg.matrix_rank(design.to_numpy(dtype=float)))
        self.n = n
        self.dof = n - self.p_full
        # block structure for within-block label permutation
        blocks = []
        for _, sub in metadata.groupby("block", sort=False):
            positions = [metadata.index.get_loc(i) for i in sub.index]
            nt = int((sub["condition"] == "FLT").sum())
            blocks.append((np.array(positions), nt))
        self.blocks = blocks
        # observed flight vector
        self.flight_obs = (metadata["condition"] == "FLT").to_numpy(dtype=float)
        # Freedman-Lane: observed residualized flight f~ = M_Z f and f~'f~ (fixed).
        self.f_res_obs = m_z @ self.flight_obs
        self.den_obs = float(self.f_res_obs @ self.f_res_obs)
        # Local orderings that put each block's observed flight units first; used to
        # map random keys onto a within-block permutation consistent with
        # permuted_flights (the nt smallest keys mark the permuted flight units).
        self._flight_first = [
            np.concatenate(
                [
                    np.flatnonzero(self.flight_obs[positions] == 1.0),
                    np.flatnonzero(self.flight_obs[positions] != 1.0),
                ]
            )
            for positions, _ in blocks
        ]
        covered = np.sort(np.concatenate([p for p, _ in blocks])) if blocks else np.array([])
        self._blocks_cover_all = bool(
            len(covered) == n and np.array_equal(covered, np.arange(n))
        )

    def coeff_var(self, flight_batch: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Vectorized flight coefficient and its (model-based) variance for a
        (batch x n) matrix of flight indicators."""
        f_res = flight_batch @ self.m_z            # (batch x n); M_Z symmetric
        den = np.einsum("bi,bi->b", f_res, f_res)  # f~'f~
        num = f_res @ self.y_res                   # f~'y~
        beta = num / den
        rss = self.yy - num**2 / den
        sigma2 = rss / self.dof
        var = sigma2 / den
        return beta, var

    def permuted_flights(self, rng: np.random.Generator, batch: int) -> np.ndarray:
        """(batch x n) within-block flight-label permutations: in each block the
        units with the nt smallest uniform keys are labelled flight."""
        out = np.zeros((batch, self.n), dtype=float)
        rows = np.arange(batch)[:, None]
        for positions, nt in self.blocks:
            m = len(positions)
            keys = rng.random((batch, m))
            if nt <= 0:
                continue
            chosen = np.argpartition(keys, nt - 1, axis=1)[:, :nt]
            out[rows, positions[chosen]] = 1.0
        return out

    def permutation_indices(self, rng: np.random.Generator, batch: int) -> np.ndarray:
        """(batch x n) uniform within-block permutations pi.

        Consumes exactly the same random keys as ``permuted_flights`` (one
        ``rng.random((batch, m))`` per block, in block order). pi maps the observed
        flight units of a block onto the units with the nt smallest keys, so for the
        same generator state ``flights_for_permutation(pi)`` equals
        ``permuted_flights`` (up to exact key ties, probability ~0).
        """
        if not self._blocks_cover_all:
            raise ValueError("blocks must partition all animals (missing block labels?)")
        perm = np.empty((batch, self.n), dtype=np.intp)
        for (positions, _), flight_first in zip(self.blocks, self._flight_first):
            m = len(positions)
            keys = rng.random((batch, m))
            order = np.argsort(keys, axis=1)
            local = np.empty((batch, m), dtype=np.intp)
            local[:, flight_first] = order  # pi(flight_first[j]) = order[j]
            perm[:, positions] = positions[local]
        return perm

    def flights_for_permutation(self, perm: np.ndarray) -> np.ndarray:
        """Permute-X flight matrix for pi: unit pi(i) receives unit i's label."""
        out = np.empty(perm.shape, dtype=float)
        out[np.arange(perm.shape[0])[:, None], perm] = self.flight_obs[None, :]
        return out

    def freedman_lane_coeff_var(self, perm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Freedman-Lane flight coefficient and model-based variance for a
        (batch x n) matrix of within-block permutations pi.

        R = M_Z e[pi]; num = R'f~; den = f~'f~; rss = R'R - num^2/den;
        var = rss / dof / den.
        """
        resid = self.y_res[perm] @ self.m_z        # (batch x n); M_Z symmetric
        num = resid @ self.f_res_obs
        beta = num / self.den_obs
        rss = np.einsum("bi,bi->b", resid, resid) - num**2 / self.den_obs
        var = rss / self.dof / self.den_obs
        return beta, var


def _pooled_t(betas: np.ndarray, variances: np.ndarray) -> float:
    """REML/mKH pooled t across missions for one permutation draw (k small)."""
    fit = random_effects_reml_mkh(betas, variances)
    return fit.t


def _validate_perm_objects(perm_objects: Mapping[str, Sequence[_MissionPerm]]) -> tuple[list[str], int]:
    models = list(perm_objects)
    if not models:
        raise ValueError("perm_objects must contain at least one family model")
    n_missions = len(perm_objects[models[0]])
    if n_missions < 2:
        raise ValueError("cross-mission permutation inference requires at least two missions")
    for model in models:
        objs = perm_objects[model]
        if len(objs) != n_missions:
            raise ValueError("every family model must have one _MissionPerm per mission")
        for mi, obj in enumerate(objs):
            ref = perm_objects[models[0]][mi]
            same = obj.n == ref.n and len(obj.blocks) == len(ref.blocks) and all(
                np.array_equal(p1, p2) and n1 == n2
                for (p1, n1), (p2, n2) in zip(obj.blocks, ref.blocks)
            ) and np.array_equal(obj.flight_obs, ref.flight_obs)
            if not same:
                raise ValueError(
                    f"mission {mi}: family models must share animals, blocks and labels "
                    "for a joint permutation"
                )
    return models, n_missions


def observed_family_t(perm_objects: Mapping[str, Sequence[_MissionPerm]]) -> pd.DataFrame:
    """Observed model-based REML/mKH pooled flight coefficient and t per family
    model, from the same closed form as the permutation null (pi = identity)."""
    rows = []
    for model, objs in perm_objects.items():
        betas = np.array([o.coeff_var(o.flight_obs[None, :])[0][0] for o in objs])
        varis = np.array([o.coeff_var(o.flight_obs[None, :])[1][0] for o in objs])
        fit = random_effects_reml_mkh(betas, varis)
        rows.append(
            {
                "model": model,
                "observed_t_model_based": fit.t,
                f"observed_{COEFFICIENT_UNITS}_model_based": fit.estimate,
            }
        )
    return pd.DataFrame(rows).set_index("model")


def family_permutation_null(
    perm_objects: Mapping[str, Sequence[_MissionPerm]],
    *,
    n_permutations: int,
    seed: int,
    chunk_size: int = 2_048,
    schemes: Sequence[str] = PERMUTATION_SCHEMES,
) -> dict[str, pd.DataFrame]:
    """Pooled REML/mKH null t (n_permutations x models) per permutation scheme.

    legacy_independent_streams reproduces the pre-correction implementation
    bit-for-bit for the same seed and chunk size: SeedSequence(seed) spawns one
    stream per (model, mission), model-major. The joint schemes spawn one stream per
    mission from a fresh SeedSequence(seed); every chunk draws one within-block
    permutation per mission and applies it to every model (joint_label as a
    permuted flight label, freedman_lane_joint as permuted reduced-model residuals).
    """
    unknown = set(schemes) - set(PERMUTATION_SCHEMES)
    if unknown:
        raise ValueError(f"unknown permutation schemes: {sorted(unknown)}")
    if n_permutations < 1:
        raise ValueError("n_permutations must be positive")
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    models, n_missions = _validate_perm_objects(perm_objects)
    null = {s: {m: np.empty(n_permutations) for m in models} for s in schemes}

    if LEGACY_SCHEME in schemes:
        # one independent rng stream per (model, mission), model-major (pre-correction)
        child_seeds = np.random.SeedSequence(seed).spawn(len(models) * n_missions)
        rngs = {}
        si = 0
        for model in models:
            for mi in range(n_missions):
                rngs[(model, mi)] = np.random.default_rng(child_seeds[si])
                si += 1
        for model in models:
            objs = perm_objects[model]
            for start in range(0, n_permutations, chunk_size):
                stop = min(start + chunk_size, n_permutations)
                batch = stop - start
                betas = np.empty((n_missions, batch))
                varis = np.empty((n_missions, batch))
                for mi, o in enumerate(objs):
                    flights = o.permuted_flights(rngs[(model, mi)], batch)
                    betas[mi], varis[mi] = o.coeff_var(flights)
                _, _, _, t_vals = _batch_reml_mkh(betas.T, varis.T)
                null[LEGACY_SCHEME][model][start:stop] = t_vals

    joint = [s for s in schemes if s != LEGACY_SCHEME]
    if joint:
        mission_rngs = [
            np.random.default_rng(child)
            for child in np.random.SeedSequence(seed).spawn(n_missions)
        ]
        for start in range(0, n_permutations, chunk_size):
            stop = min(start + chunk_size, n_permutations)
            batch = stop - start
            perms = [
                perm_objects[models[0]][mi].permutation_indices(mission_rngs[mi], batch)
                for mi in range(n_missions)
            ]
            for scheme in joint:
                for model in models:
                    betas = np.empty((n_missions, batch))
                    varis = np.empty((n_missions, batch))
                    for mi, o in enumerate(perm_objects[model]):
                        if scheme == JOINT_LABEL_SCHEME:
                            betas[mi], varis[mi] = o.coeff_var(
                                o.flights_for_permutation(perms[mi])
                            )
                        else:
                            betas[mi], varis[mi] = o.freedman_lane_coeff_var(perms[mi])
                    _, _, _, t_vals = _batch_reml_mkh(betas.T, varis.T)
                    null[scheme][model][start:stop] = t_vals

    return {s: pd.DataFrame(null[s], columns=models) for s in schemes}


def scheme_seed_record(seed: int, n_models: int, n_missions: int) -> dict[str, dict[str, object]]:
    """Manifest description of the RNG streams behind each scheme."""
    return {
        LEGACY_SCHEME: {
            "seed": int(seed),
            "rng": f"SeedSequence({int(seed)}).spawn({n_models * n_missions}); "
                   "one stream per (model, mission), model-major",
            "primary": False,
            "description": SCHEME_DESCRIPTIONS[LEGACY_SCHEME],
        },
        JOINT_LABEL_SCHEME: {
            "seed": int(seed),
            "rng": f"SeedSequence({int(seed)}).spawn({n_missions}); one stream per mission "
                   "shared by all family models (streams equal the legacy P1 streams)",
            "primary": False,
            "description": SCHEME_DESCRIPTIONS[JOINT_LABEL_SCHEME],
        },
        FREEDMAN_LANE_SCHEME: {
            "seed": int(seed),
            "rng": f"SeedSequence({int(seed)}).spawn({n_missions}); same within-block "
                   "permutations as joint_label (common random numbers)",
            "primary": True,
            "description": SCHEME_DESCRIPTIONS[FREEDMAN_LANE_SCHEME],
        },
    }


def permutation_table(
    observed: pd.DataFrame,
    null_t: Mapping[str, pd.DataFrame],
    *,
    seed: int,
) -> pd.DataFrame:
    """Per-scheme marginal and max-|T| family-wise plus-one p-values."""
    rows = []
    for scheme, null in null_t.items():
        n_perm = int(len(null))
        models = list(null.columns)
        null_max = np.max(np.abs(null.to_numpy()), axis=1)
        q95 = float(np.quantile(null_max, 0.95))
        for model in models:
            obs = float(observed.loc[model, "observed_t_model_based"])
            emp = (1.0 + np.count_nonzero(np.abs(null[model].to_numpy()) >= abs(obs))) / (n_perm + 1.0)
            fwer = (1.0 + np.count_nonzero(null_max >= abs(obs))) / (n_perm + 1.0)
            rows.append(
                {
                    "permutation_scheme": scheme,
                    "model": model,
                    "observed_t_model_based": obs,
                    "empirical_p_two_sided": emp,
                    "max_t_fwer": fwer,
                    "n_permutations": n_perm,
                    f"observed_{COEFFICIENT_UNITS}_model_based": float(
                        observed.loc[model, f"observed_{COEFFICIENT_UNITS}_model_based"]
                    ),
                    "max_t_fwer_mc_se": float(np.sqrt(fwer * (1.0 - fwer) / n_perm)),
                    "null_max_abs_t_q95": q95,
                    "family_size": len(models),
                    "seed": int(seed),
                    "primary_scheme": scheme == PRIMARY_PERMUTATION_SCHEME,
                }
            )
    return pd.DataFrame(rows)


def _sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def apply_gene_map_override(config: dict, gene_map: Path | None) -> dict:
    """Point config["gene_mapping"]["path"] at an alternative ID map, in memory only
    (the config file is never edited). Returns the config for chaining."""
    if gene_map is not None:
        resolved = Path(gene_map).expanduser().resolve()
        if not resolved.exists():
            raise FileNotFoundError(f"--gene-map {resolved} does not exist")
        config["gene_mapping"]["path"] = str(resolved)
    return config


def _with_coef_column(frame: pd.DataFrame) -> pd.DataFrame:
    """Add coef_score_units (= estimate) right after estimate; units are raw
    score-unit OLS coefficients, never Hedges g."""
    frame = frame.copy()
    frame.insert(frame.columns.get_loc("estimate") + 1, COEFFICIENT_UNITS, frame["estimate"])
    return frame


def run(args):
    config = yaml.safe_load(args.config.read_text())
    config = apply_gene_map_override(config, getattr(args, "gene_map", None))
    gene_map_path = REPO / str(config["gene_mapping"]["path"])
    missions = load_primary_missions(config, REPO)
    threshold = float(config["eligibility"]["cpm_threshold"])
    tiers = pd.read_csv(args.tiers, sep="\t")

    podocyte_all = _frozen_genes(tiers, PODOCYTE_SET)
    structural_all = _frozen_genes(tiers, STRUCTURAL_SET)
    overlap = sorted(set(podocyte_all) & set(structural_all))
    structural_disjoint = [g for g in structural_all if g not in set(podocyte_all)]

    effect_rows, coverage_rows, score_rows = [], [], []
    # per-mission closed-form permutation objects for each family model
    perm_objects: dict[str, list[_MissionPerm]] = {m: [] for m in FAMILY}
    mission_order = list(missions.keys())

    for mission, data in missions.items():
        eligible = set(cpm_eligible_genes(data, threshold))
        pod_used = [g for g in podocyte_all if g in eligible and g in data.expression.index]
        struc_used = [g for g in structural_disjoint if g in eligible and g in data.expression.index]
        if len(pod_used) < 8 or len(struc_used) < 8:
            raise RuntimeError(f"{mission}: podocyte={len(pod_used)}, structural={len(struc_used)}")
        podocyte = score_signed_axis(data.expression, {g: 1 for g in pod_used}).scores
        structural = score_signed_axis(data.expression, {g: 1 for g in struc_used}).scores
        diff = (podocyte - structural).rename("signed_axis_score")

        scored = {"podocyte": podocyte, "structural": structural, "diff": diff}
        coverage_rows.append(
            {
                "mission": mission,
                "n_podocyte_genes": len(pod_used),
                "n_structural_disjoint_genes": len(struc_used),
                "podocyte_structural_score_pearson_r": float(podocyte.corr(structural)),
            }
        )
        for label, outcome_name, adjuster_name in MODELS:
            outcome = scored[outcome_name]
            adjuster = scored[adjuster_name] if adjuster_name else None
            effect_rows.extend(
                _regression_rows(mission, data.metadata, outcome, adjuster, label)
            )
            if label in FAMILY:
                perm_objects[label].append(
                    _MissionPerm(data.metadata, outcome, adjuster)
                )
        for animal in data.metadata.index:
            score_rows.append(
                {
                    "mission": mission,
                    "animal": animal,
                    "condition": data.metadata.loc[animal, "condition"],
                    "block": data.metadata.loc[animal, "block"],
                    "podocyte_score": float(podocyte.loc[animal]),
                    "structural_score": float(structural.loc[animal]),
                }
            )

    effects = pd.DataFrame(effect_rows)

    # --- pooled REML/mKH per model (HC3 primary, model_based secondary) ---
    meta_rows, lomo_rows = [], []
    for (model, variance_type), sub in effects.groupby(["model", "variance_type"], sort=False):
        sub = sub.set_index("mission").loc[mission_order]
        fit = random_effects_reml_mkh(sub["estimate"].to_numpy(), sub["variance"].to_numpy())
        meta_rows.append(
            {
                "model": model,
                "variance_type": variance_type,
                "estimate": fit.estimate,
                "ci_low_mkh": fit.ci_low,
                "ci_high_mkh": fit.ci_high,
                "p_mkh": fit.p,
                "t_mkh": fit.t,
                "tau2": fit.tau2,
                "i_squared": fit.i_squared,
                "prediction_low": fit.prediction_low,
                "prediction_high": fit.prediction_high,
                "maximum_weight": float(fit.weights.max()),
                "k_missions": fit.k,
            }
        )
        if variance_type == "HC3" and model in LOMO_MODELS:
            lomo = leave_one_mission_out(
                sub["estimate"].to_numpy(), sub["variance"].to_numpy(), list(sub.index)
            )
            lomo.insert(0, "model", model)
            lomo_rows.append(lomo)
    meta = pd.DataFrame(meta_rows)

    # --- blocked permutation null with max-|T| over the family, three schemes ---
    n_perm = int(args.permutations)
    seed = int(args.seed)
    observed = observed_family_t(perm_objects)
    null_t = family_permutation_null(
        perm_objects,
        n_permutations=n_perm,
        seed=seed,
        chunk_size=int(args.chunk_size),
        schemes=PERMUTATION_SCHEMES,
    )
    permutation = permutation_table(observed, null_t, seed=seed)
    lomo_table = pd.concat(lomo_rows, ignore_index=True) if lomo_rows else pd.DataFrame()

    # --- write ---
    out = args.results
    out.mkdir(parents=True, exist_ok=True)
    _with_coef_column(effects).to_csv(
        out / "podocyte_scaffold_mission_effects.tsv", sep="\t", index=False
    )
    _with_coef_column(meta).to_csv(out / "podocyte_scaffold_meta.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(out / "podocyte_scaffold_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(score_rows).to_csv(out / "podocyte_scaffold_scores.tsv", sep="\t", index=False)
    permutation.to_csv(out / "podocyte_scaffold_permutation.tsv", sep="\t", index=False)
    if not lomo_table.empty:
        _with_coef_column(lomo_table).to_csv(
            out / "podocyte_scaffold_leave_one_mission.tsv", sep="\t", index=False
        )
    manifest = {
        "podocyte_set": PODOCYTE_SET,
        "structural_set": STRUCTURAL_SET,
        "n_podocyte_frozen": len(podocyte_all),
        "n_structural_frozen": len(structural_all),
        "overlap_genes_removed_from_structural": overlap,
        "n_structural_disjoint": len(structural_disjoint),
        "cpm_threshold": threshold,
        "missions": mission_order,
        "n_permutations": n_perm,
        "seed": seed,
        "chunk_size": int(args.chunk_size),
        "family_maxT": FAMILY,
        "estimate_units": (
            f"{COEFFICIENT_UNITS}: raw score-unit OLS flight coefficients (not Hedges g); "
            "'estimate' columns hold the same values for backward compatibility"
        ),
        "permutation_schemes": scheme_seed_record(seed, len(FAMILY), len(mission_order)),
        "primary_permutation_scheme": PRIMARY_PERMUTATION_SCHEME,
        "permutation_statistic": "REML/mKH pooled t of model-based per-mission flight coefficients",
        "leave_one_mission_models": list(LOMO_MODELS),
        "verdict_rule": (
            "K0 verdict rests on the HC3 REML/mKH interval for P1; no permutation scheme "
            "can upgrade it (prereg k0_rule)"
        ),
        "config": str(args.config),
        "config_sha256": _sha256(Path(args.config)),
        "tiers": str(args.tiers),
        "tiers_sha256": _sha256(Path(args.tiers)),
        "gene_map": str(gene_map_path),
        "gene_map_sha256": _sha256(gene_map_path),
        "gene_map_overridden": getattr(args, "gene_map", None) is not None,
    }
    (out / "podocyte_scaffold_manifest.json").write_text(json.dumps(manifest, indent=2))

    # --- console summary ---
    print("=== K0: podocyte vs structural-scaffold specificity ===")
    print(f"podocyte {len(podocyte_all)} genes; structural {len(structural_all)} "
          f"(minus {len(overlap)} overlap = {len(structural_disjoint)} disjoint)")
    show = meta[meta["variance_type"] == "HC3"].set_index("model")
    for label, *_ in MODELS:
        r = show.loc[label]
        print(f"  {label:26s} coef={r['estimate']:+.3f}  95%CI [{r['ci_low_mkh']:+.3f}, "
              f"{r['ci_high_mkh']:+.3f}]  p={r['p_mkh']:.4f}  I2={r['i_squared']:.0f}%")
    print("  --- blocked permutation (model-based t, maxT over family) ---")
    for _, r in permutation.iterrows():
        flag = " (primary)" if r["primary_scheme"] else ""
        print(f"  [{r['permutation_scheme']}{flag}] {r['model']:26s} "
              f"emp_p={r['empirical_p_two_sided']:.4f}  maxT_FWER={r['max_t_fwer']:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--chunk-size", type=int, default=2_048)
    parser.add_argument(
        "--gene-map",
        type=Path,
        default=None,
        help=(
            "Alternative Ensembl->symbol ID map; replaces config gene_mapping.path in "
            "memory only (the config file is not edited)."
        ),
    )
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
