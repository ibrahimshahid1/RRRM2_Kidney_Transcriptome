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

A blocked flight-label permutation (labels permuted within mission x stratum,
structural/podocyte covariates held fixed) calibrates a max-|T| family-wise test
over {P1, S1, D}. Leave-one-mission repeats the adjusted synthesis.

If P1 keeps a positive flight coefficient whose interval excludes zero and it
survives the max-|T| family, the podocyte program is separable from structural
drift and the "recurrent podocyte-associated program" title stands. If P1
collapses toward zero (as barrier-core did against the disjoint podocyte proxy,
0.507 -> 0.070), the honest headline is a broad structural-expression response.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
import yaml

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.data import cpm_eligible_genes, load_primary_missions  # noqa: E402
from src.clinical_axes.statistics import (  # noqa: E402
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
    """Precomputed per-mission quantities for a fast blocked permutation of the
    flight label under one adjusted model (nuisance columns Z held fixed)."""

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
        out = np.zeros((batch, self.n), dtype=float)
        for positions, nt in self.blocks:
            m = len(positions)
            keys = rng.random((batch, m))
            chosen = np.argpartition(keys, nt - 1, axis=1)[:, :nt] if nt > 0 else np.empty((batch, 0), int)
            for b in range(batch):
                out[b, positions[chosen[b]]] = 1.0
        return out


def _pooled_t(betas: np.ndarray, variances: np.ndarray) -> float:
    """REML/mKH pooled t across missions for one permutation draw (k small)."""
    fit = random_effects_reml_mkh(betas, variances)
    return fit.t


def run(args):
    config = yaml.safe_load(args.config.read_text())
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
        if variance_type == "HC3" and model in ("P1_podocyte_adjusted", "S1_structural_adjusted"):
            lomo = leave_one_mission_out(
                sub["estimate"].to_numpy(), sub["variance"].to_numpy(), list(sub.index)
            )
            lomo.insert(0, "model", model)
            lomo_rows.append(lomo)
    meta = pd.DataFrame(meta_rows)

    # --- blocked permutation null with max-|T| over the family ---
    n_perm = int(args.permutations)
    rng_root = np.random.SeedSequence(args.seed)
    # observed model-based pooled t per family model (from closed form, for exact
    # comparability with the null)
    observed_t = {}
    for model in FAMILY:
        objs = perm_objects[model]
        betas = np.array([o.coeff_var(o.flight_obs[None, :])[0][0] for o in objs])
        varis = np.array([o.coeff_var(o.flight_obs[None, :])[1][0] for o in objs])
        observed_t[model] = _pooled_t(betas, varis)

    null_t = {model: np.empty(n_perm) for model in FAMILY}
    # one independent rng stream per (model, mission) so blocks permute independently
    child_seeds = rng_root.spawn(len(FAMILY) * len(mission_order))
    rngs = {}
    si = 0
    for model in FAMILY:
        for mi in range(len(mission_order)):
            rngs[(model, mi)] = np.random.default_rng(child_seeds[si])
            si += 1

    chunk = int(args.chunk_size)
    for model in FAMILY:
        objs = perm_objects[model]
        for start in range(0, n_perm, chunk):
            stop = min(start + chunk, n_perm)
            batch = stop - start
            betas = np.empty((len(objs), batch))
            varis = np.empty((len(objs), batch))
            for mi, o in enumerate(objs):
                flights = o.permuted_flights(rngs[(model, mi)], batch)
                b, v = o.coeff_var(flights)
                betas[mi] = b
                varis[mi] = v
            # pool each permutation: reuse the batched REML/mKH by transposing to
            # (batch x k) then computing pooled t row-wise
            from src.clinical_axes.statistics import _batch_reml_mkh  # noqa: E402
            _, _, _, t_vals = _batch_reml_mkh(betas.T, varis.T)
            null_t[model][start:stop] = t_vals

    null_max = np.max(np.abs(np.column_stack([null_t[m] for m in FAMILY])), axis=1)
    perm_rows = []
    for model in FAMILY:
        obs = observed_t[model]
        emp = (1.0 + np.count_nonzero(np.abs(null_t[model]) >= abs(obs))) / (n_perm + 1.0)
        fwer = (1.0 + np.count_nonzero(null_max >= abs(obs))) / (n_perm + 1.0)
        perm_rows.append(
            {
                "model": model,
                "observed_t_model_based": obs,
                "empirical_p_two_sided": emp,
                "max_t_fwer": fwer,
                "n_permutations": n_perm,
            }
        )
    permutation = pd.DataFrame(perm_rows)
    lomo_table = pd.concat(lomo_rows, ignore_index=True) if lomo_rows else pd.DataFrame()

    # --- write ---
    out = args.results
    out.mkdir(parents=True, exist_ok=True)
    effects.to_csv(out / "podocyte_scaffold_mission_effects.tsv", sep="\t", index=False)
    meta.to_csv(out / "podocyte_scaffold_meta.tsv", sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(out / "podocyte_scaffold_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(score_rows).to_csv(out / "podocyte_scaffold_scores.tsv", sep="\t", index=False)
    permutation.to_csv(out / "podocyte_scaffold_permutation.tsv", sep="\t", index=False)
    if not lomo_table.empty:
        lomo_table.to_csv(out / "podocyte_scaffold_leave_one_mission.tsv", sep="\t", index=False)
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
        "seed": int(args.seed),
        "family_maxT": FAMILY,
    }
    (out / "podocyte_scaffold_manifest.json").write_text(json.dumps(manifest, indent=2))

    # --- console summary ---
    print("=== K0: podocyte vs structural-scaffold specificity ===")
    print(f"podocyte {len(podocyte_all)} genes; structural {len(structural_all)} "
          f"(minus {len(overlap)} overlap = {len(structural_disjoint)} disjoint)")
    show = meta[meta["variance_type"] == "HC3"].set_index("model")
    for label, *_ in MODELS:
        r = show.loc[label]
        print(f"  {label:26s} g={r['estimate']:+.3f}  95%CI [{r['ci_low_mkh']:+.3f}, "
              f"{r['ci_high_mkh']:+.3f}]  p={r['p_mkh']:.4f}  I2={r['i_squared']:.0f}%")
    print("  --- blocked permutation (model-based t, maxT over family) ---")
    for _, r in permutation.iterrows():
        print(f"  {r['model']:26s} emp_p={r['empirical_p_two_sided']:.4f}  "
              f"maxT_FWER={r['max_t_fwer']:.4f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--tiers", type=Path, default=DEFAULT_TIERS)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--permutations", type=int, default=100_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--chunk-size", type=int, default=2_048)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
