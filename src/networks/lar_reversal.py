"""LAR reversal and mechanism-switch analysis utilities."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from src.common import bh_fdr, id_map_lookup
from src.networks.mechanism_axis import safe_z


EPS = 1e-12

# Independent state-space coordinates used for cosine/projection geometry.
STATE_VECTOR_FEATURES = (
    "matrix_component",
    "dct_transport_component",
    "immune_context_component",
    "preservation_stress_component",
)

STATE_SCALAR_FEATURES = ("matrix_minus_dct",)
STATE_COMPONENT_FEATURES = STATE_VECTOR_FEATURES + STATE_SCALAR_FEATURES

# Backward-compatible name for callers that expect the vector feature set.
STATE_FEATURES = STATE_VECTOR_FEATURES

MECHANISM_FEATURES = (
    "ecm_organization",
    "fibrosis_tgfb_emt",
    "tlr4_innate",
    "integrin_cell_adhesion",
    "s1p_s1pr3",
    "mmp_adam_proteolysis",
    "dct_ncc_wnk_transport",
    "tubular_transport_broad",
    "oxidative_stress_nrf2",
    "macrophage_inflammation",
    "preservation_stress_response",
)

CORE_CLOCK_GENES = (
    "Arntl", "Clock", "Npas2", "Per1", "Per2", "Per3", "Cry1", "Cry2",
    "Nr1d1", "Nr1d2", "Rora", "Rorc", "Dbp", "Tef", "Hlf",
)

PER_CRY_GENES = ("Per1", "Per2", "Per3", "Cry1", "Cry2")
BMAL_CLOCK_GENES = ("Arntl", "Clock", "Npas2")

DCT_WNK_GENES = (
    "Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Klhl3", "Cul3",
    "Calb1", "Kcnj10", "Kcnj16", "Sgk1", "Nedd4l",
)

ALDOSTERONE_ENAC_GENES = ("Scnn1a", "Scnn1b", "Scnn1g", "Nr3c2", "Hsd11b2", "Sgk1")

S1P_AXIS_GENES = (
    "Sphk1", "Sphk2", "Sgpl1", "Sgpp1", "Sgpp2", "S1pr1", "S1pr2",
    "S1pr3", "S1pr4", "S1pr5", "Spns2", "Lpar1", "Lpar2",
)

RBM3_PRESERVATION_GENES = (
    "Rbm3", "Cirbp", "Ddit3", "Hspa1a", "Hspa1b", "Hsp90aa1",
    "Fos", "Jun", "Egr1", "Atf3", "Ier2", "Btg2",
)

TARGETED_PANELS: Mapping[str, tuple[str, ...]] = {
    "clock_core": CORE_CLOCK_GENES,
    "per_cry": PER_CRY_GENES,
    "bmal_clock": BMAL_CLOCK_GENES,
    "dct_wnk_expression": DCT_WNK_GENES,
    "aldosterone_enac_expression": ALDOSTERONE_ENAC_GENES,
    "s1p_axis_expression": S1P_AXIS_GENES,
    "rbm3_preservation_expression": RBM3_PRESERVATION_GENES,
}

SINGLE_GENE_FEATURES = ("Rbm3", "Per2", "S1pr3", "Slc12a3", "Wnk4")


@dataclass(frozen=True)
class VectorStats:
    """Observed vector relationship for a feature set."""

    feature_set: str
    age_scope: str
    n_features: int
    iss_norm: float
    lar_norm: float
    cos_lar_iss: float
    cos_lar_negative_iss: float
    beta_lar_to_iss: float
    lar_to_iss_magnitude_ratio: float


def finite_common(a: pd.Series, b: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Return aligned finite entries from two vectors."""
    common = sorted(set(a.index) & set(b.index))
    aa = a.loc[common].astype(float)
    bb = b.loc[common].astype(float)
    mask = np.isfinite(aa.to_numpy()) & np.isfinite(bb.to_numpy())
    return aa.loc[mask], bb.loc[mask]


def cosine_similarity(a: pd.Series | Sequence[float], b: pd.Series | Sequence[float]) -> float:
    """Cosine similarity with finite-value masking."""
    if isinstance(a, pd.Series) and isinstance(b, pd.Series):
        aa, bb = finite_common(a, b)
        av = aa.to_numpy(dtype=float)
        bv = bb.to_numpy(dtype=float)
    else:
        av = np.asarray(a, dtype=float)
        bv = np.asarray(b, dtype=float)
        mask = np.isfinite(av) & np.isfinite(bv)
        av = av[mask]
        bv = bv[mask]
    if len(av) == 0:
        return float("nan")
    denom = float(np.linalg.norm(av) * np.linalg.norm(bv))
    if denom < EPS:
        return float("nan")
    return float(np.clip(np.dot(av, bv) / denom, -1.0, 1.0))


def projection_beta(query: pd.Series, reference: pd.Series) -> float:
    """Projection coefficient of query onto reference."""
    q, r = finite_common(query, reference)
    rv = r.to_numpy(dtype=float)
    denom = float(np.dot(rv, rv))
    if denom < EPS:
        return float("nan")
    return float(np.dot(q.to_numpy(dtype=float), rv) / denom)


def vector_stats(feature_set: str, age_scope: str, iss: pd.Series, lar: pd.Series) -> VectorStats:
    """Summarize LAR versus ISS-T vector geometry."""
    iss2, lar2 = finite_common(iss, lar)
    iss_norm = float(np.linalg.norm(iss2.to_numpy(dtype=float)))
    lar_norm = float(np.linalg.norm(lar2.to_numpy(dtype=float)))
    cos_li = cosine_similarity(lar2, iss2)
    beta = projection_beta(lar2, iss2)
    return VectorStats(
        feature_set=feature_set,
        age_scope=age_scope,
        n_features=int(len(iss2)),
        iss_norm=iss_norm,
        lar_norm=lar_norm,
        cos_lar_iss=cos_li,
        cos_lar_negative_iss=-cos_li if np.isfinite(cos_li) else np.nan,
        beta_lar_to_iss=beta,
        lar_to_iss_magnitude_ratio=lar_norm / iss_norm if iss_norm > EPS else np.nan,
    )


def classify_vector_relationship(row: Mapping[str, object]) -> str:
    """Translate vector geometry into the three-model vocabulary."""
    cos_li = float(row.get("cos_lar_iss", np.nan))
    ci_high = float(row.get("bootstrap_ci_high", np.nan))
    beta = float(row.get("beta_lar_to_iss", np.nan))
    ratio = float(row.get("lar_to_iss_magnitude_ratio", np.nan))
    if np.isfinite(cos_li) and np.isfinite(ci_high) and ci_high < 0 and beta <= -0.25:
        return "model_B_reversal_candidate"
    if np.isfinite(beta) and np.isfinite(ratio) and abs(beta) < 0.25 and ratio < 0.60:
        return "model_A_attenuation_candidate"
    if np.isfinite(cos_li) and cos_li < -0.20:
        return "model_C_partial_reversal_or_switch"
    if np.isfinite(cos_li) and cos_li > 0.25 and np.isfinite(beta) and beta > 0:
        return "same_direction_candidate"
    return "mixed_or_inconclusive"


def rrrm2_flt_gc_rows(df: pd.DataFrame, *, age_scope: str = "pooled") -> pd.DataFrame:
    """Return RRRM-2 FLT/GC rows for vector effects."""
    rows = df[
        df["study"].astype(str).eq("RRRM-2")
        & df["condition"].astype(str).isin(["FLT", "GC"])
        & df["Arm"].astype(str).isin(["ISS-T", "LAR"])
    ].copy()
    if age_scope != "pooled":
        rows = rows[rows["Age"].astype(str).eq(age_scope)].copy()
    return rows


def arm_effect(
    rows: pd.DataFrame,
    feature_cols: Sequence[str],
    arm: str,
    *,
    age_scope: str = "pooled",
) -> pd.Series:
    """Age-balanced FLT-minus-GC vector for one RRRM-2 arm."""
    arm_rows = rows[rows["Arm"].astype(str).eq(arm)].copy()
    if age_scope != "pooled":
        arm_rows = arm_rows[arm_rows["Age"].astype(str).eq(age_scope)].copy()
    diffs: list[pd.Series] = []
    group_iter = arm_rows.groupby("Age", dropna=False) if age_scope == "pooled" else [(age_scope, arm_rows)]
    for _, sub in group_iter:
        flt = sub[sub["condition"].astype(str).eq("FLT")]
        gc = sub[sub["condition"].astype(str).eq("GC")]
        if flt.empty or gc.empty:
            continue
        diffs.append(flt.loc[:, feature_cols].astype(float).mean(axis=0) - gc.loc[:, feature_cols].astype(float).mean(axis=0))
    if not diffs:
        return pd.Series(np.nan, index=list(feature_cols), dtype=float)
    return pd.concat(diffs, axis=1).mean(axis=1)


def osd513_effect(df: pd.DataFrame, feature_cols: Sequence[str]) -> pd.Series:
    """OSD-513 FLT-minus-GC vector."""
    rows = df[df["study"].astype(str).eq("OSD-513") & df["condition"].astype(str).isin(["FLT", "GC"])].copy()
    flt = rows[rows["condition"].astype(str).eq("FLT")]
    gc = rows[rows["condition"].astype(str).eq("GC")]
    if flt.empty or gc.empty:
        return pd.Series(np.nan, index=list(feature_cols), dtype=float)
    return flt.loc[:, feature_cols].astype(float).mean(axis=0) - gc.loc[:, feature_cols].astype(float).mean(axis=0)


def bootstrap_rrrm2_rows(rows: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Resample RRRM-2 rows within Arm x Age x condition cells."""
    parts: list[pd.DataFrame] = []
    for _, idx in rows.groupby(["Arm", "Age", "condition"], dropna=False).groups.items():
        idx = list(idx)
        chosen = rng.choice(idx, size=len(idx), replace=True)
        parts.append(rows.loc[chosen])
    return pd.concat(parts, ignore_index=True)


def permute_rrrm2_conditions(rows: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Permute FLT/GC labels within Arm x Age strata."""
    out = rows.copy()
    for _, idx in out.groupby(["Arm", "Age"], dropna=False).groups.items():
        labels = out.loc[idx, "condition"].to_numpy(copy=True)
        rng.shuffle(labels)
        out.loc[idx, "condition"] = labels
    return out


def bootstrap_external_rows(rows: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Resample external rows within condition cells."""
    parts: list[pd.DataFrame] = []
    group_cols = [c for c in ["condition"] if c in rows.columns]
    for _, idx in rows.groupby(group_cols, dropna=False).groups.items():
        idx = list(idx)
        chosen = rng.choice(idx, size=len(idx), replace=True)
        parts.append(rows.loc[chosen])
    return pd.concat(parts, ignore_index=True)


def permute_external_conditions(rows: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    """Permute external FLT/GC labels."""
    out = rows.copy()
    labels = out["condition"].to_numpy(copy=True)
    rng.shuffle(labels)
    out["condition"] = labels
    return out


def reversal_summary_for_features(
    df: pd.DataFrame,
    feature_cols: Sequence[str],
    feature_set: str,
    *,
    age_scope: str,
    n_bootstrap: int,
    n_permutation: int,
    rng: np.random.Generator,
    include_osd513: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute LAR-vs-ISS-T reversal summary for one feature matrix."""
    feature_cols = [c for c in feature_cols if c in df.columns]
    rows = rrrm2_flt_gc_rows(df, age_scope=age_scope)
    values = rows.loc[:, feature_cols].astype(float).to_numpy()
    arms = rows["Arm"].astype(str).to_numpy()
    ages = rows["Age"].astype(str).to_numpy()
    conditions = rows["condition"].astype(str).to_numpy()

    def effect_np(
        vals: np.ndarray,
        arm_labels: np.ndarray,
        age_labels: np.ndarray,
        condition_labels: np.ndarray,
        arm: str,
    ) -> pd.Series:
        age_values = sorted(set(age_labels)) if age_scope == "pooled" else [age_scope]
        diffs: list[np.ndarray] = []
        for age in age_values:
            flt = (arm_labels == arm) & (age_labels == age) & (condition_labels == "FLT")
            gc = (arm_labels == arm) & (age_labels == age) & (condition_labels == "GC")
            if flt.any() and gc.any():
                diffs.append(np.nanmean(vals[flt], axis=0) - np.nanmean(vals[gc], axis=0))
        if not diffs:
            return pd.Series(np.nan, index=feature_cols, dtype=float)
        return pd.Series(np.nanmean(np.vstack(diffs), axis=0), index=feature_cols, dtype=float)

    def boot_indices() -> np.ndarray:
        picked: list[np.ndarray] = []
        keys = pd.DataFrame({"arm": arms, "age": ages, "condition": conditions})
        for _, idx in keys.groupby(["arm", "age", "condition"], dropna=False).groups.items():
            pos = np.asarray(list(idx), dtype=int)
            picked.append(rng.choice(pos, size=len(pos), replace=True))
        return np.concatenate(picked)

    def permuted_conditions() -> np.ndarray:
        out = conditions.copy()
        keys = pd.DataFrame({"arm": arms, "age": ages})
        for _, idx in keys.groupby(["arm", "age"], dropna=False).groups.items():
            pos = np.asarray(list(idx), dtype=int)
            labels = out[pos].copy()
            rng.shuffle(labels)
            out[pos] = labels
        return out

    iss = effect_np(values, arms, ages, conditions, "ISS-T")
    lar = effect_np(values, arms, ages, conditions, "LAR")
    obs = vector_stats(feature_set, age_scope, iss, lar).__dict__

    boot_cos: list[float] = []
    boot_beta: list[float] = []
    for _ in range(n_bootstrap):
        idx = boot_indices()
        b_iss = effect_np(values[idx], arms[idx], ages[idx], conditions[idx], "ISS-T")
        b_lar = effect_np(values[idx], arms[idx], ages[idx], conditions[idx], "LAR")
        boot_cos.append(cosine_similarity(b_lar, b_iss))
        boot_beta.append(projection_beta(b_lar, b_iss))

    perm_cos: list[float] = []
    for _ in range(n_permutation):
        perm_conditions = permuted_conditions()
        p_iss = effect_np(values, arms, ages, perm_conditions, "ISS-T")
        p_lar = effect_np(values, arms, ages, perm_conditions, "LAR")
        perm_cos.append(cosine_similarity(p_lar, p_iss))

    boot = np.asarray(boot_cos, dtype=float)
    perm = np.asarray(perm_cos, dtype=float)
    obs_cos = float(obs["cos_lar_iss"])
    obs.update({
        "bootstrap_median": float(np.nanmedian(boot)),
        "bootstrap_ci_low": float(np.nanquantile(boot, 0.025)),
        "bootstrap_ci_high": float(np.nanquantile(boot, 0.975)),
        "bootstrap_beta_ci_low": float(np.nanquantile(np.asarray(boot_beta, dtype=float), 0.025)),
        "bootstrap_beta_ci_high": float(np.nanquantile(np.asarray(boot_beta, dtype=float), 0.975)),
        "permutation_p_less": float((np.nansum(perm <= obs_cos) + 1) / (np.isfinite(perm).sum() + 1)),
        "permutation_p_two_sided_abs": float((np.nansum(np.abs(perm) >= abs(obs_cos)) + 1) / (np.isfinite(perm).sum() + 1)),
        "n_bootstrap": int(n_bootstrap),
        "n_permutation": int(n_permutation),
    })
    obs["interpretation"] = classify_vector_relationship(obs)

    ext_rows: list[dict[str, object]] = []
    if include_osd513 and df["study"].astype(str).eq("OSD-513").any():
        ext = df[df["study"].astype(str).eq("OSD-513") & df["condition"].astype(str).isin(["FLT", "GC"])].copy()
        ref = osd513_effect(ext, feature_cols)
        cos_lar_ext = cosine_similarity(lar, ref)
        beta_lar_ext = projection_beta(lar, ref)
        ext_rows.append({
            "feature_set": feature_set,
            "age_scope": age_scope,
            "comparison": "LAR_vs_OSD513",
            "cosine": cos_lar_ext,
            "cosine_vs_negative_reference": -cos_lar_ext if np.isfinite(cos_lar_ext) else np.nan,
            "beta_query_to_reference": beta_lar_ext,
            "n_features": int(len(finite_common(lar, ref)[0])),
        })
        cos_iss_ext = cosine_similarity(iss, ref)
        ext_rows.append({
            "feature_set": feature_set,
            "age_scope": age_scope,
            "comparison": "ISS-T_vs_OSD513",
            "cosine": cos_iss_ext,
            "cosine_vs_negative_reference": -cos_iss_ext if np.isfinite(cos_iss_ext) else np.nan,
            "beta_query_to_reference": projection_beta(iss, ref),
            "n_features": int(len(finite_common(iss, ref)[0])),
        })

    return pd.DataFrame([obs]), pd.DataFrame(ext_rows)


def component_effect_rows(
    df: pd.DataFrame,
    feature_cols: Sequence[str],
    *,
    feature_family: str,
) -> pd.DataFrame:
    """Arm- and age-specific FLT-minus-GC effects for sample-level scores."""
    rows = rrrm2_flt_gc_rows(df)
    out: list[dict[str, object]] = []
    scopes = ["pooled", "YNG", "OLD"]
    for feature in [c for c in feature_cols if c in rows.columns]:
        for age_scope in scopes:
            for arm in ("ISS-T", "LAR"):
                effect = arm_effect(rows, [feature], arm, age_scope=age_scope).iloc[0]
                sub = rows[rows["Arm"].astype(str).eq(arm)]
                if age_scope != "pooled":
                    sub = sub[sub["Age"].astype(str).eq(age_scope)]
                flt = sub[sub["condition"].astype(str).eq("FLT")][feature].dropna().astype(float)
                gc = sub[sub["condition"].astype(str).eq("GC")][feature].dropna().astype(float)
                p = float(stats.ttest_ind(flt, gc, equal_var=False).pvalue) if len(flt) >= 2 and len(gc) >= 2 else np.nan
                out.append({
                    "feature_family": feature_family,
                    "feature": feature,
                    "arm": arm,
                    "age_scope": age_scope,
                    "flt_n": int(len(flt)),
                    "gc_n": int(len(gc)),
                    "flt_mean": float(flt.mean()) if len(flt) else np.nan,
                    "gc_mean": float(gc.mean()) if len(gc) else np.nan,
                    "effect_flt_minus_gc": float(effect),
                    "p_welch": p,
                })
    result = pd.DataFrame(out)
    if not result.empty and result["p_welch"].notna().any():
        mask = result["p_welch"].notna()
        result.loc[mask, "q_bh_within_component_effects"] = bh_fdr(result.loc[mask, "p_welch"].to_numpy())
    return result


def interaction_table(
    df: pd.DataFrame,
    feature_cols: Sequence[str],
    *,
    feature_family: str,
    n_bootstrap: int,
    n_permutation: int,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """Arm-by-flight interaction for component scores."""
    base = rrrm2_flt_gc_rows(df)
    features = [c for c in feature_cols if c in base.columns]
    if not features:
        return pd.DataFrame()
    values = base.loc[:, features].astype(float).to_numpy()
    arms = base["Arm"].astype(str).to_numpy()
    ages = base["Age"].astype(str).to_numpy()
    conditions = base["condition"].astype(str).to_numpy()

    def effect_vec(vals: np.ndarray, arm_labels: np.ndarray, age_labels: np.ndarray, condition_labels: np.ndarray, arm: str) -> np.ndarray:
        diffs: list[np.ndarray] = []
        for age in sorted(set(age_labels)):
            flt = (arm_labels == arm) & (age_labels == age) & (condition_labels == "FLT")
            gc = (arm_labels == arm) & (age_labels == age) & (condition_labels == "GC")
            if flt.any() and gc.any():
                diffs.append(np.nanmean(vals[flt], axis=0) - np.nanmean(vals[gc], axis=0))
        if not diffs:
            return np.full(len(features), np.nan, dtype=float)
        return np.nanmean(np.vstack(diffs), axis=0)

    def interaction_vec(vals: np.ndarray, arm_labels: np.ndarray, age_labels: np.ndarray, condition_labels: np.ndarray) -> np.ndarray:
        return effect_vec(vals, arm_labels, age_labels, condition_labels, "ISS-T") - effect_vec(vals, arm_labels, age_labels, condition_labels, "LAR")

    keys_boot = pd.DataFrame({"arm": arms, "age": ages, "condition": conditions})
    boot_groups = [np.asarray(list(idx), dtype=int) for _, idx in keys_boot.groupby(["arm", "age", "condition"], dropna=False).groups.items()]
    keys_perm = pd.DataFrame({"arm": arms, "age": ages})
    perm_groups = [np.asarray(list(idx), dtype=int) for _, idx in keys_perm.groupby(["arm", "age"], dropna=False).groups.items()]

    iss_vec = effect_vec(values, arms, ages, conditions, "ISS-T")
    lar_vec = effect_vec(values, arms, ages, conditions, "LAR")
    obs_vec = iss_vec - lar_vec
    boot = np.zeros((n_bootstrap, len(features)), dtype=float)
    boot_iss = np.zeros((n_bootstrap, len(features)), dtype=float)
    boot_lar = np.zeros((n_bootstrap, len(features)), dtype=float)
    for b in range(n_bootstrap):
        idx = np.concatenate([rng.choice(pos, size=len(pos), replace=True) for pos in boot_groups])
        boot_iss[b, :] = effect_vec(values[idx], arms[idx], ages[idx], conditions[idx], "ISS-T")
        boot_lar[b, :] = effect_vec(values[idx], arms[idx], ages[idx], conditions[idx], "LAR")
        boot[b, :] = boot_iss[b, :] - boot_lar[b, :]
    perm = np.zeros((n_permutation, len(features)), dtype=float)
    for p in range(n_permutation):
        perm_conditions = conditions.copy()
        for pos in perm_groups:
            labels = perm_conditions[pos].copy()
            rng.shuffle(labels)
            perm_conditions[pos] = labels
        perm[p, :] = interaction_vec(values, arms, ages, perm_conditions)

    out: list[dict[str, object]] = []
    for j, feature in enumerate(features):
        interaction = obs_vec[j]
        perm_arr = perm[:, j]
        iss_ci_low = float(np.nanquantile(boot_iss[:, j], 0.025))
        iss_ci_high = float(np.nanquantile(boot_iss[:, j], 0.975))
        lar_ci_low = float(np.nanquantile(boot_lar[:, j], 0.025))
        lar_ci_high = float(np.nanquantile(boot_lar[:, j], 0.975))
        iss_nonzero = (iss_ci_low > 0 and iss_ci_high > 0) or (iss_ci_low < 0 and iss_ci_high < 0)
        lar_nonzero = (lar_ci_low > 0 and lar_ci_high > 0) or (lar_ci_low < 0 and lar_ci_high < 0)
        out.append({
            "feature_family": feature_family,
            "feature": feature,
            "iss_t_effect": float(iss_vec[j]),
            "lar_effect": float(lar_vec[j]),
            "interaction_iss_minus_lar": float(interaction),
            "lar_opposes_iss": bool(np.sign(iss_vec[j]) != np.sign(lar_vec[j]) and iss_nonzero and lar_nonzero),
            "iss_t_bootstrap_ci_low": iss_ci_low,
            "iss_t_bootstrap_ci_high": iss_ci_high,
            "lar_bootstrap_ci_low": lar_ci_low,
            "lar_bootstrap_ci_high": lar_ci_high,
            "bootstrap_ci_low": float(np.nanquantile(boot[:, j], 0.025)),
            "bootstrap_ci_high": float(np.nanquantile(boot[:, j], 0.975)),
            "permutation_p_two_sided": float((np.nansum(np.abs(perm_arr) >= abs(interaction)) + 1) / (np.isfinite(perm_arr).sum() + 1)),
            "n_bootstrap": int(n_bootstrap),
            "n_permutation": int(n_permutation),
        })
    result = pd.DataFrame(out)
    if not result.empty and result["permutation_p_two_sided"].notna().any():
        result["q_bh_within_interactions"] = bh_fdr(result["permutation_p_two_sided"].to_numpy(dtype=float))
    return result


def leave_one_feature_out(iss: pd.Series, lar: pd.Series, *, feature_set: str, age_scope: str) -> pd.DataFrame:
    """Leave-one-feature-out cosine sensitivity for LAR vs ISS-T."""
    ii, ll = finite_common(iss, lar)
    rows = [{
        "feature_set": feature_set,
        "age_scope": age_scope,
        "dropped_feature": "__full__",
        "n_features": int(len(ii)),
        "cos_lar_iss": cosine_similarity(ll, ii),
        "cos_lar_negative_iss": -cosine_similarity(ll, ii),
        "beta_lar_to_iss": projection_beta(ll, ii),
    }]
    for feature in ii.index:
        keep = [f for f in ii.index if f != feature]
        if len(keep) < 2:
            continue
        cos_val = cosine_similarity(ll.loc[keep], ii.loc[keep])
        rows.append({
            "feature_set": feature_set,
            "age_scope": age_scope,
            "dropped_feature": feature,
            "n_features": int(len(keep)),
            "cos_lar_iss": cos_val,
            "cos_lar_negative_iss": -cos_val if np.isfinite(cos_val) else np.nan,
            "beta_lar_to_iss": projection_beta(ll.loc[keep], ii.loc[keep]),
        })
    return pd.DataFrame(rows)


def resolve_symbol_ids(
    symbols: Sequence[str],
    id_map_path: str,
    expression_genes: set[str],
) -> dict[str, list[str]]:
    """Resolve symbols to expression IDs."""
    _ens_to_symbol, symbol_to_ens = id_map_lookup(id_map_path)
    out: dict[str, list[str]] = {}
    for symbol in symbols:
        ids = sorted(set(symbol_to_ens.get(symbol.lower(), set())) & expression_genes)
        out[symbol] = ids
    return out


def score_symbol_panels(
    expression: pd.DataFrame,
    id_map_path: str,
    panels: Mapping[str, Sequence[str]] = TARGETED_PANELS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score targeted panels as sample x score mean gene-wise z-expression."""
    expr = expression.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    expression_genes = set(expr.index)
    rows: list[dict[str, object]] = []
    scores: dict[str, np.ndarray] = {}
    for panel, symbols in panels.items():
        resolved = resolve_symbol_ids(symbols, id_map_path, expression_genes)
        ids = sorted({gid for vals in resolved.values() for gid in vals})
        if ids:
            mat = expr.loc[ids].to_numpy(dtype=float)
            mu = np.nanmean(mat, axis=1, keepdims=True)
            sd = np.nanstd(mat, axis=1, ddof=1, keepdims=True)
            sd[~np.isfinite(sd) | (sd < EPS)] = 1.0
            scores[panel] = np.nanmean((mat - mu) / sd, axis=0)
        for sym, sym_ids in resolved.items():
            rows.append({
                "panel": panel,
                "query_symbol": sym,
                "resolved_gene_ids": ",".join(sym_ids),
                "n_resolved_in_expression": int(len(sym_ids)),
            })
    score_df = pd.DataFrame(scores, index=expr.columns)

    singles = resolve_symbol_ids(SINGLE_GENE_FEATURES, id_map_path, expression_genes)
    for symbol, ids in singles.items():
        if not ids:
            continue
        vals = expr.loc[ids].mean(axis=0).to_numpy(dtype=float)
        score_df[f"{symbol}_expression"] = safe_z(vals)
    return score_df, pd.DataFrame(rows)


def sample_feature_frame(
    metadata: pd.DataFrame,
    scores: pd.DataFrame,
    *,
    sample_col: str = "Sample Name",
    study: str = "RRRM-2",
    scenario: str = "primary",
) -> pd.DataFrame:
    """Join sample-indexed scores to RRRM-2 metadata columns."""
    meta = metadata.copy()
    if sample_col not in meta.columns:
        raise ValueError(f"Metadata missing sample column {sample_col!r}")
    keep = [c for c in [sample_col, "Age", "Arm", "EnvGroup", "condition"] if c in meta.columns]
    out = meta[keep].drop_duplicates(sample_col).copy()
    if "condition" not in out.columns and "EnvGroup" in out.columns:
        out["condition"] = out["EnvGroup"].astype(str)
    out["study"] = study
    out["scenario"] = scenario
    merged = out.merge(scores, left_on=sample_col, right_index=True, how="inner")
    if sample_col != "Sample Name":
        merged = merged.rename(columns={sample_col: "Sample Name"})
    return merged


def gene_effect_scatter(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    id_map_path: str,
    *,
    sample_col: str = "Sample Name",
    priority: pd.DataFrame | None = None,
    panel_symbols: Mapping[str, Sequence[str]] | None = None,
) -> pd.DataFrame:
    """Build ISS-T vs LAR gene-level FLT-minus-GC scatter data."""
    expr = expression.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if expr.index.duplicated().any():
        expr = expr.groupby(expr.index).mean()
    meta = metadata.copy()
    meta[sample_col] = meta[sample_col].astype(str)
    common = [s for s in expr.columns if s in set(meta[sample_col])]
    feature_cols = list(expr.index)
    sample_meta = meta.drop_duplicates(sample_col).set_index(sample_col).loc[common]
    if "condition" not in sample_meta.columns and "EnvGroup" in sample_meta.columns:
        sample_meta["condition"] = sample_meta["EnvGroup"].astype(str)
    values = expr.loc[feature_cols, common].T.to_numpy(dtype=float)
    arms = sample_meta["Arm"].astype(str).to_numpy()
    ages = sample_meta["Age"].astype(str).to_numpy()
    conditions = sample_meta["condition"].astype(str).to_numpy()

    def effect_np(arm: str) -> np.ndarray:
        diffs: list[np.ndarray] = []
        for age in sorted(set(ages)):
            flt = (arms == arm) & (ages == age) & (conditions == "FLT")
            gc = (arms == arm) & (ages == age) & (conditions == "GC")
            if flt.any() and gc.any():
                diffs.append(np.nanmean(values[flt], axis=0) - np.nanmean(values[gc], axis=0))
        if not diffs:
            return np.full(len(feature_cols), np.nan, dtype=float)
        return np.nanmean(np.vstack(diffs), axis=0)

    iss = pd.Series(effect_np("ISS-T"), index=feature_cols)
    lar = pd.Series(effect_np("LAR"), index=feature_cols)
    ens_to_symbol, _symbol_to_ens = id_map_lookup(id_map_path)
    out = pd.DataFrame({
        "gene": feature_cols,
        "mgi_symbol": [ens_to_symbol.get(g, "") for g in feature_cols],
        "iss_t_effect": iss.reindex(feature_cols).to_numpy(dtype=float),
        "lar_effect": lar.reindex(feature_cols).to_numpy(dtype=float),
    })
    out["interaction_iss_minus_lar"] = out["iss_t_effect"] - out["lar_effect"]
    out["reversal_product"] = out["iss_t_effect"] * out["lar_effect"]
    out["quadrant"] = np.select(
        [
            (out["iss_t_effect"] > 0) & (out["lar_effect"] < 0),
            (out["iss_t_effect"] < 0) & (out["lar_effect"] > 0),
            (out["iss_t_effect"] > 0) & (out["lar_effect"] > 0),
            (out["iss_t_effect"] < 0) & (out["lar_effect"] < 0),
        ],
        ["ISS_up_LAR_down", "ISS_down_LAR_up", "same_up", "same_down"],
        default="near_zero_or_mixed",
    )
    panel_symbols = panel_symbols or {}
    panel_lookup: dict[str, list[str]] = {}
    for panel, symbols in panel_symbols.items():
        for symbol in symbols:
            panel_lookup.setdefault(symbol.casefold(), []).append(panel)
    out["highlight_panels"] = out["mgi_symbol"].astype(str).str.casefold().map(
        lambda s: ",".join(sorted(panel_lookup.get(s, [])))
    )
    if priority is not None and not priority.empty and "gene" in priority.columns:
        pr = priority.copy()
        pr["network_priority_rank"] = np.arange(1, len(pr) + 1)
        cols = [
            c for c in [
                "gene", "network_priority_rank", "composite_score", "silent_composite_score",
                "D_ISS", "D_LAR", "N_ISS", "N_LAR", "embedding_specificity",
                "lioness_node_specificity", "OSD513_support",
            ] if c in pr.columns
        ]
        out = out.merge(pr[cols], on="gene", how="left")
    return out.sort_values("interaction_iss_minus_lar", ascending=False, kind="mergesort")


def mechanism_axis_definitions(feature_cols: Sequence[str]) -> pd.DataFrame:
    """Hand-defined unit axes for mechanism-switch decomposition."""
    features = list(feature_cols)
    axes: dict[str, dict[str, float]] = {
        "matrix_high_dct_low": {
            "matrix_component": 1.0,
            "ecm_organization": 1.0,
            "fibrosis_tgfb_emt": 0.75,
            "integrin_cell_adhesion": 0.75,
            "mmp_adam_proteolysis": 0.75,
            "dct_transport_component": -1.0,
            "dct_ncc_wnk_transport": -0.75,
        },
        "circadian_dct_state": {
            "clock_core": 1.0,
            "per_cry": 0.75,
            "bmal_clock": 0.50,
            "dct_wnk_expression": 0.75,
            "aldosterone_enac_expression": 0.50,
            "Per2_expression": 0.50,
        },
        "s1p_matrix": {
            "s1p_s1pr3": 1.0,
            "s1p_axis_expression": 1.0,
            "S1pr3_expression": 0.75,
            "matrix_component": 0.50,
            "integrin_cell_adhesion": 0.50,
        },
        "preservation_cold_stress": {
            "preservation_stress_component": 1.0,
            "preservation_stress_response": 1.0,
            "rbm3_preservation_expression": 1.0,
            "Rbm3_expression": 0.75,
        },
        "immune_tlr4_macrophage": {
            "immune_context_component": 1.0,
            "tlr4_innate": 1.0,
            "macrophage_inflammation": 1.0,
        },
        "transport_rebound": {
            "dct_transport_component": 1.0,
            "dct_ncc_wnk_transport": 0.75,
            "dct_wnk_expression": 1.0,
            "aldosterone_enac_expression": 0.75,
            "Slc12a3_expression": 0.50,
            "Wnk4_expression": 0.50,
        },
    }
    mat = pd.DataFrame(0.0, index=features, columns=axes.keys())
    for axis, weights in axes.items():
        for feat, weight in weights.items():
            if feat in mat.index:
                mat.loc[feat, axis] = weight
        norm = float(np.linalg.norm(mat[axis].to_numpy(dtype=float)))
        if norm > EPS:
            mat[axis] = mat[axis] / norm
    return mat


def mechanism_switch_decomposition(
    df: pd.DataFrame,
    feature_cols: Sequence[str],
    *,
    age_scope: str = "pooled",
) -> pd.DataFrame:
    """Project arm effects onto candidate mechanism-switch axes."""
    rows = rrrm2_flt_gc_rows(df, age_scope=age_scope)
    features = [c for c in feature_cols if c in rows.columns]
    axes = mechanism_axis_definitions(features)
    out: list[dict[str, object]] = []
    arm_vectors = {
        "ISS-T": arm_effect(rows, features, "ISS-T", age_scope=age_scope),
        "LAR": arm_effect(rows, features, "LAR", age_scope=age_scope),
    }
    axis_matrix = axes.to_numpy(dtype=float)
    for arm, vec in arm_vectors.items():
        v = vec.reindex(features).fillna(0.0).to_numpy(dtype=float)
        multi = np.linalg.lstsq(axis_matrix, v, rcond=None)[0] if axis_matrix.size else np.asarray([])
        fitted = axis_matrix @ multi if axis_matrix.size else np.zeros_like(v)
        residual_norm = float(np.linalg.norm(v - fitted))
        effect_norm = float(np.linalg.norm(v))
        variance_explained = 1.0 - (residual_norm ** 2 / effect_norm ** 2) if effect_norm > EPS else np.nan
        for j, axis in enumerate(axes.columns):
            axis_vec = axes[axis]
            coef = projection_beta(pd.Series(v, index=features), axis_vec)
            out.append({
                "age_scope": age_scope,
                "arm": arm,
                "axis": axis,
                "simple_projection_coefficient": coef,
                "multi_axis_regression_coefficient": float(multi[j]) if len(multi) else np.nan,
                "effect_norm": effect_norm,
                "residual_norm_all_axes": residual_norm,
                "variance_explained_all_axes": variance_explained,
            })
    return pd.DataFrame(out)


def spearman_pair(df: pd.DataFrame, x: str, y: str) -> tuple[int, float, float]:
    """Spearman correlation with defensive missing-data handling."""
    sub = df[[x, y]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(sub) < 4 or sub[x].nunique() < 2 or sub[y].nunique() < 2:
        return int(len(sub)), np.nan, np.nan
    rho, p = stats.spearmanr(sub[x].astype(float), sub[y].astype(float))
    return int(len(sub)), float(rho), float(p)
