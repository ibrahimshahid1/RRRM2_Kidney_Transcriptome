"""Module convergence: Fisher's-exact overlap between"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO = Path(__file__).resolve().parents[1]
LIONESS_RUN = REPO / "data/results/run_20260517_213205_2500g"
WGCNA_RUN = REPO / "data/results/run_20260505_remediated_2500g/wgcna"
OUTDIR = LIONESS_RUN / "module_convergence"
OUTDIR.mkdir(parents=True, exist_ok=True)


def bh_qvals(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    out = np.empty_like(q)
    out[order] = q
    return out


def fisher_per_module(
    candidate_genes: set[str],
    module_assignments: pd.DataFrame,
    label: str,
    analysis_universe: set[str] | None = None,
    universe_source: str = "WGCNA",
) -> pd.DataFrame:
    wgcna_universe = set(module_assignments["gene"])
    if analysis_universe is None:
        universe = wgcna_universe
    else:
        universe = wgcna_universe & analysis_universe
    cand_in_universe = candidate_genes & universe
    N = len(universe)
    K = len(cand_in_universe)

    rows = []
    for color, sub in module_assignments.groupby("module_color"):
        module_genes = set(sub["gene"]) & universe
        M = len(module_genes)
        x = len(cand_in_universe & module_genes)
        table = [[x, K - x], [M - x, N - K - (M - x)]]
        odds, p = fisher_exact(table, alternative="greater")
        expected = K * M / N if N else float("nan")
        rows.append({
            "module": color,
            "module_size": M,
            "n_candidates_total": K,
            "n_candidates_raw": len(candidate_genes),
            "n_candidates_in_module": x,
            "expected_under_independence": expected,
            "fold_enrichment": (x / expected) if expected > 0 else np.nan,
            "odds_ratio_greater": odds,
            "p_fisher_greater": p,
        })

    df = pd.DataFrame(rows)
    df["q_BH"] = bh_qvals(df["p_fisher_greater"].values)
    df = df.sort_values("p_fisher_greater").reset_index(drop=True)
    df.insert(0, "candidate_set", label)
    df.insert(1, "candidates_dropped_off_panel", len(candidate_genes - universe))
    df.insert(2, "universe_size", N)
    df.insert(3, "universe_source", universe_source)
    return df


def load_lioness_candidates() -> set[str]:
    perm = pd.read_csv(
        LIONESS_RUN / "phase6_uncertainty" / "ISS_T_YNG_FLT_minus_GC_perm_pvals.tsv",
        sep="\t",
    )
    sig = perm[perm["q_BH_edge_sum"] < 0.05]
    return set(sig["gene"].tolist())


def load_silent_shifters(contrast: str) -> set[str]:
    path = (
        LIONESS_RUN
        / "phase5_silent_shifters_strict"
        / f"{contrast}_silent_shifters.tsv"
    )
    df = pd.read_csv(path, sep="\t")
    df = df[df["strict_silent_shifter"] == True]
    return set(df["gene"].tolist())


def load_phase6_universe(contrast: str) -> set[str]:
    path = LIONESS_RUN / "phase6_uncertainty" / f"{contrast}_perm_pvals.tsv"
    df = pd.read_csv(path, sep="\t")
    return set(df["gene"].astype(str))


def load_module_assignments() -> pd.DataFrame:
    df = pd.read_csv(WGCNA_RUN / "module_assignments.csv")
    df = df.rename(columns={"module_color": "module_color"})
    return df[["gene", "module_color", "symbol"]].copy()


def reframe_eigengene_contrasts() -> pd.DataFrame:
    """Per-module ISS-T vs LAR flight-effect comparison from existing module_trait fits.

    Coding: Age reference = Young, Arm reference = ISS-T, FlightStatus reference = Control.
      ISS-T Young flight effect = Flight
      LAR   Young flight effect = Flight + Flight:ArmLAR
      ISS-T Old   flight effect = Flight + Flight:AgeOld
      LAR   Old   flight effect = Flight + Flight:AgeOld + Flight:ArmLAR + Flight:AgeOld:ArmLAR
    """
    mt = pd.read_csv(WGCNA_RUN / "module_trait_association.csv")
    mt["term"] = mt["term"].str.strip('"')
    mt["module"] = mt["module"].str.strip('"')
    integ = pd.read_csv(WGCNA_RUN / "integrated_module_table.csv")

    me_to_color = dict(zip(integ["ME"], integ["color"]))

    pivot_est = mt.pivot(index="module", columns="term", values="estimate")
    pivot_p = mt.pivot(index="module", columns="term", values="p_value")
    pivot_q = mt.pivot(index="module", columns="term", values="q_value")

    rows = []
    for me in pivot_est.index:
        e = pivot_est.loc[me]
        p = pivot_p.loc[me]
        q = pivot_q.loc[me]
        iss_t_young = e.get("Flight", np.nan)
        lar_young = iss_t_young + e.get("Flight:ArmLAR", np.nan)
        iss_t_old = iss_t_young + e.get("Flight:AgeOld", np.nan)
        lar_old = (
            iss_t_young
            + e.get("Flight:AgeOld", np.nan)
            + e.get("Flight:ArmLAR", np.nan)
            + e.get("Flight:AgeOld:ArmLAR", np.nan)
        )
        rows.append({
            "ME": me,
            "module": me_to_color.get(me, ""),
            "ISS_T_Young_flight_effect": iss_t_young,
            "LAR_Young_flight_effect": lar_young,
            "ISS_T_Old_flight_effect": iss_t_old,
            "LAR_Old_flight_effect": lar_old,
            "p_Flight_main_eq_ISS_T_Young": p.get("Flight", np.nan),
            "q_Flight_main_eq_ISS_T_Young": q.get("Flight", np.nan),
            "p_Flight_x_ArmLAR_interaction": p.get("Flight:ArmLAR", np.nan),
            "q_Flight_x_ArmLAR_interaction": q.get("Flight:ArmLAR", np.nan),
        })
    df = pd.DataFrame(rows)
    df["ISS_T_arm_specific_pattern"] = (
        (df["q_Flight_main_eq_ISS_T_Young"] < 0.05)
        & (df["q_Flight_x_ArmLAR_interaction"] < 0.05)
        & (np.sign(df["ISS_T_Young_flight_effect"]) != np.sign(df["LAR_Young_flight_effect"]))
    )
    return df.sort_values("q_Flight_main_eq_ISS_T_Young").reset_index(drop=True)


def main() -> None:
    module_assignments = load_module_assignments()
    print(
        f"WGCNA panel: {len(module_assignments)} genes, "
        f"{module_assignments['module_color'].nunique()} modules"
    )

    # Analysis (1)
    lioness_genes = load_lioness_candidates()
    lioness_universe = load_phase6_universe("ISS_T_YNG_FLT_minus_GC")
    print(f"\n[1] LIONESS ISS-T-Young q<0.05 candidates: {len(lioness_genes)}")
    df1 = fisher_per_module(
        candidate_genes=lioness_genes,
        module_assignments=module_assignments,
        label="ISS_T_YNG_LIONESS_q<0.05",
        analysis_universe=lioness_universe,
        universe_source="WGCNA intersect ISS_T_YNG phase6 tested genes",
    )
    out1 = OUTDIR / "fisher_lioness_ISS_T_YNG_vs_modules.tsv"
    df1.to_csv(out1, sep="\t", index=False)
    print(f"  -> {out1}")
    print(df1.head(5).to_string(index=False))

    # Analysis (2)
    silent_dfs = []
    for contrast in [
        "ISS_T_YNG_FLT_minus_GC",
        "ISS_T_OLD_FLT_minus_GC",
        "LAR_YNG_FLT_minus_GC",
        "LAR_OLD_FLT_minus_GC",
    ]:
        ss = load_silent_shifters(contrast)
        ss_universe = load_phase6_universe(contrast)
        print(f"\n[2] Silent shifters in {contrast}: {len(ss)}")
        df2 = fisher_per_module(
            candidate_genes=ss,
            module_assignments=module_assignments,
            label=contrast,
            analysis_universe=ss_universe,
            universe_source=f"WGCNA intersect {contrast} phase6 tested genes",
        )
        silent_dfs.append(df2)
    df2_all = pd.concat(silent_dfs, ignore_index=True)
    out2 = OUTDIR / "fisher_silent_shifters_vs_modules.tsv"
    df2_all.to_csv(out2, sep="\t", index=False)
    print(f"  -> {out2}")
    print(
        df2_all[df2_all["q_BH"] < 0.10]
        .sort_values(["candidate_set", "q_BH"])
        .head(20)
        .to_string(index=False)
    )

    # Analysis (3)
    df3 = reframe_eigengene_contrasts()
    out3 = OUTDIR / "eigengene_stratum_specific_reframe.tsv"
    df3.to_csv(out3, sep="\t", index=False)
    print(f"\n[3] Eigengene stratum-specific reframe -> {out3}")
    print(
        df3[
            [
                "ME",
                "module",
                "ISS_T_Young_flight_effect",
                "LAR_Young_flight_effect",
                "ISS_T_Old_flight_effect",
                "LAR_Old_flight_effect",
                "q_Flight_main_eq_ISS_T_Young",
                "q_Flight_x_ArmLAR_interaction",
                "ISS_T_arm_specific_pattern",
            ]
        ]
        .head(10)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
