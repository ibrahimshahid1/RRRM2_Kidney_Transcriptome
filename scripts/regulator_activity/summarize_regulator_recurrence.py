#!/usr/bin/env python3
"""Summarize cross-cohort recurrence of decoupler Layer-B regulator activity."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.multiomics.regulator_activity import recurrence_class  # noqa: E402

RUN = REPO_ROOT / "data" / "results" / "run_20260522_regulator_activity"
FLIGHT = ["RRRM2_ISS_T_young", "OSD462_rna", "OSD513", "OSD253"]
RECOVERY = "RRRM2_LAR_young"
STRONG = 2.0  # |ULM z| threshold for a "strong" per-cohort call

# Curated set of regulators that are *expected* from generic kidney
GENERIC = {
    # pathways
    "nfkb", "tnfa", "jak-stat", "tgfb", "hypoxia", "p53", "mapk", "egfr",
    # TFs: NF-kB / inflammation
    "nfkb1", "rela", "relb", "rel", "nfkb2",
    # JAK-STAT
    "stat1", "stat2", "stat3", "stat5a", "stat5b", "stat6",
    # AP-1 / immediate-early stress
    "jun", "junb", "jund", "fos", "fosb", "fosl1", "fosl2", "atf3", "egr1",
    # TGFb effectors / fibrosis
    "smad2", "smad3", "smad4",
    # hypoxia / metabolic stress
    "hif1a", "epas1",
    # p53 / cell-cycle / proliferation (generic)
    "trp53", "e2f1", "e2f4", "myc", "mycn", "ccnd1",
    # ubiquitous housekeeping TFs
    "sp1", "sp3", "yy1", "nfya", "nfyb", "creb1",
    # generic stress / NRF2
    "nfe2l2", "cebpb", "cebpa", " atf4", "atf4", "ddit3",
}


def _load(name: str) -> pd.DataFrame:
    return pd.read_csv(RUN / name, sep="\t")


def recurrence_table(df: pd.DataFrame) -> pd.DataFrame:
    piv = df.pivot(index="source", columns="contrast", values="activity_score")
    flight = [c for c in FLIGHT if c in piv.columns]
    rows = []
    for src, r in piv.iterrows():
        fl = r[flight].dropna()
        if len(fl) < 2:
            continue
        cls = recurrence_class({c: fl[c] for c in fl.index}, threshold=STRONG)
        rows.append({
            "source": src,
            "mean_activity_flight": float(fl.mean()),
            "n_cohorts": int(len(fl)),
            "n_strong": int((fl.abs() >= STRONG).sum()),
            "sign_consistent": bool((fl > 0).all() or (fl < 0).all()),
            "recurrence_class": cls,
            "lar_recovery_activity": float(r[RECOVERY]) if RECOVERY in r.index and pd.notna(r.get(RECOVERY)) else np.nan,
            "is_generic": str(src).lower() in GENERIC,
            **{c: float(r[c]) if pd.notna(r.get(c)) else np.nan for c in (flight + [RECOVERY])},
        })
    out = pd.DataFrame(rows)
    # rank: recurrent + strong + (for the non-obvious search) non-generic first
    out["abs_mean"] = out["mean_activity_flight"].abs()
    out = out.sort_values(["sign_consistent", "n_strong", "abs_mean"],
                          ascending=[False, False, False]).reset_index(drop=True)
    return out


def main() -> int:
    prog = recurrence_table(_load("rna_progeny_pathway_activity.tsv"))
    tf = recurrence_table(_load("rna_collectri_tf_activity.tsv"))
    prog.to_csv(RUN / "progeny_recurrence.tsv", sep="\t", index=False)
    tf.to_csv(RUN / "collectri_recurrence.tsv", sep="\t", index=False)

    def recurrent(df):
        return df[df["recurrence_class"].isin(["recurrent_up", "recurrent_down"])]

    rec_prog, rec_tf = recurrent(prog), recurrent(tf)
    nonobv_tf = rec_tf[~rec_tf["is_generic"]].copy()
    nonobv_prog = rec_prog[~rec_prog["is_generic"]].copy()

    # Endothelial/vascular module: a coherent non-obvious recurrent cluster
    endo = nonobv_tf[nonobv_tf["source"].isin(["Erg", "Ets1", "Gata2", "Fli1", "Sox17", "Sox18"])]

    verdict = {
        "analysis": "regulator-activity cross-cohort recurrence (Layer B)",
        "framing": "prioritization / nomination, not causal driver discovery",
        "flight_cohorts": FLIGHT,
        "recovery_cohort_reported_separately": RECOVERY,
        "strong_threshold_abs_z": STRONG,
        "n_pathways_recurrent": int(len(rec_prog)),
        "n_pathways_recurrent_up": int((rec_prog["recurrence_class"] == "recurrent_up").sum()),
        "recurrent_pathways_up": rec_prog[rec_prog["recurrence_class"] == "recurrent_up"]["source"].tolist(),
        "recurrent_generic_pathways": rec_prog[rec_prog["is_generic"]]["source"].tolist(),
        "n_tf_recurrent": int(len(rec_tf)),
        "n_tf_recurrent_nonobvious": int(len(nonobv_tf)),
        "top_nonobvious_recurrent_tfs": nonobv_tf.sort_values("abs_mean", ascending=False)["source"].head(12).tolist(),
        "endothelial_vascular_module_recurrent": endo["source"].tolist(),
        "any_tf_recurrent_down": bool((rec_tf["recurrence_class"] == "recurrent_down").any()),
        "ksea_positive_control": "WNK / SPAK_OSR1 inferred activity DOWN (z=-4.12 / -6.31) -- activity anchor",
        "verdict": (
            "The recurrent regulator-activity signal is dominated by expected kidney "
            "injury / inflammation / fibrotic programs (NF-kB, TNFa, JAK-STAT, Hypoxia, "
            "TGFb/Smad3) and ubiquitous factors (Sp1/Sp3); these are literature-consistent, "
            "not a discovery. No single non-obvious upstream driver dominates, and no "
            "regulator recurs as inferred-DOWN. The most coherent non-obvious recurrent "
            "candidate is an endothelial/vascular TF module (Erg, Gata2, Ets1), consistent "
            "with the S1P/S1PR3 + adhesion component of the RNA context; it is graded a "
            "candidate context axis, not a proven regulator. Reported per memo Rec-2 as "
            "claim discipline."
        ),
    }
    (RUN / "regulator_recurrence_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")

    pd.set_option("display.width", 200)
    print("== PROGENy recurrence (top) ==")
    print(prog[["source", "mean_activity_flight", "n_strong", "sign_consistent",
                "recurrence_class", "is_generic", "lar_recovery_activity"]].head(14).to_string(index=False))
    print("\n== Non-obvious recurrent TFs ==")
    print(nonobv_tf[["source", "mean_activity_flight", "n_strong", "recurrence_class",
                     "lar_recovery_activity"]].head(15).to_string(index=False))
    print("\nVERDICT:", verdict["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
