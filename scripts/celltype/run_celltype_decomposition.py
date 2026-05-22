#!/usr/bin/env python3
"""Cell-type marker-panel decomposition (memo Recommendation 4).

Bounds the bulk-RNA "is this a DCT program change or a composition/dilution
change?" ambiguity without new single-cell data:

1. Cross-cohort marker-panel flight effects (DCT identity vs DCT transport vs
   PT / TAL / CNT-CD / endothelial / stromal / macrophage), reusing the
   per-cohort gene-effect tables built for regulator Layer B.
2. OSD-462 animal-matched test: per-sample compartment scores, the flight
   effect of each, and whether the DCT-transport / immune / stromal scores
   covary with the measured NCC/SPAK regulatory phospho score (dilution check).
3. Scenario decision per memo Section 8.3.
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.common import id_map_lookup  # noqa: E402
from src.multiomics.celltype_panels import (  # noqa: E402
    KIDNEY_PANELS, decide_scenario, panel_flight_effect, per_sample_panel_scores)
from src.multiomics.phenotype_anchor import (  # noqa: E402
    NCC_REGULATORY_SITES, channel_to_animal, per_sample_score, rna_sample_animal)
from scripts.regulator_activity.run_phenotype_anchor import (  # noqa: E402
    load_phospho_per_sample)

OSD462 = REPO_ROOT / "data" / "external" / "osdr" / "OSD-462"
VST = OSD462 / "GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv"
SAMPLES = OSD462 / "GLDS-462_rna_seq_SampleTable_mRNA_GLbulkRNAseq.csv"
PHOS = OSD462 / "GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx"
ID_MAP = REPO_ROOT / "data" / "processed" / "resources" / "id_map.tsv"
SPEC = REPO_ROOT / "data" / "results" / "run_20260522_regulator_activity" / "rna_effects" / "rna_effects_spec.json"
OUT = REPO_ROOT / "data" / "results" / "run_20260522_celltype_decomposition"


def log(m: str) -> None:
    print(f"{datetime.now():%H:%M:%S} [celltype] {m}")


def welch_effect(scores: pd.DataFrame, cond: pd.Series, fl: str, gc: str) -> pd.DataFrame:
    rows = []
    for panel, row in scores.iterrows():
        a = row[[c for c in scores.columns if cond.get(c) == fl]].dropna().to_numpy(float)
        b = row[[c for c in scores.columns if cond.get(c) == gc]].dropna().to_numpy(float)
        if len(a) < 2 or len(b) < 2:
            rows.append({"panel": panel, "flt_minus_gc": np.nan, "t": np.nan, "p": np.nan})
            continue
        t, p = stats.ttest_ind(a, b, equal_var=False)
        rows.append({"panel": panel, "flt_minus_gc": float(a.mean() - b.mean()),
                     "t": float(t), "p": float(p), "n_flt": len(a), "n_gc": len(b)})
    return pd.DataFrame(rows)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    ens_to_sym, sym_to_ens = id_map_lookup(str(ID_MAP))

    # ---- 1. cross-cohort panel flight effects -------------------------------
    log("Cross-cohort marker-panel flight effects")
    spec = json.loads(SPEC.read_text())
    cohort_eff = {}
    for label, rel in spec.items():
        gs = pd.read_csv(REPO_ROOT / rel, sep="\t")
        pe = panel_flight_effect(gs, KIDNEY_PANELS)
        cohort_eff[label] = pe.set_index("panel")["mean_stat"]
    cross = pd.DataFrame(cohort_eff)
    cross.to_csv(OUT / "celltype_flight_effects_by_cohort.tsv", sep="\t")
    print(cross.round(2).to_string())

    # ---- 2. OSD-462 per-sample compartment scores ---------------------------
    log("OSD-462 per-sample compartment scores")
    vst = pd.read_csv(VST, index_col=0)
    vst.index = vst.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    scores = per_sample_panel_scores(vst, KIDNEY_PANELS, sym_to_ens)
    samp = pd.read_csv(SAMPLES, index_col=0)["condition"]
    cond = samp.reindex(scores.columns)
    eff = welch_effect(scores, cond, "Space.Flight", "Ground.Control")
    eff.to_csv(OUT / "osd462_compartment_flight_effects.tsv", sep="\t", index=False)
    scores.to_csv(OUT / "osd462_compartment_scores_per_sample.tsv", sep="\t")
    print("\nOSD-462 flight effect (Space.Flight - Ground.Control), per-sample mean-z panels:")
    print(eff.round(3).to_string(index=False))

    # ---- 3. animal-matched: panels vs NCC regulatory phospho ----------------
    log("Animal-matched compartment <-> NCC phospho correlations")
    phospho = load_phospho_per_sample(PHOS)
    reg_keys = [f"{g}|{s}" for g, s in NCC_REGULATORY_SITES]
    reg_per_channel = per_sample_score(phospho, reg_keys)
    ncc = {}
    for ch, v in reg_per_channel.items():
        m = channel_to_animal(ch)
        if m:
            ncc.setdefault(f"{m[0]}|{m[1]}", []).append(v)
    ncc = pd.Series({k: float(np.mean(v)) for k, v in ncc.items()})

    def by_animal(s: pd.Series) -> pd.Series:
        rows = {}
        for sample, val in s.items():
            m = rna_sample_animal(sample)
            if m:
                rows.setdefault(f"{m[0]}|{m[1]}", []).append(val)
        return pd.Series({k: float(np.mean(v)) for k, v in rows.items()})

    corr_rows = []
    flight_ground = [a for a in ncc.index if a.split("|")[0] in ("Space Flight", "Ground Control")]
    for panel in ["dct_transport", "dct_identity", "macrophage_immune",
                  "stromal_fibroblast", "endothelial"]:
        pa = by_animal(scores.loc[panel])
        common = sorted(set(pa.index) & set(ncc.index) & set(flight_ground))
        rho, p = stats.spearmanr(pa.reindex(common), ncc.reindex(common))
        corr_rows.append({"panel": panel, "n_animals": len(common),
                          "spearman_vs_ncc_phospho": float(rho), "p": float(p)})
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(OUT / "osd462_compartment_vs_ncc_phospho.tsv", sep="\t", index=False)
    print("\nAnimal-matched Spearman vs NCC regulatory phospho (FLT+GC):")
    print(corr.round(3).to_string(index=False))

    # ---- 4. scenario decision (memo 8.3) ------------------------------------
    e = eff.set_index("panel")["flt_minus_gc"]
    decision = decide_scenario(e["dct_transport"], e["dct_identity"],
                               e["stromal_fibroblast"], e["macrophage_immune"])
    # recurrence: DCT transport down in how many cohorts; identity preserved?
    dct_t = cross.loc["dct_transport"]
    dct_i = cross.loc["dct_identity"]
    decision["dct_transport_down_cohorts"] = int((dct_t < 0).sum())
    decision["dct_identity_down_cohorts"] = int((dct_i < 0).sum())
    decision["n_cohorts"] = int(cross.shape[1])
    verdict = {
        "analysis": "cell-type marker-panel decomposition (memo Rec 4)",
        "timestamp": datetime.now().isoformat(),
        "osd462_scenario": decision,
        "dilution_check": (
            "If DCT-transport-low tracked immune-high per-animal, the DCT signal "
            "could be interstitial dilution; reported in "
            "osd462_compartment_vs_ncc_phospho.tsv."),
        "panels": {k: v for k, v in KIDNEY_PANELS.items()},
    }
    (OUT / "celltype_decomposition_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    log(f"scenario = {decision['scenario']}")
    print("\nREADING:", decision["reading"])
    log(f"Wrote outputs to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
