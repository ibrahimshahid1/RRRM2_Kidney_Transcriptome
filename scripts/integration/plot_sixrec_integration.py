#!/usr/bin/env python3
"""Consolidated six-recommendation integration figure.

Ties together the memo's Recommendation analyses into one cross-omic
decoupling / evidence-ladder panel:

  A  RNA pathway-vector recurrence (RRRM-2 ISS-T vs OSD-462)        [Rec 3 / 1]
  B  Protein-abundance concordance null (matrix & DCT)              [Rec 3]
  C  Phospho layer: NCC regulatory vs non-regulatory + total NCC    [Rec 3]
  D  Regulator activity: KSEA WNK/SPAK + recurrent pathways/TFs     [Rec 2]
  E  Cell-type decomposition: DCT down, interstitial up (cohorts)   [Rec 4]
  F  Evidence ladder + honest grades (Rec 1-6)                      [Rec 5/6]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RES = REPO / "data" / "results"
ANCHOR = RES / "run_20260522_osd462_anchor" / "osd462_anchor"
REG = RES / "run_20260522_regulator_activity"
CELL = RES / "run_20260522_celltype_decomposition"
OUTDIRS = [RES / "run_20260522_sixrec_integration", REPO / "latex_paper" / "figures_osd462"]

UP, DOWN, NEUT = "#c0392b", "#2471a3", "#7f8c8d"


def panel_a(ax):
    df = pd.read_csv(ANCHOR / "osd462_rna_pathway_effects.tsv", sep="\t")
    x = df["rrrm2_iss_t_pathway_effect"]; y = df["osd462_rna_pathway_effect"]
    ax.axhline(0, color="k", lw=.6); ax.axvline(0, color="k", lw=.6)
    ax.scatter(x, y, c=[UP if v > 0 else DOWN for v in y], s=42, zorder=3, edgecolor="k", lw=.4)
    for _, r in df.iterrows():
        if abs(r["osd462_rna_pathway_effect"]) > 0.1 or "dct" in r["pathway"]:
            ax.annotate(r["pathway"].replace("_", " ")[:14],
                        (r["rrrm2_iss_t_pathway_effect"], r["osd462_rna_pathway_effect"]),
                        fontsize=5.5, ha="left", va="bottom")
    ax.set_xlabel("RRRM-2 ISS-T pathway effect"); ax.set_ylabel("OSD-462 RNA pathway effect")
    ax.set_title("A  RNA context recurs\ncosine = 0.87 (gate PASS)", fontsize=9, loc="left", weight="bold")


def panel_b(ax):
    df = pd.read_csv(ANCHOR / "protein_concordance_summary.tsv", sep="\t")
    df = df[df["kind"] == "family"]
    names = [n.replace("_", " ")[:12] for n in df["gene_set"]]
    eff = df["signed_mean_effect"].to_numpy(float)
    lo = df["signed_mean_null_lo"].to_numpy(float); med = df["signed_mean_null_med"].to_numpy(float)
    y = np.arange(len(names))
    ax.axvline(0, color="k", lw=.6)
    ax.barh(y, eff, color=[UP if e > 0 else DOWN for e in eff], edgecolor="k", lw=.4)
    ax.scatter(med, y, color="k", marker="|", s=80, zorder=4, label="matched null median")
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=6.5)
    ax.set_xlabel("protein signed-mean flight effect")
    ax.set_title("B  Protein abundance:\nnull / opposite (H1,H2 falsified)", fontsize=9, loc="left", weight="bold")
    ax.legend(fontsize=5.5, loc="lower right")


def panel_c(ax):
    s = pd.read_csv(ANCHOR / "phospho_axis_summary.tsv", sep="\t")
    s["key"] = s["gene_symbol"] + " " + s["site_position"].astype(str)
    reg = [("Slc12a3", "53"), ("Slc12a3", "65"), ("Slc12a3", "68"),
           ("Stk39", "382"), ("Stk39", "383")]
    non = [("Slc12a3", "96"), ("Slc12a3", "120"), ("Slc12a3", "122"), ("Slc12a3", "124")]
    def eff(g, p):
        r = s[(s.gene_symbol == g) & (s.site_position.astype(str) == p)]
        return float(r["phospho_effect"].iloc[0]) if len(r) else np.nan
    labels, vals, cols = [], [], []
    for g, p in reg:
        labels.append(f"{g[:4]} {p}*"); vals.append(eff(g, p)); cols.append(DOWN)
    for g, p in non:
        labels.append(f"{g[:4]} {p}"); vals.append(eff(g, p)); cols.append(NEUT)
    labels.append("NCC total\nprotein"); vals.append(0.0889); cols.append("#27ae60")
    ax.axhline(0, color="k", lw=.6)
    ax.bar(range(len(vals)), vals, color=cols, edgecolor="k", lw=.4)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=5.8)
    ax.set_ylabel("flight effect (log2)")
    ax.set_title("C  Phospho layer: regulatory* down,\nnon-reg flat, total NCC flat (H3)", fontsize=9, loc="left", weight="bold")


def panel_d(ax):
    items, vals, cols = [], [], []
    # KSEA kinases
    k = pd.read_csv(REG / "osd462_kinase_activity_summary.tsv", sep="\t")
    for _, r in k.iterrows():
        items.append(f"{r['kinase']} (KSEA)"); vals.append(r["z_score"]); cols.append(DOWN)
    # recurrent injury pathways (mean flight activity) + endothelial TF module
    pr = pd.read_csv(REG / "progeny_recurrence.tsv", sep="\t").set_index("source")
    for p in ["NFkB", "TNFa", "Hypoxia", "JAK-STAT"]:
        if p in pr.index:
            items.append(f"{p} (path)"); vals.append(pr.loc[p, "mean_activity_flight"]); cols.append(UP)
    tf = pd.read_csv(REG / "collectri_recurrence.tsv", sep="\t").set_index("source")
    endo = [t for t in ["Erg", "Gata2", "Ets1"] if t in tf.index]
    if endo:
        items.append("Endothelial TFs\n(Erg/Gata2/Ets1)")
        vals.append(float(tf.loc[endo, "mean_activity_flight"].mean())); cols.append(UP)
    y = np.arange(len(items))
    ax.axvline(0, color="k", lw=.6)
    ax.barh(y, vals, color=cols, edgecolor="k", lw=.4)
    ax.set_yticks(y); ax.set_yticklabels(items, fontsize=6); ax.invert_yaxis()
    ax.set_xlabel("activity (KSEA z / ULM z)")
    ax.set_title("D  Regulators: WNK/SPAK down (anchor);\ninjury pathways up; endothelial TFs", fontsize=9, loc="left", weight="bold")


def panel_e(ax):
    c = pd.read_csv(CELL / "celltype_flight_effects_by_cohort.tsv", sep="\t", index_col=0)
    cohorts = [x for x in ["RRRM2_ISS_T_young", "OSD462_rna", "OSD513", "OSD253"] if x in c.columns]
    panels = ["dct_identity", "dct_transport", "endothelial", "stromal_fibroblast", "macrophage_immune"]
    sub = c.loc[panels, cohorts]
    x = np.arange(len(panels)); w = 0.8 / len(cohorts)
    for i, co in enumerate(cohorts):
        ax.bar(x + i * w, sub[co].to_numpy(float), w, label=co.replace("_", " "),
               edgecolor="k", lw=.3)
    ax.axhline(0, color="k", lw=.6)
    ax.set_xticks(x + 0.4 - w / 2)
    ax.set_xticklabels([p.replace("_", "\n") for p in panels], fontsize=6)
    ax.set_ylabel("mean flight stat")
    ax.legend(fontsize=5, ncol=2, loc="lower left")
    ax.set_title("E  Cell-type: DCT down, endothelial/stromal up\n(endo vs NCC phospho rho=-0.76, p<.001)", fontsize=9, loc="left", weight="bold")


def panel_f(ax):
    ax.axis("off")
    rows = [
        ("Recurrent RNA context (matrix-high/DCT-low)", "supported", UP),
        ("Protein-abundance concordance", "falsified (null)", DOWN),
        ("NCC/SPAK regulatory phospho suppression", "supported (anchor)", UP),
        ("Network-candidate translation", "falsified (null)", DOWN),
        ("Frozen RNA <-> NCC activity, group level", "supported", UP),
        ("...per-animal predictive link", "underpowered", NEUT),
        ("Non-obvious regulator (beyond injury)", "none; endo. axis nominated", NEUT),
        ("DCT-low = pure transcriptional suppression", "no; interstitial component", NEUT),
        ("Animal-matched physiology", "absent (Tier-3 future work)", NEUT),
    ]
    ax.set_title("F  Evidence ladder & honest grades", fontsize=9, loc="left", weight="bold")
    y = 0.94
    for label, grade, col in rows:
        ax.text(0.01, y, "●", color=col, fontsize=9, va="top")
        ax.text(0.06, y, label, fontsize=6.6, va="top")
        ax.text(0.99, y, grade, fontsize=6.4, va="top", ha="right", color=col, style="italic")
        y -= 0.105


def main() -> int:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    panel_a(axes[0, 0]); panel_b(axes[0, 1]); panel_c(axes[0, 2])
    panel_d(axes[1, 0]); panel_e(axes[1, 1]); panel_f(axes[1, 2])
    fig.suptitle("Cross-omic decoupling of the spaceflight kidney DCT/NCC phenotype "
                 "— six-recommendation integration", fontsize=12, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    for d in OUTDIRS:
        d.mkdir(parents=True, exist_ok=True)
        for ext in ("png", "pdf"):
            fig.savefig(d / f"fig_sixrec_integration.{ext}", dpi=180, bbox_inches="tight")
    print(f"wrote fig_sixrec_integration.png/.pdf to {[str(d) for d in OUTDIRS]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
