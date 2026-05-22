#!/usr/bin/env python3
"""Layer 5 - Multi-omics anchor figures.

Produces two figures from the anchor outputs:
  1. ``fig_osd462_multiomics_dashboard`` (2x2): RNA->protein scatter; DCT-axis
     protein-abundance bars; WNK-SPAK/OSR1-NCC phosphosite effects with CIs;
     network-candidate translation enrichment vs matched null.
  2. ``fig_osd462_rna_recurrence`` (1x2): pathway-vector concordance (RNA gate)
     and the cross-layer DCT/NCC summary (RNA / protein abundance / phospho).

Usage::

    python scripts/osd462/05_plot_dashboard.py --run RUN_NAME
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
try:
    from adjustText import adjust_text  # noqa: E402
except Exception:  # pragma: no cover - optional plotting polish
    adjust_text = None

from _common import REPO_ROOT, anchor_dir, default_run  # noqa: E402

C_MATRIX = "#E76F51"
C_DCT = "#264653"
C_TLR = "#E9C46A"
C_S1P = "#2A9D8F"
C_GREY = "#B8BDC4"
C_FLT = "#E63946"
C_GC = "#457B9D"
C_DOWN = "#1D3557"
C_SIG = "#E63946"
SET_COLORS = {"ecm_organization": C_MATRIX, "dct_ncc_wnk_transport": C_DCT,
              "tlr4_innate": C_TLR, "s1p_s1pr3": C_S1P,
              "tubular_transport_broad": "#8D99AE"}
FIG_DIRS = [REPO_ROOT / "latex_paper" / "figures_osd462"]


def setup():
    plt.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "font.size": 9,
        "axes.titlesize": 10, "axes.labelsize": 9, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "legend.fontsize": 7.5, "font.family": "sans-serif",
        "axes.spines.top": False, "axes.spines.right": False,
    })


def _save(fig, name, out_dir):
    dirs = [out_dir] + FIG_DIRS
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        fig.savefig(d / f"{name}.pdf", bbox_inches="tight")
        fig.savefig(d / f"{name}.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"[fig] wrote {name}.pdf/.png to {[str(d) for d in dirs]}")


def dashboard(out_dir):
    master = pd.read_csv(out_dir / "osd462_flight_effects.tsv", sep="\t")
    tgt = pd.read_csv(out_dir / "protein_concordance_target_genes.tsv", sep="\t")
    summ = pd.read_csv(out_dir / "protein_concordance_summary.tsv", sep="\t")
    phos = pd.read_csv(out_dir / "phospho_axis_summary.tsv", sep="\t")
    net = pd.read_csv(out_dir / "network_translation.tsv", sep="\t")
    bg = pd.read_csv(out_dir / "protein_concordance_background.tsv", sep="\t").iloc[0]

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))

    # ── Panel A: RNA -> protein scatter ──────────────────────────────────────
    ax = axes[0, 0]
    core = master.dropna(subset=["rrrm2_iss_t_rna_effect", "protein_flight_effect"])
    ax.scatter(core["rrrm2_iss_t_rna_effect"], core["protein_flight_effect"],
               s=5, c=C_GREY, alpha=0.25, edgecolors="none", rasterized=True, label="all genes")
    for s, color in SET_COLORS.items():
        sub = tgt[tgt["gene_set"] == s]
        ax.scatter(sub["rrrm2_iss_t_rna_effect"], sub["protein_flight_effect"],
                   s=34, c=color, edgecolors="k", linewidths=0.4, label=s.replace("_", " "))
    ax.axhline(0, color="k", lw=0.6, ls=":"); ax.axvline(0, color="k", lw=0.6, ls=":")
    ax.set_xlabel("RRRM-2 ISS-T RNA flight effect (FLT−GC, log2)")
    ax.set_ylabel("OSD-462 protein flight effect (FL−GC, log2)")
    ax.set_title("A  Transcript→protein concordance", loc="left", fontweight="bold")
    ax.text(0.03, 0.97, f"genome-wide Spearman = {bg['gw_spearman_rrrm2_rna_vs_protein']:.3f}\n"
                        f"Pearson = {bg['gw_pearson_rrrm2_rna_vs_protein']:.3f}",
            transform=ax.transAxes, va="top", fontsize=7.5,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9))
    ax.legend(loc="lower right", framealpha=0.9, markerscale=1.0)
    ax.set_ylim(-0.45, 0.45)

    # ── Panel B: DCT/NCC/WNK protein-abundance bars ──────────────────────────
    ax = axes[0, 1]
    dct_genes = ["Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Calb1",
                 "Nedd4l", "Cul3", "Kcnj10", "Kcnj16"]
    d = master[master["gene_symbol"].isin(dct_genes)].set_index("gene_symbol").reindex(dct_genes).dropna(subset=["protein_flight_effect"])
    y = np.arange(len(d))
    colors = [C_FLT if v > 0 else C_GC for v in d["protein_flight_effect"]]
    ax.barh(y, d["protein_flight_effect"], color=colors, edgecolor="k", linewidth=0.4)
    ax.set_yticks(y); ax.set_yticklabels(d.index, fontstyle="italic")
    ax.axvline(0, color="k", lw=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("protein flight effect (FL−GC, log2)")
    ax.set_title("B  DCT / NCC / WNK protein abundance", loc="left", fontweight="bold")
    ax.text(0.97, 0.04, "abundance ≈ flat\n(NCC slightly up)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7.5, color="0.3")

    # ── Panel C: phosphosite effects with CI ─────────────────────────────────
    ax = axes[1, 0]
    show = phos[phos["gene_symbol"].isin(["Slc12a3", "Stk39", "Oxsr1"])].copy()
    show["label"] = show["gene_symbol"] + " " + show["site_position"].astype(str)
    show = show.sort_values(["gene_symbol", "phospho_effect"])
    y = np.arange(len(show))
    sig = show["phospho_p_value"] < 0.05
    colors = [C_SIG if s else "0.6" for s in sig]
    ax.errorbar(show["phospho_effect"], y,
                xerr=[show["phospho_effect"] - show["phospho_ci_low"],
                      show["phospho_ci_high"] - show["phospho_effect"]],
                fmt="none", ecolor="0.5", elinewidth=0.8, capsize=2, zorder=1)
    ax.scatter(show["phospho_effect"], y, c=colors, s=26, edgecolors="k",
               linewidths=0.4, zorder=2)
    ax.set_yticks(y); ax.set_yticklabels(show["label"], fontsize=6.5)
    ax.axvline(0, color="k", lw=0.8)
    ax.invert_yaxis()
    ax.set_xlabel("phosphosite flight effect (FL−GC, log2)")
    ax.set_title("C  WNK-SPAK/OSR1-NCC phosphosites", loc="left", fontweight="bold")
    ax.scatter([], [], c=C_SIG, s=26, edgecolors="k", linewidths=0.4, label="p < 0.05")
    ax.scatter([], [], c="0.6", s=26, edgecolors="k", linewidths=0.4, label="n.s.")
    ax.legend(loc="lower left", framealpha=0.9)

    # ── Panel D: network translation enrichment ──────────────────────────────
    ax = axes[1, 1]
    k = net["top_k"].astype(str)
    x = np.arange(len(net)); w = 0.18
    ax.bar(x - 1.5 * w, net["protein_mean_abs_effect"], w, color=C_DCT, label="protein obs")
    ax.bar(x - 0.5 * w, net["protein_null_median"], w, color=C_DCT, alpha=0.4, label="protein null")
    ax.bar(x + 0.5 * w, net["phospho_mean_max_abs_effect"], w, color=C_S1P, label="phospho obs")
    ax.bar(x + 1.5 * w, net["phospho_null_median"], w, color=C_S1P, alpha=0.4, label="phospho null")
    for xi, (_, r) in zip(x, net.iterrows()):
        ax.text(xi - w, max(r["protein_mean_abs_effect"], r["protein_null_median"]) + 0.005,
                f"p={r['protein_enrichment_p']:.2f}", ha="center", fontsize=6.2, color=C_DCT)
        ax.text(xi + w, max(r["phospho_mean_max_abs_effect"], r["phospho_null_median"]) + 0.005,
                f"p={r['phospho_enrichment_p']:.2f}", ha="center", fontsize=6.2, color=C_S1P)
    ax.set_xticks(x); ax.set_xticklabels("top-" + k)
    ax.set_ylabel("mean |effect| (log2)")
    ax.set_title("D  Network-candidate translation", loc="left", fontweight="bold")
    ax.legend(loc="upper right", ncol=2, framealpha=0.9)

    fig.suptitle("OSD-462 / RR-10 multi-omics anchor — protein & phosphoprotein test of the "
                 "RRRM-2 matrix-high / DCT-low remodeling signal", fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, "fig_osd462_multiomics_dashboard", out_dir)


def recurrence_fig(out_dir):
    pv = pd.read_csv(out_dir / "osd462_rna_pathway_effects.tsv", sep="\t")
    summ = pd.read_csv(out_dir / "osd462_rna_recurrence.tsv", sep="\t").iloc[0]
    phos = pd.read_csv(out_dir / "phospho_axis_summary.tsv", sep="\t")
    master = pd.read_csv(out_dir / "osd462_flight_effects.tsv", sep="\t")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    # Panel A: pathway-vector concordance (RNA gate)
    ax = axes[0]
    lim = max(pv["osd462_rna_pathway_effect"].abs().max(),
              pv["rrrm2_iss_t_pathway_effect"].abs().max()) * 1.15
    ax.plot([-lim, lim], [-lim, lim], ls="--", color="0.6", lw=0.8)
    ax.axhline(0, color="k", lw=0.5, ls=":"); ax.axvline(0, color="k", lw=0.5, ls=":")
    texts = []
    for _, r in pv.iterrows():
        agree = r["sign_agree"]
        ax.scatter(r["rrrm2_iss_t_pathway_effect"], r["osd462_rna_pathway_effect"],
                   s=46, c=(C_DCT if agree else C_SIG), edgecolors="k", linewidths=0.4, zorder=3)
        texts.append(ax.annotate(r["pathway"].replace("_", " "),
                                 (r["rrrm2_iss_t_pathway_effect"], r["osd462_rna_pathway_effect"]),
                                 fontsize=6.1, xytext=(4, 3), textcoords="offset points"))
    if adjust_text is not None and texts:
        adjust_text(texts, ax=ax)
    ax.set_xlabel("RRRM-2 ISS-T RNA pathway effect"); ax.set_ylabel("OSD-462 RNA pathway effect (SF−GC)")
    ax.set_title("A  RNA recurrence gate (pathway vector)", loc="left", fontweight="bold")
    ax.text(0.03, 0.97, f"cosine = {summ['point_cosine']:.2f}\n95% CI [{summ['ci_low']:.2f}, "
                        f"{summ['ci_high']:.2f}]\nrecurrence = {'PASS' if summ['recurrence_pass'] else 'FAIL'}",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(boxstyle="round", fc="#EAF4EF", ec="0.6"))

    # Panel B: cross-layer DCT/NCC summary
    ax = axes[1]
    ncc_row = master.loc[master["gene_symbol"] == "Slc12a3"].iloc[0]
    ncc_rna = float(ncc_row["osd462_rna_effect"])
    ncc_prot = float(ncc_row["protein_flight_effect"])
    ncc_sites = phos[phos["gene_symbol"] == "Slc12a3"]
    reg = ncc_sites[ncc_sites["site_position"].astype(str).isin(["53", "65", "68", "89", "65;68", "58;65"])]
    ncc_phos = float(reg["phospho_effect"].mean())
    spak = phos[(phos["gene_symbol"] == "Stk39") & (phos["site_position"].astype(str).isin(["366", "382", "383", "382;383"]))]
    spak_phos = float(spak["phospho_effect"].mean())
    labels = ["NCC\nRNA", "NCC\nprotein", "NCC reg.\nphospho", "SPAK\nphospho"]
    vals = [ncc_rna, ncc_prot, ncc_phos, spak_phos]
    colors = [C_GC, "0.6", C_DOWN, C_DOWN]
    bars = ax.bar(np.arange(4), vals, color=colors, edgecolor="k", linewidth=0.5)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(np.arange(4)); ax.set_xticklabels(labels)
    ax.set_ylabel("flight effect (FL−GC, log2)")
    ax.set_title("B  DCT/NCC axis across layers", loc="left", fontweight="bold")
    ax.set_ylim(min(vals) - 0.14, max(vals) + 0.14)
    for xi, v in enumerate(vals):
        ax.text(xi, v + (0.03 if v >= 0 else -0.03), f"{v:+.2f}", ha="center",
                va="bottom" if v >= 0 else "top", fontsize=8)
    ax.text(0.5, 0.96, "RNA down; abundance slightly up; activating phosphorylation suppressed",
            transform=ax.transAxes, ha="center", va="top", fontsize=7, color="0.3")

    fig.tight_layout()
    _save(fig, "fig_osd462_rna_recurrence", out_dir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    args = ap.parse_args()
    out_dir = anchor_dir(args.run) / "figures"
    src_dir = anchor_dir(args.run)
    setup()
    # read inputs from the anchor dir, write figures to anchor/figures + latex dir
    global anchor_inputs
    dashboard(src_dir)
    recurrence_fig(src_dir)
    # move generated figs that dashboard/recurrence saved into src_dir into figures/
    print("[fig] done")


if __name__ == "__main__":
    main()
