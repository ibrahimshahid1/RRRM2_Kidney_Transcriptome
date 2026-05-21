#!/usr/bin/env python3
"""Plot the LAR attenuation/divergence and mechanism-context dashboard."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


REPO = Path(__file__).resolve().parents[1]
DEFAULT_LAR = REPO / "data/results/run_20260519_000547_2500g/contrast_vectors/mechanism_axis/tubulointerstitial_state/lar_reversal"
DEFAULT_OUT = REPO / "latex_paper/figures_v4"

C_ISS = "#264653"
C_LAR = "#2A9D8F"
C_FLT = "#E63946"
C_GC = "#457B9D"
C_GREY = "#9AA0A6"
BG = "#FAFAFA"
RNG = np.random.default_rng(20260520)


def setup() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "font.family": "sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    sns.set_style("whitegrid")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")


def clean(ax: plt.Axes) -> None:
    ax.set_facecolor(BG)
    ax.grid(True, alpha=0.25)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def read_table(root: Path, name: str) -> pd.DataFrame:
    path = root / name
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t")


def panel_gene_scatter(ax: plt.Axes, root: Path) -> None:
    df = read_table(root, "lar_reversal_gene_scatter.tsv")
    if df.empty:
        return
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["iss_t_effect", "lar_effect"])
    ax.scatter(df["iss_t_effect"], df["lar_effect"], s=8, color=C_GREY, alpha=0.22, linewidths=0)
    panel_map = {
        "clock_core": "#7B2CBF",
        "dct_wnk_expression": "#0077B6",
        "s1p_axis_expression": "#F77F00",
        "rbm3_preservation_expression": "#D62828",
    }
    for panel, color in panel_map.items():
        sub = df[df["highlight_panels"].fillna("").str.contains(panel, regex=False)]
        if not sub.empty:
            ax.scatter(sub["iss_t_effect"], sub["lar_effect"], s=34, color=color, alpha=0.85, edgecolor="white", linewidth=0.4, label=panel.replace("_expression", ""))
    if "network_priority_rank" in df.columns:
        top = df[df["network_priority_rank"].le(25)]
        ax.scatter(top["iss_t_effect"], top["lar_effect"], s=46, facecolor="none", edgecolor="#111111", linewidth=0.8, label="top 25 network priority")
    for symbol in ["Rbm3", "Per2", "S1pr3", "Slc12a3", "Wnk4"]:
        hit = df[df["mgi_symbol"].astype(str).eq(symbol)]
        if hit.empty:
            continue
        row = hit.iloc[0]
        ax.text(row["iss_t_effect"], row["lar_effect"], f" {symbol}", fontsize=7, weight="bold")
    lim = np.nanmax(np.abs(pd.concat([df["iss_t_effect"], df["lar_effect"]]).to_numpy(dtype=float)))
    lim = min(max(lim, 0.5), 3.0)
    ax.plot([-lim, lim], [-lim, lim], color="#777777", lw=0.8, ls="--")
    ax.plot([-lim, lim], [lim, -lim], color="#222222", lw=0.9, ls=":")
    ax.axhline(0, color="#999999", lw=0.7)
    ax.axvline(0, color="#999999", lw=0.7)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("ISS-T FLT - GC residual expression")
    ax.set_ylabel("LAR FLT - GC residual expression")
    ax.set_title("Gene-level direction check")
    ax.legend(frameon=True, loc="upper right", ncol=1)
    clean(ax)


def panel_cosine_heatmap(ax: plt.Axes, root: Path) -> None:
    df = read_table(root, "lar_reversal_vector_summary.tsv")
    df = df[df["feature_set"].isin(["state_space", "mechanism_scores", "special_clock_s1p_rbm3", "wgcna_module_eigengenes", "all_gene_expression"])].copy()
    if df.empty:
        return
    order = ["state_space", "mechanism_scores", "special_clock_s1p_rbm3", "wgcna_module_eigengenes", "all_gene_expression"]
    piv = df.pivot_table(index="feature_set", columns="age_scope", values="cos_lar_iss", aggfunc="first").reindex(order)
    piv = piv[[c for c in ["pooled", "YNG", "OLD"] if c in piv.columns]]
    sns.heatmap(piv, ax=ax, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", cbar_kws={"label": "cos(LAR, ISS-T)"})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Directional geometry by feature layer")


def panel_component_effects(ax: plt.Axes, root: Path) -> None:
    df = read_table(root, "lar_reversal_pathway_interactions.tsv")
    keep = [
        "matrix_component",
        "dct_transport_component",
        "matrix_minus_dct",
        "s1p_s1pr3",
        "clock_core",
        "Rbm3_expression",
        "rbm3_preservation_expression",
    ]
    df = df[df["feature"].isin(keep)].copy()
    if df.empty:
        return
    long = df.melt(
        id_vars=["feature"],
        value_vars=["iss_t_effect", "lar_effect"],
        var_name="arm",
        value_name="effect",
    )
    long["arm"] = long["arm"].replace({"iss_t_effect": "ISS-T", "lar_effect": "LAR"})
    sns.barplot(data=long, y="feature", x="effect", hue="arm", palette={"ISS-T": C_ISS, "LAR": C_LAR}, ax=ax)
    ax.axvline(0, color="#777777", lw=0.8)
    ax.set_xlabel("FLT - GC effect")
    ax.set_ylabel("")
    ax.set_title("Component effects")
    ax.legend(frameon=True, loc="lower right")
    clean(ax)


def panel_switch_axes(ax: plt.Axes, root: Path) -> None:
    df = read_table(root, "lar_mechanism_switch_scores.tsv")
    if df.empty:
        return
    piv = df.pivot_table(index="axis", columns="arm", values="multi_axis_regression_coefficient", aggfunc="first")
    sns.heatmap(piv, ax=ax, cmap="vlag", center=0, annot=True, fmt=".2f", cbar_kws={"label": "axis coefficient"})
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Exploratory mechanism-axis projections")


def draw_points(ax: plt.Axes, x: float, vals: pd.Series, color: str) -> None:
    vals = vals.dropna().astype(float).to_numpy()
    if len(vals) == 0:
        return
    ax.scatter(x + RNG.uniform(-0.10, 0.10, len(vals)), vals, s=30, color=color, alpha=0.78, edgecolor="white", linewidth=0.5)
    ax.plot([x - 0.17, x + 0.17], [np.median(vals), np.median(vals)], color=color, lw=2)


def panel_targeted_expression(ax: plt.Axes, root: Path) -> None:
    df = read_table(root, "lar_reversal_sample_scores.tsv")
    features = [f for f in ["Rbm3_expression", "Per2_expression", "S1pr3_expression"] if f in df.columns]
    if df.empty or not features:
        return
    sub = df[df["condition"].isin(["GC", "FLT"])].copy()
    x = 0
    ticks = []
    labels = []
    for feature in features:
        for arm in ["ISS-T", "LAR"]:
            for cond, color in [("GC", C_GC), ("FLT", C_FLT)]:
                vals = sub[sub["Arm"].eq(arm) & sub["condition"].eq(cond)][feature]
                draw_points(ax, x, vals, color)
                ticks.append(x)
                labels.append(f"{feature.replace('_expression','')}\n{arm}\n{cond}")
                x += 1
            x += 0.25
        x += 0.65
    ax.axhline(0, color="#999999", lw=0.8, ls="--")
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=55, ha="right")
    ax.set_ylabel("gene z-expression")
    ax.set_title("Rbm3 / Per2 / S1pr3 by mouse")
    clean(ax)


def panel_age_stratified(ax: plt.Axes, root: Path) -> None:
    df = read_table(root, "lar_reversal_component_effects.tsv")
    df = df[df["feature"].eq("matrix_minus_dct") & df["age_scope"].isin(["YNG", "OLD", "pooled"])].copy()
    if df.empty:
        return
    sns.pointplot(data=df, x="age_scope", y="effect_flt_minus_gc", hue="arm", palette={"ISS-T": C_ISS, "LAR": C_LAR}, dodge=0.35, linestyle="none", ax=ax)
    ax.axhline(0, color="#777777", lw=0.8)
    ax.set_xlabel("")
    ax.set_ylabel("matrix_minus_dct FLT - GC")
    ax.set_title("Age-stratified state score")
    ax.legend(frameon=True, loc="best")
    clean(ax)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lar-dir", default=str(DEFAULT_LAR))
    ap.add_argument("--outdir", default=str(DEFAULT_OUT))
    ap.add_argument("--name", default="fig12_lar_reversal_dashboard")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(args.lar_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    setup()
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1.0, 1.0], hspace=0.42, wspace=0.30)
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, 0]),
        fig.add_subplot(gs[2, 1]),
    ]
    panel_gene_scatter(axes[0], root)
    panel_cosine_heatmap(axes[1], root)
    panel_component_effects(axes[2], root)
    panel_switch_axes(axes[3], root)
    panel_targeted_expression(axes[4], root)
    panel_age_stratified(axes[5], root)
    for ax, label in zip(axes, "ABCDEF"):
        panel_label(ax, label)
    fig.suptitle("LAR attenuation/divergence and mechanism-context analysis", fontsize=14, fontweight="bold", y=0.995)
    fig.savefig(outdir / f"{args.name}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(outdir / f"{args.name}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[saved] {outdir / (args.name + '.pdf')}")


if __name__ == "__main__":
    main()
