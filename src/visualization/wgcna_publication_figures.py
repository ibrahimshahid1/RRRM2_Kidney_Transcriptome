#!/usr/bin/env python3
"""WGCNA Publication Figures"""
from __future__ import annotations
import argparse, os, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch
from scipy import stats as sp_stats

warnings.filterwarnings("ignore")

# Color palette
C_FLT = "#E63946"
C_GC  = "#457B9D"
C_ISS = "#264653"
C_LAR = "#2A9D8F"
BG    = "#FAFAFA"


def load_data(results_dir: str):
    """Load all WGCNA outputs needed for figures."""
    rd = Path(results_dir)
    wdir = rd / "wgcna"
    fdir = wdir / "manuscript_followup"

    data = {}
    data["wdir"] = wdir
    data["fdir"] = fdir

    # Module assignments
    data["mod_assign"] = pd.read_csv(wdir / "module_assignments.csv")

    # Eigengenes
    data["MEs"] = pd.read_csv(wdir / "module_eigengenes.csv", index_col=0)

    # Trait associations
    data["trait"] = pd.read_csv(wdir / "module_trait_association.csv")

    # Preservation
    data["pres"] = pd.read_csv(wdir / "module_preservation.csv")

    # ME → color map
    me_map = {f"ME{k}": v for k, v in
              data["mod_assign"].groupby("module_num")["module_color"].first().items()}
    data["me_map"] = me_map
    data["color_to_me"] = {v: k for k, v in me_map.items()}

    # Simple-effect contrasts
    contrasts_path = fdir / "simple_effect_contrasts.csv"
    if contrasts_path.exists():
        data["contrasts"] = pd.read_csv(contrasts_path)

    # Hub genes
    hub_path = fdir / "grey60_hub_genes_kME.csv"
    if hub_path.exists():
        data["hubs"] = pd.read_csv(hub_path)

    # External validation
    ext_dir = rd / "external_validation_wgcna"
    if ext_dir.exists():
        ext_results = {}
        for study_dir in sorted(ext_dir.iterdir()):
            proj_file = study_dir / "wgcna_module_projection.tsv"
            if proj_file.exists():
                ext_results[study_dir.name] = pd.read_csv(proj_file, sep="\t")
        data["external"] = ext_results

    # Metadata for eigengene plots
    for scan_dir in sorted((rd.parent).glob("run_*"), reverse=True):
        meta_path = scan_dir / "phase1_residuals/meta_phase1.tsv.gz"
        if meta_path.exists():
            meta = pd.read_csv(meta_path, sep="\t", compression="gzip")
            scol = [c for c in meta.columns if "Sample Name" in c or c == "sample"][0]
            meta = meta.set_index(scol, drop=False)
            meta["EnvGroup"] = meta["EnvGroup"].replace({"HGC": "GC", "VGC": "VIV"})
            meta["Age"] = meta["Age"].replace({"Young": "YNG", "young": "YNG", "Yng": "YNG", "old": "OLD"})
            meta["Arm"] = meta["Arm"].replace({"ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T"})
            keep = meta["EnvGroup"].isin(["FLT", "GC"])
            data["meta"] = meta[keep]
            break

    return data


def fig_module_trait_heatmap(data: dict, out: Path):
    """Heatmap of module-trait correlations with significance stars."""
    trait = data["trait"]
    me_map = data["me_map"]
    terms = ["Flight", "AgeOld", "ArmLAR", "Flight:AgeOld", "Flight:ArmLAR"]
    modules = sorted([m for m in me_map if me_map[m] != "grey"],
                     key=lambda x: int(x[2:]))

    mat = np.zeros((len(modules), len(terms)))
    stars = [["" for _ in terms] for _ in modules]

    for i, me in enumerate(modules):
        for j, term in enumerate(terms):
            row = trait[(trait["module"] == me) & (trait["term"] == term)]
            if len(row):
                mat[i, j] = row["estimate"].values[0]
                q = row["q_value"].values[0]
                if q < 0.001: stars[i][j] = "***"
                elif q < 0.01: stars[i][j] = "**"
                elif q < 0.05: stars[i][j] = "*"

    fig, ax = plt.subplots(figsize=(7, 10))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")

    for i in range(len(modules)):
        for j in range(len(terms)):
            txt = f"{mat[i,j]:.2f}\n{stars[i][j]}"
            color = "white" if abs(mat[i,j]) > 0.35 else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color=color)

    ylabels = [f"{me} ({me_map[me]})" for me in modules]
    ax.set_yticks(range(len(modules)))
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xticks(range(len(terms)))
    ax.set_xticklabels(terms, fontsize=9, rotation=30, ha="right")
    ax.set_title("Module–Trait Association (estimate, BH q)", fontsize=12, fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.03, label="Coefficient estimate")
    fig.tight_layout()
    fig.savefig(out / "fig1_module_trait_heatmap.pdf", dpi=300)
    fig.savefig(out / "fig1_module_trait_heatmap.png", dpi=300)
    plt.close(fig)
    print(f"  fig1_module_trait_heatmap")


def fig_eigengene_dotplot(data: dict, out: Path, mod_color: str, fig_num: int):
    """Eigengene dot plot by EnvGroup × Arm × Age."""
    me_name = data["color_to_me"].get(mod_color)
    if me_name is None or "meta" not in data:
        return
    meta = data["meta"]
    MEs = data["MEs"]
    common = [s for s in meta.index if s in MEs.index]
    meta = meta.loc[common]
    me_vals = MEs.loc[common, me_name].values

    groups = [
        ("FLT / ISS-T / YNG", "FLT", "ISS-T", "YNG"),
        ("GC / ISS-T / YNG",  "GC",  "ISS-T", "YNG"),
        ("FLT / ISS-T / OLD", "FLT", "ISS-T", "OLD"),
        ("GC / ISS-T / OLD",  "GC",  "ISS-T", "OLD"),
        ("FLT / LAR / YNG",   "FLT", "LAR",   "YNG"),
        ("GC / LAR / YNG",    "GC",  "LAR",   "YNG"),
        ("FLT / LAR / OLD",   "FLT", "LAR",   "OLD"),
        ("GC / LAR / OLD",    "GC",  "LAR",   "OLD"),
    ]

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    rng = np.random.default_rng(42)

    for i, (label, env, arm, age) in enumerate(groups):
        mask = (meta["EnvGroup"] == env) & (meta["Arm"] == arm) & (meta["Age"] == age)
        vals = me_vals[mask.values]
        if len(vals) == 0:
            continue
        jitter = rng.uniform(-0.15, 0.15, len(vals))
        c = C_FLT if env == "FLT" else C_GC
        ax.scatter([i]*len(vals) + jitter, vals, c=c, s=50, alpha=0.75,
                   edgecolors="white", linewidths=0.6, zorder=3)
        ax.plot([i-0.25, i+0.25], [np.median(vals)]*2, c=c, lw=2.5, zorder=4)

    for sep in [1.5, 3.5, 5.5]:
        ax.axvline(sep, color="#DDD", lw=0.8, ls=":")
    ax.axhline(0, color="#AAA", lw=0.5, ls="--")

    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], rotation=40, ha="right", fontsize=8)

    n_genes = len(data["mod_assign"][data["mod_assign"]["module_color"] == mod_color])
    ax.set_ylabel(f"{me_name} eigengene", fontsize=11)
    ax.set_title(f"Module eigengene: {mod_color} ({n_genes} genes)", fontsize=13, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(handles=[Patch(fc=C_FLT, label="FLT"), Patch(fc=C_GC, label="GC")],
              loc="upper right", framealpha=0.9)
    fig.tight_layout()
    fig.savefig(out / f"fig{fig_num}_eigengene_{mod_color}.pdf", dpi=300)
    fig.savefig(out / f"fig{fig_num}_eigengene_{mod_color}.png", dpi=300)
    plt.close(fig)
    print(f"  fig{fig_num}_eigengene_{mod_color}")


def fig_preservation(data: dict, out: Path):
    """Module preservation bar chart."""
    pres = data["pres"]
    pres = pres[pres["module_color"] != "grey"].sort_values("Zsummary", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 7))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    colors_bar = []
    for _, r in pres.iterrows():
        if r["Zsummary"] > 10:
            colors_bar.append("#2A9D8F")
        elif r["Zsummary"] > 2:
            colors_bar.append("#E9C46A")
        else:
            colors_bar.append("#E76F51")

    bars = ax.barh(range(len(pres)), pres["Zsummary"].values, color=colors_bar,
                   edgecolor="white", linewidth=0.5)
    ax.set_yticks(range(len(pres)))
    ax.set_yticklabels([f"{r['module_color']} (n={r['moduleSize']})"
                        for _, r in pres.iterrows()], fontsize=8)

    ax.axvline(10, color="#264653", ls="--", lw=1, label="Strong (>10)")
    ax.axvline(2, color="#E76F51", ls="--", lw=1, label="Moderate (>2)")
    ax.set_xlabel("Zsummary (preservation)", fontsize=11)
    ax.set_title("Module Preservation: GC → FLT\n(combined WGCNA modules)", fontsize=12, fontweight="bold")
    ax.legend(loc="lower right", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / "fig4_module_preservation.pdf", dpi=300)
    fig.savefig(out / "fig4_module_preservation.png", dpi=300)
    plt.close(fig)
    print(f"  fig4_module_preservation")


def fig_hub_lollipop(data: dict, out: Path):
    """Hub gene kME lollipop for grey60."""
    if "hubs" not in data:
        return
    hubs = data["hubs"].head(20).copy()
    hubs = hubs.iloc[::-1]  # reverse for bottom-up plot

    fig, ax = plt.subplots(figsize=(7, 6))
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    colors = ["#E63946" if "ecm" in str(a).lower() or "fibr" in str(a).lower()
              else "#457B9D" for a in hubs["pathway_annotation"]]

    ax.hlines(range(len(hubs)), 0, hubs["kME"].values, color=colors, lw=2, alpha=0.7)
    ax.scatter(hubs["kME"].values, range(len(hubs)), c=colors, s=60, zorder=3,
               edgecolors="white", linewidths=0.5)

    labels = [f"{r['symbol']}" if pd.notna(r['symbol']) and r['symbol'] != ''
              else r['gene'][:15] for _, r in hubs.iterrows()]
    ax.set_yticks(range(len(hubs)))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Module membership (kME)", fontsize=11)
    ax.set_title("Grey60 top 20 hub genes by |kME|", fontsize=12, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(handles=[Patch(fc="#E63946", label="ECM / fibrosis"),
                       Patch(fc="#457B9D", label="Other")],
              loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig5_grey60_hub_genes.pdf", dpi=300)
    fig.savefig(out / "fig5_grey60_hub_genes.png", dpi=300)
    plt.close(fig)
    print(f"  fig5_grey60_hub_genes")


def fig_contrast_forest(data: dict, out: Path):
    """Forest plot of simple-effect contrasts for significant modules."""
    if "contrasts" not in data:
        return

    cdf = data["contrasts"]
    sig_colors = cdf[cdf["q"] < 0.10]["color"].unique()
    if len(sig_colors) == 0:
        return

    strata_order = ["ISS-T_YNG", "ISS-T_OLD", "LAR_YNG", "LAR_OLD"]
    strata_labels = ["ISS-T Young", "ISS-T Old", "LAR Young", "LAR Old"]
    strata_colors = [C_ISS, "#1D3557", C_LAR, "#52B788"]

    fig, axes = plt.subplots(1, len(sig_colors), figsize=(4*len(sig_colors), 5),
                             sharey=True, squeeze=False)
    fig.patch.set_facecolor(BG)

    for col_idx, color in enumerate(sorted(sig_colors)):
        ax = axes[0, col_idx]
        ax.set_facecolor(BG)
        mod_data = cdf[cdf["color"] == color]

        for i, (stratum, label, sc) in enumerate(zip(strata_order, strata_labels, strata_colors)):
            row = mod_data[mod_data["stratum"] == stratum]
            if len(row):
                r = row.iloc[0]
                ax.errorbar(r["estimate"], i, xerr=1.96*r["SE"],
                           fmt="o", color=sc, capsize=4, markersize=7, lw=1.5)
                sig = "***" if r["q"] < 0.001 else "**" if r["q"] < 0.01 else "*" if r["q"] < 0.05 else ""
                ax.text(r["estimate"] + 1.96*r["SE"] + 0.02, i, sig,
                       fontsize=10, va="center", fontweight="bold", color=sc)

        ax.axvline(0, color="#AAA", ls="--", lw=0.8)
        ax.set_title(f"{color}", fontsize=11, fontweight="bold")
        ax.set_xlabel("FLT − GC (eigengene)", fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0, 0].set_yticks(range(len(strata_labels)))
    axes[0, 0].set_yticklabels(strata_labels, fontsize=9)
    fig.suptitle("Simple-Effect Contrasts: FLT vs GC per Age × Arm", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out / "fig6_contrast_forest.pdf", dpi=300)
    fig.savefig(out / "fig6_contrast_forest.png", dpi=300)
    plt.close(fig)
    print(f"  fig6_contrast_forest")


def fig_external_validation(data: dict, out: Path):
    """Cross-study module projection comparison."""
    if "external" not in data or len(data["external"]) == 0:
        return

    key_modules = ["grey60", "red", "salmon", "green", "pink", "royalblue"]
    studies = sorted(data["external"].keys())
    study_colors = ["#264653", "#2A9D8F", "#E9C46A", "#E76F51"]

    fig, axes = plt.subplots(1, len(key_modules), figsize=(3*len(key_modules), 5),
                             sharey=True, squeeze=False)
    fig.patch.set_facecolor(BG)

    for col_idx, mod in enumerate(key_modules):
        ax = axes[0, col_idx]
        ax.set_facecolor(BG)
        for i, study in enumerate(studies):
            df = data["external"][study]
            row = df[df["module"] == mod]
            if len(row):
                r = row.iloc[0]
                c = study_colors[i % len(study_colors)]
                marker = "D" if r["perm_q"] < 0.05 else "o"
                ax.scatter(r["mean_diff"], i, c=c, s=80 if r["perm_q"] < 0.05 else 50,
                          marker=marker, edgecolors="white", linewidths=0.5, zorder=3)
                if r["perm_q"] < 0.05:
                    ax.text(r["mean_diff"] + 0.05, i, f"q={r['perm_q']:.3f}",
                           fontsize=7, va="center", color=c)

        ax.axvline(0, color="#AAA", ls="--", lw=0.8)
        ax.set_title(mod, fontsize=10, fontweight="bold")
        ax.set_xlabel("FLT−GC score", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0, 0].set_yticks(range(len(studies)))
    axes[0, 0].set_yticklabels(studies, fontsize=8)
    fig.suptitle("External Validation: Module Score FLT − GC", fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out / "fig7_external_validation.pdf", dpi=300)
    fig.savefig(out / "fig7_external_validation.png", dpi=300)
    plt.close(fig)
    print(f"  fig7_external_validation")


def main():
    ap = argparse.ArgumentParser(description="WGCNA publication figures")
    ap.add_argument("--results_dir", required=True)
    args = ap.parse_args()

    results_dir = args.results_dir
    out = Path(results_dir) / "figures" / "publication"
    out.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("WGCNA Publication Figures")
    print("="*60)

    data = load_data(results_dir)

    fig_module_trait_heatmap(data, out)
    fig_eigengene_dotplot(data, out, "grey60", 2)
    fig_eigengene_dotplot(data, out, "red", 3)
    fig_preservation(data, out)
    fig_hub_lollipop(data, out)
    fig_contrast_forest(data, out)
    fig_external_validation(data, out)

    # Also generate supplementary eigengene plots for other significant modules
    supp = out / "supplementary"
    supp.mkdir(exist_ok=True)
    for i, mod in enumerate(["salmon", "green", "pink", "royalblue"], start=1):
        fig_eigengene_dotplot(data, supp, mod, f"S{i}")

    print(f"\nFigures saved to: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
