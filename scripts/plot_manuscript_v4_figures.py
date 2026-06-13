#!/usr/bin/env python3
"""Generate manuscript v4 figures from mechanism-axis outputs."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


REPO = Path(__file__).resolve().parents[1]
RUN25 = REPO / "data/results/run_20260519_000547_2500g/contrast_vectors"
RUN50 = REPO / "data/results/run_20260519_5000g_sensitivity_3seed/contrast_vectors"
RUN_CROSS = REPO / "data/results/run_20260519_cosine_perm_null/contrast_vectors"
MECH25 = RUN25 / "mechanism_axis"
MECH50 = RUN50 / "mechanism_axis"
TSTATE = MECH25 / "tubulointerstitial_state"
OUT = REPO / "latex_paper/figures_v4"

C_FLT = "#E63946"
C_GC = "#457B9D"
C_VIV = "#2A9D8F"
C_ISS = "#264653"
C_LAR = "#2A9D8F"
C_YELLOW = "#E9C46A"
C_GREY = "#8D99AE"
BG = "#FAFAFA"

COLORS = {
    "GC": C_GC,
    "FLT": C_FLT,
    "VIV": C_VIV,
    "C57BL/6J": C_FLT,
    "C3H/HeJ": C_GC,
}

RNG = np.random.default_rng(42)


def _setup() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "font.family": "sans-serif",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    sns.set_style("whitegrid")


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight", dpi=300)
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [saved] {name}")


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def clean_axis(ax: plt.Axes) -> None:
    ax.set_facecolor(BG)
    ax.grid(True, axis="y", alpha=0.28)
    ax.grid(False, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw_points(
    ax: plt.Axes,
    x: float,
    vals: pd.Series,
    color: str,
    *,
    label: str | None = None,
    median: bool = True,
    size: int = 46,
) -> None:
    vals = vals.dropna().astype(float).to_numpy()
    if len(vals) == 0:
        return
    xs = x + RNG.uniform(-0.15, 0.15, len(vals))
    ax.scatter(
        xs,
        vals,
        s=size,
        c=color,
        alpha=0.78,
        edgecolors="white",
        linewidths=0.65,
        zorder=3,
        label=label,
    )
    if median:
        med = float(np.median(vals))
        ax.plot([x - 0.24, x + 0.24], [med, med], color=color, lw=2.6, zorder=4)


def load_scores() -> pd.DataFrame:
    return pd.read_csv(MECH25 / "mechanism_score_matrix.tsv", sep="\t")


def load_state_scores() -> pd.DataFrame:
    return pd.read_csv(TSTATE / "state_space_scores.tsv", sep="\t")


def load_state_effects() -> pd.DataFrame:
    return pd.read_csv(TSTATE / "state_space_group_effects.tsv", sep="\t")


def load_state_overlay() -> pd.DataFrame:
    return pd.read_csv(TSTATE / "rrrm2_deconvolution_overlay.tsv", sep="\t")


def load_state_correlations() -> pd.DataFrame:
    return pd.read_csv(TSTATE / "state_space_component_correlations.tsv", sep="\t")


def state_space_base(ax: plt.Axes, *, ylabel: bool = True) -> None:
    ax.axhline(0, color="#999999", linewidth=0.8, linestyle="--")
    ax.axvline(0, color="#999999", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Matrix component\nhigh = ECM/adhesion/proteolysis")
    if ylabel:
        ax.set_ylabel("DCT/NCC-WNK transport component\nnative direction; high = stronger transport")
    else:
        ax.set_ylabel("")
    ax.annotate(
        "matrix-high /\nDCT-low",
        xy=(0.90, 0.12),
        xytext=(0.58, 0.35),
        xycoords="axes fraction",
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "lw": 1.5, "color": "#333333"},
        ha="center",
        va="center",
        fontsize=8,
        color="#333333",
    )
    clean_axis(ax)


def draw_state_arrow(
    ax: plt.Axes,
    sub: pd.DataFrame,
    *,
    start_condition: str,
    end_condition: str,
    color: str,
    label: str,
) -> None:
    start = sub[sub["condition"].eq(start_condition)][["matrix_component", "dct_transport_component"]].median()
    end = sub[sub["condition"].eq(end_condition)][["matrix_component", "dct_transport_component"]].median()
    if start.isna().any() or end.isna().any():
        return
    ax.annotate(
        "",
        xy=(end["matrix_component"], end["dct_transport_component"]),
        xytext=(start["matrix_component"], start["dct_transport_component"]),
        arrowprops={"arrowstyle": "->", "lw": 2.2, "color": color, "shrinkA": 4, "shrinkB": 4},
        zorder=5,
    )
    ax.text(
        end["matrix_component"],
        end["dct_transport_component"],
        f" {label}",
        color=color,
        fontsize=8,
        fontweight="bold",
        va="center",
    )


def plot_rrrm2_remodeling(ax: plt.Axes, *, compact: bool = False) -> None:
    df = load_scores()
    df = df[df["study"].eq("RRRM-2") & df["condition"].isin(["GC", "FLT"])].copy()
    groups = [
        ("ISS-T", "YNG", "GC"),
        ("ISS-T", "YNG", "FLT"),
        ("ISS-T", "OLD", "GC"),
        ("ISS-T", "OLD", "FLT"),
        ("LAR", "YNG", "GC"),
        ("LAR", "YNG", "FLT"),
        ("LAR", "OLD", "GC"),
        ("LAR", "OLD", "FLT"),
    ]
    xticks, labels = [], []
    for i, (arm, age, cond) in enumerate(groups):
        vals = df[df["Arm"].eq(arm) & df["Age"].eq(age) & df["condition"].eq(cond)]["remodeling_score"]
        draw_points(ax, i, vals, COLORS[cond], label=cond if i < 2 else None, size=40 if compact else 46)
        xticks.append(i)
        labels.append(f"{age}\n{cond}" if compact else f"{arm}\n{age}\n{cond}")
    for sep in [3.5]:
        ax.axvline(sep, color="#D0D0D0", lw=0.9, ls=":")
    ax.axhline(0, color="#999999", lw=0.8, ls="--")
    ax.set_xticks(xticks)
    ax.set_xticklabels(labels, rotation=0 if compact else 25, ha="center" if compact else "right")
    ax.set_ylabel("RNA remodeling score")
    ax.set_title("RRRM-2 terminal-flight remodeling score")
    ax.legend(frameon=True, loc="upper left")
    clean_axis(ax)


def plot_osd513_separation(ax: plt.Axes) -> None:
    df = load_scores()
    sub = df[df["study"].eq("OSD-513") & df["condition"].isin(["GC", "FLT", "VIV"])].copy()
    order = ["GC", "FLT", "VIV"]
    for x, cond in enumerate(order):
        vals = sub[sub["condition"].eq(cond)]["remodeling_score"]
        draw_points(ax, x, vals, COLORS[cond], label=cond)
    ax.axhline(0, color="#999999", lw=0.8, ls="--")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order)
    ax.set_ylabel("RNA remodeling score")
    ax.set_title("OSD-513 reproduces flight separation")
    clean_axis(ax)


def cosine(a: pd.Series, b: pd.Series) -> float:
    aa = a.to_numpy(dtype=float)
    bb = b.to_numpy(dtype=float)
    return float(np.dot(aa, bb) / (np.linalg.norm(aa) * np.linalg.norm(bb)))


def leave_one_pathway_out_table() -> pd.DataFrame:
    axis = pd.read_csv(MECH25 / "recurrent_axis_weights.tsv", sep="\t")
    r = axis.set_index("pathway")["rrrm2_iss_effect"].astype(float)
    o = axis.set_index("pathway")["osd513_effect"].astype(float)
    pathways = axis["pathway"].tolist()
    full = cosine(r, o)
    rows = [{"pathway": "Full panel", "cosine": full, "delta": 0.0}]
    for p in pathways:
        keep = [x for x in pathways if x != p]
        c = cosine(r.loc[keep], o.loc[keep])
        rows.append({"pathway": f"Drop {p}", "cosine": c, "delta": c - full})
    return pd.DataFrame(rows).sort_values("cosine")


def plot_leave_one_pathway_out(ax: plt.Axes) -> None:
    out = leave_one_pathway_out_table()
    colors = [C_ISS if x == "Full panel" else C_GC for x in out["pathway"]]
    labels = out["pathway"].str.replace("_", " ", regex=False)
    ax.barh(labels, out["cosine"], color=colors, alpha=0.88)
    full = float(out.loc[out["pathway"].eq("Full panel"), "cosine"].iloc[0])
    ax.axvline(full, color=C_FLT, linestyle="--", linewidth=1.3, label=f"Full cosine = {full:.2f}")
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_xlabel("Terminal-flight vs OSD-513 cosine")
    ax.set_title("Alignment is not driven by one pathway")
    ax.legend(frameon=True, loc="lower right")
    clean_axis(ax)


def plot_cross_alignment(ax: plt.Axes) -> None:
    cross_path = RUN_CROSS / "cross_osdr_recurrence/cross_osdr_alignment_summary.tsv"
    if not cross_path.exists():
        cross_path = RUN25 / "cross_osdr_recurrence/cross_osdr_alignment_summary.tsv"
    cross = pd.read_csv(cross_path, sep="\t")
    arm_labels = {"ISS-T": "terminal flight", "LAR": "live return"}
    scenario_labels = {
        "powered_hgc": "",
        "original_gc_blue_light": "blue-light GC",
        "rerun_gc_white_light": "white-light GC",
    }
    cross["label"] = cross.apply(
        lambda row: " / ".join(
            part
            for part in [
                str(row["study"]),
                scenario_labels.get(str(row["scenario"]), str(row["scenario"]).replace("_", " ")),
                arm_labels.get(str(row["arm"]), str(row["arm"])),
            ]
            if part
        ),
        axis=1,
    )
    cross = cross.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(cross))
    colors = np.where(
        cross["claimable_recurrence_pass"].astype(bool),
        C_ISS,
        np.where(cross["point_estimate"] < 0, C_LAR, C_GC),
    )
    ax.hlines(y, cross["ci_low"], cross["ci_high"], color=colors, linewidth=2.6, alpha=0.9)
    ax.scatter(cross["point_estimate"], y, color=colors, s=55, edgecolor="white", linewidth=0.6, zorder=3)
    if "permutation_p_directional" in cross.columns:
        for yi, row in zip(y, cross.to_dict("records")):
            if row.get("study") == "OSD-513":
                ax.text(
                    1.02,
                    yi,
                    f"perm. p={row['permutation_p_directional']:.3f}",
                    va="center",
                    ha="right",
                    fontsize=7.2,
                    color="#333333",
                )
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(cross["label"], fontsize=8)
    ax.set_xlabel("Pathway-vector cosine")
    ax.set_title("External RNA-vector alignment")
    ax.set_xlim(-0.9, 1.08)
    clean_axis(ax)


def fig_main_result_multipanel() -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.3))
    fig.patch.set_facecolor(BG)
    plot_rrrm2_remodeling(axes[0, 0], compact=True)
    plot_osd513_separation(axes[0, 1])
    plot_cross_alignment(axes[1, 0])
    plot_leave_one_pathway_out(axes[1, 1])
    for label, ax in zip(["A", "B", "C", "D"], axes.ravel()):
        panel_label(ax, label)
    fig.suptitle(
        "Terminal-flight RNA remodeling recurs across mouse kidney cohorts",
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    save(fig, "fig1_main_result_multipanel")


def fig_rrrm2_remodeling_points() -> None:
    fig, ax = plt.subplots(figsize=(10.2, 5.1))
    fig.patch.set_facecolor(BG)
    plot_rrrm2_remodeling(ax)
    panel_label(ax, "A")
    save(fig, "fig2_rrrm2_remodeling_points")


def fig_pathway_scores_by_flight() -> None:
    df = load_scores()
    panels = [
        ("RRRM-2", "ISS-T", "RRRM-2 ISS-T"),
        ("RRRM-2", "LAR", "RRRM-2 LAR"),
        ("OSD-513", None, "OSD-513"),
    ]
    scores = [
        ("ecm_organization", "ECM organization"),
        ("dct_ncc_wnk_transport", "DCT/NCC-WNK transport"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.1), sharex="col")
    fig.patch.set_facecolor(BG)
    for col, (study, arm, title) in enumerate(panels):
        sub = df[df["study"].eq(study) & df["condition"].isin(["GC", "FLT"])].copy()
        if arm is not None:
            sub = sub[sub["Arm"].eq(arm)]
        for row, (score, ylabel) in enumerate(scores):
            ax = axes[row, col]
            for x, cond in enumerate(["GC", "FLT"]):
                vals = sub[sub["condition"].eq(cond)][score]
                draw_points(ax, x, vals, COLORS[cond])
            ax.axhline(0, color="#999999", linewidth=0.8, linestyle="--")
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["GC", "FLT"])
            if row == 0:
                ax.set_title(title, fontweight="bold")
            if col == 0:
                ax.set_ylabel(ylabel)
            clean_axis(ax)
    for label, ax in zip(["A", "B", "C", "D", "E", "F"], axes.ravel()):
        panel_label(ax, label)
    fig.suptitle("Axis-defining pathway scores by flight group", fontsize=14, fontweight="bold", y=1.02)
    save(fig, "fig3_pathway_scores_by_flight")


def fig_alignment_sensitivity() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.7))
    fig.patch.set_facecolor(BG)
    plot_cross_alignment(axes[0])
    plot_leave_one_pathway_out(axes[1])
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle("Cross-OSDR alignment and pathway-level sensitivity", fontsize=14, fontweight="bold", y=1.02)
    save(fig, "fig4_cross_alignment_leave_one_out")


def fig_osd253_strain_points() -> None:
    df = load_scores()
    sub = df[df["study"].eq("OSD-253")].copy()
    scenarios = [
        ("Original GC", ["Ground Control", "Space Flight"]),
        ("White-LED rerun GC", ["Ground Control Rerun", "Space Flight"]),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8), sharey=True)
    fig.patch.set_facecolor(BG)
    for ax, (title, conditions) in zip(axes, scenarios):
        ss = sub[sub["condition"].isin(conditions)]
        xticks, labels = [], []
        x = 0
        for strain in ["C57BL/6J", "C3H/HeJ"]:
            for cond in conditions:
                vals = ss[ss["strain"].eq(strain) & ss["condition"].eq(cond)]["remodeling_score"]
                draw_points(ax, x, vals, COLORS[strain])
                xticks.append(x)
                labels.append(f"{strain}\n{'GC' if 'Ground' in cond else 'FLT'}")
                x += 1
            x += 0.45
        ax.axhline(0, color="#999999", linewidth=0.8, linestyle="--")
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(title, fontweight="bold")
        ax.set_ylabel("Remodeling score" if title == "Original GC" else "")
        clean_axis(ax)
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle("OSD-253 remodeling score by strain and control scenario", fontsize=14, fontweight="bold", y=1.02)
    save(fig, "fig5_osd253_strain_points")


def fig_mechanism_heatmap() -> None:
    assoc = pd.read_csv(MECH25 / "mechanism_score_associations.tsv", sep="\t")
    sub = assoc[assoc["comparison"].eq("score_vs_remodeling_spearman")].copy()
    order = [
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
    ]
    cols = [("RRRM-2", "primary"), ("OSD-513", "powered_hgc"), ("OSD-253", "all_samples")]
    mat = pd.DataFrame(index=order, columns=[f"{a}\n{b}" for a, b in cols], dtype=float)
    qmat = mat.copy()
    for study, scenario in cols:
        label = f"{study}\n{scenario}"
        ss = sub[sub["study"].eq(study) & sub["scenario"].eq(scenario)].set_index("mechanism")
        mat[label] = ss["estimate"].reindex(order)
        qmat[label] = ss["q_bh"].reindex(order)

    fig, ax = plt.subplots(figsize=(6.2, 5.8))
    fig.patch.set_facecolor(BG)
    labels = [x.replace("_", " ") for x in mat.index]
    sns.heatmap(
        mat.astype(float),
        ax=ax,
        cmap="RdBu_r",
        vmin=-1,
        vmax=1,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Spearman rho"},
        yticklabels=labels,
    )
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            q = qmat.iloc[i, j]
            if pd.notna(val):
                star = "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""
                ax.text(j + 0.5, i + 0.5, f"{val:.2f}{star}", ha="center", va="center", fontsize=8)
    ax.set_title("Mechanism scores associated with remodeling axis", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    save(fig, "fig6_mechanism_score_heatmap")


def fig_gene_priority_sensitivity() -> None:
    rows = []
    for label, root in [("2.5k", MECH25), ("5k", MECH50)]:
        top = pd.read_csv(root / "gene_priority_enrichment_topk.tsv", sep="\t")
        if "universe_size" in top.columns:
            keep = set()
            for n in top["universe_size"].dropna().astype(int).unique():
                keep.update(max(1, int(round(n * f))) for f in (0.01, 0.02, 0.05, 0.10))
            top = top[top["top_k"].astype(int).isin(keep)].copy()
        top["run"] = label
        rows.append(top)
    df = pd.concat(rows, ignore_index=True)
    targets = [
        "tubular_transport_broad",
        "preservation_stress_response",
        "oxidative_stress_nrf2",
        "tlr4_innate",
        "ecm_organization",
        "fibrosis_tgfb_emt",
    ]
    sub = df[df["gene_set"].isin(targets) & df["score_col"].eq("composite_score")].copy()
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    fig.patch.set_facecolor(BG)
    line_palette = sns.color_palette("husl", len(targets))
    for ax, run in zip(axes, ["2.5k", "5k"]):
        ss = sub[sub["run"].eq(run)]
        for color, gs in zip(line_palette, targets):
            gg = ss[ss["gene_set"].eq(gs)].sort_values("top_k")
            if gg.empty:
                continue
            ax.plot(
                gg["top_k"],
                -np.log10(gg["q_bh"].clip(lower=1e-12)),
                marker="o",
                linewidth=1.8,
                color=color,
                label=gs.replace("_", " "),
            )
        ax.axhline(-np.log10(0.10), color="#777777", linestyle="--", linewidth=0.9)
        ax.set_title(f"{run} network universe", fontweight="bold")
        ax.set_xlabel("Top-k threshold")
        ax.set_ylabel("-log10 BH q")
        ax.set_ylim(0, 1.5)
        clean_axis(ax)
    axes[1].legend(frameon=True, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    panel_label(axes[0], "A")
    panel_label(axes[1], "B")
    fig.suptitle("Gene-priority enrichment under proportional thresholds", fontsize=14, fontweight="bold", y=1.02)
    save(fig, "fig7_gene_priority_sensitivity")


def fig_tubulointerstitial_state_space() -> None:
    df = load_state_scores()
    fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.7), sharex=False, sharey=False)
    fig.patch.set_facecolor(BG)

    ax = axes[0]
    sub = df[df["study"].eq("RRRM-2") & df["condition"].isin(["GC", "FLT"])].copy()
    for arm, marker, arrow_color in [("ISS-T", "o", C_ISS), ("LAR", "s", C_LAR)]:
        aa = sub[sub["Arm"].eq(arm)]
        for cond in ["GC", "FLT"]:
            ss = aa[aa["condition"].eq(cond)]
            ax.scatter(
                ss["matrix_component"],
                ss["dct_transport_component"],
                s=52,
                c=COLORS[cond],
                marker=marker,
                alpha=0.78,
                edgecolors="white",
                linewidths=0.65,
                label=f"{arm} {cond}",
            )
        draw_state_arrow(ax, aa, start_condition="GC", end_condition="FLT", color=arrow_color, label=arm)
    state_space_base(ax)
    ax.set_title("RRRM-2: arm-specific state shift", fontweight="bold")
    ax.legend(frameon=True, fontsize=7, loc="upper right", ncol=2)

    ax = axes[1]
    sub = df[df["study"].eq("OSD-513") & df["condition"].isin(["GC", "FLT", "VIV"])].copy()
    for cond, marker in [("GC", "o"), ("FLT", "o"), ("VIV", "^")]:
        ss = sub[sub["condition"].eq(cond)]
        ax.scatter(
            ss["matrix_component"],
            ss["dct_transport_component"],
            s=58,
            c=COLORS[cond],
            marker=marker,
            alpha=0.80,
            edgecolors="white",
            linewidths=0.65,
            label=cond,
        )
    draw_state_arrow(ax, sub, start_condition="GC", end_condition="FLT", color=C_ISS, label="FLT")
    state_space_base(ax, ylabel=False)
    ax.set_title("OSD-513: same state direction", fontweight="bold")
    ax.legend(frameon=True, fontsize=8, loc="upper right")

    ax = axes[2]
    sub = df[
        df["study"].eq("OSD-253")
        & df["condition"].isin(["Ground Control Rerun", "Space Flight"])
        & df["strain"].isin(["C57BL/6J", "C3H/HeJ"])
    ].copy()
    markers = {"Ground Control Rerun": "o", "Space Flight": "X"}
    labels = {"Ground Control Rerun": "rerun GC", "Space Flight": "FLT"}
    for strain, color in [("C57BL/6J", C_FLT), ("C3H/HeJ", C_GC)]:
        aa = sub[sub["strain"].eq(strain)]
        for cond in ["Ground Control Rerun", "Space Flight"]:
            ss = aa[aa["condition"].eq(cond)]
            ax.scatter(
                ss["matrix_component"],
                ss["dct_transport_component"],
                s=56,
                c=color,
                marker=markers[cond],
                alpha=0.74,
                edgecolors="white",
                linewidths=0.65,
                label=f"{strain} {labels[cond]}",
            )
        draw_state_arrow(ax, aa, start_condition="Ground Control Rerun", end_condition="Space Flight", color=color, label=strain.split("/")[0])
    state_space_base(ax, ylabel=False)
    ax.set_title("OSD-253: white-LED rerun control", fontweight="bold")
    ax.legend(frameon=True, fontsize=7, loc="upper right", ncol=1)

    for label, ax in zip(["A", "B", "C"], axes):
        panel_label(ax, label)
    fig.suptitle("Matrix-high / DCT-low tubulointerstitial state space", fontsize=14, fontweight="bold", y=1.02)
    save(fig, "fig8_tubulointerstitial_state_space")


def _corr_text(corr: pd.DataFrame, component: str, cell_estimate: str) -> str:
    row = corr[corr["component"].eq(component) & corr["cell_estimate"].eq(cell_estimate)]
    if row.empty:
        return "rho = NA"
    row = row.iloc[0]
    if pd.isna(row["spearman_rho"]):
        return "rho = NA"
    return f"rho={row['spearman_rho']:.2f}, q={row['q_bh']:.2f}"


def fig_deconvolution_overlay() -> None:
    overlay = load_state_overlay()
    corr = load_state_correlations()
    pairs = [
        ("matrix_minus_dct", "prop_DCT", "State scalar vs DCT fraction"),
        ("matrix_minus_dct", "prop_Fibroblast", "State scalar vs fibroblast fraction"),
        ("dct_transport_component", "prop_DCT", "DCT component vs DCT fraction"),
        ("matrix_component", "prop_Fibroblast", "Matrix component vs fibroblast fraction"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 7.2))
    fig.patch.set_facecolor(BG)
    for ax, (component, cell, title) in zip(axes.ravel(), pairs):
        for cond in ["BSL", "VIV", "GC", "FLT"]:
            ss = overlay[overlay["condition"].eq(cond)]
            if ss.empty:
                continue
            color = {"BSL": C_GREY, "VIV": C_VIV, "GC": C_GC, "FLT": C_FLT}.get(cond, "#555555")
            ax.scatter(
                ss[cell],
                ss[component],
                s=42,
                c=color,
                alpha=0.72,
                edgecolors="white",
                linewidths=0.55,
                label=cond,
            )
        clean_axis(ax)
        ax.set_title(f"{title}\n{_corr_text(corr, component, cell)}", fontweight="bold", fontsize=10)
        ax.set_xlabel(cell.replace("prop_", "estimated ").replace("_", " "))
        ylabel = {
            "matrix_minus_dct": "matrix_minus_dct",
            "dct_transport_component": "DCT transport component",
            "matrix_component": "matrix component",
        }[component]
        ax.set_ylabel(ylabel)
    axes[0, 0].legend(frameon=True, fontsize=8, loc="best")
    for label, ax in zip(["A", "B", "C", "D"], axes.ravel()):
        panel_label(ax, label)
    fig.suptitle("RRRM-2 deconvolution overlay for the state score", fontsize=14, fontweight="bold", y=1.02)
    save(fig, "fig9_rrrm2_deconvolution_overlay")


def fig_targeted_gene_heatmap() -> None:
    panel = pd.read_csv(TSTATE / "targeted_renal_gene_panel.tsv", sep="\t")
    panel = panel[panel["status"].eq("scored") & panel["condition"].isin(["GC", "FLT"])].copy()
    panel["group"] = panel["Arm"] + " " + panel["Age"] + " " + panel["condition"]
    group_order = [
        "ISS-T YNG GC",
        "ISS-T YNG FLT",
        "ISS-T OLD GC",
        "ISS-T OLD FLT",
        "LAR YNG GC",
        "LAR YNG FLT",
        "LAR OLD GC",
        "LAR OLD FLT",
    ]
    panel["gene_key"] = panel["panel_class"] + "|" + panel["mgi_symbol"]
    panel["expr_z"] = panel.groupby("gene_key")["expression"].transform(
        lambda s: (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) > 0 else 0.0
    )
    mat = panel.pivot_table(index=["panel_class", "mgi_symbol"], columns="group", values="expr_z", aggfunc="mean")
    mat = mat.reindex(columns=group_order)
    class_order = [
        "ECM_adhesion_proteolysis",
        "TLR_innate_macrophage",
        "DCT_NCC_WNK_transport",
        "Tubular_metabolic_priority",
    ]
    row_order = []
    for cls in class_order:
        rows = [idx for idx in mat.index if idx[0] == cls]
        row_order.extend(sorted(rows, key=lambda x: x[1].casefold()))
    mat = mat.loc[row_order]

    fig, ax = plt.subplots(figsize=(8.8, 9.8))
    fig.patch.set_facecolor(BG)
    sns.heatmap(
        mat.astype(float),
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        linewidths=0.15,
        linecolor="white",
        cbar_kws={"label": "within-gene z"},
        yticklabels=[idx[1] for idx in mat.index],
        xticklabels=[x.replace(" ", "\n") for x in mat.columns],
    )
    offset = 0
    for cls in class_order:
        n = sum(1 for idx in mat.index if idx[0] == cls)
        if n == 0:
            continue
        ax.hlines(offset, *ax.get_xlim(), colors="#333333", linewidth=0.8)
        offset += n
    ax.hlines(offset, *ax.get_xlim(), colors="#333333", linewidth=0.8)
    ax.set_title("Targeted renal biology panel in RRRM-2", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("")
    save(fig, "fig10_targeted_renal_gene_heatmap")


def fig_state_component_effects() -> None:
    effects = load_state_effects()
    component_order = [
        "matrix_component",
        "dct_transport_component",
        "immune_context_component",
        "preservation_stress_component",
        "matrix_minus_dct",
    ]
    component_labels = {
        "matrix_component": "Matrix",
        "dct_transport_component": "DCT transport\n(native)",
        "immune_context_component": "Immune\ncontext",
        "preservation_stress_component": "Preservation\nstress",
        "matrix_minus_dct": "Matrix - DCT",
    }
    contrasts = [
        ("ISS-T_age_pooled_FLT_minus_GC", "RRRM-2\nISS-T"),
        ("LAR_age_pooled_FLT_minus_GC", "RRRM-2\nLAR"),
        ("OSD-513_FLT_minus_GC", "OSD-513"),
        ("OSD-253_rerun_gc_white_light_C57BL/6J_FLT_minus_GC", "OSD-253\nC57 rerun"),
        ("OSD-253_rerun_gc_white_light_C3H/HeJ_FLT_minus_GC", "OSD-253\nC3H rerun"),
    ]
    mat = pd.DataFrame(index=[label for _, label in contrasts], columns=component_order, dtype=float)
    annot = mat.copy().astype(object)
    for contrast, label in contrasts:
        ss = effects[effects["contrast"].eq(contrast)].set_index("component")
        for comp in component_order:
            if comp not in ss.index:
                continue
            val = float(ss.loc[comp, "effect_flt_minus_gc"])
            q = ss.loc[comp, "q_bh_within_component_effects"]
            star = "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""
            mat.loc[label, comp] = val
            annot.loc[label, comp] = f"{val:.2f}{star}"

    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    fig.patch.set_facecolor(BG)
    sns.heatmap(
        mat.astype(float),
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-3,
        vmax=3,
        annot=annot,
        fmt="",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "FLT - GC effect"},
        xticklabels=[component_labels[c] for c in component_order],
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Component effects defining the matrix-high / DCT-low state", fontweight="bold")
    save(fig, "fig11_state_component_effects")


def main() -> None:
    _setup()
    fig_main_result_multipanel()
    fig_rrrm2_remodeling_points()
    fig_pathway_scores_by_flight()
    fig_alignment_sensitivity()
    fig_osd253_strain_points()
    fig_mechanism_heatmap()
    fig_gene_priority_sensitivity()
    fig_tubulointerstitial_state_space()
    fig_deconvolution_overlay()
    fig_targeted_gene_heatmap()
    fig_state_component_effects()
    print(f"Wrote figures to {OUT}")


if __name__ == "__main__":
    main()
