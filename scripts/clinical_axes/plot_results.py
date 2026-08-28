#!/usr/bin/env python3
"""Render publication-sized diagnostics for the clinical renal-axis audit."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


REPO = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"
DEFAULT_FIGURES = REPO / "figures/clinical_renal_axes_cross_mission"

AXIS_LABELS = {
    "glomerular_barrier_identity_loss": "Glomerular barrier\nidentity loss",
    "tubular_epithelial_injury_induction": "Tubular epithelial\ninjury induction",
    "fibrosis_maladaptive_remodeling": "Fibrosis / maladaptive\nremodelling",
    "distal_transport_identity_loss": "Distal transport\nidentity loss",
}

MISSION_COLORS = {
    "OSD-102": "#4C78A8",
    "OSD-163": "#F58518",
    "OSD-253": "#54A24B",
    "OSD-462": "#E45756",
    "OSD-771": "#7A5195",
}


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.dpi": 160,
            "savefig.dpi": 300,
        }
    )


def _save(fig: plt.Figure, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), bbox_inches="tight", facecolor="white")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def primary_forest(results: Path, figures: Path) -> None:
    meta = pd.read_csv(results / "primary_meta_results.tsv", sep="\t").set_index("axis")
    missions = pd.read_csv(results / "primary_mission_effects.tsv", sep="\t")
    axes = list(AXIS_LABELS)
    fig, panels = plt.subplots(1, 4, figsize=(12.2, 4.3), sharex=True, sharey=True)
    y = np.arange(5)[::-1]
    for panel, axis in zip(panels, axes, strict=True):
        sub = missions[missions["axis"] == axis].set_index("mission").loc[list(MISSION_COLORS)]
        estimates = sub["estimate"].to_numpy()
        errors = 1.96 * sub["standard_error"].to_numpy()
        for index, mission in enumerate(sub.index):
            panel.errorbar(
                estimates[index],
                y[index],
                xerr=errors[index],
                fmt="o",
                color=MISSION_COLORS[mission],
                ecolor=MISSION_COLORS[mission],
                capsize=2,
                markersize=5,
                linewidth=1.2,
            )
        pooled = meta.loc[axis]
        panel.errorbar(
            pooled["estimate"],
            -1,
            xerr=np.array(
                [
                    [pooled["estimate"] - pooled["ci_low_mkh"]],
                    [pooled["ci_high_mkh"] - pooled["estimate"]],
                ]
            ),
            fmt="D",
            color="black",
            ecolor="black",
            capsize=3,
            markersize=5,
            linewidth=1.5,
        )
        panel.axvline(0, color="#777777", linewidth=0.8, linestyle="--")
        panel.axhline(-0.5, color="#D9D9D9", linewidth=0.8)
        panel.set_title(AXIS_LABELS[axis])
        panel.set_xlim(-2.7, 2.7)
        panel.set_xticks([-2, -1, 0, 1, 2])
        panel.grid(axis="x", color="#EEEEEE", linewidth=0.6)
        panel.text(
            0.02,
            0.02,
            f"maxT p={pooled['max_t_fwer']:.3g}\nI²={pooled['i_squared']:.0f}%",
            transform=panel.transAxes,
            va="bottom",
            fontsize=8,
        )
    panels[0].set_yticks(list(y) + [-1])
    panels[0].set_yticklabels(list(MISSION_COLORS) + ["Pooled"])
    fig.supxlabel("Hedges g (flight minus matched ground; positive = declared adverse direction)")
    fig.suptitle("Terminal-flight effects across five independent mouse missions", y=1.02)
    _save(fig, figures / "figure_1_primary_axes_forest")


def compartment_forest(results: Path, figures: Path) -> None:
    table = pd.read_csv(results / "compartment_context_meta_results.tsv", sep="\t")
    selected = table[table["tier"] == "high_specificity"].copy()
    control = table[table["axis"] == "broad_structural_scaffold_control__all"].copy()
    selected = pd.concat([selected, control], ignore_index=True)
    selected["label"] = selected["report_compartment"].replace(
        {
            "DCT1_context": "DCT1",
            "DCT2_CNT_context": "DCT2/CNT",
            "mesenchymal_stromal": "Fibroblast/stromal",
            "structural_scaffold_control": "Broad structural control",
        }
    )
    selected["label"] = selected["label"].str.replace("_", " ", regex=False)
    selected = selected.sort_values("estimate")
    y = np.arange(len(selected))
    significant = selected["max_t_fwer"] <= 0.05
    colors = np.where(significant, "#7A5195", "#9E9E9E")
    fig, ax = plt.subplots(figsize=(7.6, max(4.5, 0.37 * len(selected) + 1.4)))
    lower = selected["estimate"] - selected["ci_low_mkh"]
    upper = selected["ci_high_mkh"] - selected["estimate"]
    for i, (_, row) in enumerate(selected.iterrows()):
        ax.errorbar(
            row["estimate"],
            y[i],
            xerr=np.array([[lower.iloc[i]], [upper.iloc[i]]]),
            fmt="o",
            color=colors[i],
            ecolor=colors[i],
            capsize=2,
            markersize=5,
        )
        ax.text(
            1.02,
            y[i],
            f"FWER={row['max_t_fwer']:.3g}",
            transform=ax.get_yaxis_transform(),
            va="center",
            fontsize=7.5,
        )
    ax.axvline(0, color="#777777", linewidth=0.8, linestyle="--")
    ax.set_yticks(y)
    ax.set_yticklabels(selected["label"])
    ax.set_xlabel("Hedges g (higher compartment-marker transcript abundance in flight →)")
    ax.set_title("Frozen high-specificity kidney-compartment comparison")
    ax.grid(axis="x", color="#EEEEEE", linewidth=0.6)
    ax.set_xlim(
        min(-1.7, float(selected["ci_low_mkh"].min()) - 0.1),
        max(1.7, float(selected["ci_high_mkh"].max()) + 0.1),
    )
    legend = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#7A5195", label="maxT FWER ≤ 0.05", markersize=6),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#9E9E9E", label="does not pass family", markersize=6),
    ]
    ax.legend(handles=legend, loc="lower right", frameon=False)
    ax.text(
        0.0,
        -0.12,
        "Only prespecified high-specificity sets with ≥8 observable markers are shown.",
        transform=ax.transAxes,
        fontsize=7.5,
        color="#555555",
    )
    fig.subplots_adjust(right=0.80)
    _save(fig, figures / "figure_2_compartment_specificity_forest")


def podocyte_adversarial(results: Path, figures: Path) -> None:
    null = pd.read_csv(results / "barrier_matched_panel_null.tsv.gz", sep="\t")
    summary = pd.read_csv(results / "barrier_matched_panel_summary.tsv", sep="\t").iloc[0]
    proxy = pd.read_csv(results / "barrier_proxy_adjustment_meta.tsv", sep="\t")
    proxy = proxy[proxy["variance_type"] == "HC3"].copy()
    proxy["label"] = proxy["model"].map(
        {
            "unadjusted_blocked": "Unadjusted",
            "adjusted_disjoint_podocyte_proxy": "Adjusted for disjoint\npodocyte proxy",
        }
    )
    fig, panels = plt.subplots(1, 2, figsize=(10.2, 4.2))
    panels[0].hist(null["estimate"], bins=45, color="#B8C7D9", edgecolor="white")
    panels[0].axvline(summary["target_meta_estimate"], color="#7A5195", linewidth=2)
    panels[0].text(
        summary["target_meta_estimate"],
        panels[0].get_ylim()[1] * 0.92,
        f"  barrier core\n  matched p={summary['matched_two_sided_p']:.3g}",
        color="#5E3C78",
        va="top",
        fontsize=8,
    )
    panels[0].set_xlabel("Pooled expression-increase effect for matched 6-gene sets")
    panels[0].set_ylabel("Matched sets")
    panels[0].set_title("Observability- and breadth-matched null")

    y = np.arange(len(proxy))[::-1]
    for i, (_, row) in enumerate(proxy.iterrows()):
        panels[1].errorbar(
            row["estimate"],
            y[i],
            xerr=np.array(
                [[row["estimate"] - row["ci_low_mkh"]], [row["ci_high_mkh"] - row["estimate"]]]
            ),
            fmt="o",
            color="#7A5195"
            if row["model"] == "unadjusted_blocked"
            else "#E45756",
            capsize=3,
            markersize=6,
        )
    panels[1].axvline(0, color="#777777", linewidth=0.8, linestyle="--")
    panels[1].set_yticks(y)
    panels[1].set_yticklabels(proxy["label"])
    panels[1].set_xlabel("Flight coefficient (standardized score units)")
    panels[1].set_title("Barrier core is not specific beyond\nthe broader podocyte program")
    panels[1].grid(axis="x", color="#EEEEEE", linewidth=0.6)
    fig.suptitle("Adversarial interpretation tests for the glomerular signal", y=1.02)
    fig.tight_layout()
    _save(fig, figures / "figure_3_podocyte_adversarial_tests")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--figures", type=Path, default=DEFAULT_FIGURES)
    args = parser.parse_args()
    _style()
    primary_forest(args.results, args.figures)
    compartment_forest(args.results, args.figures)
    podocyte_adversarial(args.results, args.figures)
    print(f"Wrote figures to {args.figures}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
