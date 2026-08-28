#!/usr/bin/env python3
"""Compile the locked Grey60 decision object and decision-dashboard figure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve(path: str) -> Path:
    candidate = Path(path)
    return candidate if candidate.is_absolute() else REPO_ROOT / candidate


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def plot_gate_panel(
    ax: plt.Axes,
    statuses: list[tuple[str, str, bool]],
) -> None:
    labels = [f"{gate}  {label}" for gate, label, _ in statuses][::-1]
    passed = [state for _, _, state in statuses][::-1]
    colors = ["#238636" if state else "#cf222e" for state in passed]
    y = np.arange(len(labels))
    ax.barh(y, np.ones(len(y)), color=colors, height=0.68)
    ax.set_yticks(y, labels)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    for yi, state in zip(y, passed):
        ax.text(
            0.96,
            yi,
            "PASS" if state else "FAIL",
            ha="right",
            va="center",
            color="white",
            fontsize=9,
            fontweight="bold",
        )
    ax.set_title("A  Locked go/no gates", loc="left", fontweight="bold")
    ax.text(
        0,
        -1.0,
        "A–E were required for GO; F only controls distinctness wording.",
        fontsize=8.5,
        color="#57606a",
    )
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_flight_blind_panel(
    ax: plt.Axes,
    grid: pd.DataFrame,
    *,
    jaccard_threshold: float,
    p_threshold: float,
) -> None:
    y = -np.log10(grid["projected_permutation_p"].clip(lower=1e-12))
    sizes = 20 + 7 * grid["overlap"].to_numpy(dtype=float)
    colors = np.where(grid["primary_variant"], "#8250df", "#54aeff")
    ax.scatter(
        grid["jaccard"],
        y,
        s=sizes,
        c=colors,
        alpha=0.78,
        edgecolor="white",
        linewidth=0.6,
    )
    primary = grid.loc[grid["primary_variant"]].iloc[0]
    primary_y = -np.log10(max(float(primary["projected_permutation_p"]), 1e-12))
    ax.annotate(
        (
            f"primary: {int(primary['overlap'])} genes\n"
            f"Jaccard={primary['jaccard']:.3f}, p={primary['projected_permutation_p']:.3f}"
        ),
        (primary["jaccard"], primary_y),
        xytext=(18, 18),
        textcoords="offset points",
        fontsize=8,
        arrowprops={"arrowstyle": "-", "color": "#57606a"},
    )
    ax.axvline(jaccard_threshold, color="#cf222e", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(p_threshold), color="#cf222e", linestyle="--", linewidth=1)
    ax.set_xlim(left=0, right=max(jaccard_threshold * 1.08, grid["jaccard"].max() * 1.25))
    ax.set_ylim(bottom=0, top=max(-np.log10(p_threshold) * 1.12, y.max() * 1.25))
    ax.set_xlabel("Jaccard overlap with frozen Grey60")
    ax.set_ylabel(r"$-\log_{10}$ projected blocked-permutation p")
    ax.set_title(
        "B  Flight-blind reconstruction: 0/27 pass",
        loc="left",
        fontweight="bold",
    )
    ax.grid(alpha=0.18)


def plot_external_forest(
    ax: plt.Axes,
    cohort: pd.DataFrame,
    terminal_meta: pd.Series,
) -> None:
    study_order = ["OSD-102", "OSD-163", "OSD-253", "OSD-462"]
    rows = [
        cohort.loc[
            (cohort["study"] == study) & (cohort["context"] == "terminal")
        ].iloc[0]
        for study in study_order
    ]
    labels = study_order + ["Terminal meta", "OSD-513 (live return)"]
    estimates = [float(row["estimate_g"]) for row in rows]
    lows = [float(row["ci_low"]) for row in rows]
    highs = [float(row["ci_high"]) for row in rows]
    estimates.append(float(terminal_meta["estimate"]))
    lows.append(float(terminal_meta["ci_low_hartung_knapp"]))
    highs.append(float(terminal_meta["ci_high_hartung_knapp"]))
    row513 = cohort.loc[
        (cohort["study"] == "OSD-513") & (cohort["context"] == "live_return")
    ].iloc[0]
    estimates.append(float(row513["estimate_g"]))
    lows.append(float(row513["ci_low"]))
    highs.append(float(row513["ci_high"]))

    y = np.arange(len(labels))[::-1]
    colors = ["#0969da"] * 4 + ["#8250df", "#6e7781"]
    markers = ["s"] * 4 + ["D", "o"]
    for yi, estimate, low, high, color, marker in zip(
        y, estimates, lows, highs, colors, markers
    ):
        ax.errorbar(
            estimate,
            yi,
            xerr=np.array([[estimate - low], [high - estimate]]),
            fmt=marker,
            color=color,
            ecolor=color,
            capsize=3,
            markersize=6,
            linewidth=1.3,
        )
    ax.axvline(0, color="#24292f", linewidth=1)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Hedges g (flight minus control)")
    ax.set_title(
        "C  External recurrence: terminal meta unresolved",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="x", alpha=0.18)
    ax.text(
        0.01,
        0.01,
        "OSD-513 is shown as a moderator and was excluded from the terminal meta-analysis.",
        transform=ax.transAxes,
        fontsize=8,
        color="#57606a",
    )


def plot_osd462_panel(ax: plt.Axes, prep: pd.DataFrame, cohort: pd.DataFrame) -> None:
    prep_order = ["UPX", "mRNA", "totRNA"]
    rows = [prep.loc[prep["preparation"] == item].iloc[0] for item in prep_order]
    combined = cohort.loc[
        (cohort["study"] == "OSD-462") & (cohort["context"] == "terminal")
    ].iloc[0]
    labels = prep_order + ["within-animal mean"]
    estimates = [float(row["estimate_g"]) for row in rows] + [
        float(combined["estimate_g"])
    ]
    lows = [float(row["ci_low"]) for row in rows] + [float(combined["ci_low"])]
    highs = [float(row["ci_high"]) for row in rows] + [float(combined["ci_high"])]
    y = np.arange(len(labels))[::-1]
    colors = ["#54aeff", "#238636", "#cf222e", "#8250df"]
    for yi, estimate, low, high, color in zip(y, estimates, lows, highs, colors):
        ax.errorbar(
            estimate,
            yi,
            xerr=np.array([[estimate - low], [high - estimate]]),
            fmt="s",
            color=color,
            ecolor=color,
            capsize=3,
            markersize=6,
            linewidth=1.3,
        )
    ax.axvline(0, color="#24292f", linewidth=1)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Hedges g (flight minus control)")
    ax.set_title(
        "D  OSD-462 changes sign by RNA preparation",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="x", alpha=0.18)


def run(args: argparse.Namespace) -> None:
    config_path = resolve(args.config)
    cfg = yaml.safe_load(config_path.read_text())
    root = resolve(args.outdir or cfg["output_dir"])
    internal = load_json(root / "internal" / "internal_gate_status.json")
    gate_b = load_json(root / "internal" / "compactness" / "gate_B_status.json")
    gate_e = load_json(root / "external" / "gate_E_status.json")

    statuses = [
        ("A", "selection-adjusted internal evidence", bool(internal["gate_A"]["pass"])),
        ("B", "flight-blind recovery", bool(gate_b["pass"])),
        ("C", "mouse/gene influence", bool(internal["gate_C"]["pass"])),
        ("D", "technical/composition/baseline", bool(internal["gate_D"]["pass"])),
        ("E", "terminal external recurrence", bool(gate_e["pass"])),
        ("F", "generic-state distinctness", bool(internal["gate_F"]["pass"])),
    ]
    required = {gate: state for gate, _, state in statuses if gate in "ABCDE"}
    overall_pass = all(required.values())
    decision = {
        "analysis_id": cfg["analysis_id"],
        "decision": (
            "GO_STANDALONE_GREY60" if overall_pass else "NO_GO_RETIRE_STANDALONE_GREY60"
        ),
        "overall_pass": overall_pass,
        "locked_rule": "All of Gates A-E must pass; Gate F controls wording only.",
        "gates": {
            gate: {"label": label, "pass": state}
            for gate, label, state in statuses
        },
        "decisive_failures": [
            gate for gate in ("B", "E") if not dict((x[0], x[2]) for x in statuses)[gate]
        ],
        "supporting_failures": [
            gate for gate in ("A", "C", "D") if not dict((x[0], x[2]) for x in statuses)[gate]
        ],
        "interpretation": (
            "Grey60 is a coherent, influential ISS-terminal-associated 48-gene "
            "score in OSD-771, but it is not recovered as a compact flight-blind "
            "network module and it does not recur convincingly across compatible "
            "terminal kidney cohorts."
        ),
        "permitted_use": (
            "Exploratory OSD-771 context, a case study of activity shift versus "
            "module preservation, or a technical benchmark; not a recurrent "
            "spaceflight kidney module paper."
        ),
        "recommended_next_project": "OSD-462 matched-library robustness",
        "source_files": {
            "locked_config": str(config_path.relative_to(REPO_ROOT)),
            "gate_A_C_D_F": "internal/internal_gate_status.json",
            "gate_B": "internal/compactness/gate_B_status.json",
            "gate_E": "external/gate_E_status.json",
        },
    }

    decision_dir = root / "decision"
    figure_dir = root / "figures"
    decision_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    (decision_dir / "decision.json").write_text(json.dumps(decision, indent=2) + "\n")
    pd.DataFrame(
        [
            {
                "gate": gate,
                "label": label,
                "required_for_go": gate in "ABCDE",
                "pass": state,
            }
            for gate, label, state in statuses
        ]
    ).to_csv(decision_dir / "gate_summary.tsv", sep="\t", index=False)

    grid = pd.read_csv(
        root / "internal" / "compactness" / "flight_blind_grid_summary.tsv",
        sep="\t",
    )
    cohort = pd.read_csv(root / "external" / "cohort_effects.tsv", sep="\t")
    terminal_meta = pd.read_csv(
        root / "external" / "terminal_meta_summary.tsv", sep="\t"
    ).iloc[0]
    prep = pd.read_csv(
        root / "external" / "osd462_preparation_effects.tsv", sep="\t"
    )

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.2), constrained_layout=True)
    fig.suptitle(
        "Grey60 adversarial reanalysis: NO-GO for a standalone recurrent-module paper",
        fontsize=15,
        fontweight="bold",
    )
    plot_gate_panel(axes[0, 0], statuses)
    plot_flight_blind_panel(
        axes[0, 1],
        grid,
        jaccard_threshold=float(
            cfg["go_no_go"]["gate_B"]["flight_blind_jaccard_gte"]
        ),
        p_threshold=float(
            cfg["go_no_go"]["gate_B"]["flight_blind_projected_permutation_p_lte"]
        ),
    )
    plot_external_forest(axes[1, 0], cohort, terminal_meta)
    plot_osd462_panel(axes[1, 1], prep, cohort)
    for suffix in ("png", "svg", "pdf"):
        fig.savefig(
            figure_dir / f"grey60_decision_dashboard.{suffix}",
            dpi=220 if suffix == "png" else None,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(json.dumps(decision, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default="config/grey60_adversarial_reanalysis.yaml"
    )
    parser.add_argument("--outdir", default="")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
