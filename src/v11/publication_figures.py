#!/usr/bin/env python3
"""Publication-ready v11 figures.

The v11 figure labels use manuscript-facing terminology:
  * DCT1 panels show enrichment among DCT1-high parent genes in whole-kidney
    OSD-462 phosphoproteomics.
  * Composition-aware panels separate robust top-decile enrichment from weak
    continuous DCT1-gradient models.
  * Spatial panels are external IRI reference contextualization.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


TIME_ORDER = ["hour4", "hour12", "day2", "day14", "week6"]
TIME_LABELS = {
    "hour4": "4 h",
    "hour12": "12 h",
    "day2": "day 2",
    "day14": "day 14",
    "week6": "week 6",
}
MODEL_ORDER = [
    "M0_raw",
    "M1_parent_protein_abundance",
    "M2_dct_score",
    "M3_endothelial_stromal",
    "M4_full_parent_dct_endothelial_stromal",
    "M5_parent_composition_pc1",
]
MODEL_LABELS = {
    "M0_raw": "M0 raw",
    "M1_parent_protein_abundance": "M1 parent protein",
    "M2_dct_score": "M2 DCT score",
    "M3_endothelial_stromal": "M3 endothelial + stromal",
    "M4_full_parent_dct_endothelial_stromal": "M4 parent + compartments",
    "M5_parent_composition_pc1": "M5 composition PC",
}
PALETTE = {
    "red": "#C63D3D",
    "blue": "#356D9D",
    "teal": "#2A9D8F",
    "gold": "#D9A441",
    "purple": "#7A5AA6",
    "gray": "#6C757D",
    "light_gray": "#E9ECEF",
    "dark": "#1F2933",
}

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "font.family": "sans-serif",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.edgecolor": "#343A40",
            "xtick.color": "#343A40",
            "ytick.color": "#343A40",
            "text.color": "#212529",
        }
    )
    sns.set_theme(style="whitegrid", rc={"grid.color": "#DEE2E6", "grid.linewidth": 0.6})


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save(fig: plt.Figure, out_dir: Path, name: str) -> None:
    ensure(out_dir)
    fig.savefig(out_dir / f"{name}.png", bbox_inches="tight", dpi=300)
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [saved] {name}.png/.pdf")


def clean_axis(ax: plt.Axes, *, xgrid: bool = True, ygrid: bool = False) -> None:
    ax.set_facecolor("#FFFFFF")
    if xgrid:
        ax.grid(True, axis="x", alpha=0.35)
    else:
        ax.grid(False, axis="x")
    if ygrid:
        ax.grid(True, axis="y", alpha=0.25)
    else:
        ax.grid(False, axis="y")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        color=PALETTE["dark"],
    )


def axis_has_visible_data(ax: plt.Axes) -> bool:
    """Return True when an axis has at least one visible artist with data."""
    for line in ax.lines:
        if line.get_visible() and len(line.get_xdata()) and len(line.get_ydata()):
            return True
    for collection in ax.collections:
        if collection.get_visible():
            offsets = getattr(collection, "get_offsets", lambda: [])()
            if len(offsets):
                return True
    for patch in ax.patches:
        if patch.get_visible() and (patch.get_width() or patch.get_height()):
            return True
    return False


def legend_overlaps_visible_data(ax: plt.Axes) -> bool:
    """Conservative figure-QA helper: does the legend box cover plotted data?"""
    legend = ax.get_legend()
    if legend is None or not legend.get_visible():
        return False
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend_bbox = legend.get_window_extent(renderer=renderer)
    points = []
    for line in ax.lines:
        if not line.get_visible():
            continue
        xy = np.column_stack([line.get_xdata(), line.get_ydata()])
        if len(xy):
            points.append(ax.transData.transform(xy))
    for collection in ax.collections:
        if not collection.get_visible():
            continue
        offsets = getattr(collection, "get_offsets", lambda: np.empty((0, 2)))()
        if len(offsets):
            points.append(ax.transData.transform(offsets))
    if not points:
        return False
    xy = np.vstack(points)
    return bool(np.any([legend_bbox.contains(x, y) for x, y in xy]))


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def odds_ci(row: pd.Series) -> tuple[float, float]:
    count_cols = {
        "suppressed_in_flag": ["suppressed_in_flag", "table_suppressed_in_flag"],
        "suppressed_not_flag": ["suppressed_not_flag", "table_suppressed_not_flag"],
        "background_in_flag": ["background_in_flag", "table_background_in_flag"],
        "background_not_flag": ["background_not_flag", "table_background_not_flag"],
    }
    vals = {}
    for key, candidates in count_cols.items():
        for col in candidates:
            if col in row.index and pd.notna(row[col]):
                vals[key] = float(row[col])
                break
        else:
            vals[key] = np.nan
    a = vals["suppressed_in_flag"] + 0.5
    b = vals["suppressed_not_flag"] + 0.5
    c = vals["background_in_flag"] + 0.5
    d = vals["background_not_flag"] + 0.5
    if not np.all(np.isfinite([a, b, c, d])):
        return np.nan, np.nan
    log_or = np.log((a * d) / (b * c))
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return float(np.exp(log_or - 1.96 * se)), float(np.exp(log_or + 1.96 * se))


def q_label(q: float) -> str:
    if not np.isfinite(q):
        return "q=NA"
    if q < 1e-3:
        return f"q={q:.1e}"
    return f"q={q:.3f}"


def fig_cross_cohort_recurrence(run_root: Path, out_dir: Path) -> bool:
    """Clean replacement for the legacy cross-cohort figure.

    The older static Figure 1 mixed leave-one-pathway-out diagnostics with the
    headline OSD-513 cosine, which left a stale label in the rendered panel.
    This version uses the canonical TSV artifacts and labels each statistic
    explicitly.
    """
    cross_path = (
        PROJECT_ROOT
        / "data/results/run_20260518_201823_2500g/contrast_vectors/"
        / "cross_osdr_recurrence/cross_osdr_alignment_summary.tsv"
    )
    osd462_path = PROJECT_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_rna_recurrence.tsv"
    osd462_effect_path = PROJECT_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_rna_pathway_effects.tsv"
    osd462_loo_path = PROJECT_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_rna_recurrence_loo.tsv"
    required = [cross_path, osd462_path, osd462_effect_path, osd462_loo_path]
    if any(not p.exists() for p in required):
        missing = [str(p) for p in required if not p.exists()]
        print(f"  [skip] missing recurrence inputs: {missing}")
        return False

    cross = read_tsv(cross_path)
    osd462 = read_tsv(osd462_path).iloc[0]
    effects = read_tsv(osd462_effect_path)
    loo = read_tsv(osd462_loo_path)

    rows = []
    for _, row in cross[cross["study"].eq("OSD-513")].iterrows():
        rows.append(
            {
                "label": f"RRRM-2 {row['arm']} vs OSD-513",
                "cosine": float(row["point_estimate"]),
                "ci_low": float(row["ci_low"]),
                "ci_high": float(row["ci_high"]),
                "color": PALETTE["red"] if row["arm"] == "ISS-T" else PALETTE["blue"],
            }
        )
    rows.append(
        {
            "label": "RRRM-2 ISS-T vs OSD-462 RNA",
            "cosine": float(osd462["point_cosine"]),
            "ci_low": float(osd462["ci_low"]),
            "ci_high": float(osd462["ci_high"]),
            "color": PALETTE["teal"],
        }
    )
    align = pd.DataFrame(rows)

    fig = plt.figure(figsize=(10.8, 7.1))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1], height_ratios=[1, 1], wspace=0.32, hspace=0.42)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])

    y = np.arange(len(align))
    for i, row in align.iterrows():
        ax1.errorbar(
            row["cosine"],
            i,
            xerr=[[row["cosine"] - row["ci_low"]], [row["ci_high"] - row["cosine"]]],
            fmt="o",
            color=row["color"],
            ecolor=row["color"],
            capsize=3,
            ms=7,
        )
        ax1.text(row["cosine"], i - 0.16, f"{row['cosine']:.2f}", ha="center", va="bottom", fontsize=8)
    ax1.axvline(0, color="#495057", lw=1, ls="--")
    ax1.set_yticks(y)
    ax1.set_yticklabels(align["label"])
    ax1.invert_yaxis()
    ax1.set_xlim(-0.85, 1.0)
    ax1.set_xlabel("Pathway-vector cosine")
    ax1.set_title("Canonical recurrence estimates")
    clean_axis(ax1, xgrid=True)
    panel_label(ax1, "A")

    colors = np.where(effects["sign_agree"].astype(bool), PALETTE["teal"], PALETTE["gray"])
    ax2.scatter(
        effects["rrrm2_iss_t_pathway_effect"],
        effects["osd462_rna_pathway_effect"],
        s=45,
        color=colors,
        edgecolor="white",
        linewidth=0.6,
    )
    for _, row in effects.iterrows():
        label = str(row["pathway"]).replace("_", " ")
        if row["pathway"] in {"dct_ncc_wnk_transport", "ecm_organization", "mmp_adam_proteolysis"}:
            ax2.text(row["rrrm2_iss_t_pathway_effect"], row["osd462_rna_pathway_effect"], label, fontsize=7)
    ax2.axhline(0, color="#adb5bd", lw=0.8)
    ax2.axvline(0, color="#adb5bd", lw=0.8)
    ax2.set_xlabel("RRRM-2 ISS-T RNA effect")
    ax2.set_ylabel("OSD-462 RNA effect")
    ax2.set_title("OSD-462 pathway concordance")
    clean_axis(ax2, xgrid=True, ygrid=True)
    panel_label(ax2, "B")

    loo = loo.sort_values("cosine")
    labels = [str(x).replace("_", " ") for x in loo["dropped_pathway"]]
    ax3.barh(np.arange(len(loo)), loo["cosine"], color=PALETTE["teal"], alpha=0.85)
    ax3.axvline(float(osd462["point_cosine"]), color=PALETTE["red"], lw=1.2, ls="--", label=f"Full cosine {float(osd462['point_cosine']):.2f}")
    ax3.set_yticks(np.arange(len(loo)))
    ax3.set_yticklabels(labels, fontsize=7)
    ax3.set_xlim(0.75, 0.93)
    ax3.set_xlabel("Cosine after dropping one pathway")
    ax3.set_title("OSD-462 leave-one-pathway-out")
    ax3.legend(frameon=True, loc="lower right", fontsize=8)
    clean_axis(ax3, xgrid=True)
    panel_label(ax3, "C")

    fig.suptitle("Cross-cohort RNA recurrence around the spaceflight kidney remodeling state", fontsize=13, fontweight="bold", y=0.985)
    fig.text(
        0.5,
        0.01,
        "OSD-513 uses nine shared pathway features; OSD-462 uses the 11-pathway multi-omic anchor panel.",
        ha="center",
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "fig1_main_result_multipanel")
    return True


def fig_h2_primary_enrichment(run_root: Path, out_dir: Path) -> bool:
    path = run_root / "h2_enrichment/h2_dct1_sensitivity_summary.tsv"
    if not path.exists():
        print(f"  [skip] missing {path}")
        return False

    df = read_tsv(path)
    keep_tests = ["fisher_dct1_top_decile", "fisher_dct2_bottom_decile", "fisher_dct1_top_quartile"]
    keep_analyses = [
        "primary_p05",
        "exclude_anchor_genes",
        "exclude_ncc_sites",
        "composite_sites_excluded",
        "one_site_per_parent_gene",
        "single_position_one_site_per_parent_gene",
        "strict_q10",
    ]
    sub = df[df["test"].isin(keep_tests) & df["analysis"].isin(keep_analyses)].copy()
    if sub.empty:
        print("  [skip] no DCT1 enrichment rows")
        return False
    sub["analysis"] = pd.Categorical(sub["analysis"], keep_analyses, ordered=True)
    sub["test_label"] = sub["test"].map(
        {
            "fisher_dct1_top_decile": "DCT1 top decile",
            "fisher_dct2_bottom_decile": "DCT2-bottom decile",
            "fisher_dct1_top_quartile": "DCT1 top quartile",
        }
    )
    sub = sub.sort_values(["analysis", "test_label"])
    ci = sub.apply(odds_ci, axis=1, result_type="expand")
    sub["ci_low"] = ci[0]
    sub["ci_high"] = ci[1]
    y_labels = {
        "primary_p05": "Primary p<0.05",
        "exclude_anchor_genes": "Anchor genes excluded",
        "exclude_ncc_sites": "NCC sites excluded",
        "composite_sites_excluded": "Composite sites excluded",
        "one_site_per_parent_gene": "One row per parent gene",
        "single_position_one_site_per_parent_gene": "Single-position row per gene",
        "strict_q10": "Strict q<0.10",
    }

    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    offsets = {"DCT1 top decile": -0.18, "DCT2-bottom decile": 0.0, "DCT1 top quartile": 0.18}
    colors = {"DCT1 top decile": PALETTE["red"], "DCT2-bottom decile": PALETTE["teal"], "DCT1 top quartile": PALETTE["blue"]}
    y_base = np.arange(len(keep_analyses))
    for label, part in sub.groupby("test_label", sort=False):
        if label not in offsets:
            continue
        ys = np.array([keep_analyses.index(a) for a in part["analysis"].astype(str)]) + offsets[label]
        xerr = np.vstack([part["statistic"] - part["ci_low"], part["ci_high"] - part["statistic"]])
        ax.errorbar(
            part["statistic"],
            ys,
            xerr=xerr,
            fmt="o",
            color=colors[label],
            ecolor=colors[label],
            elinewidth=1.2,
            capsize=2,
            ms=5,
            label=label,
        )
    ax.axvline(1, color="#495057", lw=1, ls="--")
    ax.set_yticks(y_base)
    ax.set_yticklabels([y_labels[a] for a in keep_analyses])
    ax.invert_yaxis()
    ax.set_xlabel("Odds ratio for flight-suppressed phosphosites")
    ax.set_title("Distal-nephron subtype-prior parent-gene enrichment")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), ncol=2, frameon=True)
    clean_axis(ax)
    ax.text(
        0.01,
        -0.18,
        "External DCT1/DCT2 reference; OSD-462 is whole-kidney phosphoproteomics.",
        transform=ax.transAxes,
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "v11_dct1_parent_gene_enrichment")
    # Alias for the revised claim hierarchy while preserving the historical filename.
    fig_alias = plt.figure(figsize=(0.1, 0.1))
    plt.close(fig_alias)
    for ext in ["png", "pdf"]:
        src = out_dir / f"v11_dct1_parent_gene_enrichment.{ext}"
        dst = out_dir / f"v11_distal_nephron_prior_enrichment.{ext}"
        if src.exists():
            dst.write_bytes(src.read_bytes())
    return True


def fig_h2_composition_robustness(run_root: Path, out_dir: Path) -> bool:
    base = run_root / "h2_composition_adjusted"
    enrich_path = base / "h2_composition_adjusted_suppression_enrichment_single.tsv"
    effect_path = base / "h2_composition_effect_level_dct1_ladder_single.tsv"
    fixed_path = base / "h2_composition_site_fixed_interaction_ladder_single.tsv"
    if not enrich_path.exists() or not effect_path.exists() or not fixed_path.exists():
        print("  [skip] missing composition-aware outputs")
        return False

    enrich = read_tsv(enrich_path)
    enrich = enrich[enrich["model"].isin(MODEL_ORDER)].copy()
    enrich["model"] = pd.Categorical(enrich["model"], MODEL_ORDER, ordered=True)
    enrich["flag_label"] = enrich["flag"].map(
        {"dct1_top_decile": "DCT1 top decile", "dct1_top_quartile": "DCT1 top quartile"}
    )
    ci = enrich.apply(odds_ci, axis=1, result_type="expand")
    enrich["ci_low"] = ci[0]
    enrich["ci_high"] = ci[1]

    effects = read_tsv(effect_path)
    effects = effects[
        (effects["model"].isin(MODEL_ORDER)) & (effects["second_stage_adjustment"].eq("full_second_stage"))
    ].copy()
    effects["model"] = pd.Categorical(effects["model"], MODEL_ORDER, ordered=True)
    effects["source"] = "Two-stage adjusted effect"

    fixed = read_tsv(fixed_path)
    fixed_model_map = {
        "LM0_raw_site_fixed": "M0_raw",
        "LM1_parent_site_fixed": "M1_parent_protein_abundance",
        "LM2_dct_site_fixed": "M2_dct_score",
        "LM3_endothelial_stromal_site_fixed": "M3_endothelial_stromal",
        "LM4_full_site_fixed": "M4_full_parent_dct_endothelial_stromal",
        "LM5_parent_composition_pc1_site_fixed": "M5_parent_composition_pc1",
    }
    fixed["model_clean"] = fixed["model"].map(fixed_model_map)
    fixed = fixed[fixed["model_clean"].isin(MODEL_ORDER)].copy()
    fixed["model"] = pd.Categorical(fixed["model_clean"], MODEL_ORDER, ordered=True)
    fixed["source"] = "One-shot site-fixed interaction"
    coef = pd.concat(
        [
            effects[["model", "source", "coefficient", "se", "q_value"]],
            fixed[["model", "source", "coefficient", "se", "q_value"]],
        ],
        ignore_index=True,
    )
    coef["ci_low"] = coef["coefficient"] - 1.96 * coef["se"]
    coef["ci_high"] = coef["coefficient"] + 1.96 * coef["se"]

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 5.4), gridspec_kw={"width_ratios": [1.08, 1]})
    ax = axes[0]
    offsets = {"DCT1 top decile": -0.12, "DCT1 top quartile": 0.12}
    colors = {"DCT1 top decile": PALETTE["red"], "DCT1 top quartile": PALETTE["blue"]}
    y_base = np.arange(len(MODEL_ORDER))
    for label, part in enrich.groupby("flag_label", sort=False):
        ys = np.array([MODEL_ORDER.index(str(m)) for m in part["model"]]) + offsets[label]
        xerr = np.vstack([part["odds_ratio"] - part["ci_low"], part["ci_high"] - part["odds_ratio"]])
        ax.errorbar(
            part["odds_ratio"],
            ys,
            xerr=xerr,
            fmt="o",
            color=colors[label],
            ecolor=colors[label],
            elinewidth=1.15,
            capsize=2,
            ms=5,
            label=label,
        )
    m4 = enrich[(enrich["model"].astype(str).eq("M4_full_parent_dct_endothelial_stromal")) & (enrich["flag"].eq("dct1_top_decile"))]
    if not m4.empty:
        row = m4.iloc[0]
        y = MODEL_ORDER.index("M4_full_parent_dct_endothelial_stromal") - 0.12
        ax.annotate(
            f"M4 OR={row['odds_ratio']:.2f}\n{q_label(row['q_value'])}",
            xy=(row["odds_ratio"], y),
            xytext=(1.72, y - 0.35),
            arrowprops={"arrowstyle": "->", "lw": 1, "color": PALETTE["gray"]},
            fontsize=8,
            color=PALETTE["dark"],
        )
    ax.axvline(1, color="#495057", lw=1, ls="--")
    ax.set_yticks(y_base)
    ax.set_yticklabels([MODEL_LABELS[m] for m in MODEL_ORDER])
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted suppression enrichment odds ratio")
    ax.set_title("Top-decile enrichment persists after adjustment")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.16), frameon=True)
    clean_axis(ax)
    panel_label(ax, "A")

    ax = axes[1]
    source_offsets = {"Two-stage adjusted effect": -0.11, "One-shot site-fixed interaction": 0.11}
    source_colors = {
        "Two-stage adjusted effect": PALETTE["purple"],
        "One-shot site-fixed interaction": PALETTE["teal"],
    }
    for label, part in coef.groupby("source", sort=False):
        ys = np.array([MODEL_ORDER.index(str(m)) for m in part["model"]]) + source_offsets[label]
        ax.errorbar(
            part["coefficient"],
            ys,
            xerr=np.vstack([part["coefficient"] - part["ci_low"], part["ci_high"] - part["coefficient"]]),
            fmt="o",
            color=source_colors[label],
            ecolor=source_colors[label],
            elinewidth=1.15,
            capsize=2,
            ms=5,
            label=label,
        )
    ax.axvline(0, color="#495057", lw=1, ls="--")
    ax.set_yticks(y_base)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel("DCT1 reference coefficient\nnegative = stronger flight suppression")
    ax.set_title("Continuous gradient models remain weak")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    clean_axis(ax)
    panel_label(ax, "B")
    fig.suptitle("Parent-protein and composition-score sensitivity", fontsize=13, fontweight="bold", y=1.02)
    save(fig, out_dir, "v11_parent_protein_composition_sensitivity")
    return True


def fig_h3_mediation_forest(run_root: Path, out_dir: Path) -> bool:
    path = run_root / "h3_mediation/h3_mediation_model_summary.tsv"
    if not path.exists():
        print(f"  [skip] missing {path}")
        return False
    df = read_tsv(path)
    sub = df[df["parameter"].eq("indirect")].copy()
    if sub.empty:
        print("  [skip] no indirect mediation rows")
        return False
    label_map = {
        "endothelial": "Endothelial",
        "stromal_fibroblast": "Stromal/fibroblast",
        "dct_identity": "DCT identity",
        "matrix_endothelial_composite": "Matrix/endothelial composite",
    }
    order = ["endothelial", "matrix_endothelial_composite", "stromal_fibroblast", "dct_identity"]
    sub["mediator"] = pd.Categorical(sub["mediator"], order, ordered=True)
    sub = sub.sort_values("mediator")
    y = np.arange(len(sub))
    colors = [PALETTE["red"] if hi < 0 else PALETTE["gray"] for hi in sub["ci_high"]]

    fig, ax = plt.subplots(figsize=(7.8, 4.2))
    for yi, (_, row), color in zip(y, sub.iterrows(), colors):
        ax.errorbar(
            row["posterior_median"],
            yi,
            xerr=[[row["posterior_median"] - row["ci_low"]], [row["ci_high"] - row["posterior_median"]]],
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=1.5,
            capsize=3,
            ms=6,
        )
        ax.text(
            0.02,
            yi,
            f"P(<0)={row['p_less_than_zero']:.3f}",
            va="center",
            ha="left",
            fontsize=8,
            color=PALETTE["gray"],
        )
    ax.axvline(0, color="#495057", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([label_map.get(str(m), str(m)) for m in sub["mediator"].astype(str)])
    ax.invert_yaxis()
    ax.set_xlabel("Covariance-decomposition path summary for NCC regulatory phospho score")
    ax.set_title("Exploratory remodeling/transporter covariance decomposition")
    clean_axis(ax)
    ax.text(
        0.01,
        -0.2,
        "Approximate Bayesian OLS; cross-sectional bulk tissue; not causal mediation.",
        transform=ax.transAxes,
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "v11_exploratory_mediation_forest")
    for ext in ["png", "pdf"]:
        src = out_dir / f"v11_exploratory_mediation_forest.{ext}"
        dst = out_dir / f"v11_exploratory_covariance_decomposition.{ext}"
        if src.exists():
            dst.write_bytes(src.read_bytes())
    return True


def fig_spatial_reference(run_root: Path, out_dir: Path) -> bool:
    base = run_root / "spatial_reference"
    pathway_path = base / "visium_timepoint_pathway_cosines.tsv"
    gene_path = base / "visium_timepoint_gene_cosines.tsv"
    niche_path = base / "visium_niche_cosines.tsv"
    if not pathway_path.exists() or not gene_path.exists() or not niche_path.exists():
        print("  [skip] missing spatial reference outputs")
        return False

    pathway = read_tsv(pathway_path)
    gene = read_tsv(gene_path)
    niche = read_tsv(niche_path)

    pathway = pathway[
        pathway["comparison"].eq("visium_timepoint_minus_sham") & pathway["timepoint"].isin(TIME_ORDER)
    ].copy()
    gene = gene[gene["comparison"].eq("visium_timepoint_minus_sham") & gene["timepoint"].isin(TIME_ORDER)].copy()
    pathway["timepoint"] = pd.Categorical(pathway["timepoint"], TIME_ORDER, ordered=True)
    gene["timepoint"] = pd.Categorical(gene["timepoint"], TIME_ORDER, ordered=True)

    fig = plt.figure(figsize=(12.6, 7.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.18], hspace=0.42, wspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    vector_colors = {"osd462_rna": PALETTE["red"], "rrrm2_iss_t": PALETTE["blue"]}
    vector_labels = {"osd462_rna": "OSD-462 RNA", "rrrm2_iss_t": "RRRM-2 ISS-T"}
    for vector, part in pathway.groupby("spaceflight_pathway_vector"):
        part = part.sort_values("timepoint")
        ax1.plot(
            [TIME_LABELS[str(t)] for t in part["timepoint"]],
            part["cosine"],
            marker="o",
            lw=2,
            color=vector_colors.get(vector, PALETTE["gray"]),
            label=vector_labels.get(vector, vector),
        )
    ax1.axhline(0, color="#495057", lw=0.8)
    ax1.set_ylim(0, 0.88)
    ax1.set_ylabel("Pathway cosine vs sham")
    ax1.set_title("Visium IRI timepoint pathway alignment")
    ax1.legend(frameon=True, loc="lower right")
    clean_axis(ax1, ygrid=True)
    panel_label(ax1, "A")

    for vector, part in gene.groupby("spaceflight_vector"):
        part = part.sort_values("timepoint")
        ax2.plot(
            [TIME_LABELS[str(t)] for t in part["timepoint"]],
            part["cosine"],
            marker="o",
            lw=2,
            color=vector_colors.get(vector, PALETTE["gray"]),
            label=vector_labels.get(vector, vector),
        )
    ax2.axhline(0, color="#495057", lw=0.8)
    ax2.set_ylabel("Genome-wide gene cosine vs sham")
    ax2.set_title("Gene-level alignment is modest")
    clean_axis(ax2, ygrid=True)
    panel_label(ax2, "B")

    sub = niche[
        niche["spaceflight_vector"].eq("osd462_rna")
        & niche["comparison"].eq("visium_niche_pathway_minus_sham_same_niche")
        & niche["timepoint"].isin(TIME_ORDER)
    ].copy()
    if not sub.empty:
        niche_order = (
            sub.groupby("niche")["cosine"]
            .max()
            .sort_values(ascending=False)
            .index.tolist()
        )
        sub["timepoint"] = pd.Categorical(sub["timepoint"], TIME_ORDER, ordered=True)
        sub["niche"] = pd.Categorical(sub["niche"], niche_order, ordered=True)
        mat = sub.pivot_table(index="niche", columns="timepoint", values="cosine", observed=False)
        mat = mat.rename(index=lambda x: str(x).replace("_", " "), columns=TIME_LABELS)
        sns.heatmap(
            mat,
            ax=ax3,
            cmap=sns.light_palette(PALETTE["red"], as_cmap=True),
            vmin=0.55,
            vmax=0.85,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Pathway cosine"},
        )
        ax3.set_xlabel("IRI timepoint")
        ax3.set_ylabel("")
        ax3.set_title("OSD-462 pathway vector alignment by marker-enriched Visium niche")
        panel_label(ax3, "C")

    fig.suptitle("External IRI spatial reference contextualization", fontsize=13, fontweight="bold", y=1.01)
    fig.text(
        0.5,
        0.005,
        "GSE269622/GSE269719 contextualize known kidney injury/repair spatial programs; they do not localize RR-10 spaceflight lesions.",
        ha="center",
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "v11_spatial_reference_projection")
    return True


def fig_xenium_inventory(run_root: Path, out_dir: Path) -> bool:
    path = run_root / "spatial_reference/xenium_annotation_inventory.tsv"
    verdict_path = run_root / "spatial_reference/spatial_reference_projection_verdict.json"
    if not path.exists():
        print(f"  [skip] missing {path}")
        return False
    df = read_tsv(path)
    cell = df[df["column"].eq("celltype_plot")].copy()
    cn = df[df["column"].eq("CN")].copy()
    if cell.empty and cn.empty:
        print("  [skip] no Xenium annotation rows")
        return False
    cell = cell.sort_values("count", ascending=False).head(12)
    cn = cn.sort_values("count", ascending=False)
    n_cells = None
    n_genes = None
    if verdict_path.exists():
        verdict = json.loads(verdict_path.read_text())
        n_cells = verdict.get("xenium", {}).get("n_cells")
        n_genes = verdict.get("xenium", {}).get("n_genes")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.4))
    ax = axes[0]
    sns.barplot(data=cell, y="value", x="count", ax=ax, color=PALETTE["blue"])
    ax.set_xlabel("Cells")
    ax.set_ylabel("")
    ax.set_title("Xenium cell-type annotation inventory")
    clean_axis(ax)
    panel_label(ax, "A")

    ax = axes[1]
    sns.barplot(data=cn, y="value", x="count", ax=ax, color=PALETTE["teal"])
    ax.set_xlabel("Cells")
    ax.set_ylabel("")
    ax.set_title("Xenium cellular-neighborhood inventory")
    clean_axis(ax)
    panel_label(ax, "B")

    if n_cells and n_genes:
        fig.suptitle(
            f"GSE269719 Xenium targeted-panel context: {n_cells:,} cells, {n_genes} genes",
            fontsize=13,
            fontweight="bold",
            y=1.02,
        )
    save(fig, out_dir, "v11_xenium_annotation_inventory")
    return True


def fig_perturbation_triangulation(run_root: Path, out_dir: Path) -> bool:
    lowk_path = run_root / "perturbation/lowk_dct_alignment_summary.tsv"
    target_path = run_root / "perturbation/lowk_target_gene_table.tsv"
    occ_path = run_root / "h2_occupancy/h2_occupancy_dct1_enrichment.tsv"
    spatial_path = run_root / "spatial_reference/visium_dct_transport_timepoint_summary.tsv"
    if not all(p.exists() for p in [lowk_path, target_path, occ_path, spatial_path]):
        print("  [skip] missing perturbation-triangulation outputs")
        return False

    lowk = read_tsv(lowk_path)
    target = read_tsv(target_path)
    occ = read_tsv(occ_path)
    spatial = read_tsv(spatial_path)

    fig = plt.figure(figsize=(13.8, 10.2))
    gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.48)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    gs_d = gs[1, 1].subgridspec(2, 1, hspace=0.42)
    ax4a = fig.add_subplot(gs_d[0, 0])
    ax4b = fig.add_subplot(gs_d[1, 0], sharex=ax4a)

    sub = lowk[
        lowk["spaceflight_vector"].isin(["osd462_rna", "rrrm2_iss_t_rna", "osd513_rna_stat"])
        & lowk["gene_subset"].isin(["all_overlap", "dct1_top_decile_genes", "dct2_top_decile_genes", "transport_target_genes"])
    ].copy()
    vector_labels = {
        "osd462_rna": "OSD-462 RNA",
        "rrrm2_iss_t_rna": "RRRM-2 ISS-T RNA",
        "osd513_rna_stat": "OSD-513 RNA",
    }
    subset_order = ["all_overlap", "dct1_top_decile_genes", "dct2_top_decile_genes", "transport_target_genes"]
    subset_labels = {
        "all_overlap": "All DCT-prior genes",
        "dct1_top_decile_genes": "DCT1-prior top decile",
        "dct2_top_decile_genes": "DCT2-prior top decile",
        "transport_target_genes": "Transport targets",
    }
    if not sub.empty:
        sub["gene_subset"] = pd.Categorical(sub["gene_subset"], subset_order, ordered=True)
        pivot = sub.pivot_table(index="gene_subset", columns="spaceflight_vector", values="cosine", observed=False)
        pivot = pivot.rename(index=subset_labels, columns=vector_labels)
        sns.heatmap(
            pivot,
            ax=ax1,
            cmap="vlag",
            center=0,
            vmin=-0.45,
            vmax=0.45,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar_kws={"label": "Cosine"},
        )
        ax1.set_xlabel("")
        ax1.set_ylabel("")
        ax1.set_title("Low-K DCT response vs spaceflight RNA")
        ax1.tick_params(axis="x", labelrotation=35, labelsize=8)
        ax1.tick_params(axis="y", labelsize=8)
        for tick in ax1.get_xticklabels():
            tick.set_ha("right")
        panel_label(ax1, "A")

    available_cols = ["lowk_effect_kd_minus_nk", "osd462_rna", "rrrm2_iss_t_rna", "osd513_rna_stat"]
    heat_cols = [c for c in available_cols if c in target.columns]
    heat = target.set_index("gene_symbol")[heat_cols].copy()
    heat = heat.apply(lambda col: np.clip((col - col.mean()) / col.std(ddof=0), -2, 2), axis=0)
    heat = heat.rename(
        columns={
            "lowk_effect_kd_minus_nk": "Low-K DCT",
            "osd462_rna": "OSD-462 RNA",
            "rrrm2_iss_t_rna": "RRRM-2 ISS-T",
            "osd513_rna_stat": "OSD-513 RNA",
        }
    )
    sns.heatmap(
        heat,
        ax=ax2,
        cmap="vlag",
        center=0,
        vmin=-2,
        vmax=2,
        annot=False,
        linewidths=0.4,
        cbar_kws={"label": "Column z-score"},
    )
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.set_title("Focused DCT/WNK target-gene directions")
    ax2.tick_params(axis="x", labelrotation=35, labelsize=8)
    ax2.tick_params(axis="y", labelsize=8)
    for tick in ax2.get_xticklabels():
        tick.set_ha("right")
    panel_label(ax2, "B")

    occ_sub = occ[
        occ["analysis"].isin(["occupancy_p05", "occupancy_q10", "occupancy_one_row_per_parent_gene_p05", "occupancy_single_position_one_row_per_parent_gene_p05"])
        & occ["flag"].isin(["dct1_top_decile", "dct1_top_quartile"])
    ].copy()
    if not occ_sub.empty:
        occ_sub["analysis_label"] = occ_sub["analysis"].map(
            {
                "occupancy_p05": "Nominal p<0.05",
                "occupancy_q10": "Strict q<0.10",
                "occupancy_one_row_per_parent_gene_p05": "One row per gene",
                "occupancy_single_position_one_row_per_parent_gene_p05": "Single-position row per gene",
            }
        )
        occ_sub["flag_label"] = occ_sub["flag"].map(
            {"dct1_top_decile": "DCT1 top decile", "dct1_top_quartile": "DCT1 top quartile"}
        )
        order = ["Nominal p<0.05", "Strict q<0.10", "One row per gene", "Single-position row per gene"]
        y_base = np.arange(len(order))
        offsets = {"DCT1 top decile": -0.12, "DCT1 top quartile": 0.12}
        colors = {"DCT1 top decile": PALETTE["red"], "DCT1 top quartile": PALETTE["blue"]}
        for label, part in occ_sub.groupby("flag_label", sort=False):
            ys = np.array([order.index(x) for x in part["analysis_label"]]) + offsets[label]
            ax3.errorbar(
                part["odds_ratio"],
                ys,
                xerr=np.vstack([part["odds_ratio"] - part["ci_low"], part["ci_high"] - part["odds_ratio"]]),
                fmt="o",
                ms=5,
                color=colors[label],
                ecolor=colors[label],
                capsize=2,
                label=label,
            )
        ax3.axvline(1, color="#495057", lw=1, ls="--")
        ax3.set_yticks(y_base)
        ax3.set_yticklabels(order)
        ax3.invert_yaxis()
        ax3.set_xlabel("Odds ratio")
        ax3.set_title("Parent-protein-normalized phosphosite enrichment")
        ax3.legend(loc="upper left", bbox_to_anchor=(0.0, -0.20), ncol=2, frameon=True)
        clean_axis(ax3)
        panel_label(ax3, "C")

    sp = spatial[
        spatial["group"].isin(["all_spots", "dct_marker_top_quartile", "dct_adjacent_spots"])
        & spatial["comparison"].str.endswith("minus_sham_same_group")
        & spatial["timepoint"].isin(TIME_ORDER)
    ].copy()
    if not sp.empty:
        group_labels = {
            "all_spots": "All spots",
            "dct_marker_top_quartile": "DCT-high spots",
            "dct_adjacent_spots": "DCT-adjacent spots",
        }
        colors = {
            "All spots": PALETTE["gray"],
            "DCT-high spots": PALETTE["blue"],
            "DCT-adjacent spots": PALETTE["teal"],
        }
        sp["timepoint"] = pd.Categorical(sp["timepoint"], TIME_ORDER, ordered=True)
        sp["group_label"] = sp["group"].map(group_labels)
        for label, part in sp.sort_values("timepoint").groupby("group_label", sort=False):
            ax = ax4a if label == "DCT-high spots" else ax4b
            ax.plot(
                [TIME_LABELS[str(t)] for t in part["timepoint"]],
                part["mean_difference_vs_sham"],
                marker="o",
                lw=1.8,
                color=colors.get(label, PALETTE["gray"]),
                label=label,
            )
        for ax, title in [(ax4a, "DCT-marker-high spots"), (ax4b, "All and DCT-adjacent spots")]:
            ax.axhline(0, color="#495057", lw=1, ls="--")
            ax.set_ylabel("Difference vs sham")
            ax.set_title(title)
            ax.legend(frameon=True, loc="center left", bbox_to_anchor=(1.02, 0.5))
            clean_axis(ax, ygrid=True)
        panel_label(ax4a, "D")
        ax4b.set_xlabel("IRI timepoint")
        plt.setp(ax4a.get_xticklabels(), visible=False)

    fig.suptitle("Public perturbation-reference checks around the spaceflight DCT/phosphoproteomic signature", fontsize=13, fontweight="bold", y=0.985)
    fig.text(
        0.5,
        0.006,
        "Low-K uses DCT-enriched pseudobulk; spatial scores use external IRI Visium and are predictive context, not spaceflight spatial validation.",
        ha="center",
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "v11_perturbation_triangulation")
    return True


def generate_all(run_root: Path, out_dir: Path | None = None) -> int:
    configure_style()
    out_dir = out_dir or run_root / "figures" / "v11"
    ensure(out_dir)
    print(f"[v11-figures] output: {out_dir}")

    made = 0
    for fn in [
        fig_cross_cohort_recurrence,
        fig_h2_primary_enrichment,
        fig_h2_composition_robustness,
        fig_h3_mediation_forest,
        fig_spatial_reference,
        fig_perturbation_triangulation,
        fig_xenium_inventory,
    ]:
        try:
            made += int(fn(run_root, out_dir))
        except Exception as exc:
            print(f"  [ERROR] {fn.__name__}: {exc}")

    readme = out_dir / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# V11 Publication Figures",
                "",
                "Generated from v11 DCT1/phosphoproteome/mediation outputs.",
                "",
                "DCT1 panels show enrichment among DCT1-high parent genes in whole-kidney OSD-462 phosphoproteomics.",
                "Spatial panels use external IRI references for context.",
                "",
                f"Figures generated: {made}",
                "",
            ]
        )
    )
    print(f"[v11-figures] generated {made} figure sets")
    return made


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate v11 publication figures")
    parser.add_argument("--run-root", "--results_dir", dest="run_root", required=True)
    parser.add_argument("--out-dir", default=None)
    args = parser.parse_args()
    run_root = Path(args.run_root).resolve()
    out_dir = Path(args.out_dir).resolve() if args.out_dir else None
    generate_all(run_root, out_dir)


if __name__ == "__main__":
    main()
