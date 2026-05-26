#!/usr/bin/env python3
"""Publication-ready v11 figures.

The v11 figures are intentionally claim-disciplined:
  * H2 panels show DCT1-prior enrichment, not DCT-isolated phosphoproteomics.
  * Composition-aware panels separate robust top-decile enrichment from weak
    continuous DCT1-gradient models.
  * Spatial panels are external IRI reference contextualization, not RR-10
    spatial validation.
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
    "M4_full_parent_dct_endothelial_stromal": "M4 full conservative",
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


def fig_h2_primary_enrichment(run_root: Path, out_dir: Path) -> bool:
    path = run_root / "h2_enrichment/h2_dct1_sensitivity_summary.tsv"
    if not path.exists():
        print(f"  [skip] missing {path}")
        return False

    df = read_tsv(path)
    keep_tests = ["fisher_dct1_top_decile", "fisher_dct1_top_quartile"]
    keep_analyses = [
        "primary_p05",
        "exclude_anchor_genes",
        "exclude_ncc_sites",
        "single_site_only_p05",
        "strict_q10",
    ]
    sub = df[df["test"].isin(keep_tests) & df["analysis"].isin(keep_analyses)].copy()
    if sub.empty:
        print("  [skip] no H2 enrichment rows")
        return False
    sub["analysis"] = pd.Categorical(sub["analysis"], keep_analyses, ordered=True)
    sub["test_label"] = sub["test"].map(
        {
            "fisher_dct1_top_decile": "DCT1 top decile",
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
        "single_site_only_p05": "Single-site only",
        "strict_q10": "Strict q<0.10",
    }

    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    offsets = {"DCT1 top decile": -0.12, "DCT1 top quartile": 0.12}
    colors = {"DCT1 top decile": PALETTE["red"], "DCT1 top quartile": PALETTE["blue"]}
    y_base = np.arange(len(keep_analyses))
    for label, part in sub.groupby("test_label", sort=False):
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
    ax.set_title("H2 primary DCT1-prior enrichment")
    ax.legend(loc="lower right", frameon=True)
    clean_axis(ax)
    ax.text(
        0.01,
        -0.18,
        "Reference-prior enrichment; OSD-462 is whole-kidney phosphoproteomics.",
        transform=ax.transAxes,
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "v11_h2_primary_dct1_enrichment")
    return True


def fig_h2_composition_robustness(run_root: Path, out_dir: Path) -> bool:
    base = run_root / "h2_composition_adjusted"
    enrich_path = base / "h2_composition_adjusted_suppression_enrichment_single.tsv"
    effect_path = base / "h2_composition_effect_level_dct1_ladder_single.tsv"
    fixed_path = base / "h2_composition_site_fixed_interaction_ladder_single.tsv"
    if not enrich_path.exists() or not effect_path.exists() or not fixed_path.exists():
        print("  [skip] missing H2 composition-aware outputs")
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
    ax.set_title("Subset enrichment survives conservative adjustment")
    ax.legend(loc="lower right", frameon=True)
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
    ax.set_xlabel("DCT1-prior coefficient\nnegative = stronger flight suppression")
    ax.set_title("Continuous gradient models remain weak")
    ax.legend(loc="lower left", frameon=True)
    clean_axis(ax)
    panel_label(ax, "B")
    fig.suptitle("H2 composition-aware robustness", fontsize=13, fontweight="bold", y=1.02)
    save(fig, out_dir, "v11_h2_composition_robustness")
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
    ax.set_xlabel("Indirect effect on NCC regulatory phospho score")
    ax.set_title("H3 mediation-specified remodeling/transporter link")
    clean_axis(ax)
    ax.text(
        0.01,
        -0.2,
        "Approximate Bayesian OLS; cross-sectional bulk tissue. Intervals are hypothesis-specification, not causal proof.",
        transform=ax.transAxes,
        fontsize=8,
        color=PALETTE["gray"],
    )
    save(fig, out_dir, "v11_h3_mediation_forest")
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


def generate_all(run_root: Path, out_dir: Path | None = None) -> int:
    configure_style()
    out_dir = out_dir or run_root / "figures" / "v11"
    ensure(out_dir)
    print(f"[v11-figures] output: {out_dir}")

    made = 0
    for fn in [
        fig_h2_primary_enrichment,
        fig_h2_composition_robustness,
        fig_h3_mediation_forest,
        fig_spatial_reference,
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
                "Claim boundary: H2 panels show DCT1-prior enrichment in whole-kidney OSD-462 phosphoproteomics, not DCT-isolated phosphoproteomics.",
                "Claim boundary: spatial panels are external IRI reference contextualization, not direct RR-10 spatial validation.",
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
