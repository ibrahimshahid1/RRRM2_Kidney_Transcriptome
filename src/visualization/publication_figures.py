#!/usr/bin/env python3
"""
publication_figures.py — Generate all key publication-ready figures.

Figures produced (all saved under <results_dir>/figures/publication/):
  1. rewiring_vs_logfc.png        — 4-panel scatter: rewiring vs |log2FC|, annotated
  2. edge_sum_perm_rank.png       — Per-contrast edge-sum rank plots with candidates labeled
  3. pathway_dotplot.png          — OR dot plot across all contrasts × pathways
  4. pathway_barplot.png          — Horizontal bar chart of top pathway hits
  5. age_comparison.png           — ISS-T Young vs ISS-T Old rewiring scatter
  6. arm_comparison.png           — ISS-T vs LAR recovery (persistence) plot
  7. pipeline_schematic.png       — Text-based pipeline overview
  8. retinol_subnetwork.png       — Retinol gene rewiring across contrasts
  9. result_tables.tsv            — Table 1 (sig genes) + Table 2 (sig pathways)

Usage (standalone):
    python -m src.visualization.publication_figures --results_dir data/results/run_XYZ

Called automatically by Phase 9 of the pipeline.
"""
from __future__ import annotations

import argparse
import gzip
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

warnings.filterwarnings("ignore", category=UserWarning)

# ── Style ──────────────────────────────────────────────────────────────────
PALETTE = {
    "ISS_T_YNG": "#E63946",
    "ISS_T_OLD": "#457B9D",
    "LAR_YNG":   "#2A9D8F",
    "LAR_OLD":   "#E9C46A",
}
LABEL_MAP = {
    "ISS_T_YNG": "ISS-T Young",
    "ISS_T_OLD": "ISS-T Old",
    "LAR_YNG":   "LAR Young",
    "LAR_OLD":   "LAR Old",
}
CONTRAST_REWIRING_NAMES = {
    "ISS_T_YNG": "ISS_T_YNG_FLT_minus_GC",
    "ISS_T_OLD": "ISS_T_OLD_FLT_minus_GC",
    "LAR_YNG":   "LAR_YNG_FLT_minus_GC",
    "LAR_OLD":   "LAR_OLD_FLT_minus_GC",
}
CONTRAST_DE_NAMES = {
    "ISS_T_YNG": "ISS_T_YNG_FLT_vs_GC_gene_DE.tsv",
    "ISS_T_OLD": "ISS_T_OLD_FLT_vs_GC_gene_DE.tsv",
    "LAR_YNG":   "LAR_YNG_FLT_vs_GC_gene_DE.tsv",
    "LAR_OLD":   "LAR_OLD_FLT_vs_GC_gene_DE.tsv",
}
RETINOL_GENES = [
    "Aldh1a7", "Aldh1a2", "Cyp26b1", "Cyp2a5", "Cyp4a12b",
    "Cyp4a31", "Retsat", "Ugt1a10", "Ugt1a2", "Ugt1a7c", "Bco2",
]

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
})
sns.set_style("whitegrid")


# ── Helpers ────────────────────────────────────────────────────────────────

def _ensure(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def _save(fig: plt.Figure, path: Path) -> None:
    _ensure(path.parent)
    fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"  [saved] {path.name}")


def _load_id_map(repo_root: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Return (eid→symbol, symbol_lower→eid) dicts."""
    p = repo_root / "data/processed/resources/id_map.tsv"
    if not p.exists():
        return {}, {}
    df = pd.read_csv(p, sep="\t", comment="#")
    col_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if "ensembl" in cl:
            col_map[c] = "eid"
        elif "symbol" in cl or "mgi" in cl:
            col_map[c] = "sym"
    df = df.rename(columns=col_map)
    e2s = dict(zip(df["eid"].astype(str), df["sym"].astype(str)))
    s2e = {}
    for eid, sym in e2s.items():
        s2e[sym.lower()] = eid
    return e2s, s2e


def _load_rewiring(results_dir: Path, contrast_key: str) -> Optional[pd.DataFrame]:
    name = CONTRAST_REWIRING_NAMES.get(contrast_key)
    if name is None:
        return None
    p = results_dir / f"phase3_rewiring/{name}_rewiring_agg.tsv"
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t")


def _load_de(repo_root: Path, contrast_key: str) -> Optional[pd.DataFrame]:
    fname = CONTRAST_DE_NAMES.get(contrast_key)
    if fname is None:
        return None
    p = repo_root / f"data/processed/gene_level_DE/{fname}"
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t")


def _load_gene_perm_inference(results_dir: Path, contrast_key: str) -> Optional[pd.DataFrame]:
    """Load Phase 6 edge-sum/candidate gene inference tables."""
    for suffix in ["GC", "GND"]:
        base = CONTRAST_REWIRING_NAMES.get(contrast_key, "").replace("GC", suffix)
        for fname in [
            f"{base}_perm_candidate_genes.tsv",
            f"{base}_perm_edge_sum_pvals.tsv",
            f"{base}_perm_pvals.tsv",
        ]:
            p = results_dir / f"phase6_uncertainty/{fname}"
            if p.exists():
                return pd.read_csv(p, sep="\t")
    return None


def _q_col(df: pd.DataFrame) -> Optional[str]:
    for c in ["q_BH_candidate", "q_BH_edge_sum", "q_BH", "q_BB_two_stage"]:
        if c in df.columns:
            return c
    return None


def _p_col(df: pd.DataFrame) -> Optional[str]:
    for c in ["p_perm_edge_sum", "p_perm"]:
        if c in df.columns:
            return c
    return None


def _edge_sum_col(df: pd.DataFrame) -> Optional[str]:
    for c in ["edge_sum_node_rewiring_obs", "rewiring_abs_obs"]:
        if c in df.columns:
            return c
    return None


def _load_perm_pvals(results_dir: Path, contrast_key: str) -> Optional[pd.DataFrame]:
    for suffix in ["GC", "GND"]:
        base = CONTRAST_REWIRING_NAMES.get(contrast_key, "").replace("GC", suffix)
        p = results_dir / f"phase6_uncertainty/{base}_perm_pvals.tsv"
        if p.exists():
            return pd.read_csv(p, sep="\t")
    return None


def _load_enrichment(results_dir: Path, contrast_key: str,
                     significant_only: bool = True) -> Optional[pd.DataFrame]:
    name = CONTRAST_REWIRING_NAMES.get(contrast_key)
    if name is None:
        return None
    suffix = "gene_set_enrichment_significant.tsv" if significant_only else "gene_set_enrichment.tsv"
    p = results_dir / f"phase7_grounding/{name}/{suffix}"
    if not p.exists():
        return None
    df = pd.read_csv(p, sep="\t")
    df["contrast"] = contrast_key
    return df


# ── Figure 1: Rewiring vs |log2FC| ────────────────────────────────────────

def fig_rewiring_vs_logfc(results_dir: Path, repo_root: Path, out_dir: Path) -> None:
    """4-panel scatter: rewiring vs |log2FC|, key genes annotated."""
    e2s, _ = _load_id_map(repo_root)
    contrasts = ["ISS_T_YNG", "ISS_T_OLD", "LAR_YNG", "LAR_OLD"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for ax, ckey in zip(axes, contrasts):
        rew = _load_rewiring(results_dir, ckey)
        de  = _load_de(repo_root, ckey)
        if rew is None or de is None:
            ax.text(0.5, 0.5, "Data not available", ha="center", va="center",
                    transform=ax.transAxes)
            ax.set_title(LABEL_MAP[ckey])
            continue

        df = rew.rename(columns={"rewiring_mean": "rewiring"}).merge(de, on="gene", how="inner")
        df["absFC"] = df["log2FC"].abs()
        df["sym"]   = df["gene"].map(e2s).fillna("???")

        # Spearman rho
        rho, pval = stats.spearmanr(df["rewiring"], df["absFC"])

        # Background points — color by rewiring decile
        df["decile"] = pd.qcut(df["rewiring"], 10, labels=False)
        colors = plt.cm.RdYlBu_r(df["decile"] / 9)

        ax.scatter(df["absFC"], df["rewiring"], c=colors, alpha=0.35, s=8, rasterized=True)

        # Load valid gene-level inference for annotation.
        foc = _load_gene_perm_inference(results_dir, ckey)
        annotated = set()
        qcol = _q_col(foc) if foc is not None else None
        if foc is not None and qcol:
            sig = foc[foc[qcol] < 0.05]["gene"]
            ann_df = df[df["gene"].isin(sig)].sort_values("rewiring", ascending=False)
            ax.scatter(ann_df["absFC"], ann_df["rewiring"],
                       c=PALETTE[ckey], s=60, zorder=5, edgecolors="black", linewidths=0.5)
            for _, row in ann_df.iterrows():
                ax.annotate(row["sym"],
                            xy=(row["absFC"], row["rewiring"]),
                            xytext=(4, 2), textcoords="offset points",
                            fontsize=7, color="#222222",
                            arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
                annotated.add(row["gene"])

        # Also annotate top-5 rewired if not already annotated
        for _, row in df.sort_values("rewiring", ascending=False).head(5).iterrows():
            if row["gene"] not in annotated:
                ax.annotate(row["sym"],
                            xy=(row["absFC"], row["rewiring"]),
                            xytext=(4, 2), textcoords="offset points",
                            fontsize=6.5, color="gray",
                            arrowprops=dict(arrowstyle="-", color="#cccccc", lw=0.5))

        ax.set_xlabel("|log₂ Fold Change|")
        ax.set_ylabel("Rewiring (cosine distance)")
        ax.set_title(f"{LABEL_MAP[ckey]}\nSpearman ρ = {rho:.3f}, p = {pval:.2e}")
        ax.tick_params(labelsize=8)

        # Correlation line (loess-approximated via numpy polyfit degree 1)
        z = np.polyfit(df["absFC"], df["rewiring"], 1)
        xline = np.linspace(df["absFC"].min(), df["absFC"].max(), 100)
        ax.plot(xline, np.polyval(z, xline), "k--", lw=1, alpha=0.6)

    fig.suptitle("Rewiring vs Differential Expression: Near-Zero Correlation\n"
                 "(Spaceflight network rewiring is independent of expression fold change)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    _save(fig, out_dir / "rewiring_vs_logfc.png")


# ── Figure 2: Edge-Sum Permutation Rank Plots ────────────────────────────

def fig_focused_perm_rank(results_dir: Path, repo_root: Path, out_dir: Path) -> None:
    """Per-contrast: rank vs edge-sum rewiring, significant candidates labeled."""
    e2s, _ = _load_id_map(repo_root)
    contrasts = ["ISS_T_YNG", "ISS_T_OLD", "LAR_YNG", "LAR_OLD"]

    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    axes = axes.flatten()

    for ax, ckey in zip(axes, contrasts):
        perm = _load_perm_pvals(results_dir, ckey)
        foc  = _load_gene_perm_inference(results_dir, ckey)

        if perm is None:
            ax.text(0.5, 0.5, "Permutation data not available",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(LABEL_MAP[ckey])
            continue

        rew_col = _edge_sum_col(perm)
        if rew_col is None:
            ax.text(0.5, 0.5, "Edge-sum column not available",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(LABEL_MAP[ckey])
            continue

        perm = perm.sort_values(rew_col, ascending=False).reset_index(drop=True)
        perm["rank"] = perm.index + 1
        perm["sym"]  = perm["gene"].map(e2s).fillna("???")

        # Color: significant (FDR<0.05), near-sig (FDR<0.10), rest
        sig_genes  = set()
        near_genes = set()
        qcol = _q_col(foc) if foc is not None else None
        if foc is not None and qcol:
            sig_genes  = set(foc[foc[qcol] < 0.05]["gene"])
            near_genes = set(foc[(foc[qcol] >= 0.05) & (foc[qcol] < 0.10)]["gene"])

        color_list = []
        for g in perm["gene"]:
            if g in sig_genes:
                color_list.append(PALETTE[ckey])
            elif g in near_genes:
                color_list.append("#FFAA00")
            else:
                color_list.append("#CCCCCC")

        sizes = [60 if g in sig_genes else 25 if g in near_genes else 6 for g in perm["gene"]]

        ax.scatter(perm["rank"], perm[rew_col],
                   c=color_list, s=sizes, alpha=0.7, rasterized=True)

        # Annotate significant genes
        ann_df = perm[perm["gene"].isin(sig_genes)].head(15)
        for _, row in ann_df.iterrows():
            ax.annotate(row["sym"],
                        xy=(row["rank"], row[rew_col]),
                        xytext=(6, 0), textcoords="offset points",
                        fontsize=7, color="#222222",
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

        n_sig  = len(sig_genes)
        n_near = len(near_genes)
        ax.set_xlabel("Gene rank (by rewiring magnitude)")
        ax.set_ylabel("Edge-sum node rewiring")
        ax.set_title(f"{LABEL_MAP[ckey]}\n"
                     f"{n_sig} FDR<0.05 ●  {n_near} FDR<0.10 ●")

        # Legend patches
        legend_handles = [
            mpatches.Patch(color=PALETTE[ckey], label=f"FDR < 0.05 (n={n_sig})"),
            mpatches.Patch(color="#FFAA00",     label=f"FDR < 0.10 (n={n_near})"),
            mpatches.Patch(color="#CCCCCC",     label="Not significant"),
        ]
        ax.legend(handles=legend_handles, fontsize=7, loc="upper right")

    fig.suptitle("Phase 6 Edge-Sum Node-Rewiring Permutation Test\n"
                 "(full-domain BH, candidate-only BH only when pre-registered)",
                 fontsize=12, y=1.02)
    plt.tight_layout()
    _save(fig, out_dir / "edge_sum_perm_rank.png")


# ── Figure 3: Pathway Enrichment Dot Plot ────────────────────────────────

def _build_enrichment_table(results_dir: Path) -> pd.DataFrame:
    """Merge significant enrichment across all contrasts."""
    frames = []
    for ckey in ["ISS_T_YNG", "ISS_T_OLD", "LAR_YNG", "LAR_OLD"]:
        df = _load_enrichment(results_dir, ckey, significant_only=True)
        if df is not None:
            df["contrast_label"] = LABEL_MAP[ckey]
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def fig_pathway_dotplot(results_dir: Path, out_dir: Path) -> None:
    """Bubble plot: pathway × contrast, bubble size = hits, color = -log10(q)."""
    df = _build_enrichment_table(results_dir)
    if df.empty:
        print("  [skip] No enrichment data found for pathway dotplot")
        return

    # Shorten pathway names
    def _shorten(name: str) -> str:
        name = name.split("::")[-1]
        if len(name) > 45:
            name = name[:43] + "…"
        return name

    df["short_name"]  = df["gene_set"].apply(_shorten)
    df["neg_log10_q"] = -np.log10(df["q_BH"].clip(1e-10))
    df["log_OR"]      = np.log2(df["odds_ratio"].clip(0.01))

    # Pivot for layout
    contrasts_ordered = ["ISS-T Young", "ISS-T Old", "LAR Young", "LAR Old"]
    pathway_order = (df.groupby("short_name")["neg_log10_q"]
                       .max().sort_values(ascending=False).index.tolist())

    fig_h = max(5, len(pathway_order) * 0.55)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    x_pos = {c: i for i, c in enumerate(contrasts_ordered)}
    y_pos = {p: i for i, p in enumerate(reversed(pathway_order))}

    for _, row in df.iterrows():
        x = x_pos.get(row["contrast_label"])
        y = y_pos.get(row["short_name"])
        if x is None or y is None:
            continue
        size  = max(20, min(300, row["hits_in_top_decile"] * 18))
        color = row["neg_log10_q"]
        ax.scatter(x, y, s=size, c=color, cmap="YlOrRd",
                   vmin=1, vmax=df["neg_log10_q"].quantile(0.95),
                   edgecolors="black", linewidths=0.4, zorder=3)

    ax.set_xticks(range(len(contrasts_ordered)))
    ax.set_xticklabels(contrasts_ordered, rotation=20, ha="right")
    ax.set_yticks(range(len(pathway_order)))
    ax.set_yticklabels(list(reversed(pathway_order)), fontsize=8)
    ax.set_xlim(-0.6, len(contrasts_ordered) - 0.4)
    ax.set_ylim(-0.6, len(pathway_order) - 0.4)
    ax.grid(True, alpha=0.3)

    sm = plt.cm.ScalarMappable(cmap="YlOrRd",
                                norm=plt.Normalize(1, df["neg_log10_q"].quantile(0.95)))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("−log₁₀(q)", fontsize=9)

    # Size legend
    for ns, label in [(3, "3 genes"), (6, "6 genes"), (10, "10 genes")]:
        ax.scatter([], [], s=ns * 18, c="gray", edgecolors="black",
                   linewidths=0.4, label=label, alpha=0.7)
    ax.legend(title="Hits in top decile", fontsize=8, title_fontsize=8,
              loc="lower right", framealpha=0.8)

    ax.set_title("Gene Set Enrichment: Rewired Pathways Across Contrasts\n"
                 "(Top-decile Fisher's exact test, BH-corrected)", fontsize=11)
    plt.tight_layout()
    _save(fig, out_dir / "pathway_dotplot.png")


def fig_pathway_barplot(results_dir: Path, out_dir: Path) -> None:
    """Horizontal bar chart: all significant pathways, colored by contrast."""
    df = _build_enrichment_table(results_dir)
    if df.empty:
        print("  [skip] No enrichment data found for pathway barplot")
        return

    def _shorten(name: str) -> str:
        name = name.split("::")[-1]
        if len(name) > 50:
            name = name[:48] + "…"
        return name

    df["short_name"]  = df["gene_set"].apply(_shorten)
    df["neg_log10_q"] = -np.log10(df["q_BH"].clip(1e-10))
    df["log_OR"]      = np.log2(df["odds_ratio"].clip(0.01))

    df = df.sort_values(["contrast_label", "neg_log10_q"], ascending=[True, False])

    bar_labels   = [f"{row['short_name']} ({row['contrast_label']})"
                    for _, row in df.iterrows()]
    bar_values   = df["log_OR"].values
    bar_colors   = [PALETTE.get(k, "#888888")
                    for k in df["contrast"].values]

    fig_h = max(6, len(bar_labels) * 0.38)
    fig, ax = plt.subplots(figsize=(10, fig_h))

    y_pos = np.arange(len(bar_labels))
    bars  = ax.barh(y_pos, bar_values, color=bar_colors, edgecolor="white",
                    linewidth=0.5, alpha=0.85)

    # Annotate with q value
    for i, (val, (_, row)) in enumerate(zip(bar_values, df.iterrows())):
        q_str = f"q={row['q_BH']:.3f}"
        ax.text(val + 0.05, i, q_str, va="center", fontsize=7, color="gray")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(bar_labels, fontsize=8)
    ax.set_xlabel("log₂(Odds Ratio)", fontsize=10)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title("Significant Enriched Pathways in Rewired Gene Sets\n"
                 "(Fisher's exact test, FDR < 0.05)", fontsize=11)

    legend_handles = [mpatches.Patch(color=v, label=LABEL_MAP[k])
                      for k, v in PALETTE.items()]
    ax.legend(handles=legend_handles, fontsize=8, loc="lower right")

    plt.tight_layout()
    _save(fig, out_dir / "pathway_barplot.png")


# ── Figure 4: Age Comparison ─────────────────────────────────────────────

def fig_age_comparison(results_dir: Path, repo_root: Path, out_dir: Path) -> None:
    """Scatter: ISS-T Young rewiring vs ISS-T Old rewiring, divergent genes labeled."""
    e2s, _ = _load_id_map(repo_root)

    rew_yng = _load_rewiring(results_dir, "ISS_T_YNG")
    rew_old = _load_rewiring(results_dir, "ISS_T_OLD")
    if rew_yng is None or rew_old is None:
        print("  [skip] Age comparison requires ISS-T YNG and OLD rewiring data")
        return

    df = rew_yng[["gene", "rewiring_mean"]].rename(columns={"rewiring_mean": "rew_yng"})
    df = df.merge(rew_old[["gene", "rewiring_mean"]].rename(columns={"rewiring_mean": "rew_old"}),
                  on="gene")
    df["sym"]   = df["gene"].map(e2s).fillna("???")
    df["delta"] = (df["rew_yng"] - df["rew_old"]).abs()

    # Load ISS_T interaction file
    int_file = results_dir / "phase5_derived/ISS_T_interaction.tsv"
    if int_file.exists():
        intf = pd.read_csv(int_file, sep="\t")
        df = df.merge(intf[["gene", "delta_interaction"]], on="gene", how="left")
    else:
        df["delta_interaction"] = df["delta"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Left panel: scatter
    ax = axes[0]
    rho, pval = stats.spearmanr(df["rew_yng"], df["rew_old"])

    sc = ax.scatter(df["rew_yng"], df["rew_old"],
                    c=df["delta_interaction"], cmap="coolwarm",
                    s=12, alpha=0.5, rasterized=True)
    fig.colorbar(sc, ax=ax, label="Age-dependent rewiring (|Δold − Δyng|)")

    # Diagonal = equal rewiring in both ages
    lim_max = max(df["rew_yng"].max(), df["rew_old"].max()) * 1.05
    ax.plot([0, lim_max], [0, lim_max], "k--", lw=1, alpha=0.5)

    # Annotate top age-divergent genes
    top_div = df.nlargest(12, "delta_interaction")
    for _, row in top_div.iterrows():
        ax.annotate(row["sym"],
                    xy=(row["rew_yng"], row["rew_old"]),
                    xytext=(4, 2), textcoords="offset points",
                    fontsize=7, color="#222222",
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

    ax.set_xlabel("Rewiring — Young (ISS-T FLT−GC)")
    ax.set_ylabel("Rewiring — Old (ISS-T FLT−GC)")
    ax.set_title(f"Age-Dependent Rewiring (ISS-T)\nSpearman ρ = {rho:.3f}, p = {pval:.1e}")

    # Right panel: horizontal lollipop of top age-interaction genes
    ax2 = axes[1]
    top20 = df.nlargest(20, "delta_interaction").sort_values("delta_interaction")
    colors2 = [PALETTE["ISS_T_YNG"] if y > o else PALETTE["ISS_T_OLD"]
               for y, o in zip(top20["rew_yng"], top20["rew_old"])]
    ax2.barh(top20["sym"], top20["rew_yng"], color=PALETTE["ISS_T_YNG"],
             alpha=0.6, label="Young")
    ax2.barh(top20["sym"], -top20["rew_old"], color=PALETTE["ISS_T_OLD"],
             alpha=0.6, label="Old")
    ax2.axvline(0, color="black", lw=0.8)
    ax2.set_xlabel("← Old   Rewiring   Young →")
    ax2.set_title("Top 20 Age-Divergent Genes (ISS-T)")
    ax2.legend(fontsize=8)
    ax2.tick_params(axis="y", labelsize=8)

    fig.suptitle("Age-Dependent Network Rewiring in Spaceflight Kidney", fontsize=12)
    plt.tight_layout()
    _save(fig, out_dir / "age_comparison.png")


# ── Figure 5: Arm Comparison (ISS-T vs LAR Recovery) ────────────────────

def fig_arm_comparison(results_dir: Path, repo_root: Path, out_dir: Path) -> None:
    """ISS-T vs LAR persistence: which rewiring events resolve after return?"""
    e2s, _ = _load_id_map(repo_root)

    rew_iss_yng = _load_rewiring(results_dir, "ISS_T_YNG")
    rew_lar_yng = _load_rewiring(results_dir, "LAR_YNG")
    rew_iss_old = _load_rewiring(results_dir, "ISS_T_OLD")
    rew_lar_old = _load_rewiring(results_dir, "LAR_OLD")

    if any(x is None for x in [rew_iss_yng, rew_lar_yng, rew_iss_old, rew_lar_old]):
        print("  [skip] Arm comparison requires all 4 contrast rewiring files")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, arm_iss, arm_lar, age_label, persist_file in [
        (axes[0], rew_iss_yng, rew_lar_yng, "Young",
         results_dir / "phase5_derived/ISS_minus_LAR_YNG_persistence.tsv"),
        (axes[1], rew_iss_old, rew_lar_old, "Old",
         results_dir / "phase5_derived/ISS_minus_LAR_OLD_persistence.tsv"),
    ]:
        df = (arm_iss[["gene", "rewiring_mean"]]
              .rename(columns={"rewiring_mean": "rew_iss"})
              .merge(arm_lar[["gene", "rewiring_mean"]]
                     .rename(columns={"rewiring_mean": "rew_lar"}), on="gene"))
        df["sym"]        = df["gene"].map(e2s).fillna("???")
        df["persistent"] = df[["rew_iss", "rew_lar"]].min(axis=1)
        df["resolved"]   = (df["rew_iss"] - df["rew_lar"]).clip(0)

        rho, _ = stats.spearmanr(df["rew_iss"], df["rew_lar"])

        # Load persistence delta if available
        if persist_file.exists():
            pf = pd.read_csv(persist_file, sep="\t")
            df = df.merge(pf[["gene", "ISS_minus_LAR_delta"]].rename(
                columns={"ISS_minus_LAR_delta": "iss_minus_lar"}), on="gene", how="left")
            df["iss_minus_lar"] = df["iss_minus_lar"].fillna(0)
        else:
            df["iss_minus_lar"] = df["rew_iss"] - df["rew_lar"]

        sc = ax.scatter(df["rew_iss"], df["rew_lar"],
                        c=df["iss_minus_lar"], cmap="RdBu_r",
                        vmin=-0.2, vmax=0.2,
                        s=12, alpha=0.5, rasterized=True)
        fig.colorbar(sc, ax=ax, label="ISS-T − LAR rewiring\n(+= ISS-specific, −= LAR-specific)")

        lim_max = max(df["rew_iss"].max(), df["rew_lar"].max()) * 1.05
        ax.plot([0, lim_max], [0, lim_max], "k--", lw=1, alpha=0.5)

        # Annotate top ISS-T-specific genes (high ISS-T, low LAR)
        top_iss = df.nlargest(8, "iss_minus_lar")
        for _, row in top_iss.iterrows():
            ax.annotate(row["sym"], xy=(row["rew_iss"], row["rew_lar"]),
                        xytext=(4, 2), textcoords="offset points",
                        fontsize=7, color=PALETTE["ISS_T_YNG"],
                        arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))

        ax.set_xlabel("ISS-T rewiring (FLT−GC)")
        ax.set_ylabel("LAR rewiring (FLT−GC)")
        ax.set_title(f"Recovery After Return — {age_label}\n"
                     f"Spearman ρ = {rho:.3f} (above diagonal = ISS-T-specific)")

    fig.suptitle("Spaceflight Network Rewiring: ISS-T vs Live Animal Return (LAR)\n"
                 "Genes above diagonal = rewiring persists only during spaceflight",
                 fontsize=11)
    plt.tight_layout()
    _save(fig, out_dir / "arm_comparison.png")


# ── Figure 6: Pipeline Schematic ─────────────────────────────────────────

def fig_pipeline_schematic(out_dir: Path) -> None:
    """Clean text-flow pipeline diagram."""
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 7)
    ax.axis("off")

    phases = [
        (0.5, 6.0,  "Phase 0\nDeconvolution\n(MuSiC + TMS)",        "#AED6F1"),
        (3.0, 6.0,  "Phase 1\nResidualisation\n(DESeq2-VST + SVA)",  "#A9DFBF"),
        (5.5, 6.0,  "Phase 2\nNetwork Construction\n(LIONESS + Regression)", "#FAD7A0"),
        (8.0, 6.0,  "Phase 3\nEmbeddings\n(node2vec + Procrustes)",  "#F9E79F"),
        (10.5,6.0,  "Phase 5\nSilent Shifters\n+ Interaction Metrics","#D2B4DE"),
        (0.5, 3.0,  "Phase 6\nEdge-Sum Permutation\n+ Candidate/Hierarchical FDR", "#F1948A"),
        (3.5, 3.0,  "Phase 6\nFull Regression\n(n=80, factorial)",   "#F1948A"),
        (6.5, 3.0,  "Phase 7\nPathway Enrichment\n(KEGG/Reactome/Hallmark)", "#ABEBC6"),
        (9.5, 3.0,  "Phase 9\nFigure Generation\n(this output)",     "#D5DBDB"),
    ]
    arrows_horizontal = [
        (0.5, 3.0, 6.0), (3.0, 5.5, 6.0), (5.5, 8.0, 6.0), (8.0, 10.5, 6.0),
    ]
    arrows_down = [
        (10.5, 6.0, 3.0),
        (0.5,  6.0, 3.0),
    ]

    box_w, box_h = 2.3, 1.1
    for (x, y, label, color) in phases:
        rect = mpatches.FancyBboxPatch(
            (x - box_w / 2, y - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.08", facecolor=color,
            edgecolor="#555555", linewidth=1.2, zorder=3
        )
        ax.add_patch(rect)
        ax.text(x, y, label, ha="center", va="center",
                fontsize=8, fontweight="bold", zorder=4)

    # Horizontal arrows top row
    for x1, x2, y in arrows_horizontal:
        ax.annotate("", xy=(x2 - box_w / 2, y), xytext=(x1 + box_w / 2, y),
                    arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

    # Vertical arrows (Phase 5 → Phase 6, Phase 0 → Phase 6 context)
    ax.annotate("", xy=(0.5, 3.0 + box_h / 2), xytext=(0.5, 6.0 - box_h / 2),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))
    ax.annotate("", xy=(9.5, 3.0 + box_h / 2), xytext=(10.5, 6.0 - box_h / 2),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

    # Phase 6a → Phase 7 horizontal
    ax.annotate("", xy=(6.5 - box_w / 2, 3.0), xytext=(3.5 + box_w / 2, 3.0),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))
    ax.annotate("", xy=(9.5 - box_w / 2, 3.0), xytext=(6.5 + box_w / 2, 3.0),
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.5))

    # Row label annotations
    ax.text(0.1, 6.0, "Data\nPreperation", ha="left", va="center",
            fontsize=9, color="gray", rotation=90)
    ax.text(0.1, 3.0, "Inference\n& Grounding", ha="left", va="center",
            fontsize=9, color="gray", rotation=90)

    # Input/output labels
    ax.text(7.0, 0.8,
            "Input: 80 bulk RNA-seq kidney samples (OSD-771 RRRM-2)\n"
            "Design: 2×2 factorial (Age×Arm) × 4 conditions, n=5 per cell\n"
            "Output: Rewiring scores, permutation FDR, pathway enrichment",
            ha="center", va="center", fontsize=9,
            bbox=dict(facecolor="#F8F9FA", edgecolor="#CCCCCC",
                      boxstyle="round,pad=0.4"))

    ax.set_title("RRRM-2 Kidney Co-Expression Network Rewiring Pipeline",
                 fontsize=13, fontweight="bold", pad=10)
    plt.tight_layout()
    _save(fig, out_dir / "pipeline_schematic.png")


# ── Figure 7: Retinol Subnetwork ─────────────────────────────────────────

def fig_retinol_subnetwork(results_dir: Path, repo_root: Path, out_dir: Path) -> None:
    """Heatmap + barplot of retinol pathway gene rewiring across all contrasts."""
    e2s, s2e = _load_id_map(repo_root)

    # Collect rewiring for retinol genes across all contrasts
    records = []
    for ckey in ["ISS_T_YNG", "ISS_T_OLD", "LAR_YNG", "LAR_OLD"]:
        rew = _load_rewiring(results_dir, ckey)
        perm = _load_perm_pvals(results_dir, ckey)
        foc  = _load_gene_perm_inference(results_dir, ckey)
        if rew is None:
            continue
        rew["sym"] = rew["gene"].map(e2s).fillna("???")

        qcol = _q_col(foc) if foc is not None else None
        sig_set  = set(foc["gene"]) if foc is not None and qcol else set()
        sig05    = set(foc[foc[qcol] < 0.05]["gene"]) if foc is not None and qcol else set()
        sig10    = set(foc[foc[qcol] < 0.10]["gene"]) if foc is not None and qcol else set()

        # p_perm lookup
        p_lookup = {}
        if perm is not None:
            pcol = _p_col(perm)
            p_lookup = dict(zip(perm["gene"], perm[pcol])) if pcol else {}

        for sym in RETINOL_GENES:
            eid = s2e.get(sym.lower())
            row = rew[rew["gene"] == eid] if eid else pd.DataFrame()
            if row.empty:
                # Try by symbol match
                row = rew[rew["sym"].str.lower() == sym.lower()]
            if row.empty:
                continue
            gene = row["gene"].iloc[0]
            records.append({
                "sym": sym,
                "contrast": ckey,
                "contrast_label": LABEL_MAP[ckey],
                "rewiring": row["rewiring_mean"].iloc[0],
                "p_perm": p_lookup.get(gene, 1.0),
                "sig05": gene in sig05,
                "sig10": gene in sig10,
            })

    if not records:
        print("  [skip] No retinol gene data found")
        return

    df = pd.DataFrame(records)
    pivot = df.pivot_table(index="sym", columns="contrast_label",
                           values="rewiring", aggfunc="mean")

    # Reorder columns
    col_order = [c for c in ["ISS-T Young", "ISS-T Old", "LAR Young", "LAR Old"]
                 if c in pivot.columns]
    pivot = pivot[col_order]
    # Sort rows by mean rewiring across contrasts
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Left: heatmap
    ax_h = axes[0]
    sns.heatmap(pivot, ax=ax_h, cmap="YlOrRd", annot=True, fmt=".3f",
                linewidths=0.5, linecolor="white",
                cbar_kws={"label": "Rewiring (cosine distance)"},
                annot_kws={"size": 8})
    ax_h.set_title("Retinol Metabolism Gene Rewiring\nAcross All Contrasts", fontsize=10)
    ax_h.set_ylabel("Gene")
    ax_h.set_xlabel("")
    ax_h.tick_params(axis="x", rotation=20)

    # Right: grouped barplot showing retinol genes vs genome background per contrast
    ax_b = axes[1]
    contrast_order = ["ISS_T_YNG", "ISS_T_OLD", "LAR_YNG", "LAR_OLD"]
    retinol_means  = []
    background_means = []
    retinol_errs   = []
    background_errs = []

    for ckey in contrast_order:
        rew = _load_rewiring(results_dir, ckey)
        if rew is None:
            retinol_means.append(0)
            background_means.append(0)
            retinol_errs.append(0)
            background_errs.append(0)
            continue
        rew["sym"] = rew["gene"].map(e2s).fillna("???")
        retinol_rows = rew[rew["sym"].isin(RETINOL_GENES)]
        rest = rew[~rew["sym"].isin(RETINOL_GENES)]
        retinol_means.append(retinol_rows["rewiring_mean"].mean() if len(retinol_rows) else 0)
        retinol_errs.append(retinol_rows["rewiring_mean"].sem() if len(retinol_rows) else 0)
        background_means.append(rest["rewiring_mean"].mean())
        background_errs.append(rest["rewiring_mean"].sem())

    x   = np.arange(len(contrast_order))
    w   = 0.35
    labels = [LABEL_MAP[c] for c in contrast_order]

    bars1 = ax_b.bar(x - w/2, retinol_means, w, yerr=retinol_errs,
                     color="#E07B39", label="Retinol pathway genes", capsize=4)
    bars2 = ax_b.bar(x + w/2, background_means, w, yerr=background_errs,
                     color="#AAAAAA", label="Genome background", capsize=4, alpha=0.7)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(labels, rotation=15, ha="right")
    ax_b.set_ylabel("Mean rewiring (cosine distance)")
    ax_b.set_title("Retinol Genes vs Background Rewiring")
    ax_b.legend(fontsize=8)

    fig.suptitle("Retinol / Retinoic Acid Metabolism — Convergent Multi-Scale Signal",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    _save(fig, out_dir / "retinol_subnetwork.png")


# ── Table Generation ──────────────────────────────────────────────────────

def _build_result_tables(results_dir: Path, repo_root: Path, out_dir: Path) -> None:
    """Write Table 1 (sig genes) and Table 2 (sig pathways) as TSVs."""
    e2s, _ = _load_id_map(repo_root)

    # Table 1: valid gene-level edge-sum/candidate candidates
    rows = []
    for ckey in ["ISS_T_YNG", "ISS_T_OLD", "LAR_YNG", "LAR_OLD"]:
        foc = _load_gene_perm_inference(results_dir, ckey)
        rew = _load_rewiring(results_dir, ckey)
        de  = _load_de(repo_root, ckey)
        if foc is None or rew is None:
            continue
        qcol = _q_col(foc)
        if qcol is None:
            continue
        sig = foc[foc[qcol] < 0.05].copy()
        if sig.empty:
            continue
        sig["sym"] = sig["gene"].map(e2s).fillna("???")
        if de is not None:
            sig = sig.merge(de[["gene", "log2FC", "FDR"]].rename(
                columns={"FDR": "DE_FDR"}), on="gene", how="left")
        sig["contrast"]       = LABEL_MAP[ckey]
        rew_col = _edge_sum_col(sig)
        if rew_col:
            sig["rewiring_rank"] = sig[rew_col].rank(ascending=False).astype(int)
        pcol = _p_col(sig)
        sig["q_value"] = sig[qcol]
        sig["p_value"] = sig[pcol] if pcol else np.nan
        sig["q_value_source"] = qcol
        rows.append(sig)

    if rows:
        t1 = pd.concat(rows, ignore_index=True)
        cols = ["contrast", "sym", "gene", "edge_sum_node_rewiring_obs", "rewiring_abs_obs",
                "rewiring_rank", "p_value", "q_value", "q_value_source"]
        if "log2FC" in t1.columns:
            cols += ["log2FC", "DE_FDR"]
        t1 = t1[[c for c in cols if c in t1.columns]]
        t1 = t1.sort_values(["contrast", "q_value"])
        t1.to_csv(out_dir / "table1_edge_sum_candidates.tsv", sep="\t", index=False)
        print(f"  [saved] table1_edge_sum_candidates.tsv ({len(t1)} rows)")

    # Table 2: Significant enrichment pathways
    enr = _build_enrichment_table(results_dir)
    if not enr.empty:
        def _shorten(name: str) -> str:
            return name.split("::")[-1]
        enr["pathway_short"] = enr["gene_set"].apply(_shorten)
        cols2 = ["contrast_label", "pathway_short", "gene_set", "library",
                 "set_size_in_universe", "hits_in_top_decile",
                 "odds_ratio", "p", "q_BH", "hit_symbols"]
        enr2 = enr[[c for c in cols2 if c in enr.columns]]
        enr2 = enr2.sort_values(["contrast_label", "q_BH"])
        enr2.to_csv(out_dir / "table2_significant_pathways.tsv", sep="\t", index=False)
        print(f"  [saved] table2_significant_pathways.tsv ({len(enr2)} rows)")


# ── Main ──────────────────────────────────────────────────────────────────

def generate_all_figures(results_dir: Path, repo_root: Path) -> None:
    """Generate all publication figures for a given run."""
    out_dir = results_dir / "figures" / "publication"
    _ensure(out_dir)

    print(f"  Output: {out_dir}")

    print("\n  [1/8] Rewiring vs |log2FC| scatter...")
    try:
        fig_rewiring_vs_logfc(results_dir, repo_root, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [2/8] Edge-sum permutation rank plots...")
    try:
        fig_focused_perm_rank(results_dir, repo_root, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [3/8] Pathway enrichment dot plot...")
    try:
        fig_pathway_dotplot(results_dir, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [4/8] Pathway enrichment bar chart...")
    try:
        fig_pathway_barplot(results_dir, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [5/8] Age comparison figure...")
    try:
        fig_age_comparison(results_dir, repo_root, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [6/8] Arm comparison (ISS-T vs LAR recovery)...")
    try:
        fig_arm_comparison(results_dir, repo_root, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [7/8] Pipeline schematic...")
    try:
        fig_pipeline_schematic(out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  [8/8] Retinol subnetwork...")
    try:
        fig_retinol_subnetwork(results_dir, repo_root, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print("\n  Building result tables...")
    try:
        _build_result_tables(results_dir, repo_root, out_dir)
    except Exception as e:
        print(f"    [ERROR] {e}")

    print(f"\n  Done. All outputs in: {out_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate publication figures")
    parser.add_argument("--results_dir", required=True,
                        help="Path to a run's results directory (e.g. data/results/run_XYZ)")
    parser.add_argument("--repo_root", default=None,
                        help="Repository root (default: auto-detected from file location)")
    args = parser.parse_args()

    results_dir = Path(args.results_dir).resolve()
    repo_root   = (Path(args.repo_root).resolve() if args.repo_root
                   else Path(__file__).resolve().parents[2])

    print("=" * 70)
    print("  Publication Figure Generator")
    print(f"  Run:  {results_dir.name}")
    print(f"  Root: {repo_root}")
    print("=" * 70)

    generate_all_figures(results_dir, repo_root)


if __name__ == "__main__":
    main()
