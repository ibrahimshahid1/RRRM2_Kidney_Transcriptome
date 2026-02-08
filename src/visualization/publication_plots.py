#!/usr/bin/env python
"""
plot_output.py - Generate publication-ready plots and tables for each pipeline phase.

Organizes outputs into a clean folder structure mirroring the export script's phase organization.

Usage:
    python scripts/plot_output.py --out_root plots --tag paper_figures
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Configure matplotlib
plt.rcParams.update({
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.figsize': (8, 6),
})
sns.set_style("whitegrid")

REPO_ROOT = Path(__file__).resolve().parents[1]


# Utilities
def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8")


def load_gene_symbols(results_dir: Path) -> Dict[str, str]:
    """Load Ensembl to gene symbol mapping."""
    # Try versioned path first, then fall back to legacy path
    symbol_file = results_dir / "phase6_uncertainty" / "ensembl_to_symbol.tsv"
    if not symbol_file.exists():
        symbol_file = results_dir.parent / "phase6_uncertainty" / "ensembl_to_symbol.tsv"
    if symbol_file.exists():
        df = pd.read_csv(symbol_file, sep="\t")
        if "ensembl" in df.columns and "symbol" in df.columns:
            return dict(zip(df["ensembl"], df["symbol"]))
        elif len(df.columns) >= 2:
            return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
    return {}


def safe_savefig(fig: plt.Figure, path: Path, **kwargs) -> None:
    """Save figure with error handling."""
    ensure_dir(path.parent)
    try:
        fig.savefig(path, bbox_inches="tight", **kwargs)
        print(f"  Saved: {path.name}")
    except Exception as e:
        print(f"  [ERROR] Failed to save {path.name}: {e}")
    finally:
        plt.close(fig)


# Phase 0: Deconvolution
def plot_phase0_deconvolution(data_dir: Path, out_dir: Path, results_dir: Path = None) -> Tuple[int, int]:
    """Generate deconvolution plots and tables."""
    phase_out = out_dir / "phase0_deconvolution"
    ensure_dir(phase_out)
    plots_created = 0
    tables_created = 0

    deconv_dir = data_dir / "processed" / "deconvolution"
    
    # Try to load proportions data
    prop_files = [
        deconv_dir / "music_segment_direct_proportions_CLR.csv",
        deconv_dir / "music_cluster_proportions.csv",
        deconv_dir / "music_group_proportions.csv",
    ]
    
    for prop_file in prop_files:
        if not prop_file.exists():
            continue
            
        try:
            df = pd.read_csv(prop_file)
            if df.empty:
                continue
                
            # Determine sample and cell type columns
            id_cols = [c for c in df.columns if c.lower() in ('sample', 'sample_id', 'id', 'row.names')]
            if not id_cols:
                id_cols = [df.columns[0]]
            
            value_cols = [c for c in df.columns if c not in id_cols]
            if not value_cols:
                continue
            
            # Melt for plotting
            df_long = df.melt(id_vars=id_cols, value_vars=value_cols, 
                              var_name="Cell_Type", value_name="Proportion")
            
            # Boxplot
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.boxplot(data=df_long, x="Cell_Type", y="Proportion", ax=ax, palette="Set2")
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
            ax.set_title(f"Cell Type Proportions ({prop_file.stem})")
            ax.set_xlabel("Cell Type")
            ax.set_ylabel("Proportion")
            safe_savefig(fig, phase_out / f"{prop_file.stem}_boxplot.png")
            plots_created += 1
            
            # Heatmap (samples x cell types)
            if len(df) > 1 and len(value_cols) > 1:
                fig, ax = plt.subplots(figsize=(12, max(8, len(df) * 0.3)))
                heatmap_data = df[value_cols].values
                
                # Truncate sample names for display
                sample_labels = df[id_cols[0]].astype(str).str[:20].tolist()
                
                sns.heatmap(heatmap_data, 
                           xticklabels=value_cols,
                           yticklabels=sample_labels if len(sample_labels) <= 50 else False,
                           cmap="YlOrRd", ax=ax, cbar_kws={"label": "Proportion"})
                ax.set_title(f"Cell Type Proportions Heatmap ({prop_file.stem})")
                ax.set_xlabel("Cell Type")
                ax.set_ylabel("Sample")
                safe_savefig(fig, phase_out / f"{prop_file.stem}_heatmap.png")
                plots_created += 1
            
            # Summary table
            summary = df[value_cols].describe().T
            summary.to_csv(phase_out / f"{prop_file.stem}_summary.tsv", sep="\t")
            tables_created += 1
            print(f"  Saved: {prop_file.stem}_summary.tsv")
            
        except Exception as e:
            print(f"  [WARN] Error processing {prop_file.name}: {e}")
    
    return plots_created, tables_created


# Phase 3: Rewiring
def plot_phase3_rewiring(data_dir: Path, out_dir: Path, gene_symbols: Dict[str, str], results_dir: Path = None) -> Tuple[int, int]:
    """Generate rewiring plots and tables."""
    phase_out = out_dir / "phase3_rewiring"
    ensure_dir(phase_out)
    plots_created = 0
    tables_created = 0

    # Use versioned results_dir if provided, else fall back to data_dir
    rewiring_dir = (results_dir or data_dir / "results") / "phase3_rewiring"
    if not rewiring_dir.exists():
        print(f"  [SKIP] Phase 3 rewiring directory not found: {rewiring_dir}")
        return 0, 0

    # Find all aggregated rewiring files
    agg_files = list(rewiring_dir.glob("*_rewiring_agg.tsv"))
    if not agg_files:
        print("  [SKIP] No aggregated rewiring files found")
        return 0, 0

    all_rewiring = []
    contrast_colors = {}
    color_palette = sns.color_palette("husl", len(agg_files))

    for i, agg_file in enumerate(agg_files):
        contrast = agg_file.stem.replace("_rewiring_agg", "")
        contrast_colors[contrast] = color_palette[i]
        
        try:
            df = pd.read_csv(agg_file, sep="\t")
            if "gene" not in df.columns or "rewiring_mean" not in df.columns:
                continue
            
            df["contrast"] = contrast
            df["symbol"] = df["gene"].map(gene_symbols).fillna(df["gene"])
            all_rewiring.append(df)
            
            # Top rewired genes table for this contrast
            top_genes = df.nlargest(50, "rewiring_mean")[["gene", "symbol", "rewiring_mean", "rewiring_std", "rank_std"]]
            top_genes.to_csv(phase_out / f"{contrast}_top50_rewired.tsv", sep="\t", index=False)
            tables_created += 1
            print(f"  Saved: {contrast}_top50_rewired.tsv")
            
        except Exception as e:
            print(f"  [WARN] Error loading {agg_file.name}: {e}")

    if not all_rewiring:
        return plots_created, tables_created

    combined = pd.concat(all_rewiring, ignore_index=True)

    # 1. Distribution of rewiring scores across contrasts
    fig, ax = plt.subplots(figsize=(10, 6))
    for contrast in combined["contrast"].unique():
        subset = combined[combined["contrast"] == contrast]
        ax.hist(subset["rewiring_mean"], bins=50, alpha=0.5, label=contrast, 
                color=contrast_colors.get(contrast))
    ax.set_xlabel("Rewiring Score (mean)")
    ax.set_ylabel("Number of Genes")
    ax.set_title("Distribution of Gene Rewiring Scores by Contrast")
    ax.legend(loc="upper right", fontsize=8)
    safe_savefig(fig, phase_out / "rewiring_distribution.png")
    plots_created += 1

    # 2. Violin plot comparing contrasts
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=combined, x="contrast", y="rewiring_mean", ax=ax, palette="husl")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.set_xlabel("Contrast")
    ax.set_ylabel("Rewiring Score (mean)")
    ax.set_title("Gene Rewiring Scores by Contrast")
    safe_savefig(fig, phase_out / "rewiring_violin.png")
    plots_created += 1

    # 3. Seed variance plot (rewiring_std distribution)
    if "rewiring_std" in combined.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(data=combined, x="contrast", y="rewiring_std", ax=ax, palette="husl")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_xlabel("Contrast")
        ax.set_ylabel("Rewiring Score StdDev (across seeds)")
        ax.set_title("Seed-to-Seed Variance in Rewiring Scores")
        safe_savefig(fig, phase_out / "rewiring_seed_variance.png")
        plots_created += 1

    # 4. Mean vs StdDev scatter (stability plot)
    fig, ax = plt.subplots(figsize=(10, 8))
    for contrast in combined["contrast"].unique():
        subset = combined[combined["contrast"] == contrast]
        ax.scatter(subset["rewiring_mean"], subset["rewiring_std"], 
                   alpha=0.3, s=10, label=contrast, color=contrast_colors.get(contrast))
    ax.set_xlabel("Rewiring Score (mean)")
    ax.set_ylabel("Rewiring Score (std)")
    ax.set_title("Rewiring Score Stability: Mean vs Variance")
    ax.legend(loc="upper right", fontsize=8)
    safe_savefig(fig, phase_out / "rewiring_mean_vs_std.png")
    plots_created += 1

    # Combined summary table
    summary = combined.groupby("contrast").agg({
        "rewiring_mean": ["mean", "std", "max"],
        "rewiring_std": ["mean", "max"],
        "gene": "count"
    }).round(4)
    summary.columns = ["_".join(col) for col in summary.columns]
    summary = summary.rename(columns={"gene_count": "n_genes"})
    summary.to_csv(phase_out / "rewiring_summary.tsv", sep="\t")
    tables_created += 1
    print(f"  Saved: rewiring_summary.tsv")

    return plots_created, tables_created


# Phase 5: Silent Shifters
def plot_phase5_silent_shifters(data_dir: Path, out_dir: Path, gene_symbols: Dict[str, str], results_dir: Path = None) -> Tuple[int, int]:
    """Generate silent shifter plots and tables."""
    phase_out = out_dir / "phase5_silent_shifters"
    ensure_dir(phase_out)
    plots_created = 0
    tables_created = 0

    # Use versioned results_dir if provided
    ss_dir = (results_dir or data_dir / "results") / "phase5_silent_shifters_strict"
    if not ss_dir.exists():
        print(f"  [SKIP] Phase 5 silent shifters directory not found: {ss_dir}")
        return 0, 0

    ss_files = list(ss_dir.glob("*_silent_shifters.tsv"))
    if not ss_files:
        print("  [SKIP] No silent shifter files found")
        return 0, 0

    all_ss = []
    counts_data = []

    for ss_file in ss_files:
        contrast = ss_file.stem.replace("_silent_shifters", "")
        
        try:
            df = pd.read_csv(ss_file, sep="\t")
            if df.empty:
                continue
                
            df["contrast"] = contrast
            df["symbol"] = df["gene"].map(gene_symbols).fillna(df["gene"])
            all_ss.append(df)
            
            # Count statistics
            n_total = len(df)
            n_high_rewiring = df["high_rewiring"].sum() if "high_rewiring" in df.columns else n_total
            n_low_de = df["low_DE"].sum() if "low_DE" in df.columns else n_total
            n_supported = df["supported"].sum() if "supported" in df.columns else 0
            
            counts_data.append({
                "contrast": contrast,
                "total_candidates": n_total,
                "high_rewiring": n_high_rewiring,
                "low_DE": n_low_de,
                "supported": n_supported
            })
            
            # Per-contrast table with symbols
            out_df = df[["gene", "symbol", "rewiring_mean", "rewiring_std"]].copy()
            if "logFC" in df.columns:
                out_df["logFC"] = df["logFC"]
            if "FDR" in df.columns:
                out_df["FDR"] = df["FDR"]
            if "q_BH" in df.columns:
                out_df["q_BH"] = df["q_BH"]
            out_df.to_csv(phase_out / f"{contrast}_silent_shifters.tsv", sep="\t", index=False)
            tables_created += 1
            print(f"  Saved: {contrast}_silent_shifters.tsv")
            
        except Exception as e:
            print(f"  [WARN] Error loading {ss_file.name}: {e}")

    if not counts_data:
        return plots_created, tables_created

    counts_df = pd.DataFrame(counts_data)

    # 1. Bar chart of silent shifter counts per contrast
    fig, ax = plt.subplots(figsize=(10, 6))
    x = range(len(counts_df))
    ax.bar(x, counts_df["total_candidates"], color="steelblue", alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(counts_df["contrast"], rotation=45, ha="right")
    ax.set_xlabel("Contrast")
    ax.set_ylabel("Number of Silent Shifter Candidates")
    ax.set_title("Silent Shifter Candidates by Contrast")
    for i, v in enumerate(counts_df["total_candidates"]):
        ax.text(i, v + 1, str(v), ha="center", fontsize=9)
    safe_savefig(fig, phase_out / "silent_shifter_counts.png")
    plots_created += 1

    # 2. Quadrant plot (rewiring vs q-value if available)
    if all_ss:
        combined = pd.concat(all_ss, ignore_index=True)
        
        if "q_BH" in combined.columns and "rewiring_mean" in combined.columns:
            fig, axes = plt.subplots(2, 2, figsize=(14, 12))
            axes = axes.flatten()
            
            for i, contrast in enumerate(combined["contrast"].unique()[:4]):
                subset = combined[combined["contrast"] == contrast].copy()
                subset["neg_log_q"] = -np.log10(subset["q_BH"].clip(lower=1e-10))
                
                ax = axes[i]
                ax.scatter(subset["neg_log_q"], subset["rewiring_mean"], 
                          alpha=0.5, s=20, c="steelblue")
                
                # Add quadrant lines
                ax.axhline(y=subset["rewiring_mean"].quantile(0.75), 
                          color="red", linestyle="--", alpha=0.5, label="High rewiring (75th)")
                ax.axvline(x=-np.log10(0.05), 
                          color="green", linestyle="--", alpha=0.5, label="q=0.05")
                
                ax.set_xlabel("-log10(q-value)")
                ax.set_ylabel("Rewiring Score")
                ax.set_title(contrast)
                ax.legend(fontsize=7, loc="upper left")
            
            plt.suptitle("Quadrant Plot: Rewiring vs Statistical Significance", fontsize=14)
            plt.tight_layout()
            safe_savefig(fig, phase_out / "quadrant_plot.png")
            plots_created += 1

        # 3. Rewiring distribution for silent shifters
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.violinplot(data=combined, x="contrast", y="rewiring_mean", ax=ax, palette="Set2")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_xlabel("Contrast")
        ax.set_ylabel("Rewiring Score")
        ax.set_title("Rewiring Scores of Silent Shifter Candidates")
        safe_savefig(fig, phase_out / "silent_shifter_rewiring.png")
        plots_created += 1

    # Summary counts table
    counts_df.to_csv(phase_out / "silent_shifter_summary.tsv", sep="\t", index=False)
    tables_created += 1
    print(f"  Saved: silent_shifter_summary.tsv")

    return plots_created, tables_created


# Phase 6: Uncertainty
def plot_phase6_uncertainty(data_dir: Path, out_dir: Path, gene_symbols: Dict[str, str], results_dir: Path = None) -> Tuple[int, int]:
    """Generate uncertainty (p-value, CI) plots and tables."""
    phase_out = out_dir / "phase6_uncertainty"
    ensure_dir(phase_out)
    plots_created = 0
    tables_created = 0

    # Use versioned results_dir if provided
    unc_dir = (results_dir or data_dir / "results") / "phase6_uncertainty"
    if not unc_dir.exists():
        print(f"  [SKIP] Phase 6 uncertainty directory not found: {unc_dir}")
        return 0, 0

    # P-value files
    pval_files = list(unc_dir.glob("*_perm_pvals.tsv"))
    ci_files = list(unc_dir.glob("*_bootstrap_ci.tsv"))

    all_pvals = []
    
    for pval_file in pval_files:
        contrast = pval_file.stem.replace("_perm_pvals", "")
        
        try:
            df = pd.read_csv(pval_file, sep="\t")
            if "gene" not in df.columns:
                continue
            
            df["contrast"] = contrast
            df["symbol"] = df["gene"].map(gene_symbols).fillna(df["gene"])
            all_pvals.append(df)
            
            # Significant genes table (q < 0.05)
            if "q_BH" in df.columns:
                sig_genes = df[df["q_BH"] < 0.05].sort_values("q_BH")
                if len(sig_genes) > 0:
                    out_cols = ["gene", "symbol"]
                    if "rewiring_abs_obs" in df.columns:
                        out_cols.append("rewiring_abs_obs")
                    if "p_perm" in df.columns:
                        out_cols.append("p_perm")
                    out_cols.append("q_BH")
                    sig_genes[out_cols].to_csv(
                        phase_out / f"{contrast}_significant_genes.tsv", sep="\t", index=False)
                    tables_created += 1
                    print(f"  Saved: {contrast}_significant_genes.tsv ({len(sig_genes)} genes)")
            
        except Exception as e:
            print(f"  [WARN] Error loading {pval_file.name}: {e}")

    if not all_pvals:
        return plots_created, tables_created

    combined = pd.concat(all_pvals, ignore_index=True)

    # 1. P-value histogram
    if "p_perm" in combined.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        for contrast in combined["contrast"].unique():
            subset = combined[combined["contrast"] == contrast]
            ax.hist(subset["p_perm"], bins=50, alpha=0.4, label=contrast)
        ax.axhline(y=len(combined) / len(combined["contrast"].unique()) / 50, 
                   color="red", linestyle="--", alpha=0.5, label="Uniform expectation")
        ax.set_xlabel("Permutation P-value")
        ax.set_ylabel("Number of Genes")
        ax.set_title("Distribution of Permutation P-values")
        ax.legend(fontsize=8)
        safe_savefig(fig, phase_out / "pvalue_histogram.png")
        plots_created += 1

    # 2. Q-Q plot
    if "p_perm" in combined.columns:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        axes = axes.flatten()
        
        for i, contrast in enumerate(combined["contrast"].unique()[:4]):
            subset = combined[combined["contrast"] == contrast]
            observed = np.sort(subset["p_perm"].dropna().values)
            expected = np.linspace(0, 1, len(observed))
            
            ax = axes[i]
            ax.scatter(expected, observed, alpha=0.3, s=5)
            ax.plot([0, 1], [0, 1], "r--", label="y=x (uniform)")
            ax.set_xlabel("Expected P-value")
            ax.set_ylabel("Observed P-value")
            ax.set_title(contrast)
            ax.legend(fontsize=8)
        
        plt.suptitle("Q-Q Plots: Observed vs Expected P-values", fontsize=14)
        plt.tight_layout()
        safe_savefig(fig, phase_out / "qq_plot.png")
        plots_created += 1

    # 3. -log10(q) distribution
    if "q_BH" in combined.columns:
        fig, ax = plt.subplots(figsize=(10, 6))
        combined["neg_log10_q"] = -np.log10(combined["q_BH"].clip(lower=1e-10))
        sns.boxplot(data=combined, x="contrast", y="neg_log10_q", ax=ax, palette="Set3")
        ax.axhline(y=-np.log10(0.05), color="red", linestyle="--", label="q=0.05")
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        ax.set_xlabel("Contrast")
        ax.set_ylabel("-log10(q-value)")
        ax.set_title("Distribution of Adjusted P-values by Contrast")
        ax.legend()
        safe_savefig(fig, phase_out / "qvalue_distribution.png")
        plots_created += 1

    # Summary statistics
    summary_rows = []
    for contrast in combined["contrast"].unique():
        subset = combined[combined["contrast"] == contrast]
        n_sig = (subset["q_BH"] < 0.05).sum() if "q_BH" in subset.columns else 0
        n_sig_01 = (subset["q_BH"] < 0.01).sum() if "q_BH" in subset.columns else 0
        summary_rows.append({
            "contrast": contrast,
            "n_genes": len(subset),
            "n_significant_q05": n_sig,
            "n_significant_q01": n_sig_01,
            "median_pvalue": subset["p_perm"].median() if "p_perm" in subset.columns else np.nan,
        })
    
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(phase_out / "uncertainty_summary.tsv", sep="\t", index=False)
    tables_created += 1
    print(f"  Saved: uncertainty_summary.tsv")

    return plots_created, tables_created


# Phase 6: Edge Regression (All 80 Samples)
def plot_phase6_regression(data_dir: Path, out_dir: Path, gene_symbols: Dict[str, str], results_dir: Path = None) -> Tuple[int, int]:
    """Generate plots and tables for regression-based network analysis (all 80 samples)."""
    phase_out = out_dir / "phase6_regression"
    ensure_dir(phase_out)
    plots_created = 0
    tables_created = 0

    # Use versioned results_dir if provided
    reg_dir = (results_dir or data_dir / "results") / "phase6_regression"
    if not reg_dir.exists():
        print(f"  [SKIP] Phase 6 regression directory not found: {reg_dir}")
        return 0, 0

    # Load flight effect and age×flight interaction results
    results = {}
    for effect_type in ["flight_effect", "age_flight_interaction"]:
        file_path = reg_dir / f"gene_{effect_type}.tsv"
        if file_path.exists():
            try:
                df = pd.read_csv(file_path, sep="\t")
                df["symbol"] = df["gene"].map(gene_symbols).fillna(df["gene"])
                results[effect_type] = df
            except Exception as e:
                print(f"  [WARN] Error loading {file_path.name}: {e}")

    if not results:
        print("  [SKIP] No regression result files found")
        return 0, 0

    # Process each effect type
    for effect_type, df in results.items():
        effect_label = effect_type.replace("_", " ").title()
        
        # 1. Export top significant genes (q < 0.05)
        sig_genes = df[df["q_BH"] < 0.05].copy()
        if len(sig_genes) > 0:
            # Top 100 or all significant
            top_sig = sig_genes.nsmallest(min(100, len(sig_genes)), "p_stouffer")
            out_cols = ["gene", "symbol", "n_edges", "median_beta", "median_t", "p_stouffer", "q_BH"]
            top_sig[out_cols].to_csv(phase_out / f"{effect_type}_top_significant.tsv", sep="\t", index=False)
            tables_created += 1
            print(f"  Saved: {effect_type}_top_significant.tsv ({len(top_sig)} genes)")
        
        # 2. Full results table
        df.to_csv(phase_out / f"{effect_type}_full.tsv", sep="\t", index=False)
        tables_created += 1
        print(f"  Saved: {effect_type}_full.tsv")
        
        # 3. Volcano plot (-log10(p) vs median_beta)
        fig, ax = plt.subplots(figsize=(10, 8))
        df["neg_log10_p"] = -np.log10(df["p_stouffer"].clip(lower=1e-300))
        
        # Color by significance
        colors = np.where(df["q_BH"] < 0.05, "red", "gray")
        alpha = np.where(df["q_BH"] < 0.05, 0.7, 0.3)
        
        ax.scatter(df["median_beta"], df["neg_log10_p"], c=colors, alpha=0.5, s=10)
        
        # Significance threshold line
        if (df["q_BH"] < 0.05).any():
            sig_threshold = df[df["q_BH"] < 0.05]["neg_log10_p"].min()
            ax.axhline(y=sig_threshold, color="red", linestyle="--", alpha=0.5, label=f"q=0.05")
        
        ax.axvline(x=0, color="gray", linestyle="-", alpha=0.3)
        ax.set_xlabel("Median Beta (effect direction)")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title(f"Volcano Plot: {effect_label}\n(Red = FDR < 0.05)")
        
        # Annotate top 10 genes
        top10 = df.nsmallest(10, "p_stouffer")
        for _, row in top10.iterrows():
            ax.annotate(row["symbol"], (row["median_beta"], row["neg_log10_p"]),
                       fontsize=7, alpha=0.8)
        
        safe_savefig(fig, phase_out / f"{effect_type}_volcano.png")
        plots_created += 1
        
        # 4. Distribution of p-values (histogram)
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.hist(df["p_stouffer"].dropna(), bins=50, color="steelblue", alpha=0.7, edgecolor="black")
        ax.axhline(y=len(df) / 50, color="red", linestyle="--", label="Uniform expectation")
        ax.set_xlabel("Stouffer P-value")
        ax.set_ylabel("Number of Genes")
        ax.set_title(f"P-value Distribution: {effect_label}")
        ax.legend()
        safe_savefig(fig, phase_out / f"{effect_type}_pvalue_hist.png")
        plots_created += 1
        
        # 5. Effect size distribution for significant genes
        if len(sig_genes) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.hist(sig_genes["median_beta"], bins=30, color="coral", alpha=0.7, edgecolor="black")
            ax.axvline(x=0, color="gray", linestyle="--")
            ax.set_xlabel("Median Beta (effect size)")
            ax.set_ylabel("Number of Significant Genes")
            ax.set_title(f"Effect Size Distribution: {effect_label} (q < 0.05)")
            safe_savefig(fig, phase_out / f"{effect_type}_effect_size_dist.png")
            plots_created += 1

    # 6. Comparison plot: Flight Effect vs Age×Flight Interaction
    if "flight_effect" in results and "age_flight_interaction" in results:
        flight_df = results["flight_effect"]
        interaction_df = results["age_flight_interaction"]
        
        # Merge on gene
        merged = flight_df[["gene", "symbol", "q_BH", "median_beta"]].merge(
            interaction_df[["gene", "q_BH", "median_beta"]],
            on="gene", suffixes=("_flight", "_interaction")
        )
        
        if len(merged) > 0:
            fig, ax = plt.subplots(figsize=(10, 10))
            
            # Color by significance category
            colors = []
            for _, row in merged.iterrows():
                if row["q_BH_flight"] < 0.05 and row["q_BH_interaction"] < 0.05:
                    colors.append("purple")  # Both
                elif row["q_BH_flight"] < 0.05:
                    colors.append("blue")    # Flight only
                elif row["q_BH_interaction"] < 0.05:
                    colors.append("orange")  # Interaction only
                else:
                    colors.append("gray")
            
            merged["neg_log_q_flight"] = -np.log10(merged["q_BH_flight"].clip(lower=1e-300))
            merged["neg_log_q_interaction"] = -np.log10(merged["q_BH_interaction"].clip(lower=1e-300))
            
            ax.scatter(merged["neg_log_q_flight"], merged["neg_log_q_interaction"],
                      c=colors, alpha=0.5, s=15)
            
            # Add significance thresholds
            ax.axhline(y=-np.log10(0.05), color="orange", linestyle="--", alpha=0.5)
            ax.axvline(x=-np.log10(0.05), color="blue", linestyle="--", alpha=0.5)
            
            ax.set_xlabel("-log10(q-value) Flight Effect")
            ax.set_ylabel("-log10(q-value) Age×Flight Interaction")
            ax.set_title("Comparison: Flight Effect vs Age×Flight Interaction\n" +
                        "Blue=Flight only, Orange=Interaction only, Purple=Both")
            
            # Annotate genes significant in both
            both_sig = merged[(merged["q_BH_flight"] < 0.05) & (merged["q_BH_interaction"] < 0.05)]
            for _, row in both_sig.head(15).iterrows():
                ax.annotate(row["symbol"], 
                           (row["neg_log_q_flight"], row["neg_log_q_interaction"]),
                           fontsize=6, alpha=0.8)
            
            safe_savefig(fig, phase_out / "flight_vs_interaction_comparison.png")
            plots_created += 1
            
            # Summary table of categories
            summary = {
                "category": ["Flight Effect Only", "Age×Flight Only", "Both", "Neither"],
                "count": [
                    ((merged["q_BH_flight"] < 0.05) & (merged["q_BH_interaction"] >= 0.05)).sum(),
                    ((merged["q_BH_flight"] >= 0.05) & (merged["q_BH_interaction"] < 0.05)).sum(),
                    ((merged["q_BH_flight"] < 0.05) & (merged["q_BH_interaction"] < 0.05)).sum(),
                    ((merged["q_BH_flight"] >= 0.05) & (merged["q_BH_interaction"] >= 0.05)).sum(),
                ]
            }
            pd.DataFrame(summary).to_csv(phase_out / "significance_category_counts.tsv", sep="\t", index=False)
            tables_created += 1
            print(f"  Saved: significance_category_counts.tsv")

    # 7. Summary statistics
    summary_rows = []
    for effect_type, df in results.items():
        summary_rows.append({
            "effect": effect_type,
            "n_genes": len(df),
            "n_significant_q05": (df["q_BH"] < 0.05).sum(),
            "n_significant_q01": (df["q_BH"] < 0.01).sum(),
            "n_significant_q001": (df["q_BH"] < 0.001).sum(),
            "median_p_stouffer": df["p_stouffer"].median(),
            "min_p_stouffer": df["p_stouffer"].min(),
        })
    
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(phase_out / "regression_summary.tsv", sep="\t", index=False)
    tables_created += 1
    print(f"  Saved: regression_summary.tsv")

    return plots_created, tables_created


# Summary README
def write_summary_readme(out_dir: Path, phase_stats: Dict[str, Tuple[int, int]]) -> None:
    """Write a summary README for the output folder."""
    total_plots = sum(s[0] for s in phase_stats.values())
    total_tables = sum(s[1] for s in phase_stats.values())
    
    lines = [
        "# RRRM2 Kidney Transcriptome — Plot Outputs",
        "",
        f"Generated: {dt.datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Summary",
        "",
        f"- **Total plots generated**: {total_plots}",
        f"- **Total tables generated**: {total_tables}",
        "",
        "## Phase Breakdown",
        "",
    ]
    
    for phase, (plots, tables) in phase_stats.items():
        lines.append(f"### {phase}")
        lines.append(f"- Plots: {plots}")
        lines.append(f"- Tables: {tables}")
        lines.append("")
    
    lines.extend([
        "## File Formats",
        "",
        "- **PNG**: High-resolution (300 DPI) figures for publication",
        "- **TSV**: Tab-separated tables for downstream analysis",
        "",
        "## Usage",
        "",
        "To regenerate these plots:",
        "```bash",
        "python scripts/plot_output.py --out_root plots --tag <your_tag>",
        "```",
    ])
    
    write_text(out_dir / "README.md", "\n".join(lines))


# Main
def main() -> None:
    os.chdir(REPO_ROOT)

    ap = argparse.ArgumentParser(
        description="Generate publication-ready plots and tables for each pipeline phase."
    )
    ap.add_argument("--out_root", default="plots", 
                    help="Where to create plot folders (relative to repo root).")
    ap.add_argument("--tag", default="figures", 
                    help="Tag added to output folder name.")
    ap.add_argument("--data_dir", default="data",
                    help="Data directory relative to repo root.")
    ap.add_argument("--results_dir", default=None,
                    help="Versioned results directory (overrides data_dir/results). Use with pipeline runner.")
    ap.add_argument("--out_dir", default=None,
                    help="Direct output directory (overrides out_root/timestamp_tag).")
    args = ap.parse_args()

    # Determine output directory
    if args.out_dir:
        out_dir = Path(args.out_dir).resolve()
    else:
        timestamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = (REPO_ROOT / args.out_root / f"{timestamp}_{args.tag}").resolve()
    ensure_dir(out_dir)
    
    data_dir = REPO_ROOT / args.data_dir
    
    # Resolve versioned results directory
    results_dir = Path(args.results_dir).resolve() if args.results_dir else None

    print("=" * 72)
    print("RRRM2 Kidney Transcriptome - Plot Generator")
    print("=" * 72)
    print(f"Output folder: {out_dir}")
    print(f"Data directory: {data_dir}")
    if results_dir:
        print(f"Results directory: {results_dir}")
    print()

    # Load gene symbols
    print("Loading gene symbol mappings...")
    gene_symbols = load_gene_symbols(results_dir or data_dir / "results")
    print(f"  Loaded {len(gene_symbols)} gene symbol mappings")
    print()

    phase_stats = {}

    # Phase 0: Deconvolution
    print("[Phase 0] Deconvolution plots...")
    plots, tables = plot_phase0_deconvolution(data_dir, out_dir, results_dir)
    phase_stats["phase0_deconvolution"] = (plots, tables)
    print()

    # Phase 3: Rewiring
    print("[Phase 3] Rewiring plots...")
    plots, tables = plot_phase3_rewiring(data_dir, out_dir, gene_symbols, results_dir)
    phase_stats["phase3_rewiring"] = (plots, tables)
    print()

    # Phase 5: Silent Shifters
    print("[Phase 5] Silent shifter plots...")
    plots, tables = plot_phase5_silent_shifters(data_dir, out_dir, gene_symbols, results_dir)
    phase_stats["phase5_silent_shifters"] = (plots, tables)
    print()

    # Phase 6: Uncertainty (permutation-based)
    print("[Phase 6] Uncertainty plots (permutation-based)...")
    plots, tables = plot_phase6_uncertainty(data_dir, out_dir, gene_symbols, results_dir)
    phase_stats["phase6_uncertainty"] = (plots, tables)
    print()

    # Phase 6: Regression (all 80 samples)
    print("[Phase 6] Regression plots (all 80 samples)...")
    plots, tables = plot_phase6_regression(data_dir, out_dir, gene_symbols, results_dir)
    phase_stats["phase6_regression"] = (plots, tables)
    print()

    # Write summary README
    write_summary_readme(out_dir, phase_stats)

    # Final summary
    total_plots = sum(s[0] for s in phase_stats.values())
    total_tables = sum(s[1] for s in phase_stats.values())

    print("=" * 72)
    print("PLOT GENERATION COMPLETE")
    print("=" * 72)
    print(f"Output folder: {out_dir}")
    print(f"Total plots: {total_plots}")
    print(f"Total tables: {total_tables}")
    print()
    
    for phase, (plots, tables) in phase_stats.items():
        print(f"  {phase}: {plots} plots, {tables} tables")


if __name__ == "__main__":
    main()