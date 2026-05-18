# scripts/plot_skeleton_diagnostics.py
"""
Phase 2 Skeleton Diagnostics and Visualization

Generates publication-quality figures for the network skeleton:
1. Partial correlation distribution (edge weight histogram)
2. Degree distribution (edges per gene)
3. Top hub genes table
4. Network graph visualization (top genes by degree)
5. LIONESS z-score statistics

Usage:
    python scripts/plot_skeleton_diagnostics.py
"""
from __future__ import annotations

import argparse
from pathlib import Path
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter

# Try to import networkx for graph visualization
try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False
    print("[WARN] networkx not installed. Graph visualization will be skipped.")

REPO_ROOT = Path(__file__).resolve().parents[1]


def load_skeleton(phase2_dir: Path) -> tuple[list[str], pd.DataFrame, np.ndarray, np.ndarray]:
    """Load skeleton data."""
    genes = [g.strip() for g in (phase2_dir / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    edges = pd.read_csv(phase2_dir / "skeleton_edges.tsv", sep="\t")
    edge_i = np.load(phase2_dir / "edge_i.npy")
    edge_j = np.load(phase2_dir / "edge_j.npy")
    return genes, edges, edge_i, edge_j


def plot_pcorr_distribution(edges: pd.DataFrame, outdir: Path):
    """Plot partial correlation distribution."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Raw partial correlations
    ax1 = axes[0]
    pcorr = edges["pcorr"].values
    ax1.hist(pcorr, bins=80, color="#3498db", edgecolor="white", alpha=0.8)
    ax1.axvline(0, color="red", linestyle="--", linewidth=1.5, label="Zero")
    ax1.axvline(np.median(pcorr), color="orange", linestyle="-", linewidth=1.5, 
                label=f"Median: {np.median(pcorr):.3f}")
    ax1.set_xlabel("Partial Correlation", fontsize=12)
    ax1.set_ylabel("Edge Count", fontsize=12)
    ax1.set_title("Partial Correlation Distribution (All Edges)", fontsize=13, fontweight="bold")
    ax1.legend(loc="upper right")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # Absolute partial correlations
    ax2 = axes[1]
    abs_pcorr = edges["abs_pcorr"].values
    ax2.hist(abs_pcorr, bins=60, color="#2ecc71", edgecolor="white", alpha=0.8)
    ax2.axvline(np.median(abs_pcorr), color="orange", linestyle="-", linewidth=1.5,
                label=f"Median: {np.median(abs_pcorr):.3f}")
    ax2.axvline(np.percentile(abs_pcorr, 95), color="red", linestyle="--", linewidth=1.5,
                label=f"95th: {np.percentile(abs_pcorr, 95):.3f}")
    ax2.set_xlabel("|Partial Correlation|", fontsize=12)
    ax2.set_ylabel("Edge Count", fontsize=12)
    ax2.set_title("Absolute Partial Correlation Distribution", fontsize=13, fontweight="bold")
    ax2.legend(loc="upper right")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    plt.tight_layout()
    fig.savefig(outdir / "skeleton_pcorr_distribution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: skeleton_pcorr_distribution.png")


def plot_degree_distribution(genes: list[str], edge_i: np.ndarray, edge_j: np.ndarray, outdir: Path):
    """Plot degree distribution."""
    # Count degree per gene
    degree = Counter()
    for i, j in zip(edge_i, edge_j):
        degree[i] += 1
        degree[j] += 1
    
    degrees = [degree.get(i, 0) for i in range(len(genes))]
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    ax1 = axes[0]
    ax1.hist(degrees, bins=50, color="#9b59b6", edgecolor="white", alpha=0.8)
    ax1.axvline(np.mean(degrees), color="red", linestyle="-", linewidth=1.5,
                label=f"Mean: {np.mean(degrees):.1f}")
    ax1.axvline(np.median(degrees), color="orange", linestyle="--", linewidth=1.5,
                label=f"Median: {np.median(degrees):.0f}")
    ax1.set_xlabel("Degree (edges per gene)", fontsize=12)
    ax1.set_ylabel("Gene Count", fontsize=12)
    ax1.set_title("Degree Distribution", fontsize=13, fontweight="bold")
    ax1.legend(loc="upper right")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # Log-log plot for scale-free assessment
    ax2 = axes[1]
    degree_counts = Counter(degrees)
    x = sorted(degree_counts.keys())
    y = [degree_counts[d] for d in x]
    ax2.scatter(x, y, color="#e74c3c", alpha=0.7, s=30)
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("Degree (log scale)", fontsize=12)
    ax2.set_ylabel("Frequency (log scale)", fontsize=12)
    ax2.set_title("Degree Distribution (Log-Log)", fontsize=13, fontweight="bold")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    plt.tight_layout()
    fig.savefig(outdir / "skeleton_degree_distribution.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: skeleton_degree_distribution.png")
    
    return degrees


def save_hub_genes(genes: list[str], degrees: list[int], outdir: Path, top_n: int = 50):
    """Save top hub genes table."""
    df = pd.DataFrame({
        "gene": genes,
        "degree": degrees
    }).sort_values("degree", ascending=False).reset_index(drop=True)
    
    df["rank"] = range(1, len(df) + 1)
    df = df[["rank", "gene", "degree"]]
    
    # Save full table
    df.to_csv(outdir / "skeleton_gene_degrees.tsv", sep="\t", index=False)
    print(f"  Saved: skeleton_gene_degrees.tsv (all {len(df)} genes)")
    
    # Save top N
    top = df.head(top_n)
    top.to_csv(outdir / f"skeleton_top{top_n}_hub_genes.tsv", sep="\t", index=False)
    print(f"  Saved: skeleton_top{top_n}_hub_genes.tsv")
    
    return df


def plot_network_subgraph(genes: list[str], edges: pd.DataFrame, degrees: list[int], 
                          outdir: Path, top_n: int = 100):
    """Plot network graph for top hub genes."""
    if not HAS_NETWORKX:
        print("  [SKIP] Network graph requires networkx")
        return
    
    # Get top genes by degree
    gene_degree = pd.DataFrame({"gene": genes, "degree": degrees, "idx": range(len(genes))})
    top_genes = gene_degree.nlargest(top_n, "degree")
    top_idx = set(top_genes["idx"].values)
    top_gene_set = set(top_genes["gene"].values)
    
    # Filter edges to top genes only
    mask = edges["gene_i"].isin(top_gene_set) & edges["gene_j"].isin(top_gene_set)
    sub_edges = edges[mask]
    
    # Build graph
    G = nx.Graph()
    for _, row in top_genes.iterrows():
        G.add_node(row["gene"], degree=row["degree"])
    
    for _, row in sub_edges.iterrows():
        G.add_edge(row["gene_i"], row["gene_j"], weight=abs(row["pcorr"]))
    
    # Layout
    pos = nx.spring_layout(G, k=2/np.sqrt(len(G.nodes())), iterations=50, seed=42)
    
    # Node sizes by degree
    node_sizes = [G.nodes[n].get("degree", 10) * 3 for n in G.nodes()]
    
    # Edge widths by weight
    edge_weights = [G[u][v].get("weight", 0.1) * 2 for u, v in G.edges()]
    
    fig, ax = plt.subplots(figsize=(14, 14))
    
    # Draw
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=edge_weights, edge_color="#888888", ax=ax)
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color="#3498db", 
                           alpha=0.8, edgecolors="white", linewidths=1, ax=ax)
    
    # Label only top 20 hubs
    top20 = top_genes.head(20)["gene"].tolist()
    labels = {n: n for n in G.nodes() if n in top20}
    nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight="bold", ax=ax)
    
    ax.set_title(f"Network Graph: Top {top_n} Hub Genes\n({len(sub_edges)} edges shown)", 
                 fontsize=14, fontweight="bold")
    ax.axis("off")
    
    plt.tight_layout()
    fig.savefig(outdir / f"skeleton_network_top{top_n}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: skeleton_network_top{top_n}.png")


def plot_lioness_statistics(phase2_dir: Path, outdir: Path):
    """Plot LIONESS edge-weight statistics."""
    lioness_path = phase2_dir / "lioness_edges.npy"
    if not lioness_path.exists():
        lioness_path = phase2_dir / "lioness_z_edges.npy"
    if not lioness_path.exists():
        print("  [SKIP] No LIONESS edge weights found")
        return
    
    Z = np.load(lioness_path)
    N, E = Z.shape
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Overall edge-weight distribution
    ax1 = axes[0, 0]
    z_flat = Z.flatten()
    sample_idx = np.random.choice(len(z_flat), size=min(500000, len(z_flat)), replace=False)
    z_sample = z_flat[sample_idx]
    ax1.hist(z_sample, bins=100, color="#3498db", edgecolor="white", alpha=0.8)
    ax1.axvline(0, color="red", linestyle="--", linewidth=1)
    ax1.set_xlabel("LIONESS edge weight", fontsize=11)
    ax1.set_ylabel("Count", fontsize=11)
    ax1.set_title(f"LIONESS Edge-Weight Distribution (sampled)", fontsize=12, fontweight="bold")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    
    # 2. Per-sample mean z
    ax2 = axes[0, 1]
    sample_means = Z.mean(axis=1)
    ax2.bar(range(N), sample_means, color="#2ecc71", alpha=0.8)
    ax2.axhline(0, color="red", linestyle="--", linewidth=1)
    ax2.set_xlabel("Sample Index", fontsize=11)
    ax2.set_ylabel("Mean z-score", fontsize=11)
    ax2.set_title("Per-Sample Mean Edge Weight", fontsize=12, fontweight="bold")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    
    # 3. Per-sample variance
    ax3 = axes[1, 0]
    sample_vars = Z.var(axis=1)
    ax3.bar(range(N), sample_vars, color="#9b59b6", alpha=0.8)
    ax3.set_xlabel("Sample Index", fontsize=11)
    ax3.set_ylabel("Variance", fontsize=11)
    ax3.set_title("Per-Sample Edge Weight Variance", fontsize=12, fontweight="bold")
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    
    # 4. Edge-level variance
    ax4 = axes[1, 1]
    edge_vars = Z.var(axis=0)
    ax4.hist(edge_vars, bins=80, color="#e74c3c", edgecolor="white", alpha=0.8)
    ax4.axvline(np.median(edge_vars), color="orange", linestyle="-", linewidth=1.5,
                label=f"Median: {np.median(edge_vars):.3f}")
    ax4.set_xlabel("Edge Variance (across samples)", fontsize=11)
    ax4.set_ylabel("Edge Count", fontsize=11)
    ax4.set_title("Per-Edge Variance Distribution", fontsize=12, fontweight="bold")
    ax4.legend()
    ax4.spines["top"].set_visible(False)
    ax4.spines["right"].set_visible(False)
    
    plt.tight_layout()
    fig.savefig(outdir / "lioness_statistics.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: lioness_statistics.png")


def plot_summary_stats(genes: list[str], edges: pd.DataFrame, degrees: list[int], outdir: Path):
    """Create summary statistics figure."""
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.axis("off")
    
    stats = [
        ("Total Genes", f"{len(genes):,}"),
        ("Total Edges", f"{len(edges):,}"),
        ("Mean Degree", f"{np.mean(degrees):.1f}"),
        ("Median Degree", f"{np.median(degrees):.0f}"),
        ("Max Degree", f"{max(degrees)}"),
        ("Min Degree", f"{min(degrees)}"),
        ("Mean |pcorr|", f"{edges['abs_pcorr'].mean():.4f}"),
        ("Median |pcorr|", f"{edges['abs_pcorr'].median():.4f}"),
        ("95th %ile |pcorr|", f"{edges['abs_pcorr'].quantile(0.95):.4f}"),
        ("Positive Edges", f"{(edges['pcorr'] > 0).sum():,} ({100*(edges['pcorr'] > 0).mean():.1f}%)"),
        ("Negative Edges", f"{(edges['pcorr'] < 0).sum():,} ({100*(edges['pcorr'] < 0).mean():.1f}%)"),
    ]
    
    # Table
    table = ax.table(
        cellText=[[s[0], s[1]] for s in stats],
        colLabels=["Metric", "Value"],
        loc="center",
        cellLoc="left",
        colWidths=[0.5, 0.3]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1.2, 1.8)
    
    # Style header
    for i in range(2):
        table[(0, i)].set_facecolor("#3498db")
        table[(0, i)].set_text_props(color="white", fontweight="bold")
    
    ax.set_title("Phase 2 Skeleton Summary Statistics", fontsize=14, fontweight="bold", pad=20)
    
    plt.tight_layout()
    fig.savefig(outdir / "skeleton_summary_stats.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: skeleton_summary_stats.png")


def main():
    ap = argparse.ArgumentParser(description="Generate Phase 2 skeleton diagnostic plots")
    ap.add_argument("--phase2_dir", 
                    default=str(REPO_ROOT / "data/processed/networks/phase2"),
                    help="Phase 2 directory")
    ap.add_argument("--outdir", default="", help="Output directory (default: timestamped)")
    ap.add_argument("--top_n", type=int, default=100, help="Top N hubs for network graph")
    args = ap.parse_args()
    
    phase2 = Path(args.phase2_dir)
    
    if args.outdir:
        outdir = Path(args.outdir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = REPO_ROOT / "plots" / f"{ts}_skeleton_diagnostics"
    outdir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 60)
    print("Phase 2 Skeleton Diagnostics")
    print("=" * 60)
    print(f"Input: {phase2}")
    print(f"Output: {outdir}")
    
    # Load data
    print("\nLoading skeleton data...")
    genes, edges, edge_i, edge_j = load_skeleton(phase2)
    print(f"  Genes: {len(genes)}")
    print(f"  Edges: {len(edges)}")
    
    # Generate plots
    print("\nGenerating plots...")
    
    plot_pcorr_distribution(edges, outdir)
    degrees = plot_degree_distribution(genes, edge_i, edge_j, outdir)
    save_hub_genes(genes, degrees, outdir)
    plot_network_subgraph(genes, edges, degrees, outdir, top_n=args.top_n)
    plot_lioness_statistics(phase2, outdir)
    plot_summary_stats(genes, edges, degrees, outdir)
    
    print(f"\n{'=' * 60}")
    print(f"Done! All outputs in: {outdir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
