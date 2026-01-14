# scripts/phase3_procrustes_rewiring.py
"""
Phase 3.3: Procrustes alignment + cosine rewiring (multi-seed)

- Builds anchors from the 4 FLT–GC rewiring tables (lowest median rewiring_abs)
- Aligns embeddings per seed via orthogonal Procrustes
- Computes rewiring = 1 − cosine_similarity
- Reports seed mean/std + rank variance

Usage:
    python scripts/phase3_procrustes_rewiring.py
"""
from __future__ import annotations

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def cosine_distance_rows(A: np.ndarray, B: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """1 - cosine similarity row-wise."""
    An = np.linalg.norm(A, axis=1) + eps
    Bn = np.linalg.norm(B, axis=1) + eps
    sim = (A * B).sum(axis=1) / (An * Bn)
    sim = np.clip(sim, -1.0, 1.0)
    return 1.0 - sim


def pick_phase2_anchors(node_rewiring_dir: Path, genes: list[str], anchor_k: int) -> np.ndarray:
    """
    Anchor recipe:
      - read the 4 FLT_minus_GC rewiring tables
      - compute per gene median rewiring_abs
      - take bottom K genes as anchors (lowest rewiring = most stable)
    """
    needed = [
        "ISS_T_YNG_FLT_minus_GC_node_rewiring.tsv",
        "ISS_T_OLD_FLT_minus_GC_node_rewiring.tsv",
        "LAR_YNG_FLT_minus_GC_node_rewiring.tsv",
        "LAR_OLD_FLT_minus_GC_node_rewiring.tsv",
    ]

    dfs = []
    for fn in needed:
        p = node_rewiring_dir / fn
        if not p.exists():
            raise FileNotFoundError(f"Missing required rewiring table for anchors: {p}")
        df = pd.read_csv(p, sep="\t")
        dfs.append(df[["gene", "rewiring_abs"]].rename(columns={"rewiring_abs": fn}))

    merged = dfs[0]
    for d in dfs[1:]:
        merged = merged.merge(d, on="gene", how="inner")

    cols = [c for c in merged.columns if c != "gene"]
    merged["median_abs"] = merged[cols].median(axis=1)
    merged = merged.sort_values("median_abs", ascending=True, kind="mergesort")

    anchor_genes = merged["gene"].head(anchor_k).tolist()
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    anchor_idx = np.array([gene_to_idx[g] for g in anchor_genes if g in gene_to_idx], dtype=np.int32)

    if anchor_idx.size < min(50, anchor_k):
        raise RuntimeError(f"Too few anchors mapped to genes list: got {anchor_idx.size} / {anchor_k}")

    return anchor_idx


def orthogonal_procrustes_align(B: np.ndarray, A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotate B to best match A using orthogonal Procrustes:
      find R minimizing ||B R - A||_F
      return (B_aligned, R)
    """
    try:
        from scipy.linalg import orthogonal_procrustes
    except Exception as e:
        raise ImportError("scipy is required: pip install scipy") from e

    R, _ = orthogonal_procrustes(B, A)
    return B @ R, R


def load_embedding(emb_dir: Path, condition: str, seed: int) -> np.ndarray:
    p = emb_dir / condition / f"seed_{seed}.npy"
    if not p.exists():
        raise FileNotFoundError(f"Missing embedding: {p}")
    return np.load(p).astype(np.float64)


def main():
    os.chdir(REPO_ROOT)

    ap = argparse.ArgumentParser(description="Phase 3.3: Procrustes alignment + cosine rewiring (multi-seed)")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2")
    ap.add_argument("--node_rewiring_dir", default="data/processed/networks/phase3/node_rewiring")
    ap.add_argument("--emb_dir", default="data/processed/networks/phase3/embeddings")
    ap.add_argument("--outdir", default="data/results/phase3_rewiring")
    ap.add_argument("--anchor_k", type=int, default=200)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--num_seeds", type=int, default=10)
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    node_rewiring_dir = Path(args.node_rewiring_dir)
    emb_dir = Path(args.emb_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 3.3: Procrustes Alignment + Cosine Rewiring")
    print("=" * 60)

    genes = [g for g in (phase2 / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    G = len(genes)
    print(f"\nGenes: {G}")

    # Pick anchors from modeled deltas
    anchors = pick_phase2_anchors(node_rewiring_dir, genes, anchor_k=args.anchor_k)
    (outdir / "anchors.txt").write_text("\n".join([genes[i] for i in anchors.tolist()]) + "\n")
    print(f"Anchors selected: {anchors.size} (saved anchors.txt)")

    # Define the key FLT vs GC comparisons (matches Phase 2 contrasts)
    # Focus on FLT–GC within each Arm/Age
    pairs = [
        ("ISS_T_YNG_FLT_minus_GC", "Pred_YNG_ISS_T_FLT", "Pred_YNG_ISS_T_GC"),
        ("ISS_T_OLD_FLT_minus_GC", "Pred_OLD_ISS_T_FLT", "Pred_OLD_ISS_T_GC"),
        ("LAR_YNG_FLT_minus_GC",   "Pred_YNG_LAR_FLT",   "Pred_YNG_LAR_GC"),
        ("LAR_OLD_FLT_minus_GC",   "Pred_OLD_LAR_FLT",   "Pred_OLD_LAR_GC"),
    ]

    seeds = [args.seed0 + k for k in range(args.num_seeds)]
    print(f"Seeds: {seeds}")

    summary_rows = []

    for contrast_name, cond_FLT, cond_GC in pairs:
        per_seed_rewiring = []

        print(f"\n=== {contrast_name}: {cond_FLT} vs {cond_GC} ===")
        
        for sd in seeds:
            # Load embeddings
            emb_FLT = load_embedding(emb_dir, cond_FLT, sd)
            emb_GC = load_embedding(emb_dir, cond_GC, sd)

            if emb_FLT.shape[0] != G or emb_GC.shape[0] != G:
                raise ValueError(f"Embedding size mismatch: expected {G} rows")

            # Align GC to FLT using anchors
            FLT_anchor = emb_FLT[anchors, :]
            GC_anchor = emb_GC[anchors, :]
            
            # Compute rotation on anchors, apply to full GC
            _, R = orthogonal_procrustes_align(GC_anchor, FLT_anchor)
            GC_aligned = emb_GC @ R

            # Rewiring = 1 - cosine_sim (FLT vs aligned GC)
            rew = cosine_distance_rows(emb_FLT, GC_aligned)
            
            df = pd.DataFrame({"gene": genes, "rewiring": rew.astype(float)})
            df = df.sort_values("rewiring", ascending=False, kind="mergesort").reset_index(drop=True)

            seed_path = outdir / f"{contrast_name}_seed_{sd}.tsv"
            df.to_csv(seed_path, sep="\t", index=False)
            per_seed_rewiring.append(rew.astype(np.float64))

            print(f"  seed {sd}: saved (top gene: {df.iloc[0]['gene']} = {df.iloc[0]['rewiring']:.4f})")

        # Aggregate across seeds
        stack = np.stack(per_seed_rewiring, axis=0)  # (S, G)
        mean_rew = stack.mean(axis=0)
        std_rew = stack.std(axis=0)

        # Rank variance across seeds
        # For each seed, rank genes by rewiring (0=highest), then compute std of ranks
        ranks = np.argsort(np.argsort(-stack, axis=1), axis=1).astype(np.float64)  # 0=highest
        rank_std = ranks.std(axis=0)

        agg = pd.DataFrame(
            {
                "gene": genes,
                "rewiring_mean": mean_rew.astype(float),
                "rewiring_std": std_rew.astype(float),
                "rank_std": rank_std.astype(float),
            }
        ).sort_values("rewiring_mean", ascending=False, kind="mergesort").reset_index(drop=True)

        agg_path = outdir / f"{contrast_name}_rewiring_agg.tsv"
        agg.to_csv(agg_path, sep="\t", index=False)
        print(f"  Aggregate saved: {agg_path}")
        
        # Summary stats for this contrast
        summary_rows.append({
            "contrast": contrast_name,
            "n_genes": G,
            "n_seeds": len(seeds),
            "mean_rewiring": float(mean_rew.mean()),
            "median_rank_std": float(np.median(rank_std)),
            "top_gene": agg.iloc[0]["gene"],
            "top_rewiring": float(agg.iloc[0]["rewiring_mean"]),
        })

    # Save summary
    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "procrustes_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print(f"Done. Summary: {summary_path}")
    print("=" * 60)
    print(f"\nOutputs in: {outdir}")


if __name__ == "__main__":
    main()
