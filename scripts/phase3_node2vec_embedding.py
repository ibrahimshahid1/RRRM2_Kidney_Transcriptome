# scripts/phase3_node2vec_embedding.py
"""
Phase 3.2: PecanPy node2vec embeddings on Pred_*_z_hat.npy (multi-seed, fixed topology)

Builds biased random walks + trains embeddings using PecanPy for speed.

Graph:
  - fixed topology from Phase 2 skeleton (edge_i, edge_j)
  - condition-specific edge weights: w = abs(tanh(z_hat)) in [0,1)

Targets:
  - by default embed only FLT/GC predicted networks via --patterns
  - multi-seed for robustness; saves embeddings per seed + stability stats

Usage:
  python scripts/phase3_node2vec_embedding.py
  python scripts/phase3_node2vec_embedding.py --num_seeds 2 --num_walks 20 --walk_length 40  # quick test
  python scripts/phase3_node2vec_embedding.py --patterns "Pred_*_FLT_z_hat.npy,Pred_*_GC_z_hat.npy"
"""
from __future__ import annotations

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import random
import inspect

REPO_ROOT = Path(__file__).resolve().parents[1]


def write_edgelist_weighted(
    path: Path,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    w_edge: np.ndarray,
    node_names: list[str],
    w_min: float = 1e-8,
) -> None:
    """
    Write weighted edgelist expected by PecanPy:
      node_u <tab> node_v <tab> weight
    We use integer node IDs as strings to keep things compact.
    """
    # PecanPy accepts node ids as strings; we'll use "0..N-1"
    # Filter extremely tiny weights to reduce I/O and walk noise.
    m = w_edge > w_min
    ei = edge_i[m]
    ej = edge_j[m]
    ww = w_edge[m]

    with path.open("w", encoding="utf-8") as f:
        for a, b, w in zip(ei.tolist(), ej.tolist(), ww.tolist()):
            f.write(f"{a}\t{b}\t{w:.8g}\n")


def call_embed_version_safe(g, **kwargs):
    """
    Call PecanPy g.embed(...) but only with args supported by the installed version.
    This prevents crashes like: unexpected keyword argument 'workers' or 'seed'.
    """
    sig = inspect.signature(g.embed)
    allowed = set(sig.parameters.keys())
    safe = {k: v for k, v in kwargs.items() if k in allowed}
    return g.embed(**safe)


def main():
    os.chdir(REPO_ROOT)

    ap = argparse.ArgumentParser(description="Phase 3.2: PecanPy node2vec embeddings (multi-seed, fixed topology)")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2")
    ap.add_argument("--reg_dir", default="data/processed/networks/phase2/regression")
    ap.add_argument("--outdir", default="data/processed/networks/phase3/embeddings")

    # FLT/GC-only focus
    ap.add_argument(
        "--patterns",
        default="Pred_*_FLT_z_hat.npy,Pred_*_GC_z_hat.npy",
        help="Comma-separated glob patterns inside reg_dir (e.g., Pred_*_FLT_z_hat.npy,Pred_*_GC_z_hat.npy)",
    )

    # node2vec hyperparams (your spec)
    ap.add_argument("--dimensions", type=int, default=128)
    ap.add_argument("--walk_length", type=int, default=80)
    ap.add_argument("--num_walks", type=int, default=200)
    ap.add_argument("--window", type=int, default=10)  # skip-gram window
    ap.add_argument("--epochs", type=int, default=1)   # keep 1; increase if needed
    ap.add_argument("--p", type=float, default=0.25)
    ap.add_argument("--q", type=float, default=4.0)

    ap.add_argument("--num_seeds", type=int, default=10)
    ap.add_argument("--seed0", type=int, default=0)

    # Practical speed knobs
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--weight_min", type=float, default=1e-8, help="Drop edges with weight <= this before writing edgelist")
    ap.add_argument("--keep_edgelists", action="store_true", help="Keep generated .edgelist files (debug/audit)")
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    reg = Path(args.reg_dir)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 3.2: PecanPy Node2Vec Embeddings (Multi-seed)")
    print("=" * 60)

    # Load skeleton
    genes = [g for g in (phase2 / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    edge_i = np.load(phase2 / "edge_i.npy")
    edge_j = np.load(phase2 / "edge_j.npy")
    num_nodes = len(genes)
    E = edge_i.shape[0]

    print(f"\nNodes: {num_nodes}, Edges: {E}")
    print(f"Hyperparams: dim={args.dimensions} walks={args.num_walks} len={args.walk_length} p={args.p} q={args.q}")
    print(f"Seeds: {args.num_seeds} (starting from {args.seed0})")
    print(f"Patterns: {args.patterns}")

    # Find predicted networks (FLT/GC only by default)
    patterns = [p.strip() for p in args.patterns.split(",") if p.strip()]
    pred_files: list[Path] = []
    for pat in patterns:
        pred_files.extend(reg.glob(pat))
    pred_files = sorted(set(pred_files))

    if not pred_files:
        raise FileNotFoundError(f"No predicted networks matched patterns={patterns} in {reg}")

    print(f"\nFound {len(pred_files)} predicted networks to embed")

    # Import PecanPy
    try:
        from pecanpy import pecanpy as pecanpy_mod
    except Exception as e:
        raise ImportError("pecanpy is required. Install with: pip install pecanpy") from e

    all_stats_rows = []

    for f in pred_files:
        cond = f.stem.replace("_z_hat", "")
        cond_dir = out / cond
        cond_dir.mkdir(parents=True, exist_ok=True)

        z = np.load(f).astype(np.float64)
        if z.shape[0] != E:
            raise ValueError(f"{f.name}: expected {E} edges but got {z.shape[0]}")

        # Convert Fisher-z to correlation-ish magnitude weights
        w_edge = np.abs(np.tanh(z)).astype(np.float32)

        # Diagnostics
        w_mean = float(w_edge.mean())
        w_99 = float(np.percentile(w_edge, 99))
        n_tiny = int((w_edge <= args.weight_min).sum())

        print(f"\n=== {cond} ===")
        print(f"  Edge weight stats: mean={w_mean:.4f}, 99th={w_99:.4f}, #dropped(<= {args.weight_min:g})={n_tiny}")

        # Build a weighted edgelist file for PecanPy (fast + simple)
        edgelist_path = cond_dir / f"{cond}.edgelist"
        write_edgelist_weighted(
            edgelist_path,
            edge_i=edge_i,
            edge_j=edge_j,
            w_edge=w_edge,
            node_names=genes,
            w_min=float(args.weight_min),
        )

        # Multi-seed embeddings
        embs = []
        seeds = [args.seed0 + k for k in range(args.num_seeds)]

        for sd in seeds:
            print(f"  Embedding seed {sd}...", end=" ", flush=True)

            # Control RNG outside pecanpy (embed() may not accept seed)
            random.seed(sd)
            np.random.seed(sd)

            # PecanPy graph + embedding
            # Use SparseOTF for speed/memory on large graphs; it builds transition probs on-the-fly.
            g = pecanpy_mod.SparseOTF(
                p=float(args.p),
                q=float(args.q),
                workers=int(args.workers),
                verbose=False,
            )
            g.read_edg(str(edgelist_path), weighted=True, directed=False)

            # node2vec embedding (PecanPy includes training)
            # Use a version-safe helper to handle API variations (seed, workers, names)
            emb = call_embed_version_safe(
                g,
                dim=int(args.dimensions),
                dimensions=int(args.dimensions),   # some versions use "dimensions"
                num_walks=int(args.num_walks),
                walk_length=int(args.walk_length),
                window=int(args.window),
                window_size=int(args.window),      # some versions use "window_size"
                epochs=int(args.epochs),
                epoch=int(args.epochs),            # some versions use "epoch"
                verbose=False,
                seed=int(sd),                      # ignored if unsupported
                workers=1,                         # ignored if unsupported
            )
            

            # PecanPy returns dict-like mapping or ndarray depending on version.
            # Normalize to (N, D) array ordered by node id 0..N-1.
            if isinstance(emb, dict):
                mat = np.zeros((num_nodes, int(args.dimensions)), dtype=np.float32)
                for k, v in emb.items():
                    mat[int(k), :] = np.asarray(v, dtype=np.float32)
            else:
                mat = np.asarray(emb, dtype=np.float32)
                if mat.shape[0] != num_nodes:
                    raise RuntimeError(f"PecanPy returned {mat.shape[0]} nodes, expected {num_nodes}")

            np.save(cond_dir / f"seed_{sd}.npy", mat)
            embs.append(mat)
            print(f"saved {mat.shape}")

        # Seed stability
        stack = np.stack(embs, axis=0)  # (S, N, D)
        mean = stack.mean(axis=0)
        std = stack.std(axis=0)
        mean_norm = np.linalg.norm(mean, axis=1)
        std_norm = np.linalg.norm(std, axis=1)

        stats_df = pd.DataFrame(
            {"gene": genes, "mean_norm": mean_norm.astype(float), "std_norm": std_norm.astype(float)}
        )
        stats_path = cond_dir / "seed_stability.tsv"
        stats_df.to_csv(stats_path, sep="\t", index=False)
        print(f"  Saved seed stability: {stats_path}")

        all_stats_rows.append(
            {
                "condition": cond,
                "n_genes": num_nodes,
                "dim": int(args.dimensions),
                "num_seeds": int(args.num_seeds),
                "mean_std_norm": float(stats_df["std_norm"].mean()),
                "median_std_norm": float(stats_df["std_norm"].median()),
                "edge_weight_mean": w_mean,
                "edge_weight_99pct": w_99,
                "edges_dropped": int(n_tiny),
            }
        )

        if not args.keep_edgelists:
            try:
                edgelist_path.unlink(missing_ok=True)
            except Exception:
                pass

    summary = pd.DataFrame(all_stats_rows)
    summary_path = out / "embedding_seed_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print(f"Done. Summary: {summary_path}")
    print("=" * 60)


if __name__ == "__main__":
    main()
