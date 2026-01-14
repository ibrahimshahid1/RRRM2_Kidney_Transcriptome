# scripts/run_phase3_pipeline.py
"""
Phase 3 Unified Pipeline: Node2vec Embedding

Runs all Phase 3 steps in sequence:
  1. Node rewiring from delta-z (fast)
  2. Node2vec embeddings (multi-seed, topology fixed)
  3. Procrustes alignment + cosine rewiring

Usage:
    python scripts/run_phase3_pipeline.py
    python scripts/run_phase3_pipeline.py --num_seeds 2 --num_walks 20  # quick test
"""
from __future__ import annotations

import os
import sys
import argparse
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def run(cmd: list[str], description: str = "") -> bool:
    """Run a command with nice output."""
    print("\n" + "=" * 70)
    if description:
        print(f"STEP: {description}")
    print(" ".join(cmd))
    print("=" * 70)
    
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    
    if result.returncode != 0:
        print(f"\n[FAILED] Exit code: {result.returncode}")
        return False
    
    print("[SUCCESS]")
    return True


def main():
    os.chdir(REPO_ROOT)

    ap = argparse.ArgumentParser(
        description="Phase 3 unified pipeline runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/run_phase3_pipeline.py
  python scripts/run_phase3_pipeline.py --num_seeds 2 --num_walks 20  # quick test
  python scripts/run_phase3_pipeline.py --skip_node2vec --skip_procrustes  # rewiring only
"""
    )
    
    # Control flags
    ap.add_argument("--skip_rewiring", action="store_true", help="Skip node rewiring from deltas")
    ap.add_argument("--skip_node2vec", action="store_true", help="Skip node2vec embedding")
    ap.add_argument("--skip_procrustes", action="store_true", help="Skip Procrustes alignment")

    # Node2vec tuning
    ap.add_argument("--num_seeds", type=int, default=10, help="Number of random seeds")
    ap.add_argument("--num_walks", type=int, default=200, help="Walks per node")
    ap.add_argument("--walk_length", type=int, default=80, help="Walk length")
    ap.add_argument("--dimensions", type=int, default=128, help="Embedding dimensions")
    ap.add_argument("--p", type=float, default=0.25, help="Return parameter")
    ap.add_argument("--q", type=float, default=4.0, help="In-out parameter")

    # Procrustes anchors
    ap.add_argument("--anchor_k", type=int, default=200, help="Number of anchor genes")

    args = ap.parse_args()

    py = sys.executable

    print("=" * 70)
    print("PHASE 3 UNIFIED PIPELINE")
    print("Node2vec Embedding (Multi-seed, Topology Fixed)")
    print("=" * 70)
    print(f"\nParameters:")
    print(f"  num_seeds    = {args.num_seeds}")
    print(f"  num_walks    = {args.num_walks}")
    print(f"  walk_length  = {args.walk_length}")
    print(f"  dimensions   = {args.dimensions}")
    print(f"  p, q         = {args.p}, {args.q}")
    print(f"  anchor_k     = {args.anchor_k}")

    success = True

    # -------------------------------------------------------------------------
    # Step 1: Node rewiring from deltas
    # -------------------------------------------------------------------------
    if not args.skip_rewiring:
        if not run(
            [py, "scripts/phase3_node_rewiring_from_deltas.py"],
            "Phase 3.1: Node Rewiring from Deltas"
        ):
            print("\n[PIPELINE FAILED] Node rewiring failed")
            return 1
    else:
        print("\n[SKIPPED] Phase 3.1: Node Rewiring from Deltas")

    # -------------------------------------------------------------------------
    # Step 2: Node2vec embeddings
    # -------------------------------------------------------------------------
    if not args.skip_node2vec:
        if not run(
            [
                py, "scripts/phase3_node2vec_embedding.py",
                "--num_seeds", str(args.num_seeds),
                "--num_walks", str(args.num_walks),
                "--walk_length", str(args.walk_length),
                "--dimensions", str(args.dimensions),
                "--p", str(args.p),
                "--q", str(args.q),
            ],
            "Phase 3.2: Node2Vec Embeddings"
        ):
            print("\n[PIPELINE FAILED] Node2vec embedding failed")
            return 1
    else:
        print("\n[SKIPPED] Phase 3.2: Node2Vec Embeddings")

    # -------------------------------------------------------------------------
    # Step 3: Procrustes alignment + rewiring
    # -------------------------------------------------------------------------
    if not args.skip_procrustes:
        if not run(
            [
                py, "scripts/phase3_procrustes_rewiring.py",
                "--num_seeds", str(args.num_seeds),
                "--anchor_k", str(args.anchor_k),
            ],
            "Phase 3.3: Procrustes Alignment + Rewiring"
        ):
            print("\n[PIPELINE FAILED] Procrustes alignment failed")
            return 1
    else:
        print("\n[SKIPPED] Phase 3.3: Procrustes Alignment + Rewiring")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE" if success else "PIPELINE COMPLETE (with warnings)")
    print("=" * 70)
    
    # List outputs
    outputs = [
        ("Node rewiring", "data/processed/networks/phase3/node_rewiring"),
        ("Embeddings", "data/processed/networks/phase3/embeddings"),
        ("Procrustes rewiring", "data/results/phase3_rewiring"),
    ]
    
    print("\nOutputs:")
    for name, path in outputs:
        p = REPO_ROOT / path
        status = "✓" if p.exists() else "✗"
        n_files = len(list(p.glob("*"))) if p.exists() else 0
        print(f"  {status} {name}: {path} ({n_files} items)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
