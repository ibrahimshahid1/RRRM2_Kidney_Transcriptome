# scripts/run_phase2_pipeline.py
"""
Phase 2 Unified Pipeline: Cell-Standardized Shared Skeleton Construction

Runs all Phase 2 steps in sequence:
  1. Build skeleton E (cell-standardized partial correlation)
  2. Compute LIONESS Fisher-z weights on E
  3. Edge-wise regression + predicted networks (requires rpy2/limma)

Usage:
    python scripts/run_phase2_pipeline.py
    python scripts/run_phase2_pipeline.py --max_genes 1200 --topk 80 --skip_regression
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run_step(name: str, cmd: list[str], skip: bool = False) -> bool:
    """Run a pipeline step with nice output."""
    print("\n" + "=" * 70)
    print(f"STEP: {name}")
    print("=" * 70)
    
    if skip:
        print("  [SKIPPED]")
        return True
    
    print(f"  Command: {' '.join(cmd)}")
    print("-" * 70)
    
    result = subprocess.run(cmd, cwd=Path(__file__).parent.parent)
    
    if result.returncode != 0:
        print(f"\n  [FAILED] Exit code: {result.returncode}")
        return False
    
    print(f"  [SUCCESS]")
    return True


def main():
    ap = argparse.ArgumentParser(
        description="Run all Phase 2 pipeline steps",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/run_phase2_pipeline.py
  python scripts/run_phase2_pipeline.py --max_genes 1500 --topk 100
  python scripts/run_phase2_pipeline.py --skip_regression  # Skip if no rpy2
"""
    )
    
    # Skeleton params
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--outdir", default="data/processed/networks/phase2")
    ap.add_argument("--max_genes", type=int, default=1200)
    ap.add_argument("--topk", type=int, default=80)
    ap.add_argument("--cell_cols", default="Age,Arm,EnvGroup")
    
    # Regression params
    ap.add_argument("--add_covariates", default="LibraryBatch,SeqInstr,ReadDepth,rRNA")
    
    # Control flags
    ap.add_argument("--skip_skeleton", action="store_true", help="Skip skeleton building")
    ap.add_argument("--skip_lioness", action="store_true", help="Skip LIONESS computation")
    ap.add_argument("--skip_regression", action="store_true", help="Skip edge regression (use if no rpy2)")
    
    args = ap.parse_args()

    print("=" * 70)
    print("PHASE 2 UNIFIED PIPELINE")
    print("Cell-Standardized Shared Skeleton Construction")
    print("=" * 70)
    print(f"\nParameters:")
    print(f"  max_genes    = {args.max_genes}")
    print(f"  topk         = {args.topk}")
    print(f"  cell_cols    = {args.cell_cols}")
    print(f"  outdir       = {args.outdir}")

    scripts_dir = Path(__file__).parent
    success = True

    # -------------------------------------------------------------------------
    # Step 1: Build Skeleton E
    # -------------------------------------------------------------------------
    cmd = [
        sys.executable, str(scripts_dir / "build_phase2_skeleton.py"),
        "--rtech", args.rtech,
        "--meta", args.meta,
        "--outdir", args.outdir,
        "--max_genes", str(args.max_genes),
        "--topk", str(args.topk),
        "--cell_cols", args.cell_cols,
    ]
    if not run_step("Build Skeleton E (Step 2.1-2.3)", cmd, skip=args.skip_skeleton):
        print("\n[PIPELINE FAILED] Skeleton building failed")
        return 1

    # -------------------------------------------------------------------------
    # Step 2: Compute LIONESS on Skeleton
    # -------------------------------------------------------------------------
    cmd = [
        sys.executable, str(scripts_dir / "compute_phase2_lioness_on_skeleton.py"),
        "--rtech", args.rtech,
        "--meta", args.meta,
        "--phase2_dir", args.outdir,
    ]
    if not run_step("LIONESS Fisher-z Weights (Step A1)", cmd, skip=args.skip_lioness):
        print("\n[PIPELINE FAILED] LIONESS computation failed")
        return 1

    # -------------------------------------------------------------------------
    # Step 3: Edge-wise Regression
    # -------------------------------------------------------------------------
    cmd = [
        sys.executable, str(scripts_dir / "phase2_edge_regression.py"),
        "--meta", args.meta,
        "--phase2_dir", args.outdir,
        "--add_covariates", args.add_covariates,
    ]
    if not run_step("Edge-wise Regression (Step A2-A3)", cmd, skip=args.skip_regression):
        if not args.skip_regression:
            print("\n[WARNING] Edge regression failed (likely missing rpy2/limma)")
            print("  Run with --skip_regression or install rpy2")
        success = False

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------
    print("\n" + "=" * 70)
    print("PIPELINE COMPLETE" if success else "PIPELINE COMPLETE (with warnings)")
    print("=" * 70)
    
    outdir = Path(args.outdir)
    print(f"\nOutputs in {outdir}:")
    
    # List outputs
    expected = [
        "phase2_genes.txt",
        "skeleton_edges.tsv",
        "edge_i.npy", "edge_j.npy",
        "lioness_z_edges.npy",
        "lioness_samples.txt",
        "edge_index.tsv",
    ]
    for f in expected:
        p = outdir / f
        status = "✓" if p.exists() else "✗"
        print(f"  {status} {f}")
    
    reg_dir = outdir / "regression"
    if reg_dir.exists():
        print(f"\nRegression outputs in {reg_dir}:")
        for f in sorted(reg_dir.iterdir()):
            print(f"  ✓ {f.name}")
    elif not args.skip_regression:
        print(f"\n  ✗ regression/ (missing - run edge regression separately)")

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
