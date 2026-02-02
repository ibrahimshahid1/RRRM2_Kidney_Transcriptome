# scripts/run_phase4_7_pipeline.py
"""
Master runner for Phases 4-7 of the kidney transcriptome network analysis.

Executes in order:
  1. Phase 5: Derived metrics (interaction + persistence)
  2. Phase 4: Anchor QC report
  3. Phase 6: Permutation + bootstrap uncertainty
  4. Phase 7: Biological grounding
  5. Phase 5: Silent shifters (with Phase 6 support)

Usage:
  python scripts/run_phase4_7_pipeline.py [--quick]

Options:
  --quick   Run Phase 6 with minimal iterations for testing (10 perm, 10 boot)
"""

from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = REPO_ROOT / "scripts"


def run_script(name: str, args: list[str] = None):
    """Run a Python script and check for errors."""
    script_path = SCRIPTS / name
    cmd = [sys.executable, str(script_path)] + (args or [])
    print(f"\n{'='*60}")
    print(f"Running: {name}")
    print(f"{'='*60}")
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    if result.returncode != 0:
        print(f"[ERROR] {name} failed with exit code {result.returncode}")
        sys.exit(result.returncode)
    print(f"[OK] {name} completed successfully")


def main():
    ap = argparse.ArgumentParser(description="Run Phases 4-7 pipeline")
    ap.add_argument("--quick", action="store_true",
                    help="Quick mode: fewer permutations/bootstraps for testing")
    ap.add_argument("--K_perm", type=int, default=200,
                    help="Number of permutations for Phase 6 (default: 200)")
    ap.add_argument("--B_boot", type=int, default=200,
                    help="Number of bootstraps for Phase 6 (default: 200)")
    ap.add_argument("--skip_phase6", action="store_true",
                    help="Skip Phase 6 (useful if already run)")
    args = ap.parse_args()

    # Override K/B for quick mode
    if args.quick:
        args.K_perm = 10
        args.B_boot = 10
        print("[QUICK MODE] Using 10 permutations and 10 bootstraps")

    print(f"\n{'#'*60}")
    print("# PHASES 4-7 PIPELINE")
    print(f"{'#'*60}")

    # Phase 5: Derived metrics
    run_script("phase5_derive_interaction_persistence.py")

    # Phase 4: Anchor QC
    run_script("phase4_anchor_qc_report.py")

    # Phase 6: Uncertainty (permutation + bootstrap)
    if not args.skip_phase6:
        run_script("phase6_perm_bootstrap_node_rewiring.py", [
            f"--K_perm={args.K_perm}",
            f"--B_boot={args.B_boot}",
        ])
    else:
        print("\n[SKIP] Phase 6 (--skip_phase6 flag)")

    # Phase 7: Biological grounding
    # Run on all 4 contrasts
    rewiring_dir = REPO_ROOT / "data/results/phase3_rewiring"
    contrasts = [
        "ISS_T_YNG_FLT_minus_GC",
        "ISS_T_OLD_FLT_minus_GC",
        "LAR_YNG_FLT_minus_GC",
        "LAR_OLD_FLT_minus_GC",
    ]
    
    for contrast in contrasts:
        rewiring_file = rewiring_dir / f"{contrast}_rewiring_agg.tsv"
        if rewiring_file.exists():
            run_script("phase7_grounding_fast.py", [
                f"--rewiring={rewiring_file}",
                f"--outdir={REPO_ROOT / 'data/results/phase7_grounding' / contrast}",
            ])
        else:
            print(f"[SKIP] {contrast} - rewiring file not found")

    # Phase 5: Silent shifters (with Phase 6 support)
    perm_dir = REPO_ROOT / "data/results/phase6_uncertainty"
    silent_outdir = REPO_ROOT / "data/results/phase5_silent_shifters_strict"
    
    for contrast in contrasts:
        rewiring_file = rewiring_dir / f"{contrast}_rewiring_agg.tsv"
        perm_file = perm_dir / f"{contrast}_perm_pvals.tsv"
        
        if not rewiring_file.exists():
            print(f"[SKIP] Silent shifters for {contrast} - rewiring file not found")
            continue
            
        script_args = [f"--rewiring={rewiring_file}", f"--outdir={silent_outdir}"]
        if perm_file.exists():
            script_args.append(f"--perm={perm_file}")
        
        run_script("phase5_build_silent_shifters_strict.py", script_args)

    # Summary
    print(f"\n{'#'*60}")
    print("# PIPELINE COMPLETE")
    print(f"{'#'*60}")
    print("\nOutputs:")
    print(f"  Phase 4: {REPO_ROOT / 'data/results/phase4_anchor_qc'}")
    print(f"  Phase 5 derived: {REPO_ROOT / 'data/results/phase5_derived'}")
    print(f"  Phase 5 silent: {silent_outdir}")
    print(f"  Phase 6: {perm_dir}")
    print(f"  Phase 7: {REPO_ROOT / 'data/results/phase7_grounding'}")


if __name__ == "__main__":
    main()
