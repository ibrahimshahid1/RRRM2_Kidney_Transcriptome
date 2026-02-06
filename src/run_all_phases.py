#!/usr/bin/env python3
"""
RRRM-2 Full Pipeline Runner
===========================
Runs all phases of the RRRM-2 network rewiring analysis pipeline.

Usage:
    python scripts/run_all_phases.py                    # Run all phases
    python scripts/run_all_phases.py --phases 2 3       # Run specific phases
    python scripts/run_all_phases.py --skip-r           # Skip R-dependent steps
    python scripts/run_all_phases.py --dry-run          # Show commands without running

Requires:
    - Python environment with dependencies from requirements.txt
    - R with DESeq2, limma, sva, MuSiC (for Phase 0-1 and edge regression)
"""

from __future__ import annotations
import argparse
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

# Repository root
REPO_ROOT = Path(__file__).resolve().parent.parent


def log(msg: str, level: str = "INFO"):
    """Print timestamped log message."""
    ts = datetime.now().strftime("%H:%M:%S")
    prefix = {"INFO": "[INFO]", "OK": "[OK]", "WARN": "[WARN]", "ERROR": "[ERROR]", "RUN": "[RUN]"}
    print(f"{ts} {prefix.get(level, '*')} {msg}")


def run_python(module: str, args: list = None, dry_run: bool = False) -> bool:
    """Run a Python module."""
    cmd = [sys.executable, "-m", module] + (args or [])
    log(f"python -m {module} {' '.join(args or [])}", "RUN")
    
    if dry_run:
        return True
    
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        log(f"Failed: {module}", "ERROR")
        return False
    return True


def run_rscript(script: str, args: list = None, dry_run: bool = False) -> bool:
    """Run an R script."""
    script_path = REPO_ROOT / script
    cmd = ["Rscript", str(script_path)] + (args or [])
    log(f"Rscript {script}", "RUN")
    
    if dry_run:
        return True
    
    if not script_path.exists():
        log(f"Script not found: {script}", "WARN")
        return False
    
    result = subprocess.run(cmd, cwd=REPO_ROOT)
    if result.returncode != 0:
        log(f"Failed: {script}", "ERROR")
        return False
    return True


def phase_0(dry_run: bool = False, skip_r: bool = False) -> bool:
    """Phase 0: Deconvolution (R)"""
    log("PHASE 0: Cell-type Deconvolution")
    
    if skip_r:
        log("Skipping R-dependent deconvolution", "WARN")
        return True
    
    return run_rscript("src/preprocessing/deconvolution.R", dry_run=dry_run)


def phase_1(dry_run: bool = False, skip_r: bool = False) -> bool:
    """Phase 1: VST + Residualization (R)"""
    log("PHASE 1: VST Normalization + Residualization")
    
    if skip_r:
        log("Skipping R-dependent residualization", "WARN")
        return True
    
    success = run_rscript("src/preprocessing/residualization.R", dry_run=dry_run)
    if not success:
        return False
    
    # Export to Python format
    return run_rscript("src/preprocessing/export_phase1.R", dry_run=dry_run)


def phase_2(dry_run: bool = False, skip_r: bool = False, 
            max_genes: int = 2500, topk: int = 80) -> bool:
    """Phase 2: Network Construction"""
    log("PHASE 2: Network Skeleton + LIONESS + Edge Regression")
    
    # Step 2.1-2.3: Build shared sparse skeleton
    if not run_python("src.networks.shared_topology", 
                      [f"--max_genes={max_genes}", f"--topk={topk}"], 
                      dry_run=dry_run):
        return False
    
    # Step A1: LIONESS edge weights
    if not run_python("src.networks.lioness", dry_run=dry_run):
        return False
    
    # Step A2-A3: Edge-wise regression (requires R/limma)
    if skip_r:
        log("Skipping R-dependent edge regression", "WARN")
    else:
        if not run_python("src.networks.edge_regression", dry_run=dry_run):
            return False
    
    return True


def phase_3(dry_run: bool = False, num_seeds: int = 10) -> bool:
    """Phase 3: Embeddings + Procrustes"""
    log("PHASE 3: Node2Vec Embeddings + Procrustes Alignment")
    
    # Step 3.2: Node2Vec embeddings
    if not run_python("src.networks.embeddings", 
                      [f"--num_seeds={num_seeds}"], 
                      dry_run=dry_run):
        return False
    
    # Step 3.3: Procrustes alignment
    if not run_python("src.networks.procrustes", dry_run=dry_run):
        return False
    
    return True


def phase_5(dry_run: bool = False) -> bool:
    """Phase 5: Silent Shifters + Interaction Metrics"""
    log("PHASE 5: Silent Shifters + Interaction Metrics")
    
    # Interaction metrics
    if not run_python("src.statistics.interaction_metrics", dry_run=dry_run):
        return False
    
    # Silent shifters for each contrast
    contrasts = [
        "ISS_T_YNG_FLT_minus_GC",
        "ISS_T_OLD_FLT_minus_GC", 
        "LAR_YNG_FLT_minus_GC",
        "LAR_OLD_FLT_minus_GC",
    ]
    
    for contrast in contrasts:
        rewiring_file = REPO_ROOT / f"data/results/phase3_rewiring/{contrast}_rewiring_agg.tsv"
        if not rewiring_file.exists():
            log(f"Skipping {contrast} - rewiring file not found", "WARN")
            continue
        
        if not run_python("src.statistics.silent_shifters",
                          [f"--rewiring={rewiring_file}"],
                          dry_run=dry_run):
            log(f"Warning: silent_shifters failed for {contrast}", "WARN")
    
    return True


def phase_6(dry_run: bool = False, skip_r: bool = False) -> bool:
    """Phase 6: Uncertainty Estimation + Full Regression"""
    log("PHASE 6: Permutation + Bootstrap + Full Regression")
    
    # Permutation and bootstrap
    if not run_python("src.statistics.permutation_bootstrap", dry_run=dry_run):
        return False
    
    # Full factorial regression
    if not run_python("src.statistics.full_regression", dry_run=dry_run):
        return False
    
    # Gene-level DE (R)
    if skip_r:
        log("Skipping R-dependent DE", "WARN")
    else:
        run_rscript("src/statistics/differential_expression.R", dry_run=dry_run)
    
    return True


def phase_7(dry_run: bool = False) -> bool:
    """Phase 7: Biological Grounding"""
    log("PHASE 7: Biological Grounding + Enrichment")
    
    contrasts = [
        "ISS_T_YNG_FLT_minus_GC",
        "ISS_T_OLD_FLT_minus_GC",
        "LAR_YNG_FLT_minus_GC", 
        "LAR_OLD_FLT_minus_GC",
    ]
    
    for contrast in contrasts:
        rewiring_file = REPO_ROOT / f"data/results/phase3_rewiring/{contrast}_rewiring_agg.tsv"
        if not rewiring_file.exists():
            log(f"Skipping {contrast} - rewiring file not found", "WARN")
            continue
        
        outdir = REPO_ROOT / f"data/results/phase7_grounding/{contrast}"
        if not run_python("src.enrichment.biological_grounding",
                          [f"--rewiring={rewiring_file}", f"--outdir={outdir}"],
                          dry_run=dry_run):
            log(f"Warning: enrichment failed for {contrast}", "WARN")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Run RRRM-2 pipeline (all phases or selected)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/run_all_phases.py                  # Run all phases
  python scripts/run_all_phases.py --phases 2 3 5   # Run phases 2, 3, 5
  python scripts/run_all_phases.py --skip-r         # Skip R-dependent steps
  python scripts/run_all_phases.py --start 3        # Start from phase 3
        """
    )
    parser.add_argument("--phases", nargs="+", type=int, default=None,
                        help="Specific phases to run (default: all)")
    parser.add_argument("--start", type=int, default=0,
                        help="Start from this phase (default: 0)")
    parser.add_argument("--skip-r", action="store_true",
                        help="Skip R-dependent steps")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--max-genes", type=int, default=2500,
                        help="Max genes for skeleton (default: 2500)")
    parser.add_argument("--topk", type=int, default=80,
                        help="Top-k neighbors per gene (default: 80)")
    parser.add_argument("--num-seeds", type=int, default=10,
                        help="Number of embedding seeds (default: 10)")
    
    args = parser.parse_args()
    
    os.chdir(REPO_ROOT)
    
    # Define phase runners
    phases = {
        0: lambda: phase_0(args.dry_run, args.skip_r),
        1: lambda: phase_1(args.dry_run, args.skip_r),
        2: lambda: phase_2(args.dry_run, args.skip_r, args.max_genes, args.topk),
        3: lambda: phase_3(args.dry_run, args.num_seeds),
        5: lambda: phase_5(args.dry_run),
        6: lambda: phase_6(args.dry_run, args.skip_r),
        7: lambda: phase_7(args.dry_run),
    }
    
    # Determine which phases to run
    if args.phases:
        to_run = sorted(set(args.phases) & set(phases.keys()))
    else:
        to_run = [p for p in sorted(phases.keys()) if p >= args.start]
    
    log(f"RRRM-2 Pipeline Runner", "INFO")
    log(f"Phases to run: {to_run}")
    log(f"Skip R steps: {args.skip_r}")
    log(f"Dry run: {args.dry_run}")
    print()
    
    # Run phases
    failed = []
    for phase_num in to_run:
        if not phases[phase_num]():
            failed.append(phase_num)
            log(f"Phase {phase_num} failed", "ERROR")
    
    # Summary
    print()
    if failed:
        log(f"Pipeline completed with errors in phases: {failed}", "WARN")
        sys.exit(1)
    else:
        log("Pipeline completed successfully", "OK")
        log(f"Results in: {REPO_ROOT / 'data/results/'}")


if __name__ == "__main__":
    main()
