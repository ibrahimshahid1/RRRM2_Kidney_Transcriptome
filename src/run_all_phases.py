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
import json
from pathlib import Path
from datetime import datetime

# Repository root
REPO_ROOT = Path(__file__).resolve().parent.parent

# Global run configuration (set by init_run)
RUN_ID = None
RESULTS_DIR = None
NETWORKS_DIR = None


def init_run(run_id: str, max_genes: int, topk: int, num_seeds: int, 
             phases: list, skip_r: bool) -> tuple:
    """Initialize versioned output directories and save run metadata."""
    global RUN_ID, RESULTS_DIR, NETWORKS_DIR
    
    RUN_ID = run_id
    RESULTS_DIR = REPO_ROOT / "data/results" / run_id
    NETWORKS_DIR = REPO_ROOT / "data/processed/networks" / run_id
    
    # Create directories
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    NETWORKS_DIR.mkdir(parents=True, exist_ok=True)
    
    # Save run metadata
    metadata = {
        "run_id": run_id,
        "timestamp": datetime.now().isoformat(),
        "parameters": {
            "max_genes": max_genes,
            "topk": topk,
            "num_seeds": num_seeds,
        },
        "phases": phases,
        "skip_r": skip_r,
    }
    
    # Try to get git commit
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, cwd=REPO_ROOT
        )
        if result.returncode == 0:
            metadata["git_commit"] = result.stdout.strip()
    except:
        pass
    
    with open(RESULTS_DIR / "run_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)
    
    log(f"Run ID: {run_id}")
    log(f"Results: {RESULTS_DIR}")
    log(f"Networks: {NETWORKS_DIR}")
    
    return RESULTS_DIR, NETWORKS_DIR


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
    
    deconv_dir = os.environ.get("RRRM_DECONV_DIR",
                                str(REPO_ROOT / "data/processed/deconvolution/latest"))
    return run_rscript("src/preprocessing/deconvolution.R",
                       args=[f"--outdir={deconv_dir}"], dry_run=dry_run)


def phase_1(dry_run: bool = False, skip_r: bool = False) -> bool:
    """Phase 1: VST + Residualization (R)"""
    log("PHASE 1: VST Normalization + Residualization")
    
    if skip_r:
        log("Skipping R-dependent residualization", "WARN")
        return True
    
    # Pass CLR path from the versioned deconvolution directory
    deconv_dir = os.environ.get("RRRM_DECONV_DIR",
                                str(REPO_ROOT / "data/processed/deconvolution/latest"))
    clr_path = Path(deconv_dir) / "music_segment_direct_proportions_CLR.csv"
    resid_args = [f"--clr={clr_path}"]

    # Check for DCT marker panel to force-include
    dct_panel_path = REPO_ROOT / "data/processed/dct_markers" / "dct_marker_panel.txt"
    if dct_panel_path.exists():
        resid_args.append(f"--force_keep={dct_panel_path}")
        log(f"Passing DCT marker panel to residualization: {dct_panel_path}")
    
    success = run_rscript("src/preprocessing/residualization.R", args=resid_args, dry_run=dry_run)
    if not success:
        return False

    # Export to Python format
    return run_rscript("src/preprocessing/export_phase1.R", dry_run=dry_run)


def phase_1_5(dry_run: bool = False) -> bool:
    """Phase 1.5: Dataset-Derived DCT Marker Discovery"""
    log("PHASE 1.5: DCT Marker Discovery")

    # Use VST (NOT Rtech) to avoid the landmine where Rtech already has CLR regressed out
    vst_path = REPO_ROOT / "data/processed/vst_normalized" / "GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
    deconv_dir = os.environ.get("RRRM_DECONV_DIR",
                                str(REPO_ROOT / "data/processed/deconvolution/latest"))
    clr_path = Path(deconv_dir) / "music_segment_direct_proportions_CLR.csv"
    meta_path = REPO_ROOT / "data/processed/phase1_residuals" / "meta_phase1.tsv.gz"
    outdir = REPO_ROOT / "data/processed/dct_markers"

    if not vst_path.exists():
        log(f"VST file not found: {vst_path}", "WARN")
        return False
    if not clr_path.exists():
        log(f"CLR file not found: {clr_path}", "WARN")
        return False

    return run_python("src.markers.discover_dct", [
        f"--vst={vst_path}",
        f"--meta={meta_path}",
        f"--clr={clr_path}",
        f"--outdir={outdir}",
        "--tech_cols=LibraryBatch,ReadDepth,rRNA",
    ], dry_run=dry_run)


def phase_2(dry_run: bool = False, skip_r: bool = False,
            max_genes: int = 2500, topk: int = 80) -> bool:
    """Phase 2: Network Construction"""
    log("PHASE 2: Network Skeleton + LIONESS + Edge Regression")

    # Get versioned output directory from environment
    networks_dir = os.environ.get("RRRM_NETWORKS_DIR", str(REPO_ROOT / "data/processed/networks/phase2"))

    # Step 2.1-2.3: Build shared sparse skeleton
    # Check if DCT marker panel exists for force-inclusion
    dct_panel_path = REPO_ROOT / "data/processed/dct_markers" / "dct_marker_panel.txt"
    force_args = []
    if dct_panel_path.exists():
        force_args = [f"--force_include={dct_panel_path}"]
        log(f"Including DCT marker panel: {dct_panel_path}")
    else:
        log("No DCT marker panel found; using variance-only gene selection", "WARN")

    if not run_python("src.networks.shared_topology",
                      [f"--max_genes={max_genes}", f"--topk={topk}",
                       f"--outdir={networks_dir}"] + force_args,
                      dry_run=dry_run):
        return False

    # Step A1: LIONESS edge weights
    if not run_python("src.networks.lioness",
                      [f"--phase2_dir={networks_dir}"],
                      dry_run=dry_run):
        return False

    # Step A2-A3: Edge-wise regression (requires R/limma)
    if skip_r:
        log("Skipping R-dependent edge regression", "WARN")
    else:
        if not run_python("src.networks.edge_regression",
                          [f"--phase2_dir={networks_dir}",
                           f"--outdir={networks_dir}/regression"],
                          dry_run=dry_run):
            return False

    return True


def phase_3(dry_run: bool = False, num_seeds: int = 10) -> bool:
    """Phase 3: Embeddings + Procrustes"""
    log("PHASE 3: Node2Vec Embeddings + Procrustes Alignment")

    # Get versioned directories from environment
    networks_dir = os.environ.get("RRRM_NETWORKS_DIR", str(REPO_ROOT / "data/processed/networks/phase2"))
    results_dir = os.environ.get("RRRM_RESULTS_DIR", str(REPO_ROOT / "data/results"))

    # Step 3.2: Node2Vec embeddings
    if not run_python("src.networks.embeddings",
                      [f"--phase2_dir={networks_dir}",
                       f"--reg_dir={networks_dir}/regression",
                       f"--outdir={results_dir}/phase3_embeddings",
                       f"--num_seeds={num_seeds}"],
                      dry_run=dry_run):
        return False

    # Step 3.3: Procrustes alignment
    if not run_python("src.networks.procrustes",
                      [f"--phase2_dir={networks_dir}",
                       f"--emb_dir={results_dir}/phase3_embeddings",
                       f"--outdir={results_dir}/phase3_rewiring",
                       f"--num_seeds={num_seeds}"],
                      dry_run=dry_run):
        return False

    return True


def phase_5(dry_run: bool = False) -> bool:
    """Phase 5: Silent Shifters + Interaction Metrics"""
    log("PHASE 5: Silent Shifters + Interaction Metrics")

    # Get versioned output directory from environment
    results_dir = os.environ.get("RRRM_RESULTS_DIR", str(REPO_ROOT / "data/results"))

    # Interaction metrics
    if not run_python("src.statistics.interaction_metrics",
                      [f"--in_dir={results_dir}/phase3_rewiring",
                       f"--out_dir={results_dir}/phase5_derived"],
                      dry_run=dry_run):
        return False

    # Check for regression results to use as support
    reg_file = Path(results_dir) / "phase6_regression/gene_arm_flight_interaction.tsv"
    reg_args = []
    if reg_file.exists():
        reg_args = [f"--regression_results={reg_file}"]
        log(f"Using regression results for support: {reg_file}")
    else:
        log("Regression results not found; Silent Shifters will lack statistical support", "WARN")

    # Silent shifters for each contrast
    contrasts = [
        "ISS_T_YNG_FLT_minus_GC",
        "ISS_T_OLD_FLT_minus_GC",
        "LAR_YNG_FLT_minus_GC",
        "LAR_OLD_FLT_minus_GC",
    ]

    for contrast in contrasts:
        rewiring_file = Path(results_dir) / f"phase3_rewiring/{contrast}_rewiring_agg.tsv"
        if not rewiring_file.exists():
            log(f"Skipping {contrast} - rewiring file not found", "WARN")
            continue

        if not run_python("src.statistics.silent_shifters",
                          [f"--rewiring={rewiring_file}",
                           f"--outdir={results_dir}/phase5_silent_shifters_strict"] + reg_args,
                          dry_run=dry_run):
            log(f"Warning: silent_shifters failed for {contrast}", "WARN")

    return True


def phase_6(dry_run: bool = False, skip_r: bool = False) -> bool:
    """Phase 6: Uncertainty Estimation + Full Regression"""
    log("PHASE 6: Permutation + Bootstrap + Full Regression")

    # Get versioned output directory from environment
    results_dir = os.environ.get("RRRM_RESULTS_DIR", str(REPO_ROOT / "data/results"))
    networks_dir = os.environ.get("RRRM_NETWORKS_DIR", str(REPO_ROOT / "data/processed/networks/phase2"))

    # Permutation and bootstrap
    if not run_python("src.statistics.permutation_bootstrap",
                      [f"--phase2_dir={networks_dir}",
                       f"--outdir={results_dir}/phase6_uncertainty"],
                      dry_run=dry_run):
        return False

    # Full factorial regression
    if not run_python("src.statistics.full_regression",
                      [f"--phase2_dir={networks_dir}",
                       f"--outdir={results_dir}/phase6_regression"],
                      dry_run=dry_run):
        return False

    # Gene-level DE (R)
    if skip_r:
        log("Skipping R-dependent DE", "WARN")
    else:
        run_rscript("src/statistics/differential_expression.R", dry_run=dry_run)

    return True


def phase_9(dry_run: bool = False) -> bool:
    """Phase 9: Generate Publication-Ready Figures (runs last)"""
    log("PHASE 9: Figure Generation")

    # Get versioned results directory from environment
    results_dir = os.environ.get("RRRM_RESULTS_DIR", str(REPO_ROOT / "data/results"))
    figures_dir = Path(results_dir) / "figures"

    # Generate figures using publication_plots module
    if not run_python("src.visualization.publication_plots",
                      [f"--results_dir={results_dir}",
                       f"--out_dir={figures_dir}"],
                      dry_run=dry_run):
        return False

    return True


def phase_7(dry_run: bool = False) -> bool:
    """Phase 7: Biological Grounding"""
    log("PHASE 7: Biological Grounding + Enrichment")

    # Get versioned output directory from environment
    results_dir = os.environ.get("RRRM_RESULTS_DIR", str(REPO_ROOT / "data/results"))

    contrasts = [
        "ISS_T_YNG_FLT_minus_GC",
        "ISS_T_OLD_FLT_minus_GC",
        "LAR_YNG_FLT_minus_GC",
        "LAR_OLD_FLT_minus_GC",
    ]

    # Ensure ID map exists — use full expressed gene list (not just network HVGs)
    # so that pathway genes outside the 2500-gene skeleton are still mappable.
    map_path = REPO_ROOT / "data/processed/resources" / "id_map.tsv"
    # Rebuild if map is stale (covers fewer genes than what Rtech has)
    rtech_path = REPO_ROOT / "data/processed/phase1_residuals" / "Rtech.tsv.gz"
    needs_rebuild = not map_path.exists()
    if map_path.exists() and rtech_path.exists():
        # Check if existing map covers all expressed genes
        existing_lines = sum(1 for _ in open(map_path)) - 1  # subtract header
        if existing_lines < 5000:  # likely only covers 2500 HVGs
            log(f"ID map has only {existing_lines} genes — rebuilding with full expressed gene set")
            needs_rebuild = True
    if needs_rebuild:
        log("Building Ensembl→Symbol ID map (first run)...")
        networks_dir = os.environ.get("RRRM_NETWORKS_DIR", str(REPO_ROOT / "data/processed/networks/phase2"))
        # Prefer the full Rtech gene list (all expressed genes) over phase2_genes.txt
        rtech_path = REPO_ROOT / "data/processed/phase1_residuals" / "Rtech.tsv.gz"
        gene_list = Path(networks_dir) / "phase2_genes.txt"
        if rtech_path.exists():
            # Extract gene list from Rtech (first column = gene IDs)
            import gzip, csv
            full_gene_path = map_path.parent / "all_expressed_genes.txt"
            full_gene_path.parent.mkdir(parents=True, exist_ok=True)
            with gzip.open(rtech_path, "rt") as f:
                reader = csv.reader(f, delimiter="\t")
                header = next(reader)  # skip header
                gene_ids = [row[0] for row in reader if row]
            with open(full_gene_path, "w") as f:
                f.write("\n".join(gene_ids) + "\n")
            log(f"Extracted {len(gene_ids)} genes from Rtech for ID map")
            run_python("src.data.build_id_map", [
                f"--genes={full_gene_path}",
                f"--outdir={map_path.parent}",
            ], dry_run=dry_run)
        elif gene_list.exists():
            run_python("src.data.build_id_map", [
                f"--genes={gene_list}",
                f"--outdir={map_path.parent}",
            ], dry_run=dry_run)
        else:
            log(f"Gene list not found at {gene_list}, cannot build ID map", "WARN")

    map_args = []
    if map_path.exists():
        map_args = [f"--map={map_path}"]
    else:
        log("ID map not available; gene set enrichment may show 0 overlap", "WARN")

    # Find a gene-level DE file to expand the enrichment universe
    de_dir = REPO_ROOT / "data/processed/gene_level_DE"
    gene_de_args = []
    if de_dir.exists():
        de_files = sorted(de_dir.glob("*_gene_DE.tsv"))
        if de_files:
            # Use the first DE file (all have the same gene universe)
            gene_de_args = [f"--gene_de={de_files[0]}"]
            log(f"Expanding enrichment universe with {de_files[0].name}")

    for contrast in contrasts:
        rewiring_file = Path(results_dir) / f"phase3_rewiring/{contrast}_rewiring_agg.tsv"
        if not rewiring_file.exists():
            log(f"Skipping {contrast} - rewiring file not found", "WARN")
            continue

        outdir = Path(results_dir) / f"phase7_grounding/{contrast}"
        if not run_python("src.enrichment.biological_grounding",
                          [f"--rewiring={rewiring_file}", f"--outdir={outdir}"]
                          + map_args + gene_de_args,
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
    parser.add_argument("--phases", nargs="+", type=float, default=None,
                        help="Specific phases to run (default: all). Supports 1.5 for DCT discovery.")
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
    parser.add_argument("--run-id", type=str, default=None,
                        help="Run identifier for versioned outputs (default: auto-generated)")
    
    args = parser.parse_args()
    
    os.chdir(REPO_ROOT)
    
    # Determine which phases to run first (needed for init_run)
    # Note: Phase 9 (figures) runs last, after all data phases
    phases_available = [0, 1, 1.5, 2, 3, 5, 6, 7, 9]
    if args.phases:
        to_run = sorted(set(args.phases) & set(phases_available))
    else:
        to_run = [p for p in sorted(phases_available) if p >= args.start]
    
    # Generate run ID if not provided
    run_id = args.run_id
    if run_id is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_id = f"run_{timestamp}_{args.max_genes}g"
    
    # Initialize versioned directories (skip for dry-run)
    if not args.dry_run:
        init_run(run_id, args.max_genes, args.topk, args.num_seeds, to_run, args.skip_r)
    else:
        log(f"[DRY-RUN] Would create run: {run_id}")
    
    # Define phase runners - pass versioned paths via environment
    os.environ["RRRM_RUN_ID"] = run_id
    os.environ["RRRM_RESULTS_DIR"] = str(REPO_ROOT / "data/results" / run_id)
    os.environ["RRRM_NETWORKS_DIR"] = str(REPO_ROOT / "data/processed/networks" / run_id)
    os.environ["RRRM_DECONV_DIR"] = str(REPO_ROOT / "data/processed/deconvolution" / run_id)
    
    phases = {
        0: lambda: phase_0(args.dry_run, args.skip_r),
        1: lambda: phase_1(args.dry_run, args.skip_r),
        1.5: lambda: phase_1_5(args.dry_run),
        2: lambda: phase_2(args.dry_run, args.skip_r, args.max_genes, args.topk),
        3: lambda: phase_3(args.dry_run, args.num_seeds),
        5: lambda: phase_5(args.dry_run),
        6: lambda: phase_6(args.dry_run, args.skip_r),
        7: lambda: phase_7(args.dry_run),
        9: lambda: phase_9(args.dry_run),  # Figure generation runs last
    }
    
    log(f"RRRM-2 Pipeline Runner", "INFO")
    log(f"Phases to run: {to_run}")
    log(f"Max genes: {args.max_genes}")
    log(f"Skip R steps: {args.skip_r}")
    log(f"Dry run: {args.dry_run}")
    print()
    
    # Run phases in dependency-aware order (Phase 6 BEFORE Phase 5 for regression support)
    execution_order = [0, 1, 1.5, 2, 3, 6, 5, 7, 9]
    phases_to_execute = [p for p in execution_order if p in to_run]
    
    failed = []
    for phase_num in phases_to_execute:
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
        if RESULTS_DIR:
            log(f"Results in: {RESULTS_DIR}")
        else:
            log(f"Results would be in: data/results/{run_id}/")


if __name__ == "__main__":
    main()

