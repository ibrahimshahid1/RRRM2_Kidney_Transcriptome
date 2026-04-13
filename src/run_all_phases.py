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


# ── Artifact resolution ──────────────────────────────────────────────────────

def _find_latest_run_dir() -> Path | None:
    """Return the most recent run_* directory under data/results/ that exists."""
    results_root = REPO_ROOT / "data" / "results"
    if not results_root.exists():
        return None
    run_dirs = sorted(results_root.glob("run_*"), reverse=True)
    for rd in run_dirs:
        if rd.is_dir():
            return rd
    return None


def find_artifact(relpath: str, phase_subdir: str | None = None) -> Path | None:
    """Resolve an artifact from the current run or the most recent prior run.

    Search order:
      1. Current run's results dir  (data/results/<current_run>/<relpath>)
      2. Most recent prior run       (data/results/<latest_run>/<relpath>)
      3. Legacy global location      (data/processed/<relpath>)

    *phase_subdir* is unused but kept for future sub-folder logic.

    Returns the first existing Path, or None.
    """
    candidates: list[Path] = []

    # 1. Current run
    cur = os.environ.get("RRRM_RESULTS_DIR")
    if cur:
        candidates.append(Path(cur) / relpath)

    # 2. Most recent prior run
    latest = _find_latest_run_dir()
    if latest and (not cur or str(latest) != cur):
        candidates.append(latest / relpath)

    # 3. Legacy global locations (for backward compat)
    candidates.append(REPO_ROOT / "data" / "processed" / relpath)

    for c in candidates:
        if c.exists():
            return c
    return None


def init_run(run_id: str, max_genes: int, topk: int, num_seeds: int, 
             phases: list, skip_r: bool) -> tuple:
    """Initialize versioned output directories and save run metadata.

    All pipeline outputs live under a single run directory:
        data/results/<run_id>/
            deconvolution/     ← Phase 0
            phase1_residuals/  ← Phase 1
            dct_markers/       ← Phase 1.5
            networks/          ← Phase 2
            phase3_embeddings/ ← Phase 3
            phase3_rewiring/   ← Phase 3
            ...                ← Phase 5-9
            run_metadata.json
    """
    global RUN_ID, RESULTS_DIR, NETWORKS_DIR
    
    RUN_ID = run_id
    RESULTS_DIR = REPO_ROOT / "data/results" / run_id
    NETWORKS_DIR = RESULTS_DIR / "networks"
    
    # Create all sub-directories up front
    for subdir in ["deconvolution", "phase1_residuals", "dct_markers", "networks"]:
        (RESULTS_DIR / subdir).mkdir(parents=True, exist_ok=True)
    
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


def build_id_map_if_needed(dry_run: bool = False) -> None:
    """Build/rebuild the Ensembl→Symbol ID map if needed.

    Uses the full Rtech gene list (all expressed genes) to ensure pathway
    genes outside the network skeleton are still mappable.  Also auto-resolves
    curated gene symbols from config/gene_sets.yaml.
    """
    import gzip, csv

    map_path = REPO_ROOT / "data/processed/resources" / "id_map.tsv"
    rtech_resolved = find_artifact("phase1_residuals/Rtech.tsv.gz")

    needs_rebuild = not map_path.exists()
    if map_path.exists() and rtech_resolved and rtech_resolved.exists():
        existing_lines = sum(1 for _ in open(map_path)) - 1
        if existing_lines < 5000:
            log(f"ID map has only {existing_lines} genes — rebuilding with full gene set")
            needs_rebuild = True

    if not needs_rebuild:
        return

    log("Building Ensembl→Symbol ID map...")

    networks_dir = os.environ.get("RRRM_NETWORKS_DIR", "")
    gene_list = Path(networks_dir) / "phase2_genes.txt" if networks_dir else None

    if rtech_resolved and rtech_resolved.exists():
        full_gene_path = map_path.parent / "all_expressed_genes.txt"
        full_gene_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(rtech_resolved, "rt") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # skip header
            gene_ids = [row[0] for row in reader if row]
        with open(full_gene_path, "w") as f:
            f.write("\n".join(gene_ids) + "\n")
        log(f"Extracted {len(gene_ids)} genes from Rtech for ID map")
        run_python("src.data.build_id_map", [
            f"--genes={full_gene_path}",
            f"--outdir={map_path.parent}",
        ], dry_run=dry_run)
    elif gene_list and gene_list.exists():
        run_python("src.data.build_id_map", [
            f"--genes={gene_list}",
            f"--outdir={map_path.parent}",
        ], dry_run=dry_run)
    else:
        log("Gene list not found, cannot build ID map", "WARN")


def phase_0(dry_run: bool = False, skip_r: bool = False) -> bool:
    """Phase 0: Deconvolution (R)"""
    log("PHASE 0: Cell-type Deconvolution")
    
    if skip_r:
        log("Skipping R-dependent deconvolution", "WARN")
        return True
    
    deconv_dir = os.environ.get("RRRM_DECONV_DIR")
    return run_rscript("src/preprocessing/deconvolution.R",
                       args=[f"--outdir={deconv_dir}"], dry_run=dry_run)


def phase_1(dry_run: bool = False, skip_r: bool = False,
            preserve_dct: bool = False) -> bool:
    """Phase 1: VST + Residualization (R)"""
    log("PHASE 1: VST Normalization + Residualization")
    
    if skip_r:
        log("Skipping R-dependent residualization", "WARN")
        return True
    
    # Resolve CLR from current or latest run
    clr_path = find_artifact("deconvolution/music_segment_direct_proportions_CLR.csv")
    if clr_path is None:
        log("CLR file not found in any run directory", "ERROR")
        return False
    log(f"Using CLR: {clr_path}")

    phase1_dir = os.environ.get("RRRM_PHASE1_DIR")
    resid_args = [f"--clr={clr_path}", f"--outdir={phase1_dir}"]

    # Check for DCT marker panel (current run first, then latest)
    dct_panel = find_artifact("dct_markers/dct_marker_panel.txt")
    if dct_panel is None:
        # Also check legacy global location
        legacy = REPO_ROOT / "data/processed/dct_markers/dct_marker_panel.txt"
        if legacy.exists():
            dct_panel = legacy
    if dct_panel is not None:
        resid_args.append(f"--force_keep={dct_panel}")
        log(f"Passing DCT marker panel to residualization: {dct_panel}")

    if preserve_dct:
        resid_args.append("--preserve_dct")
        log("Preserving DCT proportions (NOT regressing out)")
    
    success = run_rscript("src/preprocessing/residualization.R", args=resid_args, dry_run=dry_run)
    if not success:
        return False

    # Export to Python format — pass the output dir so it writes into the run
    if not run_rscript("src/preprocessing/export_phase1.R",
                       args=[f"--outdir={phase1_dir}"], dry_run=dry_run):
        return False

    # Build Ensembl→Symbol ID map early (needed by gene_set_loader for YAML resolution)
    build_id_map_if_needed(dry_run)
    return True


def phase_1_5(dry_run: bool = False) -> bool:
    """Phase 1.5: Dataset-Derived DCT Marker Discovery"""
    log("PHASE 1.5: DCT Marker Discovery")

    # Static input (not pipeline-generated)
    vst_path = REPO_ROOT / "data/processed/vst_normalized" / "GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
    if not vst_path.exists():
        log(f"VST file not found: {vst_path}", "WARN")
        return False

    # Resolve CLR from current or latest run
    clr_path = find_artifact("deconvolution/music_segment_direct_proportions_CLR.csv")
    if clr_path is None:
        log("CLR file not found in any run directory", "WARN")
        return False
    log(f"Using CLR: {clr_path}")

    # Resolve metadata from current or latest run
    meta_path = find_artifact("phase1_residuals/meta_phase1.tsv.gz")
    if meta_path is None:
        log("meta_phase1.tsv.gz not found in any run directory", "WARN")
        return False
    log(f"Using metadata: {meta_path}")

    # Output goes into current run
    outdir = os.environ.get("RRRM_DCT_DIR")

    return run_python("src.markers.discover_dct", [
        f"--vst={vst_path}",
        f"--meta={meta_path}",
        f"--clr={clr_path}",
        f"--outdir={outdir}",
        "--tech_cols=LibraryBatch,ReadDepth,rRNA",
    ], dry_run=dry_run)


def phase_1_5b(dry_run: bool = False) -> bool:
    """Phase 1.5b: All-Segment Marker Discovery"""
    log("PHASE 1.5b: All-Segment Marker Discovery")

    vst_path = REPO_ROOT / "data/processed/vst_normalized" / "GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
    if not vst_path.exists():
        log(f"VST file not found: {vst_path}", "WARN")
        return False

    clr_path = find_artifact("deconvolution/music_segment_direct_proportions_CLR.csv")
    if clr_path is None:
        log("CLR file not found", "WARN")
        return False

    meta_path = find_artifact("phase1_residuals/meta_phase1.tsv.gz")
    if meta_path is None:
        log("meta_phase1.tsv.gz not found", "WARN")
        return False

    results_dir = os.environ.get("RRRM_RESULTS_DIR")
    outdir = Path(results_dir) / "segment_markers"

    return run_python("src.markers.discover_markers", [
        f"--vst={vst_path}",
        f"--meta={meta_path}",
        f"--clr={clr_path}",
        f"--outdir={outdir}",
    ], dry_run=dry_run)


def phase_2(dry_run: bool = False, skip_r: bool = False,
            max_genes: int = 2500, topk: int = 80,
            pool_controls: bool = False) -> bool:
    """Phase 2: Network Construction"""
    log("PHASE 2: Network Skeleton + LIONESS + Edge Regression")

    networks_dir = os.environ.get("RRRM_NETWORKS_DIR")

    # Check if DCT marker panel exists (current run first, then latest)
    dct_panel = find_artifact("dct_markers/dct_marker_panel.txt")
    force_args = []
    if dct_panel is not None:
        force_args = [f"--force_include={dct_panel}"]
        log(f"Including DCT marker panel: {dct_panel}")
    else:
        log("No DCT marker panel found; using variance-only gene selection", "WARN")

    if not run_python("src.networks.shared_topology",
                      [f"--max_genes={max_genes}", f"--topk={topk}",
                       f"--outdir={networks_dir}",
                       f"--id_map=data/processed/resources/id_map.tsv",
                       f"--biotype_filter=protein_coding"] + force_args,
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
        reg_args = [f"--phase2_dir={networks_dir}",
                    f"--outdir={networks_dir}/regression"]
        if pool_controls:
            reg_args.append("--pool-controls")
            log("Edge regression: pooling GC+VIV+BSL as ground controls")
        if not run_python("src.networks.edge_regression", reg_args, dry_run=dry_run):
            return False

    return True


def phase_3(dry_run: bool = False, num_seeds: int = 10,
            pool_controls: bool = False) -> bool:
    """Phase 3: Embeddings + Procrustes"""
    log("PHASE 3: Node2Vec Embeddings + Procrustes Alignment")

    networks_dir = os.environ.get("RRRM_NETWORKS_DIR")
    results_dir = os.environ.get("RRRM_RESULTS_DIR")

    # Step 3.2: Node2Vec embeddings
    # Include GND patterns if pooling controls
    emb_patterns = "Pred_*_FLT_z_hat.npy,Pred_*_GC_z_hat.npy"
    if pool_controls:
        emb_patterns = "Pred_*_FLT_z_hat.npy,Pred_*_GND_z_hat.npy,Pred_*_GC_z_hat.npy"
        log("Including GND (pooled ground) networks for embedding")

    if not run_python("src.networks.embeddings",
                      [f"--phase2_dir={networks_dir}",
                       f"--reg_dir={networks_dir}/regression",
                       f"--outdir={results_dir}/phase3_embeddings",
                       f"--num_seeds={num_seeds}",
                       f"--patterns={emb_patterns}"],
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

    results_dir = os.environ.get("RRRM_RESULTS_DIR")

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

    # Auto-detect available contrasts (supports both _GC and _GND suffixes)
    rewiring_dir = Path(results_dir) / "phase3_rewiring"
    contrast_files = sorted(rewiring_dir.glob("*_FLT_minus_*_rewiring_agg.tsv"))
    if not contrast_files:
        log("No rewiring agg files found for silent shifter analysis", "WARN")
        return True

    for rewiring_file in contrast_files:
        contrast = rewiring_file.stem.replace("_rewiring_agg", "")
        if not run_python("src.statistics.silent_shifters",
                          [f"--rewiring={rewiring_file}",
                           f"--outdir={results_dir}/phase5_silent_shifters_strict"] + reg_args,
                          dry_run=dry_run):
            log(f"Warning: silent_shifters failed for {contrast}", "WARN")

    return True


def phase_6(dry_run: bool = False, skip_r: bool = False,
            pool_controls: bool = False,
            focused_permutation: bool = False,
            pathway_permutation: bool = False) -> bool:
    """Phase 6: Uncertainty Estimation + Full Regression"""
    log("PHASE 6: Permutation + Bootstrap + Full Regression")

    results_dir = os.environ.get("RRRM_RESULTS_DIR")
    networks_dir = os.environ.get("RRRM_NETWORKS_DIR")

    # Permutation and bootstrap
    perm_args = [f"--phase2_dir={networks_dir}",
                 f"--outdir={results_dir}/phase6_uncertainty"]
    if pool_controls:
        perm_args.append("--pool-controls")
        log("Permutation test: pooling GC+VIV+BSL as ground controls")
    if focused_permutation:
        perm_args.append("--focused-permutation")
        # Increase permutation budget for focused testing to push p-values
        # below the BH threshold for the top-decile (~250 genes).
        # With K=50,000, min achievable p = 1/50001 ≈ 2×10⁻⁵,
        # vs BH threshold for rank-1: 0.05/250 = 2×10⁻⁴ — comfortably reachable.
        perm_args.append("--K_perm=50000")
        log("Permutation test: focused BH on top-decile genes (K=50,000)")
    if pathway_permutation:
        perm_args.append("--pathway-permutation")
        # Auto-detect gene map for pathway symbol resolution
        map_path = REPO_ROOT / "data/processed/resources" / "id_map.tsv"
        if map_path.exists():
            perm_args.append(f"--gene-map={map_path}")
        log("Permutation test: pathway-level competitive enrichment")
    if not run_python("src.statistics.permutation_bootstrap", perm_args, dry_run=dry_run):
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

    results_dir = os.environ.get("RRRM_RESULTS_DIR")
    figures_dir = Path(results_dir) / "figures"

    # 9a — Existing diagnostic / QC plots (publication_plots)
    if not run_python("src.visualization.publication_plots",
                      [f"--results_dir={results_dir}",
                       f"--out_dir={figures_dir}"],
                      dry_run=dry_run):
        log("publication_plots module failed (non-fatal, continuing)", "WARN")

    # 9b — New key publication figures (rewiring-vs-DE, focused perm, pathway,
    #       age comparison, arm comparison, pipeline schematic, retinol subnetwork)
    log("PHASE 9b: Key Publication Figures")
    if not run_python("src.visualization.publication_figures",
                      [f"--results_dir={results_dir}"],
                      dry_run=dry_run):
        log("publication_figures module failed (non-fatal)", "WARN")

    pub_dir = Path(results_dir) / "figures" / "publication"
    log(f"Publication figures → {pub_dir}")
    return True


def phase_8(dry_run: bool = False, max_genes: int = 2500, topk: int = 80) -> bool:
    """Phase 8: Leakage-Safe Cross-Validation"""
    log("PHASE 8: Leakage-Safe Cross-Validation")

    results_dir = os.environ.get("RRRM_RESULTS_DIR")
    networks_dir = os.environ.get("RRRM_NETWORKS_DIR")

    meta_path = find_artifact("phase1_residuals/meta_phase1.tsv.gz")
    rtech_path = find_artifact("phase1_residuals/Rtech.tsv.gz")

    if meta_path is None or rtech_path is None:
        log("Phase 1 outputs not found; cannot run validation", "WARN")
        return False

    return run_python("src.validation.cross_validation", [
        f"--phase2_dir={networks_dir}",
        f"--rtech={rtech_path}",
        f"--meta={meta_path}",
        f"--outdir={results_dir}/phase8_validation",
        f"--max_genes={max_genes}",
        f"--topk={topk}",
    ], dry_run=dry_run)


def phase_8b(dry_run: bool = False, max_genes: int = 2500, topk: int = 80) -> bool:
    """Phase 8b: Enhanced Predictive Validation

    Three improvements over Phase 8:
      1. Stratum-specific LOO-CV (n=10 per stratum, 5 FLT + 5 GC)
      2. Expression-based classification as control baseline
      3. Permutation null distribution (1,000 shuffles) with p-values
    """
    log("PHASE 8b: Enhanced Predictive Validation")

    results_dir = os.environ.get("RRRM_RESULTS_DIR")
    networks_dir = os.environ.get("RRRM_NETWORKS_DIR")

    meta_path = find_artifact("phase1_residuals/meta_phase1.tsv.gz")
    rtech_path = find_artifact("phase1_residuals/Rtech.tsv.gz")

    if meta_path is None or rtech_path is None:
        log("Phase 1 outputs not found; cannot run enhanced validation", "WARN")
        return False

    return run_python("src.validation.enhanced_cv", [
        f"--phase2_dir={networks_dir}",
        f"--rtech={rtech_path}",
        f"--meta={meta_path}",
        f"--outdir={results_dir}/phase8_validation",
        f"--max_genes={max_genes}",
        f"--topk={topk}",
        "--n_perms=1000",
    ], dry_run=dry_run)


def phase_7(dry_run: bool = False) -> bool:
    """Phase 7: Biological Grounding"""
    log("PHASE 7: Biological Grounding + Enrichment")

    results_dir = os.environ.get("RRRM_RESULTS_DIR")

    # Auto-detect available contrasts (supports both _GC and _GND suffixes)
    rewiring_dir = Path(results_dir) / "phase3_rewiring"
    contrast_files = sorted(rewiring_dir.glob("*_FLT_minus_*_rewiring_agg.tsv"))
    if not contrast_files:
        log("No rewiring agg files found for Phase 7", "WARN")
        return True

    contrasts = [f.stem.replace("_rewiring_agg", "") for f in contrast_files]
    log(f"Detected contrasts: {contrasts}")

    # Ensure ID map exists (may have been built by Phase 1, but rebuild if stale)
    build_id_map_if_needed(dry_run)

    map_path = REPO_ROOT / "data/processed/resources" / "id_map.tsv"

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
    parser.add_argument("--pool-controls", action="store_true",
                        help="Pool GC+VIV+BSL as ground reference (triples control n)")
    parser.add_argument("--preserve-dct", action="store_true",
                        help="Do NOT regress out DCT proportions (preserves DCT signal)")
    parser.add_argument("--focused-permutation", action="store_true",
                        help="Run focused BH correction on top-decile genes only (~250). "
                             "Reduces the multiple-testing burden from 2500 to ~250 hypotheses.")
    parser.add_argument("--pathway-permutation", action="store_true",
                        help="Run pathway-level permutation testing (competitive GSEA-style). "
                             "Tests ~100 pathways instead of ~2500 individual genes.")
    
    args = parser.parse_args()
    
    os.chdir(REPO_ROOT)
    
    # Determine which phases to run first (needed for init_run)
    # Note: Phase 9 (figures) runs last, after all data phases
    phases_available = [0, 1, 1.5, 2, 3, 5, 6, 7, 8, 9]
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
    # All outputs live under data/results/<run_id>/
    run_root = str(REPO_ROOT / "data/results" / run_id)
    os.environ["RRRM_RUN_ID"] = run_id
    os.environ["RRRM_RESULTS_DIR"] = run_root
    os.environ["RRRM_DECONV_DIR"] = str(Path(run_root) / "deconvolution")
    os.environ["RRRM_PHASE1_DIR"] = str(Path(run_root) / "phase1_residuals")
    os.environ["RRRM_DCT_DIR"] = str(Path(run_root) / "dct_markers")
    os.environ["RRRM_NETWORKS_DIR"] = str(Path(run_root) / "networks")
    
    phases = {
        0: lambda: phase_0(args.dry_run, args.skip_r),
        1: lambda: phase_1(args.dry_run, args.skip_r, args.preserve_dct),
        1.5: lambda: phase_1_5(args.dry_run),
        1.6: lambda: phase_1_5b(args.dry_run),
        2: lambda: phase_2(args.dry_run, args.skip_r, args.max_genes, args.topk, args.pool_controls),
        3: lambda: phase_3(args.dry_run, args.num_seeds, args.pool_controls),
        5: lambda: phase_5(args.dry_run),
        6: lambda: phase_6(args.dry_run, args.skip_r, args.pool_controls,
                           args.focused_permutation, args.pathway_permutation),
        7: lambda: phase_7(args.dry_run),
        8: lambda: phase_8(args.dry_run, args.max_genes, args.topk),
        8.5: lambda: phase_8b(args.dry_run, args.max_genes, args.topk),
        9: lambda: phase_9(args.dry_run),  # Figure generation runs last
    }
    
    log(f"RRRM-2 Pipeline Runner", "INFO")
    log(f"Phases to run: {to_run}")
    log(f"Max genes: {args.max_genes}")
    log(f"Skip R steps: {args.skip_r}")
    log(f"Dry run: {args.dry_run}")
    print()
    
    # Run phases in dependency-aware order (Phase 6 BEFORE Phase 5 for regression support)
    execution_order = [0, 1, 1.5, 2, 3, 6, 5, 7, 8, 9]
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

