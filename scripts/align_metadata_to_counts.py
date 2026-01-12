#!/usr/bin/env python3
"""
Metadata Alignment Script for GLDS-674 / OSD-771
Aligns metadata to count matrices (STAR, RSEM, VST) and creates aligned datasets.

This script:
1. Reads count files (STAR unnormalized, RSEM unnormalized, VST normalized)
2. Reads ISA metadata
3. Aligns metadata to samples across all count matrices
4. Optionally parses design factors from sample names
5. Saves aligned outputs to data/processed/aligned_outputs/

Usage:
    python scripts/align_metadata_to_counts.py
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import pandas as pd
import sys


# --- Paths Configuration ---
# All paths are relative to the repository root
REPO_ROOT = Path(__file__).parent.parent

STAR_COUNTS_PATH = REPO_ROOT / "data/raw/counts/GLDS-674_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv"
RSEM_COUNTS_PATH = REPO_ROOT / "data/raw/counts/GLDS-674_rna_seq_RSEM_Unnormalized_Counts_rRNArm_GLbulkRNAseq.csv"
VST_COUNTS_PATH  = REPO_ROOT / "data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"
META_PATH        = REPO_ROOT / "data/raw/metadata/a_OSD-771_transcription-profiling_rna-sequencing-(rna-seq)_Illumina.txt"

OUTPUT_DIR = REPO_ROOT / "data/processed/aligned_outputs"


# --- Helper Functions ---
def read_counts_csv(path: Path) -> pd.DataFrame:
    """
    Reads a GeneLab counts CSV shaped like:
      gene_id, sample1, sample2, ...
    Returns a DataFrame indexed by gene_id with sample columns.
    """
    if not path.exists():
        raise FileNotFoundError(f"Count file not found: {path}")
    
    df = pd.read_csv(path, dtype=str)
    gene_col = df.columns[0]  # usually "Unnamed: 0" or first column
    df = df.set_index(gene_col)

    # Enforce sample IDs as strings, keep original exact names
    df.columns = df.columns.astype(str)

    # Convert all entries to numeric (counts / VST can be float)
    df = df.apply(pd.to_numeric, errors="raise")
    return df


def read_metadata_txt(path: Path) -> pd.DataFrame:
    """
    Reads the ISA-like metadata TSV. Uses 'Sample Name' as the index.
    """
    if not path.exists():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    
    meta = pd.read_csv(path, sep="\t", dtype=str)
    if "Sample Name" not in meta.columns:
        raise ValueError("Expected a 'Sample Name' column in metadata file.")
    meta = meta.set_index("Sample Name")
    meta.index = meta.index.astype(str)
    return meta


def align_metadata_to_samples(meta: pd.DataFrame, canonical_samples: list[str]) -> pd.DataFrame:
    """
    Align metadata to the canonical sample ordering.
    Raises error if any samples are missing from metadata.
    """
    missing = [s for s in canonical_samples if s not in meta.index]
    if missing:
        raise ValueError(f"Metadata missing {len(missing)} canonical samples. Example: {missing[:5]}")
    meta_aligned = meta.loc[canonical_samples].copy()  # subset + reorder
    assert list(meta_aligned.index) == canonical_samples
    return meta_aligned


def align_counts_to_samples(counts: pd.DataFrame, canonical_samples: list[str], name: str) -> pd.DataFrame:
    """
    Align count matrix columns to the canonical sample ordering.
    Raises error if any samples are missing from count matrix.
    """
    missing = [s for s in canonical_samples if s not in counts.columns]
    if missing:
        raise ValueError(f"{name} counts missing {len(missing)} canonical samples. Example: {missing[:5]}")

    # Subset + reorder columns
    aligned = counts.loc[:, canonical_samples].copy()
    assert list(aligned.columns) == canonical_samples
    return aligned


def parse_design_from_sample_name(sample_name: str) -> dict:
    """
    Optional: extract simple design factors from the RRRM2 sample name.
    Example: RRRM2_R-KDN_BSL_ISS-T_YNG_BY1
    
    Returns dict with keys: EnvGroup, Arm, Age, AnimalCode
    """
    parts = sample_name.split("_")
    env_levels = {"FLT", "HGC", "VIV", "BSL", "GC", "VGC"}
    arm_levels = {"ISS-T", "LAR"}
    age_levels = {"YNG", "OLD"}

    env = next((p for p in parts if p in env_levels), None)
    arm = next((p for p in parts if p in arm_levels), None)
    age = next((p for p in parts if p in age_levels), None)
    animal = parts[-1] if parts else None

    return {"EnvGroup": env, "Arm": arm, "Age": age, "AnimalCode": animal}


@dataclass
class DatasetBundle:
    """Container for aligned count matrix and metadata."""
    counts: pd.DataFrame        # genes × samples
    meta: pd.DataFrame          # samples × metadata


# --- Main Alignment Logic ---
def main():
    print("=" * 70)
    print("GLDS-674 Metadata Alignment to Count Matrices")
    print("=" * 70)
    
    # 1. Read RSEM counts (canonical reference)
    print(f"\n[1/5] Reading RSEM counts (canonical reference)...")
    print(f"      Path: {RSEM_COUNTS_PATH}")
    rsem = read_counts_csv(RSEM_COUNTS_PATH)
    canonical_samples = list(rsem.columns)
    print(f"      ✓ Loaded {rsem.shape[0]} genes × {rsem.shape[1]} samples")
    
    # 2. Read metadata
    print(f"\n[2/5] Reading ISA metadata...")
    print(f"      Path: {META_PATH}")
    meta = read_metadata_txt(META_PATH)
    print(f"      ✓ Loaded metadata for {len(meta)} samples with {len(meta.columns)} columns")
    
    # 3. Align metadata to canonical samples
    print(f"\n[3/5] Aligning metadata to canonical sample ordering...")
    meta_aligned = align_metadata_to_samples(meta, canonical_samples)
    
    # Optional: add parsed design columns (useful later in Python modeling)
    print(f"      Parsing design factors from sample names...")
    parsed = pd.DataFrame([parse_design_from_sample_name(s) for s in meta_aligned.index], 
                          index=meta_aligned.index)
    meta_aligned = pd.concat([meta_aligned, parsed], axis=1)
    print(f"      ✓ Metadata aligned with {len(parsed.columns)} additional design columns")
    
    # 4. Read and align other count matrices
    print(f"\n[4/5] Reading and aligning count matrices...")
    
    print(f"      - STAR unnormalized: {STAR_COUNTS_PATH}")
    star = read_counts_csv(STAR_COUNTS_PATH)
    star_aligned = align_counts_to_samples(star, canonical_samples, name="STAR")
    print(f"        ✓ {star_aligned.shape[0]} genes × {star_aligned.shape[1]} samples")
    
    print(f"      - VST normalized: {VST_COUNTS_PATH}")
    vst = read_counts_csv(VST_COUNTS_PATH)
    vst_aligned = align_counts_to_samples(vst, canonical_samples, name="VST_rRNArm")
    print(f"        ✓ {vst_aligned.shape[0]} genes × {vst_aligned.shape[1]} samples")
    
    # RSEM is already canonical, but align for consistency
    rsem_aligned = align_counts_to_samples(rsem, canonical_samples, name="RSEM_rRNArm")
    
    # 5. Verify alignment consistency
    print(f"\n[5/5] Verifying alignment consistency...")
    assert list(star_aligned.columns) == list(rsem_aligned.columns) == list(vst_aligned.columns), \
        "Sample ordering mismatch between count matrices!"
    assert list(meta_aligned.index) == list(rsem_aligned.columns), \
        "Sample ordering mismatch between metadata and counts!"
    print(f" All datasets aligned to {len(canonical_samples)} canonical samples")
    
    # Create dataset bundles
    bundles = {
        "star_raw": DatasetBundle(counts=star_aligned, meta=meta_aligned),
        "rsem_rRNArm_raw": DatasetBundle(counts=rsem_aligned, meta=meta_aligned),
        "vst_rRNArm": DatasetBundle(counts=vst_aligned, meta=meta_aligned),
    }
    
    print(f"\n      Created {len(bundles)} aligned dataset bundles:")
    for name in bundles.keys():
        print(f"    - {name}")
    
    # --- Save aligned outputs ---
    print(f"\n{'=' * 70}")
    print("Saving Aligned Outputs")
    print(f"{'=' * 70}")
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Save metadata once (since it's shared across all datasets)
    meta_out = OUTPUT_DIR / "metadata_aligned.tsv"
    meta_aligned.to_csv(meta_out, sep="\t")
    print(f"  ✓ Metadata: {meta_out}")
    
    # Save each count matrix separately
    for name, bundle in bundles.items():
        counts_out = OUTPUT_DIR / f"{name}_counts.csv"
        bundle.counts.to_csv(counts_out)
        print(f"  ✓ {name:20s}: {counts_out}")
    
    print(f"\n{'=' * 70}")
    print("✅ Alignment Complete!")
    print(f"{'=' * 70}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Canonical sample count: {len(canonical_samples)}")
    print()


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n❌ Error: {e}", file=sys.stderr)
        sys.exit(1)
