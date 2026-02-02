# scripts/phase4_anchor_qc_report.py
"""
Phase 4 (prereg documentation): Anchor QC report.

Goal: produce a small report showing anchors are "stable" across groups
and optionally not-DE (if you provide a gene-level DE table).

Inputs:
  - anchors.txt from Phase 3.3
  - meta_phase1.tsv.gz
  - OPTIONAL: gene_DE.tsv with columns: gene, log2FC(or logFC), FDR(or adj.P.Val)

Outputs:
  data/results/phase4_anchor_qc/
    anchors_preregistered.txt
    anchor_qc.tsv
    anchor_qc_summary.json
"""

from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]


def find_sample_col(meta: pd.DataFrame) -> str:
    """Find the sample identifier column in metadata."""
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            return col
    return meta.columns[0]


def normalize_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Normalize Age, Arm, EnvGroup labels for consistency."""
    meta = meta.copy()
    meta["Age"] = meta["Age"].astype(str).replace({
        "Young": "YNG", "Yng": "YNG", "young": "YNG",
        "Old": "OLD", "old": "OLD"
    })
    meta["Arm"] = meta["Arm"].astype(str).replace({
        "ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T", "ISS T": "ISS-T"
    })
    meta["Arm"] = meta["Arm"].astype(str).replace({
        "LAR_T": "LAR", "LAR-T": "LAR", "LAR T": "LAR"
    })
    meta["EnvGroup"] = meta["EnvGroup"].astype(str).replace({
        "HGC": "GC", "VGC": "VIV", "HGC/GC": "GC", "VIV/VGC": "VIV"
    })
    return meta


def main():
    ap = argparse.ArgumentParser(description="Generate anchor QC report for prereg")
    ap.add_argument("--anchors", 
                    default=str(REPO_ROOT / "data/results/phase3_rewiring/anchors.txt"),
                    help="Path to anchors.txt from Phase 3.3")
    ap.add_argument("--meta", 
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"),
                    help="Metadata file")
    ap.add_argument("--expr_cell_means", default="", 
                    help="Optional TSV with columns: gene + 8 cell-mean columns")
    ap.add_argument("--gene_de", default="", 
                    help="Optional TSV with columns: gene, logFC/log2FC, FDR/adj.P.Val")
    ap.add_argument("--outdir", 
                    default=str(REPO_ROOT / "data/results/phase4_anchor_qc"),
                    help="Output directory")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load anchors
    anchors_path = Path(args.anchors)
    if not anchors_path.exists():
        raise FileNotFoundError(f"Anchors file not found: {anchors_path}")
    
    anchors = [x.strip() for x in anchors_path.read_text().splitlines() if x.strip()]
    print(f"Loaded {len(anchors)} anchors from {anchors_path}")

    # Write preregistered anchors (canonical copy)
    prereg_path = outdir / "anchors_preregistered.txt"
    prereg_path.write_text("\n".join(anchors) + "\n")
    print(f"Wrote {prereg_path}")

    # Load metadata
    meta_path = Path(args.meta)
    if meta_path.exists():
        meta = pd.read_csv(meta_path, sep="\t", compression="gzip")
        meta = meta.set_index(find_sample_col(meta), drop=False)
        meta = normalize_labels(meta)
        meta["cell"] = (meta["Age"].astype(str) + "_" + 
                        meta["Arm"].astype(str) + "_" + 
                        meta["EnvGroup"].astype(str))
        print(f"Loaded metadata with {len(meta)} samples")
    else:
        print(f"[WARN] Metadata not found: {meta_path}")
        meta = None

    # Build QC table
    qc = pd.DataFrame({"gene": anchors})

    # --- Stability metrics from expression cell means ---
    if args.expr_cell_means and Path(args.expr_cell_means).exists():
        cm = pd.read_csv(args.expr_cell_means, sep="\t")
        if "gene" not in cm.columns:
            raise ValueError("expr_cell_means must have a 'gene' column")
        cm = cm.set_index("gene")
        cell_cols = [c for c in cm.columns if c.startswith(("YNG_", "OLD_"))]
        if len(cell_cols) < 4:
            raise ValueError(f"expr_cell_means doesn't look like cell means. Found cols: {cm.columns[:10]}")
        
        # Filter to anchors present in cell means
        valid_anchors = [a for a in anchors if a in cm.index]
        stab = cm.loc[valid_anchors, cell_cols].copy()
        
        # Compute stability metrics
        between_cell_sd = stab.std(axis=1)
        between_cell_range = stab.max(axis=1) - stab.min(axis=1)
        between_cell_absmean = stab.abs().mean(axis=1)
        
        qc = qc.merge(pd.DataFrame({
            "gene": valid_anchors,
            "between_cell_sd": between_cell_sd.values,
            "between_cell_range": between_cell_range.values,
            "between_cell_absmean": between_cell_absmean.values,
        }), on="gene", how="left")
        print(f"Added stability metrics for {len(valid_anchors)} anchors")
    else:
        qc["between_cell_sd"] = np.nan
        qc["between_cell_range"] = np.nan
        qc["between_cell_absmean"] = np.nan
        print("[INFO] No expr_cell_means provided, stability metrics set to NA")

    # --- Optional: DE filters ---
    if args.gene_de and Path(args.gene_de).exists():
        de = pd.read_csv(args.gene_de, sep=None, engine="python")
        if "gene" not in de.columns:
            raise ValueError("gene_de must have a 'gene' column")
        
        # Normalize column names
        de = de.rename(columns={
            "adj.P.Val": "FDR", "padj": "FDR", "qvalue": "FDR",
            "log2FC": "logFC"
        })
        
        if "FDR" not in de.columns:
            raise ValueError(f"gene_de missing FDR/adj.P.Val column. Columns: {list(de.columns)}")
        if "logFC" not in de.columns:
            raise ValueError(f"gene_de missing logFC/log2FC column. Columns: {list(de.columns)}")
        
        de = de[["gene", "logFC", "FDR"]].drop_duplicates("gene")
        qc = qc.merge(de, on="gene", how="left")
        qc["passes_relaxed_DE"] = (qc["FDR"] > 0.2) & (qc["logFC"].abs() < 0.3)
        print(f"Added DE metrics from {args.gene_de}")
    else:
        qc["logFC"] = np.nan
        qc["FDR"] = np.nan
        qc["passes_relaxed_DE"] = np.nan
        print("[INFO] No gene_de provided, DE metrics set to NA")

    # Write QC table
    qc_path = outdir / "anchor_qc.tsv"
    qc.to_csv(qc_path, sep="\t", index=False)
    print(f"Wrote {qc_path}")

    # Write summary JSON
    summary = {
        "n_anchors": int(len(anchors)),
        "n_with_DE": int(qc["FDR"].notna().sum()),
        "n_pass_relaxed_DE": int((qc["passes_relaxed_DE"] == True).sum()) if qc["passes_relaxed_DE"].notna().any() else None,
        "n_with_stability_metrics": int(qc["between_cell_sd"].notna().sum()),
        "note": "If expr_cell_means was not provided, stability metrics are NA. Provide a gene x cellMeans TSV to fill them."
    }
    summary_path = outdir / "anchor_qc_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path}")

    print(f"\n[OK] Anchor QC outputs written to: {outdir}")


if __name__ == "__main__":
    main()
