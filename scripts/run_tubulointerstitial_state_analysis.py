#!/usr/bin/env python3
"""Run tubulointerstitial state-space analysis for manuscript v4."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

SCRIPT_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(SCRIPT_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(SCRIPT_REPO_ROOT))

from src.common import REPO_ROOT
from src.networks.tubulointerstitial_state import (
    MATRIX_COMPONENTS,
    DCT_COMPONENT,
    IMMUNE_COMPONENTS,
    PRESERVATION_COMPONENT,
    TARGET_GENE_PANELS,
    build_component_correlations,
    build_group_effects,
    build_rrrm2_deconvolution_overlay,
    build_state_space_scores,
    build_targeted_gene_panel,
)


DEFAULT_MECH = REPO_ROOT / "data/results/run_20260519_000547_2500g/contrast_vectors/mechanism_axis"
DEFAULT_DECONV = REPO_ROOT / "data/results/run_20260517_213205_2500g/deconvolution"


def file_sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def resolve(path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else REPO_ROOT / p


def read_expression(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    expr = pd.read_csv(path, sep="\t", compression=compression, index_col=0)
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    return expr


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Tubulointerstitial state-space analysis")
    ap.add_argument(
        "--mechanism-score-matrix",
        default=str(DEFAULT_MECH / "mechanism_score_matrix.tsv"),
        help="Mechanism score matrix from the mechanism-axis phase.",
    )
    ap.add_argument(
        "--deconv-dir",
        default=str(DEFAULT_DECONV),
        help="RRRM-2 deconvolution output directory.",
    )
    ap.add_argument(
        "--rtech",
        default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"),
        help="RRRM-2 residualized expression matrix.",
    )
    ap.add_argument(
        "--id-map",
        default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"),
        help="Ensembl-to-MGI symbol map.",
    )
    ap.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Defaults to mechanism_score_matrix parent / tubulointerstitial_state.",
    )
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    mechanism_path = resolve(args.mechanism_score_matrix)
    deconv_dir = resolve(args.deconv_dir)
    rtech_path = resolve(args.rtech)
    id_map_path = resolve(args.id_map)
    outdir = resolve(args.outdir) if args.outdir else mechanism_path.parent / "tubulointerstitial_state"
    outdir.mkdir(parents=True, exist_ok=True)

    mechanism = pd.read_csv(mechanism_path, sep="\t")
    state_scores = build_state_space_scores(mechanism)
    state_scores.to_csv(outdir / "state_space_scores.tsv", sep="\t", index=False)

    group_effects = build_group_effects(state_scores)
    group_effects.to_csv(outdir / "state_space_group_effects.tsv", sep="\t", index=False)

    cluster_path = deconv_dir / "music_cluster_proportions.csv"
    clr_path = deconv_dir / "music_segment_direct_proportions_CLR.csv"
    cluster = pd.read_csv(cluster_path)
    clr = pd.read_csv(clr_path) if clr_path.exists() else pd.DataFrame()
    overlay = build_rrrm2_deconvolution_overlay(state_scores, cluster, clr)
    overlay.to_csv(outdir / "rrrm2_deconvolution_overlay.tsv", sep="\t", index=False)

    correlations = build_component_correlations(overlay)
    correlations.to_csv(outdir / "state_space_component_correlations.tsv", sep="\t", index=False)

    expression = read_expression(rtech_path)
    rrrm2_meta = state_scores[state_scores["study"].eq("RRRM-2")].copy()
    gene_panel = build_targeted_gene_panel(expression, str(id_map_path), rrrm2_meta, TARGET_GENE_PANELS)
    gene_panel.to_csv(outdir / "targeted_renal_gene_panel.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "tubulointerstitial state-space analysis",
        "timestamp": datetime.now().isoformat(),
        "git_commit": git_commit(),
        "component_definitions": {
            "matrix_component": {
                "source_scores": list(MATRIX_COMPONENTS),
                "definition": "mean within-study/scenario z-score",
                "direction": "higher = ECM/adhesion/proteolysis higher",
            },
            "dct_transport_component": {
                "source_scores": [DCT_COMPONENT],
                "definition": "within-study/scenario z-score",
                "direction": "higher = DCT/NCC-WNK transport higher; native direction preserved",
            },
            "immune_context_component": {
                "source_scores": list(IMMUNE_COMPONENTS),
                "definition": "mean within-study/scenario z-score",
                "direction": "higher = TLR4/innate and macrophage/inflammatory context higher",
            },
            "preservation_stress_component": {
                "source_scores": [PRESERVATION_COMPONENT],
                "definition": "within-study/scenario z-score",
                "direction": "higher = preservation/sample-stress signature higher",
            },
            "matrix_minus_dct": {
                "definition": "matrix_component - dct_transport_component",
                "direction": "higher = matrix-high / DCT-low state",
            },
        },
        "state_space_interpretation": {
            "x_axis": "matrix_component",
            "y_axis": "dct_transport_component",
            "remodeling_direction": "rightward and downward",
        },
        "inputs": {
            "mechanism_score_matrix": str(mechanism_path),
            "mechanism_score_matrix_sha256": file_sha256(mechanism_path),
            "deconvolution_cluster_proportions": str(cluster_path),
            "deconvolution_cluster_proportions_sha256": file_sha256(cluster_path),
            "deconvolution_clr_proportions": str(clr_path),
            "deconvolution_clr_proportions_sha256": file_sha256(clr_path),
            "rtech": str(rtech_path),
            "rtech_sha256": file_sha256(rtech_path),
            "id_map": str(id_map_path),
            "id_map_sha256": file_sha256(id_map_path),
        },
        "outputs": {
            "state_space_scores": str(outdir / "state_space_scores.tsv"),
            "state_space_group_effects": str(outdir / "state_space_group_effects.tsv"),
            "state_space_component_correlations": str(outdir / "state_space_component_correlations.tsv"),
            "rrrm2_deconvolution_overlay": str(outdir / "rrrm2_deconvolution_overlay.tsv"),
            "targeted_renal_gene_panel": str(outdir / "targeted_renal_gene_panel.tsv"),
        },
        "counts": {
            "state_space_rows": int(len(state_scores)),
            "group_effect_rows": int(len(group_effects)),
            "rrrm2_deconvolution_rows": int(len(overlay)),
            "component_correlation_rows": int(len(correlations)),
            "targeted_gene_panel_rows": int(len(gene_panel)),
            "targeted_gene_panel_scored_genes": int(gene_panel.loc[gene_panel["status"].eq("scored"), "gene"].nunique()) if not gene_panel.empty else 0,
        },
    }
    (outdir / "state_space_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"[state] wrote {outdir}")


if __name__ == "__main__":
    main()
