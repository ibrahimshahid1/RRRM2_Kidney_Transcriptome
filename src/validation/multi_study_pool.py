# src/validation/multi_study_pool.py
"""
Pre-registered OSD-102 + OSD-771 LAR-Young pooling.

This is separate from independent external replication. Pooling is restricted to
RRRM-2/OSD-771 LAR-Young and individual OSD-102 FLT/GC mice. ComBat-seq with
study as batch and treatment preserved is required before any pooled-network
claim. PCA checks must show that study separation is reduced while treatment
signal is not erased; otherwise the module recommends fixed-effects or
meta-analysis and forbids pooled network validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression

from src.common import REPO_ROOT, find_sample_col, normalize_labels


def enforce_lar_young_scope(meta: pd.DataFrame, study_col: str = "study") -> pd.DataFrame:
    meta = meta.copy()
    if study_col not in meta.columns:
        raise ValueError(f"Metadata must contain '{study_col}' for multi-study pooling")
    if "OSD-568" in set(meta[study_col].astype(str)):
        raise ValueError("OSD-568 is excluded and cannot enter multi-study pooling")
    allowed = meta[study_col].isin(["OSD-771", "OSD-102"])
    if {"Age", "Arm"}.issubset(meta.columns):
        lar_young = (meta["Age"].isin(["YNG", "Young"])) & (meta["Arm"].isin(["LAR", "LAR_T", "LAR-T"]))
        osd102 = meta[study_col] == "OSD-102"
        allowed = allowed & (lar_young | osd102)
    if not allowed.all():
        dropped = int((~allowed).sum())
        print(f"[WARN] Dropping {dropped} samples outside OSD-771 LAR-Young + OSD-102 scope")
    return meta.loc[allowed].copy()


def pca_batch_treatment_check(counts: pd.DataFrame, meta: pd.DataFrame, study_col: str, treatment_col: str) -> dict:
    common = [s for s in counts.columns if s in meta.index]
    X = np.log1p(counts[common].T.values.astype(float))
    pcs = PCA(n_components=min(5, X.shape[0] - 1, X.shape[1])).fit_transform(X)

    def pc1_accuracy(label_col: str) -> float:
        y = meta.loc[common, label_col].astype(str).values
        if len(set(y)) < 2:
            return float("nan")
        clf = LogisticRegression(max_iter=1000)
        clf.fit(pcs[:, : min(3, pcs.shape[1])], y)
        return float((clf.predict(pcs[:, : min(3, pcs.shape[1])]) == y).mean())

    return {
        "n_samples": len(common),
        "study_predictability_from_pcs": pc1_accuracy(study_col),
        "treatment_predictability_from_pcs": pc1_accuracy(treatment_col),
        "passes_pca_gate": False,
        "gate_note": "Run on post-ComBat-seq counts and compare to pre-correction baseline before enabling pooled claims.",
    }


def run_combat_seq_guard(counts_path: Path, meta_path: Path, outdir: Path) -> Path:
    """Guarded placeholder for ComBat-seq execution via R/sva."""
    combat_out = outdir / "combat_seq_counts.tsv"
    manifest = {
        "required_method": "ComBat-seq",
        "batch": "study",
        "preserve": "treatment/EnvGroup",
        "counts_path": str(counts_path),
        "meta_path": str(meta_path),
        "status": "not_run_by_python_guard",
        "message": "Run an R/sva ComBat-seq command or supply --combat_counts for PCA gating.",
    }
    (outdir / "combat_seq_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return combat_out


def main() -> None:
    ap = argparse.ArgumentParser(description="OSD-102 + OSD-771 LAR-Young pooling guard")
    ap.add_argument("--counts", required=True, help="Raw counts genes × samples for candidate pooling")
    ap.add_argument("--meta", required=True, help="Metadata with study and treatment columns")
    ap.add_argument("--combat_counts", default="", help="Post-ComBat-seq counts; if absent, only guard manifest is written")
    ap.add_argument("--study_col", default="study")
    ap.add_argument("--treatment_col", default="EnvGroup")
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/multi_study_pool"))
    args = ap.parse_args()

    counts = pd.read_csv(args.counts, sep=None, engine="python", index_col=0)
    meta = pd.read_csv(args.meta, sep=None, engine="python")
    sample_col = find_sample_col(meta)
    meta = normalize_labels(meta.set_index(sample_col, drop=False))
    meta = enforce_lar_young_scope(meta, args.study_col)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    meta.to_csv(outdir / "pooling_scope_samples.tsv", sep="\t", index=True)

    if args.combat_counts:
        combat_counts = pd.read_csv(args.combat_counts, sep=None, engine="python", index_col=0)
        check = pca_batch_treatment_check(combat_counts, meta, args.study_col, args.treatment_col)
        check["pooled_network_claim_allowed"] = bool(check["passes_pca_gate"])
        check["fallback_if_failed"] = "fixed_effects_or_meta_analysis"
        pd.DataFrame([check]).to_csv(outdir / "pca_pooling_gate.tsv", sep="\t", index=False)
        print(f"[OK] PCA pooling gate written to {outdir}")
    else:
        run_combat_seq_guard(Path(args.counts), Path(args.meta), outdir)
        result = {
            "pooled_network_claim_allowed": False,
            "reason": "ComBat-seq/PCA checks have not been run",
            "fallback": "fixed_effects_or_meta_analysis",
        }
        (outdir / "pooling_guard_result.json").write_text(json.dumps(result, indent=2) + "\n")
        print(f"[OK] Pooling guard written to {outdir}; no pooled claim allowed yet")


if __name__ == "__main__":
    main()
