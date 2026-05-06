# src/validation/continuous_target.py
"""
Fold-safe continuous-target validation for kidney stress/injury scores.

Targets are built from independent markers such as Havcr1/KIM-1 and Lcn2/NGAL.
Those marker genes are excluded from predictor features to avoid circularity.
Network and expression baselines are evaluated inside cross-validation folds by
Pearson and Spearman correlation.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

from src.common import REPO_ROOT, find_sample_col, id_map_lookup, normalize_labels
from src.validation.sample_features import node_strength

DEFAULT_MARKERS = ["Havcr1", "Lcn2", "Spp1", "Timp1", "Vcam1"]


def resolve_markers(markers: list[str], id_map: Path, genes: set[str]) -> list[str]:
    _, symbol_to_ens = id_map_lookup(id_map)
    resolved: list[str] = []
    for marker in markers:
        if marker in genes:
            resolved.append(marker)
        else:
            resolved.extend(sorted(symbol_to_ens.get(marker.lower(), set()) & genes))
    return sorted(set(resolved))


def zscore_rows(X: np.ndarray) -> np.ndarray:
    mu = X.mean(axis=1, keepdims=True)
    sd = X.std(axis=1, keepdims=True) + 1e-8
    return (X - mu) / sd


def build_stress_score(rtech: pd.DataFrame, marker_ids: list[str]) -> np.ndarray:
    if not marker_ids:
        raise ValueError("No stress marker IDs resolved in expression matrix")
    marker_expr = rtech.loc[marker_ids].values.astype(float)
    return zscore_rows(marker_expr).mean(axis=0)


def fold_safe_regression_metrics(X: np.ndarray, y: np.ndarray, n_splits: int, seed: int) -> dict:
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    pred = np.empty_like(y, dtype=float)
    for train_idx, test_idx in kf.split(X):
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X[train_idx])
        X_test = scaler.transform(X[test_idx])
        model = Ridge(alpha=1.0)
        model.fit(X_train, y[train_idx])
        pred[test_idx] = model.predict(X_test)
    pear = pearsonr(y, pred)[0] if np.std(pred) > 0 and np.std(y) > 0 else np.nan
    spear = spearmanr(y, pred).correlation if np.std(pred) > 0 and np.std(y) > 0 else np.nan
    return {"pearson": float(pear), "spearman": float(spear)}


def build_fold_features(
    feature_type: str,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    strength: np.ndarray,
    lioness: np.ndarray,
    expr: np.ndarray,
    n_strength: int = 50,
    n_expr: int = 100,
    n_pcs: int = 10,
) -> tuple[np.ndarray, np.ndarray]:
    """Select/scale/PCA inputs using training samples only."""
    from sklearn.decomposition import PCA

    parts_train: list[np.ndarray] = []
    parts_test: list[np.ndarray] = []

    if feature_type in {"network", "combined"}:
        k = min(n_strength, strength.shape[1])
        node_rank = strength[train_idx].var(axis=0).argsort()[-k:]
        parts_train.append(strength[train_idx][:, node_rank])
        parts_test.append(strength[test_idx][:, node_rank])
        n_pc = min(n_pcs, len(train_idx) - 1, lioness.shape[1])
        if n_pc > 0:
            pca = PCA(n_components=n_pc)
            parts_train.append(pca.fit_transform(lioness[train_idx]))
            parts_test.append(pca.transform(lioness[test_idx]))

    if feature_type in {"expression_baseline", "combined"}:
        k = min(n_expr, expr.shape[1])
        expr_rank = expr[train_idx].var(axis=0).argsort()[-k:]
        parts_train.append(expr[train_idx][:, expr_rank])
        parts_test.append(expr[test_idx][:, expr_rank])

    return np.hstack(parts_train), np.hstack(parts_test)


def fold_safe_source_regression_metrics(
    feature_type: str,
    strength: np.ndarray,
    lioness: np.ndarray,
    expr: np.ndarray,
    y: np.ndarray,
    n_splits: int,
    seed: int,
) -> tuple[dict, int]:
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    pred = np.empty_like(y, dtype=float)
    n_features_seen = []
    for train_idx, test_idx in kf.split(y):
        X_train, X_test = build_fold_features(feature_type, train_idx, test_idx, strength, lioness, expr)
        n_features_seen.append(X_train.shape[1])
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)
        model = Ridge(alpha=1.0)
        model.fit(X_train, y[train_idx])
        pred[test_idx] = model.predict(X_test)
    pear = pearsonr(y, pred)[0] if np.std(pred) > 0 and np.std(y) > 0 else np.nan
    spear = spearmanr(y, pred).correlation if np.std(pred) > 0 and np.std(y) > 0 else np.nan
    return {"pearson": float(pear), "spearman": float(spear)}, int(np.median(n_features_seen))


def main() -> None:
    ap = argparse.ArgumentParser(description="Continuous kidney-stress target validation")
    ap.add_argument("--rtech", default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"))
    ap.add_argument("--meta", default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"))
    ap.add_argument("--phase2_dir", default=str(REPO_ROOT / "data/processed/networks/phase2"))
    ap.add_argument("--id_map", default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"))
    ap.add_argument("--markers", default=",".join(DEFAULT_MARKERS))
    ap.add_argument("--n_splits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/continuous_target"))
    args = ap.parse_args()

    rtech = pd.read_csv(args.rtech, sep="\t", compression="gzip", index_col=0)
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = normalize_labels(meta.set_index(sample_col, drop=False))
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]

    markers = [m.strip() for m in args.markers.split(",") if m.strip()]
    marker_ids = resolve_markers(markers, Path(args.id_map), set(rtech.index))
    y = build_stress_score(rtech, marker_ids)

    phase2 = Path(args.phase2_dir)
    genes = [g.strip() for g in (phase2 / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    edge_i = np.load(phase2 / "edge_i.npy")
    edge_j = np.load(phase2 / "edge_j.npy")
    lioness = np.load(phase2 / "lioness_z_edges.npy")
    samples = [s.strip() for s in (phase2 / "lioness_samples.txt").read_text().splitlines() if s.strip()]
    sample_pos = [samples.index(s) for s in common if s in samples]
    common_lioness_samples = [samples[i] for i in sample_pos]
    lioness = lioness[sample_pos]
    y_lioness = pd.Series(y, index=common).loc[common_lioness_samples].values

    strength = node_strength(lioness, edge_i, edge_j, len(genes))
    marker_gene_idx = {genes.index(g) for g in marker_ids if g in genes}
    eligible_nodes = [i for i in range(len(genes)) if i not in marker_gene_idx]
    strength_source = strength[:, eligible_nodes]

    expr = rtech.loc[[g for g in rtech.index if g not in marker_ids]].T.loc[common_lioness_samples]
    expr_source = expr.values.astype(float)

    rows = []
    for feature_type in ["network", "expression_baseline", "combined"]:
        metrics, n_features = fold_safe_source_regression_metrics(
            feature_type,
            strength_source,
            lioness,
            expr_source,
            y_lioness,
            args.n_splits,
            args.seed,
        )
        rows.append({
            "feature_type": feature_type,
            "n_samples": len(y_lioness),
            "n_features": n_features,
            "target": "kidney_stress_score",
            "markers_requested": ",".join(markers),
            "markers_resolved": ",".join(marker_ids),
            **metrics,
        })

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(outdir / "continuous_target_summary.tsv", sep="\t", index=False)
    print(f"[OK] Continuous target validation written to {outdir}")


if __name__ == "__main__":
    main()
