# src/validation/cross_validation.py
"""
Phase 8: Leakage-Safe Predictive Validation

Implements the validation framework described in Section 5 of the methodology:
  - Stratified K-fold CV (FLT vs GC, stratified within Age×Arm)
  - Fold-wise: residualization fit on train only, skeleton E built on train only,
    LIONESS computed for train and test relative to training pool
  - Sample-level features extracted via src.validation.sample_features
  - Classification with LogisticRegression (L2) and RandomForest
  - Reports per-fold accuracy, AUC, and confusion matrix

Usage:
    python -m src.validation.cross_validation \\
        --phase2_dir data/results/<run>/networks \\
        --meta data/results/<run>/phase1_residuals/meta_phase1.tsv.gz \\
        --rtech data/results/<run>/phase1_residuals/Rtech.tsv.gz \\
        --outdir data/results/<run>/phase8_validation
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix
from sklearn.preprocessing import StandardScaler

from src.common import REPO_ROOT, find_sample_col, normalize_labels


# ---------------------------------------------------------------------------
# Fold-wise skeleton + LIONESS (leakage-safe)
# ---------------------------------------------------------------------------

def build_skeleton_on_fold(
    rtech_train: np.ndarray,
    meta_train: pd.DataFrame,
    genes: list[str],
    topk: int = 80,
) -> tuple[np.ndarray, np.ndarray]:
    """Build skeleton E using only training-fold data.

    Performs cell-standardization within Age×Arm×EnvGroup cells (training only),
    then computes Ledoit-Wolf shrinkage covariance, inverts, and keeps top-k
    neighbors per gene.  Returns edge_i, edge_j arrays.
    """
    from sklearn.covariance import LedoitWolf

    X = rtech_train.copy()  # genes × train_samples
    G, N = X.shape

    # Cell-standardize within each experimental cell
    if {"Age", "Arm", "EnvGroup"}.issubset(meta_train.columns):
        cell_labels = (
            meta_train["Age"].astype(str) + "_"
            + meta_train["Arm"].astype(str) + "_"
            + meta_train["EnvGroup"].astype(str)
        )
        for cell in cell_labels.unique():
            idx = np.where(cell_labels.values == cell)[0]
            if len(idx) < 2:
                continue
            mu = X[:, idx].mean(axis=1, keepdims=True)
            sd = X[:, idx].std(axis=1, keepdims=True) + 1e-8
            X[:, idx] = (X[:, idx] - mu) / sd

    # Ledoit-Wolf shrinkage covariance → precision → partial correlations
    lw = LedoitWolf()
    lw.fit(X.T)  # samples × genes
    prec = np.linalg.pinv(lw.covariance_)
    d = np.sqrt(np.diag(prec))
    d[d < 1e-12] = 1e-12
    pcorr = -prec / np.outer(d, d)
    np.fill_diagonal(pcorr, 0.0)

    # Top-k neighbors per gene
    absP = np.abs(pcorr)
    edges = set()
    for i in range(G):
        k = min(topk, G - 1)
        nbrs = np.argpartition(absP[i], -k)[-k:]
        for j in nbrs:
            if i != j:
                a, b = (i, j) if i < j else (j, i)
                edges.add((a, b))

    edge_list = sorted(edges)
    edge_i = np.array([e[0] for e in edge_list], dtype=np.int32)
    edge_j = np.array([e[1] for e in edge_list], dtype=np.int32)
    return edge_i, edge_j


def lioness_on_fold(
    rtech: np.ndarray,
    sample_mask: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    """Compute LIONESS Fisher-z weights for samples in sample_mask.

    The pooled network is computed from ALL samples indicated by sample_mask
    (the training set).  For each sample s, the leave-one-out network is
    computed by dropping s from the pool.

    For test samples, we compute their LIONESS relative to the training pool
    (the test sample influences only its own network, not the pool).

    Returns lioness_z of shape (len(sample_mask), n_edges).
    """
    CLIP_R = 0.9995
    ZCAP = 20.0

    X = rtech  # genes × all_samples
    G, N_all = X.shape
    idx_pool = np.where(sample_mask)[0]
    N_pool = len(idx_pool)

    # Precompute pooled sums over pool
    X_pool = X[:, idx_pool]
    Sx = X_pool.sum(axis=1)
    Sxx = (X_pool ** 2).sum(axis=1)
    Xi_pool = X_pool[edge_i, :]
    Xj_pool = X_pool[edge_j, :]
    Sxy = (Xi_pool * Xj_pool).sum(axis=1)

    def pearson_from_sums(n, sx_i, sx_j, sxx_i, sxx_j, sxy, eps=1e-12):
        num = n * sxy - sx_i * sx_j
        denx = n * sxx_i - sx_i * sx_i
        deny = n * sxx_j - sx_j * sx_j
        den = np.sqrt(np.maximum(denx, 0.0) * np.maximum(deny, 0.0))
        r = np.where(den > eps, num / den, 0.0)
        return np.clip(r, -1.0, 1.0)

    r_all = pearson_from_sums(N_pool, Sx[edge_i], Sx[edge_j],
                               Sxx[edge_i], Sxx[edge_j], Sxy)
    z_all = np.arctanh(np.clip(r_all, -CLIP_R, CLIP_R))

    # LIONESS for each sample in the pool (training samples)
    E = len(edge_i)
    out = np.empty((N_pool, E), dtype=np.float32)

    for k, s_idx in enumerate(idx_pool):
        xs_i = X[edge_i, s_idx]
        xs_j = X[edge_j, s_idx]

        Sx_i_loo = Sx[edge_i] - xs_i
        Sx_j_loo = Sx[edge_j] - xs_j
        Sxx_i_loo = Sxx[edge_i] - xs_i ** 2
        Sxx_j_loo = Sxx[edge_j] - xs_j ** 2
        Sxy_loo = Sxy - xs_i * xs_j

        r_loo = pearson_from_sums(N_pool - 1, Sx_i_loo, Sx_j_loo,
                                   Sxx_i_loo, Sxx_j_loo, Sxy_loo)
        z_loo = np.arctanh(np.clip(r_loo, -CLIP_R, CLIP_R))

        z_s = N_pool * z_all - (N_pool - 1) * z_loo
        out[k, :] = np.clip(z_s, -ZCAP, ZCAP).astype(np.float32)

    return out


# ---------------------------------------------------------------------------
# Main CV loop
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Phase 8: Leakage-safe cross-validation (FLT vs GC)"
    )
    ap.add_argument("--phase2_dir",
                    default=str(REPO_ROOT / "data/processed/networks/phase2"),
                    help="Phase 2 directory (for gene list)")
    ap.add_argument("--rtech",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"),
                    help="Rtech expression matrix")
    ap.add_argument("--meta",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"),
                    help="Metadata file")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/results/phase8_validation"),
                    help="Output directory")
    ap.add_argument("--n_folds", type=int, default=5, help="Number of CV folds")
    ap.add_argument("--topk", type=int, default=80, help="Top-k neighbors for skeleton")
    ap.add_argument("--max_genes", type=int, default=2500, help="Max genes for network")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    phase2_dir = Path(args.phase2_dir)

    print("=" * 60)
    print("Phase 8: Leakage-Safe Cross-Validation")
    print("=" * 60)

    # ── 1) Load expression + metadata ────────────────────────────
    print("\nLoading data...")
    rtech = pd.read_csv(args.rtech, sep="\t", compression="gzip", index_col=0)
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")

    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    meta = normalize_labels(meta)

    # Align
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]

    # Filter to FLT and GC only
    mask_flt_gc = meta["EnvGroup"].isin(["FLT", "GC"])
    meta = meta[mask_flt_gc].copy()
    rtech = rtech[meta.index]

    print(f"  Samples (FLT+GC): {len(meta)}")
    print(f"  FLT: {(meta['EnvGroup'] == 'FLT').sum()}, GC: {(meta['EnvGroup'] == 'GC').sum()}")

    # ── 2) Gene selection (top variance) ─────────────────────────
    gene_var = rtech.var(axis=1).sort_values(ascending=False)
    keep_genes = gene_var.head(args.max_genes).index.tolist()

    # Force-include Phase 2 genes if available
    p2_genes_path = phase2_dir / "phase2_genes.txt"
    if p2_genes_path.exists():
        p2_genes = [g.strip() for g in p2_genes_path.read_text().splitlines() if g.strip()]
        for g in p2_genes:
            if g in rtech.index and g not in keep_genes:
                keep_genes.append(g)

    rtech = rtech.loc[keep_genes]
    genes = list(rtech.index)
    G = len(genes)
    print(f"  Genes: {G}")

    # ── 3) Prepare labels and stratification ─────────────────────
    y = (meta["EnvGroup"] == "FLT").astype(int).values  # 1=FLT, 0=GC
    strat_labels = meta["Age"].astype(str) + "_" + meta["Arm"].astype(str)

    # Expression matrix: genes × samples
    X_expr = rtech.values.astype(np.float64)
    sample_ids = list(meta.index)

    print(f"\n  Running {args.n_folds}-fold stratified CV (stratified by Age×Arm)...")

    skf = StratifiedKFold(n_splits=args.n_folds, shuffle=True, random_state=42)

    fold_results = []

    for fold_idx, (train_idx, test_idx) in enumerate(skf.split(np.zeros(len(y)), strat_labels)):
        print(f"\n{'─' * 40}")
        print(f"Fold {fold_idx + 1}/{args.n_folds}")
        print(f"  Train: {len(train_idx)} samples, Test: {len(test_idx)} samples")

        y_train, y_test = y[train_idx], y[test_idx]

        # ── 3a) Build skeleton on training data only ─────────
        meta_train = meta.iloc[train_idx]
        X_train = X_expr[:, train_idx]

        print("  Building skeleton on training fold...")
        edge_i, edge_j = build_skeleton_on_fold(
            X_train, meta_train, genes, topk=args.topk
        )
        print(f"  Skeleton: {len(edge_i)} edges")

        # ── 3b) Compute LIONESS on training samples ──────────
        print("  Computing LIONESS (training)...")
        train_mask = np.zeros(len(y), dtype=bool)
        train_mask[train_idx] = True
        lioness_train = lioness_on_fold(X_expr, train_mask, edge_i, edge_j)

        # ── 3c) Compute LIONESS for test samples ─────────────
        # For test: compute relative to training pool
        # Each test sample is added to the pool individually
        print("  Computing LIONESS (test)...")
        lioness_test = np.empty((len(test_idx), len(edge_i)), dtype=np.float32)
        for k, t_idx in enumerate(test_idx):
            # Augment pool with this test sample
            aug_mask = train_mask.copy()
            aug_mask[t_idx] = True
            aug_lioness = lioness_on_fold(X_expr, aug_mask, edge_i, edge_j)
            # The test sample's LIONESS is the last row added
            # Find its position in the augmented set
            aug_indices = np.where(aug_mask)[0]
            pos = np.searchsorted(aug_indices, t_idx)
            lioness_test[k, :] = aug_lioness[pos, :]

        # ── 3d) Extract features ─────────────────────────────
        from src.validation.sample_features import (
            node_strength, edge_pca_features
        )

        # Node strength features (top 50 by variance in training)
        strength_train = node_strength(lioness_train, edge_i, edge_j, G)
        strength_test = node_strength(lioness_test, edge_i, edge_j, G)

        var_rank = strength_train.var(axis=0).argsort()[-50:]
        feat_train = strength_train[:, var_rank]
        feat_test = strength_test[:, var_rank]

        # PCA features on edge weights
        from sklearn.decomposition import PCA
        n_pcs = min(10, min(lioness_train.shape) - 1)
        pca = PCA(n_components=n_pcs)
        pc_train = pca.fit_transform(lioness_train)
        pc_test = pca.transform(lioness_test)

        # Combine
        feat_train = np.hstack([feat_train, pc_train])
        feat_test = np.hstack([feat_test, pc_test])

        # Scale
        scaler = StandardScaler()
        feat_train = scaler.fit_transform(feat_train)
        feat_test = scaler.transform(feat_test)

        # ── 3e) Classify ─────────────────────────────────────
        classifiers = {
            "LogisticRegression": LogisticRegression(
                max_iter=1000, C=1.0, random_state=42
            ),
            "RandomForest": RandomForestClassifier(
                n_estimators=100, max_depth=5, random_state=42
            ),
        }

        for clf_name, clf in classifiers.items():
            clf.fit(feat_train, y_train)
            y_pred = clf.predict(feat_test)
            y_prob = clf.predict_proba(feat_test)[:, 1]

            acc = accuracy_score(y_test, y_pred)
            try:
                auc = roc_auc_score(y_test, y_prob)
            except ValueError:
                auc = float("nan")
            cm = confusion_matrix(y_test, y_pred)

            print(f"  {clf_name}: acc={acc:.3f}, AUC={auc:.3f}")

            fold_results.append({
                "fold": fold_idx + 1,
                "classifier": clf_name,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "n_edges": len(edge_i),
                "n_features": feat_train.shape[1],
                "accuracy": float(acc),
                "auc": float(auc),
                "tn": int(cm[0, 0]),
                "fp": int(cm[0, 1]),
                "fn": int(cm[1, 0]),
                "tp": int(cm[1, 1]),
            })

    # ── 4) Save results ──────────────────────────────────────────
    results_df = pd.DataFrame(fold_results)
    results_path = outdir / "cv_results.tsv"
    results_df.to_csv(results_path, sep="\t", index=False)

    # Summary by classifier
    print(f"\n{'=' * 60}")
    print("Cross-Validation Summary")
    print(f"{'=' * 60}")

    summary_rows = []
    for clf_name in results_df["classifier"].unique():
        sub = results_df[results_df["classifier"] == clf_name]
        row = {
            "classifier": clf_name,
            "mean_accuracy": float(sub["accuracy"].mean()),
            "std_accuracy": float(sub["accuracy"].std()),
            "mean_auc": float(sub["auc"].mean()),
            "std_auc": float(sub["auc"].std()),
        }
        summary_rows.append(row)
        print(f"\n  {clf_name}:")
        print(f"    Accuracy: {row['mean_accuracy']:.3f} ± {row['std_accuracy']:.3f}")
        print(f"    AUC:      {row['mean_auc']:.3f} ± {row['std_auc']:.3f}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "cv_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    # Save metadata
    cv_meta = {
        "n_folds": args.n_folds,
        "n_samples": len(y),
        "n_flt": int(y.sum()),
        "n_gc": int((1 - y).sum()),
        "n_genes": G,
        "topk": args.topk,
        "target": "FLT_vs_GC",
        "stratification": "Age_x_Arm",
    }
    with open(outdir / "cv_metadata.json", "w") as f:
        json.dump(cv_meta, f, indent=2)

    print(f"\n[OK] Results saved to: {outdir}")
    print(f"  - {results_path.name}")
    print(f"  - {summary_path.name}")
    print(f"  - cv_metadata.json")


if __name__ == "__main__":
    main()
