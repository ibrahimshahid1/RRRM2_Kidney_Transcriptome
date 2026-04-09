# src/validation/enhanced_cv.py
"""
Phase 8b: Enhanced Predictive Validation

Three key improvements over the original Phase 8:

  1. STRATUM-SPECIFIC LOO-CV: Run classification within each Age×Arm stratum
     (n=10: 5 FLT + 5 GC) using leave-one-out CV.  This avoids asking the
     classifier to learn one decision boundary across biologically divergent
     subgroups.

  2. EXPRESSION BASELINE: Run the classifier on raw expression features
     (top-variance genes) as a control.  If network features match or exceed
     expression features, the network representation contains non-redundant
     biological information.

  3. PERMUTATION BASELINE: Shuffle FLT/GC labels 1,000 times and re-run the
     full classification pipeline each time.  Report the permutation p-value
     for each stratum's accuracy, not just the raw accuracy.

Feature types compared:
  A) Network topology: node strength (top-50 most variable) + PCA on edge weights
  B) Raw expression: top-50 most variable genes' expression levels
  C) Combined: A + B

Outputs:
  enhanced_cv_results.tsv      – per-fold, per-stratum, per-feature-type results
  enhanced_cv_summary.tsv      – summary statistics
  enhanced_cv_permutation.tsv  – permutation null distributions
  enhanced_cv_metadata.json    – run configuration

Usage:
    python -m src.validation.enhanced_cv \\
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
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from src.common import REPO_ROOT, find_sample_col, normalize_labels


# ---------------------------------------------------------------------------
# Lightweight LIONESS for small strata (reuses skeleton from Phase 2)
# ---------------------------------------------------------------------------

def compute_lioness_for_stratum(
    rtech_sub: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    """Compute LIONESS Fisher-z weights for a small stratum.

    Args:
        rtech_sub: genes × samples expression matrix (already residualized)
        edge_i, edge_j: shared skeleton edges
    Returns:
        lioness_z: samples × edges
    """
    CLIP_R = 0.9995
    ZCAP = 20.0
    G, N = rtech_sub.shape
    E = len(edge_i)

    # Pooled sums
    Sx = rtech_sub.sum(axis=1)
    Sxx = (rtech_sub ** 2).sum(axis=1)
    Xi = rtech_sub[edge_i, :]
    Xj = rtech_sub[edge_j, :]
    Sxy = (Xi * Xj).sum(axis=1)

    def pearson_from_sums(n, sx_i, sx_j, sxx_i, sxx_j, sxy_val, eps=1e-12):
        num = n * sxy_val - sx_i * sx_j
        denx = n * sxx_i - sx_i * sx_i
        deny = n * sxx_j - sx_j * sx_j
        den = np.sqrt(np.maximum(denx, 0.0) * np.maximum(deny, 0.0))
        r = np.where(den > eps, num / den, 0.0)
        return np.clip(r, -1.0, 1.0)

    r_all = pearson_from_sums(N, Sx[edge_i], Sx[edge_j],
                               Sxx[edge_i], Sxx[edge_j], Sxy)
    z_all = np.arctanh(np.clip(r_all, -CLIP_R, CLIP_R))

    out = np.empty((N, E), dtype=np.float32)
    for s in range(N):
        xs_i = rtech_sub[edge_i, s]
        xs_j = rtech_sub[edge_j, s]

        Sx_i_loo = Sx[edge_i] - xs_i
        Sx_j_loo = Sx[edge_j] - xs_j
        Sxx_i_loo = Sxx[edge_i] - xs_i ** 2
        Sxx_j_loo = Sxx[edge_j] - xs_j ** 2
        Sxy_loo = Sxy - xs_i * xs_j

        r_loo = pearson_from_sums(N - 1, Sx_i_loo, Sx_j_loo,
                                   Sxx_i_loo, Sxx_j_loo, Sxy_loo)
        z_loo = np.arctanh(np.clip(r_loo, -CLIP_R, CLIP_R))

        z_s = N * z_all - (N - 1) * z_loo
        out[s, :] = np.clip(z_s, -ZCAP, ZCAP).astype(np.float32)

    return out


def extract_network_features(
    lioness_z: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    n_genes: int,
    n_strength: int = 50,
    n_pcs: int = 10,
) -> np.ndarray:
    """Extract network topology features from LIONESS weights.

    Features:
      - Top-k node strengths (sum of incident edge weights)
      - PCA on edge weights
    """
    n_samples = lioness_z.shape[0]

    # Node strength
    strength = np.zeros((n_samples, n_genes), dtype=np.float32)
    for s in range(n_samples):
        w = lioness_z[s]
        np.add.at(strength[s], edge_i, w)
        np.add.at(strength[s], edge_j, w)

    # Top-k most variable nodes
    k = min(n_strength, n_genes)
    var_rank = strength.var(axis=0).argsort()[-k:]
    feat_strength = strength[:, var_rank]

    # PCA on edge weights
    n_pc = min(n_pcs, min(lioness_z.shape) - 1, 10)
    if n_pc > 0:
        pca = PCA(n_components=n_pc)
        pc = pca.fit_transform(lioness_z)
        return np.hstack([feat_strength, pc])
    return feat_strength


def run_classification(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    y_test: np.ndarray,
    clf_name: str = "RandomForest",
) -> dict:
    """Run a single train/test classification and return metrics."""
    # Scale
    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)

    if clf_name == "RandomForest":
        clf = RandomForestClassifier(
            n_estimators=200, max_depth=3, min_samples_leaf=2,
            random_state=42
        )
    else:
        clf = LogisticRegression(
            max_iter=2000, C=0.1, penalty='l2', random_state=42, solver='lbfgs'
        )

    clf.fit(X_tr, y_train)
    y_pred = clf.predict(X_te)

    try:
        y_prob = clf.predict_proba(X_te)[:, 1]
    except Exception:
        y_prob = np.full(len(y_test), 0.5)

    acc = accuracy_score(y_test, y_pred)
    try:
        auc = roc_auc_score(y_test, y_prob)
    except ValueError:
        auc = float("nan")

    return {
        "accuracy": float(acc),
        "auc": float(auc),
        "y_pred": int(y_pred[0]) if len(y_pred) == 1 else y_pred.tolist(),
        "y_true": int(y_test[0]) if len(y_test) == 1 else y_test.tolist(),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Phase 8b: Enhanced Predictive Validation"
    )
    ap.add_argument("--phase2_dir",
                    default=str(REPO_ROOT / "data/processed/networks/phase2"),
                    help="Phase 2 directory (for gene list and edges)")
    ap.add_argument("--rtech",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"),
                    help="Rtech expression matrix")
    ap.add_argument("--meta",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"),
                    help="Metadata file")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/results/phase8_validation"),
                    help="Output directory")
    ap.add_argument("--topk", type=int, default=80, help="Top-k neighbors for skeleton")
    ap.add_argument("--max_genes", type=int, default=2500, help="Max genes for network")
    ap.add_argument("--n_perms", type=int, default=1000,
                    help="Number of permutations for null distribution")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    phase2_dir = Path(args.phase2_dir)

    print("=" * 70)
    print("Phase 8b: Enhanced Predictive Validation")
    print("  - Stratum-specific LOO-CV")
    print("  - Expression baseline comparison")
    print("  - Permutation null distribution")
    print("=" * 70)

    # ── 1) Load data ─────────────────────────────────────────────────
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
    print(f"  FLT: {(meta['EnvGroup'] == 'FLT').sum()}, "
          f"GC: {(meta['EnvGroup'] == 'GC').sum()}")

    # ── 2) Gene selection ────────────────────────────────────────────
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
    X_expr_full = rtech.values.astype(np.float64)  # genes × samples
    print(f"  Genes: {G}")

    # ── 3) Load Phase 2 skeleton (if available) ──────────────────────
    edge_i_path = phase2_dir / "edge_i.npy"
    edge_j_path = phase2_dir / "edge_j.npy"
    has_skeleton = edge_i_path.exists() and edge_j_path.exists()

    if has_skeleton:
        edge_i = np.load(edge_i_path)
        edge_j = np.load(edge_j_path)
        print(f"  Loaded skeleton: {len(edge_i)} edges")
    else:
        print("  WARNING: No Phase 2 skeleton found. Building from scratch.")
        # Fall back to building from the full dataset
        from src.validation.cross_validation import build_skeleton_on_fold
        edge_i, edge_j = build_skeleton_on_fold(
            X_expr_full, meta, genes, topk=args.topk
        )
        print(f"  Built skeleton: {len(edge_i)} edges")

    # ── 4) Define strata ─────────────────────────────────────────────
    strata = {
        "ISS-T_Young": {"Age": "YNG", "Arm": "ISS-T"},
        "ISS-T_Old": {"Age": "OLD", "Arm": "ISS-T"},
        "LAR_Young": {"Age": "YNG", "Arm": "LAR"},
        "LAR_Old": {"Age": "OLD", "Arm": "LAR"},
    }

    y_all = (meta["EnvGroup"] == "FLT").astype(int).values
    all_results = []
    perm_results = []

    # ── 5) Run stratum-specific LOO-CV ───────────────────────────────
    for stratum_name, filt in strata.items():
        mask = (meta["Age"] == filt["Age"]) & (meta["Arm"] == filt["Arm"])
        sub_idx = np.where(mask.values)[0]

        if len(sub_idx) < 4:
            print(f"\n  [SKIP] {stratum_name}: only {len(sub_idx)} samples")
            continue

        sub_meta = meta.iloc[sub_idx]
        sub_y = y_all[sub_idx]
        sub_expr = X_expr_full[:, sub_idx]  # genes × stratum_samples
        n_sub = len(sub_idx)
        n_flt = sub_y.sum()
        n_gc = n_sub - n_flt

        print(f"\n{'━' * 60}")
        print(f"Stratum: {stratum_name} (n={n_sub}: {n_flt} FLT, {n_gc} GC)")
        print(f"{'━' * 60}")

        # ── 5a) Compute LIONESS for this stratum ─────────────────
        print("  Computing LIONESS for stratum...")
        lioness_z = compute_lioness_for_stratum(sub_expr, edge_i, edge_j)
        print(f"  LIONESS shape: {lioness_z.shape}")

        # ── 5b) Extract features ─────────────────────────────────
        # Network topology features
        net_features = extract_network_features(
            lioness_z, edge_i, edge_j, G,
            n_strength=50, n_pcs=min(10, n_sub - 2)
        )

        # Expression features (top 50 most variable within stratum)
        expr_var = sub_expr.var(axis=1)
        top_expr_idx = expr_var.argsort()[-50:]
        expr_features = sub_expr[top_expr_idx, :].T  # samples × 50

        # Combined features
        combined_features = np.hstack([net_features, expr_features])

        feature_sets = {
            "network_topology": net_features,
            "expression_baseline": expr_features,
            "combined": combined_features,
        }

        # ── 5c) LOO-CV for each feature set and classifier ───────
        for feat_name, feat_matrix in feature_sets.items():
            for clf_name in ["RandomForest", "LogisticRegression"]:
                loo = LeaveOneOut()
                fold_preds = []
                fold_trues = []

                for train_idx_loo, test_idx_loo in loo.split(feat_matrix):
                    X_train = feat_matrix[train_idx_loo]
                    y_train = sub_y[train_idx_loo]
                    X_test = feat_matrix[test_idx_loo]
                    y_test = sub_y[test_idx_loo]

                    result = run_classification(
                        X_train, y_train, X_test, y_test, clf_name
                    )
                    fold_preds.append(result["y_pred"])
                    fold_trues.append(result["y_true"])

                # Aggregate LOO results
                all_preds = np.array(fold_preds)
                all_trues = np.array(fold_trues)
                loo_acc = accuracy_score(all_trues, all_preds)

                try:
                    loo_auc = roc_auc_score(all_trues, all_preds)
                except ValueError:
                    loo_auc = float("nan")

                print(f"  {feat_name:25s} | {clf_name:20s} | "
                      f"LOO acc={loo_acc:.3f} | AUC={loo_auc:.3f}")

                all_results.append({
                    "stratum": stratum_name,
                    "feature_type": feat_name,
                    "classifier": clf_name,
                    "n_samples": n_sub,
                    "n_flt": int(n_flt),
                    "n_gc": int(n_gc),
                    "n_features": feat_matrix.shape[1],
                    "loo_accuracy": float(loo_acc),
                    "loo_auc": float(loo_auc),
                    "n_correct": int((all_preds == all_trues).sum()),
                })

        # ── 5d) Permutation null (network features + RF only) ────
        print(f"\n  Running {args.n_perms} permutations (network + RF)...")
        rng = np.random.default_rng(42)
        perm_accs = []

        for p_idx in range(args.n_perms):
            # Shuffle labels
            perm_y = sub_y.copy()
            rng.shuffle(perm_y)

            # LOO-CV with shuffled labels
            loo = LeaveOneOut()
            p_preds = []
            p_trues = []

            for train_idx_loo, test_idx_loo in loo.split(net_features):
                X_train = net_features[train_idx_loo]
                y_train = perm_y[train_idx_loo]
                X_test = net_features[test_idx_loo]
                y_test = perm_y[test_idx_loo]

                res = run_classification(X_train, y_train, X_test, y_test, "RandomForest")
                p_preds.append(res["y_pred"])
                p_trues.append(res["y_true"])

            perm_acc = accuracy_score(np.array(p_trues), np.array(p_preds))
            perm_accs.append(perm_acc)

            if (p_idx + 1) % 200 == 0:
                print(f"    Permutation {p_idx + 1}/{args.n_perms}...")

        perm_accs = np.array(perm_accs)

        # Get observed accuracy for network + RF
        obs_acc = None
        for r in all_results:
            if (r["stratum"] == stratum_name and
                r["feature_type"] == "network_topology" and
                r["classifier"] == "RandomForest"):
                obs_acc = r["loo_accuracy"]
                break

        if obs_acc is not None:
            perm_p = (1 + (perm_accs >= obs_acc).sum()) / (args.n_perms + 1)
        else:
            perm_p = 1.0

        perm_results.append({
            "stratum": stratum_name,
            "observed_accuracy": float(obs_acc) if obs_acc else float("nan"),
            "perm_mean_accuracy": float(perm_accs.mean()),
            "perm_std_accuracy": float(perm_accs.std()),
            "perm_p_value": float(perm_p),
            "perm_95th": float(np.quantile(perm_accs, 0.95)),
            "n_perms": args.n_perms,
        })

        print(f"  Permutation null: mean={perm_accs.mean():.3f} ± {perm_accs.std():.3f}")
        print(f"  Observed accuracy: {obs_acc:.3f}")
        print(f"  Permutation p-value: {perm_p:.4f}")

    # ── 6) Also run pooled analysis (all strata combined) ────────────
    print(f"\n{'━' * 60}")
    print(f"Pooled Analysis (all strata, 5-fold stratified CV)")
    print(f"{'━' * 60}")

    strat_labels = meta["Age"].astype(str) + "_" + meta["Arm"].astype(str)
    n_total = len(y_all)

    # Compute LIONESS on full FLT+GC set
    print("  Computing LIONESS for full dataset...")
    lioness_full = compute_lioness_for_stratum(X_expr_full, edge_i, edge_j)

    net_features_full = extract_network_features(
        lioness_full, edge_i, edge_j, G,
        n_strength=50, n_pcs=min(20, n_total - 2)
    )

    # Expression baseline
    expr_var_full = X_expr_full.var(axis=1)
    top_expr_idx_full = expr_var_full.argsort()[-50:]
    expr_features_full = X_expr_full[top_expr_idx_full, :].T

    combined_full = np.hstack([net_features_full, expr_features_full])

    feature_sets_full = {
        "network_topology": net_features_full,
        "expression_baseline": expr_features_full,
        "combined": combined_full,
    }

    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    for feat_name, feat_matrix in feature_sets_full.items():
        for clf_name in ["RandomForest", "LogisticRegression"]:
            fold_accs = []
            fold_aucs = []

            for train_idx, test_idx in skf.split(feat_matrix, strat_labels):
                result = run_classification(
                    feat_matrix[train_idx], y_all[train_idx],
                    feat_matrix[test_idx], y_all[test_idx],
                    clf_name
                )
                fold_accs.append(result["accuracy"])
                fold_aucs.append(result["auc"])

            mean_acc = np.mean(fold_accs)
            mean_auc = np.nanmean(fold_aucs)

            print(f"  {feat_name:25s} | {clf_name:20s} | "
                  f"acc={mean_acc:.3f}±{np.std(fold_accs):.3f} | "
                  f"AUC={mean_auc:.3f}±{np.nanstd(fold_aucs):.3f}")

            all_results.append({
                "stratum": "POOLED",
                "feature_type": feat_name,
                "classifier": clf_name,
                "n_samples": n_total,
                "n_flt": int(y_all.sum()),
                "n_gc": int((1 - y_all).sum()),
                "n_features": feat_matrix.shape[1],
                "loo_accuracy": float(mean_acc),
                "loo_auc": float(mean_auc),
                "n_correct": -1,  # Not applicable for 5-fold
            })

    # ── 7) Save results ──────────────────────────────────────────────
    results_df = pd.DataFrame(all_results)
    results_path = outdir / "enhanced_cv_results.tsv"
    results_df.to_csv(results_path, sep="\t", index=False)

    perm_df = pd.DataFrame(perm_results)
    perm_path = outdir / "enhanced_cv_permutation.tsv"
    perm_df.to_csv(perm_path, sep="\t", index=False)

    # Summary: best configuration per stratum
    print(f"\n{'=' * 70}")
    print("ENHANCED VALIDATION SUMMARY")
    print(f"{'=' * 70}")

    summary_rows = []
    for stratum in results_df["stratum"].unique():
        sub = results_df[results_df["stratum"] == stratum]
        best = sub.loc[sub["loo_accuracy"].idxmax()]

        # Find permutation p-value if available
        perm_row = perm_df[perm_df["stratum"] == stratum] if stratum != "POOLED" else pd.DataFrame()
        perm_p = float(perm_row["perm_p_value"].iloc[0]) if len(perm_row) > 0 else float("nan")

        # Compare network vs expression
        net_acc = sub[(sub["feature_type"] == "network_topology") &
                      (sub["classifier"] == "RandomForest")]["loo_accuracy"].values
        expr_acc = sub[(sub["feature_type"] == "expression_baseline") &
                       (sub["classifier"] == "RandomForest")]["loo_accuracy"].values

        net_acc = float(net_acc[0]) if len(net_acc) > 0 else float("nan")
        expr_acc = float(expr_acc[0]) if len(expr_acc) > 0 else float("nan")

        row = {
            "stratum": stratum,
            "best_feature_type": best["feature_type"],
            "best_classifier": best["classifier"],
            "best_accuracy": float(best["loo_accuracy"]),
            "network_rf_accuracy": net_acc,
            "expression_rf_accuracy": expr_acc,
            "network_advantage": net_acc - expr_acc if not (np.isnan(net_acc) or np.isnan(expr_acc)) else float("nan"),
            "permutation_p_value": perm_p,
        }
        summary_rows.append(row)

        print(f"\n  {stratum}:")
        print(f"    Best: {best['feature_type']} + {best['classifier']} = {best['loo_accuracy']:.3f}")
        print(f"    Network (RF): {net_acc:.3f}  |  Expression (RF): {expr_acc:.3f}")
        if not np.isnan(perm_p):
            sig = "***" if perm_p < 0.001 else "**" if perm_p < 0.01 else "*" if perm_p < 0.05 else "n.s."
            print(f"    Permutation p-value: {perm_p:.4f} {sig}")

    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "enhanced_cv_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    # Save metadata
    cv_meta = {
        "analysis": "Enhanced Predictive Validation (Phase 8b)",
        "improvements": [
            "Stratum-specific LOO-CV (n=10 per stratum)",
            "Expression baseline comparison",
            "Permutation null distribution (n_perms={})".format(args.n_perms),
        ],
        "n_samples_total": int(len(y_all)),
        "n_flt": int(y_all.sum()),
        "n_gc": int((1 - y_all).sum()),
        "n_genes": G,
        "topk": args.topk,
        "strata": list(strata.keys()),
        "feature_types": ["network_topology", "expression_baseline", "combined"],
        "classifiers": ["RandomForest", "LogisticRegression"],
        "n_permutations": args.n_perms,
    }
    with open(outdir / "enhanced_cv_metadata.json", "w") as f:
        json.dump(cv_meta, f, indent=2)

    print(f"\n[OK] Enhanced validation results saved to: {outdir}")
    print(f"  - {results_path.name}")
    print(f"  - {perm_path.name}")
    print(f"  - {summary_path.name}")
    print(f"  - enhanced_cv_metadata.json")


if __name__ == "__main__":
    main()
