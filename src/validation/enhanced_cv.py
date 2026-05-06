# src/validation/enhanced_cv.py
"""
Phase 8b: leakage-safe enhanced predictive validation.

For every fold, this module computes the network pool, skeleton, LIONESS or
alternative sample-specific weights, feature selection, scaling, PCA, and model
fit inside the fold. It preserves expression-only baselines and evaluates
multiple LIONESS pooling modes and feature sets.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import LeaveOneOut, StratifiedKFold
from sklearn.preprocessing import StandardScaler

from src.common import REPO_ROOT, find_sample_col, id_map_lookup, normalize_labels
from src.networks.alternative_methods import compute_alternative_network
from src.validation.cross_validation import (
    DEFAULT_RF_MAX_DEPTH,
    DEFAULT_TOPK,
    build_skeleton_on_fold,
    lioness_on_fold,
)
from src.validation.sample_features import node_strength

POOL_MODES = ["age_arm_envgroup", "arm", "age", "global"]
NETWORK_FEATURE_SETS = [
    "node_strength",
    "pathway_strength",
    "sparse_edges",
    "edge_pca",
    "network_combined",
]


def load_curated_pathway_masks(gene_set_yaml: Path, id_map: Path, genes: list[str]) -> dict[str, np.ndarray]:
    try:
        import yaml
    except Exception:
        return {}
    if not gene_set_yaml.exists() or not id_map.exists():
        return {}
    _, symbol_to_ens = id_map_lookup(id_map)
    cfg = yaml.safe_load(gene_set_yaml.read_text()) or {}
    gene_universe = set(genes)
    masks: dict[str, np.ndarray] = {}
    for name, val in cfg.items():
        if not isinstance(val, dict) or "genes" not in val:
            continue
        symbols = []
        for item in val["genes"]:
            if isinstance(item, str):
                symbols.append(item)
            elif isinstance(item, list):
                symbols.extend(item)
            elif isinstance(item, dict):
                for sublist in item.values():
                    if isinstance(sublist, list):
                        symbols.extend(sublist)
        matched = set()
        for symbol in symbols:
            matched |= (symbol_to_ens.get(str(symbol).lower(), set()) & gene_universe)
        if len(matched) >= 3:
            masks[name] = np.array([g in matched for g in genes], dtype=bool)
    return masks


def run_classifier(X_train: np.ndarray, y_train: np.ndarray, X_test: np.ndarray, classifier: str) -> tuple[int, float]:
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)
    if classifier == "RandomForest":
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=DEFAULT_RF_MAX_DEPTH,
            min_samples_split=5,
            min_samples_leaf=2,
            random_state=42,
        )
    else:
        model = LogisticRegression(max_iter=2000, C=0.1, penalty="l2", random_state=42)
    model.fit(X_train, y_train)
    pred = int(model.predict(X_test)[0])
    try:
        prob = float(model.predict_proba(X_test)[0, 1])
    except Exception:
        prob = 0.5
    return pred, prob


def pool_mask_for_mode(meta: pd.DataFrame, eval_indices: np.ndarray, train_indices: np.ndarray, mode: str) -> np.ndarray:
    """Return training-only samples allowed to define the network pool."""
    train_mask = np.zeros(len(meta), dtype=bool)
    train_mask[train_indices] = True
    if mode == "global":
        return train_mask

    eval_meta = meta.iloc[eval_indices]
    if mode == "age_arm_envgroup":
        keys = set(zip(eval_meta["Age"].astype(str), eval_meta["Arm"].astype(str)))
        return train_mask & np.array([
            (str(row.Age), str(row.Arm)) in keys for row in meta.itertuples()
        ])
    if mode == "arm":
        arms = set(eval_meta["Arm"].astype(str))
        return train_mask & meta["Arm"].astype(str).isin(arms).values
    if mode == "age":
        ages = set(eval_meta["Age"].astype(str))
        return train_mask & meta["Age"].astype(str).isin(ages).values
    raise ValueError(f"Unknown pooling mode: {mode}")


def sample_specific_weights(
    method: str,
    X_expr: np.ndarray,
    pool_mask: np.ndarray,
    sample_indices: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    """Compute sample-specific edge weights relative to a training-only pool."""
    pool_indices = np.where(pool_mask)[0]
    if len(pool_indices) < 4:
        raise ValueError("Network pool has fewer than 4 training samples")

    rows = []
    in_pool = {idx: pos for pos, idx in enumerate(pool_indices)}
    pool_weights = None
    for sample_idx in sample_indices:
        if sample_idx in in_pool:
            if pool_weights is None:
                if method == "lioness":
                    pool_weights = lioness_on_fold(X_expr, pool_mask, edge_i, edge_j)
                else:
                    pool_weights = compute_alternative_network(
                        method, X_expr[:, pool_indices], edge_i, edge_j
                    )
            rows.append(pool_weights[in_pool[sample_idx]])
        else:
            aug_mask = pool_mask.copy()
            aug_mask[sample_idx] = True
            aug_indices = np.where(aug_mask)[0]
            pos = int(np.where(aug_indices == sample_idx)[0][0])
            if method == "lioness":
                aug_weights = lioness_on_fold(X_expr, aug_mask, edge_i, edge_j)
            else:
                aug_weights = compute_alternative_network(method, X_expr[:, aug_indices], edge_i, edge_j)
            rows.append(aug_weights[pos])
    return np.vstack(rows).astype(np.float32)


def pathway_strength_features(
    weights: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    pathway_masks: dict[str, np.ndarray],
) -> np.ndarray:
    feats = []
    for mask in pathway_masks.values():
        edge_mask = mask[edge_i] | mask[edge_j]
        if edge_mask.any():
            feats.append(np.abs(weights[:, edge_mask]).mean(axis=1))
    if not feats:
        return np.empty((weights.shape[0], 0), dtype=float)
    return np.vstack(feats).T


def build_network_features_fold(
    feature_set: str,
    train_w: np.ndarray,
    test_w: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    n_genes: int,
    pathway_masks: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    parts_train: list[np.ndarray] = []
    parts_test: list[np.ndarray] = []

    if feature_set in {"node_strength", "network_combined"}:
        strength_train = node_strength(train_w, edge_i, edge_j, n_genes)
        strength_test = node_strength(test_w, edge_i, edge_j, n_genes)
        k = min(50, n_genes)
        rank = strength_train.var(axis=0).argsort()[-k:]
        parts_train.append(strength_train[:, rank])
        parts_test.append(strength_test[:, rank])

    if feature_set in {"pathway_strength", "network_combined"}:
        p_train = pathway_strength_features(train_w, edge_i, edge_j, pathway_masks)
        p_test = pathway_strength_features(test_w, edge_i, edge_j, pathway_masks)
        if p_train.shape[1] > 0:
            parts_train.append(p_train)
            parts_test.append(p_test)

    if feature_set in {"sparse_edges", "network_combined"}:
        k = min(100, train_w.shape[1])
        rank = train_w.var(axis=0).argsort()[-k:]
        parts_train.append(train_w[:, rank])
        parts_test.append(test_w[:, rank])

    if feature_set in {"edge_pca", "network_combined"}:
        n_pc = min(10, min(train_w.shape) - 1)
        if n_pc > 0:
            pca = PCA(n_components=n_pc)
            parts_train.append(pca.fit_transform(train_w))
            parts_test.append(pca.transform(test_w))

    if not parts_train:
        raise ValueError(f"Feature set {feature_set} produced no features")
    return np.hstack(parts_train), np.hstack(parts_test)


def build_expression_features_fold(
    X_expr: np.ndarray,
    train_indices: np.ndarray,
    test_indices: np.ndarray,
    max_features: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    train_expr = X_expr[:, train_indices].T
    test_expr = X_expr[:, test_indices].T
    k = min(max_features, train_expr.shape[1])
    rank = train_expr.var(axis=0).argsort()[-k:]
    return train_expr[:, rank], test_expr[:, rank]


def evaluation_groups(meta: pd.DataFrame) -> dict[str, np.ndarray]:
    groups = {
        "ISS-T_Young": np.where((meta["Age"] == "YNG") & (meta["Arm"] == "ISS-T"))[0],
        "ISS-T_Old": np.where((meta["Age"] == "OLD") & (meta["Arm"] == "ISS-T"))[0],
        "LAR_Young": np.where((meta["Age"] == "YNG") & (meta["Arm"] == "LAR"))[0],
        "LAR_Old": np.where((meta["Age"] == "OLD") & (meta["Arm"] == "LAR"))[0],
        "POOLED": np.arange(len(meta)),
    }
    return {k: v for k, v in groups.items() if len(v) >= 4}


def group_folds(group_indices: np.ndarray, y_all: np.ndarray, meta: pd.DataFrame, n_splits: int):
    if len(group_indices) <= 12:
        loo = LeaveOneOut()
        for local_train, local_test in loo.split(group_indices):
            yield group_indices[local_train], group_indices[local_test]
    else:
        strat = meta.iloc[group_indices]["Age"].astype(str) + "_" + meta.iloc[group_indices]["Arm"].astype(str)
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        for local_train, local_test in skf.split(np.zeros(len(group_indices)), strat):
            yield group_indices[local_train], group_indices[local_test]


def main() -> None:
    ap = argparse.ArgumentParser(description="Phase 8b: leakage-safe enhanced predictive validation")
    ap.add_argument("--phase2_dir", default=str(REPO_ROOT / "data/processed/networks/phase2"))
    ap.add_argument("--rtech", default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"))
    ap.add_argument("--meta", default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"))
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/phase8_validation"))
    ap.add_argument("--topk", type=int, default=DEFAULT_TOPK)
    ap.add_argument("--max_genes", type=int, default=2500)
    ap.add_argument("--n_perms", type=int, default=1000)
    ap.add_argument("--n_splits", type=int, default=5)
    ap.add_argument("--pool_modes", default=",".join(POOL_MODES))
    ap.add_argument("--network_methods", default="lioness")
    ap.add_argument("--network_feature_sets", default=",".join(NETWORK_FEATURE_SETS))
    ap.add_argument("--id_map", default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"))
    ap.add_argument("--gene_sets", default=str(REPO_ROOT / "config/gene_sets.yaml"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    phase2_dir = Path(args.phase2_dir)

    rtech = pd.read_csv(args.rtech, sep="\t", compression="gzip", index_col=0)
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = normalize_labels(meta.set_index(sample_col, drop=False))
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]
    mask = meta["EnvGroup"].isin(["FLT", "GC"])
    meta = meta[mask].copy()
    rtech = rtech[meta.index]

    gene_var = rtech.var(axis=1).sort_values(ascending=False)
    keep_genes = gene_var.head(args.max_genes).index.tolist()
    p2_genes_path = phase2_dir / "phase2_genes.txt"
    if p2_genes_path.exists():
        for g in [x.strip() for x in p2_genes_path.read_text().splitlines() if x.strip()]:
            if g in rtech.index and g not in keep_genes:
                keep_genes.append(g)
    rtech = rtech.loc[keep_genes]
    genes = list(rtech.index)
    X_expr = rtech.values.astype(np.float64)
    y_all = (meta["EnvGroup"] == "FLT").astype(int).values

    pool_modes = [m.strip() for m in args.pool_modes.split(",") if m.strip()]
    network_methods = [m.strip() for m in args.network_methods.split(",") if m.strip()]
    feature_sets = [m.strip() for m in args.network_feature_sets.split(",") if m.strip()]
    pathway_masks = load_curated_pathway_masks(Path(args.gene_sets), Path(args.id_map), genes)

    results = []
    rng = np.random.default_rng(42)
    groups = evaluation_groups(meta)

    for pool_mode in pool_modes:
        for method in network_methods:
            for group_name, group_idx in groups.items():
                if len(set(y_all[group_idx])) < 2:
                    continue
                print(f"\n[{pool_mode} | {method}] {group_name}: n={len(group_idx)}")

                fold_predictions: dict[tuple[str, str], list[int]] = {}
                fold_probabilities: dict[tuple[str, str], list[float]] = {}
                fold_truth: list[int] = []
                fold_feature_counts: dict[tuple[str, str], list[int]] = {}

                for train_idx, test_idx in group_folds(group_idx, y_all, meta, args.n_splits):
                    y_train = y_all[train_idx]
                    y_test = y_all[test_idx]
                    if len(set(y_train)) < 2:
                        continue

                    pool_mask = pool_mask_for_mode(meta, group_idx, train_idx, pool_mode)
                    pool_indices = np.where(pool_mask)[0]
                    if len(pool_indices) < 4:
                        continue
                    edge_i, edge_j = build_skeleton_on_fold(
                        X_expr[:, pool_indices],
                        meta.iloc[pool_indices],
                        genes,
                        topk=args.topk,
                    )

                    try:
                        train_w = sample_specific_weights(method, X_expr, pool_mask, train_idx, edge_i, edge_j)
                        test_w = sample_specific_weights(method, X_expr, pool_mask, test_idx, edge_i, edge_j)
                    except Exception as exc:
                        results.append({
                            "pool_mode": pool_mode,
                            "network_method": method,
                            "stratum": group_name,
                            "feature_type": "network_unavailable",
                            "classifier": "",
                            "loo_accuracy": np.nan,
                            "loo_auc": np.nan,
                            "status": f"unavailable: {exc}",
                        })
                        break

                    feature_mats: dict[str, tuple[np.ndarray, np.ndarray]] = {}
                    for feature_set in feature_sets:
                        if feature_set == "expression_baseline":
                            continue
                        feature_mats[feature_set] = build_network_features_fold(
                            feature_set,
                            train_w,
                            test_w,
                            edge_i,
                            edge_j,
                            len(genes),
                            pathway_masks,
                        )
                    expr_train, expr_test = build_expression_features_fold(X_expr, train_idx, test_idx)
                    feature_mats["expression_baseline"] = (expr_train, expr_test)
                    if "combined" in feature_sets:
                        net_train, net_test = feature_mats.get("network_combined") or build_network_features_fold(
                            "network_combined", train_w, test_w, edge_i, edge_j, len(genes), pathway_masks
                        )
                        feature_mats["combined"] = (np.hstack([net_train, expr_train]), np.hstack([net_test, expr_test]))

                    for feature_name, (X_train, X_test) in feature_mats.items():
                        for classifier in ["RandomForest", "LogisticRegression"]:
                            key = (feature_name, classifier)
                            pred, prob = run_classifier(X_train, y_train, X_test, classifier)
                            fold_predictions.setdefault(key, []).append(pred)
                            fold_probabilities.setdefault(key, []).append(prob)
                            fold_feature_counts.setdefault(key, []).append(X_train.shape[1])
                    fold_truth.extend(y_test.tolist())

                truth = np.array(fold_truth)
                for key, preds in fold_predictions.items():
                    feature_name, classifier = key
                    pred_arr = np.array(preds)
                    prob_arr = np.array(fold_probabilities[key])
                    acc = accuracy_score(truth, pred_arr) if len(truth) else np.nan
                    try:
                        auc = roc_auc_score(truth, prob_arr)
                    except ValueError:
                        auc = np.nan
                    results.append({
                        "pool_mode": pool_mode,
                        "network_method": method,
                        "stratum": group_name,
                        "feature_type": feature_name,
                        "classifier": classifier,
                        "n_samples": int(len(truth)),
                        "n_flt": int(truth.sum()) if len(truth) else 0,
                        "n_gc": int((1 - truth).sum()) if len(truth) else 0,
                        "n_features": int(np.median(fold_feature_counts[key])),
                        "loo_accuracy": float(acc),
                        "loo_auc": float(auc),
                        "status": "ok",
                    })

                # Conditional permutation null for network_combined + RF.
                if args.n_perms > 0 and len(truth) and ("network_combined", "RandomForest") in fold_predictions:
                    obs = accuracy_score(truth, np.array(fold_predictions[("network_combined", "RandomForest")]))
                    null_acc = []
                    pred_fixed = np.array(fold_predictions[("network_combined", "RandomForest")])
                    for _ in range(args.n_perms):
                        perm_truth = truth.copy()
                        rng.shuffle(perm_truth)
                        null_acc.append(accuracy_score(perm_truth, pred_fixed))
                    perm_p = (1 + np.sum(np.array(null_acc) >= obs)) / (args.n_perms + 1)
                    results.append({
                        "pool_mode": pool_mode,
                        "network_method": method,
                        "stratum": group_name,
                        "feature_type": "network_combined_permutation_null",
                        "classifier": "RandomForest",
                        "n_samples": int(len(truth)),
                        "n_features": 0,
                        "loo_accuracy": float(obs),
                        "loo_auc": np.nan,
                        "permutation_p_value": float(perm_p),
                        "status": "conditional_label_shuffle",
                    })

    results_df = pd.DataFrame(results)
    results_path = outdir / "enhanced_cv_results.tsv"
    results_df.to_csv(results_path, sep="\t", index=False)

    summary_rows = []
    ok = results_df[results_df["status"] == "ok"].copy()
    for (pool_mode, method, stratum), sub in ok.groupby(["pool_mode", "network_method", "stratum"]):
        best = sub.loc[sub["loo_accuracy"].idxmax()]
        net = sub[
            (sub["feature_type"] == "network_combined") &
            (sub["classifier"] == "RandomForest")
        ]["loo_accuracy"]
        expr = sub[
            (sub["feature_type"] == "expression_baseline") &
            (sub["classifier"] == "RandomForest")
        ]["loo_accuracy"]
        net_acc = float(net.iloc[0]) if len(net) else np.nan
        expr_acc = float(expr.iloc[0]) if len(expr) else np.nan
        summary_rows.append({
            "pool_mode": pool_mode,
            "network_method": method,
            "stratum": stratum,
            "best_feature_type": best["feature_type"],
            "best_classifier": best["classifier"],
            "best_accuracy": float(best["loo_accuracy"]),
            "network_rf_accuracy": net_acc,
            "expression_rf_accuracy": expr_acc,
            "network_advantage": net_acc - expr_acc if not (np.isnan(net_acc) or np.isnan(expr_acc)) else np.nan,
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "enhanced_cv_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    metadata = {
        "analysis": "Phase 8b leakage-safe enhanced CV",
        "fold_safe_transforms": [
            "network pool",
            "skeleton",
            "sample-specific network weights",
            "feature selection",
            "PCA",
            "scaling",
            "classifier fit",
        ],
        "pool_modes": pool_modes,
        "network_methods": network_methods,
        "network_feature_sets": feature_sets,
        "expression_baseline_preserved": True,
        "rf_max_depth": DEFAULT_RF_MAX_DEPTH,
        "topk": args.topk,
        "n_permutations": args.n_perms,
    }
    (outdir / "enhanced_cv_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    print(f"\n[OK] Enhanced leakage-safe CV outputs written to {outdir}")
    print(f"  - {results_path.name}")
    print(f"  - {summary_path.name}")
    print("  - enhanced_cv_metadata.json")


if __name__ == "__main__":
    main()
