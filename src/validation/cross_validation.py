"""
Leakage-Safe Cross-Validation Framework

Implements fold-wise operations that prevent test samples from influencing
feature space definition. Critical for honest predictive validation.

Key Principle:
    All steps that define features (residualization, topology E, anchors, 
    candidate selection) must be recomputed within each training fold.

Key Functions:
    - create_cv_folds: Stratified fold creation
    - cv_fold_pipeline: Run full pipeline on one fold
    - evaluate_classification: Sample-level prediction
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Callable
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, roc_auc_score, classification_report


def create_cv_folds(
    metadata: pd.DataFrame,
    stratify_col: str = 'EnvironmentGroup',
    n_folds: int = 5,
    random_state: int = 42
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """
    Create stratified CV folds.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Sample metadata
    stratify_col : str
        Column to stratify on
    n_folds : int
        Number of folds
    random_state : int
        Random seed
        
    Returns
    -------
    folds : list of (train_idx, test_idx)
        Fold indices
    """
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_state)
    
    y_stratify = metadata[stratify_col].values
    
    folds = []
    for train_idx, test_idx in skf.split(np.zeros(len(metadata)), y_stratify):
        folds.append((train_idx, test_idx))
    
    return folds


def cv_fold_pipeline(
    fold_idx: int,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    Y_raw: np.ndarray,
    metadata: pd.DataFrame,
    pipeline_functions: Dict[str, Callable]
) -> Dict:
    """
    Run full pipeline on one CV fold (leakage-safe).
    
    Parameters
    ----------
    fold_idx : int
        Fold number
    train_idx, test_idx : np.ndarray
        Sample indices
    Y_raw : np.ndarray
        Raw expression (samples × genes)
    metadata : pd.DataFrame
        Sample metadata
    pipeline_functions : dict
        {
            'residualize': function(Y_train, meta_train) -> coefficients,
            'build_topology': function(R_train, meta_train) -> E,
            'select_anchors': function(R_train, meta_train) -> anchor_indices,
            'extract_features': function(R_sample, E, anchors) -> features
        }
        
    Returns
    -------
    results : dict
        {
            'fold': fold_idx,
            'train_features': np.ndarray,
            'test_features': np.ndarray,
            'train_labels': np.ndarray,
            'test_labels': np.ndarray,
            'edge_list': E,
            'anchors': anchor_indices
        }
        
    Notes
    -----
    This is a TEMPLATE. User must provide actual pipeline functions that:
    - Fit residualization on train only
    - Build E from train only (with cell-standardization)
    - Select anchors from train only
    - Extract LIONESS-based features per sample
    """
    print(f"\n=== Fold {fold_idx} ===")
    
    # Split data
    Y_train = Y_raw[train_idx, :]
    Y_test = Y_raw[test_idx, :]
    meta_train = metadata.iloc[train_idx]
    meta_test = metadata.iloc[test_idx]
    
    # 1. Residualize (fit on train, apply to both)
    print("  1. Residualization...")
    residualize_fn = pipeline_functions['residualize']
    R_train, R_test = residualize_fn(Y_train, Y_test, meta_train, meta_test)
    
    # 2. Build topology E (from train only)
    print("  2. Building shared topology E...")
    build_topology_fn = pipeline_functions['build_topology']
    E = build_topology_fn(R_train, meta_train)
    
    # 3. Select anchors (from train only)
    print("  3. Selecting anchor genes...")
    select_anchors_fn = pipeline_functions['select_anchors']
    anchor_indices = select_anchors_fn(R_train, meta_train)
    
    # 4. Extract sample-level features via LIONESS
    print("  4. Extracting LIONESS features...")
    extract_features_fn = pipeline_functions['extract_features']
    train_features = extract_features_fn(R_train, E, anchor_indices)
    test_features = extract_features_fn(R_test, E, anchor_indices)
    
    # Labels
    train_labels = meta_train['EnvironmentGroup'].values
    test_labels = meta_test['EnvironmentGroup'].values
    
    results = {
        'fold': fold_idx,
        'train_features': train_features,
        'test_features': test_features,
        'train_labels': train_labels,
        'test_labels': test_labels,
        'edge_list': E,
        'anchors': anchor_indices,
        'train_idx': train_idx,
        'test_idx': test_idx
    }
    
    return results


def evaluate_classification(
    cv_results: List[Dict],
    classifier_class = RandomForestClassifier,
    classifier_params: Dict = None
) -> pd.DataFrame:
    """
    Train classifiers on CV folds and evaluate.
    
    Parameters
    ----------
    cv_results : list of dict
        Output from cv_fold_pipeline for each fold
    classifier_class : sklearn classifier
        Classifier class to use
    classifier_params : dict
        Parameters for classifier
        
    Returns
    -------
    metrics : pd.DataFrame
        Per-fold performance metrics
    """
    if classifier_params is None:
        classifier_params = {'n_estimators': 100, 'random_state': 42}
    
    fold_metrics = []
    
    for fold_result in cv_results:
        fold_idx = fold_result['fold']
        
        X_train = fold_result['train_features']
        y_train = fold_result['train_labels']
        X_test = fold_result['test_features']
        y_test = fold_result['test_labels']
        
        # Train classifier
        clf = classifier_class(**classifier_params)
        clf.fit(X_train, y_train)
        
        # Predict
        y_pred = clf.predict(X_test)
        y_proba = clf.predict_proba(X_test)
        
        # Metrics
        acc = accuracy_score(y_test, y_pred)
        
        # Multiclass AUC (one-vs-rest)
        try:
            auc = roc_auc_score(y_test, y_proba, multi_class='ovr', average='weighted')
        except:
            auc = np.nan
        
        fold_metrics.append({
            'fold': fold_idx,
            'accuracy': acc,
            'auc': auc,
            'n_train': len(y_train),
            'n_test': len(y_test)
        })
        
        print(f"Fold {fold_idx}: Accuracy = {acc:.3f}, AUC = {auc:.3f}")
    
    metrics_df = pd.DataFrame(fold_metrics)
    
    print("\n=== Cross-Validation Summary ===")
    print(f"Mean Accuracy: {metrics_df['accuracy'].mean():.3f} ± {metrics_df['accuracy'].std():.3f}")
    print(f"Mean AUC: {metrics_df['auc'].mean():.3f} ± {metrics_df['auc'].std():.3f}")
    
    return metrics_df


if __name__ == "__main__":
    print("Cross-validation module loaded successfully")
    print("Key functions: create_cv_folds, cv_fold_pipeline, evaluate_classification")
