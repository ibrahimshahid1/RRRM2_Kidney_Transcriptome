"""
Leakage-Safe Cross-Validation Framework for RRRM-2

Implements fold-wise procedures where all feature engineering
is recomputed within training folds to prevent information leakage.

Key Principle (Non-negotiable):
    Test samples may NOT influence:
    - Regression coefficients for residualization
    - Skeleton edge set E (including cell-wise standardization)
    - Anchor gene selection
    - Candidate feature selection (silent shifter lists)
    - Pooled networks used for LIONESS
    - Edge-wise models for predicted networks
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict, Optional, Callable
from sklearn.model_selection import StratifiedKFold
from dataclasses import dataclass


@dataclass
class CVFold:
    """Represents a single cross-validation fold."""
    fold_idx: int
    train_samples: np.ndarray
    test_samples: np.ndarray
    train_meta: pd.DataFrame
    test_meta: pd.DataFrame


class LeakageSafeCV:
    """
    Leakage-safe cross-validation for RRRM-2 network rewiring analysis.
    
    Ensures all feature engineering is performed within training folds only.
    """
    
    def __init__(
        self,
        n_folds: int = 5,
        stratify_cols: List[str] = None,
        random_state: int = 42
    ):
        """
        Initialize CV framework.
        
        Args:
            n_folds: Number of CV folds
            stratify_cols: Columns to stratify by (e.g., ['Age', 'Arm', 'EnvGroup'])
            random_state: Random seed for reproducibility
        """
        self.n_folds = n_folds
        self.stratify_cols = stratify_cols or ['Age', 'EnvGroup']
        self.random_state = random_state
    
    def create_folds(
        self,
        meta: pd.DataFrame,
        sample_col: str = 'Sample Name'
    ) -> List[CVFold]:
        """
        Create stratified CV folds.
        
        Args:
            meta: Metadata DataFrame
            sample_col: Column containing sample identifiers
            
        Returns:
            List of CVFold objects
        """
        # Create stratification key
        strat_key = meta[self.stratify_cols].astype(str).agg('_'.join, axis=1)
        
        skf = StratifiedKFold(
            n_splits=self.n_folds,
            shuffle=True,
            random_state=self.random_state
        )
        
        samples = meta[sample_col].values
        folds = []
        
        for fold_idx, (train_idx, test_idx) in enumerate(skf.split(samples, strat_key)):
            train_samples = samples[train_idx]
            test_samples = samples[test_idx]
            
            fold = CVFold(
                fold_idx=fold_idx,
                train_samples=train_samples,
                test_samples=test_samples,
                train_meta=meta.iloc[train_idx].copy(),
                test_meta=meta.iloc[test_idx].copy()
            )
            folds.append(fold)
        
        return folds
    
    def run_fold(
        self,
        fold: CVFold,
        rtech: pd.DataFrame,
        pipeline_fn: Callable,
        **pipeline_kwargs
    ) -> Dict:
        """
        Run pipeline on a single fold with leakage-safe feature engineering.
        
        Args:
            fold: CVFold object
            rtech: Residualized expression (genes x samples)
            pipeline_fn: Function to run on fold (receives train/test data)
            **pipeline_kwargs: Additional arguments for pipeline
            
        Returns:
            Dictionary with fold results
        """
        # Extract train/test data
        train_rtech = rtech[[s for s in fold.train_samples if s in rtech.columns]]
        test_rtech = rtech[[s for s in fold.test_samples if s in rtech.columns]]
        
        # Run pipeline (all feature engineering done within)
        results = pipeline_fn(
            train_rtech=train_rtech,
            test_rtech=test_rtech,
            train_meta=fold.train_meta,
            test_meta=fold.test_meta,
            fold_idx=fold.fold_idx,
            **pipeline_kwargs
        )
        
        return results
    
    def run_cv(
        self,
        meta: pd.DataFrame,
        rtech: pd.DataFrame,
        pipeline_fn: Callable,
        sample_col: str = 'Sample Name',
        **pipeline_kwargs
    ) -> List[Dict]:
        """
        Run full cross-validation.
        
        Args:
            meta: Metadata DataFrame
            rtech: Residualized expression (genes x samples)
            pipeline_fn: Function to run on each fold
            sample_col: Column containing sample identifiers
            **pipeline_kwargs: Additional arguments for pipeline
            
        Returns:
            List of fold results
        """
        folds = self.create_folds(meta, sample_col)
        
        all_results = []
        for fold in folds:
            print(f"\n{'='*60}")
            print(f"Running Fold {fold.fold_idx + 1}/{self.n_folds}")
            print(f"  Train: {len(fold.train_samples)} samples")
            print(f"  Test: {len(fold.test_samples)} samples")
            print(f"{'='*60}")
            
            results = self.run_fold(fold, rtech, pipeline_fn, **pipeline_kwargs)
            results['fold_idx'] = fold.fold_idx
            all_results.append(results)
        
        return all_results


def fold_wise_skeleton_construction(
    train_rtech: pd.DataFrame,
    train_meta: pd.DataFrame,
    cell_cols: List[str],
    max_genes: int = 1200,
    topk: int = 80
) -> Tuple[List[str], pd.DataFrame]:
    """
    Build skeleton E using only training data.
    
    This is the key leakage-safe operation: topology is learned
    from training samples only.
    
    Args:
        train_rtech: Training expression (genes x samples)
        train_meta: Training metadata
        cell_cols: Columns defining experimental cells
        max_genes: Maximum genes for network
        topk: Top-k neighbors per gene
        
    Returns:
        Tuple of (selected_genes, skeleton_edges)
    """
    from sklearn.covariance import LedoitWolf
    
    # Select top variable genes
    gene_var = train_rtech.var(axis=1)
    selected_genes = gene_var.nlargest(max_genes).index.tolist()
    
    # Cell-standardize within training
    X = train_rtech.loc[selected_genes].T  # samples x genes
    samples = X.index.tolist()
    
    # Create cell ID for each sample
    sample_to_meta = train_meta.set_index(
        train_meta.columns[0] if 'Sample Name' not in train_meta.columns 
        else 'Sample Name'
    )
    
    Z = np.zeros_like(X.values, dtype=np.float64)
    
    # Standardize within each cell
    for cell_vals, cell_df in train_meta.groupby(cell_cols):
        cell_samples = [s for s in cell_df.iloc[:, 0].values if s in samples]
        if len(cell_samples) < 2:
            continue
        cell_idx = [samples.index(s) for s in cell_samples]
        cell_data = X.iloc[cell_idx].values
        mu = cell_data.mean(axis=0)
        sd = cell_data.std(axis=0) + 1e-8
        Z[cell_idx, :] = (cell_data - mu) / sd
    
    # Compute partial correlations using Ledoit-Wolf
    lw = LedoitWolf()
    lw.fit(Z)
    precision = lw.precision_
    
    # Convert to partial correlations
    d = np.sqrt(np.diag(precision))
    d[d < 1e-8] = 1e-8
    pcorr = -precision / np.outer(d, d)
    np.fill_diagonal(pcorr, 0)
    
    # Top-k skeleton
    abs_pc = np.abs(pcorr)
    edges = set()
    G = len(selected_genes)
    
    for i in range(G):
        top_idx = np.argpartition(abs_pc[i], -topk)[-topk:]
        for j in top_idx:
            if i != j:
                a, b = (i, j) if i < j else (j, i)
                edges.add((a, b))
    
    # Create edge DataFrame
    edge_list = []
    for i, j in sorted(edges):
        edge_list.append({
            'gene_i': selected_genes[i],
            'gene_j': selected_genes[j],
            'pcorr': pcorr[i, j]
        })
    
    skeleton_edges = pd.DataFrame(edge_list)
    
    return selected_genes, skeleton_edges


def fold_wise_lioness(
    train_rtech: pd.DataFrame,
    test_rtech: pd.DataFrame,
    genes: List[str],
    skeleton_edges: pd.DataFrame
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute LIONESS for test samples relative to training pool.
    
    Args:
        train_rtech: Training expression (genes x samples)
        test_rtech: Test expression (genes x samples)
        genes: Gene list for network
        skeleton_edges: Edge DataFrame with gene_i, gene_j
        
    Returns:
        Tuple of (train_lioness, test_lioness) matrices (samples x edges)
    """
    # Align to genes
    train_X = train_rtech.loc[genes].values  # genes x train_samples
    test_X = test_rtech.loc[genes].values    # genes x test_samples
    
    N_train = train_X.shape[1]
    N_test = test_X.shape[1]
    
    # Build edge indices
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    ii = skeleton_edges['gene_i'].map(gene_to_idx).values.astype(int)
    jj = skeleton_edges['gene_j'].map(gene_to_idx).values.astype(int)
    E = len(ii)
    
    # Pooled training correlations
    def edge_corrs(X):
        """Compute edge correlations (X is genes x samples)."""
        X_std = (X - X.mean(axis=1, keepdims=True)) / (X.std(axis=1, keepdims=True) + 1e-8)
        corrs = np.array([
            (X_std[i] * X_std[j]).mean() 
            for i, j in zip(ii, jj)
        ])
        return corrs
    
    r_all = edge_corrs(train_X)
    
    # Training LIONESS
    train_lioness = np.zeros((N_train, E), dtype=np.float32)
    for s in range(N_train):
        X_loo = np.delete(train_X, s, axis=1)
        r_loo = edge_corrs(X_loo)
        train_lioness[s, :] = N_train * r_all - (N_train - 1) * r_loo
    
    # Test LIONESS: each test sample influences only its own network
    # relative to the training pool
    test_lioness = np.zeros((N_test, E), dtype=np.float32)
    for s in range(N_test):
        # Add test sample to pool, compute correlation, then LIONESS
        X_pool = np.hstack([train_X, test_X[:, s:s+1]])
        r_pool = edge_corrs(X_pool)
        test_lioness[s, :] = (N_train + 1) * r_pool - N_train * r_all
    
    return train_lioness, test_lioness


if __name__ == "__main__":
    print("Leakage-Safe Cross-Validation Module")
    print("Available classes:")
    print("  - LeakageSafeCV: Main CV framework")
    print("  - CVFold: Single fold container")
    print("Available functions:")
    print("  - fold_wise_skeleton_construction: Build skeleton from training only")
    print("  - fold_wise_lioness: Compute LIONESS with test relative to train pool")
