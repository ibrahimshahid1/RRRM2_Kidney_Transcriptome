"""
Shared Topology Construction with Cell-Standardization

Implements Phase 2 of the pipeline: learning a single global edge list E
using cell-standardized pooled expression to prevent Simpson's paradox.

Key Innovation:
    Standardize within each experimental cell (Age × Arm × Group) before
    pooling for topology selection. This prevents edges from being selected
    due to between-cell mean shifts.

Key Functions:
    - standardize_within_cells: Cell-wise standardization per gene
    - standardize_global: Global standardization for weight computation
    - rank_transform: Spearman-like robustness
    - build_global_skeleton_topk: Top-k edge selection
    - build_skeleton_partial_corr: Partial correlation skeleton
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Optional
from sklearn.covariance import LedoitWolf
from scipy import linalg


def standardize_within_cells(
    R: np.ndarray,
    cell_ids: np.ndarray
) -> np.ndarray:
    """
    Cell-wise standardization per gene for TOPOLOGY SELECTION ONLY.
    
    Parameters
    ----------
    R : np.ndarray, shape (n_samples, n_genes)
        Residualized expression matrix
    cell_ids : np.ndarray, shape (n_samples,)
        Integer cell label for each sample (encode Age×Arm×Group)
        
    Returns
    -------
    R_star : np.ndarray
        Cell-standardized matrix for topology selection
        
    Notes
    -----
    This is the KEY innovation to prevent Simpson's paradox in edge selection.
    After standardization, pool all samples to learn edges on standardized data.
    """
    R = R.astype(np.float32)
    out = np.empty_like(R, dtype=np.float32)
    
    for cid in np.unique(cell_ids):
        idx = np.where(cell_ids == cid)[0]
        Rc = R[idx, :]
        mu = Rc.mean(axis=0, keepdims=True)
        sd = Rc.std(axis=0, keepdims=True)
        out[idx, :] = (Rc - mu) / (sd + 1e-8)
    
    return out


def standardize_global(R: np.ndarray) -> np.ndarray:
    """
    Standardize per gene across all samples (for weight computations).
    
    Parameters
    ----------
    R : np.ndarray
        Expression matrix
        
    Returns
    -------
    R_std : np.ndarray
        Globally standardized matrix
    """
    R = R.astype(np.float32)
    R -= R.mean(axis=0, keepdims=True)
    R /= (R.std(axis=0, keepdims=True) + 1e-8)
    return R


def rank_transform_inplace(R: np.ndarray) -> np.ndarray:
    """
    Spearman-like: rank-transform per gene, then Pearson on ranks.
    
    Parameters
    ----------
    R : np.ndarray
        Expression matrix
        
    Returns
    -------
    R_ranked : np.ndarray
        Rank-transformed matrix
    """
    out = np.empty_like(R, dtype=np.float32)
    for g in range(R.shape[1]):
        order = np.argsort(R[:, g])
        ranks = np.empty_like(order, dtype=np.float32)
        ranks[order] = np.arange(1, R.shape[0] + 1, dtype=np.float32)
        out[:, g] = ranks
    return out


def build_global_skeleton_topk(
    R_for_edges: np.ndarray,
    k: int = 50,
    use_abs: bool = True
) -> List[Tuple[int, int]]:
    """
    Build global edge list via top-k neighbors per gene using pooled data.
    
    Parameters
    ----------
    R_for_edges : np.ndarray
        Should be cell-standardized (and optionally rank-transformed)
    k : int
        Number of top neighbors per gene
    use_abs : bool
        If True, use |correlation| for ranking
        
    Returns
    -------
    edges : list of (i, j) tuples
        Undirected edge list (i < j)
        
    Notes
    -----
    For production, consider using partial correlation or graphical lasso
    instead of simple correlation.
    """
    # Correlation on standardized data: C = (X^T X) / (n-1)
    C = (R_for_edges.T @ R_for_edges) / (R_for_edges.shape[0] - 1)
    np.fill_diagonal(C, 0.0)
    
    if use_abs:
        absC = np.abs(C)
    else:
        absC = C
    
    edges = set()
    for i in range(absC.shape[0]):
        # Get top-k neighbors (excluding self)
        nbrs = np.argpartition(absC[i], -k)[-k:]
        for j in nbrs:
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            edges.add((a, b))
    
    return list(edges)


def build_skeleton_partial_corr(
    R_for_edges: np.ndarray,
    threshold: float = 0.1,
    max_edges: Optional[int] = None
) -> List[Tuple[int, int]]:
    """
    Build skeleton using partial correlations from precision matrix.
    
    Parameters
    ----------
    R_for_edges : np.ndarray
        Cell-standardized expression
    threshold : float
        Minimum |partial correlation| to keep edge
    max_edges : int, optional
        Maximum edges to return (keeps strongest)
        
    Returns
    -------
    edges : list of (i, j)
        Edge list
        
    Notes
    -----
    Uses Ledoit-Wolf shrinkage for well-conditioned covariance estimation.
    For large gene sets, consider graphical lasso for sparsity.
    """
    # Estimate shrunk covariance
    lw = LedoitWolf()
    cov = lw.fit(R_for_edges).covariance_
    
    # Invert to precision
    try:
        precision = linalg.inv(cov)
    except linalg.LinAlgError:
        # Fallback to pseudoinverse
        precision = linalg.pinv(cov)
    
    # Convert to partial correlations
    D = np.sqrt(np.diag(precision))
    D_outer = np.outer(D, D)
    partial_corr = -precision / D_outer
    np.fill_diagonal(partial_corr, 0.0)
    
    # Threshold
    abs_pcor = np.abs(partial_corr)
    edges = []
    for i in range(abs_pcor.shape[0]):
        for j in range(i + 1, abs_pcor.shape[1]):
            if abs_pcor[i, j] >= threshold:
                edges.append((i, j))
    
    # Sort by strength and limit
    if max_edges is not None and len(edges) > max_edges:
        edge_strengths = [(i, j, abs_pcor[i, j]) for i, j in edges]
        edge_strengths.sort(key=lambda x: x[2], reverse=True)
        edges = [(i, j) for i, j, _ in edge_strengths[:max_edges]]
    
    return edges


def save_edge_list(
    edges: List[Tuple[int, int]],
    gene_names: List[str],
    output_path: str
) -> None:
    """
    Save edge list to file with gene names.
    
    Parameters
    ----------
    edges : list of (int, int)
        Edge indices
    gene_names : list of str
        Gene identifiers
    output_path : str
        Output file path
    """
    edge_df = pd.DataFrame([
        {'gene1': gene_names[i], 'gene2': gene_names[j]}
        for i, j in edges
    ])
    edge_df.to_csv(output_path, index=False)
    print(f"Saved {len(edges)} edges to {output_path}")


def load_edge_list(
    input_path: str,
    gene_names: List[str]
) -> List[Tuple[int, int]]:
    """
    Load edge list from file and convert to indices.
    
    Parameters
    ----------
    input_path : str
        Input file path
    gene_names : list of str
        Gene identifiers for index mapping
        
    Returns
    -------
    edges : list of (int, int)
        Edge indices
    """
    edge_df = pd.read_csv(input_path)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    
    edges = []
    for _, row in edge_df.iterrows():
        i = gene_to_idx.get(row['gene1'])
        j = gene_to_idx.get(row['gene2'])
        if i is not None and j is not None:
            edges.append((i, j))
    
    return edges


if __name__ == "__main__":
    print("Shared topology module loaded successfully")
    print("Key functions: standardize_within_cells, build_global_skeleton_topk")
