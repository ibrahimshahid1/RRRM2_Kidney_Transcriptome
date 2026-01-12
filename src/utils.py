"""
Utility Functions for RRRM-2 Kidney Network Analysis

Consolidated starter code and helper functions from the research proposal.
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Tuple, Optional
from sklearn.linear_model import LinearRegression


def residualize_expr(
    Y_samples_by_genes: np.ndarray,
    X_covariates: np.ndarray
) -> np.ndarray:
    """
    Residualize each gene on covariates (global regression).
    
    See src/preprocessing/residualization.py for full implementation.
    """
    lr = LinearRegression()
    R = np.empty_like(Y_samples_by_genes, dtype=np.float32)
    
    for g in range(Y_samples_by_genes.shape[1]):
        lr.fit(X_covariates, Y_samples_by_genes[:, g])
        R[:, g] = (Y_samples_by_genes[:, g] - lr.predict(X_covariates)).astype(np.float32)
    
    return R


def standardize_global(R: np.ndarray) -> np.ndarray:
    """
    Standardize per gene across all samples (for weight computations).
    """
    R = R.astype(np.float32)
    R -= R.mean(axis=0, keepdims=True)
    R /= (R.std(axis=0, keepdims=True) + 1e-8)
    return R


def standardize_within_cells(
    R: np.ndarray,
    cell_ids: np.ndarray
) -> np.ndarray:
    """
    Cell-wise standardization per gene for TOPOLOGY SELECTION ONLY.
    
    cell_ids: length n_samples, integer cell label for each sample
              (e.g., encode Age×Arm×Group).
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


def rank_transform_inplace(R: np.ndarray) -> np.ndarray:
    """
    Spearman-like: rank-transform per gene, then Pearson on ranks.
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
    k: int = 50
) -> List[Tuple[int, int]]:
    """
    Build global edge list via top-k neighbors per gene using pooled data.
    
    R_for_edges should already be standardized appropriately (e.g., within-cell standardized),
    and optionally rank-transformed if you want Spearman-like robustness.
    """
    # Correlation on standardized data: C = (X^T X) / (n-1)
    C = (R_for_edges.T @ R_for_edges) / (R_for_edges.shape[0] - 1)
    np.fill_diagonal(C, 0.0)
    absC = np.abs(C)
    
    edges = set()
    for i in range(absC.shape[0]):
        nbrs = np.argpartition(absC[i], -k)[-k:]
        for j in nbrs:
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            edges.add((a, b))
    
    return list(edges)


def edge_corrs(
    R_std: np.ndarray,
    edges: List[Tuple[int, int]]
) -> np.ndarray:
    """
    Pearson corr on standardized columns for each edge in edges.
    """
    W = np.empty(len(edges), dtype=np.float32)
    for idx, (i, j) in enumerate(edges):
        W[idx] = (R_std[:, i] * R_std[:, j]).mean()
    return W


def lioness_sparse(
    R_std: np.ndarray,
    edges: List[Tuple[int, int]]
) -> np.ndarray:
    """
    Naive LIONESS on sparse edges.
    
    NOTE: optimize in real code (fast LOO formulas / vectorized ops).
    Returns: W_samp (n_samples, |E|) sample-specific correlations.
    """
    N = R_std.shape[0]
    w_all = edge_corrs(R_std, edges)
    W = np.empty((N, len(edges)), dtype=np.float32)
    
    for s in range(N):
        R_loo = np.delete(R_std, s, axis=0)
        w_minus = edge_corrs(R_loo, edges)
        W[s, :] = N * w_all - (N - 1) * w_minus
    
    return W


def fisher_z(W: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """
    Fisher transform of correlations: z = atanh(r).
    Clips to avoid inf at |r|=1.
    """
    Wc = np.clip(W, -1 + eps, 1 - eps)
    return np.arctanh(Wc).astype(np.float32)


def weights_to_graph(
    n_genes: int,
    edges: List[Tuple[int, int]],
    weights: np.ndarray,
    mode: str = "pos"
) -> nx.Graph:
    """
    Convert edge list + weights to NetworkX graph.
    
    mode:
      - "pos": keep only positive weights
      - "neg": keep only negative weights as positive magnitude
      - "abs_nonneg": use abs(weights)
    """
    G = nx.Graph()
    G.add_nodes_from(range(n_genes))
    
    for (i, j), w in zip(edges, weights):
        if mode == "pos":
            if w <= 0:
                continue
            G.add_edge(i, j, weight=float(w))
        elif mode == "neg":
            if w >= 0:
                continue
            G.add_edge(i, j, weight=float(-w))
        else:
            G.add_edge(i, j, weight=float(abs(w)))
    
    return G


# Additional helper functions

def load_config(config_path: str) -> dict:
    """Load YAML configuration file."""
    import yaml
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def save_matrix(matrix: np.ndarray, path: str, format: str = 'npy') -> None:
    """Save matrix to file."""
    if format == 'npy':
        np.save(path, matrix)
    elif format == 'csv':
        pd.DataFrame(matrix).to_csv(path, index=False)
    else:
        raise ValueError(f"Unknown format: {format}")


def load_matrix(path: str, format: str = 'npy') -> np.ndarray:
    """Load matrix from file."""
    if format == 'npy':
        return np.load(path)
    elif format == 'csv':
        return pd.read_csv(path).values
    else:
        raise ValueError(f"Unknown format: {format}")


if __name__ == "__main__":
    print("Utilities module loaded successfully")
    print("Available functions:")
    print("  - residualize_expr")
    print("  - standardize_global, standardize_within_cells")
    print("  - build_global_skeleton_topk")
    print("  - lioness_sparse")
    print("  - fisher_z")
    print("  - weights_to_graph")