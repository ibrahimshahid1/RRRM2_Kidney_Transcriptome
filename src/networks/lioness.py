"""
LIONESS: Sample-Specific Network Inference

Implements sample-specific edge weight estimation on a fixed edge skeleton E.
Uses the LIONESS formula to quantify how each sample influences pooled networks.

Reference:
    Kuijjer et al. (2019) Estimating Sample-Specific Regulatory Networks. iScience.

Key Functions:
    - edge_corrs: Compute correlations for edges in E
    - lioness_sparse: Sample-specific weights via LIONESS formula
    - fisher_z: Fisher z-transform for edge weights
"""

import numpy as np
from typing import List, Tuple


def edge_corrs(
    R_std: np.ndarray,
    edges: List[Tuple[int, int]]
) -> np.ndarray:
    """
    Pearson correlation on standardized columns for each edge in E.
    
    Parameters
    ----------
    R_std : np.ndarray, shape (n_samples, n_genes)
        Standardized expression (mean=0, sd=1 per gene)
    edges : list of (int, int)
        Edge list E
        
    Returns
    -------
    W : np.ndarray, shape (|E|,)
        Correlation weight for each edge
        
    Notes
    -----
    For standardized data, Pearson(i,j) = mean(x_i * x_j)
    """
    W = np.empty(len(edges), dtype=np.float32)
    for idx, (i, j) in enumerate(edges):
        W[idx] = (R_std[:, i] * R_std[:, j]).mean()
    return W


def lioness_sparse(
    R_std: np.ndarray,
    edges: List[Tuple[int, int]],
    show_progress: bool = True
) -> np.ndarray:
    """
    Naive LIONESS on sparse edges.
    
    Parameters
    ----------
    R_std : np.ndarray, shape (n_samples, n_genes)
        Standardized expression
    edges : list of (int, int)
        Fixed edge list E
    show_progress : bool
        If True, print progress
        
    Returns
    -------
    W_samp : np.ndarray, shape (n_samples, |E|)
        Sample-specific correlations for each edge
        
    Notes
    -----
    LIONESS formula: w_e(s) = N * w_e(all) - (N-1) * w_e(-s)
    
    OPTIMIZATION NOTE: For production, vectorize LOO operations or use
    fast update formulas instead of recomputing from scratch.
    """
    N = R_std.shape[0]
    n_edges = len(edges)
    
    # Pooled weights
    w_all = edge_corrs(R_std, edges)
    
    # Per-sample weights
    W = np.empty((N, n_edges), dtype=np.float32)
    
    for s in range(N):
        if show_progress and s % 10 == 0:
            print(f"LIONESS: processing sample {s+1}/{N}")
        
        # Leave-one-out
        R_loo = np.delete(R_std, s, axis=0)
        
        # Re-standardize LOO (IMPORTANT for correct correlations)
        R_loo -= R_loo.mean(axis=0, keepdims=True)
        R_loo /= (R_loo.std(axis=0, keepdims=True) + 1e-8)
        
        w_minus = edge_corrs(R_loo, edges)
        
        # LIONESS formula
        W[s, :] = N * w_all - (N - 1) * w_minus
    
    return W


def lioness_sparse_fast(
    R_std: np.ndarray,
    edges: List[Tuple[int, int]]
) -> np.ndarray:
    """
    Fast LIONESS using update formulas (advanced).
    
    Parameters
    ----------
    R_std : np.ndarray
        Standardized expression
    edges : list of (int, int)
        Edge list
        
    Returns
    -------
    W_samp : np.ndarray
        Sample-specific weights
        
    Notes
    -----
    TODO: Implement fast LOO update formulas that avoid full recomputation.
    For now, use lioness_sparse (slower but correct).
    """
    raise NotImplementedError("Fast LIONESS not yet implemented. Use lioness_sparse.")


def fisher_z(W: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """
    Fisher transform of correlations: z = atanh(r).
    
    Parameters
    ----------
    W : np.ndarray
        Correlation weights
    eps : float
        Small constant to clip correlations away from ±1
        
    Returns
    -------
    Z : np.ndarray
        Fisher z-transformed weights
        
    Notes
    -----
    Fisher z approximates normal distribution, improving regression modeling.
    Clips to avoid inf at |r|=1.
    """
    Wc = np.clip(W, -1 + eps, 1 - eps)
    return np.arctanh(Wc).astype(np.float32)


def inverse_fisher_z(Z: np.ndarray) -> np.ndarray:
    """
    Inverse Fisher transform: r = tanh(z).
    
    Parameters
    ----------
    Z : np.ndarray
        Fisher z-transformed weights
        
    Returns
    -------
    W : np.ndarray
        Correlation weights
    """
    return np.tanh(Z).astype(np.float32)


if __name__ == "__main__":
    # Example usage with toy data
    np.random.seed(42)
    R = np.random.randn(20, 10).astype(np.float32)
    
    # Standardize
    R -= R.mean(axis=0, keepdims=True)
    R /= (R.std(axis=0, keepdims=True) + 1e-8)
    
    # Toy edge list
    edges = [(0, 1), (1, 2), (2, 3), (0, 5)]
    
    print(f"Computing LIONESS for {len(edges)} edges on {R.shape[0]} samples...")
    W_samp = lioness_sparse(R, edges, show_progress=False)
    Z_samp = fisher_z(W_samp)
    
    print(f"Output shape: {W_samp.shape}")
    print(f"Sample 0 edge weights (raw): {W_samp[0, :]}")
    print(f"Sample 0 edge weights (Fisher-z): {Z_samp[0, :]}")
    print("\nLIONESS module loaded successfully")
