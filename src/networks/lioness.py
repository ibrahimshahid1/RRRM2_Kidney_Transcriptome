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
from dataclasses import dataclass
from pathlib import Path


# ============================================================================
# DENSE LIONESS IMPLEMENTATION (full correlation networks)
# ============================================================================

@dataclass(frozen=True)
class LionessIndex:
    """Index structure for reconstructing matrices from edge vectors."""
    gene_names: list
    triu_i: np.ndarray  # (n_edges,) row indices of upper triangle
    triu_j: np.ndarray  # (n_edges,) column indices of upper triangle


def _corr_from_summaries(sum_x: np.ndarray, cross: np.ndarray, n: int, eps: float = 1e-12) -> np.ndarray:
    """
    Compute Pearson correlation matrix from sufficient statistics.
    
    This allows efficient leave-one-out computation without recomputing
    the full correlation matrix from scratch.
    
    Parameters
    ----------
    sum_x : np.ndarray
        Column sums of the data matrix, shape (p,)
    cross : np.ndarray
        Cross-product matrix X.T @ X, shape (p, p)
    n : int
        Number of samples
    eps : float
        Small constant to prevent division by zero
        
    Returns
    -------
    np.ndarray
        Correlation matrix, shape (p, p)
    """
    outer = np.outer(sum_x, sum_x) / float(n)
    S = cross - outer
    cov = S / float(n - 1)

    var = np.clip(np.diag(cov), eps, None)
    sd = np.sqrt(var)
    denom = sd[:, None] * sd[None, :]
    corr = cov / denom

    corr = np.clip(corr, -1.0, 1.0)
    np.fill_diagonal(corr, 0.0)
    return corr


def lioness_correlation_edges(
    X: np.ndarray,
    gene_names: list,
    dtype=np.float32,
    verbose: bool = True,
) -> tuple:
    """
    Compute LIONESS sample-specific networks using Pearson correlation.
    
    Uses full dense correlation matrices. Suitable for ~800-1500 genes.
    
    Parameters
    ----------
    X : np.ndarray
        Expression matrix, shape (n_samples, p_genes)
    gene_names : list
        Gene names corresponding to columns of X
    dtype : np.dtype
        Output dtype for edge weights
    verbose : bool
        Print progress
        
    Returns
    -------
    edges : np.ndarray, shape (n_samples, n_edges)
    index : LionessIndex
    """
    if X.ndim != 2:
        raise ValueError(f"X must be 2D (n_samples x p_genes). Got {X.shape}")

    n, p = X.shape
    if n < 4:
        raise ValueError("Need at least 4 samples for stable correlation / leave-one-out.")

    if len(gene_names) != p:
        raise ValueError("gene_names length must match number of columns in X.")

    if verbose:
        n_edges = p * (p - 1) // 2
        print(f"LIONESS: {n} samples × {p} genes → {n_edges:,} edges per network")

    X = np.asarray(X, dtype=np.float64)
    sum_x = X.sum(axis=0)
    cross = X.T @ X

    W_all = _corr_from_summaries(sum_x, cross, n)

    triu_i, triu_j = np.triu_indices(p, k=1)
    W_all_edges = W_all[triu_i, triu_j]

    n_edges = W_all_edges.shape[0]
    edges = np.empty((n, n_edges), dtype=dtype)

    for i in range(n):
        if verbose and i % 10 == 0:
            print(f"  Processing sample {i+1}/{n}", end="\r")
        
        xi = X[i, :]
        sum_x_m = sum_x - xi
        cross_m = cross - np.outer(xi, xi)
        W_m = _corr_from_summaries(sum_x_m, cross_m, n - 1)
        W_m_edges = W_m[triu_i, triu_j]

        Wi_edges = n * (W_all_edges - W_m_edges) + W_m_edges
        edges[i, :] = Wi_edges.astype(dtype, copy=False)

    if verbose:
        print(f"  Completed all {n} samples" + " " * 20)

    index = LionessIndex(gene_names=gene_names, triu_i=triu_i, triu_j=triu_j)
    return edges, index


def edges_to_matrix(edges: np.ndarray, index: LionessIndex) -> np.ndarray:
    """Reconstruct symmetric adjacency matrix from edge vector."""
    p = len(index.gene_names)
    mat = np.zeros((p, p), dtype=edges.dtype)
    mat[index.triu_i, index.triu_j] = edges
    mat[index.triu_j, index.triu_i] = edges
    return mat


def save_lioness_bundle(outdir, edges: np.ndarray, index: LionessIndex) -> None:
    """Save LIONESS results: lioness_edges.npy, triu_i.npy, triu_j.npy, genes.txt"""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    np.save(outdir / "lioness_edges.npy", edges)
    np.save(outdir / "triu_i.npy", index.triu_i)
    np.save(outdir / "triu_j.npy", index.triu_j)
    (outdir / "genes.txt").write_text("\n".join(index.gene_names) + "\n", encoding="utf-8")
    
    print(f"Saved LIONESS bundle to {outdir}:")
    print(f"  - lioness_edges.npy: {edges.shape}")
    print(f"  - triu_i.npy, triu_j.npy: {len(index.triu_i):,} edges")
    print(f"  - genes.txt: {len(index.gene_names)} genes")


def load_lioness_bundle(indir):
    """Load LIONESS results from disk."""
    indir = Path(indir)
    edges = np.load(indir / "lioness_edges.npy")
    triu_i = np.load(indir / "triu_i.npy")
    triu_j = np.load(indir / "triu_j.npy")
    gene_names = (indir / "genes.txt").read_text(encoding="utf-8").strip().split("\n")
    index = LionessIndex(gene_names=gene_names, triu_i=triu_i, triu_j=triu_j)
    return edges, index


# ============================================================================
# SPARSE LIONESS IMPLEMENTATION (original - for fixed edge skeletons)
# ============================================================================


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
