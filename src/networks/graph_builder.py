# src/networks/graph_builder.py
"""
Graph Builder

Builds undirected weighted graphs from upper-triangle edge vectors.

Supports:
- Global threshold (top X% by |w|)
- Top-k per node (keeps each node connected)
- Hybrid mode (union of global + topk)
- Conversion to dense adjacency or edge list
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional, Tuple
import numpy as np
import pandas as pd


ThresholdMode = Literal["global", "topk", "hybrid"]


@dataclass(frozen=True)
class GraphData:
    """
    Graph representation with dense adjacency and edge list.
    
    Attributes
    ----------
    gene_names : list[str]
        Node names in order [0..p-1]
    W : np.ndarray
        Dense symmetric adjacency (p x p), float32, diag=0
    edges : pd.DataFrame
        Edge list with columns [u, v, w, abs_w]
    """
    gene_names: list
    W: np.ndarray
    edges: pd.DataFrame


def edge_vector_to_dense(
    weights: np.ndarray,
    triu_i: np.ndarray,
    triu_j: np.ndarray,
    p: int,
    dtype=np.float32,
) -> np.ndarray:
    """
    Build dense symmetric matrix W from an upper-triangle edge vector.

    Parameters
    ----------
    weights : np.ndarray
        Edge weights, shape (n_edges,)
    triu_i, triu_j : np.ndarray
        Node indices for upper triangle (k=1), shape (n_edges,)
    p : int
        Number of genes/nodes
    dtype : np.dtype
        Output dtype
        
    Returns
    -------
    np.ndarray
        Symmetric adjacency matrix, shape (p, p)
    """
    W = np.zeros((p, p), dtype=dtype)
    W[triu_i, triu_j] = weights.astype(dtype, copy=False)
    W[triu_j, triu_i] = weights.astype(dtype, copy=False)
    np.fill_diagonal(W, 0.0)
    return W


def _global_threshold_mask(
    weights: np.ndarray,
    keep_fraction: float,
    use_abs: bool = True,
) -> np.ndarray:
    """
    Returns a boolean mask over edges for global thresholding.
    Keeps the top `keep_fraction` edges by |w| (or w if use_abs=False).
    """
    if not (0.0 < keep_fraction <= 1.0):
        raise ValueError(f"keep_fraction must be in (0,1]. Got {keep_fraction}")

    score = np.abs(weights) if use_abs else weights
    n_edges = score.shape[0]
    k = int(np.ceil(n_edges * keep_fraction))
    if k <= 0:
        return np.zeros(n_edges, dtype=bool)
    if k >= n_edges:
        return np.ones(n_edges, dtype=bool)

    # Find threshold value for top-k without full sort (partition)
    kth = np.partition(score, n_edges - k)[n_edges - k]
    mask = score >= kth
    return mask


def _topk_mask_from_dense(
    W_full: np.ndarray,
    triu_i: np.ndarray,
    triu_j: np.ndarray,
    topk: int,
) -> np.ndarray:
    """
    Create a mask over the edge vector that keeps top-k edges per node by |w|.
    
    Symmetry rule: if either endpoint selects the edge in its top-k, we keep it.
    """
    if topk <= 0:
        return np.zeros_like(triu_i, dtype=bool)

    p = W_full.shape[0]
    absW = np.abs(W_full)

    keep_mat = np.zeros((p, p), dtype=bool)

    for u in range(p):
        row = absW[u].copy()
        row[u] = 0.0
        if topk >= p - 1:
            nbrs = np.where(row > 0)[0]
        else:
            idx = np.argpartition(row, -(topk))[-topk:]
            nbrs = idx[row[idx] > 0]
        keep_mat[u, nbrs] = True

    # Symmetric union
    keep_mat = keep_mat | keep_mat.T

    # Map to edge vector mask
    mask = keep_mat[triu_i, triu_j]
    return mask


def threshold_edge_vector(
    weights: np.ndarray,
    triu_i: np.ndarray,
    triu_j: np.ndarray,
    p: int,
    gene_names: Optional[list] = None,
    mode: ThresholdMode = "topk",
    keep_fraction: float = 0.01,
    topk: int = 30,
    min_abs_weight: Optional[float] = None,
    dtype=np.float32,
) -> GraphData:
    """
    Build a graph from an edge vector using thresholding.

    Parameters
    ----------
    weights : np.ndarray
        Edge weights, shape (n_edges,)
    triu_i, triu_j : np.ndarray
        Upper triangle indices
    p : int
        Number of nodes
    gene_names : list, optional
        Node names. If None, uses g0, g1, ...
    mode : str
        "global": keep top `keep_fraction` edges by |w|
        "topk": keep top `topk` edges per node by |w| (symmetric union)
        "hybrid": union of global + topk
    keep_fraction : float
        Fraction of edges to keep in global mode
    topk : int
        Number of edges per node in topk mode
    min_abs_weight : float, optional
        If set, edges with |w| < min_abs_weight are dropped
    dtype : np.dtype
        Output dtype for adjacency

    Returns
    -------
    GraphData
        Graph with dense W and edge list df
    """
    weights = np.asarray(weights)
    triu_i = np.asarray(triu_i)
    triu_j = np.asarray(triu_j)

    n_edges = weights.shape[0]
    if triu_i.shape[0] != n_edges or triu_j.shape[0] != n_edges:
        raise ValueError("weights, triu_i, triu_j must have the same length")

    if gene_names is None:
        gene_names = [f"g{i}" for i in range(p)]

    # Build masks
    if mode == "global":
        mask = _global_threshold_mask(weights, keep_fraction=keep_fraction, use_abs=True)

    elif mode == "topk":
        W_full = edge_vector_to_dense(weights, triu_i, triu_j, p, dtype=np.float32)
        mask = _topk_mask_from_dense(W_full, triu_i, triu_j, topk=topk)

    elif mode == "hybrid":
        W_full = edge_vector_to_dense(weights, triu_i, triu_j, p, dtype=np.float32)
        mask_global = _global_threshold_mask(weights, keep_fraction=keep_fraction, use_abs=True)
        mask_topk = _topk_mask_from_dense(W_full, triu_i, triu_j, topk=topk)
        mask = mask_global | mask_topk

    else:
        raise ValueError(f"Unknown mode: {mode}")

    # Optional min_abs filter
    if min_abs_weight is not None:
        mask = mask & (np.abs(weights) >= float(min_abs_weight))

    # Build adjacency from kept edges
    kept_w = weights[mask].astype(dtype, copy=False)
    kept_i = triu_i[mask]
    kept_j = triu_j[mask]

    W = np.zeros((p, p), dtype=dtype)
    W[kept_i, kept_j] = kept_w
    W[kept_j, kept_i] = kept_w
    np.fill_diagonal(W, 0.0)

    edges_df = pd.DataFrame({
        "u": kept_i.astype(np.int32),
        "v": kept_j.astype(np.int32),
        "w": kept_w.astype(np.float32),
    })
    edges_df["abs_w"] = np.abs(edges_df["w"].values)
    edges_df = edges_df.sort_values("abs_w", ascending=False).reset_index(drop=True)

    return GraphData(gene_names=gene_names, W=W, edges=edges_df)


def ensure_connected_by_backfill(
    W: np.ndarray,
    original_W: np.ndarray,
    min_degree: int = 1,
    backfill_topk: int = 5,
) -> np.ndarray:
    """
    If thresholding produces isolated nodes (degree 0), backfill a few strongest
    edges from original_W for those nodes.

    Parameters
    ----------
    W : np.ndarray
        Thresholded adjacency
    original_W : np.ndarray
        Full/less-thresholded adjacency to source backfill edges
    min_degree : int
        Minimum degree to ensure
    backfill_topk : int
        Number of edges to backfill for isolated nodes
        
    Returns
    -------
    np.ndarray
        Adjacency with backfilled connections
    """
    p = W.shape[0]
    out = W.copy()
    abs_orig = np.abs(original_W)

    for u in range(p):
        deg = int(np.count_nonzero(out[u]))
        if deg >= min_degree:
            continue

        row = abs_orig[u].copy()
        row[u] = 0.0
        if backfill_topk >= p - 1:
            nbrs = np.where(row > 0)[0]
        else:
            idx = np.argpartition(row, -(backfill_topk))[-backfill_topk:]
            nbrs = idx[row[idx] > 0]

        for v in nbrs:
            out[u, v] = original_W[u, v]
            out[v, u] = original_W[v, u]

    np.fill_diagonal(out, 0.0)
    return out


def adjacency_to_edgelist(W: np.ndarray, gene_names: list) -> pd.DataFrame:
    """
    Convert symmetric adjacency W into an undirected edge list.
    
    Parameters
    ----------
    W : np.ndarray
        Symmetric adjacency matrix
    gene_names : list
        Node names
        
    Returns
    -------
    pd.DataFrame
        Edge list with columns [u, v, gene_u, gene_v, w, abs_w]
    """
    p = W.shape[0]
    triu_i, triu_j = np.triu_indices(p, k=1)
    w = W[triu_i, triu_j]
    mask = w != 0

    triu_i = triu_i[mask]
    triu_j = triu_j[mask]
    w = w[mask].astype(np.float32, copy=False)

    df = pd.DataFrame({
        "u": triu_i.astype(np.int32),
        "v": triu_j.astype(np.int32),
        "gene_u": [gene_names[i] for i in triu_i],
        "gene_v": [gene_names[j] for j in triu_j],
        "w": w,
        "abs_w": np.abs(w),
    }).sort_values("abs_w", ascending=False).reset_index(drop=True)
    return df


if __name__ == "__main__":
    # Quick test
    np.random.seed(42)
    p = 100
    n_edges = p * (p - 1) // 2
    weights = np.random.randn(n_edges).astype(np.float32)
    triu_i, triu_j = np.triu_indices(p, k=1)
    
    g = threshold_edge_vector(weights, triu_i, triu_j, p, mode="topk", topk=10)
    print(f"Graph: {len(g.gene_names)} nodes, {len(g.edges)} edges")
    print(f"Adjacency shape: {g.W.shape}")
    print("graph_builder module loaded successfully")
