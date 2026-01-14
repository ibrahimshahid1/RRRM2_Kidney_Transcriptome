# src/statistics/rewiring_metrics.py
"""
Rewiring Metrics and Silent Shifter Detection

Quantifies network rewiring from LIONESS edge differences and embedding shifts.
Includes permutation testing for significance.

Key Functions:
    - node_rewiring_from_delta_edges: Per-gene rewiring from edge differences
    - top_differential_edges: Top edges by delta
    - delta_edges_from_samples: Compute deltas from LIONESS matrix
    - permutation_test_node_rewiring_abs: Permutation test for significance
    - embedding_shift_table: Per-gene embedding shifts
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class EdgeIndex:
    """Index structure for edge vectors."""
    gene_names: list
    triu_i: np.ndarray
    triu_j: np.ndarray


# ============================================================================
# EDGE-BASED REWIRING (for LIONESS networks)
# ============================================================================

def node_rewiring_from_delta_edges(
    delta_edges: np.ndarray,
    edge_index: EdgeIndex,
) -> pd.DataFrame:
    """
    Compute per-gene rewiring scores from edge-level differences.
    
    Given delta_edges = mean(groupA edges) - mean(groupB edges), compute
    how much each node's connectivity changed.
    
    Parameters
    ----------
    delta_edges : np.ndarray, shape (n_edges,)
        Difference in edge weights between two groups
    edge_index : EdgeIndex
        Edge index structure with triu_i, triu_j, gene_names
        
    Returns
    -------
    pd.DataFrame
        Columns: gene, rewiring_abs, rewiring_signed
        Sorted by rewiring_abs descending
    """
    delta_edges = np.asarray(delta_edges, dtype=np.float64)
    triu_i = edge_index.triu_i
    triu_j = edge_index.triu_j
    genes = edge_index.gene_names
    p = len(genes)

    abs_delta = np.abs(delta_edges)

    # Absolute rewiring: sum of |Δw| for all edges incident to node
    abs_score = (
        np.bincount(triu_i, weights=abs_delta, minlength=p) +
        np.bincount(triu_j, weights=abs_delta, minlength=p)
    )

    # Signed rewiring: sum of Δw (directional, can cancel)
    signed_score = (
        np.bincount(triu_i, weights=delta_edges, minlength=p) +
        np.bincount(triu_j, weights=delta_edges, minlength=p)
    )

    df = pd.DataFrame({
        "gene": genes,
        "rewiring_abs": abs_score,
        "rewiring_signed": signed_score,
    })
    return df.sort_values("rewiring_abs", ascending=False).reset_index(drop=True)


def top_differential_edges(
    delta_edges: np.ndarray,
    edge_index: EdgeIndex,
    top_n: int = 200,
    sort_by: Literal["abs", "signed"] = "abs",
) -> pd.DataFrame:
    """
    Return top differential edges from delta vector.
    
    Parameters
    ----------
    delta_edges : np.ndarray
        Edge-level differences
    edge_index : EdgeIndex
        Edge index structure
    top_n : int
        Number of top edges to return
    sort_by : str
        "abs": top by |delta|
        "signed": top by delta (positive first)
        
    Returns
    -------
    pd.DataFrame
        Columns: gene_u, gene_v, delta, abs_delta, u, v
    """
    delta_edges = np.asarray(delta_edges, dtype=np.float64)
    triu_i = edge_index.triu_i
    triu_j = edge_index.triu_j
    genes = edge_index.gene_names

    absd = np.abs(delta_edges)
    if top_n <= 0:
        return pd.DataFrame(columns=["gene_u", "gene_v", "delta", "abs_delta", "u", "v"])

    if sort_by == "abs":
        idx = np.argsort(absd)[::-1][:top_n]
    elif sort_by == "signed":
        idx = np.argsort(delta_edges)[::-1][:top_n]
    else:
        raise ValueError("sort_by must be 'abs' or 'signed'")

    uu = triu_i[idx].astype(np.int32)
    vv = triu_j[idx].astype(np.int32)

    df = pd.DataFrame({
        "u": uu,
        "v": vv,
        "gene_u": [genes[i] for i in uu],
        "gene_v": [genes[j] for j in vv],
        "delta": delta_edges[idx].astype(np.float32),
        "abs_delta": absd[idx].astype(np.float32),
    })
    if sort_by == "abs":
        df = df.sort_values("abs_delta", ascending=False)
    else:
        df = df.sort_values("delta", ascending=False)
    return df.reset_index(drop=True)


def delta_edges_from_samples(
    lioness_edges: np.ndarray,
    idx_A: np.ndarray,
    idx_B: np.ndarray,
) -> np.ndarray:
    """
    Compute delta_edges = mean(A) - mean(B) from LIONESS edge matrix.
    
    Parameters
    ----------
    lioness_edges : np.ndarray, shape (n_samples, n_edges)
        LIONESS edge matrix
    idx_A, idx_B : np.ndarray
        Sample indices for groups A and B
        
    Returns
    -------
    np.ndarray
        Delta edges, shape (n_edges,)
    """
    A = lioness_edges[idx_A].mean(axis=0)
    B = lioness_edges[idx_B].mean(axis=0)
    return (A - B).astype(np.float32)


def embedding_shift_table(
    gene_names: list,
    shift_values: np.ndarray,
    label: str = "shift_l2",
) -> pd.DataFrame:
    """
    Create a shift table from pre-computed shift values.
    
    Parameters
    ----------
    gene_names : list
        Gene names
    shift_values : np.ndarray
        Per-gene shift values
    label : str
        Column name for the shift values
        
    Returns
    -------
    pd.DataFrame
        Sorted by shift descending
    """
    df = pd.DataFrame({"gene": gene_names, label: np.asarray(shift_values, dtype=np.float64)})
    return df.sort_values(label, ascending=False).reset_index(drop=True)


def permutation_test_node_rewiring_abs(
    lioness_edges: np.ndarray,
    edge_index: EdgeIndex,
    idx_A: np.ndarray,
    idx_B: np.ndarray,
    n_perm: int = 200,
    seed: int = 0,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Permutation test for node rewiring_abs significance.
    
    This is simple (and computationally intensive for large n_perm),
    but works well for p~1200 and n_perm~200.
    
    Parameters
    ----------
    lioness_edges : np.ndarray, shape (n_samples, n_edges)
        LIONESS edge matrix
    edge_index : EdgeIndex
        Edge index structure
    idx_A, idx_B : np.ndarray
        Sample indices for groups A and B
    n_perm : int
        Number of permutations
    seed : int
        Random seed
    verbose : bool
        Print progress
        
    Returns
    -------
    pd.DataFrame
        Columns: gene, rewiring_abs_observed, p_value
    """
    rng = np.random.default_rng(seed)
    n = lioness_edges.shape[0]
    all_idx = np.arange(n)

    idx_A = np.asarray(idx_A, dtype=int)
    idx_B = np.asarray(idx_B, dtype=int)
    nA = len(idx_A)
    nB = len(idx_B)

    # Observed statistic
    obs_delta = delta_edges_from_samples(lioness_edges, idx_A, idx_B).astype(np.float64)
    obs = node_rewiring_from_delta_edges(obs_delta, edge_index)["rewiring_abs"].values

    perm_ge = np.zeros_like(obs, dtype=np.int32)

    for i in range(n_perm):
        if verbose and i % 50 == 0:
            print(f"  Permutation {i+1}/{n_perm}", end="\r")
        
        perm = rng.permutation(all_idx)
        pA = perm[:nA]
        pB = perm[nA:nA+nB]
        d = delta_edges_from_samples(lioness_edges, pA, pB).astype(np.float64)
        s = node_rewiring_from_delta_edges(d, edge_index)["rewiring_abs"].values
        perm_ge += (s >= obs).astype(np.int32)

    if verbose:
        print(f"  Completed {n_perm} permutations" + " " * 20)

    pvals = (perm_ge + 1) / (n_perm + 1)

    out = pd.DataFrame({
        "gene": edge_index.gene_names,
        "rewiring_abs_observed": obs,
        "p_value": pvals,
    }).sort_values("rewiring_abs_observed", ascending=False).reset_index(drop=True)

    return out


# ============================================================================
# LEGACY COMPATIBILITY - keep old function signatures working
# ============================================================================

def node_rewiring_from_delta_edges_legacy(
    delta_edges: np.ndarray,
    triu_i: np.ndarray,
    triu_j: np.ndarray,
    gene_names: list,
) -> pd.DataFrame:
    """
    Legacy interface - wraps new EdgeIndex-based function.
    """
    edge_index = EdgeIndex(gene_names=gene_names, triu_i=triu_i, triu_j=triu_j)
    return node_rewiring_from_delta_edges(delta_edges, edge_index)


if __name__ == "__main__":
    # Quick test
    np.random.seed(42)
    p = 100
    n_edges = p * (p - 1) // 2
    
    delta = np.random.randn(n_edges).astype(np.float32)
    triu_i, triu_j = np.triu_indices(p, k=1)
    genes = [f"gene_{i}" for i in range(p)]
    
    edge_index = EdgeIndex(gene_names=genes, triu_i=triu_i, triu_j=triu_j)
    
    df = node_rewiring_from_delta_edges(delta, edge_index)
    print(f"Rewiring table: {len(df)} genes")
    print(df.head())
    
    top = top_differential_edges(delta, edge_index, top_n=10)
    print(f"\nTop differential edges: {len(top)} edges")
    print(top)
    
    print("\nrewiring_metrics module loaded successfully")
