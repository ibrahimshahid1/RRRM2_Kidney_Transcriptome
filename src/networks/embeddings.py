"""
node2vec Embeddings with Procrustes Alignment

Implements multi-seed node2vec embedding and orthogonal Procrustes alignment
using pre-registered stable anchor genes.

Key Functions:
    - run_node2vec_multi_seed: Run node2vec with multiple random seeds
    - procrustes_align: Align embedding spaces using anchors
    - consensus_embedding: Average across seeds
"""

import numpy as np
import networkx as nx
from typing import List, Tuple, Dict, Optional
import warnings

# Try importing node2vec libraries
try:
    from pecanpy import node2vec as pecan_node2vec
    PECANPY_AVAILABLE = True
except ImportError:
    PECANPY_AVAILABLE = False
    warnings.warn("PecanPy not available, will use fallback")

try:
    from node2vec import Node2Vec
    NODE2VEC_AVAILABLE = True
except ImportError:
    NODE2VEC_AVAILABLE = False


def weights_to_graph(
    n_genes: int,
    edges: List[Tuple[int, int]],
    weights: np.ndarray,
    mode: str = "pos"
) -> nx.Graph:
    """
    Convert edge list + weights to NetworkX graph.
    
    Parameters
    ----------
    n_genes : int
        Number of nodes (genes)
    edges : list of (int, int)
        Edge list
    weights : np.ndarray
        Edge weights (same length as edges)
    mode : str
        - "pos": keep only positive weights
        - "neg": keep only negative weights as positive magnitude
        - "abs": use abs(weights)
        
    Returns
    -------
    G : nx.Graph
        Weighted graph
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
        else:  # abs
            G.add_edge(i, j, weight=float(abs(w)))
    
    return G


def run_node2vec_single(
    G: nx.Graph,
    dimensions: int = 128,
    walk_length: int = 80,
    num_walks: int = 200,
    p: float = 0.25,
    q: float = 4.0,
    seed: int = 42
) -> np.ndarray:
    """
    Run node2vec on a single graph with one random seed.
    
    Parameters
    ----------
    G : nx.Graph
        Input graph
    dimensions : int
        Embedding dimension
    walk_length : int
        Length of random walks
    num_walks : int
        Number of walks per node
    p : float
        Return parameter
    q : float
        In-out parameter
    seed : int
        Random seed
        
    Returns
    -------
    embeddings : np.ndarray, shape (n_nodes, dimensions)
        Node embeddings
    """
    if PECANPY_AVAILABLE:
        # Use PecanPy (faster)
        # Note: PecanPy requires specific graph format
        # This is a placeholder - adapt to actual PecanPy API
        raise NotImplementedError("PecanPy integration pending")
    
    elif NODE2VEC_AVAILABLE:
        # Use node2vec package
        node2vec = Node2Vec(
            G,
            dimensions=dimensions,
            walk_length=walk_length,
            num_walks=num_walks,
            p=p,
            q=q,
            workers=1
        )
        
        model = node2vec.fit(window=10, min_count=1, batch_words=4, seed=seed)
        
        # Extract embeddings in node order
        n_nodes = G.number_of_nodes()
        embeddings = np.zeros((n_nodes, dimensions), dtype=np.float32)
        for node_id in range(n_nodes):
            if str(node_id) in model.wv:
                embeddings[node_id, :] = model.wv[str(node_id)]
        
        return embeddings
    
    else:
        raise ImportError("No node2vec library available. Install pecanpy or node2vec.")


def run_node2vec_multi_seed(
    G: nx.Graph,
    n_seeds: int = 10,
    **node2vec_params
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run node2vec with multiple random seeds.
    
    Parameters
    ----------
    G : nx.Graph
        Input graph
    n_seeds : int
        Number of random seeds
    **node2vec_params
        Parameters for node2vec (dimensions, walk_length, etc.)
        
    Returns
    -------
    embeddings_list : np.ndarray, shape (n_seeds, n_nodes, dimensions)
        Embeddings from each seed
    variance : np.ndarray, shape (n_nodes,)
        Variance across seeds (for stability assessment)
    """
    n_nodes = G.number_of_nodes()
    dimensions = node2vec_params.get('dimensions', 128)
    
    embeddings_list = np.zeros((n_seeds, n_nodes, dimensions), dtype=np.float32)
    
    for s in range(n_seeds):
        print(f"Running node2vec seed {s+1}/{n_seeds}...")
        seed = 42 + s
        emb = run_node2vec_single(G, seed=seed, **node2vec_params)
        embeddings_list[s, :, :] = emb
    
    # Compute variance across seeds
    variance = embeddings_list.var(axis=0).mean(axis=1)
    
    return embeddings_list, variance


def procrustes_align(
    embeddings: np.ndarray,
    reference: np.ndarray,
    anchor_indices: List[int]
) -> Tuple[np.ndarray, float]:
    """
    Orthogonal Procrustes alignment using anchor genes.
    
    Parameters
    ----------
    embeddings : np.ndarray, shape (n_nodes, dimensions)
        Embedding to align (target)
    reference : np.ndarray, shape (n_nodes, dimensions)
        Reference embedding
    anchor_indices : list of int
        Indices of anchor genes for alignment
        
    Returns
    -------
    aligned : np.ndarray
        Aligned embedding
    alignment_error : float
        Frobenius norm of alignment residual
        
    Notes
    -----
    Finds rotation matrix R that minimizes ||Y - X @ R||_F for anchors.
    Then applies R to all nodes.
    """
    from scipy.linalg import orthogonal_procrustes
    
    # Extract anchor points
    X_anchor = embeddings[anchor_indices, :]
    Y_anchor = reference[anchor_indices, :]
    
    # Find rotation
    R, scale = orthogonal_procrustes(X_anchor, Y_anchor)
    
    # Apply to all nodes
    aligned = embeddings @ R
    
    # Compute alignment error on anchors
    alignment_error = np.linalg.norm(Y_anchor - X_anchor @ R, 'fro')
    
    return aligned.astype(np.float32), alignment_error


def consensus_embedding(
    embeddings_list: np.ndarray,
    reference_idx: int = 0,
    anchor_indices: Optional[List[int]] = None
) -> np.ndarray:
    """
    Average embeddings across seeds after alignment.
    
    Parameters
    ----------
    embeddings_list : np.ndarray, shape (n_seeds, n_nodes, dimensions)
        Embeddings from multiple seeds
    reference_idx : int
        Which seed to use as reference for alignment
    anchor_indices : list of int, optional
        Anchor genes for alignment. If None, align using all nodes.
        
    Returns
    -------
    consensus : np.ndarray, shape (n_nodes, dimensions)
        Average embedding across seeds
    """
    n_seeds, n_nodes, dimensions = embeddings_list.shape
    
    reference = embeddings_list[reference_idx, :, :]
    
    if anchor_indices is None:
        anchor_indices = list(range(n_nodes))
    
    aligned_list = np.zeros_like(embeddings_list)
    aligned_list[reference_idx, :, :] = reference
    
    for s in range(n_seeds):
        if s == reference_idx:
            continue
        aligned, err = procrustes_align(
            embeddings_list[s, :, :],
            reference,
            anchor_indices
        )
        aligned_list[s, :, :] = aligned
        print(f"Seed {s} alignment error: {err:.4f}")
    
    consensus = aligned_list.mean(axis=0)
    
    return consensus.astype(np.float32)


if __name__ == "__main__":
    print("Embeddings module loaded successfully")
    print("Key functions: run_node2vec_multi_seed, procrustes_align, consensus_embedding")
    print(f"PecanPy available: {PECANPY_AVAILABLE}")
    print(f"node2vec available: {NODE2VEC_AVAILABLE}")
