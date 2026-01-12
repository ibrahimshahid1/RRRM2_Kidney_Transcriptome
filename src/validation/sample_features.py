"""
Sample-Level Feature Extraction from LIONESS Networks

Extracts mouse-level features from sample-specific networks for classification
and validation. Features are biologically interpretable and pathway-focused.

Key Functions:
    - pathway_edge_strength: Mean/median edge weights in pathway subnetworks
    - node_strength_features: Sum of incident weights for candidate genes
    - shifter_connectivity: Connectivity scores for silent shifters
    - pca_features: Low-dimensional network summaries
"""

import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Set
from sklearn.decomposition import PCA


def pathway_edge_strength(
    W_samp: np.ndarray,
    edges: List[Tuple[int, int]],
    gene_names: List[str],
    pathway_genes: Set[str],
    aggregation: str = 'mean'
) -> np.ndarray:
    """
    Compute mean/median edge weight within pathway subnetwork.
    
    Parameters
    ----------
    W_samp : np.ndarray, shape (n_samples, n_edges)
        Sample-specific edge weights
    edges : list of (int, int)
        Edge list E
    gene_names : list of str
        Gene identifiers
    pathway_genes : set of str
        Pathway gene set (e.g., NCC-WNK genes)
    aggregation : str
        'mean' or 'median'
        
    Returns
    -------
    features : np.ndarray, shape (n_samples,)
        Pathway edge strength per sample
    """
    # Find pathway gene indices
    pathway_idx = {i for i, g in enumerate(gene_names) if g in pathway_genes}
    
    # Find edges within pathway
    pathway_edge_mask = np.array([
        (i in pathway_idx and j in pathway_idx)
        for i, j in edges
    ])
    
    if pathway_edge_mask.sum() == 0:
        # No pathway edges
        return np.zeros(W_samp.shape[0])
    
    W_pathway = W_samp[:, pathway_edge_mask]
    
    if aggregation == 'mean':
        features = W_pathway.mean(axis=1)
    elif aggregation == 'median':
        features = np.median(W_pathway, axis=1)
    else:
        raise ValueError(f"Unknown aggregation: {aggregation}")
    
    return features


def node_strength_features(
    W_samp: np.ndarray,
    edges: List[Tuple[int, int]],
    gene_names: List[str],
    candidate_genes: List[str]
) -> np.ndarray:
    """
    Node strength (sum of incident weights) for candidate genes.
    
    Parameters
    ----------
    W_samp : np.ndarray, shape (n_samples, n_edges)
        Sample-specific edge weights
    edges : list of (int, int)
        Edge list
    gene_names : list of str
        Gene identifiers
    candidate_genes : list of str
        Genes to compute strength for
        
    Returns
    -------
    features : np.ndarray, shape (n_samples, len(candidate_genes))
        Node strength per candidate gene per sample
    """
    n_samples = W_samp.shape[0]
    n_candidates = len(candidate_genes)
    
    features = np.zeros((n_samples, n_candidates))
    
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    
    for c, gene in enumerate(candidate_genes):
        if gene not in gene_to_idx:
            continue
        
        node_idx = gene_to_idx[gene]
        
        # Find edges incident to this node
        incident_mask = np.array([
            (i == node_idx or j == node_idx)
            for i, j in edges
        ])
        
        if incident_mask.sum() > 0:
            # Sum absolute weights (undirected strength)
            features[:, c] = np.abs(W_samp[:, incident_mask]).sum(axis=1)
    
    return features


def shifter_connectivity(
    W_samp: np.ndarray,
    edges: List[Tuple[int, int]],
    gene_names: List[str],
    shifter_genes: List[str],
    top_k_neighbors: int = 10
) -> np.ndarray:
    """
    For each silent shifter, aggregate weights to its top-k neighbors in E.
    
    Parameters
    ----------
    W_samp : np.ndarray
        Sample-specific edge weights
    edges : list of (int, int)
        Edge list
    gene_names : list of str
        Gene identifiers
    shifter_genes : list of str
        Silent shifter genes
    top_k_neighbors : int
        Number of top neighbors to consider
        
    Returns
    -------
    features : np.ndarray, shape (n_samples, len(shifter_genes))
        Shifter connectivity scores
    """
    n_samples = W_samp.shape[0]
    n_shifters = len(shifter_genes)
    
    features = np.zeros((n_samples, n_shifters))
    
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    
    for s, gene in enumerate(shifter_genes):
        if gene not in gene_to_idx:
            continue
        
        node_idx = gene_to_idx[gene]
        
        # Find all neighbors
        neighbor_edges = []
        for e_idx, (i, j) in enumerate(edges):
            if i == node_idx:
                neighbor_edges.append((e_idx, j))
            elif j == node_idx:
                neighbor_edges.append((e_idx, i))
        
        if len(neighbor_edges) == 0:
            continue
        
        # For each sample, get top-k neighbor weights
        for samp in range(n_samples):
            edge_weights = [(W_samp[samp, e_idx], nbr) for e_idx, nbr in neighbor_edges]
            edge_weights.sort(reverse=True, key=lambda x: abs(x[0]))
            
            top_weights = [w for w, _ in edge_weights[:top_k_neighbors]]
            features[samp, s] = np.mean(top_weights)
    
    return features


def pca_features(
    W_samp: np.ndarray,
    n_components: int = 10,
    stable_edge_mask: np.ndarray = None
) -> np.ndarray:
    """
    Low-dimensional PCA features from edge weight matrix.
    
    Parameters
    ----------
    W_samp : np.ndarray, shape (n_samples, n_edges)
        Sample-specific edge weights
    n_components : int
        Number of PCA components
    stable_edge_mask : np.ndarray, optional
        Boolean mask for stable edges (to reduce noise)
        
    Returns
    -------
    features : np.ndarray, shape (n_samples, n_components)
        PCA features
    """
    if stable_edge_mask is not None:
        W_stable = W_samp[:, stable_edge_mask]
    else:
        W_stable = W_samp
    
    pca = PCA(n_components=n_components)
    features = pca.fit_transform(W_stable)
    
    return features.astype(np.float32)


def combine_features(
    feature_dict: Dict[str, np.ndarray]
) -> Tuple[np.ndarray, List[str]]:
    """
    Combine multiple feature types into one matrix.
    
    Parameters
    ----------
    feature_dict : dict
        {feature_name: feature_array}
        Each array should be (n_samples, n_features_i)
        
    Returns
    -------
    X : np.ndarray, shape (n_samples, total_features)
        Combined feature matrix
    feature_names : list of str
        Feature labels
    """
    feature_names = []
    feature_list = []
    
    for name, arr in feature_dict.items():
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
            cols = [name]
        else:
            cols = [f"{name}_{i}" for i in range(arr.shape[1])]
        
        feature_list.append(arr)
        feature_names.extend(cols)
    
    X = np.hstack(feature_list)
    
    return X, feature_names


if __name__ == "__main__":
    print("Sample features module loaded successfully")
    print("Key functions: pathway_edge_strength, node_strength_features, pca_features")
