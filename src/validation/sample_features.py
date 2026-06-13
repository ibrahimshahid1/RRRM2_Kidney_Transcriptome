"""Sample-Level Feature Extraction for RRRM-2"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass


@dataclass
class SampleFeatures:
    """Container for sample-level network features."""
    sample_ids: List[str]
    feature_matrix: np.ndarray
    feature_names: List[str]
    
    def to_dataframe(self) -> pd.DataFrame:
        """Convert to DataFrame."""
        return pd.DataFrame(
            self.feature_matrix,
            index=self.sample_ids,
            columns=self.feature_names
        )


def node_strength(
    lioness_z: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    n_genes: int,
    genes: List[str] = None
) -> np.ndarray:
    """
    Compute node strength (sum of incident edge weights) per sample.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        edge_i: Source gene indices for each edge
        edge_j: Target gene indices for each edge
        n_genes: Total number of genes
        genes: Optional gene list for column names
        
    Returns:
        Node strength matrix (samples x genes)
    """
    n_samples = lioness_z.shape[0]
    strength = np.zeros((n_samples, n_genes), dtype=np.float32)
    
    for s in range(n_samples):
        w = lioness_z[s]
        np.add.at(strength[s], edge_i, w)
        np.add.at(strength[s], edge_j, w)
    
    return strength


def pathway_edge_summary(
    lioness_z: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    genes: List[str],
    pathway_genes: List[str],
    agg_func: str = 'mean'
) -> np.ndarray:
    """
    Compute summary statistics for edges within a pathway.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        edge_i: Source gene indices
        edge_j: Target gene indices
        genes: Full gene list
        pathway_genes: Genes in the pathway of interest
        agg_func: Aggregation function ('mean', 'median', 'sum', 'std')
        
    Returns:
        Pathway summary per sample (n_samples,)
    """
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    pathway_idx = set(gene_to_idx.get(g, -1) for g in pathway_genes)
    pathway_idx.discard(-1)
    
    # Find edges where both genes are in pathway
    within_mask = np.array([
        (i in pathway_idx) and (j in pathway_idx)
        for i, j in zip(edge_i, edge_j)
    ])
    
    if not within_mask.any():
        return np.zeros(lioness_z.shape[0], dtype=np.float32)
    
    pathway_weights = lioness_z[:, within_mask]
    
    if agg_func == 'mean':
        return pathway_weights.mean(axis=1)
    elif agg_func == 'median':
        return np.median(pathway_weights, axis=1)
    elif agg_func == 'sum':
        return pathway_weights.sum(axis=1)
    elif agg_func == 'std':
        return pathway_weights.std(axis=1)
    else:
        raise ValueError(f"Unknown agg_func: {agg_func}")


def shifter_connectivity(
    lioness_z: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    genes: List[str],
    shifter_genes: List[str],
    topk_neighbors: int = 10
) -> np.ndarray:
    """
    Compute shifter-centered connectivity scores.
    
    For each silent shifter, aggregate weights to its top neighbors.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        edge_i: Source gene indices
        edge_j: Target gene indices
        genes: Full gene list
        shifter_genes: Silent shifter genes
        topk_neighbors: Number of top neighbors to consider
        
    Returns:
        Shifter connectivity features (samples x len(shifter_genes))
    """
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    n_samples = lioness_z.shape[0]
    
    features = []
    feature_names = []
    
    for shifter in shifter_genes:
        if shifter not in gene_to_idx:
            continue
        
        sidx = gene_to_idx[shifter]
        
        # Find edges incident to this shifter
        incident_mask = (edge_i == sidx) | (edge_j == sidx)
        
        if not incident_mask.any():
            features.append(np.zeros(n_samples, dtype=np.float32))
            feature_names.append(f"{shifter}_connectivity")
            continue
        
        # Get mean weight across incident edges
        incident_weights = lioness_z[:, incident_mask]
        features.append(incident_weights.mean(axis=1))
        feature_names.append(f"{shifter}_connectivity")
    
    if not features:
        return np.zeros((n_samples, 0), dtype=np.float32)
    
    return np.column_stack(features)


def edge_pca_features(
    lioness_z: np.ndarray,
    n_components: int = 10,
    stable_edge_idx: np.ndarray = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract PC features from LIONESS edge weights.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        n_components: Number of PCs to extract
        stable_edge_idx: Optional subset of stable edges to use
        
    Returns:
        Tuple of (PC scores, explained variance ratios)
    """
    from sklearn.decomposition import PCA
    
    if stable_edge_idx is not None:
        X = lioness_z[:, stable_edge_idx]
    else:
        X = lioness_z
    
    n_components = min(n_components, min(X.shape) - 1)
    
    pca = PCA(n_components=n_components)
    pc_scores = pca.fit_transform(X)
    
    return pc_scores, pca.explained_variance_ratio_


def extract_all_features(
    lioness_z: np.ndarray,
    sample_ids: List[str],
    genes: List[str],
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    pathway_dict: Dict[str, List[str]] = None,
    shifter_genes: List[str] = None,
    n_pcs: int = 10
) -> SampleFeatures:
    """
    Extract comprehensive sample-level features.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        sample_ids: Sample identifiers
        genes: Gene list
        edge_i: Source gene indices
        edge_j: Target gene indices
        pathway_dict: Dictionary of pathway name -> gene list
        shifter_genes: Silent shifter genes
        n_pcs: Number of PCs to extract
        
    Returns:
        SampleFeatures container
    """
    all_features = []
    all_names = []
    
    # 1. Node strengths for all genes
    strength = node_strength(lioness_z, edge_i, edge_j, len(genes))
    # Just use top-variance genes as features
    strength_var = strength.var(axis=0)
    top_strength_idx = np.argsort(strength_var)[-50:]  # Top 50 variable
    for idx in top_strength_idx:
        all_features.append(strength[:, idx])
        all_names.append(f"strength_{genes[idx]}")
    
    # 2. Pathway summaries
    if pathway_dict:
        for pathway_name, pathway_genes in pathway_dict.items():
            for agg in ['mean', 'std']:
                feat = pathway_edge_summary(
                    lioness_z, edge_i, edge_j, genes, pathway_genes, agg
                )
                all_features.append(feat)
                all_names.append(f"pathway_{pathway_name}_{agg}")
    
    # 3. Shifter connectivity
    if shifter_genes:
        shifter_feats = shifter_connectivity(
            lioness_z, edge_i, edge_j, genes, shifter_genes
        )
        for i in range(shifter_feats.shape[1]):
            all_features.append(shifter_feats[:, i])
            all_names.append(f"shifter_{shifter_genes[i]}")
    
    # 4. Edge PCA
    pc_scores, var_ratios = edge_pca_features(lioness_z, n_pcs)
    for i in range(pc_scores.shape[1]):
        all_features.append(pc_scores[:, i])
        all_names.append(f"edge_PC{i+1}")
    
    # Stack all features
    feature_matrix = np.column_stack(all_features) if all_features else np.zeros((len(sample_ids), 0))
    
    return SampleFeatures(
        sample_ids=sample_ids,
        feature_matrix=feature_matrix.astype(np.float32),
        feature_names=all_names
    )



