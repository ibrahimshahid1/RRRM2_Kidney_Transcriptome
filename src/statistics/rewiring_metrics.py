"""
Rewiring Metrics and Silent Shifter Detection

Quantifies network rewiring as cosine distance shifts between aligned embeddings.
Defines "silent shifters" as genes with high rewiring but low differential expression.

Key Functions:
    - cosine_distance: Compute cosine distance between vectors
    - compute_rewiring_metrics: Primary factorial design metrics
    - compute_legacy_metrics: 4-bin pooled metrics
    - identify_silent_shifters: High Δ, low DE genes
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from scipy.spatial.distance import cosine


def cosine_distance(v1: np.ndarray, v2: np.ndarray) -> float:
    """
    Compute cosine distance: 1 - cos(v1, v2).
    
    Parameters
    ----------
    v1, v2 : np.ndarray
        Embedding vectors
        
    Returns
    -------
    distance : float
        Cosine distance in [0, 2]
    """
    return cosine(v1, v2)


def compute_rewiring_metrics(
    embeddings: Dict[str, np.ndarray],
    gene_names: List[str]
) -> pd.DataFrame:
    """
    Compute primary rewiring metrics over factorial design.
    
    Parameters
    ----------
    embeddings : dict
        {condition_name: embedding_array}
        Expected keys (ISS-T):
            - 'Young_ISST_FLT', 'Young_ISST_HGC'
            - 'Old_ISST_FLT', 'Old_ISST_HGC'
        And similarly for LAR
    gene_names : list of str
        Gene identifiers
        
    Returns
    -------
    metrics : pd.DataFrame
        Columns: gene, Δ_flt_young_ISST, Δ_flt_old_ISST, Δ_interaction_ISST, ...
        
    Notes
    -----
    Primary metrics (model-based):
        Δ_ISS-T_flt,young = cosine_dist(Young_ISST_FLT, Young_ISST_HGC)
        Δ_ISS-T_flt,old = cosine_dist(Old_ISST_FLT, Old_ISST_HGC)
        Δ_interaction = |Δ_old - Δ_young|
    """
    n_genes = len(gene_names)
    
    metrics = {
        'gene': gene_names
    }
    
    # ISS-T flight effects
    if all(k in embeddings for k in ['Young_ISST_FLT', 'Young_ISST_HGC']):
        delta_young_isst = np.zeros(n_genes)
        for g in range(n_genes):
            v_flt = embeddings['Young_ISST_FLT'][g, :]
            v_hgc = embeddings['Young_ISST_HGC'][g, :]
            delta_young_isst[g] = cosine_distance(v_flt, v_hgc)
        metrics['delta_flt_young_ISST'] = delta_young_isst
    
    if all(k in embeddings for k in ['Old_ISST_FLT', 'Old_ISST_HGC']):
        delta_old_isst = np.zeros(n_genes)
        for g in range(n_genes):
            v_flt = embeddings['Old_ISST_FLT'][g, :]
            v_hgc = embeddings['Old_ISST_HGC'][g, :]
            delta_old_isst[g] = cosine_distance(v_flt, v_hgc)
        metrics['delta_flt_old_ISST'] = delta_old_isst
    
    # Interaction
    if 'delta_flt_young_ISST' in metrics and 'delta_flt_old_ISST' in metrics:
        metrics['delta_interaction_ISST'] = np.abs(
            metrics['delta_flt_old_ISST'] - metrics['delta_flt_young_ISST']
        )
    
    # Repeat for LAR (if available)
    # TODO: Add LAR metrics similarly
    
    return pd.DataFrame(metrics)


def compute_legacy_metrics(
    embeddings: Dict[str, np.ndarray],
    gene_names: List[str]
) -> pd.DataFrame:
    """
    Compute legacy 4-bin pooled metrics for sensitivity/summary.
    
    Parameters
    ----------
    embeddings : dict
        Expected keys: YC, YF, OC, OF (pooled across arms)
    gene_names : list of str
        
    Returns
    -------
    metrics : pd.DataFrame
        Columns: gene, delta_age_ctrl, delta_age_flt, delta_flt_young, delta_flt_old
        
    Notes
    -----
    Legacy metrics (4-bin):
        Δ_age,ctrl = cosine_dist(OC, YC)
        Δ_age,flt = cosine_dist(OF, YF)
        Δ_flt,young = cosine_dist(YF, YC)
        Δ_flt,old = cosine_dist(OF, OC)
    """
    n_genes = len(gene_names)
    
    metrics = {
        'gene': gene_names
    }
    
    # Age effect in controls
    if 'OC' in embeddings and 'YC' in embeddings:
        delta_age_ctrl = np.array([
            cosine_distance(embeddings['OC'][g, :], embeddings['YC'][g, :])
            for g in range(n_genes)
        ])
        metrics['delta_age_ctrl'] = delta_age_ctrl
    
    # Age effect in flight
    if 'OF' in embeddings and 'YF' in embeddings:
        delta_age_flt = np.array([
            cosine_distance(embeddings['OF'][g, :], embeddings['YF'][g, :])
            for g in range(n_genes)
        ])
        metrics['delta_age_flt'] = delta_age_flt
    
    # Flight effect in young
    if 'YF' in embeddings and 'YC' in embeddings:
        delta_flt_young = np.array([
            cosine_distance(embeddings['YF'][g, :], embeddings['YC'][g, :])
            for g in range(n_genes)
        ])
        metrics['delta_flt_young'] = delta_flt_young
    
    # Flight effect in old
    if 'OF' in embeddings and 'OC' in embeddings:
        delta_flt_old = np.array([
            cosine_distance(embeddings['OF'][g, :], embeddings['OC'][g, :])
            for g in range(n_genes)
        ])
        metrics['delta_flt_old'] = delta_flt_old
    
    return pd.DataFrame(metrics)


def identify_silent_shifters(
    rewiring_df: pd.DataFrame,
    de_results: pd.DataFrame,
    delta_col: str = 'delta_flt_young_ISST',
    delta_percentile: float = 90,
    log2fc_threshold: float = 0.3,
    de_fdr_threshold: float = 0.2,
    delta_fdr_col: Optional[str] = None,
    delta_fdr_threshold: float = 0.1
) -> pd.DataFrame:
    """
    Identify silent shifters: high rewiring + low differential expression.
    
    Parameters
    ----------
    rewiring_df : pd.DataFrame
        Output from compute_rewiring_metrics
    de_results : pd.DataFrame
        Differential expression results (must have: gene, log2FC, FDR)
    delta_col : str
        Which Δ metric to use
    delta_percentile : float
        Percentile threshold for high rewiring (e.g., 90 = top 10%)
    log2fc_threshold : float
        Maximum |log2FC| for "low DE"
    de_fdr_threshold : float
        Minimum FDR for "not significantly DE"
    delta_fdr_col : str, optional
        Column name for rewiring FDR (if available)
    delta_fdr_threshold : float
        Maximum FDR for significant rewiring
        
    Returns
    -------
    silent_shifters : pd.DataFrame
        Genes meeting silent shifter criteria
    """
    # Merge rewiring and DE
    merged = rewiring_df.merge(de_results, on='gene', how='inner')
    
    # High rewiring criterion
    delta_threshold = np.percentile(merged[delta_col], delta_percentile)
    high_rewiring = merged[delta_col] >= delta_threshold
    
    # Low DE criterion
    low_de = (np.abs(merged['log2FC']) < log2fc_threshold) & (merged['FDR'] > de_fdr_threshold)
    
    # Significant rewiring (if FDR available)
    if delta_fdr_col is not None and delta_fdr_col in merged.columns:
        sig_rewiring = merged[delta_fdr_col] < delta_fdr_threshold
    else:
        sig_rewiring = high_rewiring  # Use percentile only
    
    # Combine criteria
    silent_shifters = merged[high_rewiring & low_de & sig_rewiring].copy()
    
    # Sort by rewiring magnitude
    silent_shifters = silent_shifters.sort_values(delta_col, ascending=False)
    
    return silent_shifters


if __name__ == "__main__":
    print("Rewiring metrics module loaded successfully")
    print("Key functions: compute_rewiring_metrics, identify_silent_shifters")
