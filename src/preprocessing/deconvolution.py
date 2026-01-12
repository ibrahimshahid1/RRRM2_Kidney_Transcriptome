"""
Cell-Type Deconvolution for Bulk RNA-Seq

Estimates nephron segment proportions per sample using murine kidney
single-cell reference atlases. Handles compositional transformation (CLR)
for downstream regression.

Key Functions:
    - load_reference_atlas: Load single-cell kidney reference
    - deconvolve_samples: Estimate cell-type proportions
    - clr_transform: Centered log-ratio transformation
"""

import numpy as np
import pandas as pd
from typing import Tuple, Optional, Dict
from scipy.stats import pearsonr
from sklearn.linear_model import Ridge


def load_reference_atlas(
    atlas_path: Optional[str] = None,
    atlas_source: str = "tabula_muris"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load murine kidney single-cell reference atlas.
    
    Parameters
    ----------
    atlas_path : str, optional
        Path to pre-processed atlas file
    atlas_source : str
        Reference atlas to use: 'tabula_muris', 'park2018', etc.
        
    Returns
    -------
    ref_expr : pd.DataFrame
        Gene expression profiles (genes × cell types)
    cell_metadata : pd.DataFrame
        Cell type annotations
        
    Notes
    -----
    User must provide actual reference data. This is a placeholder
    that defines expected structure.
    """
    if atlas_path is not None:
        # Load user-provided atlas
        # Expected format: genes as rows, cell types as columns
        ref_expr = pd.read_csv(atlas_path, index_col=0)
        cell_metadata = pd.DataFrame({
            'cell_type': ref_expr.columns,
            'segment': ref_expr.columns  # Map to nephron segments
        })
    else:
        # Placeholder - user must implement
        raise NotImplementedError(
            f"Please provide atlas_path or implement {atlas_source} loader"
        )
    
    return ref_expr, cell_metadata


def deconvolve_samples(
    bulk_expr: pd.DataFrame,
    ref_expr: pd.DataFrame,
    method: str = "ridge",
    alpha: float = 1.0
) -> pd.DataFrame:
    """
    Estimate cell-type proportions from bulk RNA-seq.
    
    Parameters
    ----------
    bulk_expr : pd.DataFrame, shape (n_samples, n_genes)
        Bulk expression matrix (VST-transformed)
    ref_expr : pd.DataFrame, shape (n_genes, n_cell_types)
        Reference expression profiles per cell type
    method : str
        Deconvolution method: 'ridge', 'nnls', 'cor'
    alpha : float
        Regularization parameter for ridge regression
        
    Returns
    -------
    proportions : pd.DataFrame, shape (n_samples, n_cell_types)
        Estimated cell-type proportions (rows sum to 1)
        
    Notes
    -----
    For production use, consider dedicated tools:
    - CIBERSORTx (requires web service)
    - MuSiC (R package)
    - Scaden (deep learning)
    This implements a simple regression-based approach.
    """
    # Align genes
    common_genes = bulk_expr.columns.intersection(ref_expr.index)
    bulk_aligned = bulk_expr[common_genes].values
    ref_aligned = ref_expr.loc[common_genes].values
    
    n_samples = bulk_aligned.shape[0]
    n_cell_types = ref_aligned.shape[1]
    
    proportions = np.zeros((n_samples, n_cell_types))
    
    if method == "ridge":
        # Ridge regression with non-negativity constraint approximation
        ridge = Ridge(alpha=alpha, fit_intercept=False)
        
        for i in range(n_samples):
            ridge.fit(ref_aligned, bulk_aligned[i, :])
            props = ridge.coef_
            props = np.maximum(props, 0)  # Enforce non-negativity
            props = props / props.sum()    # Normalize to sum to 1
            proportions[i, :] = props
            
    elif method == "nnls":
        from scipy.optimize import nnls
        for i in range(n_samples):
            props, _ = nnls(ref_aligned, bulk_aligned[i, :])
            props = props / props.sum()
            proportions[i, :] = props
            
    elif method == "cor":
        # Correlation-based (simple, fast, but less accurate)
        for i in range(n_samples):
            cors = np.array([
                pearsonr(bulk_aligned[i, :], ref_aligned[:, j])[0]
                for j in range(n_cell_types)
            ])
            cors = np.maximum(cors, 0)
            proportions[i, :] = cors / cors.sum()
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    prop_df = pd.DataFrame(
        proportions,
        index=bulk_expr.index,
        columns=ref_expr.columns
    )
    
    return prop_df


def clr_transform(proportions: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    """
    Centered log-ratio (CLR) transformation for compositional data.
    
    Parameters
    ----------
    proportions : pd.DataFrame
        Cell-type proportions (rows sum to 1)
    pseudocount : float
        Small constant to avoid log(0)
        
    Returns
    -------
    clr_props : pd.DataFrame
        CLR-transformed proportions
        
    Notes
    -----
    CLR avoids perfect collinearity in regression. Alternative: drop one category.
    """
    props = proportions.values + pseudocount
    
    # Geometric mean per sample
    geom_mean = np.exp(np.log(props).mean(axis=1, keepdims=True))
    
    # CLR
    clr = np.log(props / geom_mean)
    
    clr_df = pd.DataFrame(
        clr,
        index=proportions.index,
        columns=[f'CLR_{col}' for col in proportions.columns]
    )
    
    return clr_df


def aggregate_to_segments(
    proportions: pd.DataFrame,
    cell_to_segment_map: Dict[str, str]
) -> pd.DataFrame:
    """
    Aggregate cell types to nephron segments.
    
    Parameters
    ----------
    proportions : pd.DataFrame
        Cell-type proportions
    cell_to_segment_map : dict
        {cell_type: segment} mapping
        
    Returns
    -------
    segment_props : pd.DataFrame
        Aggregated segment proportions
        
    Example
    -------
    >>> mapping = {
    ...     'DCT': 'DCT',
    ...     'PT_S1': 'PT', 'PT_S2': 'PT', 'PT_S3': 'PT',
    ...     'Podocyte': 'Glom', 'Endothelial': 'Glom',
    ...     'CD_PC': 'CD', 'CD_IC': 'CD'
    ... }
    """
    segment_props = pd.DataFrame(index=proportions.index)
    
    segments = set(cell_to_segment_map.values())
    for segment in segments:
        cell_types = [ct for ct, seg in cell_to_segment_map.items() if seg == segment]
        cols_present = [ct for ct in cell_types if ct in proportions.columns]
        segment_props[segment] = proportions[cols_present].sum(axis=1)
    
    return segment_props


if __name__ == "__main__":
    # Example usage
    print("Deconvolution module loaded successfully")
    print("Key functions: load_reference_atlas, deconvolve_samples, clr_transform")
