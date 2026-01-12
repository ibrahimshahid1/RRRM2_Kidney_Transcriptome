"""
Global Residualization for Confounder Removal

Implements Phase 1 of the pipeline: global regression to remove technical and
compositional confounding while preserving biological signal for network inference.

Key Functions:
    - residualize_expr: Global regression across all samples
    - fit_sva: Surrogate variable analysis
    - create_design_matrix: Build covariate matrix from metadata
"""

import numpy as np
import pandas as pd
from typing import Tuple, Dict, Optional, List
from sklearn.linear_model import LinearRegression
from sklearn.decomposition import PCA


def residualize_expr(
    Y_samples_by_genes: np.ndarray,
    X_covariates: np.ndarray
) -> np.ndarray:
    """
    Residualize each gene on covariates (global regression).
    
    Parameters
    ----------
    Y : np.ndarray, shape (n_samples, n_genes)
        VST-transformed expression matrix
    X : np.ndarray, shape (n_samples, n_covariates)
        Covariate matrix: batch + SVs + transformed cell proportions
        
    Returns
    -------
    R : np.ndarray, shape (n_samples, n_genes)
        Residual matrix preserving biology not included in X
        
    Notes
    -----
    Do NOT include Age/Arm/Group labels in X if you want biology preserved
    for rewiring comparisons. Use separate R_tech and R_all matrices.
    """
    lr = LinearRegression()
    R = np.empty_like(Y_samples_by_genes, dtype=np.float32)
    
    for g in range(Y_samples_by_genes.shape[1]):
        lr.fit(X_covariates, Y_samples_by_genes[:, g])
        R[:, g] = (Y_samples_by_genes[:, g] - lr.predict(X_covariates)).astype(np.float32)
    
    return R


def fit_sva(
    Y: np.ndarray,
    design_matrix: pd.DataFrame,
    biological_vars: List[str],
    n_sv: Optional[int] = None
) -> Tuple[np.ndarray, int]:
    """
    Surrogate Variable Analysis to capture hidden confounding.
    
    Parameters
    ----------
    Y : np.ndarray
        Expression matrix (samples × genes)
    design_matrix : pd.DataFrame
        Known covariates (batch, etc.)
    biological_vars : list of str
        Column names of biological variables to preserve
    n_sv : int, optional
        Number of surrogate variables. If None, estimated automatically.
        
    Returns
    -------
    sv_matrix : np.ndarray, shape (n_samples, n_sv)
        Surrogate variables
    n_sv : int
        Number of SVs computed
        
    Notes
    -----
    This is a simplified PCA-based implementation. For production use,
    consider sva R package via rpy2 or python-sva.
    """
    # Remove biological signal first
    bio_cols = [col for col in design_matrix.columns if col in biological_vars]
    X_bio = design_matrix[bio_cols].values
    
    lr = LinearRegression()
    residuals = np.empty_like(Y)
    for g in range(Y.shape[1]):
        lr.fit(X_bio, Y[:, g])
        residuals[:, g] = Y[:, g] - lr.predict(X_bio)
    
    # PCA on residuals to find hidden structure
    if n_sv is None:
        # Estimate via scree plot elbow or percentage variance
        pca_full = PCA()
        pca_full.fit(residuals)
        cumvar = np.cumsum(pca_full.explained_variance_ratio_)
        n_sv = int(np.where(cumvar > 0.5)[0][0] + 1)  # 50% variance threshold
        n_sv = min(n_sv, 10)  # Cap at 10
    
    pca = PCA(n_components=n_sv)
    sv_matrix = pca.fit_transform(residuals)
    
    return sv_matrix.astype(np.float32), n_sv


def create_design_matrix(
    metadata: pd.DataFrame,
    cell_props: pd.DataFrame,
    sv_matrix: Optional[np.ndarray] = None,
    include_bio: bool = False
) -> pd.DataFrame:
    """
    Create covariate design matrix for residualization.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Sample metadata with batch variables
    cell_props : pd.DataFrame
        Cell-type proportions (CLR-transformed)
    sv_matrix : np.ndarray, optional
        Surrogate variables from SVA
    include_bio : bool, default=False
        If True, include Age/Arm/Group for R_all matrix
        If False, exclude for R_tech matrix (biology preserved)
        
    Returns
    -------
    design : pd.DataFrame
        Design matrix for regression
    """
    design_parts = []
    
    # Technical covariates
    tech_cols = ['batch', 'lane', 'shipment', 'RIN']
    tech_cols = [c for c in tech_cols if c in metadata.columns]
    if tech_cols:
        design_parts.append(metadata[tech_cols])
    
    # Cell-type proportions (CLR-transformed)
    design_parts.append(cell_props)
    
    # Surrogate variables
    if sv_matrix is not None:
        sv_df = pd.DataFrame(
            sv_matrix,
            index=metadata.index,
            columns=[f'SV{i+1}' for i in range(sv_matrix.shape[1])]
        )
        design_parts.append(sv_df)
    
    # Biological variables (optional)
    if include_bio:
        bio_cols = ['Age', 'Arm', 'EnvironmentGroup']
        bio_cols = [c for c in bio_cols if c in metadata.columns]
        if bio_cols:
            design_parts.append(pd.get_dummies(metadata[bio_cols], drop_first=True))
    
    design = pd.concat(design_parts, axis=1)
    return design


def variance_partition(
    Y: np.ndarray,
    design_matrix: pd.DataFrame,
    covariate_groups: Dict[str, List[str]]
) -> pd.DataFrame:
    """
    Partition variance explained by different covariate groups.
    
    Parameters
    ----------
    Y : np.ndarray
        Expression matrix
    design_matrix : pd.DataFrame
        Full design matrix
    covariate_groups : dict
        {group_name: [column_names]}
        
    Returns
    -------
    var_explained : pd.DataFrame
        Variance partition results per gene
    """
    n_genes = Y.shape[1]
    results = {group: np.zeros(n_genes) for group in covariate_groups}
    
    lr = LinearRegression()
    
    for g in range(n_genes):
        y = Y[:, g]
        total_var = np.var(y)
        
        for group_name, cols in covariate_groups.items():
            X_group = design_matrix[cols].values
            lr.fit(X_group, y)
            residuals = y - lr.predict(X_group)
            var_explained = 1 - (np.var(residuals) / total_var)
            results[group_name][g] = var_explained
    
    return pd.DataFrame(results)


if __name__ == "__main__":
    # Example usage
    print("Residualization module loaded successfully")
    print("Key functions: residualize_expr, fit_sva, create_design_matrix")
