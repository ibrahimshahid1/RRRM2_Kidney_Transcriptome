"""
Edge-Wise Regression and Predicted Networks

Implements edge-wise modeling over the full factorial design (Age × Arm × Group)
with empirical Bayes variance moderation. Generates predicted contrast-specific
networks by predicting edge weights under target condition profiles.

Key Functions:
    - fit_edge_wise_model: Fit regression for each edge
    - predict_network_contrast: Generate predicted network for a contrast
    - create_contrast_matrix: Build contrast vectors for comparisons
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from sklearn.linear_model import LinearRegression
import warnings


def create_design_matrix_full(
    metadata: pd.DataFrame
) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """
    Create full factorial design matrix: Age + Arm + Group + interactions.
    
    Parameters
    ----------
    metadata : pd.DataFrame
        Must contain: Age, Arm, EnvironmentGroup
        Plus technical covariates: batch, SVs, cell_props
        
    Returns
    -------
    design : pd.DataFrame
        Full model matrix
    term_columns : dict
        {term_name: [column_names]} for variance partitioning
        
    Notes
    -----
    Model: z_e ~ Age + Arm + Group + Age:Group + Arm:Group + Age:Arm + 
                   Age:Arm:Group + Tech + Deconv + SVs
    """
    # Placeholder implementation
    # User must provide actual metadata with proper factor levels
    
    # Main effects (dummy coding)
    design_parts = []
    term_columns = {}
    
    # Age (Young=0, Old=1)
    if 'Age' in metadata.columns:
        age_dummy = pd.get_dummies(metadata['Age'], prefix='Age', drop_first=True)
        design_parts.append(age_dummy)
        term_columns['Age'] = list(age_dummy.columns)
    
    # Arm (ISS-T=0, LAR=1)
    if 'Arm' in metadata.columns:
        arm_dummy = pd.get_dummies(metadata['Arm'], prefix='Arm', drop_first=True)
        design_parts.append(arm_dummy)
        term_columns['Arm'] = list(arm_dummy.columns)
    
    # EnvironmentGroup (reference=BSL)
    if 'EnvironmentGroup' in metadata.columns:
        group_dummy = pd.get_dummies(metadata['EnvironmentGroup'], prefix='Group', drop_first=True)
        design_parts.append(group_dummy)
        term_columns['Group'] = list(group_dummy.columns)
    
    # TODO: Add interaction terms (Age:Group, Arm:Group, etc.)
    # TODO: Add technical covariates from metadata
    
    design = pd.concat(design_parts, axis=1)
    
    return design, term_columns


def fit_edge_wise_model(
    Z_samp: np.ndarray,
    design: pd.DataFrame,
    use_empirical_bayes: bool = True
) -> Dict:
    """
    Fit edge-wise regression with optional empirical Bayes moderation.
    
    Parameters
    ----------
    Z_samp : np.ndarray, shape (n_samples, n_edges)
        Fisher z-transformed sample-specific edge weights
    design : pd.DataFrame, shape (n_samples, n_predictors)
        Design matrix
    use_empirical_bayes : bool
        If True, apply variance moderation across edges (limma-style)
        
    Returns
    -------
    results : dict
        {
            'coefficients': np.ndarray (n_edges, n_predictors),
            'se': np.ndarray (n_edges, n_predictors),
            'residuals': np.ndarray (n_samples, n_edges),
            'sigma2': np.ndarray (n_edges,)  # residual variance per edge
        }
        
    Notes
    -----
    Simplified implementation. For production, use statsmodels or
    implement full limma empirical Bayes with moderated t-statistics.
    """
    X = design.values
    n_samples, n_edges = Z_samp.shape
    n_pred = X.shape[1]
    
    lr = LinearRegression()
    
    coefficients = np.zeros((n_edges, n_pred), dtype=np.float32)
    residuals = np.zeros((n_samples, n_edges), dtype=np.float32)
    sigma2 = np.zeros(n_edges, dtype=np.float32)
    
    print(f"Fitting {n_edges} edge-wise models...")
    for e in range(n_edges):
        if e % 1000 == 0 and e > 0:
            print(f"  ... fitted {e}/{n_edges} edges")
        
        y = Z_samp[:, e]
        lr.fit(X, y)
        
        coefficients[e, :] = lr.coef_
        resid = y - lr.predict(X)
        residuals[:, e] = resid
        sigma2[e] = np.var(resid)
    
    # Empirical Bayes variance moderation (simplified)
    if use_empirical_bayes:
        # Fit inverse-gamma to sigma2 distribution
        d0 = 2.0  # Prior degrees of freedom (hyperparameter)
        s0_sq = np.median(sigma2)  # Prior variance estimate
        
        # Moderated variance
        df_resid = n_samples - n_pred
        sigma2_mod = (d0 * s0_sq + df_resid * sigma2) / (d0 + df_resid)
    else:
        sigma2_mod = sigma2
    
    # Standard errors
    # SE = sqrt(sigma2 * (X^T X)^{-1}_{jj})
    XtX_inv = np.linalg.inv(X.T @ X)
    se = np.zeros((n_edges, n_pred), dtype=np.float32)
    for e in range(n_edges):
        se[e, :] = np.sqrt(sigma2_mod[e] * np.diag(XtX_inv))
    
    results = {
        'coefficients': coefficients,
        'se': se,
        'residuals': residuals,
        'sigma2': sigma2,
        'sigma2_moderated': sigma2_mod,
        'design_matrix': design
    }
    
    return results


def predict_network_contrast(
    model_results: Dict,
    condition_profile: pd.Series,
    edges: List[Tuple[int, int]]
) -> np.ndarray:
    """
    Predict edge weights under a target condition profile.
    
    Parameters
    ----------
    model_results : dict
        Output from fit_edge_wise_model
    condition_profile : pd.Series
        Target condition values (must match design matrix columns)
    edges : list of (int, int)
        Edge list E
        
    Returns
    -------
    z_predicted : np.ndarray, shape (|E|,)
        Predicted Fisher-z edge weights
        
    Example
    -------
    >>> profile = pd.Series({
    ...     'Age_Old': 1, 'Arm_LAR': 0, 'Group_FLT': 1, 'Group_HGC': 0, ...
    ... })
    >>> z_hat = predict_network_contrast(results, profile, edges)
    """
    design_cols = model_results['design_matrix'].columns
    X_pred = condition_profile[design_cols].values.reshape(1, -1)
    
    coefficients = model_results['coefficients']
    z_predicted = (X_pred @ coefficients.T).flatten()
    
    return z_predicted.astype(np.float32)


def compute_contrast(
    model_results: Dict,
    profile_1: pd.Series,
    profile_2: pd.Series
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute edge-wise contrast: profile_1 - profile_2.
    
    Parameters
    ----------
    model_results : dict
        Fitted model
    profile_1, profile_2 : pd.Series
        Condition profiles to compare
        
    Returns
    -------
    delta_z : np.ndarray
        Difference in predicted z for each edge
    se_delta : np.ndarray
        Standard error of difference
    """
    design_cols = model_results['design_matrix'].columns
    
    X1 = profile_1[design_cols].values
    X2 = profile_2[design_cols].values
    X_contrast = X1 - X2
    
    coefficients = model_results['coefficients']
    se = model_results['se']
    
    delta_z = X_contrast @ coefficients.T
    
    # SE(contrast) = sqrt(contrast^T * Var(beta) * contrast)
    # Simplified: assume independence and add variances
    se_delta = np.sqrt((X_contrast**2 @ se.T**2).T)
    
    return delta_z, se_delta


if __name__ == "__main__":
    print("Edge regression module loaded successfully")
    print("Key functions: fit_edge_wise_model, predict_network_contrast")
