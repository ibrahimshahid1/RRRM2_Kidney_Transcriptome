# scripts/phase6_edge_regression_full.py
"""
Phase 6: Edge-level regression using ALL 80 samples (full factorial model).

Fits a FULL 2×2×2 factorial model per edge to properly control for all
structure (Age, Arm, Flight) and extract clean effect estimates.

Model per edge (FLT vs GC only):
  w_e ~ Age + Arm + Flight + Age:Flight + Arm:Flight + Age:Arm + Age:Arm:Flight

Key outputs:
  - Flight main effect (FLT vs GC, controlling for Age and Arm)
  - Age×Flight interaction (does flight effect differ by age?)
  - Arm×Flight interaction (does flight effect differ by arm?)
  - Age×Arm×Flight (3-way interaction)

Then aggregate edge-level p-values → gene-level p-values using Stouffer's Z,
and apply BH-FDR across genes.

Inputs:
  - data/processed/networks/phase2/lioness_z_edges.npy (N × E edges)
  - data/processed/networks/phase2/edge_i.npy, edge_j.npy
  - data/processed/networks/phase2/phase2_genes.txt
  - data/processed/networks/phase2/lioness_samples.txt
  - data/processed/phase1_residuals/meta_phase1.tsv.gz

Outputs:
  - data/results/phase6_regression/gene_<effect>.tsv for each effect
  - data/results/phase6_regression/edge_regression_stats.tsv (optional)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

REPO_ROOT = Path(__file__).resolve().parents[1]


def find_sample_col(meta: pd.DataFrame) -> str:
    """Find sample identifier column."""
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            return col
    return meta.columns[0]


def normalize_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Normalize Age, Arm, EnvGroup labels."""
    meta = meta.copy()
    meta["Age"] = meta["Age"].astype(str).str.upper().replace({
        "YOUNG": "YNG", "YNG": "YNG", "OLD": "OLD"
    })
    meta["Arm"] = meta["Arm"].astype(str).str.upper().replace({
        "ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T", "ISS-T": "ISS-T",
        "LAR_T": "LAR", "LAR-T": "LAR", "LAR": "LAR"
    })
    meta["EnvGroup"] = meta["EnvGroup"].astype(str).str.upper().replace({
        "FLIGHT": "FLT", "FLT": "FLT",
        "GROUND CONTROL": "GC", "GC": "GC",
        "HGC": "HGC", "VIV": "VIV", "VGC": "VIV"
    })
    return meta


def stouffer_combine(pvals: np.ndarray) -> float:
    """Combine p-values using Stouffer's Z method."""
    pvals = np.asarray(pvals)
    pvals = pvals[~np.isnan(pvals)]
    pvals = np.clip(pvals, 1e-300, 1 - 1e-10)  # Avoid exactly 0 or 1
    
    if len(pvals) == 0:
        return np.nan
    
    # Convert p-values to z-scores (one-tailed)
    z = stats.norm.ppf(1 - pvals)
    
    # Stouffer's Z
    combined_z = np.sum(z) / np.sqrt(len(z))
    
    # Convert back to p-value (one-tailed)
    combined_p = 1 - stats.norm.cdf(combined_z)
    
    return combined_p


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def fit_edge_regression(W: np.ndarray, X: np.ndarray, term_indices: dict) -> dict:
    """
    Fit OLS regression for all edges simultaneously.
    
    Args:
        W: Edge weight matrix (N_samples × E_edges)
        X: Design matrix (N_samples × P_terms)
        term_indices: Dict mapping term names to column indices in X
    
    Returns:
        Dict with keys like 'p_flight', 'p_age_flight', 'beta_flight', etc.
    """
    N, E = W.shape
    P = X.shape[1]
    
    # Solve OLS: β = (X'X)^-1 X'W for all edges at once
    # β is P x E
    XtX_inv = np.linalg.pinv(X.T @ X)
    beta = XtX_inv @ X.T @ W  # P x E
    
    # Residuals and residual variance
    predicted = X @ beta  # N x E
    residuals = W - predicted
    df_resid = N - P
    
    # Residual variance per edge
    var_resid = np.sum(residuals ** 2, axis=0) / df_resid  # E
    
    # Standard errors of coefficients
    # SE[β_j] = sqrt(var_resid * (X'X)^-1[j,j])
    diag_XtX_inv = np.diag(XtX_inv)  # P
    
    results = {}
    
    for term_name, term_idx in term_indices.items():
        b = beta[term_idx, :]  # E
        se = np.sqrt(var_resid * diag_XtX_inv[term_idx])  # E
        
        # t-statistic
        with np.errstate(divide='ignore', invalid='ignore'):
            t_stat = b / se
        
        # Two-tailed p-value
        p_vals = 2 * (1 - stats.t.cdf(np.abs(t_stat), df_resid))
        
        results[f"beta_{term_name}"] = b
        results[f"se_{term_name}"] = se
        results[f"t_{term_name}"] = t_stat
        results[f"p_{term_name}"] = p_vals
    
    return results


def main():
    ap = argparse.ArgumentParser(description="Edge-level regression using all 80 samples")
    ap.add_argument("--phase2_dir",
                    default=str(REPO_ROOT / "data/processed/networks/phase2"),
                    help="Directory with LIONESS edge weights")
    ap.add_argument("--meta",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"),
                    help="Metadata file")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/results/phase6_regression"),
                    help="Output directory")
    ap.add_argument("--save_edges", action="store_true",
                    help="Save edge-level statistics (large file)")
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load edge weight matrix
    print("Loading edge weight matrix...")
    
    z_path = phase2 / "lioness_z_edges.npy"
    if not z_path.exists():
        raise FileNotFoundError(f"Edge weight matrix not found: {z_path}")
    
    W = np.load(z_path)  # N_samples x E_edges
    N, E = W.shape
    print(f"  Edge weights: {N} samples × {E} edges")

    # Load edge indices
    edge_i = np.load(phase2 / "edge_i.npy")
    edge_j = np.load(phase2 / "edge_j.npy")
    print(f"  Edge index arrays loaded")

    # Load genes
    genes_path = phase2 / "phase2_genes.txt"
    genes = [g.strip() for g in genes_path.read_text().splitlines() if g.strip()]
    G = len(genes)
    print(f"  {G} genes")

    # Load sample order
    samples_path = phase2 / "lioness_samples.txt"
    samples = [s.strip() for s in samples_path.read_text().splitlines() if s.strip()]
    if len(samples) != N:
        raise ValueError(f"Sample count mismatch: {len(samples)} vs {N}")

    # 2. Load and align metadata
    print("\nLoading metadata...")
    
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    meta = normalize_labels(meta)
    
    # Align to sample order
    meta = meta.loc[samples].copy()
    print(f"  Aligned {len(meta)} samples")

    # Filter to FLT and GC samples only (exclude VIV, HGC if present)
    mask_flt_gc = meta["EnvGroup"].isin(["FLT", "GC"])
    print(f"  FLT+GC samples: {mask_flt_gc.sum()}")
    
    meta_sub = meta[mask_flt_gc].copy()
    W_sub = W[mask_flt_gc.values, :]
    
    N_sub = len(meta_sub)
    print(f"  Using {N_sub} samples for regression")

    # Print sample breakdown
    print("\nSample breakdown:")
    for (age, arm, env), count in meta_sub.groupby(["Age", "Arm", "EnvGroup"]).size().items():
        print(f"  {age} × {arm} × {env}: {count}")

    # 3. Build design matrix (FULL 2×2×2 factorial)
    print("\nBuilding FULL factorial design matrix...")
    print("  Model: w ~ Age + Arm + Flight + Age:Flight + Arm:Flight + Age:Arm + Age:Arm:Flight")
    
    # Binary encodings (centered not needed for factorial, just 0/1)
    flight = (meta_sub["EnvGroup"] == "FLT").astype(float).values  # 1=FLT, 0=GC
    age = (meta_sub["Age"] == "OLD").astype(float).values  # 1=OLD, 0=YNG
    arm = (meta_sub["Arm"] == "ISS-T").astype(float).values  # 1=ISS-T, 0=LAR
    
    # All interaction terms
    age_x_flight = age * flight
    arm_x_flight = arm * flight
    age_x_arm = age * arm
    age_x_arm_x_flight = age * arm * flight  # 3-way interaction
    
    # Full factorial design matrix (8 terms including intercept)
    X = np.column_stack([
        np.ones(N_sub),    # 0: intercept
        age,               # 1: age main effect
        arm,               # 2: arm main effect  
        flight,            # 3: flight main effect
        age_x_flight,      # 4: age × flight
        arm_x_flight,      # 5: arm × flight
        age_x_arm,         # 6: age × arm
        age_x_arm_x_flight # 7: age × arm × flight (3-way)
    ])
    
    term_indices = {
        "age": 1,
        "arm": 2,
        "flight": 3,
        "age_flight": 4,
        "arm_flight": 5,
        "age_arm": 6,
        "age_arm_flight": 7,
    }
    
    print(f"  Design matrix: {X.shape} ({X.shape[1]} terms including intercept)")
    print(f"  Terms: {list(term_indices.keys())}")

    # 4. Fit edge-level regressions
    print("\nFitting edge-level regressions...")
    
    edge_results = fit_edge_regression(W_sub, X, term_indices)
    
    print(f"  Flight effect: median p = {np.nanmedian(edge_results['p_flight']):.4f}")
    print(f"  Age×Flight: median p = {np.nanmedian(edge_results['p_age_flight']):.4f}")
    print(f"  Arm×Flight: median p = {np.nanmedian(edge_results['p_arm_flight']):.4f}")
    print(f"  Age×Arm×Flight: median p = {np.nanmedian(edge_results['p_age_arm_flight']):.4f}")

    # 5. Aggregate to gene-level using Stouffer's Z
    print("\nAggregating to gene-level (Stouffer's Z)...")
    
    # Key terms to report at gene level
    report_terms = ["flight", "age_flight", "arm_flight", "age_arm_flight"]
    
    gene_results = {term: {"p_stouffer": [], "n_edges": [], "median_beta": [], "median_t": []} 
                    for term in report_terms}
    
    for g_idx, gene in enumerate(genes):
        # Find edges incident to this gene
        incident = (edge_i == g_idx) | (edge_j == g_idx)
        n_edges = incident.sum()
        
        for term in report_terms:
            if n_edges == 0:
                gene_results[term]["p_stouffer"].append(np.nan)
                gene_results[term]["n_edges"].append(0)
                gene_results[term]["median_beta"].append(np.nan)
                gene_results[term]["median_t"].append(np.nan)
            else:
                p_edges = edge_results[f"p_{term}"][incident]
                beta_edges = edge_results[f"beta_{term}"][incident]
                t_edges = edge_results[f"t_{term}"][incident]
                
                p_combined = stouffer_combine(p_edges)
                
                gene_results[term]["p_stouffer"].append(p_combined)
                gene_results[term]["n_edges"].append(n_edges)
                gene_results[term]["median_beta"].append(np.nanmedian(beta_edges))
                gene_results[term]["median_t"].append(np.nanmedian(t_edges))
    
    # 6. Apply BH-FDR and save
    print("\nApplying BH-FDR and saving results...")
    
    # Term name to output file label mapping
    term_labels = {
        "flight": "flight_effect",
        "age_flight": "age_flight_interaction",
        "arm_flight": "arm_flight_interaction",
        "age_arm_flight": "age_arm_flight_3way",
    }
    
    for term in report_terms:
        p_vals = np.array(gene_results[term]["p_stouffer"])
        
        # Handle NaN for FDR
        valid_mask = ~np.isnan(p_vals)
        q_vals = np.full_like(p_vals, np.nan)
        if valid_mask.sum() > 0:
            q_vals[valid_mask] = bh_fdr(p_vals[valid_mask])
        
        df = pd.DataFrame({
            "gene": genes,
            "n_edges": gene_results[term]["n_edges"],
            "median_beta": gene_results[term]["median_beta"],
            "median_t": gene_results[term]["median_t"],
            "p_stouffer": p_vals,
            "q_BH": q_vals,
        })
        df = df.sort_values("p_stouffer").reset_index(drop=True)
        
        term_label = term_labels[term]
        out_path = outdir / f"gene_{term_label}.tsv"
        df.to_csv(out_path, sep="\t", index=False)
        
        n_sig = (q_vals < 0.05).sum()
        n_sig_01 = (q_vals < 0.01).sum()
        print(f"\n  {term_label}:")
        print(f"    Wrote: {out_path}")
        print(f"    Significant (q < 0.05): {n_sig}")
        print(f"    Significant (q < 0.01): {n_sig_01}")

    # 7. Optional: save edge-level stats
    if args.save_edges:
        print("\nSaving edge-level statistics...")
        edge_df = pd.DataFrame({
            "edge_i": edge_i,
            "edge_j": edge_j,
            "gene_i": [genes[i] for i in edge_i],
            "gene_j": [genes[j] for j in edge_j],
            "beta_flight": edge_results["beta_flight"],
            "p_flight": edge_results["p_flight"],
            "beta_age_flight": edge_results["beta_age_flight"],
            "p_age_flight": edge_results["p_age_flight"],
        })
        edge_path = outdir / "edge_regression_stats.tsv"
        edge_df.to_csv(edge_path, sep="\t", index=False)
        print(f"  Wrote: {edge_path}")

    print(f"\n[OK] All outputs written to: {outdir}")


if __name__ == "__main__":
    main()
