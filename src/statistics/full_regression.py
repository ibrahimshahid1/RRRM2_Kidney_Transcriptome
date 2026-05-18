# src/statistics/full_regression.py
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

Then aggregate incident edge t-statistics into signed gene-level statistics and
calibrate them by within-Age×Arm label permutation. Stouffer/Brown-style
combination is not used for primary inference because unsigned edge p-values
discard direction and edge dependence is substantial on a shared skeleton.

Inputs:
  - data/processed/networks/phase2/lioness_edges.npy (N × E edge weights)
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

from src.common import REPO_ROOT, find_sample_col, normalize_labels, bh_fdr

warnings.filterwarnings("ignore", category=RuntimeWarning)


DEFAULT_GENE_PERMUTATIONS = 1000


def stouffer_combine(pvals: np.ndarray) -> float:
    """Legacy unsigned Stouffer diagnostic; never used as primary inference."""
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


def build_factorial_design(meta_sub: pd.DataFrame, flight_override: np.ndarray | None = None) -> tuple[np.ndarray, dict]:
    """Build the full 2×2×2 design matrix for FLT vs GC samples."""
    if flight_override is None:
        flight = (meta_sub["EnvGroup"] == "FLT").astype(float).values
    else:
        flight = np.asarray(flight_override, dtype=float)
        if flight.shape[0] != len(meta_sub):
            raise ValueError("flight_override length does not match metadata")

    age = (meta_sub["Age"] == "OLD").astype(float).values
    arm = (meta_sub["Arm"] == "ISS-T").astype(float).values

    age_x_flight = age * flight
    arm_x_flight = arm * flight
    age_x_arm = age * arm
    age_x_arm_x_flight = age * arm * flight

    X = np.column_stack([
        np.ones(len(meta_sub)),
        age,
        arm,
        flight,
        age_x_flight,
        arm_x_flight,
        age_x_arm,
        age_x_arm_x_flight,
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
    return X, term_indices


def signed_incident_t_aggregate(
    edge_results: dict,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    n_genes: int,
    term: str,
) -> dict[str, np.ndarray]:
    """Aggregate incident edge t-statistics as sum(t)/sqrt(degree)."""
    t = np.asarray(edge_results[f"t_{term}"], dtype=float)
    beta = np.asarray(edge_results[f"beta_{term}"], dtype=float)
    valid_t = np.where(np.isfinite(t), t, 0.0)

    signed_sum = np.zeros(n_genes, dtype=float)
    degree = np.zeros(n_genes, dtype=int)
    np.add.at(signed_sum, edge_i, valid_t)
    np.add.at(signed_sum, edge_j, valid_t)
    np.add.at(degree, edge_i, 1)
    np.add.at(degree, edge_j, 1)
    stat = np.full(n_genes, np.nan, dtype=float)
    nonzero = degree > 0
    stat[nonzero] = signed_sum[nonzero] / np.sqrt(degree[nonzero])

    median_beta = np.full(n_genes, np.nan, dtype=float)
    median_t = np.full(n_genes, np.nan, dtype=float)
    for g_idx in range(n_genes):
        incident = (edge_i == g_idx) | (edge_j == g_idx)
        if incident.any():
            median_beta[g_idx] = np.nanmedian(beta[incident])
            median_t[g_idx] = np.nanmedian(t[incident])
    return {
        "signed_t_sum_sqrt_degree": stat,
        "n_edges": degree,
        "median_beta": median_beta,
        "median_t": median_t,
    }


def permute_flight_within_age_arm(meta_sub: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    """Permute FLT/GC labels within each Age×Arm stratum."""
    labels = (meta_sub["EnvGroup"] == "FLT").astype(int).values.copy()
    out = labels.copy()
    strata = meta_sub["Age"].astype(str) + "|" + meta_sub["Arm"].astype(str)
    for stratum in strata.unique():
        idx = np.where(strata.values == stratum)[0]
        perm = labels[idx].copy()
        rng.shuffle(perm)
        out[idx] = perm
    return out.astype(float)


def empirical_signed_gene_pvalues(
    W_sub: np.ndarray,
    meta_sub: pd.DataFrame,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    n_genes: int,
    report_terms: list[str],
    observed_stats: dict[str, np.ndarray],
    n_perm: int,
    seed: int,
) -> dict[str, np.ndarray]:
    """Calibrate signed gene-level statistics by label permutation."""
    rng = np.random.default_rng(seed)
    exceed = {term: np.zeros(n_genes, dtype=int) for term in report_terms}

    for k in range(n_perm):
        flight_perm = permute_flight_within_age_arm(meta_sub, rng)
        Xp, term_indices = build_factorial_design(meta_sub, flight_override=flight_perm)
        edge_perm = fit_edge_regression(W_sub, Xp, term_indices)
        for term in report_terms:
            agg = signed_incident_t_aggregate(edge_perm, edge_i, edge_j, n_genes, term)
            null_stat = agg["signed_t_sum_sqrt_degree"]
            obs = observed_stats[term]
            valid = np.isfinite(obs) & np.isfinite(null_stat)
            exceed[term][valid] += (np.abs(null_stat[valid]) >= np.abs(obs[valid]))
        if (k + 1) % max(1, n_perm // 10) == 0:
            print(f"    gene aggregation permutation {k + 1}/{n_perm}")

    pvals = {}
    for term in report_terms:
        p = (1.0 + exceed[term].astype(float)) / (n_perm + 1.0)
        p[~np.isfinite(observed_stats[term])] = np.nan
        pvals[term] = p
    return pvals


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
    ap.add_argument("--edge-weights", "--z", dest="edge_weights", default="lioness_edges.npy",
                    help="LIONESS edge-weight matrix filename under phase2_dir")
    ap.add_argument("--meta",
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"),
                    help="Metadata file")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/results/phase6_regression"),
                    help="Output directory")
    ap.add_argument("--save_edges", action="store_true",
                    help="Save edge-level statistics (large file)")
    ap.add_argument("--K_gene_perm", type=int, default=DEFAULT_GENE_PERMUTATIONS,
                    help="Within-stratum label permutations for primary signed gene aggregation")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--diagnostic_stouffer", action="store_true",
                    help="Also write legacy unsigned Stouffer p-values as a diagnostic column")
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # 1. Load edge weight matrix
    print("Loading edge weight matrix...")
    
    z_path = phase2 / args.edge_weights
    if not z_path.exists() and args.edge_weights == "lioness_edges.npy":
        legacy = phase2 / "lioness_z_edges.npy"
        if legacy.exists():
            z_path = legacy
            print("  WARNING: using legacy lioness_z_edges.npy. Prefer lioness_edges.npy for corrected defaults.")
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
    X, term_indices = build_factorial_design(meta_sub)
    
    print(f"  Design matrix: {X.shape} ({X.shape[1]} terms including intercept)")
    print(f"  Terms: {list(term_indices.keys())}")

    # 4. Fit edge-level regressions
    print("\nFitting edge-level regressions...")
    
    edge_results = fit_edge_regression(W_sub, X, term_indices)
    
    print(f"  Flight effect: median p = {np.nanmedian(edge_results['p_flight']):.4f}")
    print(f"  Age×Flight: median p = {np.nanmedian(edge_results['p_age_flight']):.4f}")
    print(f"  Arm×Flight: median p = {np.nanmedian(edge_results['p_arm_flight']):.4f}")
    print(f"  Age×Arm×Flight: median p = {np.nanmedian(edge_results['p_age_arm_flight']):.4f}")

    # 5. Aggregate to gene-level using signed empirical calibration
    print("\nAggregating to gene-level (signed incident-edge t statistic)...")
    
    # Key terms to report at gene level
    report_terms = ["flight", "age_flight", "arm_flight", "age_arm_flight"]
    
    gene_results = {}
    observed_stats = {}
    for term in report_terms:
        agg = signed_incident_t_aggregate(edge_results, edge_i, edge_j, G, term)
        gene_results[term] = agg
        observed_stats[term] = agg["signed_t_sum_sqrt_degree"]

    print(f"\nCalibrating signed gene aggregation with {args.K_gene_perm} within-stratum label permutations...")
    p_empirical = empirical_signed_gene_pvalues(
        W_sub=W_sub,
        meta_sub=meta_sub,
        edge_i=edge_i,
        edge_j=edge_j,
        n_genes=G,
        report_terms=report_terms,
        observed_stats=observed_stats,
        n_perm=args.K_gene_perm,
        seed=args.seed,
    )
    
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
        p_vals = np.array(p_empirical[term])
        
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
            "signed_t_sum_sqrt_degree": gene_results[term]["signed_t_sum_sqrt_degree"],
            "p_empirical_signed": p_vals,
            "q_BH_empirical_signed": q_vals,
        })
        if args.diagnostic_stouffer:
            legacy = np.full(G, np.nan, dtype=float)
            for g_idx in range(G):
                incident = (edge_i == g_idx) | (edge_j == g_idx)
                if incident.any():
                    legacy[g_idx] = stouffer_combine(edge_results[f"p_{term}"][incident])
            df["p_stouffer_legacy_unsigned_diagnostic"] = legacy
        df = df.sort_values("p_empirical_signed").reset_index(drop=True)
        
        term_label = term_labels[term]
        out_path = outdir / f"gene_{term_label}.tsv"
        df.to_csv(out_path, sep="\t", index=False)
        
        n_sig = int(np.nansum(q_vals < 0.05))
        n_sig_01 = int(np.nansum(q_vals < 0.01))
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
