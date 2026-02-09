# src/statistics/permutation_bootstrap.py
"""
Phase 6: Uncertainty / null models for rewiring (FAST).

We avoid rerunning node2vec by working in LIONESS edge-weight space:
  - Compute delta_z per edge for FLT-GC contrasts within each (Age, Arm)
  - Convert to per-gene rewiring_abs by summing |delta| over incident edges
  - Permutation test: shuffle FLT/GC labels within each (Age, Arm) stratum
  - Bootstrap CI: resample samples with replacement within FLT and within GC (stratified)

Outputs:
  data/results/phase6_uncertainty/
    <contrast>_perm_pvals.tsv
    <contrast>_bootstrap_ci.tsv

Contrast names match your pipeline:
  ISS_T_YNG_FLT_minus_GC, ISS_T_OLD_FLT_minus_GC, LAR_YNG_FLT_minus_GC, LAR_OLD_FLT_minus_GC
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# Repository root (2 levels up from src/statistics/)
REPO_ROOT = Path(__file__).resolve().parents[2]


def find_sample_col(meta: pd.DataFrame) -> str:
    """Find the sample identifier column in metadata."""
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            return col
    return meta.columns[0]


def normalize_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Normalize Age, Arm, EnvGroup labels for consistency."""
    meta = meta.copy()
    meta["Age"] = meta["Age"].astype(str).replace({
        "Young": "YNG", "Yng": "YNG", "young": "YNG",
        "Old": "OLD", "old": "OLD"
    })
    meta["Arm"] = meta["Arm"].astype(str).replace({
        "ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T", "ISS T": "ISS-T"
    })
    meta["Arm"] = meta["Arm"].astype(str).replace({
        "LAR_T": "LAR", "LAR-T": "LAR", "LAR T": "LAR"
    })
    meta["EnvGroup"] = meta["EnvGroup"].astype(str).replace({
        "HGC": "GC", "VGC": "VIV", "HGC/GC": "GC", "VIV/VGC": "VIV"
    })
    return meta


def node_rewiring_abs(delta: np.ndarray, edge_i: np.ndarray, edge_j: np.ndarray, G: int) -> np.ndarray:
    """Sum |delta| over incident edges for each gene."""
    abs_sum = np.zeros(G, dtype=np.float64)
    dabs = np.abs(delta).astype(np.float64)
    np.add.at(abs_sum, edge_i, dabs)
    np.add.at(abs_sum, edge_j, dabs)
    return abs_sum


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


def main():
    ap = argparse.ArgumentParser(description="Permutation + bootstrap uncertainty for rewiring")
    ap.add_argument("--phase2_dir", 
                    default=str(REPO_ROOT / "data/processed/networks/phase2"),
                    help="Directory with LIONESS z-scores and edge indices")
    ap.add_argument("--meta", 
                    default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"),
                    help="Metadata file")
    ap.add_argument("--z", default="lioness_z_edges.npy",
                    help="Filename of LIONESS z-scores matrix (N x E)")
    ap.add_argument("--outdir", 
                    default=str(REPO_ROOT / "data/results/phase6_uncertainty"),
                    help="Outpußt directory")
    ap.add_argument("--K_perm", type=int, default=2000,
                    help="Number of permutations")
    ap.add_argument("--B_boot", type=int, default=2000,
                    help="Number of bootstrap resamples")
    ap.add_argument("--seed", type=int, default=0,
                    help="Random seed")
    ap.add_argument("--covariates", default="",
                    help="Optional: comma-separated covariates to regress out from Z")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    phase2 = Path(args.phase2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load edges + genes
    genes_path = phase2 / "phase2_genes.txt"
    if not genes_path.exists():
        raise FileNotFoundError(f"Missing genes file: {genes_path}")
    genes = [g.strip() for g in genes_path.read_text().splitlines() if g.strip()]
    G = len(genes)
    print(f"Loaded {G} genes")

    edge_i = np.load(phase2 / "edge_i.npy")
    edge_j = np.load(phase2 / "edge_j.npy")
    E = len(edge_i)
    print(f"Loaded {E} edges")

    # Load Z (N x E) and sample order
    z_path = phase2 / args.z
    if not z_path.exists():
        raise FileNotFoundError(f"Missing Z matrix: {z_path}")
    Z = np.load(z_path)
    N = Z.shape[0]
    print(f"Loaded Z matrix: {Z.shape}")

    samples_path = phase2 / "lioness_samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(f"Missing samples file: {samples_path}")
    samples = [s.strip() for s in samples_path.read_text().splitlines() if s.strip()]
    if len(samples) != N:
        raise ValueError(f"lioness_samples.txt length ({len(samples)}) != Z rows ({N})")

    # Load + align metadata
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    meta = normalize_labels(meta)
    meta = meta.loc[samples].copy()
    print(f"Aligned metadata: {len(meta)} samples")

    # Optional: regress out covariates from Z once (fast)
    covs = [c.strip() for c in args.covariates.split(",") if c.strip() and c.strip() in meta.columns]
    if covs:
        X_parts = [np.ones((N, 1), dtype=float)]  # intercept
        for c in covs:
            v = meta[c]
            if pd.api.types.is_numeric_dtype(v):
                vv = (v.values.astype(float) - np.nanmean(v.values.astype(float))) / (np.nanstd(v.values.astype(float)) + 1e-8)
                X_parts.append(vv.reshape(-1, 1))
            else:
                d = pd.get_dummies(v.astype(str), drop_first=True)
                if d.shape[1] > 0:
                    X_parts.append(d.values.astype(float))
        X = np.concatenate(X_parts, axis=1)
        B, *_ = np.linalg.lstsq(X, Z, rcond=None)
        Z = Z - X @ B
        print(f"[OK] Residualized Z on covariates: {covs} (P={X.shape[1]})")

    # Define the four contrasts (within each Age×Arm, compare FLT vs GC)
    contrasts = [
        ("ISS_T_YNG_FLT_minus_GC", {"Age": "YNG", "Arm": "ISS-T"}),
        ("ISS_T_OLD_FLT_minus_GC", {"Age": "OLD", "Arm": "ISS-T"}),
        ("LAR_YNG_FLT_minus_GC",   {"Age": "YNG", "Arm": "LAR"}),
        ("LAR_OLD_FLT_minus_GC",   {"Age": "OLD", "Arm": "LAR"}),
    ]

    for cname, filt in contrasts:
        print(f"\n--- Processing contrast: {cname} ---")
        mask = (meta["Age"] == filt["Age"]) & (meta["Arm"] == filt["Arm"]) & (meta["EnvGroup"].isin(["FLT", "GC"]))
        sub_idx = np.where(mask.values)[0]
        sub = meta.iloc[sub_idx].copy()
        zsub = Z[sub_idx, :]  # n x E

        flt_idx = np.where(sub["EnvGroup"].values == "FLT")[0]
        gc_idx = np.where(sub["EnvGroup"].values == "GC")[0]
        
        if len(flt_idx) < 2 or len(gc_idx) < 2:
            print(f"  [SKIP] Not enough samples (FLT={len(flt_idx)} GC={len(gc_idx)})")
            continue

        print(f"  Samples: {len(sub)} total (FLT={len(flt_idx)}, GC={len(gc_idx)})")

        # Observed delta per edge (mean difference)
        delta_obs = zsub[flt_idx].mean(axis=0) - zsub[gc_idx].mean(axis=0)
        rew_obs = node_rewiring_abs(delta_obs, edge_i, edge_j, G)

        # Permutation null
        print(f"  Running {args.K_perm} permutations...")
        labels = sub["EnvGroup"].values.copy()
        null = np.zeros((args.K_perm, G), dtype=np.float32)

        for k in range(args.K_perm):
            perm = labels.copy()
            rng.shuffle(perm)
            p_flt = np.where(perm == "FLT")[0]
            p_gc = np.where(perm == "GC")[0]
            d = zsub[p_flt].mean(axis=0) - zsub[p_gc].mean(axis=0)
            null[k] = node_rewiring_abs(d, edge_i, edge_j, G).astype(np.float32)

        # One-sided p-value: P(null >= obs)
        pvals = (1.0 + (null >= rew_obs[None, :]).sum(axis=0)) / (args.K_perm + 1.0)
        qvals = bh_fdr(pvals)

        perm_out = pd.DataFrame({
            "gene": genes,
            "rewiring_abs_obs": rew_obs,
            "p_perm": pvals,
            "q_BH": qvals,
        }).sort_values("p_perm").reset_index(drop=True)
        perm_path = outdir / f"{cname}_perm_pvals.tsv"
        perm_out.to_csv(perm_path, sep="\t", index=False)
        print(f"  Wrote {perm_path}")

        # Bootstrap CI (percentile)
        print(f"  Running {args.B_boot} bootstrap resamples...")
        boot = np.zeros((args.B_boot, G), dtype=np.float32)

        for b in range(args.B_boot):
            b_flt = rng.choice(flt_idx, size=len(flt_idx), replace=True)
            b_gc = rng.choice(gc_idx, size=len(gc_idx), replace=True)
            d = zsub[b_flt].mean(axis=0) - zsub[b_gc].mean(axis=0)
            boot[b] = node_rewiring_abs(d, edge_i, edge_j, G).astype(np.float32)

        lo = np.quantile(boot, 0.025, axis=0)
        hi = np.quantile(boot, 0.975, axis=0)

        boot_out = pd.DataFrame({
            "gene": genes,
            "rewiring_abs_obs": rew_obs,
            "ci_2p5": lo,
            "ci_97p5": hi,
        }).sort_values("rewiring_abs_obs", ascending=False).reset_index(drop=True)
        boot_path = outdir / f"{cname}_bootstrap_ci.tsv"
        boot_out.to_csv(boot_path, sep="\t", index=False)
        print(f"  Wrote {boot_path}")

    print(f"\n[OK] All outputs written to: {outdir}")


if __name__ == "__main__":
    main()
