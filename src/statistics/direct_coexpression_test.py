#!/usr/bin/env python3
"""
Direct Differential Co-expression Diagnostic (Phase 10)
========================================================

Permutation-based aggregate diagnostic to determine whether LIONESS-derived
rewiring rankings reflect group-level covariance changes or are primarily
artifacts of the sample-specific network pipeline.

For each Age×Arm stratum:
  1. Compute per-edge Δz = atanh(r_FLT) − atanh(r_GC) on the shared skeleton.
     Individual edge p-values are NOT trusted (n=5 gives SE≈1).
  2. Aggregate to gene-level: S_g = mean(|Δz_e|) over incident edges.
  3. Build a label-permutation null (K shuffles of FLT/GC within stratum)
     recomputing the entire gene score each time.
  4. Compare direct-correlation gene ranking to LIONESS/node2vec rewiring
     ranking (Spearman ρ), with a permutation null for the ρ itself.
  5. Run pathway enrichment on the direct ranking and compare with
     LIONESS-derived pathway results.

The skeleton was built using all samples (including FLT/GC), so this is
labeled as a robustness diagnostic, not primary inference.

Usage:
    python -m src.statistics.direct_coexpression_test \\
        --rtech  data/results/run_XYZ/phase1_residuals/Rtech.tsv.gz \\
        --meta   data/results/run_XYZ/phase1_residuals/meta_phase1.tsv.gz \\
        --phase2_dir data/results/run_XYZ/networks \\
        --rewiring_dir data/results/run_XYZ/phase3_rewiring \\
        --outdir data/results/run_XYZ/phase10_direct_coexpr \\
        --n_perms 1000
"""
from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats as sp_stats

from src.common import REPO_ROOT, find_sample_col, normalize_labels, bh_fdr

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ── Strata ────────────────────────────────────────────────────────────────

STRATA = {
    "ISS_T_YNG": {"Age": "YNG", "Arm": "ISS-T"},
    "ISS_T_OLD": {"Age": "OLD", "Arm": "ISS-T"},
    "LAR_YNG":   {"Age": "YNG", "Arm": "LAR"},
    "LAR_OLD":   {"Age": "OLD", "Arm": "LAR"},
}

CONTRAST_REWIRING = {
    "ISS_T_YNG": "ISS_T_YNG_FLT_minus_GC",
    "ISS_T_OLD": "ISS_T_OLD_FLT_minus_GC",
    "LAR_YNG":   "LAR_YNG_FLT_minus_GC",
    "LAR_OLD":   "LAR_OLD_FLT_minus_GC",
}


# ── Core computation ──────────────────────────────────────────────────────

def _pearson_edges(X: np.ndarray, ei: np.ndarray, ej: np.ndarray) -> np.ndarray:
    """Pearson r for each edge across samples.  X is (genes, samples)."""
    Xi = X[ei, :]
    Xj = X[ej, :]
    Xi = Xi - Xi.mean(axis=1, keepdims=True)
    Xj = Xj - Xj.mean(axis=1, keepdims=True)
    num = (Xi * Xj).sum(axis=1)
    den = np.sqrt((Xi**2).sum(axis=1) * (Xj**2).sum(axis=1))
    r = np.divide(num, den, out=np.zeros(len(ei), dtype=np.float64), where=den > 1e-12)
    return np.clip(r, -0.9999, 0.9999)


def _gene_scores(delta_z: np.ndarray, ei: np.ndarray, ej: np.ndarray,
                 n_genes: int) -> np.ndarray:
    """Aggregate edge-level |Δz| to gene-level mean across incident edges."""
    abs_dz = np.abs(delta_z)
    scores = np.zeros(n_genes, dtype=np.float64)
    counts = np.zeros(n_genes, dtype=np.int64)
    # Accumulate via numpy add.at for speed
    np.add.at(scores, ei, abs_dz)
    np.add.at(scores, ej, abs_dz)
    np.add.at(counts, ei, 1)
    np.add.at(counts, ej, 1)
    with np.errstate(invalid="ignore"):
        scores = np.where(counts > 0, scores / counts, 0.0)
    return scores


def compute_observed(X: np.ndarray, flt_idx: np.ndarray, gc_idx: np.ndarray,
                     ei: np.ndarray, ej: np.ndarray,
                     n_genes: int) -> tuple[np.ndarray, np.ndarray]:
    """Compute observed Δz per edge and gene-level scores."""
    r_flt = _pearson_edges(X[:, flt_idx], ei, ej)
    r_gc  = _pearson_edges(X[:, gc_idx],  ei, ej)
    delta_z = np.arctanh(r_flt) - np.arctanh(r_gc)
    scores = _gene_scores(delta_z, ei, ej, n_genes)
    return delta_z, scores


def permutation_null(X: np.ndarray, all_idx: np.ndarray,
                     n_flt: int, ei: np.ndarray, ej: np.ndarray,
                     n_genes: int, n_perms: int,
                     rng: np.random.Generator) -> np.ndarray:
    """Build gene-score null distribution by shuffling FLT/GC labels.

    Returns (n_perms, n_genes) array of null gene scores.
    """
    null_scores = np.empty((n_perms, n_genes), dtype=np.float64)
    for k in range(n_perms):
        perm = rng.permutation(all_idx)
        flt_p = perm[:n_flt]
        gc_p  = perm[n_flt:]
        _, scores_k = compute_observed(X, flt_p, gc_p, ei, ej, n_genes)
        null_scores[k, :] = scores_k
        if (k + 1) % 200 == 0:
            print(f"    Permutation {k+1}/{n_perms}")
    return null_scores


# ── Stratum runner ────────────────────────────────────────────────────────

def run_stratum(
    stratum_key: str,
    X: np.ndarray,
    sample_labels: np.ndarray,
    ei: np.ndarray,
    ej: np.ndarray,
    genes: list[str],
    n_perms: int,
    rewiring_dir: Path | None,
    outdir: Path,
    id_map: dict[str, str] | None = None,
    rng: np.random.Generator | None = None,
) -> dict:
    """Run direct differential coexpression diagnostic for one stratum."""
    if rng is None:
        rng = np.random.default_rng(42)

    flt_mask = sample_labels == "FLT"
    gc_mask  = sample_labels == "GC"
    n_flt = flt_mask.sum()
    n_gc  = gc_mask.sum()
    n_genes = len(genes)

    print(f"\n{'─'*60}")
    print(f"  Stratum: {stratum_key}  (FLT n={n_flt}, GC n={n_gc})")
    print(f"{'─'*60}")

    if n_flt < 3 or n_gc < 3:
        print(f"  SKIP: insufficient samples")
        return {"stratum": stratum_key, "status": "skipped"}

    flt_idx = np.where(flt_mask)[0]
    gc_idx  = np.where(gc_mask)[0]
    all_idx = np.concatenate([flt_idx, gc_idx])

    # ── Observed ──────────────────────────────────────────────────────
    delta_z_obs, scores_obs = compute_observed(X, flt_idx, gc_idx, ei, ej, n_genes)

    print(f"  Observed mean gene score (mean |Δz|): {scores_obs.mean():.4f}")
    print(f"  Observed max  gene score:             {scores_obs.max():.4f}")

    # ── Permutation null ──────────────────────────────────────────────
    print(f"  Running {n_perms} label permutations...")
    null_scores = permutation_null(X, all_idx, n_flt, ei, ej, n_genes, n_perms, rng)

    # Gene-level permutation p-values (one-sided: is observed score extreme?)
    perm_p = np.mean(null_scores >= scores_obs[np.newaxis, :], axis=0)
    perm_p = np.clip(perm_p, 1.0 / (n_perms + 1), 1.0)
    perm_q = bh_fdr(perm_p)

    n_sig_005 = (perm_q < 0.05).sum()
    n_sig_010 = (perm_q < 0.10).sum()
    print(f"  Genes with perm q<0.05: {n_sig_005}")
    print(f"  Genes with perm q<0.10: {n_sig_010}")

    # Global test: is the mean observed gene score higher than expected?
    mean_obs = scores_obs.mean()
    null_means = null_scores.mean(axis=1)
    global_p = np.mean(null_means >= mean_obs)
    print(f"  Global test (mean gene score): obs={mean_obs:.4f}, "
          f"null mean={null_means.mean():.4f}, p={global_p:.4f}")

    # ── LIONESS concordance ───────────────────────────────────────────
    concordance = {}
    contrast_name = CONTRAST_REWIRING.get(stratum_key)
    if rewiring_dir and contrast_name:
        rew_file = rewiring_dir / f"{contrast_name}_rewiring_agg.tsv"
        if rew_file.exists():
            rew = pd.read_csv(rew_file, sep="\t")
            rew_map = dict(zip(rew["gene"], rew["rewiring_mean"]))
            # Build paired arrays
            direct_vals = []
            lioness_vals = []
            for i, g in enumerate(genes):
                if g in rew_map:
                    direct_vals.append(scores_obs[i])
                    lioness_vals.append(rew_map[g])
            direct_arr = np.array(direct_vals)
            lioness_arr = np.array(lioness_vals)
            n_paired = len(direct_arr)

            if n_paired > 10:
                rho_obs, rho_p = sp_stats.spearmanr(direct_arr, lioness_arr)
                concordance["spearman_rho"] = float(rho_obs)
                concordance["spearman_p_parametric"] = float(rho_p)
                concordance["n_genes_compared"] = n_paired

                # Permutation null for the concordance ρ itself:
                # Shuffle direct scores and recompute ρ with LIONESS
                null_rhos = np.empty(n_perms)
                for k in range(n_perms):
                    perm_scores = null_scores[k, :]
                    perm_direct = np.array([perm_scores[i] for i, g in enumerate(genes) if g in rew_map])
                    null_rhos[k], _ = sp_stats.spearmanr(perm_direct, lioness_arr)

                concordance["rho_perm_p"] = float(np.mean(np.abs(null_rhos) >= abs(rho_obs)))
                concordance["null_rho_mean"] = float(null_rhos.mean())
                concordance["null_rho_std"] = float(null_rhos.std())

                print(f"  LIONESS concordance: ρ = {rho_obs:.4f}")
                print(f"    Permutation p(|ρ_null| ≥ |ρ_obs|) = {concordance['rho_perm_p']:.4f}")
                print(f"    Null ρ: mean={null_rhos.mean():.4f} ± {null_rhos.std():.4f}")

                # Top-decile overlap
                n_top = max(1, n_paired // 10)
                top_direct = set(np.array([g for g in genes if g in rew_map])[
                    np.argsort(direct_arr)[-n_top:]])
                top_lioness = set(np.array([g for g in genes if g in rew_map])[
                    np.argsort(lioness_arr)[-n_top:]])
                overlap = len(top_direct & top_lioness)
                jaccard = overlap / len(top_direct | top_lioness) if top_direct | top_lioness else 0

                # Null for Jaccard overlap
                null_overlaps = []
                for k in range(min(n_perms, 500)):
                    perm_scores_k = null_scores[k, :]
                    perm_direct_k = np.array([perm_scores_k[i] for i, g in enumerate(genes) if g in rew_map])
                    top_perm = set(np.array([g for g in genes if g in rew_map])[
                        np.argsort(perm_direct_k)[-n_top:]])
                    null_overlaps.append(len(top_perm & top_lioness))
                null_overlap_mean = np.mean(null_overlaps)

                concordance["top_decile_overlap"] = overlap
                concordance["top_decile_jaccard"] = float(jaccard)
                concordance["top_decile_n"] = n_top
                concordance["top_decile_null_overlap_mean"] = float(null_overlap_mean)

                print(f"    Top-decile overlap: {overlap}/{n_top} "
                      f"(null mean: {null_overlap_mean:.1f})")

    # ── Save outputs ──────────────────────────────────────────────────
    sdir = outdir / stratum_key
    sdir.mkdir(parents=True, exist_ok=True)

    gene_df = pd.DataFrame({
        "gene": genes,
        "direct_score": scores_obs,
        "perm_p": perm_p,
        "perm_q": perm_q,
        "direct_rank": sp_stats.rankdata(-scores_obs).astype(int),
    })
    if id_map:
        gene_df["symbol"] = gene_df["gene"].map(id_map).fillna("")
    gene_df = gene_df.sort_values("direct_score", ascending=False).reset_index(drop=True)
    gene_df.to_csv(sdir / "gene_level_direct_coexpr.tsv", sep="\t", index=False)

    # Edge-level observed Δz (compressed: save only top-magnitude edges)
    top_edge_n = min(5000, len(delta_z_obs))
    top_edge_idx = np.argsort(np.abs(delta_z_obs))[-top_edge_n:]
    edge_df = pd.DataFrame({
        "gene_i": [genes[ei[k]] for k in top_edge_idx],
        "gene_j": [genes[ej[k]] for k in top_edge_idx],
        "delta_z": delta_z_obs[top_edge_idx],
        "abs_delta_z": np.abs(delta_z_obs[top_edge_idx]),
    }).sort_values("abs_delta_z", ascending=False)
    edge_df.to_csv(sdir / "top_edges_by_delta_z.tsv", sep="\t", index=False)

    result = {
        "stratum": stratum_key,
        "status": "ok",
        "n_flt": int(n_flt),
        "n_gc": int(n_gc),
        "n_perms": n_perms,
        "mean_gene_score_obs": float(mean_obs),
        "mean_gene_score_null": float(null_means.mean()),
        "global_p": float(global_p),
        "n_genes_perm_q005": int(n_sig_005),
        "n_genes_perm_q010": int(n_sig_010),
    }
    result.update({f"concordance_{k}": v for k, v in concordance.items()})
    return result


# ── Main ──────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Direct differential co-expression diagnostic (permutation-based)")
    ap.add_argument("--rtech", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--phase2_dir", required=True)
    ap.add_argument("--rewiring_dir", default="")
    ap.add_argument("--id_map", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--n_perms", type=int, default=1000,
                    help="Label permutations per stratum (default: 1000)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    phase2 = Path(args.phase2_dir)
    rewiring_dir = Path(args.rewiring_dir) if args.rewiring_dir else None
    rng = np.random.default_rng(args.seed)

    print("=" * 60)
    print("Phase 10: Direct Differential Co-expression Diagnostic")
    print(f"  Permutations: {args.n_perms}")
    print("=" * 60)

    # Load skeleton
    genes_file = phase2 / "phase2_genes.txt"
    if not genes_file.exists():
        genes_file = phase2 / "phase2_genes_effective.txt"
    genes = [g.strip() for g in genes_file.read_text().splitlines() if g.strip()]
    gene_to_idx = {g: i for i, g in enumerate(genes)}

    edges = pd.read_csv(phase2 / "skeleton_edges.tsv", sep="\t")
    valid = edges["gene_i"].isin(gene_to_idx) & edges["gene_j"].isin(gene_to_idx)
    ei = edges.loc[valid, "gene_i"].map(gene_to_idx).astype(int).values
    ej = edges.loc[valid, "gene_j"].map(gene_to_idx).astype(int).values
    print(f"  Genes: {len(genes)}, Edges: {len(ei)}")

    # Load expression + metadata
    print(f"\nLoading expression: {args.rtech}")
    rtech = pd.read_csv(args.rtech, sep="\t", compression="gzip", index_col=0)
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    meta = normalize_labels(meta)
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]
    print(f"  Aligned: {len(common)} samples")

    # Filter genes
    missing = [g for g in genes if g not in rtech.index]
    if missing:
        print(f"  WARNING: {len(missing)} genes missing, filtering")
        genes = [g for g in genes if g in rtech.index]
        gene_to_idx = {g: i for i, g in enumerate(genes)}
        valid = edges["gene_i"].isin(gene_to_idx) & edges["gene_j"].isin(gene_to_idx)
        ei = edges.loc[valid, "gene_i"].map(gene_to_idx).astype(int).values
        ej = edges.loc[valid, "gene_j"].map(gene_to_idx).astype(int).values

    X = rtech.loc[genes].values.astype(np.float64)
    sample_names = list(rtech.columns)
    print(f"  Matrix: {X.shape[0]} genes × {X.shape[1]} samples")

    # ID map
    id_map = None
    if args.id_map and Path(args.id_map).exists():
        id_df = pd.read_csv(args.id_map, sep="\t", comment="#")
        eid_col = [c for c in id_df.columns if "ensembl" in c.lower()][0]
        sym_col = [c for c in id_df.columns if "symbol" in c.lower() or "mgi" in c.lower()][0]
        id_map = dict(zip(id_df[eid_col].astype(str), id_df[sym_col].astype(str)))

    # Run strata
    all_results = []
    for stratum_key, filters in STRATA.items():
        mask = np.ones(len(meta), dtype=bool)
        for col, val in filters.items():
            mask &= (meta[col].values == val)
        stratum_meta = meta[mask]
        stratum_samples = list(stratum_meta.index)
        sample_indices = [sample_names.index(s) for s in stratum_samples if s in sample_names]
        if not sample_indices:
            print(f"\n  SKIP {stratum_key}: no samples")
            continue
        X_stratum = X[:, sample_indices]
        labels = stratum_meta.loc[[sample_names[i] for i in sample_indices], "EnvGroup"].values

        result = run_stratum(
            stratum_key, X_stratum, labels, ei, ej, genes,
            args.n_perms, rewiring_dir, outdir, id_map, rng,
        )
        all_results.append(result)

    # Summary
    summary = pd.DataFrame(all_results)
    summary.to_csv(outdir / "direct_coexpr_summary.tsv", sep="\t", index=False)

    meta_out = {
        "timestamp": datetime.now().isoformat(),
        "n_genes": len(genes), "n_edges": int(len(ei)),
        "n_perms": args.n_perms, "seed": args.seed,
        "method": "Permutation-based aggregate diagnostic: gene-level mean|Δz| "
                  "with label-shuffled null. Concordance with LIONESS ranking "
                  "tested via permutation ρ.",
        "note": "Skeleton built from all samples (diagnostic, not primary inference).",
    }
    with open(outdir / "direct_coexpr_metadata.json", "w") as f:
        json.dump(meta_out, f, indent=2)

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for r in all_results:
        if r["status"] == "skipped":
            print(f"  {r['stratum']}: SKIPPED")
            continue
        rho_str = ""
        if "concordance_spearman_rho" in r:
            rho_str = (f", LIONESS ρ={r['concordance_spearman_rho']:.4f} "
                       f"(perm p={r.get('concordance_rho_perm_p', 'N/A')})")
        print(f"  {r['stratum']}: global p={r['global_p']:.4f}, "
              f"{r['n_genes_perm_q005']} genes q<.05{rho_str}")
    print(f"\nOutputs: {outdir}")


if __name__ == "__main__":
    main()
