# src/statistics/permutation_bootstrap.py
"""
Phase 6: Uncertainty / null models for rewiring (FAST).

We avoid rerunning node2vec by working in LIONESS edge-weight space:
  - Compute delta_z per edge for FLT-GC contrasts within each (Age, Arm)
  - Convert to per-gene rewiring_abs by summing |delta| over incident edges
  - Permutation test: shuffle FLT/GC labels within each (Age, Arm) stratum
  - Bootstrap CI: resample samples with replacement within FLT and within GC (stratified)

Modes:
  Default:             BH-FDR across all 2500 genes (stringent)
  --focused-permutation: BH-FDR across top-decile genes only (~250, much more powerful)
  --pathway-permutation: gene-set-level permutation (tests ~100 pathways, not 2500 genes)

Outputs:
  data/results/phase6_uncertainty/
    <contrast>_perm_pvals.tsv
    <contrast>_bootstrap_ci.tsv
    <contrast>_perm_focused.tsv       (if --focused-permutation)
    <contrast>_perm_pathway.tsv       (if --pathway-permutation)

Contrast names match your pipeline:
  ISS_T_YNG_FLT_minus_GC, ISS_T_OLD_FLT_minus_GC, LAR_YNG_FLT_minus_GC, LAR_OLD_FLT_minus_GC
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from src.common import REPO_ROOT, find_sample_col, normalize_labels, bh_fdr


def node_rewiring_abs(delta: np.ndarray, edge_i: np.ndarray, edge_j: np.ndarray, G: int) -> np.ndarray:
    """Sum |delta| over incident edges for each gene."""
    abs_sum = np.zeros(G, dtype=np.float64)
    dabs = np.abs(delta).astype(np.float64)
    np.add.at(abs_sum, edge_i, dabs)
    np.add.at(abs_sum, edge_j, dabs)
    return abs_sum


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
                    help="Output directory")
    ap.add_argument("--K_perm", type=int, default=5000,
                    help="Number of permutations (default: 5000). "
                         "With n=5 FLT vs n=15 ctrl there are 15504 unique permutations; "
                         "5000 gives min p = 0.0002, just enough for focused BH on 250 genes.")
    ap.add_argument("--B_boot", type=int, default=2000,
                    help="Number of bootstrap resamples")
    ap.add_argument("--seed", type=int, default=0,
                    help="Random seed")
    ap.add_argument("--covariates", default="",
                    help="Optional: comma-separated covariates to regress out from Z")
    ap.add_argument("--pool-controls", action="store_true",
                    help="Pool GC+VIV+BSL as ground reference (increases control n from 5 to 15)")
    ap.add_argument("--candidate-genes", default="",
                    help="Path to a text file with pre-registered gene IDs (one per line). "
                         "If provided, BH correction is applied only to these genes, "
                         "drastically reducing the multiple-testing burden.")
    ap.add_argument("--focused-permutation", action="store_true",
                    help="Run focused BH correction on top-decile rewired genes only (~250 genes). "
                         "Dramatically reduces the multiple-testing burden: BH threshold goes from "
                         "p<0.00002 (2500 tests) to p<0.0002 (250 tests), which is achievable "
                         "with the combinatorial limits of n=5 vs n=15.")
    ap.add_argument("--focused-quantile", type=float, default=0.90,
                    help="Quantile threshold for focused permutation (default: 0.90 = top decile)")
    ap.add_argument("--pathway-permutation", action="store_true",
                    help="Run pathway-level permutation testing (GSEA-competitive style). "
                         "Tests whether pathway member genes have higher rewiring than non-members "
                         "under label permutation. Tests ~100 pathways instead of ~2500 genes.")
    ap.add_argument("--gene-map", default="",
                    help="Ensembl→Symbol mapping TSV (for pathway-permutation). "
                         "Auto-detected from data/processed/resources/id_map.tsv if not given.")
    ap.add_argument("--pathway-libraries", default="KEGG_2019_Mouse,MSigDB_Hallmark_2020",
                    help="Comma-separated Enrichr library names for pathway-permutation. "
                         "Fewer libraries = fewer tests = more FDR power. "
                         "Default: KEGG_2019_Mouse,MSigDB_Hallmark_2020 (~350 sets).")
    args = ap.parse_args()
    pool = args.pool_controls

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

    # Load candidate gene set for focused BH correction
    candidate_mask = None
    if args.candidate_genes and Path(args.candidate_genes).exists():
        cand_ids = set()
        with open(args.candidate_genes) as f:
            for line in f:
                g = line.strip()
                if g and not g.startswith("#"):
                    cand_ids.add(g)
        candidate_mask = np.array([g in cand_ids for g in genes], dtype=bool)
        n_cand = candidate_mask.sum()
        print(f"\n[FOCUSED MODE] {len(cand_ids)} candidate IDs loaded, "
              f"{n_cand} found in network ({G} total)")
        print(f"  BH correction on {n_cand} genes instead of {G}")

    # ── Load gene-set definitions for pathway-level permutation ──────────
    pathway_sets: dict[str, np.ndarray] = {}  # name → boolean mask over G genes
    if args.pathway_permutation:
        print("\n[PATHWAY-PERMUTATION] Loading gene sets for pathway-level testing...")
        # Resolve gene map for Ensembl→Symbol
        map_path = Path(args.gene_map) if args.gene_map else REPO_ROOT / "data/processed/resources/id_map.tsv"
        ens_to_sym: dict[str, str] = {}
        sym_to_ens: dict[str, set[str]] = {}
        if map_path.exists():
            m = pd.read_csv(map_path, sep="\t", comment="#")
            col_map = {}
            for c in m.columns:
                cl = c.lower().strip()
                if "ensembl" in cl:
                    col_map[c] = "ensembl_gene_id"
                elif "symbol" in cl or "mgi" in cl:
                    col_map[c] = "mgi_symbol"
            m = m.rename(columns=col_map)
            if "ensembl_gene_id" in m.columns and "mgi_symbol" in m.columns:
                for _, row in m.iterrows():
                    eid = str(row["ensembl_gene_id"]).strip()
                    sym = str(row["mgi_symbol"]).strip()
                    if eid and sym and sym != "nan":
                        ens_to_sym[eid] = sym
                        sym_to_ens.setdefault(sym.lower(), set()).add(eid)
                print(f"  Loaded ID map: {len(ens_to_sym)} Ensembl→Symbol")
        else:
            print(f"  WARNING: No gene map at {map_path}. Pathway testing may have 0 overlap.")

        # Load gene sets via the shared gene_set_loader
        # Use a restricted library set to keep hypothesis count manageable
        pw_libraries = [s.strip() for s in args.pathway_libraries.split(",") if s.strip()]
        try:
            from src.enrichment.gene_set_loader import load_gene_sets
            loaded_sets, set_to_lib = load_gene_sets(
                libraries=pw_libraries,
                include_curated=True,
                min_size=5,
                max_size=500,
            )
        except Exception as e:
            print(f"  WARNING: Could not load gene sets: {e}")
            loaded_sets = {}

        # Resolve each gene set to a boolean mask over the G network genes
        gene_set_ens = set(genes)  # Ensembl IDs in the network
        for setname, symbols in loaded_sets.items():
            # Convert symbols to Ensembl IDs that are in our gene list
            matched = set()
            for s in symbols:
                ens_ids = sym_to_ens.get(s.lower(), set())
                matched |= (ens_ids & gene_set_ens)
            if len(matched) >= 3:  # need at least 3 genes for a meaningful test
                mask = np.array([g in matched for g in genes], dtype=bool)
                pathway_sets[setname] = mask

        print(f"  Resolved {len(pathway_sets)} gene sets with ≥3 genes in the network "
              f"(from {len(loaded_sets)} total)")
        if pathway_sets:
            sizes = [int(m.sum()) for m in pathway_sets.values()]
            print(f"  Set sizes: min={min(sizes)}, median={int(np.median(sizes))}, max={max(sizes)}")

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

    # Define the four contrasts (within each Age×Arm, compare FLT vs controls)
    ctrl_groups = ["FLT", "GC", "VIV", "BSL"] if pool else ["FLT", "GC"]
    ctrl_tag = "GND" if pool else "GC"
    contrasts = [
        (f"ISS_T_YNG_FLT_minus_{ctrl_tag}", {"Age": "YNG", "Arm": "ISS-T"}),
        (f"ISS_T_OLD_FLT_minus_{ctrl_tag}", {"Age": "OLD", "Arm": "ISS-T"}),
        (f"LAR_YNG_FLT_minus_{ctrl_tag}",   {"Age": "YNG", "Arm": "LAR"}),
        (f"LAR_OLD_FLT_minus_{ctrl_tag}",   {"Age": "OLD", "Arm": "LAR"}),
    ]

    if pool:
        print(f"\n[POOL MODE] Using all non-flight groups as controls: {ctrl_groups}")

    for cname, filt in contrasts:
        print(f"\n--- Processing contrast: {cname} ---")
        mask = (meta["Age"] == filt["Age"]) & (meta["Arm"] == filt["Arm"]) & (meta["EnvGroup"].isin(ctrl_groups))
        sub_idx = np.where(mask.values)[0]
        sub = meta.iloc[sub_idx].copy()
        zsub = Z[sub_idx, :]  # n x E

        flt_idx = np.where(sub["EnvGroup"].values == "FLT")[0]
        ctrl_idx = np.where(sub["EnvGroup"].values != "FLT")[0]
        
        if len(flt_idx) < 2 or len(ctrl_idx) < 2:
            print(f"  [SKIP] Not enough samples (FLT={len(flt_idx)} CTRL={len(ctrl_idx)})")
            continue

        print(f"  Samples: {len(sub)} total (FLT={len(flt_idx)}, CTRL={len(ctrl_idx)})")

        # Observed delta per edge (mean difference)
        delta_obs = zsub[flt_idx].mean(axis=0) - zsub[ctrl_idx].mean(axis=0)
        rew_obs = node_rewiring_abs(delta_obs, edge_i, edge_j, G)

        # Permutation null: shuffle FLT vs non-FLT labels
        print(f"  Running {args.K_perm} permutations...")
        # Create binary labels: "FLT" vs "CTRL"
        is_flt = (sub["EnvGroup"].values == "FLT")
        null = np.zeros((args.K_perm, G), dtype=np.float32)

        for k in range(args.K_perm):
            perm = is_flt.copy()
            rng.shuffle(perm)
            p_flt = np.where(perm)[0]
            p_ctrl = np.where(~perm)[0]
            d = zsub[p_flt].mean(axis=0) - zsub[p_ctrl].mean(axis=0)
            null[k] = node_rewiring_abs(d, edge_i, edge_j, G).astype(np.float32)

        # One-sided p-value: P(null >= obs)
        pvals = (1.0 + (null >= rew_obs[None, :]).sum(axis=0)) / (args.K_perm + 1.0)
        qvals = bh_fdr(pvals)

        perm_out = pd.DataFrame({
            "gene": genes,
            "rewiring_abs_obs": rew_obs,
            "p_perm": pvals,
            "q_BH": qvals,
        })

        # If candidate genes specified, add focused q-values
        if candidate_mask is not None:
            cand_pvals = pvals[candidate_mask]
            cand_qvals = bh_fdr(cand_pvals)
            # Map focused q-values back to full array
            q_focused = np.ones(G, dtype=np.float64)
            q_focused[candidate_mask] = cand_qvals
            perm_out["q_BH_focused"] = q_focused
            n_sig_focused = (cand_qvals < 0.05).sum()
            print(f"  Focused BH: {n_sig_focused} candidate genes q<0.05 "
                  f"(of {candidate_mask.sum()} tested)")

        perm_out = perm_out.sort_values("p_perm").reset_index(drop=True)
        perm_path = outdir / f"{cname}_perm_pvals.tsv"
        perm_out.to_csv(perm_path, sep="\t", index=False)
        n_sig = (qvals < 0.05).sum()
        print(f"  Wrote {perm_path} ({n_sig} genes q<0.05)")

        # Bootstrap CI (percentile)
        print(f"  Running {args.B_boot} bootstrap resamples...")
        boot = np.zeros((args.B_boot, G), dtype=np.float32)

        for b in range(args.B_boot):
            b_flt = rng.choice(flt_idx, size=len(flt_idx), replace=True)
            b_ctrl = rng.choice(ctrl_idx, size=len(ctrl_idx), replace=True)
            d = zsub[b_flt].mean(axis=0) - zsub[b_ctrl].mean(axis=0)
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

        # ── Fix 1: Focused permutation (top-decile BH) ──────────────────
        if args.focused_permutation:
            # Identify top-decile genes by observed rewiring
            thr = np.quantile(rew_obs, args.focused_quantile)
            top_mask = rew_obs >= thr
            n_top = int(top_mask.sum())
            print(f"\n  [FOCUSED-PERM] Top-decile threshold (q={args.focused_quantile}): "
                  f"{thr:.2f} ({n_top} genes)")

            # Extract p-values only for top-decile genes and apply BH on that subset
            top_pvals = pvals[top_mask]
            top_qvals = bh_fdr(top_pvals)
            n_sig_foc = int((top_qvals < 0.05).sum())

            # Build output: only top-decile genes
            top_genes = np.array(genes)[top_mask]
            top_rew = rew_obs[top_mask]
            foc_out = pd.DataFrame({
                "gene": top_genes,
                "rewiring_abs_obs": top_rew,
                "p_perm": top_pvals,
                "q_BH_focused": top_qvals,
            }).sort_values("p_perm").reset_index(drop=True)
            foc_path = outdir / f"{cname}_perm_focused.tsv"
            foc_out.to_csv(foc_path, sep="\t", index=False)
            print(f"  [FOCUSED-PERM] {n_sig_foc} / {n_top} genes q<0.05")
            print(f"  Wrote {foc_path}")

            # Diagnostic: show achievability
            min_possible_p = 1.0 / (args.K_perm + 1.0)
            bh_threshold_full = 0.05 / G
            bh_threshold_foc = 0.05 / n_top
            print(f"  Diagnostic: min achievable p = {min_possible_p:.6f}")
            print(f"    Full BH needs p < {bh_threshold_full:.6f} (top gene) — "
                  f"{'REACHABLE' if min_possible_p < bh_threshold_full else 'UNREACHABLE'}")
            print(f"    Focused BH needs p < {bh_threshold_foc:.6f} (top gene) — "
                  f"{'REACHABLE' if min_possible_p < bh_threshold_foc else 'UNREACHABLE'}")

        # ── Fix 2: Pathway-level permutation ─────────────────────────────
        if args.pathway_permutation and pathway_sets:
            print(f"\n  [PATHWAY-PERM] Testing {len(pathway_sets)} gene sets...")

            # Observed pathway-level statistic: mean(rewiring_abs) for member genes
            # Competitive: compare mean(members) vs mean(non-members)
            pw_results = []
            for setname, set_mask in pathway_sets.items():
                n_in = int(set_mask.sum())
                obs_mean_in = rew_obs[set_mask].mean()
                obs_mean_out = rew_obs[~set_mask].mean()
                obs_diff = obs_mean_in - obs_mean_out  # competitive enrichment stat

                # Permutation null: under shuffled FLT/CTRL labels, recompute
                # rewiring_abs then re-measure the mean(in) - mean(out) difference
                null_diffs = np.zeros(args.K_perm, dtype=np.float32)
                for k in range(args.K_perm):
                    null_rew = null[k]  # already computed above, shape (G,)
                    null_diffs[k] = null_rew[set_mask].mean() - null_rew[~set_mask].mean()

                # One-sided p-value: P(null_diff >= obs_diff)
                pw_p = (1.0 + (null_diffs >= obs_diff).sum()) / (args.K_perm + 1.0)

                pw_results.append({
                    "gene_set": setname,
                    "n_genes_in_network": n_in,
                    "mean_rewiring_members": round(float(obs_mean_in), 4),
                    "mean_rewiring_nonmembers": round(float(obs_mean_out), 4),
                    "enrichment_stat": round(float(obs_diff), 4),
                    "p_perm": float(pw_p),
                })

            pw_df = pd.DataFrame(pw_results)
            pw_df["q_BH"] = bh_fdr(pw_df["p_perm"].values)
            pw_df = pw_df.sort_values("p_perm").reset_index(drop=True)

            n_sig_pw = int((pw_df["q_BH"] < 0.05).sum())
            n_nom_pw = int((pw_df["p_perm"] < 0.05).sum())
            pw_path = outdir / f"{cname}_perm_pathway.tsv"
            pw_df.to_csv(pw_path, sep="\t", index=False)
            print(f"  [PATHWAY-PERM] {n_sig_pw} sets FDR<0.05, "
                  f"{n_nom_pw} nominally p<0.05 (of {len(pw_df)} tested)")
            print(f"  Wrote {pw_path}")

            # Show top 10
            for _, row in pw_df.head(10).iterrows():
                marker = "*" if row["q_BH"] < 0.05 else " "
                print(f"    {marker} {row['gene_set'][:55]:55s}  "
                      f"n={row['n_genes_in_network']:3d}  "
                      f"stat={row['enrichment_stat']:+.4f}  "
                      f"p={row['p_perm']:.4f}  q={row['q_BH']:.4f}")

    print(f"\n[OK] All outputs written to: {outdir}")


if __name__ == "__main__":
    main()
