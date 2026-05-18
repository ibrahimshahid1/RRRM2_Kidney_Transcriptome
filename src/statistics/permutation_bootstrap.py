# src/statistics/permutation_bootstrap.py
"""
Phase 6: fast uncertainty tests for edge-sum node rewiring.

This module works in LIONESS edge-weight space.  For each FLT-control contrast
it computes a per-edge mean difference, then a per-gene statistic equal to the
sum of absolute incident edge differences.  These p-values are therefore
edge-sum node-rewiring tests. They are not direct inference for the Phase 3
node2vec/Procrustes cosine-distance rewiring statistic.

Focused testing is limited to statistically valid pre-specified modes:
  * --candidate-genes: BH only over an external, pre-registered candidate file.
  * --hierarchical-fdr: Benjamini-Bogomolov-style two-stage FDR over configured
    gene families/pathways. Overlapping families are additionally reported with
    a BY column because dependence is not fully eliminated by this procedure.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import (
    REPO_ROOT,
    bh_fdr,
    by_fdr,
    find_sample_col,
    id_map_lookup,
    normalize_labels,
)

DEFAULT_K_PERM = 5000
DEFAULT_B_BOOT = 2000


def edge_sum_node_rewiring_abs(
    delta: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
    n_genes: int,
) -> np.ndarray:
    """Sum |edge delta| over incident edges for each gene."""
    abs_sum = np.zeros(n_genes, dtype=np.float64)
    dabs = np.abs(delta).astype(np.float64)
    np.add.at(abs_sum, edge_i, dabs)
    np.add.at(abs_sum, edge_j, dabs)
    return abs_sum


def load_candidate_mask(path: str | Path, genes: list[str]) -> tuple[np.ndarray, set[str]]:
    """Load a pre-registered candidate gene list; no observed statistic is used."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Candidate-gene file not found: {path}")
    candidates = {
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    mask = np.array([g in candidates for g in genes], dtype=bool)
    if mask.sum() == 0:
        raise ValueError(f"No candidate genes from {path} are present in the Phase 2 panel")
    return mask, candidates


def _flatten_gene_sets_yaml(yaml_path: Path) -> dict[str, list[str]]:
    try:
        import yaml
    except Exception as exc:  # pragma: no cover
        raise ImportError("PyYAML is required for hierarchical FDR gene sets") from exc
    cfg = yaml.safe_load(yaml_path.read_text()) or {}
    sets: dict[str, list[str]] = {}
    for key, val in cfg.items():
        if not isinstance(val, dict) or "genes" not in val:
            continue
        genes: list[str] = []
        for item in val.get("genes", []):
            if isinstance(item, str):
                genes.append(item.strip())
            elif isinstance(item, list):
                genes.extend(str(x).strip() for x in item if str(x).strip())
            elif isinstance(item, dict):
                for sublist in item.values():
                    if isinstance(sublist, list):
                        genes.extend(str(x).strip() for x in sublist if str(x).strip())
        if genes:
            sets[f"curated::{key}"] = genes
    return sets


def load_hierarchical_gene_sets(
    gene_set_yaml: str | Path,
    id_map: str | Path,
    genes: list[str],
    min_size: int = 3,
) -> dict[str, np.ndarray]:
    """Resolve pre-specified YAML gene sets to masks over Phase 2 genes."""
    gene_set_yaml = Path(gene_set_yaml)
    if not gene_set_yaml.exists():
        raise FileNotFoundError(f"Hierarchical gene-set YAML not found: {gene_set_yaml}")
    _, symbol_to_ens = id_map_lookup(id_map)
    gene_universe = set(genes)
    raw_sets = _flatten_gene_sets_yaml(gene_set_yaml)

    masks: dict[str, np.ndarray] = {}
    for set_name, symbols in raw_sets.items():
        matched: set[str] = set()
        for symbol in symbols:
            if symbol in gene_universe:
                matched.add(symbol)
            matched |= (symbol_to_ens.get(symbol.lower(), set()) & gene_universe)
        if len(matched) >= min_size:
            masks[set_name] = np.array([g in matched for g in genes], dtype=bool)
    if not masks:
        raise ValueError(
            f"No hierarchical gene sets from {gene_set_yaml} overlap the Phase 2 panel "
            f"with min_size={min_size}"
        )
    return masks


def hierarchical_fdr_tables(
    rew_obs: np.ndarray,
    null: np.ndarray,
    p_gene: np.ndarray,
    genes: list[str],
    gene_sets: dict[str, np.ndarray],
    alpha: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Two-stage family/gene testing with conservative overlap diagnostics."""
    family_rows = []
    for set_name, mask in gene_sets.items():
        obs = rew_obs[mask].mean() - rew_obs[~mask].mean()
        null_stat = null[:, mask].mean(axis=1) - null[:, ~mask].mean(axis=1)
        p_family = (1.0 + (null_stat >= obs).sum()) / (null.shape[0] + 1.0)
        family_rows.append({
            "gene_set": set_name,
            "n_genes_in_panel": int(mask.sum()),
            "edge_sum_enrichment_stat": float(obs),
            "p_family_edge_sum": float(p_family),
        })

    fam = pd.DataFrame(family_rows)
    fam["q_family_BH"] = bh_fdr(fam["p_family_edge_sum"].values)
    fam["q_family_BY_overlap_conservative"] = by_fdr(fam["p_family_edge_sum"].values)
    selected = set(fam.loc[fam["q_family_BH"] < alpha, "gene_set"])

    n_families = len(fam)
    n_selected = max(len(selected), 1)
    gene_rows = []
    for set_name, mask in gene_sets.items():
        set_gene_idx = np.where(mask)[0]
        q_within = bh_fdr(p_gene[set_gene_idx])
        family_selected = set_name in selected
        if family_selected:
            q_bb = np.clip(q_within * n_families / n_selected, 0, 1)
        else:
            q_bb = np.full_like(q_within, np.nan, dtype=float)
        for local_i, gene_i in enumerate(set_gene_idx):
            gene_rows.append({
                "gene_set": set_name,
                "family_selected_q05": family_selected,
                "confirmatory_gene_tested": family_selected,
                "gene": genes[gene_i],
                "edge_sum_node_rewiring_obs": float(rew_obs[gene_i]),
                "p_perm_edge_sum": float(p_gene[gene_i]),
                "q_within_family_BH": float(q_within[local_i]),
                "q_BB_two_stage": float(q_bb[local_i]),
            })

    gene_df = pd.DataFrame(gene_rows)
    if not gene_df.empty:
        gene_df["q_BY_all_family_gene_rows_overlap_conservative"] = by_fdr(
            gene_df["p_perm_edge_sum"].values
        )
        gene_df = gene_df.sort_values(
            ["confirmatory_gene_tested", "q_BB_two_stage", "p_perm_edge_sum"],
            ascending=[False, True, True],
            na_position="last",
        ).reset_index(drop=True)
    fam = fam.sort_values("p_family_edge_sum").reset_index(drop=True)
    return fam, gene_df


def main() -> None:
    ap = argparse.ArgumentParser(description="Permutation/bootstrap tests for edge-sum node rewiring")
    ap.add_argument("--phase2_dir", default=str(REPO_ROOT / "data/processed/networks/phase2"))
    ap.add_argument("--meta", default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"))
    ap.add_argument("--edge-weights", "--z", dest="edge_weights", default="lioness_edges.npy")
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/phase6_uncertainty"))
    ap.add_argument("--K_perm", type=int, default=DEFAULT_K_PERM)
    ap.add_argument("--B_boot", type=int, default=DEFAULT_B_BOOT)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--covariates", default="")
    ap.add_argument("--pool-controls", action="store_true")
    ap.add_argument("--candidate-genes", default="",
                    help="Pre-registered Ensembl-ID file. BH is applied only to listed genes.")
    ap.add_argument("--hierarchical-fdr", action="store_true",
                    help="Run two-stage FDR over pre-specified gene families/pathways.")
    ap.add_argument("--hierarchical-gene-sets", default=str(REPO_ROOT / "config/gene_sets.yaml"))
    ap.add_argument("--hierarchical-alpha", type=float, default=0.05)
    ap.add_argument("--gene-map", default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"))
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    phase2 = Path(args.phase2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes_path = phase2 / "phase2_genes.txt"
    if not genes_path.exists():
        raise FileNotFoundError(f"Missing genes file: {genes_path}")
    genes = [g.strip() for g in genes_path.read_text().splitlines() if g.strip()]
    n_genes = len(genes)
    print(f"Loaded {n_genes} genes")

    candidate_mask = None
    if args.candidate_genes:
        candidate_mask, candidates = load_candidate_mask(args.candidate_genes, genes)
        print(
            f"\n[CANDIDATE MODE] {len(candidates)} pre-registered IDs loaded; "
            f"{int(candidate_mask.sum())} present in Phase 2 panel."
        )

    hierarchical_sets: dict[str, np.ndarray] = {}
    if args.hierarchical_fdr:
        print("\n[HIERARCHICAL FDR] Loading pre-specified gene sets...")
        hierarchical_sets = load_hierarchical_gene_sets(
            args.hierarchical_gene_sets,
            args.gene_map,
            genes,
            min_size=3,
        )
        memberships = np.zeros(n_genes, dtype=int)
        for mask in hierarchical_sets.values():
            memberships += mask.astype(int)
        overlap_count = int((memberships > 1).sum())
        print(f"  Loaded {len(hierarchical_sets)} gene sets; {overlap_count} genes are in >1 set")
        if overlap_count:
            print("  Overlap note: reporting BY-adjusted conservative columns in addition to two-stage q-values.")

    edge_i = np.load(phase2 / "edge_i.npy")
    edge_j = np.load(phase2 / "edge_j.npy")
    print(f"Loaded {len(edge_i)} edges")

    z_path = phase2 / args.edge_weights
    if not z_path.exists() and args.edge_weights == "lioness_edges.npy":
        legacy = phase2 / "lioness_z_edges.npy"
        if legacy.exists():
            z_path = legacy
            print("  WARNING: using legacy lioness_z_edges.npy. Prefer lioness_edges.npy for corrected defaults.")
    if not z_path.exists():
        raise FileNotFoundError(f"Missing Z matrix: {z_path}")
    z = np.load(z_path)
    n_samples = z.shape[0]
    print(f"Loaded Z matrix: {z.shape}")

    samples_path = phase2 / "lioness_samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(f"Missing samples file: {samples_path}")
    samples = [s.strip() for s in samples_path.read_text().splitlines() if s.strip()]
    if len(samples) != n_samples:
        raise ValueError(f"lioness_samples.txt length ({len(samples)}) != Z rows ({n_samples})")

    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)
    meta = normalize_labels(meta)
    meta = meta.loc[samples].copy()
    print(f"Aligned metadata: {len(meta)} samples")

    covs = [c.strip() for c in args.covariates.split(",") if c.strip() and c.strip() in meta.columns]
    if covs:
        x_parts = [np.ones((n_samples, 1), dtype=float)]
        for c in covs:
            v = meta[c]
            if pd.api.types.is_numeric_dtype(v):
                vv = (v.values.astype(float) - np.nanmean(v.values.astype(float))) / (
                    np.nanstd(v.values.astype(float)) + 1e-8
                )
                x_parts.append(vv.reshape(-1, 1))
            else:
                dummies = pd.get_dummies(v.astype(str), drop_first=True)
                if dummies.shape[1] > 0:
                    x_parts.append(dummies.values.astype(float))
        x = np.concatenate(x_parts, axis=1)
        beta, *_ = np.linalg.lstsq(x, z, rcond=None)
        z = z - x @ beta
        print(f"[OK] Residualized Z on covariates: {covs}")

    pool = args.pool_controls
    ctrl_groups = ["FLT", "GC", "VIV", "BSL"] if pool else ["FLT", "GC"]
    ctrl_tag = "GND" if pool else "GC"
    contrasts = [
        (f"ISS_T_YNG_FLT_minus_{ctrl_tag}", {"Age": "YNG", "Arm": "ISS-T"}),
        (f"ISS_T_OLD_FLT_minus_{ctrl_tag}", {"Age": "OLD", "Arm": "ISS-T"}),
        (f"LAR_YNG_FLT_minus_{ctrl_tag}", {"Age": "YNG", "Arm": "LAR"}),
        (f"LAR_OLD_FLT_minus_{ctrl_tag}", {"Age": "OLD", "Arm": "LAR"}),
    ]

    for cname, filt in contrasts:
        print(f"\n--- Processing contrast: {cname} ---")
        mask = (
            (meta["Age"] == filt["Age"]) &
            (meta["Arm"] == filt["Arm"]) &
            (meta["EnvGroup"].isin(ctrl_groups))
        )
        sub_idx = np.where(mask.values)[0]
        sub = meta.iloc[sub_idx].copy()
        zsub = z[sub_idx, :]

        flt_idx = np.where(sub["EnvGroup"].values == "FLT")[0]
        ctrl_idx = np.where(sub["EnvGroup"].values != "FLT")[0]
        if len(flt_idx) < 2 or len(ctrl_idx) < 2:
            print(f"  [SKIP] Not enough samples (FLT={len(flt_idx)} CTRL={len(ctrl_idx)})")
            continue

        print(f"  Samples: {len(sub)} total (FLT={len(flt_idx)}, CTRL={len(ctrl_idx)})")
        delta_obs = zsub[flt_idx].mean(axis=0) - zsub[ctrl_idx].mean(axis=0)
        rew_obs = edge_sum_node_rewiring_abs(delta_obs, edge_i, edge_j, n_genes)

        print(f"  Running {args.K_perm} edge-sum permutations...")
        is_flt = sub["EnvGroup"].values == "FLT"
        null = np.zeros((args.K_perm, n_genes), dtype=np.float32)
        for k in range(args.K_perm):
            perm = is_flt.copy()
            rng.shuffle(perm)
            p_flt = np.where(perm)[0]
            p_ctrl = np.where(~perm)[0]
            d = zsub[p_flt].mean(axis=0) - zsub[p_ctrl].mean(axis=0)
            null[k] = edge_sum_node_rewiring_abs(d, edge_i, edge_j, n_genes).astype(np.float32)

        pvals = (1.0 + (null >= rew_obs[None, :]).sum(axis=0)) / (args.K_perm + 1.0)
        qvals = bh_fdr(pvals)

        perm_out = pd.DataFrame({
            "gene": genes,
            "edge_sum_node_rewiring_obs": rew_obs,
            "p_perm_edge_sum": pvals,
            "q_BH_edge_sum": qvals,
        })

        if candidate_mask is not None:
            cand_q = bh_fdr(pvals[candidate_mask])
            q_candidate = np.ones(n_genes, dtype=float)
            q_candidate[candidate_mask] = cand_q
            perm_out["tested_in_preregistered_candidate_set"] = candidate_mask
            perm_out["q_BH_candidate"] = q_candidate
            cand_out = perm_out[perm_out["tested_in_preregistered_candidate_set"]].copy()
            cand_out = cand_out.sort_values("p_perm_edge_sum").reset_index(drop=True)
            cand_path = outdir / f"{cname}_perm_candidate_genes.tsv"
            cand_out.to_csv(cand_path, sep="\t", index=False)
            print(f"  Candidate-only BH: {(cand_q < 0.05).sum()} / {candidate_mask.sum()} q<0.05")
            print(f"  Wrote {cand_path}")

        perm_out = perm_out.sort_values("p_perm_edge_sum").reset_index(drop=True)
        perm_path = outdir / f"{cname}_perm_edge_sum_pvals.tsv"
        perm_out.to_csv(perm_path, sep="\t", index=False)

        # Backward-compatible filename, but with corrected columns.
        legacy_path = outdir / f"{cname}_perm_pvals.tsv"
        perm_out.to_csv(legacy_path, sep="\t", index=False)
        print(f"  Wrote {perm_path} ({int((qvals < 0.05).sum())} genes q<0.05)")

        print(f"  Running {args.B_boot} bootstrap resamples...")
        boot = np.zeros((args.B_boot, n_genes), dtype=np.float32)
        for b in range(args.B_boot):
            b_flt = rng.choice(flt_idx, size=len(flt_idx), replace=True)
            b_ctrl = rng.choice(ctrl_idx, size=len(ctrl_idx), replace=True)
            d = zsub[b_flt].mean(axis=0) - zsub[b_ctrl].mean(axis=0)
            boot[b] = edge_sum_node_rewiring_abs(d, edge_i, edge_j, n_genes).astype(np.float32)

        lo = np.quantile(boot, 0.025, axis=0)
        hi = np.quantile(boot, 0.975, axis=0)
        boot_out = pd.DataFrame({
            "gene": genes,
            "edge_sum_node_rewiring_obs": rew_obs,
            "edge_sum_ci_2p5": lo,
            "edge_sum_ci_97p5": hi,
        }).sort_values("edge_sum_node_rewiring_obs", ascending=False).reset_index(drop=True)
        boot_path = outdir / f"{cname}_bootstrap_edge_sum_ci.tsv"
        boot_out.to_csv(boot_path, sep="\t", index=False)
        legacy_boot_path = outdir / f"{cname}_bootstrap_ci.tsv"
        boot_out.to_csv(legacy_boot_path, sep="\t", index=False)
        print(f"  Wrote {boot_path}")

        if hierarchical_sets:
            fam_df, gene_df = hierarchical_fdr_tables(
                rew_obs=rew_obs,
                null=null,
                p_gene=pvals,
                genes=genes,
                gene_sets=hierarchical_sets,
                alpha=args.hierarchical_alpha,
            )
            fam_path = outdir / f"{cname}_hierarchical_fdr_families.tsv"
            gene_path = outdir / f"{cname}_hierarchical_fdr_genes.tsv"
            fam_df.to_csv(fam_path, sep="\t", index=False)
            gene_df.to_csv(gene_path, sep="\t", index=False)
            print(
                f"  Hierarchical FDR: {(fam_df['q_family_BH'] < args.hierarchical_alpha).sum()} "
                f"families selected at q<{args.hierarchical_alpha}"
            )
            print(f"  Wrote {fam_path}")
            print(f"  Wrote {gene_path}")

    print(f"\n[OK] All edge-sum node-rewiring outputs written to: {outdir}")


if __name__ == "__main__":
    main()
