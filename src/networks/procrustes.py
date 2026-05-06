# src/networks/procrustes.py
"""
Phase 3.3: Procrustes alignment + cosine-distance rewiring (multi-seed)

- Loads configured anchors from config/anchor_genes.yaml
- Resolves symbols to Ensembl IDs through id_map.tsv
- Aligns embeddings per seed via orthogonal Procrustes
- Computes rewiring = 1 − cosine_similarity
- Reports seed mean/std + rank variance

Usage:
    python -m src.networks.procrustes
"""
from __future__ import annotations

import os
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from src.common import REPO_ROOT, load_anchor_config, resolve_configured_genes


def cosine_distance_rows(A: np.ndarray, B: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """1 - cosine similarity row-wise."""
    An = np.linalg.norm(A, axis=1) + eps
    Bn = np.linalg.norm(B, axis=1) + eps
    sim = (A * B).sum(axis=1) / (An * Bn)
    sim = np.clip(sim, -1.0, 1.0)
    return 1.0 - sim


def load_configured_anchor_indices(
    genes: list[str],
    anchor_config: Path,
    id_map: Path,
    outdir: Path,
    min_anchors: int | None = None,
    warn_threshold: int | None = None,
) -> np.ndarray:
    """Resolve configured anchors and return their Phase 2 gene indices."""
    cfg, records = load_anchor_config(anchor_config)
    validation = cfg.get("validation", {}) if isinstance(cfg, dict) else {}
    min_anchors = int(min_anchors if min_anchors is not None else validation.get("minimum_anchors", 20))
    warn_threshold = int(warn_threshold if warn_threshold is not None else validation.get("warn_threshold", 50))

    genes_set = set(genes)
    symbols = [r["symbol"] for r in records]
    resolved = resolve_configured_genes(symbols, id_map, panel_genes=genes_set)
    rec_df = pd.DataFrame(records).rename(columns={"symbol": "query"})
    report = resolved.merge(rec_df, on="query", how="left")
    report_path = outdir / "anchors_report.tsv"
    report.to_csv(report_path, sep="\t", index=False)

    present = report[
        (report["status"] == "mapped") &
        (report["ensembl_gene_id"] != "") &
        (report["in_panel"] == True)
    ].drop_duplicates("ensembl_gene_id")

    missing = report[report["in_panel"] != True]
    if not missing.empty:
        missing_path = outdir / "anchors_missing.tsv"
        missing.to_csv(missing_path, sep="\t", index=False)
        print(f"  Missing/unavailable configured anchors: {missing['query'].nunique()} (see {missing_path})")

    if len(present) < min_anchors:
        raise RuntimeError(
            f"Configured Procrustes anchors below minimum: {len(present)} present in Phase 2 panel, "
            f"minimum_anchors={min_anchors}. Rerun Phase 2 with anchor force-inclusion and inspect "
            f"{report_path}."
        )
    if len(present) < warn_threshold:
        print(
            f"  WARNING: only {len(present)} configured anchors present; "
            f"warn_threshold={warn_threshold}."
        )

    gene_to_idx = {g: i for i, g in enumerate(genes)}
    anchor_ids = present["ensembl_gene_id"].tolist()
    anchors = np.array([gene_to_idx[g] for g in anchor_ids], dtype=np.int32)
    (outdir / "anchors.txt").write_text("\n".join(anchor_ids) + "\n")
    return anchors


def orthogonal_procrustes_align(B: np.ndarray, A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Rotate B to best match A using orthogonal Procrustes:
      find R minimizing ||B R - A||_F
      return (B_aligned, R)
    """
    try:
        from scipy.linalg import orthogonal_procrustes
    except Exception as e:
        raise ImportError("scipy is required: pip install scipy") from e

    R, _ = orthogonal_procrustes(B, A)
    return B @ R, R


def load_embedding(emb_dir: Path, condition: str, seed: int) -> np.ndarray:
    p = emb_dir / condition / f"seed_{seed}.npy"
    if not p.exists():
        raise FileNotFoundError(f"Missing embedding: {p}")
    return np.load(p).astype(np.float64)


def main():
    os.chdir(REPO_ROOT)

    ap = argparse.ArgumentParser(description="Phase 3.3: Procrustes alignment + cosine rewiring (multi-seed)")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2")
    ap.add_argument("--emb_dir", default="data/processed/networks/phase3/embeddings")
    ap.add_argument("--outdir", default="data/results/phase3_rewiring")
    ap.add_argument("--anchor_config", default="config/anchor_genes.yaml")
    ap.add_argument("--id_map", default="data/processed/resources/id_map.tsv")
    ap.add_argument("--min_anchors", type=int, default=None)
    ap.add_argument("--warn_anchors", type=int, default=None)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--num_seeds", type=int, default=10)
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    emb_dir = Path(args.emb_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 3.3: Procrustes Alignment + Cosine Rewiring")
    print("=" * 60)

    genes = [g for g in (phase2 / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    G = len(genes)
    print(f"\nGenes: {G}")

    anchor_config = Path(args.anchor_config)
    if not anchor_config.is_absolute():
        anchor_config = REPO_ROOT / anchor_config
    id_map = Path(args.id_map)
    if not id_map.is_absolute():
        id_map = REPO_ROOT / id_map

    anchors = load_configured_anchor_indices(
        genes=genes,
        anchor_config=anchor_config,
        id_map=id_map,
        outdir=outdir,
        min_anchors=args.min_anchors,
        warn_threshold=args.warn_anchors,
    )
    print(f"Configured anchors loaded: {anchors.size} (saved anchors.txt and anchors_report.tsv)")

    # Define the key FLT vs control comparisons
    # Auto-detect: prefer GND (pooled controls) if available, else GC
    gnd_files = sorted(emb_dir.glob("Pred_*_GND"))
    gc_files  = sorted(emb_dir.glob("Pred_*_GC"))
    ctrl_tag = "GND" if gnd_files else "GC"
    print(f"  Control tag: {ctrl_tag} ({'pooled GC+VIV+BSL' if ctrl_tag == 'GND' else 'GC only'})")

    pairs = [
        (f"ISS_T_YNG_FLT_minus_{ctrl_tag}", "Pred_YNG_ISS_T_FLT", f"Pred_YNG_ISS_T_{ctrl_tag}"),
        (f"ISS_T_OLD_FLT_minus_{ctrl_tag}", "Pred_OLD_ISS_T_FLT", f"Pred_OLD_ISS_T_{ctrl_tag}"),
        (f"LAR_YNG_FLT_minus_{ctrl_tag}",   "Pred_YNG_LAR_FLT",   f"Pred_YNG_LAR_{ctrl_tag}"),
        (f"LAR_OLD_FLT_minus_{ctrl_tag}",   "Pred_OLD_LAR_FLT",   f"Pred_OLD_LAR_{ctrl_tag}"),
    ]

    seeds = [args.seed0 + k for k in range(args.num_seeds)]
    print(f"Seeds: {seeds}")

    summary_rows = []

    for contrast_name, cond_FLT, cond_GC in pairs:
        per_seed_rewiring = []

        print(f"\n=== {contrast_name}: {cond_FLT} vs {cond_GC} ===")
        
        for sd in seeds:
            # Load embeddings
            emb_FLT = load_embedding(emb_dir, cond_FLT, sd)
            emb_GC = load_embedding(emb_dir, cond_GC, sd)

            if emb_FLT.shape[0] != G or emb_GC.shape[0] != G:
                raise ValueError(f"Embedding size mismatch: expected {G} rows")

            # Align GC to FLT using anchors
            FLT_anchor = emb_FLT[anchors, :]
            GC_anchor = emb_GC[anchors, :]
            
            # Compute rotation on anchors, apply to full GC
            _, R = orthogonal_procrustes_align(GC_anchor, FLT_anchor)
            GC_aligned = emb_GC @ R

            # Rewiring = 1 - cosine_sim (FLT vs aligned GC)
            rew = cosine_distance_rows(emb_FLT, GC_aligned)
            
            df = pd.DataFrame({"gene": genes, "rewiring": rew.astype(float)})
            df = df.sort_values("rewiring", ascending=False, kind="mergesort").reset_index(drop=True)

            seed_path = outdir / f"{contrast_name}_seed_{sd}.tsv"
            df.to_csv(seed_path, sep="\t", index=False)
            per_seed_rewiring.append(rew.astype(np.float64))

            print(f"  seed {sd}: saved (top gene: {df.iloc[0]['gene']} = {df.iloc[0]['rewiring']:.4f})")

        # Aggregate across seeds
        stack = np.stack(per_seed_rewiring, axis=0)  # (S, G)
        mean_rew = stack.mean(axis=0)
        std_rew = stack.std(axis=0)

        # Rank variance across seeds
        # For each seed, rank genes by rewiring (0=highest), then compute std of ranks
        ranks = np.argsort(np.argsort(-stack, axis=1), axis=1).astype(np.float64)  # 0=highest
        rank_std = ranks.std(axis=0)

        agg = pd.DataFrame(
            {
                "gene": genes,
                "rewiring_mean": mean_rew.astype(float),
                "rewiring_std": std_rew.astype(float),
                "rank_std": rank_std.astype(float),
            }
        ).sort_values("rewiring_mean", ascending=False, kind="mergesort").reset_index(drop=True)

        agg_path = outdir / f"{contrast_name}_rewiring_agg.tsv"
        agg.to_csv(agg_path, sep="\t", index=False)
        print(f"  Aggregate saved: {agg_path}")
        
        # Summary stats for this contrast
        summary_rows.append({
            "contrast": contrast_name,
            "n_genes": G,
            "n_seeds": len(seeds),
            "mean_rewiring": float(mean_rew.mean()),
            "median_rank_std": float(np.median(rank_std)),
            "top_gene": agg.iloc[0]["gene"],
            "top_rewiring": float(agg.iloc[0]["rewiring_mean"]),
        })

    # Save summary
    summary_df = pd.DataFrame(summary_rows)
    summary_path = outdir / "procrustes_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print(f"Done. Summary: {summary_path}")
    print("=" * 60)
    print(f"\nOutputs in: {outdir}")


if __name__ == "__main__":
    main()
