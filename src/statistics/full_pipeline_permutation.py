# src/statistics/full_pipeline_permutation.py
"""
Appendix-tier full-pipeline permutation for selected top-hit genes.

This is deliberately separate from src.statistics.permutation_bootstrap. The
fast Phase 6 test calibrates edge-sum node rewiring; this driver is reserved for
the expensive statistic that matches Phase 3 node2vec/Procrustes cosine-distance
rewiring. It permutes labels, reruns the requested pipeline commands into
isolated directories, and records a manifest for auditable execution.

By default the module writes a dry-run manifest only. Use --execute after a
pre-registered top-hit file and command template have been reviewed.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT

DEFAULT_COMMAND_TEMPLATE = (
    "{python} -m src.networks.shared_topology "
    "--rtech={rtech} --meta={meta} --outdir={phase2_dir} "
    "--max_genes={max_genes} --topk={topk} --id_map={id_map} "
    "--anchor_config={anchor_config} --biotype_filter=protein_coding && "
    "{python} -m src.networks.lioness "
    "--rtech={rtech} --meta={meta} --phase2_dir={phase2_dir} && "
    "{python} -m src.networks.edge_regression "
    "--meta={meta} --phase2_dir={phase2_dir} --outdir={reg_dir} && "
    "{python} -m src.networks.embeddings "
    "--phase2_dir={phase2_dir} --reg_dir={reg_dir} --outdir={emb_dir} "
    "--num_seeds={num_seeds} --seed0={seed0} && "
    "{python} -m src.networks.procrustes "
    "--phase2_dir={phase2_dir} --emb_dir={emb_dir} --outdir={rewiring_dir} "
    "--anchor_config={anchor_config} --id_map={id_map} "
    "--num_seeds={num_seeds} --seed0={seed0}"
)


def load_top_hits(path: Path, gene_col: str = "gene", max_hits: int | None = None) -> list[str]:
    if not path.exists():
        raise FileNotFoundError(f"Top-hit file not found: {path}")
    if path.suffix.lower() in {".tsv", ".csv"}:
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        df = pd.read_csv(path, sep=sep)
        if gene_col not in df.columns:
            raise ValueError(f"{path} missing gene column '{gene_col}'")
        genes = [str(g) for g in df[gene_col].dropna().drop_duplicates()]
    else:
        genes = [
            line.strip()
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
    return genes[:max_hits] if max_hits else genes


def read_metadata(path: Path) -> pd.DataFrame:
    compression = "gzip" if path.suffix == ".gz" else None
    return pd.read_csv(path, sep="\t", compression=compression)


def find_sample_col(meta: pd.DataFrame) -> str:
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            return col
    return meta.columns[0]


def permute_flight_labels_within_strata(
    meta: pd.DataFrame,
    seed: int,
    strata_cols: list[str],
    env_col: str = "EnvGroup",
) -> pd.DataFrame:
    """Permute FLT/GC labels within strata while leaving other controls unchanged."""
    missing = [c for c in strata_cols + [env_col] if c not in meta.columns]
    if missing:
        raise ValueError(f"Metadata missing required columns for permutation: {missing}")

    rng = np.random.default_rng(seed)
    out = meta.copy()
    env = out[env_col].astype(str).copy()
    strata = out[strata_cols].astype(str).agg("|".join, axis=1)
    for stratum in strata.unique():
        idx = np.where((strata == stratum).values & env.isin(["FLT", "GC"]).values)[0]
        labels = env.iloc[idx].to_numpy(copy=True)
        if len(set(labels)) < 2:
            continue
        rng.shuffle(labels)
        env.iloc[idx] = labels
    out[env_col] = env
    return out


def write_permuted_metadata(meta: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    meta.to_csv(path, sep="\t", index=False, compression=compression)


def shell_quote(value: str | Path) -> str:
    return shlex.quote(str(value))


def validate_command_template(command_template: str) -> None:
    if "{meta}" not in command_template:
        raise ValueError("Full-pipeline permutation command_template must include the {meta} placeholder")
    if "--dry-run" in command_template:
        raise ValueError("Full-pipeline permutation command_template must not include --dry-run")


def write_permutation_manifest(
    outdir: Path,
    top_hits: list[str],
    k_perm: int,
    seed: int,
    command_template: str,
    meta_path: Path,
    rtech_path: Path,
    id_map_path: Path,
    anchor_config: Path,
    max_genes: int,
    topk: int,
    num_seeds: int,
    seed0: int,
    strata_cols: list[str],
) -> pd.DataFrame:
    validate_command_template(command_template)
    meta = read_metadata(meta_path)
    rng = np.random.default_rng(seed)
    rows = []
    for perm_id in range(k_perm):
        perm_seed = int(rng.integers(0, np.iinfo(np.int32).max))
        perm_dir = outdir / f"perm_{perm_id:04d}"
        perm_phase2 = perm_dir / "networks"
        perm_reg = perm_phase2 / "regression"
        perm_emb = perm_dir / "phase3_embeddings"
        perm_rewiring = perm_dir / "phase3_rewiring"
        perm_meta = perm_dir / "meta_permuted.tsv.gz"
        permuted = permute_flight_labels_within_strata(meta, perm_seed, strata_cols=strata_cols)
        write_permuted_metadata(permuted, perm_meta)
        command = command_template.format(
            perm_id=perm_id,
            seed=perm_seed,
            outdir=shell_quote(perm_dir),
            phase2_dir=shell_quote(perm_phase2),
            reg_dir=shell_quote(perm_reg),
            emb_dir=shell_quote(perm_emb),
            rewiring_dir=shell_quote(perm_rewiring),
            meta=shell_quote(perm_meta),
            rtech=shell_quote(rtech_path),
            id_map=shell_quote(id_map_path),
            anchor_config=shell_quote(anchor_config),
            repo_root=shell_quote(REPO_ROOT),
            python=shell_quote(sys.executable),
            max_genes=max_genes,
            topk=topk,
            num_seeds=num_seeds,
            seed0=seed0,
        )
        rows.append({
            "perm_id": perm_id,
            "seed": perm_seed,
            "outdir": str(perm_dir),
            "meta_permuted": str(perm_meta),
            "phase2_dir": str(perm_phase2),
            "rewiring_dir": str(perm_rewiring),
            "command": command,
        })
    manifest = pd.DataFrame(rows)
    manifest.to_csv(outdir / "full_pipeline_permutation_manifest.tsv", sep="\t", index=False)
    (outdir / "top_hits.txt").write_text("\n".join(top_hits) + "\n")
    metadata = {
        "statistic": "Phase 3 node2vec/Procrustes cosine-distance rewiring",
        "scope": "appendix-tier selected top hits only",
        "n_top_hits": len(top_hits),
        "k_perm": k_perm,
        "seed": seed,
        "meta_path": str(meta_path),
        "rtech_path": str(rtech_path),
        "permutation": f"FLT/GC labels shuffled within {' x '.join(strata_cols)} strata",
        "requires_execute_flag": True,
    }
    (outdir / "full_pipeline_permutation_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return manifest


def collect_null_statistics(
    manifest: pd.DataFrame,
    top_hits: list[str],
    contrast: str,
    outdir: Path,
) -> pd.DataFrame:
    rows = []
    for _, row in manifest.iterrows():
        stats_path = Path(row["rewiring_dir"]) / f"{contrast}_rewiring_agg.tsv"
        if not stats_path.exists():
            rows.append({
                "perm_id": row["perm_id"],
                "gene": "",
                "contrast": contrast,
                "rewiring_mean": np.nan,
                "status": "missing_rewiring_output",
            })
            continue
        df = pd.read_csv(stats_path, sep="\t")
        if "gene" not in df.columns or "rewiring_mean" not in df.columns:
            raise ValueError(f"{stats_path} must contain gene and rewiring_mean columns")
        sub = df[df["gene"].isin(top_hits)][["gene", "rewiring_mean"]].copy()
        found = set(sub["gene"])
        for missing_gene in sorted(set(top_hits) - found):
            sub.loc[len(sub)] = {"gene": missing_gene, "rewiring_mean": np.nan}
        sub["perm_id"] = row["perm_id"]
        sub["contrast"] = contrast
        sub["status"] = np.where(sub["rewiring_mean"].isna(), "gene_missing", "ok")
        rows.extend(sub[["perm_id", "contrast", "gene", "rewiring_mean", "status"]].to_dict("records"))

    result = pd.DataFrame(rows)
    result.to_csv(outdir / "full_pipeline_permutation_null_stats.tsv", sep="\t", index=False)
    return result


def execute_manifest(manifest: pd.DataFrame, stop_on_error: bool = True) -> int:
    failures = 0
    for _, row in manifest.iterrows():
        outdir = Path(row["outdir"])
        outdir.mkdir(parents=True, exist_ok=True)
        cmd = str(row["command"])
        print(f"[RUN] {cmd}")
        result = subprocess.run(cmd, shell=True, cwd=REPO_ROOT)
        if result.returncode != 0:
            failures += 1
            print(f"[ERROR] permutation {row['perm_id']} failed with status {result.returncode}")
            if stop_on_error:
                return failures
    return failures


def main() -> None:
    ap = argparse.ArgumentParser(description="Appendix-tier full-pipeline cosine-distance permutation")
    ap.add_argument("--top_hits", required=True,
                    help="Pre-selected top-hit list/table. Must be fixed before permutations are run.")
    ap.add_argument("--gene_col", default="gene")
    ap.add_argument("--max_hits", type=int, default=None)
    ap.add_argument("--K_perm", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/full_pipeline_permutation"))
    ap.add_argument("--meta", default=str(REPO_ROOT / "data/processed/phase1_residuals/meta_phase1.tsv.gz"))
    ap.add_argument("--rtech", default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"))
    ap.add_argument("--id_map", default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"))
    ap.add_argument("--anchor_config", default=str(REPO_ROOT / "config/anchor_genes.yaml"))
    ap.add_argument("--strata_cols", default="Age,Arm")
    ap.add_argument("--contrast", default="ISS_T_YNG_FLT_minus_GC")
    ap.add_argument("--max_genes", type=int, default=2500)
    ap.add_argument("--topk", type=int, default=80)
    ap.add_argument("--num_seeds", type=int, default=10)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--command_template", default=DEFAULT_COMMAND_TEMPLATE)
    ap.add_argument("--execute", action="store_true",
                    help="Actually run command_template for each permutation.")
    ap.add_argument("--continue_on_error", action="store_true")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    top_hits = load_top_hits(Path(args.top_hits), gene_col=args.gene_col, max_hits=args.max_hits)
    if not top_hits:
        raise ValueError("No top-hit genes loaded; refusing to create a full-pipeline permutation run")

    strata_cols = [c.strip() for c in args.strata_cols.split(",") if c.strip()]
    if not strata_cols:
        raise ValueError("--strata_cols must contain at least one metadata column")

    manifest = write_permutation_manifest(
        outdir=outdir,
        top_hits=top_hits,
        k_perm=args.K_perm,
        seed=args.seed,
        command_template=args.command_template,
        meta_path=Path(args.meta),
        rtech_path=Path(args.rtech),
        id_map_path=Path(args.id_map),
        anchor_config=Path(args.anchor_config),
        max_genes=args.max_genes,
        topk=args.topk,
        num_seeds=args.num_seeds,
        seed0=args.seed0,
        strata_cols=strata_cols,
    )
    print(f"[OK] Wrote appendix-tier full-pipeline permutation manifest to {outdir}")

    if args.execute:
        failures = execute_manifest(manifest, stop_on_error=not args.continue_on_error)
        if failures:
            raise SystemExit(failures)
        collect_null_statistics(manifest, top_hits, args.contrast, outdir)
    else:
        print("[DRY] No commands executed. Review the manifest, then rerun with --execute.")


if __name__ == "__main__":
    main()
