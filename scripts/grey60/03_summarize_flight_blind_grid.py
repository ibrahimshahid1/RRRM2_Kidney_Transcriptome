#!/usr/bin/env python3
"""Summarize flight-blind WGCNA recovery and projected Grey60 effects."""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import yaml

from src.grey60.adversarial import mean_z_score, pooled_iss_effect


REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else REPO_ROOT / p


def bh(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    order = np.argsort(x)
    ranked = x[order]
    adjusted = ranked * len(x) / np.arange(1, len(x) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def blocked_iss_permutation_p(
    score: pd.Series,
    meta: pd.DataFrame,
    *,
    n_perm: int,
    seed: int,
) -> float:
    observed = pooled_iss_effect(score, meta)
    rng = np.random.default_rng(seed)
    combos = np.array(list(combinations(range(10), 5)), dtype=int)
    permuted_differences = []
    for age in ("YNG", "OLD"):
        idx = meta.index[(meta["Arm"] == "ISS-T") & (meta["Age"] == age)]
        values = score.loc[idx].to_numpy(dtype=float)
        chosen = combos[rng.integers(0, len(combos), size=n_perm)]
        labels = np.zeros((n_perm, 10), dtype=float)
        labels[np.arange(n_perm)[:, None], chosen] = 1.0
        permuted_differences.append(
            labels @ values / 5.0 - (1.0 - labels) @ values / 5.0
        )
    effects = (permuted_differences[0] + permuted_differences[1]) / 2.0
    exceed = int(np.count_nonzero(np.abs(effects) >= abs(observed)))
    return float((exceed + 1) / (n_perm + 1))


def run(args: argparse.Namespace) -> None:
    cfg = yaml.safe_load(resolve(args.config).read_text())
    inputs = {k: resolve(v) for k, v in cfg["inputs"].items()}
    root = resolve(args.outdir or cfg["output_dir"])
    grid_root = root / "flight_blind_wgcna"
    outdir = root / "internal" / "compactness"
    outdir.mkdir(parents=True, exist_ok=True)

    genes = {
        g.strip()
        for g in inputs["grey60_genes"].read_text().splitlines()
        if g.strip()
    }
    rtech = pd.read_csv(inputs["rrrm2_rtech"], sep="\t", index_col=0)
    meta = pd.read_csv(inputs["rrrm2_meta"], sep="\t")
    sample_col = next(
        c
        for c in ("Sample Name (raw_counts_colname)", "Sample Name", "sample")
        if c in meta
    )
    meta = meta.set_index(sample_col, drop=False)
    meta["EnvGroup"] = meta["EnvGroup"].replace({"HGC": "GC", "VGC": "VIV"})
    fg_samples = [
        s for s in rtech.columns if s in meta.index and meta.loc[s, "EnvGroup"] in ("FLT", "GC")
    ]
    meta_fg = meta.loc[fg_samples]
    expr_fg = rtech.loc[:, fg_samples]

    rows = []
    detail_rows = []
    files = sorted((grid_root / "grid_assignments").glob("*.tsv.gz"))
    if len(files) != 27:
        raise RuntimeError(f"Expected 27 WGCNA grid assignments; found {len(files)}")
    for i, path in enumerate(files):
        table = pd.read_csv(path, sep="\t")
        table["gene"] = table["gene"].astype(str)
        universe = set(table["gene"])
        target_in_universe = genes & universe
        module_rows = []
        for module, sub in table[table["module_color"] != "grey"].groupby("module_color"):
            members = set(sub["gene"])
            overlap = len(members & genes)
            union = len(members | genes)
            population = len(universe)
            successes = len(target_in_universe)
            p = stats.hypergeom.sf(overlap - 1, population, successes, len(members))
            module_rows.append(
                {
                    "variant": path.stem.replace(".tsv", ""),
                    "module": module,
                    "module_size": len(members),
                    "target_in_universe": successes,
                    "overlap": overlap,
                    "jaccard": overlap / union if union else 0.0,
                    "overlap_p": float(p),
                    "genes": sorted(members),
                }
            )
        if not module_rows:
            raise RuntimeError(f"No non-grey modules in {path}")
        pvals = np.array([row["overlap_p"] for row in module_rows])
        qvals = bh(pvals)
        for row, q in zip(module_rows, qvals):
            row["overlap_q"] = float(q)
        best = max(module_rows, key=lambda x: (x["overlap"], x["jaccard"]))
        score = mean_z_score(expr_fg, best["genes"])
        effect = pooled_iss_effect(score, meta_fg)
        p_perm = blocked_iss_permutation_p(
            score,
            meta_fg,
            n_perm=int(cfg["resampling"]["external_permutations"]),
            seed=int(cfg["resampling"]["seed"]) + 100 + i,
        )
        primary_variant = best["variant"] == "g05000_m30_c25"
        recovery = (
            best["overlap"] >= cfg["go_no_go"]["gate_B"]["flight_blind_overlap_genes_gte"]
            and best["overlap_q"] <= cfg["go_no_go"]["gate_B"]["flight_blind_overlap_bh_q_lte"]
            and best["jaccard"] >= cfg["go_no_go"]["gate_B"]["flight_blind_jaccard_gte"]
        )
        rows.append(
            {
                "variant": best["variant"],
                "primary_variant": primary_variant,
                "best_module": best["module"],
                "module_size": best["module_size"],
                "target_in_universe": best["target_in_universe"],
                "overlap": best["overlap"],
                "jaccard": best["jaccard"],
                "overlap_p": best["overlap_p"],
                "overlap_q": best["overlap_q"],
                "recovery_pass": recovery,
                "projected_iss_effect": effect,
                "projected_permutation_p": p_perm,
                "positive_projected_effect": effect > 0,
            }
        )
        for row in module_rows:
            detail_rows.append({k: v for k, v in row.items() if k != "genes"})

    summary = pd.DataFrame(rows).sort_values("variant")
    detail = pd.DataFrame(detail_rows).sort_values(["variant", "overlap"], ascending=[True, False])
    summary.to_csv(
        outdir / "flight_blind_grid_summary.tsv", sep="\t", index=False
    )
    detail.to_csv(
        outdir / "flight_blind_all_module_overlaps.tsv", sep="\t", index=False
    )
    primary = summary[summary["primary_variant"]]
    if len(primary) != 1:
        raise RuntimeError("Could not uniquely resolve primary flight-blind variant")
    primary = primary.iloc[0]
    recovery_count = int(summary["recovery_pass"].sum())
    positive_count = int(summary["positive_projected_effect"].sum())
    primary_pass = (
        bool(primary["recovery_pass"])
        and primary["projected_iss_effect"] > 0
        and primary["projected_permutation_p"]
        <= cfg["go_no_go"]["gate_B"]["flight_blind_projected_permutation_p_lte"]
    )
    grid_pass = (
        recovery_count >= cfg["go_no_go"]["gate_B"]["grid_recovery_count_gte"]
        and positive_count
        >= cfg["go_no_go"]["gate_B"]["grid_positive_effect_count_gte"]
    )

    existing_gate = json.loads(
        (root / "internal" / "internal_gate_status.json").read_text()
    )
    partial = existing_gate["gate_B_partial"]
    gate_b = {
        "pass": (
            primary_pass
            and grid_pass
            and bool(partial["components"]["gc_subset_exceeds_remainder"])
            and bool(partial["components"]["gc_matched_percentile"])
        ),
        "components": {
            "primary_flight_blind_recovery": primary_pass,
            "grid_robustness": grid_pass,
            "gc_subset_exceeds_remainder": bool(
                partial["components"]["gc_subset_exceeds_remainder"]
            ),
            "gc_matched_percentile": bool(
                partial["components"]["gc_matched_percentile"]
            ),
        },
        "primary_variant": primary.to_dict(),
        "grid_recovery_count": recovery_count,
        "grid_positive_effect_count": positive_count,
        "n_grid_variants": len(summary),
    }
    (outdir / "gate_B_status.json").write_text(
        json.dumps(gate_b, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else x)
        + "\n"
    )
    print(json.dumps(gate_b, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else x))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default="config/grey60_adversarial_reanalysis.yaml"
    )
    parser.add_argument("--outdir", default="")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
