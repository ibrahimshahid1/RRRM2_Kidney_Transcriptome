"""
WGCNA module-projection external validation.

Projects RRRM-2 WGCNA module gene sets into external OSD cohorts and tests
FLT vs GC module score shifts via permutation. Does NOT rebuild WGCNA in
external cohorts (they're too small).

For each external cohort × module:
  1. Compute module score = mean z-scored expression of module genes present
  2. Welch t-test (FLT vs GC) + label-permutation p-value
  3. BH correction across modules within each study
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT, bh_fdr

# Reuse the data-loading helpers from the LIONESS external validation
from src.validation.osd_external_validation import (
    STUDY_SPECS,
    load_vst_matrix,
    resolve_study_file,
    infer_flt_gc_samples,
)


def load_module_gene_lists(wgcna_dir: Path) -> dict[str, list[str]]:
    """Load WGCNA module gene lists from module_gene_lists/ directory."""
    list_dir = wgcna_dir / "module_gene_lists"
    modules = {}
    for f in sorted(list_dir.glob("*.txt")):
        color = f.stem
        if color == "grey":
            continue  # skip unassigned
        genes = [g.strip() for g in f.read_text().splitlines() if g.strip()]
        if len(genes) >= 5:
            modules[color] = genes
    return modules


def compute_module_scores(
    vst: pd.DataFrame,
    samples: list[str],
    module_genes: dict[str, list[str]],
) -> pd.DataFrame:
    """Compute module scores (mean z-scored expression) for each sample.

    Returns DataFrame: modules × samples.
    """
    available_genes = set(vst.index.astype(str))
    scores = {}
    coverage = {}
    for mod, genes in module_genes.items():
        present = [g for g in genes if g in available_genes]
        if len(present) < 3:
            continue
        coverage[mod] = len(present) / len(genes)
        expr = vst.loc[present, samples].values.astype(float)
        # Z-score each gene across samples, then mean across genes
        mu = expr.mean(axis=1, keepdims=True)
        sd = expr.std(axis=1, keepdims=True, ddof=1)
        sd[sd < 1e-12] = 1.0
        z = (expr - mu) / sd
        scores[mod] = z.mean(axis=0)
    return pd.DataFrame(scores, index=samples), coverage


def test_module_shift(
    scores: pd.DataFrame,
    flt_samples: list[str],
    gc_samples: list[str],
    k_perm: int = 5000,
    seed: int = 42,
) -> pd.DataFrame:
    """Permutation-calibrated Welch t-test for each module's FLT vs GC shift."""
    rng = np.random.default_rng(seed)
    all_samples = flt_samples + gc_samples
    labels = np.array([1] * len(flt_samples) + [0] * len(gc_samples))
    scores_aligned = scores.loc[all_samples]

    results = []
    for mod in scores_aligned.columns:
        vals = scores_aligned[mod].values
        flt_vals = vals[labels == 1]
        gc_vals = vals[labels == 0]

        mean_diff = flt_vals.mean() - gc_vals.mean()
        pooled_se = np.sqrt(
            flt_vals.var(ddof=1) / len(flt_vals) +
            gc_vals.var(ddof=1) / len(gc_vals)
        )
        t_obs = mean_diff / pooled_se if pooled_se > 1e-12 else 0.0

        # Permutation null
        exceed = 0
        for _ in range(k_perm):
            perm = labels.copy()
            rng.shuffle(perm)
            diff_p = vals[perm == 1].mean() - vals[perm == 0].mean()
            se_p = np.sqrt(
                vals[perm == 1].var(ddof=1) / (perm == 1).sum() +
                vals[perm == 0].var(ddof=1) / (perm == 0).sum()
            )
            t_p = diff_p / se_p if se_p > 1e-12 else 0.0
            if abs(t_p) >= abs(t_obs):
                exceed += 1

        p_perm = (1 + exceed) / (1 + k_perm)

        results.append({
            "module": mod,
            "mean_flt": float(flt_vals.mean()),
            "mean_gc": float(gc_vals.mean()),
            "mean_diff": float(mean_diff),
            "t_stat": float(t_obs),
            "perm_p": float(p_perm),
            "n_flt": len(flt_samples),
            "n_gc": len(gc_samples),
        })

    df = pd.DataFrame(results)
    df["perm_q"] = bh_fdr(df["perm_p"].values)
    return df.sort_values("perm_p")


def run_study(
    study: str,
    external_root: Path,
    wgcna_dir: Path,
    outdir: Path,
    k_perm: int = 5000,
    seed: int = 42,
) -> pd.DataFrame:
    """Run module-projection validation for one external cohort."""
    spec = STUDY_SPECS[study]
    study_dir = external_root / study
    vst_path = resolve_study_file(study_dir, spec["vst"], spec.get("vst_glob", "*VST_Counts*.csv"))
    vst = load_vst_matrix(vst_path)

    flt_samples, gc_samples = infer_flt_gc_samples(
        list(vst.columns),
        spec["flt_pattern"],
        spec["gc_pattern"],
        include_pattern=spec.get("include_pattern"),
        exclude_pattern=spec.get("exclude_pattern"),
    )

    module_genes = load_module_gene_lists(wgcna_dir)
    print(f"  {study}: {len(flt_samples)} FLT, {len(gc_samples)} GC, "
          f"{len(module_genes)} modules")

    scores, coverage = compute_module_scores(vst, flt_samples + gc_samples, module_genes)
    print(f"  Gene coverage: {', '.join(f'{m}={c:.0%}' for m, c in sorted(coverage.items()) if m in ['grey60', 'red', 'midnightblue', 'green', 'salmon', 'pink'])}")

    results = test_module_shift(scores, flt_samples, gc_samples, k_perm, seed)
    results["study"] = study
    results["analysis"] = spec["analysis"]

    study_out = outdir / study
    study_out.mkdir(parents=True, exist_ok=True)
    results.to_csv(study_out / "wgcna_module_projection.tsv", sep="\t", index=False)

    # Coverage table
    cov_df = pd.DataFrame([{"module": m, "coverage": c, "n_genes_total": len(module_genes[m]),
                            "n_genes_matched": int(c * len(module_genes[m]))}
                           for m, c in coverage.items()])
    cov_df.to_csv(study_out / "wgcna_gene_coverage.tsv", sep="\t", index=False)

    manifest = {
        "study": study,
        "analysis": spec["analysis"],
        "method": "WGCNA module projection (mean z-scored expression)",
        "vst": str(vst_path),
        "n_flt": len(flt_samples),
        "n_gc": len(gc_samples),
        "k_permutations": k_perm,
        "n_modules_tested": len(scores.columns),
        "claim_boundary": "Module activity comparison; modules defined in RRRM-2 FLT+GC only",
    }
    (study_out / "wgcna_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return results


def main():
    ap = argparse.ArgumentParser(
        description="WGCNA module-projection external validation")
    ap.add_argument("--external_root", default=str(REPO_ROOT / "data/external/osdr"))
    ap.add_argument("--wgcna_dir", required=True,
                    help="Path to WGCNA output directory (contains module_gene_lists/)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--studies", default="auto",
                    help="Comma-separated studies, or auto for all available")
    ap.add_argument("--k_perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    external_root = Path(args.external_root)
    wgcna_dir = Path(args.wgcna_dir)

    if args.studies.strip().lower() == "auto":
        studies = [s for s in STUDY_SPECS if (external_root / s).exists()]
    else:
        studies = [s.strip() for s in args.studies.split(",") if s.strip()]

    if not studies:
        raise ValueError(f"No supported external studies found under {external_root}")

    print("=" * 60)
    print("WGCNA Module-Projection External Validation")
    print("=" * 60)

    outdir = Path(args.outdir)
    all_results = []
    for i, study in enumerate(studies):
        print(f"\n[{study}]")
        results = run_study(
            study, external_root, wgcna_dir, outdir,
            k_perm=args.k_perm, seed=args.seed + i,
        )
        all_results.append(results)

        # Print top results
        sig = results[results["perm_q"] < 0.10]
        if len(sig):
            print(f"  Significant modules (q < 0.10):")
            for _, r in sig.iterrows():
                print(f"    {r['module']:<15} diff={r['mean_diff']:+.3f} "
                      f"t={r['t_stat']:.2f} p={r['perm_p']:.4f} q={r['perm_q']:.4f}")
        else:
            print(f"  No modules at q < 0.10")
            # Show top 3
            for _, r in results.head(3).iterrows():
                print(f"    {r['module']:<15} diff={r['mean_diff']:+.3f} "
                      f"t={r['t_stat']:.2f} p={r['perm_p']:.4f}")

    # Combined summary
    if all_results:
        summary = pd.concat(all_results, ignore_index=True)
        summary.to_csv(outdir / "wgcna_external_summary.tsv", sep="\t", index=False)

    print(f"\nOutputs: {outdir}")


if __name__ == "__main__":
    main()
