#!/usr/bin/env python3
"""Layer 4 - Same-modality RNA recurrence gate."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _common import (ID_MAP, N_BOOTSTRAP, RNA_SAMPLES, RNA_VST, RRRM2_GENE_SCATTER,
                     SEED, anchor_dir, build_symbol_to_ensembl, default_run,
                     load_mechanism_sets, resolve_symbols_to_ensembl, write_manifest)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.multiomics.osd462_anchor import (aligned_pathway_cosine,  # noqa: E402
                                          pathway_effect_vector)
from src.networks.contrast_vectors import cosine, percentile_ci  # noqa: E402

RECUR_POINT_MIN = 0.20   # same thresholds as cross-OSDR framework
RECUR_CI_MIN = 0.0


def gene_sets_to_ensembl() -> dict[str, list[str]]:
    cfg = load_mechanism_sets()
    bridge = build_symbol_to_ensembl()
    bridge.pop("_n_collisions", None)
    out: dict[str, list[str]] = {}
    for name, spec in cfg.items():
        genes = spec.get("genes", []) if isinstance(spec, dict) else spec
        mapping = resolve_symbols_to_ensembl(genes, bridge)
        out[name] = sorted(set(mapping.values()))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    ap.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    args = ap.parse_args()
    out_dir = anchor_dir(args.run)
    rng = np.random.default_rng(SEED)
    print(f"[layer4] run = {args.run}")

    gene_sets = gene_sets_to_ensembl()

    # RRRM-2 reference pathway vector (fixed input).
    ref = pd.read_csv(RRRM2_GENE_SCATTER, sep="\t", usecols=["gene", "iss_t_effect"])
    ref_effect = ref.set_index("gene")["iss_t_effect"].astype(float)
    ref_vec, ref_cov = pathway_effect_vector(ref_effect, gene_sets, min_genes=3)

    # OSD-462 RNA: VST matrix + sample labels, restricted to pathway genes.
    print("[layer4] loading OSD-462 mRNA VST matrix ...")
    samples = pd.read_csv(RNA_SAMPLES)
    samples.columns = ["sample", "condition"]
    sf = samples.loc[samples["condition"] == "Space.Flight", "sample"].tolist()
    gc = samples.loc[samples["condition"] == "Ground.Control", "sample"].tolist()

    pathway_genes = sorted({g for gs in gene_sets.values() for g in gs})
    vst = pd.read_csv(RNA_VST, index_col=0)
    vst.index = vst.index.astype(str).str.strip().str.strip('"')
    keep = [g for g in pathway_genes if g in vst.index]
    vst = vst.loc[keep]
    sf_mat = vst[sf].to_numpy(dtype=float)
    gc_mat = vst[gc].to_numpy(dtype=float)
    print(f"[layer4] pathway genes in VST: {len(keep)}; SF={len(sf)} GC={len(gc)}")

    def osd462_vector_from_indices(sf_idx, gc_idx) -> pd.Series:
        eff = sf_mat[:, sf_idx].mean(axis=1) - gc_mat[:, gc_idx].mean(axis=1)
        eff_series = pd.Series(eff, index=keep)
        vec, _ = pathway_effect_vector(eff_series, gene_sets, min_genes=3)
        return vec

    # Point estimate.
    point_vec = osd462_vector_from_indices(np.arange(len(sf)), np.arange(len(gc)))
    point_cos, shared = aligned_pathway_cosine(point_vec, ref_vec)
    print(f"[layer4] point cosine = {point_cos:.4f} over {len(shared)} pathways")

    # Bootstrap (resample SF and GC samples with replacement).
    boot = np.full(args.n_bootstrap, np.nan)
    for b in range(args.n_bootstrap):
        sf_idx = rng.integers(0, len(sf), len(sf))
        gc_idx = rng.integers(0, len(gc), len(gc))
        vec = osd462_vector_from_indices(sf_idx, gc_idx)
        c, _ = aligned_pathway_cosine(vec, ref_vec)
        boot[b] = c
    med, lo, hi = percentile_ci(boot, alpha=0.05)
    recurrence_pass = bool(np.isfinite(point_cos) and point_cos >= RECUR_POINT_MIN and lo > RECUR_CI_MIN)
    print(f"[layer4] bootstrap CI = [{lo:.4f}, {hi:.4f}]  recurrence_pass = {recurrence_pass}")

    # Leave-one-pathway-out.
    loo_rows = []
    for drop in shared:
        a = point_vec.drop(index=drop, errors="ignore")
        b_ = ref_vec.drop(index=drop, errors="ignore")
        c, sh = aligned_pathway_cosine(a, b_)
        loo_rows.append({"dropped_pathway": drop, "cosine": c, "n_pathways": len(sh)})
    loo = pd.DataFrame(loo_rows)

    # Pathway-level effects table (matrix-high / DCT-low readout).
    pv = pd.DataFrame({
        "pathway": shared,
        "osd462_rna_pathway_effect": [point_vec.get(p, np.nan) for p in shared],
        "rrrm2_iss_t_pathway_effect": [ref_vec.get(p, np.nan) for p in shared],
    })
    pv["sign_agree"] = np.sign(pv["osd462_rna_pathway_effect"]) == np.sign(pv["rrrm2_iss_t_pathway_effect"])

    # Outputs.
    summary = pd.DataFrame([{
        "comparison": "OSD-462_RNA_SF_minus_GC_vs_RRRM2_ISS-T",
        "resolution": "pathway",
        "n_pathways": len(shared),
        "point_cosine": point_cos,
        "boot_median": med, "ci_low": lo, "ci_high": hi,
        "n_bootstrap": int(np.isfinite(boot).sum()),
        "recurrence_point_min": RECUR_POINT_MIN,
        "recurrence_ci_min": RECUR_CI_MIN,
        "recurrence_pass": recurrence_pass,
        "loo_cosine_min": float(loo["cosine"].min()) if len(loo) else np.nan,
        "loo_cosine_max": float(loo["cosine"].max()) if len(loo) else np.nan,
    }])
    summary.to_csv(out_dir / "osd462_rna_recurrence.tsv", sep="\t", index=False)
    pv.to_csv(out_dir / "osd462_rna_pathway_effects.tsv", sep="\t", index=False)
    loo.to_csv(out_dir / "osd462_rna_recurrence_loo.tsv", sep="\t", index=False)
    ref_cov.to_csv(out_dir / "osd462_pathway_coverage.tsv", sep="\t", index=False)
    pd.DataFrame({"cosine_bootstrap": boot}).to_csv(
        out_dir / "osd462_rna_recurrence_bootstrap.tsv", sep="\t", index=False)

    print("[layer4] pathway effects (matrix-high / DCT-low readout):")
    print(pv.to_string(index=False))

    write_manifest(
        out_dir,
        analysis="OSD-462 multi-omics anchor - Layer 4 RNA recurrence gate",
        inputs={"rna_vst": RNA_VST, "rna_samples": RNA_SAMPLES,
                "rrrm2_gene_scatter": RRRM2_GENE_SCATTER, "id_map": ID_MAP},
        outputs={"osd462_rna_recurrence": out_dir / "osd462_rna_recurrence.tsv",
                 "osd462_rna_pathway_effects": out_dir / "osd462_rna_pathway_effects.tsv",
                 "osd462_rna_recurrence_loo": out_dir / "osd462_rna_recurrence_loo.tsv"},
        parameters={"n_bootstrap": args.n_bootstrap, "seed": SEED,
                    "recurrence_point_min": RECUR_POINT_MIN,
                    "recurrence_ci_min": RECUR_CI_MIN,
                    "gate": "Layer 4 must pass for cross-modality claims to be interpretable"},
        name="manifest_layer4.json",
    )
    gate = "PASS" if recurrence_pass else "FAIL"
    print(f"[layer4] GATE: {gate}")


if __name__ == "__main__":
    main()
