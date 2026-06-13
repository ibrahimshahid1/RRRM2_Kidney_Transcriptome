#!/usr/bin/env python3
"""Layer 3 - Network-candidate translation check."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _common import (ID_MAP, N_NULL, RRRM2_GENE_PRIORITY, SEED, anchor_dir,
                     build_symbol_to_ensembl, default_run, write_manifest)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.multiomics.osd462_anchor import assign_match_strata, matched_null_test  # noqa: E402

TOP_K = [25, 50, 100]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    ap.add_argument("--n-null", type=int, default=N_NULL)
    args = ap.parse_args()
    out_dir = anchor_dir(args.run)
    rng = np.random.default_rng(SEED)
    print(f"[layer3] run = {args.run}; n_null = {args.n_null}")

    # RRRM-2 network candidates ranked by composite score.
    prio = pd.read_csv(RRRM2_GENE_PRIORITY, sep="\t")
    prio = prio.dropna(subset=["composite_score"]).sort_values("composite_score", ascending=False)
    prio["ENSEMBL"] = prio["gene"].astype(str).str.strip()

    master = pd.read_csv(out_dir / "osd462_flight_effects.tsv", sep="\t")

    # Protein pool
    prot = master.dropna(subset=["protein_flight_effect"]).copy()
    prot["abs_protein_effect"] = prot["protein_flight_effect"].abs()
    prot["match_stratum"] = assign_match_strata(prot).values
    prot_thresh = prot["abs_protein_effect"].quantile(0.90)

    # Phospho gene-level pool
    phos_sites = pd.read_csv(out_dir / "phospho_all_sites.tsv", sep="\t")
    bridge = build_symbol_to_ensembl()
    bridge.pop("_n_collisions", None)
    phos_sites["ENSEMBL"] = phos_sites["gene_symbol"].astype(str).str.lower().map(bridge)
    phos_sites = phos_sites.dropna(subset=["ENSEMBL", "phospho_effect"])
    phos_gene = phos_sites.groupby("ENSEMBL").agg(
        max_abs_phospho=("phospho_effect", lambda s: float(np.nanmax(np.abs(s)))),
        min_p=("phospho_p_value", "min"),
        n_sites=("phospho_effect", "size"),
    ).reset_index()
    # match phospho genes on detectability (number of quantified sites)
    phos_gene["match_stratum"] = assign_match_strata(
        phos_gene.rename(columns={"max_abs_phospho": "abundance_log2", "n_sites": "n_peptides"}),
        abundance_col="abundance_log2", peptide_col="n_peptides",
        n_abundance_bins=1, n_peptide_bins=5).values
    phos_thresh = phos_gene["max_abs_phospho"].quantile(0.90)

    rows = []
    for k in TOP_K:
        cand = set(prio["ENSEMBL"].head(k))

        # protein enrichment: mean |protein effect| of candidates vs matched null
        pmask = prot["ENSEMBL"].isin(cand).to_numpy()
        n_prot = int(pmask.sum())
        if n_prot >= 3:
            r_prot = matched_null_test(
                prot, pmask, lambda df: float(df["abs_protein_effect"].mean()),
                prot["match_stratum"], "mean_abs_protein_effect",
                n_null=args.n_null, rng=rng)
            prot_hits = int((prot[pmask]["abs_protein_effect"] >= prot_thresh).sum())
            prot_obs = r_prot.observed
            prot_p = r_prot.p_greater
            prot_med = r_prot.null_median
        else:
            prot_obs = prot_p = prot_med = np.nan
            prot_hits = 0

        # phospho enrichment: mean max|phospho effect| of candidates vs matched null
        qmask = phos_gene["ENSEMBL"].isin(cand).to_numpy()
        n_phos = int(qmask.sum())
        if n_phos >= 3:
            r_phos = matched_null_test(
                phos_gene, qmask, lambda df: float(df["max_abs_phospho"].mean()),
                phos_gene["match_stratum"], "mean_max_abs_phospho",
                n_null=args.n_null, rng=rng)
            phos_hits = int((phos_gene[qmask]["max_abs_phospho"] >= phos_thresh).sum())
            phos_obs = r_phos.observed
            phos_p = r_phos.p_greater
            phos_med = r_phos.null_median
        else:
            phos_obs = phos_p = phos_med = np.nan
            phos_hits = 0

        rows.append({
            "top_k": k,
            "n_candidates_with_protein": n_prot,
            "protein_mean_abs_effect": prot_obs,
            "protein_null_median": prot_med,
            "protein_enrichment_p": prot_p,
            "protein_top_decile_hits": prot_hits,
            "n_candidates_with_phospho": n_phos,
            "phospho_mean_max_abs_effect": phos_obs,
            "phospho_null_median": phos_med,
            "phospho_enrichment_p": phos_p,
            "phospho_top_decile_hits": phos_hits,
        })
        print(f"[layer3] top{k:3d}: protein n={n_prot} obs={prot_obs:.3f} "
              f"null={prot_med:.3f} p={prot_p:.4f} | "
              f"phospho n={n_phos} obs={phos_obs:.3f} null={phos_med:.3f} p={phos_p:.4f}")

    res = pd.DataFrame(rows)
    res.to_csv(out_dir / "network_translation.tsv", sep="\t", index=False)

    # candidate-level detail table (top-100): which candidates have protein / phospho movement
    top100 = prio.head(100)[["ENSEMBL", "mgi_symbol", "composite_score"]].copy()
    top100 = top100.merge(prot[["ENSEMBL", "protein_flight_effect", "abs_protein_effect"]],
                          on="ENSEMBL", how="left")
    top100 = top100.merge(phos_gene[["ENSEMBL", "max_abs_phospho", "min_p", "n_sites"]],
                          on="ENSEMBL", how="left")
    top100["protein_top_decile"] = top100["abs_protein_effect"] >= prot_thresh
    top100["phospho_top_decile"] = top100["max_abs_phospho"] >= phos_thresh
    top100.to_csv(out_dir / "network_translation_candidates.tsv", sep="\t", index=False)

    verdict = "enriched" if (res["protein_enrichment_p"].min() < 0.05
                             or res["phospho_enrichment_p"].min() < 0.05) else "not_enriched"
    print(f"[layer3] verdict: network candidates are {verdict.upper()} among "
          f"protein/phospho-changing genes")

    write_manifest(
        out_dir,
        analysis="OSD-462 multi-omics anchor - Layer 3 network translation",
        inputs={"gene_axis_priority": RRRM2_GENE_PRIORITY,
                "osd462_flight_effects": out_dir / "osd462_flight_effects.tsv",
                "phospho_all_sites": out_dir / "phospho_all_sites.tsv", "id_map": ID_MAP},
        outputs={"network_translation": out_dir / "network_translation.tsv",
                 "network_translation_candidates": out_dir / "network_translation_candidates.tsv"},
        parameters={"n_null": args.n_null, "seed": SEED, "top_k": TOP_K,
                    "candidate_ranking": "composite_score (LIONESS/node2vec)",
                    "protein_statistic": "mean |protein flight effect|",
                    "phospho_statistic": "mean per-gene max |phospho effect|",
                    "verdict": verdict},
        name="manifest_layer3.json",
    )
    print(f"[layer3] wrote {out_dir / 'network_translation.tsv'}")


if __name__ == "__main__":
    main()
