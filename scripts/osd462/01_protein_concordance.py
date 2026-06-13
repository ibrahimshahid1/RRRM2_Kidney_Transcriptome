#!/usr/bin/env python3
"""Layer 1 - Protein-level concordance with abundance/peptide-matched null."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _common import (ID_MAP, N_NULL, PROTEOMICS_XLSX, RRRM2_GENE_SCATTER, SEED,
                     anchor_dir, build_symbol_to_ensembl, default_run,
                     load_mechanism_sets, resolve_symbols_to_ensembl, write_manifest)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common import bh_fdr  # noqa: E402
from src.multiomics.osd462_anchor import (assign_match_strata, matched_null_test,  # noqa: E402
                                          sign_agreement, spearman)

# Targeted sets (config name) + whether each is in the BH family or a control.
TARGET_SETS = {
    "dct_ncc_wnk_transport": "family",
    "ecm_organization": "family",
    "tlr4_innate": "family",
    "s1p_s1pr3": "family",
    "tubular_transport_broad": "control",
}


def gene_sets_to_ensembl(names) -> dict[str, list[str]]:
    cfg = load_mechanism_sets()
    bridge = build_symbol_to_ensembl()
    bridge.pop("_n_collisions", None)
    out: dict[str, list[str]] = {}
    for name in names:
        spec = cfg[name]
        genes = spec.get("genes", []) if isinstance(spec, dict) else spec
        out[name] = sorted(set(resolve_symbols_to_ensembl(genes, bridge).values()))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    ap.add_argument("--n-null", type=int, default=N_NULL)
    ap.add_argument("--peptide-filter", type=int, default=2,
                    help="minimum total peptides for inclusion (primary=2)")
    args = ap.parse_args()
    out_dir = anchor_dir(args.run)
    rng = np.random.default_rng(SEED)
    print(f"[layer1] run = {args.run}; n_null = {args.n_null}; peptide_filter >= {args.peptide_filter}")

    master = pd.read_csv(out_dir / "osd462_flight_effects.tsv", sep="\t")
    pool = master.dropna(subset=["protein_flight_effect", "rrrm2_iss_t_rna_effect"]).copy()
    pool = pool[pool["n_peptides"] >= args.peptide_filter].reset_index(drop=True)
    # rebuild strata on this exact pool so matched sampling is internally consistent
    pool["match_stratum"] = assign_match_strata(pool).values
    print(f"[layer1] background pool (protein+ref, >={args.peptide_filter} pep): {len(pool)} genes")

    # genome-wide background context
    gw_spearman_ref = spearman(pool["rrrm2_iss_t_rna_effect"], pool["protein_flight_effect"])
    gw_spearman_osd = spearman(pool["osd462_rna_effect"], pool["protein_flight_effect"]) \
        if "osd462_rna_effect" in pool else np.nan
    gw_pearson = float(np.corrcoef(pool["rrrm2_iss_t_rna_effect"], pool["protein_flight_effect"])[0, 1])
    print(f"[layer1] genome-wide RRRM2-RNA<->protein: Spearman={gw_spearman_ref:.3f} Pearson={gw_pearson:.3f}")

    gene_sets = gene_sets_to_ensembl(TARGET_SETS.keys())
    ens_index = pool.set_index("ENSEMBL")

    rows = []
    for name, kind in TARGET_SETS.items():
        members = [g for g in gene_sets[name] if g in ens_index.index]
        mask = pool["ENSEMBL"].isin(members).to_numpy()
        n_q = int(mask.sum())
        if n_q < 3:
            print(f"[layer1] {name}: only {n_q} quantified members -> skipped")
            continue
        sub = pool[mask]
        mean_eff = float(sub["protein_flight_effect"].mean())
        # predicted direction = sign of RRRM-2 RNA pathway effect
        pred_sign = float(np.sign(sub["rrrm2_iss_t_rna_effect"].mean())) or 1.0
        conc = spearman(sub["rrrm2_iss_t_rna_effect"], sub["protein_flight_effect"])
        sgn = sign_agreement(sub["rrrm2_iss_t_rna_effect"], sub["protein_flight_effect"])

        # within-study (same-cohort) concordance isolates transcript-protein
        # uncoupling from cross-strain divergence
        conc_osd = spearman(sub["osd462_rna_effect"], sub["protein_flight_effect"]) \
            if "osd462_rna_effect" in sub else np.nan

        def f_signed_mean(df, s=pred_sign):
            return s * float(df["protein_flight_effect"].mean())

        def f_conc(df):
            return spearman(df["rrrm2_iss_t_rna_effect"], df["protein_flight_effect"])

        def f_conc_osd(df):
            return spearman(df["osd462_rna_effect"], df["protein_flight_effect"])

        def f_sign(df):
            return sign_agreement(df["rrrm2_iss_t_rna_effect"], df["protein_flight_effect"])

        r_mean = matched_null_test(pool, mask, f_signed_mean, pool["match_stratum"],
                                   "signed_mean_protein_effect", n_null=args.n_null, rng=rng)
        r_conc = matched_null_test(pool, mask, f_conc, pool["match_stratum"],
                                   "spearman_concordance", n_null=args.n_null, rng=rng)
        r_conc_osd = matched_null_test(pool, mask, f_conc_osd, pool["match_stratum"],
                                       "spearman_concordance_within_osd462", n_null=args.n_null, rng=rng)
        r_sign = matched_null_test(pool, mask, f_sign, pool["match_stratum"],
                                   "sign_agreement", n_null=args.n_null, rng=rng)

        rows.append({
            "gene_set": name, "kind": kind, "n_quantified": n_q,
            "predicted_direction": "up" if pred_sign > 0 else "down",
            "mean_protein_effect": mean_eff,
            "signed_mean_effect": pred_sign * mean_eff,
            "signed_mean_null_med": r_mean.null_median,
            "signed_mean_null_lo": r_mean.null_ci_low,
            "signed_mean_null_hi": r_mean.null_ci_high,
            "signed_mean_p": r_mean.p_greater,
            "signed_mean_p_two_sided": r_mean.p_two_sided,
            "concordance": conc,
            "concordance_null_med": r_conc.null_median,
            "concordance_null_lo": r_conc.null_ci_low,
            "concordance_null_hi": r_conc.null_ci_high,
            "concordance_p": r_conc.p_greater,
            "concordance_within_osd462": conc_osd,
            "concordance_within_osd462_null_med": r_conc_osd.null_median,
            "concordance_within_osd462_p": r_conc_osd.p_greater,
            "sign_agreement": sgn,
            "sign_agreement_null_med": r_sign.null_median,
            "sign_agreement_p": r_sign.p_greater,
        })
        print(f"[layer1] {name:24s} n={n_q:2d} meanEff={mean_eff:+.3f} "
              f"(pred {('up' if pred_sign>0 else 'down')}, p={r_mean.p_greater:.4f}) "
              f"conc_xstudy={conc:+.3f}(p={r_conc.p_greater:.4f}) "
              f"conc_within={conc_osd:+.3f}(p={r_conc_osd.p_greater:.4f}) "
              f"signAgree={sgn:.2f}(p={r_sign.p_greater:.4f})")

    res = pd.DataFrame(rows)
    # BH within the targeted family (exclude control) per statistic
    fam = res["kind"] == "family"
    for pcol, qcol in [("signed_mean_p", "signed_mean_q"),
                       ("concordance_p", "concordance_q"),
                       ("sign_agreement_p", "sign_agreement_q")]:
        res[qcol] = np.nan
        if fam.any():
            res.loc[fam, qcol] = bh_fdr(res.loc[fam, pcol].to_numpy())

    res.attrs["genome_wide"] = {}
    out_path = out_dir / "protein_concordance_summary.tsv"
    res.to_csv(out_path, sep="\t", index=False)

    # genome-wide context as a one-row side table
    pd.DataFrame([{
        "n_background": len(pool),
        "gw_spearman_rrrm2_rna_vs_protein": gw_spearman_ref,
        "gw_pearson_rrrm2_rna_vs_protein": gw_pearson,
        "gw_spearman_osd462_rna_vs_protein": gw_spearman_osd,
        "peptide_filter": args.peptide_filter,
    }]).to_csv(out_dir / "protein_concordance_background.tsv", sep="\t", index=False)

    # per-gene target table for the figure
    tgt_all = sorted({g for gs in gene_sets.values() for g in gs})
    keep_cols = ["ENSEMBL", "gene_symbol", "rrrm2_iss_t_rna_effect", "protein_flight_effect",
                 "n_peptides", "abundance_log2"]
    if "osd462_rna_effect" in pool:
        keep_cols.insert(3, "osd462_rna_effect")
    tgt_df = pool[pool["ENSEMBL"].isin(tgt_all)][keep_cols].copy()
    # focused sets win over the broad-transport control for labeling
    setmap: dict[str, str] = {}
    for name in TARGET_SETS:  # iteration order: focused sets before control
        for g in gene_sets[name]:
            setmap.setdefault(g, name)
    tgt_df["gene_set"] = tgt_df["ENSEMBL"].map(setmap)
    tgt_df.to_csv(out_dir / "protein_concordance_target_genes.tsv", sep="\t", index=False)

    write_manifest(
        out_dir,
        analysis="OSD-462 multi-omics anchor - Layer 1 protein concordance",
        inputs={"osd462_flight_effects": out_dir / "osd462_flight_effects.tsv",
                "proteomics_xlsx": PROTEOMICS_XLSX,
                "rrrm2_gene_scatter": RRRM2_GENE_SCATTER, "id_map": ID_MAP},
        outputs={"protein_concordance_summary": out_path,
                 "protein_concordance_background": out_dir / "protein_concordance_background.tsv",
                 "protein_concordance_target_genes": out_dir / "protein_concordance_target_genes.tsv"},
        parameters={"n_null": args.n_null, "seed": SEED, "peptide_filter": args.peptide_filter,
                    "null_model": "abundance(5) x peptide(4) matched random gene sets",
                    "bh_family": "4 targeted sets (control excluded)",
                    "predicted_direction": "sign of RRRM-2 ISS-T RNA pathway effect"},
        name="manifest_layer1.json",
    )
    print(f"[layer1] wrote {out_path}")


if __name__ == "__main__":
    main()
