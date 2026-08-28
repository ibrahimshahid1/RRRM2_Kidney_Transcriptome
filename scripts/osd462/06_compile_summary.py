#!/usr/bin/env python3
"""Layer 6 - Consolidate anchor results into a single summary + run manifest."""
from __future__ import annotations

import argparse
import json

import numpy as np
import pandas as pd

from _common import (ID_MAP, PHOSPHO_XLSX, PROTEOMICS_XLSX, RNA_DE, RNA_VST,
                     RRRM2_GENE_PRIORITY, RRRM2_GENE_SCATTER, anchor_dir,
                     default_run, write_manifest)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    args = ap.parse_args()
    out = anchor_dir(args.run)

    rec = pd.read_csv(out / "osd462_rna_recurrence.tsv", sep="\t").iloc[0]
    summ = pd.read_csv(out / "protein_concordance_summary.tsv", sep="\t").set_index("gene_set")
    bg = pd.read_csv(out / "protein_concordance_background.tsv", sep="\t").iloc[0]
    verdict = pd.read_csv(out / "phospho_axis_verdict.tsv", sep="\t").iloc[0]
    net = pd.read_csv(out / "network_translation.tsv", sep="\t")
    master = pd.read_csv(out / "osd462_flight_effects.tsv", sep="\t")
    phos = pd.read_csv(out / "phospho_axis_summary.tsv", sep="\t")

    ncc = master.loc[master["gene_symbol"] == "Slc12a3"].iloc[0]

    def setrow(name, col):
        return float(summ.loc[name, col]) if name in summ.index else None

    matrix_conc = (setrow("ecm_organization", "signed_mean_p") or 1) < 0.05
    dct_conc = (setrow("dct_ncc_wnk_transport", "signed_mean_p") or 1) < 0.05
    phospho_down = bool(verdict["activity_reduced_supported"])
    net_enriched = bool(net["protein_enrichment_p"].min() < 0.05
                        or net["phospho_enrichment_p"].min() < 0.05)

    # decision-table row (plan section 7)
    if matrix_conc and dct_conc and phospho_down:
        row = "Strongest: transcript+protein+signaling concordant remodeling"
    elif matrix_conc and dct_conc:
        row = "Protein-level concordant remodeling; canonical phospho activity unresolved"
    elif matrix_conc and not dct_conc:
        row = "Matrix protein-anchored; DCT transcript-level"
    else:
        row = ("Protein-abundance concordance NULL for both matrix and DCT; "
               "co-modified phosphosite features provide context, but isolated "
               "canonical NCC/SPAK activity is unresolved")

    summary = {
        "run": args.run,
        "headline": ("RNA matrix-high/DCT-low direction recurs in OSD-462 (gate PASS); "
                     "it does NOT translate to protein-abundance concordance. "
                     "T53/Y65 and S382/S383 co-modified phosphoforms decrease, "
                     "but no isolated canonical NCC/SPAK feature qualifies."),
        "layer4_rna_gate": {
            "point_cosine": float(rec["point_cosine"]),
            "ci": [float(rec["ci_low"]), float(rec["ci_high"])],
            "n_pathways": int(rec["n_pathways"]),
            "loo_cosine_range": [float(rec["loo_cosine_min"]), float(rec["loo_cosine_max"])],
            "recurrence_pass": bool(rec["recurrence_pass"]),
        },
        "layer1_protein": {
            "genome_wide_spearman_rrrm2_rna_vs_protein": float(bg["gw_spearman_rrrm2_rna_vs_protein"]),
            "genome_wide_pearson": float(bg["gw_pearson_rrrm2_rna_vs_protein"]),
            "ecm_signed_mean_effect": setrow("ecm_organization", "signed_mean_effect"),
            "ecm_signed_mean_p_two_sided": setrow("ecm_organization", "signed_mean_p_two_sided"),
            "ecm_concordance": setrow("ecm_organization", "concordance"),
            "dct_signed_mean_effect": setrow("dct_ncc_wnk_transport", "signed_mean_effect"),
            "dct_concordance": setrow("dct_ncc_wnk_transport", "concordance"),
            "any_set_concordant_in_predicted_direction": bool(matrix_conc or dct_conc),
        },
        "layer2_phospho": {
            "checkpoint_pass": bool(verdict["checkpoint_pass"]),
            "canonical_provenance_checkpoint_pass": bool(
                verdict.get("canonical_provenance_checkpoint_pass", False)
            ),
            "ncc_total_protein_effect": float(ncc["protein_flight_effect"]),
            "ncc_rna_effect_osd462": float(ncc["osd462_rna_effect"]),
            "ncc_rna_effect_rrrm2": float(ncc["rrrm2_iss_t_rna_effect"]),
            "ncc_regulatory_mean_phospho": float(verdict["ncc_regulatory_mean_phospho_effect"]),
            "ncc_comodified_canonical_index_mean_phospho": float(
                verdict.get(
                    "ncc_comodified_canonical_index_mean_phospho_effect",
                    np.nan,
                )
            ),
            "spak_osr1_mean_phospho": float(verdict["spak_osr1_mean_phospho_effect"]),
            "chain_mean_phospho": float(verdict["wnk_spak_osr1_ncc_chain_mean"]),
            "n_ncc_regulatory_sites_sig_down": int(verdict["n_ncc_regulatory_sites_sig_down"]),
            "n_spak_osr1_sites_sig_down": int(verdict["n_spak_osr1_sites_sig_down"]),
            "activity_reduced_supported": phospho_down,
            "inference_boundary": (
                "co-modified canonical-indexed context only; no isolated "
                "canonical assay feature"
            ),
        },
        "layer3_network": {
            "protein_enrichment_p_min": float(net["protein_enrichment_p"].min()),
            "phospho_enrichment_p_min": float(net["phospho_enrichment_p"].min()),
            "enriched": net_enriched,
        },
        "hypotheses": {
            "H1_matrix_protein_concordance": "FALSIFIED" if not matrix_conc else "SUPPORTED",
            "H2_dct_protein_concordance": "FALSIFIED" if not dct_conc else "SUPPORTED",
            "H3_phospho_activity": (
                "SUPPORTED" if phospho_down
                else "UNRESOLVED/no-qualified-isolated-feature"
            ),
            "H4_network_translation": "SUPPORTED" if net_enriched else "FALSIFIED",
        },
        "decision_table_row": row,
    }
    (out / "results_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

    write_manifest(
        out,
        analysis="OSD-462 multi-omics anchor - consolidated run manifest",
        inputs={"proteomics_xlsx": PROTEOMICS_XLSX, "phospho_xlsx": PHOSPHO_XLSX,
                "rna_de": RNA_DE, "rna_vst": RNA_VST,
                "rrrm2_gene_scatter": RRRM2_GENE_SCATTER,
                "rrrm2_gene_priority": RRRM2_GENE_PRIORITY, "id_map": ID_MAP},
        outputs={"results_summary": out / "results_summary.json",
                 "osd462_flight_effects": out / "osd462_flight_effects.tsv",
                 "osd462_rna_recurrence": out / "osd462_rna_recurrence.tsv",
                 "protein_concordance_summary": out / "protein_concordance_summary.tsv",
                 "phospho_axis_summary": out / "phospho_axis_summary.tsv",
                 "phospho_axis_verdict": out / "phospho_axis_verdict.tsv",
                 "network_translation": out / "network_translation.tsv",
                 "dashboard_figure": out / "fig_osd462_multiomics_dashboard.pdf",
                 "recurrence_figure": out / "fig_osd462_rna_recurrence.pdf"},
        parameters={"layers": ["00_harmonize", "01_protein_concordance", "02_phospho_axis",
                               "03_network_translation", "04_rna_recurrence", "05_plot_dashboard"],
                    "decision_table_row": row},
        name="manifest.json",
    )
    print(f"\n[summary] wrote {out / 'results_summary.json'} and manifest.json")


if __name__ == "__main__":
    main()
