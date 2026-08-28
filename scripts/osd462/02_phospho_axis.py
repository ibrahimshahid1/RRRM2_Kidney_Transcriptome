#!/usr/bin/env python3
"""Layer 2 - Phosphoprotein activity of the WNK-SPAK/OSR1-NCC axis (conditional)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from _common import (ID_MAP, PHOSPHO_XLSX, SEED, anchor_dir, default_run,
                     write_manifest)

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.common import bh_fdr  # noqa: E402
from src.multiomics.osd462_anchor import (compute_site_flight_effect_lm,  # noqa: E402
                                          parse_tmt_sheet)
from src.multiomics.osd462_stage0 import (  # noqa: E402
    audit_ncc_spak_phosphosites,
)

# DCT / WNK-SPAK/OSR1-NCC chain plus regulators.
AXIS_GENES = ["Slc12a3", "Stk39", "Oxsr1", "Wnk1", "Wnk4", "Klhl3", "Cul3",
              "Nedd4l", "Sgk1", "Calb1"]
# The genes that gate the broad axis-coverage checkpoint (must be quantified).
GATE_GENES = {"Slc12a3", "Stk39", "Oxsr1"}  # NCC + SPAK + OSR1
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    args = ap.parse_args()
    out_dir = anchor_dir(args.run)
    print(f"[layer2] run = {args.run}")

    sheet_cfg = [
        ("siteQuant_360", "gene_symbol", "Site Position", "Protein Id", "single"),
        ("siteQuant_360_compositeSite", "geneSymbol", "sitePosStr", "proteinID", "composite"),
    ]

    all_sites = []
    coverage = []
    for sheet, gcol, scol, idcol, kind in sheet_cfg:
        tab = parse_tmt_sheet(
            PHOSPHO_XLSX, sheet, gene_col=gcol,
            peptide_cols={"Samp1-5": "Samp1-5~num_quant", "Samp6-10": "Samp6-10~num_quant"},
            id_col=idcol, extra_meta_cols=[scol])
        eff = compute_site_flight_effect_lm(tab, min_per_condition=3, channel_center=True)
        eff = eff.rename(columns={scol: "site_position", gcol: "gene_symbol"})
        eff["sheet"] = sheet
        eff["site_kind"] = kind
        all_sites.append(eff)
        for g in AXIS_GENES:
            n = int((eff["gene_symbol"] == g).sum())
            n_quant = int(((eff["gene_symbol"] == g) & np.isfinite(eff["phospho_effect"])).sum())
            coverage.append({"sheet": sheet, "gene": g, "n_sites": n, "n_quantified": n_quant})

    sites = pd.concat(all_sites, ignore_index=True)
    cov = pd.DataFrame(coverage)
    cov_g = cov.groupby("gene", as_index=False)[["n_sites", "n_quantified"]].sum()
    print("[layer2] phosphosite coverage (target chain):")
    print(cov_g.to_string(index=False))

    # Hard checkpoint
    gate_ok = {}
    for g in GATE_GENES:
        gate_ok[g] = bool(cov_g.loc[cov_g["gene"] == g, "n_quantified"].sum() > 0
                          if (cov_g["gene"] == g).any() else False)
    checkpoint_pass = all(gate_ok.values())
    print(
        f"[layer2] BROAD GENE-COVERAGE checkpoint (not canonical provenance) "
        f"{gate_ok} -> {'PASS' if checkpoint_pass else 'FAIL'}"
    )

    # Target-axis site table.  Join residue/phosphoform provenance before any
    # canonical-site claim; a position-only match is not sufficient.
    axis = sites[sites["gene_symbol"].isin(AXIS_GENES) & np.isfinite(sites["phospho_effect"])].copy()
    axis["site_position"] = (
        axis["site_position"].astype(str).str.replace(r"\.0$", "", regex=True)
    )
    provenance = audit_ncc_spak_phosphosites()
    provenance = provenance.rename(
        columns={
            "source_sheet": "sheet",
            "workbook_site_position": "site_position",
        }
    )
    provenance_cols = [
        "gene_symbol",
        "sheet",
        "site_position",
        "residue_resolved_site_labels",
        "sequence",
        "n_hash_modifications_in_sequence",
        "isolated_canonical_assay_feature",
        "residue_indexed_canonical_but_comodified",
        "canonical_claim_class",
    ]
    axis = axis.merge(
        provenance[provenance_cols],
        on=["gene_symbol", "sheet", "site_position"],
        how="left",
        validate="many_to_one",
    )

    # protein-abundance normalization (phospho-occupancy = phospho - protein)
    prot_path = out_dir / "osd462_flight_effects.tsv"
    if prot_path.exists():
        master = pd.read_csv(prot_path, sep="\t")
        prot_map = master.dropna(subset=["protein_flight_effect"]).set_index(
            "gene_symbol")["protein_flight_effect"].to_dict()
        axis["protein_flight_effect"] = axis["gene_symbol"].map(prot_map)
        axis["phospho_occupancy_effect"] = axis["phospho_effect"] - axis["protein_flight_effect"]

    axis = axis.sort_values(["gene_symbol", "site_position"]).reset_index(drop=True)
    axis["phospho_q_value"] = np.nan
    fin = np.isfinite(axis["phospho_p_value"])
    if fin.any():
        axis.loc[fin, "phospho_q_value"] = bh_fdr(axis.loc[fin, "phospho_p_value"].to_numpy())

    keep = ["gene_symbol", "site_position", "site_kind", "sheet", "n_fl", "n_gc",
            "phospho_effect", "phospho_se", "phospho_ci_low", "phospho_ci_high",
            "phospho_p_value", "phospho_q_value", "protein_flight_effect",
            "phospho_occupancy_effect", "residue_resolved_site_labels", "sequence",
            "n_hash_modifications_in_sequence", "isolated_canonical_assay_feature",
            "residue_indexed_canonical_but_comodified", "canonical_claim_class"]
    keep = [c for c in keep if c in axis.columns]
    axis_out = axis[keep]
    axis_out.to_csv(out_dir / "phospho_axis_summary.tsv", sep="\t", index=False)
    sites[np.isfinite(sites["phospho_effect"])][
        ["gene_symbol", "site_position", "site_kind", "phospho_effect", "phospho_se",
         "phospho_p_value", "n_fl", "n_gc"]].to_csv(
        out_dir / "phospho_all_sites.tsv", sep="\t", index=False)
    cov_g.to_csv(out_dir / "phospho_axis_coverage.tsv", sep="\t", index=False)

    # Residue-indexed NCC display.  The strict isolated-canonical set is empty
    # for the current workbook because T53 is carried on a T53/Y65 phosphoform.
    ncc = axis_out[axis_out["gene_symbol"] == "Slc12a3"]
    print("\n[layer2] NCC (Slc12a3) phosphosites:")
    if len(ncc):
        print(ncc[["site_position", "phospho_effect", "phospho_ci_low",
                   "phospho_ci_high", "phospho_p_value"]].to_string(index=False))
    spak = axis_out[axis_out["gene_symbol"].isin(["Stk39", "Oxsr1"])]
    print("\n[layer2] SPAK/OSR1 phosphosites:")
    if len(spak):
        print(spak[["gene_symbol", "site_position", "phospho_effect",
                    "phospho_p_value"]].to_string(index=False))

    # Cross-layer integration verdict
    def axis_mean(genes):
        s = axis_out[axis_out["gene_symbol"].isin(genes)]
        return float(s["phospho_effect"].mean()) if len(s) else np.nan

    ncc_mean = axis_mean(["Slc12a3"])
    spak_osr_mean = axis_mean(["Stk39", "Oxsr1"])
    chain_mean = axis_mean(["Slc12a3", "Stk39", "Oxsr1", "Wnk1", "Wnk4"])
    isolated = axis_out["isolated_canonical_assay_feature"].fillna(False).astype(bool)
    comodified = (
        axis_out["residue_indexed_canonical_but_comodified"]
        .fillna(False)
        .astype(bool)
    )
    ncc_reg = axis_out[(axis_out["gene_symbol"] == "Slc12a3") & isolated]
    ncc_comodified = axis_out[
        (axis_out["gene_symbol"] == "Slc12a3") & comodified
    ]
    ncc_reg_mean = float(ncc_reg["phospho_effect"].mean()) if len(ncc_reg) else np.nan

    # Strict activity support requires isolated canonical NCC and upstream
    # features.  Broad gate-gene and co-modified rows are reported separately.
    sig_down = axis_out[
        (axis_out["gene_symbol"].isin(GATE_GENES))
        & (axis_out["phospho_effect"] < 0)
        & (axis_out["phospho_p_value"] < 0.05)
    ]
    ncc_reg_sig_down = ncc_reg[(ncc_reg["phospho_effect"] < 0)
                               & (ncc_reg["phospho_p_value"] < 0.05)]
    isolated_upstream = axis_out[
        axis_out["gene_symbol"].isin(["Stk39", "Oxsr1"]) & isolated
    ]
    spak_osr_sig_down = isolated_upstream[
        (isolated_upstream["phospho_effect"] < 0)
        & (isolated_upstream["phospho_p_value"] < 0.05)
    ]
    comodified_sig_down = axis_out[
        comodified
        & (axis_out["phospho_effect"] < 0)
        & (axis_out["phospho_p_value"] < 0.05)
    ]
    canonical_checkpoint_pass = bool(len(ncc_reg) and len(isolated_upstream))
    activity_reduced = bool(
        canonical_checkpoint_pass
        and chain_mean < 0
        and np.isfinite(ncc_reg_mean)
        and ncc_reg_mean < 0
        and spak_osr_mean < 0
        and len(ncc_reg_sig_down) > 0
        and len(spak_osr_sig_down) > 0
    )

    verdict = pd.DataFrame([{
        "checkpoint_pass": checkpoint_pass,
        "broad_gene_coverage_checkpoint_pass": checkpoint_pass,
        "canonical_provenance_checkpoint_pass": canonical_checkpoint_pass,
        "ncc_mean_phospho_effect": ncc_mean,
        "ncc_comodified_canonical_index_mean_phospho_effect": (
            float(ncc_comodified["phospho_effect"].mean())
            if len(ncc_comodified) else np.nan
        ),
        "ncc_regulatory_mean_phospho_effect": ncc_reg_mean,
        "spak_osr1_mean_phospho_effect": spak_osr_mean,
        "wnk_spak_osr1_ncc_chain_mean": chain_mean,
        "n_gate_sites_sig_down": int(len(sig_down)),
        "n_isolated_canonical_assay_features": int(isolated.sum()),
        "n_comodified_canonical_index_features": int(comodified.sum()),
        "n_comodified_canonical_index_features_sig_down": int(len(comodified_sig_down)),
        "n_ncc_regulatory_sites_sig_down": int(len(ncc_reg_sig_down)),
        "n_spak_osr1_sites_sig_down": int(len(spak_osr_sig_down)),
        "activity_reduced_supported": activity_reduced,
        "interpretation": (
            "isolated canonical phosphoforms support lower NCC transport activity"
            if activity_reduced else
            "strict canonical activity claim unresolved: no isolated canonical "
            "NCC/SPAK phosphoform qualifies; co-modified site features are context only"
        ),
    }])
    verdict.to_csv(out_dir / "phospho_axis_verdict.tsv", sep="\t", index=False)
    print("\n[layer2] chain mean phospho FL-GC = "
          f"{chain_mean:+.3f}; gate sites sig down = {len(sig_down)}; "
          f"isolated canonical features = {int(isolated.sum())}; "
          f"activity_reduced_supported = {activity_reduced}")

    write_manifest(
        out_dir,
        analysis="OSD-462 multi-omics anchor - Layer 2 phospho axis",
        inputs={"phospho_xlsx": PHOSPHO_XLSX,
                "osd462_flight_effects": prot_path, "id_map": ID_MAP},
        outputs={"phospho_axis_summary": out_dir / "phospho_axis_summary.tsv",
                 "phospho_axis_verdict": out_dir / "phospho_axis_verdict.tsv",
                 "phospho_axis_coverage": out_dir / "phospho_axis_coverage.tsv"},
        parameters={"seed": SEED, "gate_genes": sorted(GATE_GENES),
                    "checkpoint_pass": checkpoint_pass,
                    "site_effect_model": "log2 scaled ~ flight + plex (per site), channel-centered",
                    "occupancy_def": "phospho_effect - protein_flight_effect"},
        name="manifest_layer2.json",
    )
    print(f"[layer2] wrote {out_dir / 'phospho_axis_summary.tsv'}")


if __name__ == "__main__":
    main()
