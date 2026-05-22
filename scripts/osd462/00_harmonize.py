#!/usr/bin/env python3
"""Layer 0 - Harmonization.

Builds the master ``osd462_flight_effects.tsv`` joining, per gene:
  * RRRM-2 ISS-T RNA flight effect (reference; from lar_reversal_gene_scatter)
  * OSD-462 RNA flight effect (Space Flight - Ground Control; OSDR DE table)
  * OSD-462 protein flight effect (FL - GC, plex-corrected; TMT workbook)
  * peptide count, plex coverage, protein abundance, matched-null strata bins

Usage::

    python scripts/osd462/00_harmonize.py [--run RUN_NAME]
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from _common import (ID_MAP, PROTEOMICS_XLSX, RNA_DE, RNA_FLIGHT_ADJP,
                     RNA_FLIGHT_COL, RRRM2_GENE_SCATTER, anchor_dir,
                     build_symbol_to_ensembl, default_run, write_manifest)

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.multiomics.osd462_anchor import (assign_match_strata, collapse_to_gene,  # noqa: E402
                                          compute_flight_effect, parse_tmt_sheet)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    args = ap.parse_args()
    out_dir = anchor_dir(args.run)
    print(f"[layer0] run = {args.run}")
    print(f"[layer0] output dir = {out_dir}")

    # 1. OSD-462 proteomics flight effect -> gene level ----------------------
    print("[layer0] parsing proteomics workbook ...")
    tab = parse_tmt_sheet(
        PROTEOMICS_XLSX, "protein_quant_2721", gene_col="Gene Symbol",
        peptide_cols={"Samp1-5": "Samp1-5 Peptides", "Samp6-10": "Samp6-10 Peptides"},
        id_col="Protein Id",
    )
    prot_eff = compute_flight_effect(tab)
    prot_eff.to_csv(out_dir / "protein_effects_by_row.tsv", sep="\t", index=False)
    gene_prot = collapse_to_gene(prot_eff, require_both_plex=True)
    gene_prot_1plex = collapse_to_gene(prot_eff, require_both_plex=False)
    print(f"[layer0] proteins quantified: {tab.n_rows}; "
          f"gene-level both-plex: {len(gene_prot)}; any-plex: {len(gene_prot_1plex)}")

    # 2. OSD-462 RNA flight effect (SF - GC) ---------------------------------
    print("[layer0] loading OSD-462 RNA differential-expression table ...")
    rna = pd.read_csv(RNA_DE, usecols=["ENSEMBL", "SYMBOL", RNA_FLIGHT_COL, RNA_FLIGHT_ADJP],
                      dtype={"ENSEMBL": str, "SYMBOL": str})
    rna = rna.rename(columns={RNA_FLIGHT_COL: "osd462_rna_effect",
                              RNA_FLIGHT_ADJP: "osd462_rna_adj_p"})
    rna["ENSEMBL"] = rna["ENSEMBL"].str.strip()
    rna = rna.dropna(subset=["ENSEMBL"]).drop_duplicates("ENSEMBL")

    # 3. RRRM-2 ISS-T reference RNA effect -----------------------------------
    print("[layer0] loading RRRM-2 ISS-T reference gene effect ...")
    ref = pd.read_csv(RRRM2_GENE_SCATTER, sep="\t",
                      usecols=["gene", "mgi_symbol", "iss_t_effect"])
    ref = ref.rename(columns={"gene": "ENSEMBL", "iss_t_effect": "rrrm2_iss_t_rna_effect"})
    ref["ENSEMBL"] = ref["ENSEMBL"].astype(str).str.strip()
    ref = ref.drop_duplicates("ENSEMBL")

    # 4. ID bridge: protein symbol -> ENSMUSG --------------------------------
    bridge = build_symbol_to_ensembl()
    n_coll = bridge.pop("_n_collisions", 0)
    sym = gene_prot["gene_symbol"].astype(str)
    gene_prot = gene_prot.assign(
        ENSEMBL=sym.str.lower().map(lambda k: bridge.get(k)))
    unresolved = gene_prot[gene_prot["ENSEMBL"].isna()]
    print(f"[layer0] symbol->ENSMUSG: {gene_prot['ENSEMBL'].notna().sum()} resolved, "
          f"{len(unresolved)} unresolved; DE-table symbol collisions skipped: {n_coll}")
    gene_prot = gene_prot.dropna(subset=["ENSEMBL"]).copy()

    # collapse protein rows that resolve to the same ENSEMBL (peptide-weighted)
    if gene_prot["ENSEMBL"].duplicated().any():
        def _agg(g):
            w = g["n_peptides"].to_numpy(float)
            w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
            if w.sum() <= 0:
                w = np.ones(len(g))
            return pd.Series({
                "gene_symbol": g.sort_values("n_peptides", ascending=False)["gene_symbol"].iloc[0],
                "protein_flight_effect": float(np.average(g["flight_effect"], weights=w)),
                "n_peptides": float(g["n_peptides"].sum()),
                "abundance_log2": float(g["abundance_log2"].mean()),
                "plex_coverage": int(g["plex_coverage"].max()),
                "n_protein_rows": int(g["n_protein_rows"].sum()),
            })
        gene_prot = gene_prot.groupby("ENSEMBL", sort=False).apply(_agg).reset_index()
    else:
        gene_prot = gene_prot.rename(columns={"flight_effect": "protein_flight_effect"})

    # 5. Master join ----------------------------------------------------------
    master = (rna.merge(ref[["ENSEMBL", "rrrm2_iss_t_rna_effect"]], on="ENSEMBL", how="outer")
                 .merge(gene_prot, on="ENSEMBL", how="outer"))
    # prefer the proteomics symbol, else RNA SYMBOL
    master["gene_symbol"] = master["gene_symbol"].fillna(master["SYMBOL"])
    master = master.drop(columns=["SYMBOL"])

    # matched-null strata over the protein-quantified background
    prot_pool = master[master["protein_flight_effect"].notna()].copy()
    strata = assign_match_strata(prot_pool)
    master["match_stratum"] = pd.Series(pd.NA, index=master.index, dtype="object")
    master.loc[prot_pool.index, "match_stratum"] = strata.values

    cols = ["ENSEMBL", "gene_symbol", "rrrm2_iss_t_rna_effect", "osd462_rna_effect",
            "osd462_rna_adj_p", "protein_flight_effect", "n_peptides", "plex_coverage",
            "n_protein_rows", "abundance_log2", "match_stratum"]
    master = master[[c for c in cols if c in master.columns]]
    master = master.sort_values("ENSEMBL").reset_index(drop=True)

    out_path = out_dir / "osd462_flight_effects.tsv"
    master.to_csv(out_path, sep="\t", index=False)

    # coverage summary
    has_rna = master["osd462_rna_effect"].notna()
    has_ref = master["rrrm2_iss_t_rna_effect"].notna()
    has_prot = master["protein_flight_effect"].notna()
    print(f"[layer0] master rows: {len(master)}")
    print(f"[layer0]   OSD-462 RNA effect: {has_rna.sum()}")
    print(f"[layer0]   RRRM-2 ref effect : {has_ref.sum()}")
    print(f"[layer0]   protein effect    : {has_prot.sum()}")
    print(f"[layer0]   RNA+ref+protein   : {(has_rna & has_ref & has_prot).sum()}")
    print(f"[layer0]   ref+protein       : {(has_ref & has_prot).sum()}")

    # also drop the 1-plex sensitivity table (sensitivity per plan §5/§6)
    gene_prot_1plex.to_csv(out_dir / "protein_effects_gene_anyplex.tsv", sep="\t", index=False)

    write_manifest(
        out_dir,
        analysis="OSD-462 multi-omics anchor - Layer 0 harmonization",
        inputs={"proteomics_xlsx": PROTEOMICS_XLSX, "rna_de": RNA_DE,
                "rrrm2_gene_scatter": RRRM2_GENE_SCATTER, "id_map": ID_MAP},
        outputs={"osd462_flight_effects": out_path,
                 "protein_effects_by_row": out_dir / "protein_effects_by_row.tsv",
                 "protein_effects_gene_anyplex": out_dir / "protein_effects_gene_anyplex.tsv"},
        parameters={"flight_effect_def": "FL - GC within plex, plex-averaged, "
                                          "log2 scaled S/N, per-channel median centered",
                    "rna_flight_def": RNA_FLIGHT_COL,
                    "primary_inclusion": "quantified in both plexes, >=2 peptides total",
                    "de_symbol_collisions_skipped": int(n_coll),
                    "n_unresolved_symbols": int(len(unresolved))},
        name="manifest_layer0.json",
    )
    print(f"[layer0] wrote {out_path}")


if __name__ == "__main__":
    main()
