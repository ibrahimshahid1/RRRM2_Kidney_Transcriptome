"""Cell-type marker-panel decomposition (memo Recommendation 4)."""
from __future__ import annotations

import numpy as np
import pandas as pd

EPS = 1e-12

# Curated mouse kidney marker panels. Sources: Chen et al. 2021 (GSE150338),
KIDNEY_PANELS: dict[str, list[str]] = {
    # --- distal convoluted tubule, split into identity vs transport program ---
    "dct_identity": ["Pvalb", "Trpm6", "Calb1", "Egf", "Tmem52b", "Klhl3"],
    "dct_transport": ["Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Sgk1",
                      "Nedd4l", "Kcnj10", "Kcnj16"],
    # --- other nephron segments (composition references) ---
    "proximal_tubule": ["Slc34a1", "Lrp2", "Slc5a2", "Slc7a13", "Gatm", "Miox"],
    "thick_ascending_limb": ["Umod", "Slc12a1", "Cldn16", "Kcnj1", "Egf"],
    "connecting_tubule_cd": ["Aqp2", "Scnn1a", "Scnn1b", "Scnn1g", "Hsd11b2",
                             "Aqp3", "Fxyd4"],
    "podocyte": ["Nphs1", "Nphs2", "Podxl", "Wt1"],
    # --- interstitial / immune / vascular compartments ---
    "endothelial": ["Pecam1", "Kdr", "Emcn", "Flt1", "Ehd3", "Cdh5"],
    "stromal_fibroblast": ["Pdgfrb", "Pdgfra", "Col1a1", "Col3a1", "Acta2",
                           "Dcn", "Lum"],
    "macrophage_immune": ["Lyz2", "C1qa", "C1qb", "Cd68", "Ptprc", "Adgre1",
                          "Itgam", "Csf1r"],
}

# Compartment scores used by the memo Section 8.2.3 split (re-uses project
# mechanism sets where available; symbols below are the marker proxies here).
COMPARTMENT_PANELS = ("dct_transport", "dct_identity", "stromal_fibroblast",
                      "macrophage_immune", "endothelial")


def panel_flight_effect(gene_stat: pd.DataFrame, panels: dict[str, list[str]],
                        *, gene_col: str = "gene", stat_col: str = "stat",
                        min_genes: int = 2) -> pd.DataFrame:
    """Mean flight stat per panel from a (gene symbol, stat) table.

    Returns one row per panel with the mean stat, n mapped genes, and the genes
    found. Panels with < ``min_genes`` mapped are still reported with NaN mean.
    """
    s = gene_stat.set_index(gene_col)[stat_col]
    s = s[~s.index.duplicated()].astype(float)
    rows = []
    for name, members in panels.items():
        present = [g for g in members if g in s.index]
        vals = s.reindex(present).dropna()
        rows.append({
            "panel": name,
            "n_members": len(members),
            "n_mapped": int(vals.size),
            "mean_stat": float(vals.mean()) if vals.size >= min_genes else np.nan,
            "genes_found": ",".join(present),
        })
    return pd.DataFrame(rows)


def zscore_rows(mat: pd.DataFrame) -> pd.DataFrame:
    mu = mat.mean(axis=1)
    sd = mat.std(axis=1, ddof=1)
    keep = sd > EPS
    return mat.loc[keep].sub(mu[keep], axis=0).div(sd[keep], axis=0)


def per_sample_panel_scores(vst: pd.DataFrame, panels: dict[str, list[str]],
                            sym_to_ens: dict[str, set[str]],
                            *, min_genes: int = 2) -> pd.DataFrame:
    """Per-sample mean-z panel scores from an ENSMUSG-indexed VST matrix.

    Returns panels x samples. ``vst`` rows are ENSMUSG ids; ``sym_to_ens`` maps
    a lowercase symbol to its ENSMUSG id set.
    """
    out = {}
    for name, members in panels.items():
        ids = sorted({e for g in members for e in sym_to_ens.get(g.lower(), set())}
                     & set(vst.index))
        if len(ids) < min_genes:
            continue
        out[name] = zscore_rows(vst.loc[ids]).mean(axis=0)
    return pd.DataFrame(out).T


def decide_scenario(dct_transport_eff: float, dct_identity_eff: float,
                    stromal_eff: float, immune_eff: float,
                    *, drop: float = -0.15, stable: float = 0.15) -> dict:
    """Classify the bulk-RNA ambiguity per memo Section 8.3.

    ``drop``/``stable`` are flight-effect thresholds on the mean-z panel score.
    """
    transport_down = dct_transport_eff <= drop
    identity_down = dct_identity_eff <= drop
    identity_stable = dct_identity_eff >= stable * -1  # i.e. not strongly down
    interstitial_up = (stromal_eff >= stable) or (immune_eff >= stable)

    if transport_down and identity_stable:
        scenario = "transport_program_suppression"
        reading = ("DCT transport program falls while DCT identity markers are "
                   "preserved -> transcriptional suppression of transport more "
                   "likely than loss of DCT cells (strong support for DCT-program "
                   "interpretation).")
    elif transport_down and identity_down:
        scenario = "segment_loss_or_composition_shift"
        reading = ("Both DCT transport and DCT identity fall -> possible segment "
                   "loss / dedifferentiation / composition shift; call it "
                   "DCT/NCC-low, not pure transcriptional suppression.")
    elif interstitial_up and not transport_down:
        scenario = "interstitial_immune_context"
        reading = ("Interstitial/immune markers rise without a clear DCT-transport "
                   "drop -> tissue-context remodeling / possible dilution.")
    else:
        scenario = "indeterminate"
        reading = "No clean scenario; report descriptively."
    return {
        "scenario": scenario, "reading": reading,
        "dct_transport_effect": dct_transport_eff,
        "dct_identity_effect": dct_identity_eff,
        "stromal_effect": stromal_eff, "immune_effect": immune_eff,
    }


__all__ = ["KIDNEY_PANELS", "COMPARTMENT_PANELS", "panel_flight_effect",
           "zscore_rows", "per_sample_panel_scores", "decide_scenario"]
