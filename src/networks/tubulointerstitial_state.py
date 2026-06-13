"""Tubulointerstitial state-space analysis utilities."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from src.common import bh_fdr, load_id_map
from src.networks.mechanism_axis import safe_z


MATRIX_COMPONENTS = ("ecm_organization", "integrin_cell_adhesion", "mmp_adam_proteolysis")
DCT_COMPONENT = "dct_ncc_wnk_transport"
IMMUNE_COMPONENTS = ("tlr4_innate", "macrophage_inflammation")
PRESERVATION_COMPONENT = "preservation_stress_response"
STATE_COMPONENTS = (
    "matrix_component",
    "dct_transport_component",
    "immune_context_component",
    "preservation_stress_component",
    "matrix_minus_dct",
)

TARGET_GENE_PANELS: Mapping[str, tuple[str, ...]] = {
    "ECM_adhesion_proteolysis": (
        "Fbn1", "Itga9", "Cdh11", "Adam12", "Mmp14", "Col1a1", "Col1a2",
        "Col3a1", "Col4a1", "Fn1", "Sparc", "Lox", "Postn", "Tnc",
    ),
    "TLR_innate_macrophage": (
        "Tril", "Tlr4", "Ly96", "Cd14", "Myd88", "Ticam1", "Ccl2",
        "Tnf", "Il1b", "Adgre1", "Aif1", "Lyz2", "Csf1r", "Cd68",
    ),
    "DCT_NCC_WNK_transport": (
        "Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Klhl3", "Cul3",
        "Calb1", "Kcnj10", "Kcnj16", "Nedd4l", "Sgk1",
    ),
    "Tubular_metabolic_priority": (
        "Rbm3", "Cyp4a14", "Bhmt", "Per2", "Hmgcr", "Hacl1", "Gpt",
        "Prodh", "Slc22a6", "Ugt1a10", "Plin5", "Lpin1", "Foxo3",
        "Herpud1", "Car4",
    ),
}


@dataclass(frozen=True)
class ContrastSpec:
    """One FLT-minus-control state-space contrast."""

    study: str
    scenario: str
    contrast: str
    arm: str = "NA"
    age: str = "pooled"
    strain: str = "NA"


def context_z(
    df: pd.DataFrame,
    columns: Sequence[str],
    group_cols: Sequence[str] = ("study", "scenario"),
) -> pd.DataFrame:
    """Z-score mechanism columns within each study/scenario context."""
    out = pd.DataFrame(index=df.index)
    grouped = df.groupby(list(group_cols), dropna=False).groups if all(c in df.columns for c in group_cols) else {}
    for col in columns:
        if col not in df.columns:
            out[col] = np.nan
            continue
        z = pd.Series(np.nan, index=df.index, dtype=float)
        if grouped:
            for _, idx in grouped.items():
                idx_list = list(idx)
                z.loc[idx_list] = safe_z(df.loc[idx_list, col].to_numpy(dtype=float))
        else:
            z.loc[:] = safe_z(df[col].to_numpy(dtype=float))
        out[col] = z
    return out


def build_state_space_scores(mechanism_matrix: pd.DataFrame) -> pd.DataFrame:
    """Add tubulointerstitial state-space components to a mechanism matrix."""
    required = set(MATRIX_COMPONENTS + (DCT_COMPONENT,) + IMMUNE_COMPONENTS + (PRESERVATION_COMPONENT,))
    missing = sorted(required - set(mechanism_matrix.columns))
    if missing:
        raise ValueError(f"Mechanism matrix missing required columns: {missing}")

    out = mechanism_matrix.copy()
    z = context_z(out, sorted(required))
    out["matrix_component"] = z.loc[:, MATRIX_COMPONENTS].mean(axis=1, skipna=True)
    out["dct_transport_component"] = z[DCT_COMPONENT]
    out["immune_context_component"] = z.loc[:, IMMUNE_COMPONENTS].mean(axis=1, skipna=True)
    out["preservation_stress_component"] = z[PRESERVATION_COMPONENT]
    out["matrix_minus_dct"] = out["matrix_component"] - out["dct_transport_component"]
    out["state_space_z_scope"] = "within_study_scenario"
    return out


def _welch_p(x: pd.Series, y: pd.Series) -> float:
    x = x.dropna().astype(float)
    y = y.dropna().astype(float)
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    return float(stats.ttest_ind(x, y, equal_var=False).pvalue)


def _effect_rows(
    rows: pd.DataFrame,
    spec: ContrastSpec,
    *,
    flt_mask: pd.Series,
    gc_mask: pd.Series,
    components: Sequence[str] = STATE_COMPONENTS,
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    flt = rows.loc[flt_mask]
    gc = rows.loc[gc_mask]
    for comp in components:
        if comp not in rows.columns:
            continue
        flt_vals = flt[comp].dropna().astype(float)
        gc_vals = gc[comp].dropna().astype(float)
        result.append({
            "study": spec.study,
            "scenario": spec.scenario,
            "contrast": spec.contrast,
            "arm": spec.arm,
            "age": spec.age,
            "strain": spec.strain,
            "component": comp,
            "flt_n": int(len(flt_vals)),
            "gc_n": int(len(gc_vals)),
            "flt_mean": float(flt_vals.mean()) if len(flt_vals) else np.nan,
            "gc_mean": float(gc_vals.mean()) if len(gc_vals) else np.nan,
            "effect_flt_minus_gc": float(flt_vals.mean() - gc_vals.mean()) if len(flt_vals) and len(gc_vals) else np.nan,
            "p_welch": _welch_p(flt_vals, gc_vals),
        })
    return result


def osd253_scenario_mask(rows: pd.DataFrame, scenario: str, strain: str) -> tuple[pd.Series, pd.Series]:
    """Return FLT and control masks for an OSD-253 strain/scenario."""
    strain_mask = rows["strain"].astype(str).eq(strain)
    space = rows["Factor Value[Spaceflight]"].astype(str)
    treatment = rows["Factor Value[Treatment]"].astype(str).str.casefold()
    flt = strain_mask & space.eq("Space Flight") & treatment.eq("white light")
    if scenario == "original_gc_blue_light":
        gc = strain_mask & space.eq("Ground Control") & treatment.eq("blue light")
    elif scenario == "rerun_gc_white_light":
        gc = strain_mask & space.eq("Ground Control Rerun") & treatment.eq("white light")
    else:
        raise ValueError(f"Unsupported OSD-253 scenario: {scenario}")
    return flt, gc


def build_group_effects(state_scores: pd.DataFrame) -> pd.DataFrame:
    """Compute FLT-minus-GC component effects across RRRM-2, OSD-513, and OSD-253."""
    rows: list[dict[str, object]] = []

    rrrm2 = state_scores[state_scores["study"].eq("RRRM-2")].copy()
    if not rrrm2.empty and {"Arm", "Age", "condition"} <= set(rrrm2.columns):
        for arm in ("ISS-T", "LAR"):
            arm_rows = rrrm2[rrrm2["Arm"].astype(str).eq(arm)]
            flt = arm_rows["condition"].astype(str).eq("FLT")
            gc = arm_rows["condition"].astype(str).eq("GC")
            rows.extend(_effect_rows(
                arm_rows,
                ContrastSpec("RRRM-2", "primary", f"{arm}_age_pooled_FLT_minus_GC", arm=arm),
                flt_mask=flt,
                gc_mask=gc,
            ))
            for age in ("YNG", "OLD"):
                age_rows = arm_rows[arm_rows["Age"].astype(str).eq(age)]
                flt = age_rows["condition"].astype(str).eq("FLT")
                gc = age_rows["condition"].astype(str).eq("GC")
                rows.extend(_effect_rows(
                    age_rows,
                    ContrastSpec("RRRM-2", "primary", f"{arm}_{age}_FLT_minus_GC", arm=arm, age=age),
                    flt_mask=flt,
                    gc_mask=gc,
                ))

    osd513 = state_scores[state_scores["study"].eq("OSD-513")].copy()
    if not osd513.empty:
        rows.extend(_effect_rows(
            osd513,
            ContrastSpec("OSD-513", "powered_hgc", "OSD-513_FLT_minus_GC"),
            flt_mask=osd513["condition"].astype(str).eq("FLT"),
            gc_mask=osd513["condition"].astype(str).eq("GC"),
        ))

    osd253 = state_scores[state_scores["study"].eq("OSD-253")].copy()
    for scenario in ("original_gc_blue_light", "rerun_gc_white_light"):
        for strain in ("C57BL/6J", "C3H/HeJ"):
            if osd253.empty:
                continue
            flt, gc = osd253_scenario_mask(osd253, scenario, strain)
            rows.extend(_effect_rows(
                osd253,
                ContrastSpec("OSD-253", scenario, f"OSD-253_{scenario}_{strain}_FLT_minus_GC", strain=strain),
                flt_mask=flt,
                gc_mask=gc,
            ))

    out = pd.DataFrame(rows)
    if not out.empty and out["p_welch"].notna().any():
        out["q_bh_within_component_effects"] = np.nan
        mask = out["p_welch"].notna()
        out.loc[mask, "q_bh_within_component_effects"] = bh_fdr(out.loc[mask, "p_welch"].to_numpy(dtype=float))
    return out


def _sample_col_first(df: pd.DataFrame) -> str:
    if "Sample Name" in df.columns:
        return "Sample Name"
    if "Unnamed: 0" in df.columns:
        return "Unnamed: 0"
    return df.columns[0]


def build_rrrm2_deconvolution_overlay(
    state_scores: pd.DataFrame,
    cluster_props: pd.DataFrame,
    clr_props: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Join RRRM-2 state scores to existing RRRM-2 deconvolution estimates."""
    rrrm2 = state_scores[state_scores["study"].eq("RRRM-2")].copy()
    if rrrm2.empty:
        return pd.DataFrame()

    c = cluster_props.copy()
    sample_col = _sample_col_first(c)
    c = c.rename(columns={sample_col: "Sample Name"})
    value_cols = [col for col in c.columns if col != "Sample Name"]
    c = c.rename(columns={col: f"prop_{col}" for col in value_cols})

    merged = rrrm2.merge(c, on="Sample Name", how="left")
    if clr_props is not None and not clr_props.empty:
        z = clr_props.copy()
        sample_col = _sample_col_first(z)
        z = z.rename(columns={sample_col: "Sample Name"})
        value_cols = [col for col in z.columns if col != "Sample Name"]
        z = z.rename(columns={col: f"clr_{col}" for col in value_cols})
        merged = merged.merge(z, on="Sample Name", how="left")
    return merged


def build_component_correlations(
    overlay: pd.DataFrame,
    *,
    components: Sequence[str] = ("matrix_component", "dct_transport_component", "immune_context_component", "preservation_stress_component", "matrix_minus_dct"),
    cell_columns: Sequence[str] = ("prop_DCT", "prop_Fibroblast", "prop_Immune", "prop_Endothelial", "prop_PT", "prop_TAL_LOH", "clr_DCT", "clr_Fibroblast", "clr_Immune"),
) -> pd.DataFrame:
    """Spearman correlations between RRRM-2 state components and deconvolution estimates."""
    rows: list[dict[str, object]] = []
    if overlay.empty:
        return pd.DataFrame()
    for comp in components:
        if comp not in overlay.columns:
            continue
        for cell in cell_columns:
            if cell not in overlay.columns:
                continue
            sub = overlay[[comp, cell]].dropna()
            if (
                len(sub) < 3
                or sub[comp].astype(float).nunique(dropna=True) < 2
                or sub[cell].astype(float).nunique(dropna=True) < 2
            ):
                rho = p = np.nan
            else:
                rho, p = stats.spearmanr(sub[comp].astype(float), sub[cell].astype(float))
            rows.append({
                "study": "RRRM-2",
                "scope": "all_samples",
                "component": comp,
                "cell_estimate": cell,
                "n": int(len(sub)),
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "p_value": float(p) if np.isfinite(p) else np.nan,
            })
    out = pd.DataFrame(rows)
    if not out.empty and out["p_value"].notna().any():
        mask = out["p_value"].notna()
        out["q_bh"] = np.nan
        out.loc[mask, "q_bh"] = bh_fdr(out.loc[mask, "p_value"].to_numpy(dtype=float))
    return out


def build_targeted_gene_panel(
    expression: pd.DataFrame,
    id_map_path: str,
    sample_metadata: pd.DataFrame,
    panels: Mapping[str, Sequence[str]] = TARGET_GENE_PANELS,
) -> pd.DataFrame:
    """Extract a long targeted residualized-expression table for renal biology genes."""
    expr = expression.copy()
    expr.index = expr.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    id_map = load_id_map(id_map_path)
    symbol_to_ids: dict[str, set[str]] = {}
    id_to_symbol: dict[str, str] = {}
    for row in id_map.itertuples():
        gid = str(row.ensembl_gene_id)
        sym = str(row.mgi_symbol)
        id_to_symbol[gid] = sym
        symbol_to_ids.setdefault(sym.casefold(), set()).add(gid)

    meta = sample_metadata.copy()
    meta = meta[meta["Sample Name"].astype(str).isin(expr.columns)].copy()
    meta_cols = [c for c in ["Sample Name", "Age", "Arm", "EnvGroup", "condition", "study", "scenario"] if c in meta.columns]
    meta = meta[meta_cols].drop_duplicates("Sample Name")

    rows: list[pd.DataFrame] = []
    expression_genes = set(expr.index)
    for panel_class, symbols in panels.items():
        for query in symbols:
            ids = symbol_to_ids.get(str(query).casefold(), set())
            if not ids:
                rows.append(pd.DataFrame([{
                    "panel_class": panel_class,
                    "query_symbol": query,
                    "gene": "",
                    "mgi_symbol": query,
                    "status": "unmapped",
                    "Sample Name": "",
                    "expression": np.nan,
                }]))
                continue
            present = sorted(ids & expression_genes)
            if not present:
                rows.append(pd.DataFrame([{
                    "panel_class": panel_class,
                    "query_symbol": query,
                    "gene": sorted(ids)[0],
                    "mgi_symbol": query,
                    "status": "not_in_expression",
                    "Sample Name": "",
                    "expression": np.nan,
                }]))
                continue
            for gid in present:
                vals = expr.loc[gid, meta["Sample Name"].tolist()].rename("expression").reset_index()
                vals = vals.rename(columns={"index": "Sample Name"})
                vals["panel_class"] = panel_class
                vals["query_symbol"] = query
                vals["gene"] = gid
                vals["mgi_symbol"] = id_to_symbol.get(gid, query)
                vals["status"] = "scored"
                vals = vals.merge(meta, on="Sample Name", how="left")
                rows.append(vals)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True, sort=False)
    leading = ["panel_class", "query_symbol", "gene", "mgi_symbol", "status", "Sample Name"]
    other = [c for c in out.columns if c not in leading]
    return out[leading + other]
