from pathlib import Path

import numpy as np
import pandas as pd

from src.networks.tubulointerstitial_state import (
    build_group_effects,
    build_rrrm2_deconvolution_overlay,
    build_state_space_scores,
    build_targeted_gene_panel,
)


def _mechanism_matrix() -> pd.DataFrame:
    return pd.DataFrame({
        "Sample Name": ["r1", "r2", "r3", "r4"],
        "study": ["RRRM-2"] * 4,
        "scenario": ["primary"] * 4,
        "condition": ["GC", "GC", "FLT", "FLT"],
        "Arm": ["ISS-T"] * 4,
        "Age": ["YNG", "OLD", "YNG", "OLD"],
        "ecm_organization": [0.0, 0.0, 2.0, 2.0],
        "integrin_cell_adhesion": [0.0, 0.0, 2.0, 2.0],
        "mmp_adam_proteolysis": [0.0, 0.0, 2.0, 2.0],
        "dct_ncc_wnk_transport": [2.0, 2.0, 0.0, 0.0],
        "tlr4_innate": [0.0, 1.0, 0.0, 1.0],
        "macrophage_inflammation": [0.0, 1.0, 0.0, 1.0],
        "preservation_stress_response": [0.0, 0.0, 0.0, 0.0],
    })


def test_state_space_preserves_native_dct_direction_and_scalar():
    state = build_state_space_scores(_mechanism_matrix())
    flt = state[state["condition"].eq("FLT")]
    gc = state[state["condition"].eq("GC")]
    assert flt["matrix_component"].mean() > gc["matrix_component"].mean()
    assert flt["dct_transport_component"].mean() < gc["dct_transport_component"].mean()
    assert np.allclose(state["matrix_minus_dct"], state["matrix_component"] - state["dct_transport_component"])


def test_rrrm2_deconvolution_join_matches_samples():
    state = build_state_space_scores(_mechanism_matrix())
    cluster = pd.DataFrame({
        "Unnamed: 0": ["r1", "r2", "r3", "r4"],
        "DCT": [0.1, 0.2, 0.3, 0.4],
        "Fibroblast": [0.01, 0.02, 0.03, 0.04],
        "Immune": [0.0, 0.0, 0.1, 0.1],
        "Endothelial": [0.2, 0.2, 0.2, 0.2],
        "PT": [0.5, 0.5, 0.4, 0.4],
        "TAL_LOH": [0.1, 0.1, 0.1, 0.1],
    })
    overlay = build_rrrm2_deconvolution_overlay(state, cluster)
    assert len(overlay) == 4
    assert overlay["prop_DCT"].notna().all()
    assert overlay["prop_Fibroblast"].notna().all()


def test_osd253_scenario_group_effects_split_by_strain():
    rows = []
    for strain in ["C57BL/6J", "C3H/HeJ"]:
        for condition, treatment in [
            ("Space Flight", "White light"),
            ("Ground Control", "Blue light"),
            ("Ground Control Rerun", "White light"),
        ]:
            for rep in range(2):
                rows.append({
                    "Sample Name": f"{strain}_{condition}_{rep}",
                    "study": "OSD-253",
                    "scenario": "all_samples",
                    "condition": condition,
                    "strain": strain,
                    "Factor Value[Spaceflight]": condition,
                    "Factor Value[Treatment]": treatment,
                    "ecm_organization": 2.0 if condition == "Space Flight" else 0.0,
                    "integrin_cell_adhesion": 2.0 if condition == "Space Flight" else 0.0,
                    "mmp_adam_proteolysis": 2.0 if condition == "Space Flight" else 0.0,
                    "dct_ncc_wnk_transport": 0.0 if condition == "Space Flight" else 2.0,
                    "tlr4_innate": 1.0,
                    "macrophage_inflammation": 1.0,
                    "preservation_stress_response": 0.0,
                })
    state = build_state_space_scores(pd.DataFrame(rows))
    effects = build_group_effects(state)
    scenarios = set(effects["scenario"])
    assert {"original_gc_blue_light", "rerun_gc_white_light"} <= scenarios
    strains = set(effects.loc[effects["study"].eq("OSD-253"), "strain"])
    assert {"C57BL/6J", "C3H/HeJ"} <= strains


def test_targeted_gene_panel_handles_missing_genes(tmp_path: Path):
    expr = pd.DataFrame({"s1": [1.0], "s2": [2.0]}, index=["g1"])
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("ensembl_gene_id\tmgi_symbol\ng1\tRbm3\ng2\tMissingInExpr\n")
    meta = pd.DataFrame({
        "Sample Name": ["s1", "s2"],
        "Age": ["YNG", "OLD"],
        "Arm": ["ISS-T", "ISS-T"],
        "EnvGroup": ["GC", "FLT"],
        "condition": ["GC", "FLT"],
        "study": ["RRRM-2", "RRRM-2"],
        "scenario": ["primary", "primary"],
    })
    panel = build_targeted_gene_panel(
        expr,
        str(id_map),
        meta,
        panels={"test": ("Rbm3", "MissingInExpr", "UnmappedGene")},
    )
    assert set(panel["status"]) == {"scored", "not_in_expression", "unmapped"}
    assert panel.loc[panel["status"].eq("scored"), "expression"].notna().all()
