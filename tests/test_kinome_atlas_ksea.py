"""Unit tests for Repair A -- kinome-wide atlas KSEA."""
import numpy as np
import pandas as pd
import pytest

from src.v11.kinome_atlas_ksea import (
    AtlasPSSM,
    ATLAS_POSITIONS,
    ATLAS_PSSM_PATH,
    SPAK_OSR1_KINASES,
    WNK_KINASES,
    _central_residue,
    add_rank,
    build_kinome_net,
    load_atlas_pssms,
    over_representation,
    percentile_by_kinase,
    score_sites,
)
from src.multiomics.regulator_activity import ksea


# #

_RESIDUES = ["R", "D", "E", "A", "S", "T", "G"]
_POS_IX = {p: i for i, p in enumerate(ATLAS_POSITIONS)}
_RES_IX = {r: i for i, r in enumerate(_RESIDUES)}


def _toy_atlas() -> AtlasPSSM:
    """BASO prefers R at -3/-2; ACID prefers D/E at +1/+2/+3. log-weights."""
    K, P, R = 2, len(ATLAS_POSITIONS), len(_RESIDUES)
    logmat = np.zeros((K, P, R), dtype=float)          # neutral baseline log2(1)=0
    # BASO (index 0)
    logmat[0, _POS_IX[-3], _RES_IX["R"]] = 3.0
    logmat[0, _POS_IX[-2], _RES_IX["R"]] = 3.0
    # ACID (index 1)
    logmat[1, _POS_IX[1], _RES_IX["D"]] = 3.0
    logmat[1, _POS_IX[2], _RES_IX["E"]] = 3.0
    logmat[1, _POS_IX[3], _RES_IX["D"]] = 3.0
    return AtlasPSSM(["BASO", "ACID"], list(ATLAS_POSITIONS), _RESIDUES, logmat)


# 13-mers, centre index 6 = S
_BASO_MOTIF = "GGGRRGSGGGGGG"   # idx3=R(-3), idx4=R(-2)
_ACID_MOTIF = "GGGGGGSDEDGGG"   # idx7=D(+1), idx8=E(+2), idx9=D(+3)
_NULL_MOTIF = "GGGGGGSGGGGGG"


def _toy_sites(seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    # 6 BASO substrates, suppressed in flight
    for i in range(6):
        rows.append(("BASOg%d" % i, str(100 + i), _BASO_MOTIF,
                     -2.0 + rng.normal(0, 0.05), 0.001))
    # 6 ACID substrates, up in flight
    for i in range(6):
        rows.append(("ACIDg%d" % i, str(200 + i), _ACID_MOTIF,
                     +2.0 + rng.normal(0, 0.05), 0.001))
    # 48 neutral background, ~null
    for i in range(48):
        rows.append(("BGg%d" % i, str(300 + i), _NULL_MOTIF,
                     rng.normal(0, 0.10), 0.8))
    return pd.DataFrame(
        rows, columns=["gene_symbol", "site_position", "motif",
                       "phospho_effect", "phospho_p_value"]
    )


# #

def test_score_sites_is_consensus_faithful():
    atlas = _toy_atlas()
    scores = score_sites([_BASO_MOTIF, _ACID_MOTIF, _NULL_MOTIF], atlas)
    # rows: motif; cols: [BASO, ACID]
    assert scores[0, 0] > scores[0, 1]          # BASO motif scores higher for BASO
    assert scores[1, 1] > scores[1, 0]          # ACID motif scores higher for ACID
    assert np.isclose(scores[0, 0], 6.0)        # 3 + 3 summed log-weights
    assert np.isclose(scores[2, 0], 0.0)        # null motif: neutral everywhere


def test_central_residue_filters_non_st():
    assert _central_residue(_BASO_MOTIF) == "S"
    assert _central_residue("GGGGGGYGGGGGG") == "Y"   # would be dropped downstream
    assert _central_residue(None) is None


# #

def test_build_kinome_net_assigns_top_percentile_sites():
    atlas = _toy_atlas()
    sites = _toy_sites()
    scores = score_sites(sites["motif"].tolist(), atlas)
    pct = percentile_by_kinase(scores)
    net = build_kinome_net(
        pct, atlas, sites["gene_symbol"].tolist(), sites["site_position"].tolist(),
        threshold=90.0,
    )
    assert set(net.columns) == {"kinase", "substrate_gene", "substrate_site"}
    baso = net[net["kinase"] == "BASO"]["substrate_gene"].tolist()
    assert len(baso) >= 3
    assert all(g.startswith("BASOg") for g in baso)     # only the planted basophilic sites


# #

def test_ksea_recovers_suppressed_kinase_negative_and_ranks_it_first():
    atlas = _toy_atlas()
    sites = _toy_sites()
    scores = score_sites(sites["motif"].tolist(), atlas)
    pct = percentile_by_kinase(scores)
    net = build_kinome_net(
        pct, atlas, sites["gene_symbol"].tolist(), sites["site_position"].tolist(),
        threshold=90.0,
    )
    table = add_rank(ksea(sites[["gene_symbol", "site_position", "phospho_effect"]], net,
                          min_substrates=3))
    row = table.set_index("kinase")
    assert row.loc["BASO", "status"] == "scored"
    assert row.loc["BASO", "z_score"] < 0           # suppressed -> inferred activity down
    assert row.loc["BASO", "p_value"] < 0.05
    assert row.loc["ACID", "z_score"] > 0           # the up kinase separates cleanly
    # z ascending rank: the most-suppressed kinase ranks first
    assert row.loc["BASO", "rank"] == 1.0


# #

def test_over_representation_flags_the_down_kinase():
    atlas = _toy_atlas()
    sites = _toy_sites()
    scores = score_sites(sites["motif"].tolist(), atlas)
    pct = percentile_by_kinase(scores)
    net = build_kinome_net(
        pct, atlas, sites["gene_symbol"].tolist(), sites["site_position"].tolist(),
        threshold=90.0,
    )
    over = over_representation(sites, net, min_substrate_genes=3)
    baso = over[over["kinase"] == "BASO"]
    assert not baso.empty
    assert baso["frac_down"].iloc[0] == 1.0         # every BASO substrate gene is down
    assert baso["p_value"].iloc[0] < 0.05


# #

def test_real_atlas_parses_and_contains_control_kinases():
    pytest.importorskip("openpyxl")
    if not ATLAS_PSSM_PATH.exists():
        pytest.skip("Johnson 2023 atlas workbook not staged")
    atlas = load_atlas_pssms()
    assert len(atlas) == 303
    assert list(atlas.positions) == list(ATLAS_POSITIONS)
    labels = set(atlas.kinases)
    for k in SPAK_OSR1_KINASES + WNK_KINASES:
        assert k in labels, f"{k} missing from atlas"
    # log-weights are finite (no log2(0) leaks)
    assert np.isfinite(atlas.logmat).all()
