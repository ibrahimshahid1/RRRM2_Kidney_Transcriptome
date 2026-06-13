"""Unit tests for the OSD-462 multi-omics anchor module."""
import numpy as np
import pandas as pd

from src.multiomics.osd462_anchor import (TmtTable, _classify_header,
                                          aligned_pathway_cosine,
                                          assign_match_strata,
                                          collapse_to_gene, compute_flight_effect,
                                          compute_site_flight_effect_lm,
                                          matched_null_test, pathway_effect_vector,
                                          sign_agreement, spearman)


def _synthetic_table(effects, n_per=3):
    """Build a TmtTable where gene g has FL = GC + effects[g] (log2) in both plexes.

    Scaled values are 2**(base + condition shift); BL is set equal to GC.
    """
    plexes = ["Samp1-5", "Samp6-10"]
    conds = ["BL", "FL", "GC"]
    channels = []
    colnames = []
    for plex in plexes:
        for cond in conds:
            for r in range(n_per):
                col = f"{plex}|{cond}|{cond}-{r}"
                channels.append({"column": col, "plex": plex, "condition": cond, "sample": f"{cond}-{r}"})
                colnames.append(col)
    channel_df = pd.DataFrame(channels)

    rows = []
    meta_rows = []
    base = 4.0
    for g, eff in effects.items():
        vals = []
        for ch in channels:
            shift = eff if ch["condition"] == "FL" else 0.0
            vals.append(2.0 ** (base + shift))
        rows.append(vals)
        meta_rows.append({"gene_symbol": g, "n_pep_Samp1-5": 5, "n_pep_Samp6-10": 5})
    scaled = pd.DataFrame(rows, columns=colnames)
    meta = pd.DataFrame(meta_rows)
    return TmtTable(meta=meta, scaled=scaled, channels=channel_df, sheet="synthetic")


def test_classify_header():
    info = _classify_header("Samp1-5~rq_129n_sn scaled", "FL-01")
    assert info == {"plex": "Samp1-5", "value_type": "scaled", "condition": "FL", "sample": "FL-01"}
    assert _classify_header("Gene Symbol", "FL-01") is None
    assert _classify_header("Samp6-10~126_sn_sum", "GC-12")["value_type"] == "sum"


def test_flight_effect_recovers_known_shift():
    tab = _synthetic_table({"Up": 1.0, "Down": -0.5, "Flat": 0.0})
    eff = compute_flight_effect(tab, channel_center=False).set_index("gene_symbol")
    assert np.isclose(eff.loc["Up", "flight_effect"], 1.0, atol=1e-6)
    assert np.isclose(eff.loc["Down", "flight_effect"], -0.5, atol=1e-6)
    assert np.isclose(eff.loc["Flat", "flight_effect"], 0.0, atol=1e-6)
    assert (eff["plex_coverage"] == 2).all()
    assert np.isclose(eff.loc["Up", "n_peptides"], 10.0)


def test_channel_centering_is_invariant_for_balanced_design():
    # When up- and down-movers balance, per-channel medians are equal across
    # conditions so loading normalization must not change individual effects.
    tab = _synthetic_table({"Up": 1.0, "Up2": 1.0, "Down": -1.0, "Down2": -1.0})
    e_raw = compute_flight_effect(tab, channel_center=False).set_index("gene_symbol")
    e_cen = compute_flight_effect(tab, channel_center=True).set_index("gene_symbol")
    for g in ("Up", "Down"):
        assert np.isclose(e_raw.loc[g, "flight_effect"], e_cen.loc[g, "flight_effect"], atol=1e-6)


def test_collapse_to_gene_peptide_weighted():
    tab = _synthetic_table({"A": 1.0})
    eff = compute_flight_effect(tab, channel_center=False)
    # duplicate the gene with a different effect and peptide weight
    extra = eff.iloc[[0]].copy()
    extra["flight_effect"] = -1.0
    extra["n_peptides"] = 30.0
    eff.loc[0, "n_peptides"] = 10.0
    combined = pd.concat([eff, extra], ignore_index=True)
    collapsed = collapse_to_gene(combined).set_index("gene_symbol")
    # weighted mean of (1.0 w=10) and (-1.0 w=30) = -0.5
    assert np.isclose(collapsed.loc["A", "flight_effect"], -0.5, atol=1e-6)
    assert collapsed.loc["A", "n_protein_rows"] == 2


def test_site_lm_effect_and_significance():
    tab = _synthetic_table({"Sig": 1.0, "Null": 0.0}, n_per=3)
    site = compute_site_flight_effect_lm(tab, min_per_condition=3, channel_center=False)
    site = site.set_index("gene_symbol")
    assert np.isclose(site.loc["Sig", "phospho_effect"], 1.0, atol=1e-6)
    assert site.loc["Sig", "phospho_p_value"] < 0.05
    assert np.isclose(site.loc["Null", "phospho_effect"], 0.0, atol=1e-6)


def test_spearman_and_sign_agreement():
    a = np.array([1.0, 2.0, 3.0, 4.0])
    assert np.isclose(spearman(a, a), 1.0)
    assert np.isclose(spearman(a, -a), -1.0)
    assert np.isclose(sign_agreement(np.array([1, -1, 1]), np.array([1, -1, -1])), 2 / 3)


def test_pathway_cosine_direction():
    eff = pd.Series({"g1": 1.0, "g2": 1.0, "g3": -1.0, "g4": -1.0})
    sets = {"up": ["g1", "g2", "x"], "down": ["g3", "g4", "y"]}
    vec, cov = pathway_effect_vector(eff, sets, min_genes=2)
    assert vec["up"] > 0 and vec["down"] < 0
    # identical reference -> cosine 1
    c, shared = aligned_pathway_cosine(vec, vec)
    assert np.isclose(c, 1.0)
    # opposite reference -> cosine -1
    c2, _ = aligned_pathway_cosine(vec, -vec)
    assert np.isclose(c2, -1.0)


def test_matched_null_observed_and_pvalue_range():
    rng = np.random.default_rng(0)
    n = 400
    pool = pd.DataFrame({
        "value": rng.normal(size=n),
        "abundance_log2": rng.normal(size=n),
        "n_peptides": rng.integers(2, 40, size=n),
    })
    pool["match_stratum"] = assign_match_strata(pool).values
    # target = 20 genes with an injected positive shift
    mask = np.zeros(n, dtype=bool)
    mask[:20] = True
    pool.loc[:19, "value"] = pool.loc[:19, "value"] + 5.0
    res = matched_null_test(pool, mask, lambda df: float(df["value"].mean()),
                            pool["match_stratum"], "mean", n_null=2000, rng=rng)
    assert np.isclose(res.observed, pool[mask]["value"].mean())
    assert 0.0 <= res.p_greater <= 1.0
    assert res.p_greater < 0.05  # injected positive shift -> enriched
    assert res.n_target == 20
