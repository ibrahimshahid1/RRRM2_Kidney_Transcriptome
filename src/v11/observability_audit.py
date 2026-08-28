"""Module 3 - proteome observability-bias audit."""
from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy as np
import pandas as pd

from src.multiomics.phenotype_anchor import (
    NCC_NONREGULATORY_SITES,
    OSD462_COMODIFIED_CANONICAL_INDEX_FEATURES,
)
from src.v11.matched_null import (
    assign_match_strata,
    ols_slope,
    run_matched_null,
    sign_agreement_rate,
)



def collapse_observability_to_gene(
    by_row: pd.DataFrame,
    gene_col: str = "gene_symbol",
    n_channels_col: str = "n_channels_used",
    peptide_col: str = "n_peptides",
    abundance_col: str = "abundance_log2",
) -> pd.DataFrame:
    """Peptide-weighted per-gene collapse of per-protein observability.

    Inputs come from ``osd462_anchor/protein_effects_by_row.tsv`` which is
    one row per protein with ``n_channels_used`` (count of finite scaled
    S/N channels) and ``n_peptides``.  Output is one row per gene_symbol.
    """
    df = by_row.copy()
    df[gene_col] = df[gene_col].astype(str).str.strip()
    df = df[df[gene_col].ne("")]
    for col in (n_channels_col, peptide_col, abundance_col):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    n_channels_total = float(df[n_channels_col].max(skipna=True) or 1.0)

    def _agg(g: pd.DataFrame) -> pd.Series:
        w = g[peptide_col].fillna(0).to_numpy(dtype=float)
        w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
        if w.sum() <= 0:
            w = np.ones(len(g))
        nch = g[n_channels_col].to_numpy(dtype=float)
        return pd.Series({
            "n_channels_used_mean": float(np.average(nch, weights=w)),
            "n_channels_total": n_channels_total,
            "missing_fraction": float(1.0 - np.average(nch, weights=w) / n_channels_total),
            "abundance_log2_mean": float(np.average(g[abundance_col].to_numpy(dtype=float),
                                                    weights=w)),
            "n_peptides_total": float(g[peptide_col].fillna(0).sum()),
            "n_protein_rows": int(len(g)),
        })

    return df.groupby(gene_col, sort=False).apply(_agg).reset_index().rename(
        columns={gene_col: "gene_symbol"}
    )


def merge_observability_into_pool(
    pool: pd.DataFrame,
    obs: pd.DataFrame,
    pool_gene_col: str = "gene_symbol",
) -> pd.DataFrame:
    """Add per-gene ``missing_fraction``, ``n_channels_used_mean`` to a Module-2 pool."""
    out = pool.copy()
    out["gene_symbol_upper"] = out[pool_gene_col].astype(str).str.upper()
    obs2 = obs.copy()
    obs2["gene_symbol_upper"] = obs2["gene_symbol"].astype(str).str.upper()
    keep = ["gene_symbol_upper", "missing_fraction", "n_channels_used_mean",
            "n_channels_total"]
    out = out.merge(obs2[keep].drop_duplicates("gene_symbol_upper"),
                    on="gene_symbol_upper", how="left")
    out = out.drop(columns=["gene_symbol_upper"])
    return out


def assign_observability_strata(
    pool: pd.DataFrame,
    *,
    n_missing_bins: int = 3,
    missing_col: str = "missing_fraction",
    abundance_col: str = "abundance_log2",
    peptide_col: str = "n_peptides",
    n_abundance_bins: int = 5,
    n_peptide_bins: int = 4,
) -> pd.Series:
    """Extended (abundance × peptide × missing-fraction) joint stratum label.

    Falls back to the standard 5×4 strata if every row has the same
    ``missing_fraction`` (most TMT 2-plex scaled-S/N pools).  In that case
    the extra dimension is uninformative and the audit's Module-2 vs
    Module-3 q-value comparison should show no difference — itself a
    reportable finding ("the proteome has effectively zero missingness;
    detectability bias cannot explain the mismatch").
    """
    if missing_col not in pool.columns:
        raise KeyError(f"observability strata require '{missing_col}' column")
    miss = pd.to_numeric(pool[missing_col], errors="coerce").fillna(0.0)
    base = assign_match_strata(pool,
                               abundance_col=abundance_col,
                               peptide_col=peptide_col,
                               n_abundance_bins=n_abundance_bins,
                               n_peptide_bins=n_peptide_bins)
    if miss.nunique(dropna=True) <= 1:
        # No information in missingness — fall back to standard 5×4 strata so
        # the Module 2 / Module 3 q-value comparison is well-defined.
        return base
    try:
        miss_codes = pd.qcut(miss.rank(method="first"), q=n_missing_bins,
                             labels=False, duplicates="drop")
    except ValueError:
        miss_codes = pd.Series(0, index=miss.index)
    miss_suffix = "x0" + miss_codes.astype("Int64").astype(str)
    return (base.astype(str) + "_" + miss_suffix.astype(str)).reset_index(drop=True)



def detectability_gradient(
    rna_table: pd.DataFrame,
    protein_pool: pd.DataFrame,
    rna_col: str = "rrrm2_iss_t_rna_effect",
    n_bins: int = 5,
) -> pd.DataFrame:
    """Per-RNA-effect-magnitude-decile fraction of genes also protein-quantified.

    ``rna_table`` is the RNA universe (every gene with a finite RNA effect);
    ``protein_pool`` is the subset additionally protein-quantified.

    A monotone decline across bins would suggest large-RNA-effect genes
    are systematically detection-limited at the protein level (a real
    confounder).  A flat profile rules that out at the RNA-effect level.
    """
    df = rna_table.dropna(subset=[rna_col]).copy()
    df["abs_rna_effect"] = df[rna_col].abs()
    df["rna_bin"] = pd.qcut(df["abs_rna_effect"].rank(method="first"),
                            q=n_bins, labels=False, duplicates="drop")
    protein_ids = set(protein_pool["ENSEMBL"].astype(str).dropna())
    df["protein_quantified"] = df["ENSEMBL"].astype(str).isin(protein_ids)
    grouped = df.groupby("rna_bin", as_index=False).agg(
        n_rna=("ENSEMBL", "size"),
        n_protein=("protein_quantified", "sum"),
        mean_abs_rna_effect=("abs_rna_effect", "mean"),
        rna_effect_lo=("abs_rna_effect", "min"),
        rna_effect_hi=("abs_rna_effect", "max"),
    )
    grouped["fraction_protein_quantified"] = grouped["n_protein"] / grouped["n_rna"]
    return grouped



def high_coverage_subset(
    pool: pd.DataFrame,
    min_peptides: int = 3,
    max_missing_fraction: float = 0.20,
    peptide_col: str = "n_peptides",
    missing_col: str = "missing_fraction",
) -> pd.DataFrame:
    """Restrict to high-confidence quantification.

    Defaults (``min_peptides=3``, ``max_missing_fraction=0.2``) mirror the
    informal "high-coverage subset" used by reviewer-prep audits in the
    proteomics literature; both can be loosened or tightened for a
    sensitivity sweep.
    """
    df = pool.copy()
    df[peptide_col] = pd.to_numeric(df[peptide_col], errors="coerce")
    if missing_col not in df.columns:
        df[missing_col] = 0.0
    df[missing_col] = pd.to_numeric(df[missing_col], errors="coerce").fillna(0.0)
    keep = (df[peptide_col].fillna(0) >= min_peptides) & (df[missing_col] <= max_missing_fraction)
    sub = df.loc[keep].reset_index(drop=True)
    return sub



def propagation_with_strata(
    pool: pd.DataFrame,
    gene_sets: Mapping[str, list[str]],
    *,
    stratum_col: str,
    n_null: int = 10000,
    seed: int = 20260607,
    rna_col: str = "osd462_rna_effect",
    protein_col: str = "protein_flight_effect",
    gene_col: str = "gene_upper",
    min_members: int = 3,
) -> pd.DataFrame:
    """Per-pathway matched-null propagation test using ``pool[stratum_col]``.

    Mirrors :func:`src.v11.rna_protein_propagation.compute_propagation_per_pathway`
    but draws strata from an arbitrary column on ``pool`` (this is how
    Module 3 swaps in the observability-extended stratum).
    """
    rng = np.random.default_rng(seed)
    pool = pool.reset_index(drop=True)
    if stratum_col != "match_stratum":
        pool["_audit_stratum"] = pool[stratum_col].values
        pool["match_stratum"] = pool["_audit_stratum"]
    rows: list[dict] = []
    for pathway, members in gene_sets.items():
        members_upper = {str(g).upper() for g in members}
        mask = pool[gene_col].isin(members_upper).to_numpy()
        n_q = int(mask.sum())
        base: dict = {"pathway": pathway, "n_quantified": n_q, "stratum": stratum_col}
        if n_q < min_members:
            base["skipped"] = True
            rows.append(base)
            continue
        sub = pool[mask]
        rna_mean = float(pd.to_numeric(sub[rna_col], errors="coerce").mean())
        pred_sign = float(np.sign(rna_mean)) or 1.0

        def f_slope(df: pd.DataFrame) -> float:
            return ols_slope(df[rna_col], df[protein_col])

        def f_signed_mean(df: pd.DataFrame, s: float = pred_sign) -> float:
            vals = pd.to_numeric(df[protein_col], errors="coerce").to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            return s * float(vals.mean()) if vals.size else float("nan")

        def f_sign_agree(df: pd.DataFrame) -> float:
            return sign_agreement_rate(df[rna_col], df[protein_col])

        slope_res = run_matched_null(pool, mask, f_slope, "protein_slope",
                                     n_null=n_null, rng=rng)
        signed_res = run_matched_null(pool, mask, f_signed_mean,
                                      "signed_mean_protein_by_rna",
                                      n_null=n_null, rng=rng)
        sign_res = run_matched_null(pool, mask, f_sign_agree, "protein_sign_agreement",
                                    n_null=n_null, rng=rng)

        base.update({
            "skipped": False,
            "rna_pathway_mean": rna_mean,
            "predicted_direction": "up" if pred_sign > 0 else "down",
            "protein_slope": slope_res.observed,
            "protein_slope_p_two_sided": slope_res.p_two_sided,
            "protein_slope_p_greater": slope_res.p_greater,
            "signed_mean_protein_by_rna": signed_res.observed,
            "signed_mean_p_greater": signed_res.p_greater,
            "signed_mean_p_two_sided": signed_res.p_two_sided,
            "protein_sign_agreement": sign_res.observed,
            "protein_sign_agreement_p_greater": sign_res.p_greater,
        })
        rows.append(base)
    return pd.DataFrame(rows)



def ncc_site_observability(
    phospho_all_sites: pd.DataFrame,
    *,
    site_position_col: str = "site_position",
    n_fl_col: str = "n_fl",
    n_gc_col: str = "n_gc",
    effect_col: str = "phospho_effect",
) -> pd.DataFrame:
    """Per-NCC/SPAK site: observability metrics + percentile vs the phosphoproteome.

    For each residue-indexed co-modified context feature (and the
    non-regulatory sentinels), report:
      - n_fl + n_gc (channels with a finite, positive scaled S/N),
      - missing_fraction percentile within the full phospho table,
      - intensity (effect magnitude) percentile,
      - an explicit role tag that does not imply isolated canonical occupancy.

    This is an observability audit only. It cannot upgrade position-only effect
    rows into isolated canonical-site evidence.
    """
    sites = phospho_all_sites.copy()
    sites["gene_upper"] = sites["gene_symbol"].astype(str).str.upper()
    sites[n_fl_col] = pd.to_numeric(sites[n_fl_col], errors="coerce")
    sites[n_gc_col] = pd.to_numeric(sites[n_gc_col], errors="coerce")
    sites[effect_col] = pd.to_numeric(sites[effect_col], errors="coerce")
    n_total = sites[n_fl_col].fillna(0) + sites[n_gc_col].fillna(0)
    sites["n_total_used"] = n_total
    sites["abs_effect"] = sites[effect_col].abs()
    sites["coverage_percentile"] = sites["n_total_used"].rank(pct=True) * 100
    sites["abs_effect_percentile"] = sites["abs_effect"].rank(pct=True) * 100

    def _key(s: object) -> str:
        return str(s).strip()

    site_rows = []
    for gene, pos in OSD462_COMODIFIED_CANONICAL_INDEX_FEATURES:
        sel = sites[(sites["gene_upper"].eq(str(gene).upper()))
                    & (sites[site_position_col].astype(str).map(_key).eq(_key(pos)))]
        role = "canonical-indexed co-modified context"
        for _, row in sel.iterrows():
            site_rows.append({
                "gene_symbol": gene, "site_position": pos, "role": role,
                "n_fl": row[n_fl_col], "n_gc": row[n_gc_col],
                "n_total_used": row["n_total_used"],
                "coverage_percentile": row["coverage_percentile"],
                "phospho_effect": row[effect_col],
                "abs_effect_percentile": row["abs_effect_percentile"],
            })
    for gene, pos in NCC_NONREGULATORY_SITES:
        sel = sites[(sites["gene_upper"].eq(str(gene).upper()))
                    & (sites[site_position_col].astype(str).map(_key).eq(_key(pos)))]
        role = "non-regulatory control"
        for _, row in sel.iterrows():
            site_rows.append({
                "gene_symbol": gene, "site_position": pos, "role": role,
                "n_fl": row[n_fl_col], "n_gc": row[n_gc_col],
                "n_total_used": row["n_total_used"],
                "coverage_percentile": row["coverage_percentile"],
                "phospho_effect": row[effect_col],
                "abs_effect_percentile": row["abs_effect_percentile"],
            })
    return pd.DataFrame(site_rows)


__all__ = [
    "assign_observability_strata",
    "collapse_observability_to_gene",
    "detectability_gradient",
    "high_coverage_subset",
    "merge_observability_into_pool",
    "ncc_site_observability",
    "propagation_with_strata",
]
