"""RNA-to-protein propagation scoring for the v11 layer-specificity extension."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from src.common import bh_fdr
from src.v11.matched_null import (
    ols_slope,
    prepare_matched_pool,
    run_matched_null,
    sign_agreement_rate,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RUN_ROOT = REPO_ROOT / "data/results/run_20260606_v11_layer_specificity"
DEFAULT_ANCHOR_DIR = REPO_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor"
DEFAULT_GENE_SETS = REPO_ROOT / "config/mechanism_gene_sets.yaml"
DEFAULT_SEED = 20260606


@dataclass(frozen=True)
class PropagationConfig:
    """Configuration for pathway-level RNA-to-protein propagation scoring."""

    run_root: Path = DEFAULT_RUN_ROOT
    anchor_dir: Path = DEFAULT_ANCHOR_DIR
    gene_sets_path: Path = DEFAULT_GENE_SETS
    n_null: int = 10000
    n_bootstrap: int = 2000
    seed: int = DEFAULT_SEED
    peptide_filter: int = 2
    rna_threshold: float = 0.04
    protein_threshold: float = 0.02
    phospho_threshold: float = 0.02
    min_coobserved: int = 3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def upper_symbol(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.upper()


def direction(value: float, threshold: float) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        value = float("nan")
    if not np.isfinite(value):
        return "unobserved"
    if abs(value) < threshold:
        return "flat"
    return "up" if value > 0 else "down"


def load_gene_sets(path: Path) -> dict[str, list[str]]:
    """Load curated mechanism sets as pathway -> upper-case symbols."""
    cfg = yaml.safe_load(path.read_text())
    out: dict[str, list[str]] = {}
    for name, spec in cfg.items():
        genes = spec.get("genes", []) if isinstance(spec, dict) else spec
        out[str(name)] = sorted({str(g).strip().upper() for g in genes if str(g).strip()})
    return out


def phospho_parent_effects(phospho: pd.DataFrame) -> pd.DataFrame:
    """Collapse phosphosite effects to parent-gene means for layer assignment."""
    df = phospho.copy()
    df["gene_upper"] = upper_symbol(df["gene_symbol"])
    df["phospho_effect"] = pd.to_numeric(df["phospho_effect"], errors="coerce")
    if "phospho_p_value" in df.columns:
        df["phospho_p_value"] = pd.to_numeric(df["phospho_p_value"], errors="coerce")
    else:
        df["phospho_p_value"] = np.nan
    return (
        df.dropna(subset=["gene_upper"])
        .groupby("gene_upper", as_index=False)
        .agg(
            phospho_parent_effect=("phospho_effect", "mean"),
            n_phosphosites=("phospho_effect", "count"),
            min_phospho_p_value=("phospho_p_value", "min"),
        )
    )


def build_layer_table(master: pd.DataFrame, phospho: pd.DataFrame) -> pd.DataFrame:
    """Join harmonized gene-level RNA/protein effects to parent-gene phospho effects."""
    df = master.copy()
    df["gene_upper"] = upper_symbol(df["gene_symbol"])
    numeric = [
        "rrrm2_iss_t_rna_effect",
        "osd462_rna_effect",
        "osd462_rna_adj_p",
        "protein_flight_effect",
        "n_peptides",
        "plex_coverage",
        "n_protein_rows",
        "abundance_log2",
    ]
    for col in numeric:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    parent = phospho_parent_effects(phospho)
    df = df.merge(parent, on="gene_upper", how="left")
    return df


def classify_gene(
    row: pd.Series,
    rna_threshold: float = 0.04,
    protein_threshold: float = 0.02,
    phospho_threshold: float = 0.02,
) -> tuple[str, str, str, str]:
    """Return ternary class, detail class, and layer directions for one gene."""
    rna_dir = direction(row.get("osd462_rna_effect", np.nan), rna_threshold)
    protein_dir = direction(row.get("protein_flight_effect", np.nan), protein_threshold)
    phospho_dir = direction(row.get("phospho_parent_effect", np.nan), phospho_threshold)

    if rna_dir == "unobserved":
        return "rna_unobserved", "rna_unobserved", protein_dir, phospho_dir
    if rna_dir == "flat":
        return "rna_flat", "rna_flat", protein_dir, phospho_dir

    protein_concordant = protein_dir in {"up", "down"} and protein_dir == rna_dir
    protein_discordant = protein_dir in {"up", "down"} and protein_dir != rna_dir
    phospho_concordant = phospho_dir in {"up", "down"} and phospho_dir == rna_dir

    if protein_concordant:
        return "rna_to_protein", "rna_to_protein", protein_dir, phospho_dir
    if phospho_concordant:
        detail = "rna_to_phospho_with_protein_discordance" if protein_discordant else "rna_to_phospho"
        return "rna_to_phospho", detail, protein_dir, phospho_dir
    if protein_discordant:
        return "rna_only", "rna_protein_discordant", protein_dir, phospho_dir
    return "rna_only", "rna_only_or_protein_unobserved", protein_dir, phospho_dir


def pathway_gene_table(
    layer: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    cfg: PropagationConfig,
) -> pd.DataFrame:
    """One row per pathway member resolved through the harmonized OSD-462 table."""
    rows = []
    by_symbol = layer.dropna(subset=["gene_upper"]).copy()
    for pathway, symbols in gene_sets.items():
        seen: set[tuple[str, str]] = set()
        for symbol in symbols:
            hits = by_symbol[by_symbol["gene_upper"].eq(symbol)]
            if hits.empty:
                rows.append(
                    {
                        "pathway": pathway,
                        "query_symbol": symbol,
                        "ENSEMBL": pd.NA,
                        "gene_symbol": pd.NA,
                        "gene_upper": symbol,
                        "resolved": False,
                    }
                )
                continue
            for _, hit in hits.iterrows():
                key = (pathway, str(hit["ENSEMBL"]))
                if key in seen:
                    continue
                seen.add(key)
                row = {
                    "pathway": pathway,
                    "query_symbol": symbol,
                    "ENSEMBL": hit["ENSEMBL"],
                    "gene_symbol": hit["gene_symbol"],
                    "gene_upper": hit["gene_upper"],
                    "resolved": True,
                    "osd462_rna_effect": hit.get("osd462_rna_effect", np.nan),
                    "rrrm2_iss_t_rna_effect": hit.get("rrrm2_iss_t_rna_effect", np.nan),
                    "protein_flight_effect": hit.get("protein_flight_effect", np.nan),
                    "phospho_parent_effect": hit.get("phospho_parent_effect", np.nan),
                    "n_peptides": hit.get("n_peptides", np.nan),
                    "abundance_log2": hit.get("abundance_log2", np.nan),
                    "n_phosphosites": hit.get("n_phosphosites", np.nan),
                }
                ternary, detail, protein_dir, phospho_dir = classify_gene(
                    pd.Series(row),
                    rna_threshold=cfg.rna_threshold,
                    protein_threshold=cfg.protein_threshold,
                    phospho_threshold=cfg.phospho_threshold,
                )
                row["rna_dir"] = direction(row["osd462_rna_effect"], cfg.rna_threshold)
                row["protein_dir"] = protein_dir
                row["phospho_dir"] = phospho_dir
                row["ternary_class"] = ternary
                row["detail_class"] = detail
                rows.append(row)
    out = pd.DataFrame(rows)
    for col in [
        "osd462_rna_effect",
        "rrrm2_iss_t_rna_effect",
        "protein_flight_effect",
        "phospho_parent_effect",
        "n_peptides",
        "abundance_log2",
        "n_phosphosites",
    ]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out


def bootstrap_fraction_ci(
    classes: pd.Series,
    label: str,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    values = classes.astype(str).to_numpy()
    if values.size == 0 or n_bootstrap <= 0:
        return float("nan"), float("nan")
    boot = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        draw = rng.choice(values, size=values.size, replace=True)
        boot[i] = np.mean(draw == label)
    return float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))


def _safe_ratio(numer: int, denom: int) -> float:
    return float(numer / denom) if denom else float("nan")


def pathway_assignment(row: pd.Series) -> str:
    """Conservative pathway-level layer assignment from calibrated quantities."""
    signed_mean = float(row.get("signed_mean_protein_by_rna", np.nan))
    inverse_mean_q = float(row.get("protein_inverse_signed_mean_q_greater", np.nan))
    slope = float(row.get("protein_slope", np.nan))
    slope_q = float(row.get("protein_slope_q_greater", np.nan))
    inverse_q = float(row.get("protein_inverse_slope_q_greater", np.nan))
    if np.isfinite(signed_mean) and signed_mean < 0 and np.isfinite(inverse_mean_q) and inverse_mean_q < 0.10:
        return "protein_inverted_calibrated"
    if np.isfinite(slope) and slope > 0 and np.isfinite(slope_q) and slope_q < 0.10:
        return "RNA_to_protein_calibrated"
    if np.isfinite(slope) and slope < 0 and np.isfinite(inverse_q) and inverse_q < 0.10:
        return "protein_inverted_calibrated"
    if row.get("n_rna_phospho", 0) >= 2 and float(row.get("fraction_rna_to_phospho", 0.0)) >= 0.25:
        return "RNA_to_phospho_candidate"
    return "RNA_not_propagated"


def summarize_pathways(
    layer: pd.DataFrame,
    genes: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    cfg: PropagationConfig,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    pool = prepare_matched_pool(
        layer,
        required_cols=["osd462_rna_effect", "protein_flight_effect"],
        peptide_filter=cfg.peptide_filter,
    )
    rng = np.random.default_rng(cfg.seed)
    rows = []
    coverage_rows = []

    for pathway, symbols in gene_sets.items():
        part = genes[genes["pathway"].eq(pathway)].copy()
        resolved = part[part["resolved"].fillna(False)].copy()
        rna_directional = resolved[resolved["rna_dir"].isin(["up", "down"])]
        denom = int(len(rna_directional))

        n_rna_to_protein = int((rna_directional["ternary_class"] == "rna_to_protein").sum())
        n_rna_to_phospho = int((rna_directional["ternary_class"] == "rna_to_phospho").sum())
        n_rna_only = int((rna_directional["ternary_class"] == "rna_only").sum())
        n_discordant = int((rna_directional["detail_class"] == "rna_protein_discordant").sum())

        pathway_seed = int(hashlib.sha256(pathway.encode("utf-8")).hexdigest()[:8], 16)
        ci_rng = np.random.default_rng(cfg.seed + pathway_seed % 100000)
        protein_ci = bootstrap_fraction_ci(rna_directional["ternary_class"], "rna_to_protein", cfg.n_bootstrap, ci_rng)
        phospho_ci = bootstrap_fraction_ci(rna_directional["ternary_class"], "rna_to_phospho", cfg.n_bootstrap, ci_rng)
        rna_only_ci = bootstrap_fraction_ci(rna_directional["ternary_class"], "rna_only", cfg.n_bootstrap, ci_rng)

        coobserved_genes = set(
            resolved.loc[
                np.isfinite(resolved["osd462_rna_effect"]) & np.isfinite(resolved["protein_flight_effect"]),
                "ENSEMBL",
            ].astype(str)
        )
        target_mask = pool["ENSEMBL"].astype(str).isin(coobserved_genes).to_numpy() if len(pool) else np.array([])
        n_coobserved = int(target_mask.sum()) if len(pool) else 0

        row = {
            "pathway": pathway,
            "n_set": int(len(symbols)),
            "n_resolved": int(resolved["ENSEMBL"].nunique()),
            "n_rna": int(np.isfinite(resolved["osd462_rna_effect"]).sum()),
            "n_protein": int(np.isfinite(resolved["protein_flight_effect"]).sum()),
            "n_rna_protein": n_coobserved,
            "n_phospho_parent": int(np.isfinite(resolved["phospho_parent_effect"]).sum()),
            "n_rna_phospho": int(
                (
                    np.isfinite(resolved["osd462_rna_effect"])
                    & np.isfinite(resolved["phospho_parent_effect"])
                ).sum()
            ),
            "n_rna_directional": denom,
            "n_rna_to_protein": n_rna_to_protein,
            "n_rna_to_phospho": n_rna_to_phospho,
            "n_rna_only": n_rna_only,
            "n_rna_protein_discordant": n_discordant,
            "fraction_rna_to_protein": _safe_ratio(n_rna_to_protein, denom),
            "fraction_rna_to_protein_ci_low": protein_ci[0],
            "fraction_rna_to_protein_ci_high": protein_ci[1],
            "fraction_rna_to_phospho": _safe_ratio(n_rna_to_phospho, denom),
            "fraction_rna_to_phospho_ci_low": phospho_ci[0],
            "fraction_rna_to_phospho_ci_high": phospho_ci[1],
            "fraction_rna_only": _safe_ratio(n_rna_only, denom),
            "fraction_rna_only_ci_low": rna_only_ci[0],
            "fraction_rna_only_ci_high": rna_only_ci[1],
            "fraction_rna_protein_discordant": _safe_ratio(n_discordant, denom),
        }

        if n_coobserved >= cfg.min_coobserved:
            observed_pool = pool[target_mask]
            mean_rna_effect = float(observed_pool["osd462_rna_effect"].mean())
            mean_protein_effect = float(observed_pool["protein_flight_effect"].mean())
            pathway_rna_sign = float(np.sign(mean_rna_effect)) or 1.0

            def slope_stat(df: pd.DataFrame) -> float:
                return ols_slope(df["osd462_rna_effect"], df["protein_flight_effect"], min_n=cfg.min_coobserved)

            def inverse_slope_stat(df: pd.DataFrame) -> float:
                slope = slope_stat(df)
                return -slope if np.isfinite(slope) else float("nan")

            def signed_mean_stat(df: pd.DataFrame, sign: float = pathway_rna_sign) -> float:
                return sign * float(df["protein_flight_effect"].mean())

            def inverse_signed_mean_stat(df: pd.DataFrame) -> float:
                val = signed_mean_stat(df)
                return -val if np.isfinite(val) else float("nan")

            def sign_stat(df: pd.DataFrame) -> float:
                return sign_agreement_rate(df["osd462_rna_effect"], df["protein_flight_effect"], min_n=cfg.min_coobserved)

            slope_null = run_matched_null(
                pool,
                target_mask,
                slope_stat,
                "protein_slope",
                n_null=cfg.n_null,
                rng=rng,
            )
            inverse_null = run_matched_null(
                pool,
                target_mask,
                inverse_slope_stat,
                "protein_inverse_slope",
                n_null=cfg.n_null,
                rng=rng,
            )
            signed_mean_null = run_matched_null(
                pool,
                target_mask,
                signed_mean_stat,
                "protein_signed_mean_by_rna_direction",
                n_null=cfg.n_null,
                rng=rng,
            )
            inverse_signed_mean_null = run_matched_null(
                pool,
                target_mask,
                inverse_signed_mean_stat,
                "protein_inverse_signed_mean_by_rna_direction",
                n_null=cfg.n_null,
                rng=rng,
            )
            sign_null = run_matched_null(
                pool,
                target_mask,
                sign_stat,
                "protein_sign_agreement",
                n_null=cfg.n_null,
                rng=rng,
            )
            row.update(
                {
                    "mean_rna_effect_coobserved": mean_rna_effect,
                    "mean_protein_effect_coobserved": mean_protein_effect,
                    "signed_mean_protein_by_rna": signed_mean_null.observed,
                    "signed_mean_protein_null_median": signed_mean_null.null_median,
                    "signed_mean_protein_null_ci_low": signed_mean_null.null_ci_low,
                    "signed_mean_protein_null_ci_high": signed_mean_null.null_ci_high,
                    "protein_signed_mean_p_greater": signed_mean_null.p_greater,
                    "protein_signed_mean_p_two_sided": signed_mean_null.p_two_sided,
                    "protein_inverse_signed_mean_p_greater": inverse_signed_mean_null.p_greater,
                    "protein_slope": slope_null.observed,
                    "protein_slope_null_median": slope_null.null_median,
                    "protein_slope_null_ci_low": slope_null.null_ci_low,
                    "protein_slope_null_ci_high": slope_null.null_ci_high,
                    "protein_slope_p_greater": slope_null.p_greater,
                    "protein_slope_p_two_sided": slope_null.p_two_sided,
                    "protein_inverse_slope_p_greater": inverse_null.p_greater,
                    "protein_sign_agreement": sign_null.observed,
                    "protein_sign_agreement_null_median": sign_null.null_median,
                    "protein_sign_agreement_p_greater": sign_null.p_greater,
                    "n_null_valid": slope_null.n_null_valid,
                }
            )
        else:
            row.update(
                {
                    "mean_rna_effect_coobserved": np.nan,
                    "mean_protein_effect_coobserved": np.nan,
                    "signed_mean_protein_by_rna": np.nan,
                    "signed_mean_protein_null_median": np.nan,
                    "signed_mean_protein_null_ci_low": np.nan,
                    "signed_mean_protein_null_ci_high": np.nan,
                    "protein_signed_mean_p_greater": np.nan,
                    "protein_signed_mean_p_two_sided": np.nan,
                    "protein_inverse_signed_mean_p_greater": np.nan,
                    "protein_slope": np.nan,
                    "protein_slope_null_median": np.nan,
                    "protein_slope_null_ci_low": np.nan,
                    "protein_slope_null_ci_high": np.nan,
                    "protein_slope_p_greater": np.nan,
                    "protein_slope_p_two_sided": np.nan,
                    "protein_inverse_slope_p_greater": np.nan,
                    "protein_sign_agreement": np.nan,
                    "protein_sign_agreement_null_median": np.nan,
                    "protein_sign_agreement_p_greater": np.nan,
                    "n_null_valid": 0,
                }
            )
        rows.append(row)
        coverage_rows.append({k: row[k] for k in [
            "pathway",
            "n_set",
            "n_resolved",
            "n_rna",
            "n_protein",
            "n_rna_protein",
            "n_phospho_parent",
            "n_rna_phospho",
        ]})

    summary = pd.DataFrame(rows)
    summary["protein_propagation_score"] = summary["protein_slope"]
    for p_col, q_col in [
        ("protein_slope_p_greater", "protein_slope_q_greater"),
        ("protein_inverse_slope_p_greater", "protein_inverse_slope_q_greater"),
        ("protein_signed_mean_p_greater", "protein_signed_mean_q_greater"),
        ("protein_inverse_signed_mean_p_greater", "protein_inverse_signed_mean_q_greater"),
        ("protein_sign_agreement_p_greater", "protein_sign_agreement_q_greater"),
    ]:
        summary[q_col] = np.nan
        ok = np.isfinite(pd.to_numeric(summary[p_col], errors="coerce"))
        if ok.any():
            summary.loc[ok, q_col] = bh_fdr(summary.loc[ok, p_col].to_numpy(dtype=float))
    summary["layer_assignment"] = summary.apply(pathway_assignment, axis=1)
    summary = summary.sort_values(["protein_slope", "pathway"], ascending=[True, True]).reset_index(drop=True)
    coverage = pd.DataFrame(coverage_rows).sort_values("pathway").reset_index(drop=True)
    return summary, coverage


def plot_propagation(summary: pd.DataFrame, out_dir: Path) -> None:
    """Write a compact propagation figure for manuscript/supplement review."""
    mpl_cache = out_dir.parent / "logs" / "matplotlib"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))
    os.environ.setdefault("XDG_CACHE_HOME", str(out_dir.parent / "logs" / "cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir.mkdir(parents=True, exist_ok=True)
    df = summary.sort_values("signed_mean_protein_by_rna", ascending=True).reset_index(drop=True)
    y = np.arange(len(df))

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.8), gridspec_kw={"width_ratios": [1.15, 1.0]})
    ax = axes[0]
    ax.axvline(0, color="#555555", lw=1)
    lo = df["signed_mean_protein_null_ci_low"].to_numpy(dtype=float)
    hi = df["signed_mean_protein_null_ci_high"].to_numpy(dtype=float)
    signed_mean = df["signed_mean_protein_by_rna"].to_numpy(dtype=float)
    ax.hlines(y, lo, hi, color="#b8b8b8", lw=3, label="matched null 95% interval")
    colors = np.where(signed_mean < 0, "#b24a3b", "#2f6f9f")
    ax.scatter(signed_mean, y, c=colors, s=48, zorder=3, edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(df["pathway"].str.replace("_", " "), fontsize=8)
    ax.set_xlabel("Protein effect signed by RNA direction")
    ax.set_title("Directional Protein Carry-Through", fontsize=11)
    ax.grid(axis="x", color="#e5e5e5", lw=0.8)

    ax2 = axes[1]
    stack_cols = [
        ("fraction_rna_to_protein", "RNA to protein", "#2f6f9f"),
        ("fraction_rna_to_phospho", "RNA to phospho", "#5c8f48"),
        ("fraction_rna_protein_discordant", "RNA/protein inverted", "#b24a3b"),
        ("fraction_rna_only", "RNA only/other", "#b7a36a"),
    ]
    left = np.zeros(len(df), dtype=float)
    for col, label, color in stack_cols:
        vals = df[col].fillna(0).to_numpy(dtype=float)
        ax2.barh(y, vals, left=left, color=color, height=0.72, label=label)
        left += vals
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xlim(0, 1)
    ax2.set_xlabel("Fraction of RNA-directional pathway genes")
    ax2.set_title("Per-Gene Layer Class", fontsize=11)
    ax2.grid(axis="x", color="#e5e5e5", lw=0.8)
    ax2.legend(loc="lower center", bbox_to_anchor=(0.5, -0.25), ncol=2, fontsize=8, frameon=False)

    fig.tight_layout()
    fig.savefig(out_dir / "v11_rna_protein_propagation.png", dpi=220, bbox_inches="tight")
    fig.savefig(out_dir / "v11_rna_protein_propagation.pdf", bbox_inches="tight")
    plt.close(fig)


def write_manifest(cfg: PropagationConfig, outputs: dict[str, Path], row_counts: dict[str, int]) -> Path:
    manifest_dir = cfg.run_root / "manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    inputs = {
        "osd462_flight_effects": cfg.anchor_dir / "osd462_flight_effects.tsv",
        "phospho_all_sites": cfg.anchor_dir / "phospho_all_sites.tsv",
        "mechanism_gene_sets": cfg.gene_sets_path,
    }
    manifest = {
        "analysis": "v11 RNA-to-protein propagation score",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "parameters": {
            "n_null": cfg.n_null,
            "n_bootstrap": cfg.n_bootstrap,
            "seed": cfg.seed,
            "peptide_filter": cfg.peptide_filter,
            "rna_threshold": cfg.rna_threshold,
            "protein_threshold": cfg.protein_threshold,
            "phospho_threshold": cfg.phospho_threshold,
            "min_coobserved": cfg.min_coobserved,
            "null_model": "abundance x peptide-count matched random gene sets",
        },
        "inputs": {},
        "outputs": {name: str(path) for name, path in outputs.items()},
        "row_counts": row_counts,
    }
    for name, path in inputs.items():
        manifest["inputs"][name] = str(path)
        if path.exists():
            manifest["inputs"][f"{name}_sha256"] = sha256(path)
    path = manifest_dir / "v11_rna_protein_propagation_manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path


def run_propagation(cfg: PropagationConfig) -> dict[str, Path]:
    prop_dir = cfg.run_root / "propagation"
    fig_dir = cfg.run_root / "figures"
    prop_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    master_path = cfg.anchor_dir / "osd462_flight_effects.tsv"
    phospho_path = cfg.anchor_dir / "phospho_all_sites.tsv"
    master = read_tsv(master_path)
    phospho = read_tsv(phospho_path)
    gene_sets = load_gene_sets(cfg.gene_sets_path)
    layer = build_layer_table(master, phospho)
    genes = pathway_gene_table(layer, gene_sets, cfg)
    summary, coverage = summarize_pathways(layer, genes, gene_sets, cfg)

    summary_path = prop_dir / "rna_protein_propagation_summary.tsv"
    genes_path = prop_dir / "rna_protein_propagation_gene_classes.tsv"
    coverage_path = prop_dir / "rna_protein_propagation_coverage.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    genes.to_csv(genes_path, sep="\t", index=False)
    coverage.to_csv(coverage_path, sep="\t", index=False)
    plot_propagation(summary, fig_dir)

    outputs = {
        "summary": summary_path,
        "gene_classes": genes_path,
        "coverage": coverage_path,
        "figure_png": fig_dir / "v11_rna_protein_propagation.png",
        "figure_pdf": fig_dir / "v11_rna_protein_propagation.pdf",
    }
    write_manifest(
        cfg,
        outputs=outputs,
        row_counts={
            "master_rows": int(len(master)),
            "phosphosite_rows": int(len(phospho)),
            "pathway_member_rows": int(len(genes)),
            "summary_rows": int(len(summary)),
        },
    )
    return outputs


__all__ = [
    "PropagationConfig",
    "build_layer_table",
    "classify_gene",
    "direction",
    "load_gene_sets",
    "pathway_gene_table",
    "phospho_parent_effects",
    "run_propagation",
    "summarize_pathways",
]
