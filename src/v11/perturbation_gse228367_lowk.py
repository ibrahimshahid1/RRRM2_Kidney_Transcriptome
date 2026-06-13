#!/usr/bin/env python3
"""GSE228367 low-potassium DCT perturbation alignment."""

from __future__ import annotations

import argparse
import gc
import json
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_ROOT = REPO_ROOT / "data/results/run_20260526_v11_dct1_phospho_mediation"
GSE_DIR = REPO_ROOT / "data/external/dct_reference/GSE228367"
RAW_TAR = GSE_DIR / "GSE228367_RAW.tar"
OSD462_EFFECTS = REPO_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_flight_effects.tsv"
RRRM2_GENE = (
    REPO_ROOT
    / "data/results/run_20260519_000547_2500g/contrast_vectors/mechanism_axis/"
    "tubulointerstitial_state/lar_reversal/lar_reversal_gene_scatter.tsv"
)
OSD513_GENE = REPO_ROOT / "data/results/run_20260522_regulator_activity/rna_effects/OSD513_gene_stat.tsv"

TARGET_GENES = [
    "Slc12a3",
    "Stk39",
    "Oxsr1",
    "Wnk1",
    "Wnk4",
    "Kcnj10",
    "Kcnj16",
    "Clcnkb",
    "Bsnd",
    "Klhl3",
    "Cul3",
    "Ppp1ca",
    "Ppp1r1a",
    "Calb1",
]


def bh(pvals) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.ones_like(p)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(vals)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[ok] = restored
    return out


def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    mu = np.nanmean(v)
    sd = np.nanstd(v)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(v) * np.nan
    return (v - mu) / sd


def cosine(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 3:
        return float("nan")
    x = x[keep]
    y = y[keep]
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    if denom <= 0:
        return float("nan")
    return float(np.dot(x, y) / denom)


def extract_filtered_h5(out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    if not RAW_TAR.exists():
        raise FileNotFoundError(
            f"Missing {RAW_TAR}. Download from "
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228367/suppl/GSE228367_RAW.tar"
        )
    h5_paths: list[Path] = []
    with tarfile.open(RAW_TAR, "r") as tf:
        members = [
            m
            for m in tf.getmembers()
            if m.name.endswith("_filtered_feature_bc_matrix.h5")
        ]
        for member in members:
            out = out_dir / Path(member.name).name
            if not out.exists() or out.stat().st_size == 0:
                src = tf.extractfile(member)
                if src is None:
                    raise FileNotFoundError(member.name)
                with out.open("wb") as fh:
                    while True:
                        chunk = src.read(1 << 20)
                        if not chunk:
                            break
                        fh.write(chunk)
            h5_paths.append(out)
    return sorted(h5_paths)


def parse_sample(path: Path) -> dict:
    m = re.match(r"^(GSM\d+)_(NK|KD)(\d+)_filtered_feature_bc_matrix\.h5$", path.name)
    if not m:
        raise ValueError(f"Cannot parse GSE228367 sample name: {path.name}")
    diet = "normal_K" if m.group(2) == "NK" else "low_K"
    return {"gsm": m.group(1), "sample": f"{m.group(2)}{m.group(3)}", "diet": diet, "replicate": int(m.group(3))}


def sample_mean_expression(h5_path: Path) -> pd.Series:
    import scanpy as sc

    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    mean = np.asarray(adata.X.mean(axis=0)).ravel()
    out = pd.Series(mean, index=adata.var_names.astype(str), name=h5_path.name)
    del adata
    gc.collect()
    return out


def build_lowk_effects(out_dir: Path) -> pd.DataFrame:
    mean_path = out_dir / "lowk_replicate_mean_expression.tsv.gz"
    if mean_path.exists():
        means = pd.read_csv(mean_path, sep="\t", index_col=0)
    else:
        h5_paths = extract_filtered_h5(out_dir / "gse228367_filtered_h5")
        means = []
        meta_rows = []
        for h5 in h5_paths:
            info = parse_sample(h5)
            print(f"[lowK] reading {h5.name}")
            s = sample_mean_expression(h5)
            s.name = info["sample"]
            means.append(s)
            meta_rows.append(info | {"file": str(h5)})
        means = pd.concat(means, axis=1)
        means.to_csv(mean_path, sep="\t", compression="gzip")
        pd.DataFrame(meta_rows).to_csv(out_dir / "lowk_sample_manifest.tsv", sep="\t", index=False)

    nk_cols = [c for c in means.columns if str(c).startswith("NK")]
    kd_cols = [c for c in means.columns if str(c).startswith("KD")]
    if len(nk_cols) < 2 or len(kd_cols) < 2:
        raise ValueError(f"Expected at least two NK and KD columns; got NK={nk_cols}, KD={kd_cols}")

    nk = means[nk_cols]
    kd = means[kd_cols]
    eff = kd.mean(axis=1) - nk.mean(axis=1)
    stat, p = stats.ttest_ind(kd.T, nk.T, equal_var=False, nan_policy="omit")
    out = pd.DataFrame(
        {
            "gene_symbol": means.index.astype(str),
            "lowk_effect_kd_minus_nk": eff.to_numpy(dtype=float),
            "lowk_t_stat": stat,
            "lowk_p_value": p,
            "lowk_mean_kd": kd.mean(axis=1).to_numpy(dtype=float),
            "lowk_mean_nk": nk.mean(axis=1).to_numpy(dtype=float),
            "lowk_n_kd_reps": len(kd_cols),
            "lowk_n_nk_reps": len(nk_cols),
        }
    )
    out["lowk_q_value"] = bh(out["lowk_p_value"])
    out.to_csv(out_dir / "lowk_dct_pseudobulk_effects.tsv", sep="\t", index=False)
    return out


def load_spaceflight_vectors() -> dict[str, pd.Series]:
    vectors: dict[str, pd.Series] = {}
    if OSD462_EFFECTS.exists():
        osd = pd.read_csv(OSD462_EFFECTS, sep="\t")
        vectors["osd462_rna"] = (
            osd.dropna(subset=["gene_symbol", "osd462_rna_effect"])
            .drop_duplicates("gene_symbol")
            .set_index("gene_symbol")["osd462_rna_effect"]
        )
    if RRRM2_GENE.exists():
        rr = pd.read_csv(RRRM2_GENE, sep="\t")
        vectors["rrrm2_iss_t_rna"] = (
            rr.dropna(subset=["mgi_symbol", "iss_t_effect"])
            .drop_duplicates("mgi_symbol")
            .set_index("mgi_symbol")["iss_t_effect"]
        )
    if OSD513_GENE.exists():
        osd513 = pd.read_csv(OSD513_GENE, sep="\t")
        vectors["osd513_rna_stat"] = (
            osd513.dropna(subset=["gene", "stat"])
            .drop_duplicates("gene")
            .set_index("gene")["stat"]
        )
    return vectors


def load_gene_sets(run_root: Path) -> dict[str, set[str]]:
    prior = pd.read_csv(run_root / "dct_prior/gse228367_dct1_vs_dct2_de.tsv", sep="\t")
    score = pd.to_numeric(prior["dct1_enrichment_score"], errors="coerce")
    top_decile = score.quantile(0.90)
    bottom_decile = score.quantile(0.10)
    top_quartile = score.quantile(0.75)
    bottom_quartile = score.quantile(0.25)
    sets = {
        "all_overlap": set(prior["gene_symbol"].astype(str)),
        "dct1_core_genes": set(prior.loc[prior["dct_expression_class"].eq("DCT1_core"), "gene_symbol"].astype(str)),
        "dct2_core_genes": set(prior.loc[prior["dct_expression_class"].eq("DCT2_core"), "gene_symbol"].astype(str)),
        "dct_shared_genes": set(prior.loc[prior["dct_expression_class"].eq("DCT_shared"), "gene_symbol"].astype(str)),
        "dct1_top_decile_genes": set(prior.loc[score >= top_decile, "gene_symbol"].astype(str)),
        "dct2_top_decile_genes": set(prior.loc[score <= bottom_decile, "gene_symbol"].astype(str)),
        "dct1_top_quartile_genes": set(prior.loc[score >= top_quartile, "gene_symbol"].astype(str)),
        "dct2_top_quartile_genes": set(prior.loc[score <= bottom_quartile, "gene_symbol"].astype(str)),
        "transport_target_genes": set(TARGET_GENES),
    }
    return sets


def alignment_rows(lowk: pd.DataFrame, vectors: dict[str, pd.Series], gene_sets: dict[str, set[str]], n_boot: int) -> pd.DataFrame:
    rng = np.random.default_rng(20260527)
    low = lowk.dropna(subset=["lowk_effect_kd_minus_nk"]).drop_duplicates("gene_symbol").set_index("gene_symbol")
    rows = []
    for vector_name, vec in vectors.items():
        joined = low[["lowk_effect_kd_minus_nk"]].join(vec.rename("flight_effect"), how="inner").dropna()
        for set_name, genes in gene_sets.items():
            sub = joined.loc[joined.index.intersection(genes)].copy()
            if len(sub) < 10:
                continue
            x = sub["lowk_effect_kd_minus_nk"].to_numpy(dtype=float)
            y = sub["flight_effect"].to_numpy(dtype=float)
            c = cosine(x, y)
            rho, rho_p = stats.spearmanr(x, y, nan_policy="omit")
            boot_cos = []
            boot_rho = []
            for _ in range(n_boot):
                idx = rng.integers(0, len(sub), len(sub))
                bx = x[idx]
                by = y[idx]
                boot_cos.append(cosine(bx, by))
                try:
                    boot_rho.append(stats.spearmanr(bx, by, nan_policy="omit").statistic)
                except Exception:
                    boot_rho.append(np.nan)
            rows.append(
                {
                    "spaceflight_vector": vector_name,
                    "gene_subset": set_name,
                    "n_genes": len(sub),
                    "cosine": c,
                    "cosine_ci_low": float(np.nanpercentile(boot_cos, 2.5)),
                    "cosine_ci_high": float(np.nanpercentile(boot_cos, 97.5)),
                    "spearman_rho": float(rho),
                    "spearman_p_value": float(rho_p),
                    "spearman_ci_low": float(np.nanpercentile(boot_rho, 2.5)),
                    "spearman_ci_high": float(np.nanpercentile(boot_rho, 97.5)),
                    "interpretation": (
                        "anti_aligned" if c < -0.1 else "aligned" if c > 0.1 else "weak_or_null"
                    ),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["spearman_q_value"] = bh(out["spearman_p_value"])
    return out


def target_gene_table(lowk: pd.DataFrame, vectors: dict[str, pd.Series], run_root: Path) -> pd.DataFrame:
    prior = pd.read_csv(run_root / "dct_prior/gse228367_dct1_vs_dct2_de.tsv", sep="\t")
    prior = prior.drop_duplicates("gene_symbol").set_index("gene_symbol")
    low = lowk.drop_duplicates("gene_symbol").set_index("gene_symbol")
    rows = []
    for gene in TARGET_GENES:
        row = {
            "gene_symbol": gene,
            "lowk_effect_kd_minus_nk": low.get("lowk_effect_kd_minus_nk", pd.Series(dtype=float)).get(gene, np.nan),
            "lowk_p_value": low.get("lowk_p_value", pd.Series(dtype=float)).get(gene, np.nan),
            "lowk_q_value": low.get("lowk_q_value", pd.Series(dtype=float)).get(gene, np.nan),
            "dct_expression_class": prior.get("dct_expression_class", pd.Series(dtype=str)).get(gene, "not_in_prior"),
            "dct1_enrichment_score": prior.get("dct1_enrichment_score", pd.Series(dtype=float)).get(gene, np.nan),
        }
        for name, vec in vectors.items():
            row[name] = vec.get(gene, np.nan)
            if np.isfinite(row["lowk_effect_kd_minus_nk"]) and np.isfinite(row[name]):
                prod = row["lowk_effect_kd_minus_nk"] * row[name]
                row[f"{name}_vs_lowk_direction"] = "opposite" if prod < 0 else "same" if prod > 0 else "zero"
            else:
                row[f"{name}_vs_lowk_direction"] = "missing"
        rows.append(row)
    return pd.DataFrame(rows)


def write_verdict(out_dir: Path, align: pd.DataFrame, targets: pd.DataFrame) -> None:
    primary = align[
        align["spaceflight_vector"].eq("osd462_rna") & align["gene_subset"].eq("transport_target_genes")
    ]
    if primary.empty:
        primary = align[align["spaceflight_vector"].eq("osd462_rna") & align["gene_subset"].eq("dct1_core_genes")]
    hit = False
    label = "not_evaluable"
    primary_record = {}
    if not primary.empty:
        row = primary.iloc[0]
        target_dirs = targets["osd462_rna_vs_lowk_direction"].value_counts().to_dict() if "osd462_rna_vs_lowk_direction" in targets else {}
        opposite = target_dirs.get("opposite", 0)
        observed = sum(v for k, v in target_dirs.items() if k in {"opposite", "same"})
        opposite_fraction = float(opposite / observed) if observed else float("nan")
        hit = bool(row["cosine"] < -0.3 and row["cosine_ci_high"] < 0 and opposite_fraction >= 0.5)
        label = "strong_lowk_antialignment" if hit else "weak_or_bounded_lowk_alignment"
        primary_record = row.to_dict() | {"target_opposite_fraction": opposite_fraction}
    verdict = {
        "status": "complete",
        "classification": label,
        "promote_to_primary_result": hit,
        "primary_record": primary_record,
        "boundary": (
            "GSE228367 raw filtered matrices are used as DCT-enriched pseudobulk KD-NK response. "
            "DCT1/DCT2 specificity is gene-prior restricted, not cell-isolated low-K DCT1/DCT2 modeling."
        ),
    }
    (out_dir / "lowk_alignment_verdict.json").write_text(json.dumps(verdict, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=RUN_ROOT)
    ap.add_argument("--n-bootstrap", type=int, default=1000)
    args = ap.parse_args()

    out_dir = args.run_root / "perturbation"
    out_dir.mkdir(parents=True, exist_ok=True)
    lowk = build_lowk_effects(out_dir)
    vectors = load_spaceflight_vectors()
    gene_sets = load_gene_sets(args.run_root)
    align = alignment_rows(lowk, vectors, gene_sets, args.n_bootstrap)
    align.to_csv(out_dir / "lowk_dct_alignment_summary.tsv", sep="\t", index=False)
    target = target_gene_table(lowk, vectors, args.run_root)
    target.to_csv(out_dir / "lowk_target_gene_table.tsv", sep="\t", index=False)
    spec = align[
        align["gene_subset"].isin(
            [
                "dct1_core_genes",
                "dct2_core_genes",
                "dct_shared_genes",
                "dct1_top_decile_genes",
                "dct2_top_decile_genes",
                "dct1_top_quartile_genes",
                "dct2_top_quartile_genes",
            ]
        )
    ].copy()
    spec.to_csv(out_dir / "lowk_dct1_dct2_specificity.tsv", sep="\t", index=False)
    write_verdict(out_dir, align, target)
    print(f"[lowK] complete: {out_dir}")


if __name__ == "__main__":
    main()
