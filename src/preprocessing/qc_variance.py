"""Pre-residualization sample QC and approximate variance partitioning."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression

from src.common import REPO_ROOT, find_sample_col, normalize_labels


def _clean_id(x: object) -> str:
    value = str(x).strip()
    value = Path(value).name
    for suffix in [".bam", ".sam", ".fastq", ".fq", ".gz", ".txt", ".csv"]:
        if value.lower().endswith(suffix):
            value = value[: -len(suffix)]
    return value.replace(".", "_").replace("-", "_").replace("/", "_")


def load_counts(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    counts = pd.read_csv(path, sep=sep)
    gene_col = counts.columns[0]
    counts = counts.set_index(gene_col)
    counts.index = counts.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    counts = counts.groupby(level=0).sum()
    return counts.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def align_counts_meta(counts: pd.DataFrame, meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample_col = find_sample_col(meta)
    meta = normalize_labels(meta.copy())
    meta["_clean_sample"] = meta[sample_col].map(_clean_id)
    count_clean = {_clean_id(c): c for c in counts.columns}
    common_clean = [s for s in meta["_clean_sample"] if s in count_clean]
    if not common_clean:
        raise ValueError("Could not align counts to metadata for QC")
    count_cols = [count_clean[s] for s in common_clean]
    meta = meta.set_index("_clean_sample").loc[common_clean].copy()
    counts = counts[count_cols].copy()
    counts.columns = common_clean
    return counts, meta


def log_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.divide(lib, axis=1) * 1_000_000.0
    return np.log2(cpm.fillna(0.0) + 1.0)


def sample_qc_table(counts: pd.DataFrame, expr: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    detected = (counts > 0).sum(axis=0)
    med_expr = expr.median(axis=0)
    qc = pd.DataFrame({
        "sample": counts.columns,
        "library_size": lib.values.astype(float),
        "detected_genes": detected.values.astype(int),
        "median_log2cpm": med_expr.values.astype(float),
    })
    numeric = qc[["library_size", "detected_genes", "median_log2cpm"]]
    z = (numeric - numeric.mean()) / (numeric.std(ddof=0).replace(0, np.nan))
    qc["max_abs_qc_z"] = z.abs().max(axis=1).fillna(0.0).values
    qc["qc_outlier_flag"] = qc["max_abs_qc_z"] >= 3.0
    return qc


def pca_outliers(expr: pd.DataFrame, n_components: int = 5) -> tuple[pd.DataFrame, np.ndarray]:
    top_genes = expr.var(axis=1).sort_values(ascending=False).head(min(2000, expr.shape[0])).index
    X = expr.loc[top_genes].T.values.astype(float)
    n_components = min(n_components, X.shape[0] - 1, X.shape[1])
    if n_components < 2:
        return pd.DataFrame({"sample": expr.columns, "pca_outlier_flag": False}), np.array([])
    pca = PCA(n_components=n_components, random_state=0)
    scores = pca.fit_transform(X)
    dist = np.sqrt(((scores - scores.mean(axis=0)) ** 2).sum(axis=1))
    z = (dist - dist.mean()) / (dist.std(ddof=0) + 1e-8)
    out = pd.DataFrame({
        "sample": expr.columns,
        "pca_distance": dist,
        "pca_distance_z": z,
        "pca_outlier_flag": np.abs(z) >= 3.0,
    })
    for i in range(n_components):
        out[f"PC{i + 1}"] = scores[:, i]
    return out, pca.explained_variance_ratio_


def variance_partition(expr: pd.DataFrame, meta: pd.DataFrame, covariates: list[str], max_genes: int = 1000) -> pd.DataFrame:
    genes = expr.var(axis=1).sort_values(ascending=False).head(min(max_genes, expr.shape[0])).index
    y = expr.loc[genes].T.values.astype(float)
    rows = []
    for cov in covariates:
        if cov not in meta.columns:
            continue
        series = meta[cov]
        if pd.api.types.is_numeric_dtype(series):
            x = pd.to_numeric(series, errors="coerce").fillna(series.median()).values.reshape(-1, 1)
        else:
            dummies = pd.get_dummies(series.astype(str), drop_first=False)
            if dummies.shape[1] < 2:
                continue
            x = dummies.values.astype(float)
        model = LinearRegression()
        model.fit(x, y)
        pred = model.predict(x)
        ss_res = ((y - pred) ** 2).sum(axis=0)
        ss_tot = ((y - y.mean(axis=0)) ** 2).sum(axis=0)
        r2 = np.divide(ss_res, ss_tot, out=np.ones_like(ss_tot), where=ss_tot > 0)
        r2 = 1.0 - r2
        rows.append({
            "covariate": cov,
            "n_genes": int(len(genes)),
            "median_gene_r2": float(np.nanmedian(r2)),
            "mean_gene_r2": float(np.nanmean(r2)),
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Pre-residualization QC and approximate variance partitioning")
    ap.add_argument("--counts", default=str(REPO_ROOT / "data/processed/aligned_outputs/rsem_rRNArm_raw_counts.csv"))
    ap.add_argument("--meta", default=str(REPO_ROOT / "data/processed/aligned_outputs/metadata_aligned.tsv"))
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/pre_residualization_qc"))
    ap.add_argument("--covariates", default="Age,Arm,EnvGroup,LibraryBatch,SeqInstr,ReadDepth,rRNA,RIN")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    counts = load_counts(Path(args.counts))
    meta = pd.read_csv(args.meta, sep="\t")
    counts, meta = align_counts_meta(counts, meta)
    expr = log_cpm(counts)

    sample_qc = sample_qc_table(counts, expr)
    pca_qc, explained = pca_outliers(expr)
    covariates = [c.strip() for c in args.covariates.split(",") if c.strip()]
    var_part = variance_partition(expr, meta, covariates)

    sample_qc.merge(pca_qc, on="sample", how="left").to_csv(
        outdir / "sample_qc_outliers.tsv", sep="\t", index=False
    )
    var_part.to_csv(outdir / "variance_partition_approx.tsv", sep="\t", index=False)
    metadata = {
        "counts": str(args.counts),
        "meta": str(args.meta),
        "n_samples": int(counts.shape[1]),
        "n_genes": int(counts.shape[0]),
        "pca_explained_variance_ratio": explained.tolist(),
        "note": "QC only; no samples are removed automatically.",
    }
    (outdir / "qc_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"[OK] Pre-residualization QC outputs written to {outdir}")


if __name__ == "__main__":
    main()
