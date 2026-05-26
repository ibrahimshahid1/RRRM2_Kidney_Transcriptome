#!/usr/bin/env python3
"""External kidney injury/repair spatial reference projection for v11.

Primary spatial source:
  * GSE269622 Visium, whole-transcriptome Space Ranger archives.

Secondary spatial source:
  * GSE269719 processed Xenium AnnData, used only for annotation/neighborhood
    inventory because the panel is targeted and not a genome-wide projection
    source.

This script does not localize the RR-10 spaceflight lesion. It asks which IRI
timepoints or marker-enriched Visium spatial niches most resemble the bulk
spaceflight RNA remodeling vector.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.spatial import cKDTree


REPO_ROOT = Path(__file__).resolve().parents[2]
RUN_ROOT = REPO_ROOT / "data/results/run_20260526_v11_dct1_phospho_mediation"
VISIUM_DIR = REPO_ROOT / "data/external/spatial_reference/GSE269622_Visium"
XENIUM_H5AD = REPO_ROOT / "data/external/spatial_reference/GSE269719_Xenium/Xenium.h5ad"
OUT_DIR = RUN_ROOT / "spatial_reference"
MECHANISM_SETS = REPO_ROOT / "config/mechanism_gene_sets.yaml"
OSD462_EFFECTS = REPO_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_flight_effects.tsv"
OSD462_PATHWAY = REPO_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor/osd462_rna_pathway_effects.tsv"
RRRM2_GENE = (
    REPO_ROOT
    / "data/results/run_20260519_000547_2500g/contrast_vectors/mechanism_axis/tubulointerstitial_state/lar_reversal/lar_reversal_gene_scatter.tsv"
)


TIME_ORDER = {
    "sham": 0,
    "hour4": 1,
    "hour12": 2,
    "day2": 3,
    "day14": 4,
    "week6": 5,
}

NICHE_MARKERS = {
    "injured_tubule_enriched": ["Havcr1", "Vcam1", "Krt8", "Krt18", "Lcn2", "Sox9"],
    "fibroblast_interstitial_enriched": ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdgfra", "Fn1", "Sparc"],
    "endothelial_enriched": ["Pecam1", "Cdh5", "Kdr", "Flt1", "Emcn", "Klf2"],
    "immune_enriched": ["Ptprc", "Lyz2", "Adgre1", "C1qa", "Itgam", "Cd68"],
    "dct_enriched": ["Slc12a3", "Pvalb", "Trpm6", "Wnk4", "Stk39", "Klhl3"],
    "fibro_inflammatory_repair_enriched": [
        "Vcam1",
        "Ccl2",
        "Il34",
        "Serpine1",
        "Cd44",
        "Klf5",
        "Col1a1",
        "Lyz2",
    ],
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cosine(a: pd.Series, b: pd.Series) -> tuple[float, int]:
    joined = pd.concat([a, b], axis=1, join="inner").dropna()
    if joined.empty:
        return np.nan, 0
    x = joined.iloc[:, 0].to_numpy(dtype=float)
    y = joined.iloc[:, 1].to_numpy(dtype=float)
    denom = np.linalg.norm(x) * np.linalg.norm(y)
    if denom <= 0:
        return np.nan, len(joined)
    return float(np.dot(x, y) / denom), int(len(joined))


def parse_timepoint(path: Path) -> str:
    name = path.name
    for key in TIME_ORDER:
        if f"_{key}_" in name:
            return key
    raise ValueError(f"Cannot infer timepoint from {path}")


def extract_member(archive: Path, member_suffix: str, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    out = dest / Path(member_suffix).name
    if out.exists() and out.stat().st_size > 0:
        return out
    with tarfile.open(archive, "r:gz") as tf:
        matches = [m for m in tf.getmembers() if m.name.endswith(member_suffix)]
        if not matches:
            raise FileNotFoundError(f"{member_suffix} not found in {archive}")
        member = matches[0]
        src = tf.extractfile(member)
        if src is None:
            raise FileNotFoundError(member.name)
        with out.open("wb") as fh:
            while True:
                chunk = src.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
    return out


def load_gene_sets() -> dict[str, list[str]]:
    raw = yaml.safe_load(MECHANISM_SETS.read_text())
    return {name: [str(g).title() for g in cfg.get("genes", [])] for name, cfg in raw.items()}


def load_spaceflight_vectors() -> tuple[dict[str, pd.Series], dict[str, pd.Series]]:
    osd = pd.read_csv(OSD462_EFFECTS, sep="\t")
    osd_vec = (
        osd.dropna(subset=["gene_symbol", "osd462_rna_effect"])
        .drop_duplicates("gene_symbol")
        .set_index("gene_symbol")["osd462_rna_effect"]
    )
    rr = pd.read_csv(RRRM2_GENE, sep="\t")
    rr_vec = (
        rr.dropna(subset=["mgi_symbol", "iss_t_effect"])
        .drop_duplicates("mgi_symbol")
        .set_index("mgi_symbol")["iss_t_effect"]
    )
    pathway = pd.read_csv(OSD462_PATHWAY, sep="\t").set_index("pathway")
    pathway_vectors = {
        "osd462_rna": pathway["osd462_rna_pathway_effect"],
        "rrrm2_iss_t": pathway["rrrm2_iss_t_pathway_effect"],
    }
    return {"osd462_rna": osd_vec, "rrrm2_iss_t": rr_vec}, pathway_vectors


def read_positions(path: Path) -> pd.DataFrame:
    pos = pd.read_csv(path)
    if "barcode" not in pos.columns:
        pos = pd.read_csv(path, header=None)
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]
    return pos


def expression_mean(adata) -> pd.Series:
    x = adata.X
    arr = np.asarray(x.mean(axis=0)).ravel()
    return pd.Series(arr, index=adata.var_names.astype(str))


def score_genes(adata, genes: list[str]) -> pd.Series:
    gene_to_idx = {str(g).lower(): i for i, g in enumerate(adata.var_names)}
    idx = [gene_to_idx[g.lower()] for g in genes if g.lower() in gene_to_idx]
    if not idx:
        return pd.Series(np.nan, index=adata.obs_names)
    x = adata.X[:, idx]
    return pd.Series(np.asarray(x.mean(axis=1)).ravel(), index=adata.obs_names)


def dct_adjacent_mask(adata, dct_score: pd.Series, positions: pd.DataFrame) -> pd.Series:
    pos = positions.set_index("barcode")
    common = adata.obs_names.intersection(pos.index)
    mask = pd.Series(False, index=adata.obs_names)
    if len(common) < 10 or dct_score.dropna().empty:
        return mask
    coords = pos.loc[common, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)
    dct_common = dct_score.loc[common]
    dct_high = dct_common >= dct_common.quantile(0.75)
    if dct_high.sum() == 0:
        return mask
    tree_all = cKDTree(coords)
    dists, _ = tree_all.query(coords, k=2)
    radius = float(np.nanmedian(dists[:, 1]) * 1.6)
    tree_dct = cKDTree(coords[dct_high.to_numpy()])
    nearest, _ = tree_dct.query(coords, k=1)
    adjacent = (nearest <= radius) & (~dct_high.to_numpy())
    mask.loc[common] = adjacent
    return mask


def load_visium_sample(archive: Path, extracted_root: Path):
    import scanpy as sc

    timepoint = parse_timepoint(archive)
    sample_dir = extracted_root / archive.stem.replace(".tar", "")
    h5 = extract_member(archive, "outs/filtered_feature_bc_matrix.h5", sample_dir)
    pos_path = extract_member(archive, "outs/spatial/tissue_positions.csv", sample_dir)
    adata = sc.read_10x_h5(h5)
    adata.var_names_make_unique()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    positions = read_positions(pos_path)
    return timepoint, adata, positions


def pathway_effects(vec: pd.Series, gene_sets: dict[str, list[str]]) -> pd.Series:
    out = {}
    lower = {str(g).lower(): g for g in vec.index}
    for name, genes in gene_sets.items():
        present = [lower[g.lower()] for g in genes if g.lower() in lower]
        out[name] = float(vec.loc[present].mean()) if present else np.nan
    return pd.Series(out)


def run_visium_projection() -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    extracted_root = OUT_DIR / "visium_extracted_minimal"
    archives = sorted(VISIUM_DIR.glob("GSM*_visium_*_spaceranger.tar.gz"))
    if not archives:
        raise FileNotFoundError(f"No Visium archives found in {VISIUM_DIR}")

    gene_sets = load_gene_sets()
    sf_gene_vectors, sf_pathway_vectors = load_spaceflight_vectors()
    sample_means = {}
    spot_score_rows = []
    niche_means = {}
    inventory_rows = []

    for archive in archives:
        timepoint, adata, positions = load_visium_sample(archive, extracted_root)
        sample_means[timepoint] = expression_mean(adata)
        score_frame = pd.DataFrame(index=adata.obs_names)
        for niche, genes in NICHE_MARKERS.items():
            score_frame[niche] = score_genes(adata, genes)
        score_frame["dct_adjacent_spots"] = dct_adjacent_mask(adata, score_frame["dct_enriched"], positions)

        for niche in NICHE_MARKERS:
            score = score_frame[niche]
            keep = score >= score.quantile(0.75)
            if keep.sum() >= 10:
                niche_means[(timepoint, niche)] = expression_mean(adata[keep.to_numpy(), :])
            spot_score_rows.append(
                {
                    "timepoint": timepoint,
                    "niche": niche,
                    "n_selected_spots": int(keep.sum()),
                    "n_total_spots": int(adata.n_obs),
                    "selection_rule": "top_quartile_marker_score",
                }
            )
        adj = score_frame["dct_adjacent_spots"].astype(bool)
        if adj.sum() >= 10:
            niche_means[(timepoint, "dct_adjacent_spots")] = expression_mean(adata[adj.to_numpy(), :])
        spot_score_rows.append(
            {
                "timepoint": timepoint,
                "niche": "dct_adjacent_spots",
                "n_selected_spots": int(adj.sum()),
                "n_total_spots": int(adata.n_obs),
                "selection_rule": "within_1.6_median_spot_spacing_of_dct_high_spots_not_dct_high",
            }
        )
        inventory_rows.append(
            {
                "archive": archive.name,
                "timepoint": timepoint,
                "n_spots": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "sha256": sha256(archive),
            }
        )

    sample_df = pd.DataFrame(sample_means).T.sort_index(key=lambda x: x.map(TIME_ORDER))
    grand = sample_df.mean(axis=0)
    sham = sample_df.loc["sham"]

    gene_rows = []
    pathway_rows = []
    for timepoint in sample_df.index:
        state_vec = sample_df.loc[timepoint] - grand
        effect_vec = sample_df.loc[timepoint] - sham
        for sf_name, sf_vec in sf_gene_vectors.items():
            c_state, n_state = cosine(sf_vec, state_vec)
            c_effect, n_effect = cosine(sf_vec, effect_vec)
            gene_rows.append(
                {
                    "timepoint": timepoint,
                    "spaceflight_vector": sf_name,
                    "comparison": "visium_timepoint_state_minus_grand_mean",
                    "cosine": c_state,
                    "n_genes": n_state,
                }
            )
            gene_rows.append(
                {
                    "timepoint": timepoint,
                    "spaceflight_vector": sf_name,
                    "comparison": "visium_timepoint_minus_sham",
                    "cosine": c_effect,
                    "n_genes": n_effect,
                }
            )
        spatial_pathway_state = pathway_effects(state_vec, gene_sets)
        spatial_pathway_effect = pathway_effects(effect_vec, gene_sets)
        for sf_name, sf_path in sf_pathway_vectors.items():
            c_state, n_state = cosine(sf_path, spatial_pathway_state)
            c_effect, n_effect = cosine(sf_path, spatial_pathway_effect)
            pathway_rows.append(
                {
                    "timepoint": timepoint,
                    "spaceflight_pathway_vector": sf_name,
                    "comparison": "visium_timepoint_state_minus_grand_mean",
                    "cosine": c_state,
                    "n_pathways": n_state,
                }
            )
            pathway_rows.append(
                {
                    "timepoint": timepoint,
                    "spaceflight_pathway_vector": sf_name,
                    "comparison": "visium_timepoint_minus_sham",
                    "cosine": c_effect,
                    "n_pathways": n_effect,
                }
            )

    niche_rows = []
    for (timepoint, niche), vec in niche_means.items():
        if ("sham", niche) not in niche_means:
            continue
        effect_vec = vec - niche_means[("sham", niche)]
        spatial_path = pathway_effects(effect_vec, gene_sets)
        for sf_name, sf_vec in sf_gene_vectors.items():
            c, n = cosine(sf_vec, effect_vec)
            niche_rows.append(
                {
                    "timepoint": timepoint,
                    "niche": niche,
                    "spaceflight_vector": sf_name,
                    "comparison": "visium_niche_mean_minus_sham_same_niche",
                    "cosine": c,
                    "n_genes": n,
                }
            )
        for sf_name, sf_path in sf_pathway_vectors.items():
            c, n = cosine(sf_path, spatial_path)
            niche_rows.append(
                {
                    "timepoint": timepoint,
                    "niche": niche,
                    "spaceflight_vector": sf_name,
                    "comparison": "visium_niche_pathway_minus_sham_same_niche",
                    "cosine": c,
                    "n_pathways": n,
                }
            )

    pd.DataFrame(inventory_rows).to_csv(OUT_DIR / "spatial_reference_dataset_inventory.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(OUT_DIR / "visium_timepoint_gene_cosines.tsv", sep="\t", index=False)
    pd.DataFrame(pathway_rows).to_csv(OUT_DIR / "visium_timepoint_pathway_cosines.tsv", sep="\t", index=False)
    pd.DataFrame(niche_rows).to_csv(OUT_DIR / "visium_niche_cosines.tsv", sep="\t", index=False)
    pd.DataFrame(spot_score_rows).to_csv(OUT_DIR / "visium_spot_niche_counts.tsv", sep="\t", index=False)

    best = pd.DataFrame(pathway_rows)
    best = best[best["comparison"] == "visium_timepoint_minus_sham"].copy()
    best = best[best["timepoint"] != "sham"].sort_values("cosine", ascending=False)
    return {
        "n_visium_archives": len(archives),
        "best_timepoint_pathway_match": best.head(1).to_dict(orient="records"),
        "claim_label": "external injury/repair spatial contextualization, not spaceflight spatial validation",
    }


def run_xenium_inventory() -> dict:
    if not XENIUM_H5AD.exists():
        return {"status": "not_downloaded", "path": str(XENIUM_H5AD)}
    h5 = ad.read_h5ad(XENIUM_H5AD, backed="r")
    obs = h5.obs
    rows = []
    candidate_cols = [
        c
        for c in obs.columns
        if any(k in c.lower() for k in ["cell", "type", "annot", "time", "sample", "condition", "cn", "neighbor"])
    ]
    for col in candidate_cols[:30]:
        vc = obs[col].astype(str).value_counts(dropna=False).head(50)
        for value, count in vc.items():
            rows.append({"column": col, "value": value, "count": int(count)})
    pd.DataFrame(rows).to_csv(OUT_DIR / "xenium_annotation_inventory.tsv", sep="\t", index=False)
    summary = {
        "status": "loaded_backed",
        "path": str(XENIUM_H5AD),
        "sha256": sha256(XENIUM_H5AD),
        "n_cells": int(h5.n_obs),
        "n_genes": int(h5.n_vars),
        "obs_columns": list(obs.columns),
        "obsm_keys": list(h5.obsm.keys()),
        "claim_label": "secondary targeted-panel annotation/neighborhood context only",
    }
    h5.file.close()
    return summary


def main() -> None:
    global RUN_ROOT, OUT_DIR

    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=RUN_ROOT)
    ap.add_argument("--skip-visium", action="store_true")
    ap.add_argument("--skip-xenium", action="store_true")
    args = ap.parse_args()
    RUN_ROOT = args.run_root
    OUT_DIR = RUN_ROOT / "spatial_reference"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.skip_visium and (OUT_DIR / "spatial_reference_projection_verdict.json").exists():
        previous = json.loads((OUT_DIR / "spatial_reference_projection_verdict.json").read_text())
        visium = previous.get("visium", {"status": "skipped"})
    elif args.skip_visium:
        visium = {"status": "skipped"}
    else:
        visium = run_visium_projection()
    xenium = {"status": "skipped"} if args.skip_xenium else run_xenium_inventory()
    verdict = {
        "status": "complete",
        "visium": visium,
        "xenium": xenium,
        "interpretation_boundary": (
            "GSE269622/GSE269719 are external IRI injury/repair references. "
            "They spatially contextualize the bulk spaceflight RNA vector but do not validate or localize "
            "spaceflight lesions in RR-10 tissue."
        ),
    }
    (OUT_DIR / "spatial_reference_projection_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(f"[spatial] complete: {OUT_DIR}")


if __name__ == "__main__":
    main()
