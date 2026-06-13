#!/usr/bin/env python3
"""Phenotype-anchoring layer -- animal-matched RNA state vs NCC activity."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.common import id_map_lookup  # noqa: E402
from src.multiomics.phenotype_anchor import (  # noqa: E402
    NCC_NONREGULATORY_SITES,
    NCC_REGULATORY_SITES,
    NCC_REGULATORY_SITES_SENS,
    channel_to_animal,
    compare_scores,
    per_sample_score,
    result_to_dict,
    rna_sample_animal,
)

PHOS_XLSX = "data/external/osdr/OSD-462/GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx"
RNA_VST = "data/external/osdr/OSD-462/GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv"
ID_MAP = "data/processed/resources/id_map.tsv"
# DCT/NCC-WNK transport gene set (project-standard panel)
DCT_GENES = ("Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Klhl3", "Cul3",
             "Calb1", "Kcnj10", "Kcnj16", "Sgk1", "Nedd4l")


def log(m: str) -> None:
    print(f"{datetime.now():%H:%M:%S} [phenotype] {m}")


def sha256(p: Path) -> str | None:
    if not p.exists():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def load_phospho_per_sample(xlsx: Path) -> pd.DataFrame:
    """Return a (gene|site) x channel matrix of log2 phosphosite intensities."""
    raw = pd.read_excel(xlsx, "siteQuant_360", header=None)
    labels = raw.iloc[1].tolist()  # row 1 holds per-channel labels
    df = pd.read_excel(xlsx, "siteQuant_360", header=2)
    df.columns = [str(c) for c in df.columns]
    sn_cols = [c for c in df.columns if "sn_sum" in c.lower()]
    # align each sn_sum column to its channel label by raw-sheet position
    pos = {c: list(df.columns).index(c) for c in sn_cols}
    chan = {c: str(labels[pos[c]]).strip() for c in sn_cols}
    mat = df[sn_cols].apply(pd.to_numeric, errors="coerce")
    mat.index = [f"{g}|{s}" for g, s in zip(df["gene_symbol"].astype(str).str.strip(),
                                            df["Site Position"].astype(str).str.strip())]
    mat.columns = [chan[c] for c in sn_cols]
    mat = mat[mat > 0]
    return np.log2(mat)


def load_rna_dct(vst_path: Path, id_map_path: Path, *, drop_slc12a3: bool = False) -> pd.DataFrame:
    """Return a DCT/NCC-WNK gene x RNA-sample VST matrix (ENSMUSG-indexed)."""
    vst = pd.read_csv(vst_path, index_col=0)
    vst.index = vst.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    _ens_to_sym, sym_to_ens = id_map_lookup(str(id_map_path))
    genes = [g for g in DCT_GENES if not (drop_slc12a3 and g == "Slc12a3")]
    ids = sorted({e for g in genes for e in sym_to_ens.get(g.lower(), set())} & set(vst.index))
    sub = vst.loc[ids].copy()
    sub.index = [f"{i}" for i in sub.index]
    return sub


def collapse_rna_by_animal(score: pd.Series) -> tuple[pd.Series, pd.Series]:
    """Collapse techrep RNA columns to one value per animal; return (score, condition)."""
    rows = {}
    for sample, val in score.items():
        m = rna_sample_animal(sample)
        if m is None:
            continue
        cond, animal = m
        rows.setdefault((cond, animal), []).append(val)
    idx = sorted(rows)
    s = pd.Series({f"{c}|{a}": float(np.mean(v)) for (c, a), v in rows.items()})
    cond = pd.Series({f"{c}|{a}": c for (c, a) in rows})
    return s.reindex([f"{c}|{a}" for c, a in idx]), cond.reindex([f"{c}|{a}" for c, a in idx])


def phospho_score_by_animal(phospho: pd.DataFrame, sites: tuple) -> tuple[pd.Series, pd.Series]:
    """Per-animal phospho score and condition, keyed identically to the RNA score."""
    keys = [f"{g}|{s}" for g, s in sites]
    per_channel = per_sample_score(phospho, keys)
    rows, cond = {}, {}
    for channel, val in per_channel.items():
        m = channel_to_animal(channel)
        if m is None:
            continue
        c, a = m
        rows.setdefault(f"{c}|{a}", []).append(val)
        cond[f"{c}|{a}"] = c
    s = pd.Series({k: float(np.mean(v)) for k, v in rows.items()})
    return s, pd.Series(cond)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phospho", default=PHOS_XLSX)
    ap.add_argument("--rna-vst", default=RNA_VST)
    ap.add_argument("--id-map", default=ID_MAP)
    ap.add_argument("--outdir", default="")
    args = ap.parse_args()

    phos_path = REPO_ROOT / args.phospho
    vst_path = REPO_ROOT / args.rna_vst
    idmap_path = REPO_ROOT / args.id_map
    outdir = (REPO_ROOT / args.outdir if args.outdir else
              REPO_ROOT / f"data/results/run_{datetime.now():%Y%m%d}_phenotype_anchor")
    outdir.mkdir(parents=True, exist_ok=True)

    log("Loading OSD-462 phosphoproteomics (per-sample siteQuant)")
    phospho = load_phospho_per_sample(phos_path)
    log(f"  phosphosites x channels: {phospho.shape}")

    log("Building per-animal NCC activity scores")
    reg_score, reg_cond = phospho_score_by_animal(phospho, NCC_REGULATORY_SITES)
    sens_score, _ = phospho_score_by_animal(phospho, NCC_REGULATORY_SITES_SENS)
    ctrl_score, _ = phospho_score_by_animal(phospho, NCC_NONREGULATORY_SITES)

    log("Building per-animal RNA DCT transport state scores")
    rna_full = load_rna_dct(vst_path, idmap_path, drop_slc12a3=False)
    rna_noslc = load_rna_dct(vst_path, idmap_path, drop_slc12a3=True)
    # RNA score = mean-z of DCT/NCC-WNK transport genes (higher = more DCT
    rna_score_full, rna_cond = collapse_rna_by_animal(per_sample_score(rna_full, list(rna_full.index)))
    rna_score_noslc, _ = collapse_rna_by_animal(per_sample_score(rna_noslc, list(rna_noslc.index)))

    common = sorted(set(reg_score.index) & set(rna_score_full.index))
    log(f"  animals matched (RNA + phospho): {len(common)}")
    cond = reg_cond.reindex(common)

    comparisons = {
        "regulatory_vs_rna": compare_scores(rna_score_full.reindex(common),
                                            reg_score.reindex(common), cond),
        "regulatory_sens_T89_vs_rna": compare_scores(rna_score_full.reindex(common),
                                                     sens_score.reindex(common), cond),
        "control_nonregulatory_vs_rna": compare_scores(rna_score_full.reindex(common),
                                                       ctrl_score.reindex(common), cond),
        "regulatory_vs_rna_no_Slc12a3": compare_scores(rna_score_noslc.reindex(common),
                                                       reg_score.reindex(common), cond),
    }

    rows = []
    for name, res in comparisons.items():
        d = result_to_dict(res)
        d["comparison"] = name
        rows.append(d)
    summary = pd.DataFrame(rows)[
        ["comparison", "n_animals", "group_rna_flt_minus_gc", "group_phospho_flt_minus_gc",
         "spearman_all", "spearman_all_p", "spearman_condition_adjusted",
         "spearman_condition_adjusted_p", "interpretation"]]
    summary.to_csv(outdir / "phenotype_anchor_summary.tsv", sep="\t", index=False)

    per_animal = pd.DataFrame({
        "condition": cond,
        "rna_dct_transport_score": rna_score_full.reindex(common),
        "ncc_activity_score_regulatory": reg_score.reindex(common),
        "ncc_activity_score_nonregulatory_control": ctrl_score.reindex(common),
    }).reset_index(names="animal")
    per_animal.to_csv(outdir / "phenotype_anchor_per_animal.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "phenotype-anchoring: animal-matched RNA state vs NCC activity",
        "framing": "molecular activity-phenotype anchoring (Tier 1); not tissue/physiology",
        "timestamp": datetime.now().isoformat(),
        "inputs": {
            "phospho_xlsx": str(phos_path), "phospho_sha256": sha256(phos_path),
            "rna_vst": str(vst_path), "rna_vst_sha256": sha256(vst_path),
            "id_map": str(idmap_path), "id_map_sha256": sha256(idmap_path),
        },
        "regulatory_sites": [f"{g}:{s}" for g, s in NCC_REGULATORY_SITES],
        "nonregulatory_control_sites": [f"{g}:{s}" for g, s in NCC_NONREGULATORY_SITES],
        "dct_rna_gene_set": list(DCT_GENES),
        "n_animals_matched": len(common),
        "comparisons": {k: result_to_dict(v) for k, v in comparisons.items()},
        "outputs": {
            "summary": str(outdir / "phenotype_anchor_summary.tsv"),
            "per_animal": str(outdir / "phenotype_anchor_per_animal.tsv"),
        },
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"Wrote outputs to {outdir}")
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
