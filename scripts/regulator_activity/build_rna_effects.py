#!/usr/bin/env python3
"""Build per-cohort gene-level flight-effect tables for regulator Layer B.

decoupler PROGENy / CollecTRI priors are keyed by *gene symbol*, so each cohort
table is emitted as ``gene`` (mgi symbol) + ``stat`` (a signed flight-effect
statistic, Space Flight - Ground Control). The per-cohort effect definition is:

* RRRM-2 ISS-T / LAR (Young): Wald-style z = log2FC / lfcSE from the project
  limma gene-level DE tables (``data/processed/gene_level_DE``).
* OSD-462 RNA: the OSDR-provided moderated ``Stat_(Space Flight)v(Ground
  Control)`` (symbol bridge from the same DE table).
* OSD-513 / OSD-253: Welch t per gene (Space Flight vs Ground Control) computed
  directly from the GeneLab VST matrix (no DE table downloaded). OSD-253 pools
  strains/durations into a single cohort-level SF-vs-GC contrast and is treated
  as a context cohort.

All effects are *within cohort*; decoupler Layer B is run per cohort and only
the cross-cohort recurrence of activity sign is interpreted (never raw pooling).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.common import id_map_lookup  # noqa: E402

OSD = REPO_ROOT / "data" / "external" / "osdr"
GENE_DE = REPO_ROOT / "data" / "processed" / "gene_level_DE"
ID_MAP = REPO_ROOT / "data" / "processed" / "resources" / "id_map.tsv"
EPS = 1e-12


def _strip_version(idx: pd.Index) -> pd.Index:
    return idx.astype(str).str.replace(r"\.\d+$", "", regex=True)


def _collapse_symbol(df: pd.DataFrame) -> pd.DataFrame:
    """Mean-collapse duplicate gene symbols; drop missing/blank symbols."""
    df = df[df["gene"].notna() & (df["gene"].astype(str).str.strip() != "")]
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["stat"])
    return df.groupby("gene", as_index=False)["stat"].mean()


def from_limma_de(path: Path, ens_to_sym: dict) -> pd.DataFrame:
    de = pd.read_csv(path, sep="\t")
    de["gene"] = _strip_version(de["gene"]).map(ens_to_sym)
    z = de["log2FC"].astype(float) / de["lfcSE"].astype(float).replace(0, np.nan)
    out = pd.DataFrame({"gene": de["gene"], "stat": z})
    return _collapse_symbol(out)


def from_osd462_de() -> pd.DataFrame:
    path = OSD / "OSD-462" / "GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv"
    de = pd.read_csv(path, usecols=["SYMBOL", "Stat_(Space Flight)v(Ground Control)"])
    out = pd.DataFrame({"gene": de["SYMBOL"],
                        "stat": de["Stat_(Space Flight)v(Ground Control)"].astype(float)})
    return _collapse_symbol(out)


def _runsheet_groups(runsheet: Path) -> pd.Series:
    rs = pd.read_csv(runsheet)
    name_col = "Sample Name" if "Sample Name" in rs.columns else rs.columns[0]
    sf_col = next(c for c in rs.columns if "Spaceflight" in c or "Space Flight" in c)
    grp = {}
    for name, val in zip(rs[name_col], rs[sf_col]):
        v = str(val).lower()
        if "flight" in v:
            grp[str(name)] = "FLT"
        elif "ground" in v:
            grp[str(name)] = "GC"
    return pd.Series(grp)


def welch_t_from_vst(vst_path: Path, groups: pd.Series, ens_to_sym: dict) -> pd.DataFrame:
    vst = pd.read_csv(vst_path, index_col=0)
    vst.index = _strip_version(vst.index)
    fl = [c for c in vst.columns if groups.get(c) == "FLT"]
    gc = [c for c in vst.columns if groups.get(c) == "GC"]
    if not fl or not gc:
        # fall back to substring matching on column names
        fl = [c for c in vst.columns if "FLT" in c or "_FL" in c]
        gc = [c for c in vst.columns if "_GC" in c or "Ground" in c]
    print(f"    VST {vst_path.name}: n_FLT={len(fl)} n_GC={len(gc)} genes={vst.shape[0]}")
    a = vst[fl].to_numpy(float)
    b = vst[gc].to_numpy(float)
    ma, mb = a.mean(1), b.mean(1)
    va, vb = a.var(1, ddof=1), b.var(1, ddof=1)
    se = np.sqrt(va / a.shape[1] + vb / b.shape[1])
    t = np.where(se > EPS, (ma - mb) / se, np.nan)
    out = pd.DataFrame({"gene": pd.Index(vst.index).map(ens_to_sym), "stat": t})
    return _collapse_symbol(out)


def main() -> int:
    ens_to_sym, _ = id_map_lookup(str(ID_MAP))
    outdir = REPO_ROOT / "data" / "results" / "run_20260522_regulator_activity" / "rna_effects"
    outdir.mkdir(parents=True, exist_ok=True)

    cohorts: dict[str, pd.DataFrame] = {}
    print("[build] RRRM-2 ISS-T Young (limma z)")
    cohorts["RRRM2_ISS_T_young"] = from_limma_de(GENE_DE / "ISS_T_YNG_FLT_vs_GC_gene_DE.tsv", ens_to_sym)
    print("[build] RRRM-2 LAR Young (limma z)")
    cohorts["RRRM2_LAR_young"] = from_limma_de(GENE_DE / "LAR_YNG_FLT_vs_GC_gene_DE.tsv", ens_to_sym)
    print("[build] OSD-462 RNA (OSDR moderated Stat)")
    cohorts["OSD462_rna"] = from_osd462_de()
    print("[build] OSD-513 (Welch t from VST)")
    cohorts["OSD513"] = welch_t_from_vst(
        OSD / "OSD-513" / "GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        _runsheet_groups(OSD / "OSD-513" / "GLDS-513_rna_seq_bulkRNASeq_v2_runsheet.csv"),
        ens_to_sym)
    print("[build] OSD-253 (Welch t from VST; strains pooled, context cohort)")
    cohorts["OSD253"] = welch_t_from_vst(
        OSD / "OSD-253" / "GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        _runsheet_groups(OSD / "OSD-253" / "GLDS-253_rna_seq_bulkRNASeq_v2_runsheet.csv"),
        ens_to_sym)

    spec = {}
    for label, df in cohorts.items():
        path = outdir / f"{label}_gene_stat.tsv"
        df.to_csv(path, sep="\t", index=False)
        spec[label] = str(path.relative_to(REPO_ROOT))
        print(f"  wrote {label}: {len(df)} genes -> {path.name}")

    (outdir / "rna_effects_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    print(f"[build] spec -> {outdir / 'rna_effects_spec.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
