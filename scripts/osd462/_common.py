"""Shared paths and helpers for the OSD-462 multi-omics anchor layer scripts."""
from __future__ import annotations

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml

import sys
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.common import id_map_lookup  # noqa: E402

# Reference RRRM-2 run (fixed input; not re-derived)
RRRM2_RUN = "run_20260519_000547_2500g"
RRRM2_BASE = REPO_ROOT / "data" / "results" / RRRM2_RUN / "contrast_vectors" / "mechanism_axis"

# Input paths
OSD462_DIR = REPO_ROOT / "data" / "external" / "osdr" / "OSD-462"
PROTEOMICS_XLSX = OSD462_DIR / "GLDS-462_proteomics_2021-12-31_tc884-885_Protein_WorkUp.xlsx"
PHOSPHO_XLSX = OSD462_DIR / "GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx"
RNA_DE = OSD462_DIR / "GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv"
RNA_VST = OSD462_DIR / "GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv"
RNA_SAMPLES = OSD462_DIR / "GLDS-462_rna_seq_SampleTable_mRNA_GLbulkRNAseq.csv"

RRRM2_GENE_SCATTER = (RRRM2_BASE / "tubulointerstitial_state" / "lar_reversal"
                      / "lar_reversal_gene_scatter.tsv")
RRRM2_GENE_PRIORITY = RRRM2_BASE / "gene_axis_priority.tsv"
RRRM2_STATE_EFFECTS = (RRRM2_BASE / "tubulointerstitial_state"
                       / "state_space_group_effects.tsv")

ID_MAP = REPO_ROOT / "data" / "processed" / "resources" / "id_map.tsv"
MECHANISM_SETS = REPO_ROOT / "config" / "mechanism_gene_sets.yaml"

# OSDR RNA differential-expression contrast: Log2fc_(A)v(B) = A - B, so
# Space Flight - Ground Control is the flight effect.
RNA_FLIGHT_COL = "Log2fc_(Space Flight)v(Ground Control)"
RNA_FLIGHT_ADJP = "Adj.p.value_(Space Flight)v(Ground Control)"

SEED = 20260521
N_NULL = 10000
N_BOOTSTRAP = 2000


def anchor_dir(run: str) -> Path:
    d = REPO_ROOT / "data" / "results" / run / "osd462_anchor"
    d.mkdir(parents=True, exist_ok=True)
    return d


def default_run() -> str:
    return f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}_osd462_anchor"


def sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git_commit() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True
        ).strip()
    except Exception:
        return "unknown"


def write_manifest(out_dir: Path, analysis: str, inputs: dict, outputs: dict,
                   parameters: dict, name: str = "manifest.json") -> Path:
    manifest = {
        "analysis": analysis,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "parameters": parameters,
        "inputs": {},
        "outputs": {k: str(v) for k, v in outputs.items()},
    }
    for key, path in inputs.items():
        p = Path(path)
        manifest["inputs"][key] = str(p)
        if p.exists() and p.is_file():
            manifest["inputs"][f"{key}_sha256"] = sha256(p)
    path = out_dir / name
    path.write_text(json.dumps(manifest, indent=2))
    return path


def load_mechanism_sets(symbols: bool = False) -> dict[str, dict]:
    """Return the mechanism gene sets config (raw)."""
    return yaml.safe_load(MECHANISM_SETS.read_text())


def build_symbol_to_ensembl() -> dict[str, str]:
    """Symbol (lowercase) -> ENSMUSG, OSD-462 DE table preferred, id_map fallback.

    The OSD-462 differential-expression table provides an authoritative
    SYMBOL<->ENSEMBL bridge for this organism/build; we prefer it and fall back
    to the project id_map for symbols it does not contain.  One-to-many
    collisions keep the first occurrence (logged by the caller).
    """
    de = pd.read_csv(RNA_DE, usecols=["ENSEMBL", "SYMBOL"], dtype=str)
    de = de.dropna(subset=["ENSEMBL", "SYMBOL"])
    bridge: dict[str, str] = {}
    collisions = 0
    for sym, ens in zip(de["SYMBOL"], de["ENSEMBL"]):
        key = sym.strip().lower()
        if key in bridge:
            collisions += 1
            continue
        bridge[key] = ens.strip()
    # id_map fallback
    _, symbol_to_ens = id_map_lookup(ID_MAP)
    for key, ens_set in symbol_to_ens.items():
        if key not in bridge and len(ens_set) == 1:
            bridge[key] = next(iter(ens_set))
    bridge["_n_collisions"] = collisions  # type: ignore[assignment]
    return bridge


def resolve_symbols_to_ensembl(symbols, bridge: dict[str, str]) -> dict[str, str]:
    """Map an iterable of gene symbols to ENSMUSG using ``bridge``."""
    out = {}
    for s in symbols:
        key = str(s).strip().lower()
        if key in bridge and not key.startswith("_"):
            out[str(s)] = bridge[key]
    return out
