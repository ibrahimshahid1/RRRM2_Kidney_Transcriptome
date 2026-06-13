#!/usr/bin/env python3
"""v10 regulator-activity prioritization -- orchestrator."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.multiomics.regulator_activity import (  # noqa: E402
    grade_axis,
    ksea,
    ksea_positive_control_passes,
    load_kinase_substrate_net,
    recurrence_class,
    run_ulm_activity,
)

DEF_PHOSPHO = "data/results/run_20260521_osd462_anchor/osd462_anchor/phospho_all_sites.tsv"
DEF_KSNET = "data/external/kinase_substrate/renal_kinase_substrate_core.tsv"
CONTROL_KINASES = ("WNK", "SPAK_OSR1")


def log(msg: str) -> None:
    print(f"{datetime.now():%H:%M:%S} [v10] {msg}")


def sha256(path: Path) -> str | None:
    if not path.exists():
        return None
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def resolve(p: str | Path) -> Path:
    p = Path(p)
    return p if p.is_absolute() else REPO_ROOT / p


# Layer A -- KSEA kinase activity

def layer_a_kinase(phospho_path: Path, ksnet_path: Path, outdir: Path) -> dict:
    log("Layer A: KSEA kinase-activity inference on OSD-462 phosphoproteomics")
    sites = pd.read_csv(phospho_path, sep="\t")
    net = load_kinase_substrate_net(str(ksnet_path))
    table = ksea(sites, net, min_substrates=3)
    table.to_csv(outdir / "osd462_kinase_activity_summary.tsv", sep="\t", index=False)
    control_pass = ksea_positive_control_passes(table, control_kinases=CONTROL_KINASES)
    log(f"  kinases scored: {int(table['status'].eq('scored').sum())}; "
        f"WNK-SPAK/OSR1 positive control passed: {control_pass}")
    if not control_pass:
        log("  WARNING: positive control did not pass -- inspect the kinase-substrate net "
            "before trusting any kinase result.")
    return {
        "n_kinases_scored": int(table["status"].eq("scored").sum()),
        "positive_control_passed": bool(control_pass),
        "output": str(outdir / "osd462_kinase_activity_summary.tsv"),
    }


# Layer B -- TF / pathway activity (decoupler, network-dependent priors)

def _load_rna_effects(spec: dict[str, str]) -> pd.DataFrame:
    """Build a contrasts x genes matrix from per-cohort gene-effect tables.

    Each file must be a TSV with a gene-id column and a flight-effect column;
    the loader auto-detects common column names.
    """
    gene_cands = ["gene", "gene_id", "ensembl_gene_id", "ensembl"]
    eff_cands = ["effect", "logFC", "log2fc", "flight_effect", "iss_t_effect",
                 "stat", "t"]
    rows = {}
    for label, path in spec.items():
        df = pd.read_csv(resolve(path), sep=None, engine="python")
        gcol = next((c for c in df.columns if c.lower() in gene_cands), df.columns[0])
        ecol = next((c for c in df.columns if c.lower() in eff_cands), None)
        if ecol is None:
            raise ValueError(f"{label}: no recognizable effect column in {path}")
        s = df.set_index(gcol)[ecol].astype(float)
        s = s[~s.index.duplicated()]
        rows[label] = s
    mat = pd.DataFrame(rows).T  # contrasts x genes
    return mat.fillna(0.0)


def layer_b_tf_pathway(rna_spec: dict[str, str], outdir: Path) -> dict:
    log("Layer B: TF / pathway activity inference (decoupler ULM)")
    try:
        import decoupler as dc
    except Exception as exc:  # pragma: no cover
        log(f"  SKIPPED: decoupler not importable ({exc!r}).")
        return {"status": "skipped_no_decoupler"}
    try:
        progeny = dc.op.progeny(organism="mouse")
        collectri = dc.op.collectri(organism="mouse")
    except Exception as exc:
        log(f"  SKIPPED: could not fetch PROGENy/CollecTRI priors ({type(exc).__name__}). "
            "This step needs network access to omnipath/zenodo.")
        log("  Run recipe (on a machine with network):")
        log("    python3 -c \"import decoupler as dc; "
            "dc.op.progeny(organism='mouse'); dc.op.collectri(organism='mouse')\"  # warms cache")
        log("    then re-run this script with --rna-effects pointing at the cohort tables.")
        return {"status": "skipped_no_network_priors"}

    rna = _load_rna_effects(rna_spec)
    results = {}
    for name, net in (("progeny_pathway", progeny), ("collectri_tf", collectri)):
        long = run_ulm_activity(rna, net)
        long.to_csv(outdir / f"rna_{name}_activity.tsv", sep="\t", index=False)
        results[name] = str(outdir / f"rna_{name}_activity.tsv")
        log(f"  wrote {name}: {len(long)} rows")
    return {"status": "completed", "outputs": results}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phospho-sites", default=DEF_PHOSPHO)
    ap.add_argument("--ks-net", default=DEF_KSNET)
    ap.add_argument("--rna-effects", default="",
                    help="JSON: {cohort_label: gene-effect-table-path}")
    ap.add_argument("--outdir", default="")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    phospho = resolve(args.phospho_sites)
    ksnet = resolve(args.ks_net)
    outdir = resolve(args.outdir) if args.outdir else (
        REPO_ROOT / f"data/results/run_{datetime.now():%Y%m%d}_regulator_activity")
    outdir.mkdir(parents=True, exist_ok=True)

    layer_a = layer_a_kinase(phospho, ksnet, outdir)

    rna_spec = json.loads(args.rna_effects) if args.rna_effects else {}
    layer_b = (layer_b_tf_pathway(rna_spec, outdir) if rna_spec
               else {"status": "skipped_no_rna_effects_provided"})

    manifest = {
        "analysis": "v10 regulator-activity prioritization",
        "framing": "computational mechanism prioritization, not causal discovery",
        "timestamp": datetime.now().isoformat(),
        "inputs": {
            "phospho_sites": str(phospho), "phospho_sites_sha256": sha256(phospho),
            "kinase_substrate_net": str(ksnet), "kinase_substrate_net_sha256": sha256(ksnet),
            "rna_effects_spec": rna_spec,
        },
        "layer_a_kinase": layer_a,
        "layer_b_tf_pathway": layer_b,
        "positive_control": {
            "kinases": list(CONTROL_KINASES),
            "rule": "WNK and SPAK_OSR1 must return inferred_activity_down at p<0.05",
            "passed": layer_a.get("positive_control_passed"),
        },
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"Wrote outputs and manifest to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
