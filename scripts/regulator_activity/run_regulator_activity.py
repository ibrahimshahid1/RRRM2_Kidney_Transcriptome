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
from src.multiomics.osd462_stage0 import (  # noqa: E402
    audit_ncc_spak_phosphosites,
    isolated_canonical_assay_features,
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
    provenance = isolated_canonical_assay_features(
        audit_ncc_spak_phosphosites()
    )
    eligible_keys = set(
        zip(
            provenance["gene_symbol"].astype(str),
            provenance["workbook_site_position"].astype(str),
        )
    )
    assay_net = net[
        [
            (str(gene), str(site)) in eligible_keys
            for gene, site in zip(net["substrate_gene"], net["substrate_site"])
        ]
    ].copy()
    if assay_net.empty:
        # Fail closed. A position match to T53 or S383 is insufficient because
        # both OSD-462 rows are co-modified phosphoforms.
        table = pd.DataFrame(
            {
                "kinase": sorted(net["kinase"].unique()),
                "n_substrates_quantified": 0,
                "mean_substrate_effect": np.nan,
                "z_score": np.nan,
                "p_value": np.nan,
                "status": "no_qualified_isolated_substrates",
                "q_value": np.nan,
                "direction": "not_inferred",
                "background_n_sites": int(len(sites)),
            }
        )
    else:
        table = ksea(sites, assay_net, min_substrates=3)
    table.to_csv(outdir / "osd462_kinase_activity_summary.tsv", sep="\t", index=False)
    control_pass = ksea_positive_control_passes(table, control_kinases=CONTROL_KINASES)
    log(f"  kinases scored: {int(table['status'].eq('scored').sum())}; "
        f"WNK-SPAK/OSR1 positive control passed: {control_pass}")
    if not control_pass:
        log(
            "  WARNING: kinase activity not inferred -- no isolated canonical "
            "OSD-462 substrate phosphoform passes provenance."
        )
    return {
        "n_kinases_scored": int(table["status"].eq("scored").sum()),
        "n_qualified_isolated_substrate_features": int(len(assay_net)),
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
            "rule": (
                "Only residue- and phosphoform-qualified isolated substrates "
                "may enter KSEA; otherwise activity is not inferred."
            ),
            "passed": layer_a.get("positive_control_passed"),
        },
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    log(f"Wrote outputs and manifest to {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
