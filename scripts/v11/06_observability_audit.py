#!/usr/bin/env python3
"""v11 Module 3 — proteome observability-bias audit.

Pre-empts the first-line reviewer objection to the v11 RNA→protein
discordance: "the mismatch is just proteome detectability."

Pipeline:

  1. Build per-gene observability features from
     ``protein_effects_by_row.tsv`` (peptide-weighted ``n_channels_used``,
     ``missing_fraction``).
  2. Detectability gradient — fraction of RNA-detected genes also
     protein-quantified, per RNA-effect-magnitude decile.
  3. Observability-matched re-test — extend the matched-null stratum with
     a missing-fraction bin; rerun the Module 2 per-pathway propagation
     tests on the extended strata and report q-values alongside Module
     2's.
  4. High-coverage subset re-test — restrict to ``n_peptides >= 3`` and
     ``missing_fraction <= 0.2``; rerun the per-pathway propagation
     tests.
  5. NCC/SPAK phosphosite observability check — confirm the
     pre-specified regulatory sites are NOT in the low-observability
     tail of the phosphoproteome.

Usage::

    python scripts/v11/06_observability_audit.py [--n-null 10000]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.common import bh_fdr  # noqa: E402
from src.v11.matched_null import assign_match_strata  # noqa: E402
from src.v11.observability_audit import (  # noqa: E402
    assign_observability_strata,
    collapse_observability_to_gene,
    detectability_gradient,
    high_coverage_subset,
    merge_observability_into_pool,
    ncc_site_observability,
    propagation_with_strata,
)
from src.v11.rna_protein_propagation import (  # noqa: E402
    build_layer_table,
    load_gene_sets,
    pathway_gene_table,
    PropagationConfig,
)


DEFAULT_RUN_ROOT = REPO_ROOT / "data/results/run_20260606_v11_layer_specificity"
OSD462_ANCHOR = REPO_ROOT / "data/results/run_20260522_osd462_anchor/osd462_anchor"
DEFAULT_GENE_SETS = REPO_ROOT / "config/mechanism_gene_sets.yaml"
SEED = 20260607
N_NULL = 10000


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _bh_within_family(df: pd.DataFrame, p_col: str, q_col: str) -> pd.DataFrame:
    out = df.copy()
    out[q_col] = float("nan")
    ok = ~out["skipped"].fillna(False).astype(bool) & out[p_col].notna()
    if ok.any():
        out.loc[ok, q_col] = bh_fdr(out.loc[ok, p_col].astype(float).to_numpy())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    ap.add_argument("--anchor-dir", type=Path, default=OSD462_ANCHOR)
    ap.add_argument("--gene-sets", type=Path, default=DEFAULT_GENE_SETS)
    ap.add_argument("--n-null", type=int, default=N_NULL)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--peptide-filter", type=int, default=2)
    ap.add_argument("--high-coverage-min-peptides", type=int, default=3)
    ap.add_argument("--high-coverage-max-missing", type=float, default=0.20)
    args = ap.parse_args()

    out_dir = args.run_root / "observability"
    out_dir.mkdir(parents=True, exist_ok=True)
    manifests_dir = args.run_root / "manifests"
    manifests_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. inputs ---------------------------------------------------------
    flight_path = args.anchor_dir / "osd462_flight_effects.tsv"
    by_row_path = args.anchor_dir / "protein_effects_by_row.tsv"
    phospho_path = args.anchor_dir / "phospho_all_sites.tsv"
    print(f"[observability] anchor: {args.anchor_dir}")
    flight = pd.read_csv(flight_path, sep="\t")
    by_row = pd.read_csv(by_row_path, sep="\t")
    phospho = pd.read_csv(phospho_path, sep="\t")
    print(f"[observability] flight rows: {len(flight)}; by-row rows: {len(by_row)};"
          f" phospho sites: {len(phospho)}")

    # ---- 2. per-gene observability features --------------------------------
    obs = collapse_observability_to_gene(by_row)
    obs.to_csv(out_dir / "v11_observability_per_gene.tsv", sep="\t", index=False)
    print(f"[observability] per-gene observability: n_genes={len(obs)};"
          f" median_missing_fraction={obs['missing_fraction'].median():.4f}")

    # ---- 3. build Module 2's pool with observability merged in -------------
    cfg = PropagationConfig(
        run_root=args.run_root, anchor_dir=args.anchor_dir,
        gene_sets_path=args.gene_sets, n_null=args.n_null, seed=args.seed,
        peptide_filter=args.peptide_filter,
    )
    layer = build_layer_table(flight, phospho)
    gene_sets = load_gene_sets(args.gene_sets)
    genes = pathway_gene_table(layer, gene_sets, cfg)

    # Pool with observability columns added; keep Module 2's match_stratum and
    # rebuild a Module 3 (abundance × peptide × missing-fraction) stratum.
    pool = flight.dropna(subset=["osd462_rna_effect", "protein_flight_effect"]).copy()
    pool["n_peptides"] = pd.to_numeric(pool["n_peptides"], errors="coerce")
    pool = pool[pool["n_peptides"].fillna(0) >= args.peptide_filter].copy()
    pool = pool.reset_index(drop=True)
    pool["gene_upper"] = pool["gene_symbol"].astype(str).str.upper()
    pool = merge_observability_into_pool(pool, obs)
    pool["missing_fraction"] = pool["missing_fraction"].fillna(0.0)
    # Rebuild Module 2-equivalent strata on the trimmed pool for a fair comparison.
    pool["match_stratum"] = assign_match_strata(pool).values
    pool["match_stratum_observability"] = assign_observability_strata(pool).values
    print(f"[observability] pool size: {len(pool)};"
          f" unique standard strata: {pool['match_stratum'].nunique()};"
          f" unique extended strata: {pool['match_stratum_observability'].nunique()}")

    # ---- 4. detectability gradient -----------------------------------------
    rna_table = flight[["ENSEMBL", "rrrm2_iss_t_rna_effect"]].dropna(
        subset=["rrrm2_iss_t_rna_effect"]
    ).copy()
    grad = detectability_gradient(rna_table, pool)
    grad.to_csv(out_dir / "v11_observability_detectability_gradient.tsv", sep="\t",
                index=False)
    print("[observability] detectability gradient (RNA |effect| decile → fraction protein-quantified):")
    print(grad.to_string(index=False))

    # ---- 5. observability-matched per-pathway re-test ----------------------
    standard = propagation_with_strata(
        pool, gene_sets, stratum_col="match_stratum",
        n_null=args.n_null, seed=args.seed,
    )
    standard = _bh_within_family(standard, "protein_slope_p_two_sided",
                                 "protein_slope_q_two_sided")
    standard = _bh_within_family(standard, "signed_mean_p_greater",
                                 "signed_mean_q_greater")
    standard["stratum"] = "abundance_x_peptide"

    observ = propagation_with_strata(
        pool, gene_sets, stratum_col="match_stratum_observability",
        n_null=args.n_null, seed=args.seed,
    )
    observ = _bh_within_family(observ, "protein_slope_p_two_sided",
                               "protein_slope_q_two_sided")
    observ = _bh_within_family(observ, "signed_mean_p_greater",
                               "signed_mean_q_greater")
    observ["stratum"] = "abundance_x_peptide_x_missingness"

    audit = pd.concat([standard, observ], ignore_index=True)
    audit_path = out_dir / "v11_observability_matched_propagation.tsv"
    audit.to_csv(audit_path, sep="\t", index=False)

    # Delta table — q-value shift per pathway between the two stratifications.
    delta_rows = []
    for pathway in standard["pathway"]:
        s_row = standard[standard["pathway"].eq(pathway)].iloc[0]
        o_row = observ[observ["pathway"].eq(pathway)].iloc[0]
        delta_rows.append({
            "pathway": pathway,
            "n_quantified_standard": s_row["n_quantified"],
            "n_quantified_observability": o_row["n_quantified"],
            "protein_slope_standard": s_row.get("protein_slope", float("nan")),
            "protein_slope_observability": o_row.get("protein_slope", float("nan")),
            "protein_slope_q_standard": s_row.get("protein_slope_q_two_sided", float("nan")),
            "protein_slope_q_observability": o_row.get("protein_slope_q_two_sided", float("nan")),
            "signed_mean_standard": s_row.get("signed_mean_protein_by_rna", float("nan")),
            "signed_mean_observability": o_row.get("signed_mean_protein_by_rna", float("nan")),
            "signed_mean_q_standard": s_row.get("signed_mean_q_greater", float("nan")),
            "signed_mean_q_observability": o_row.get("signed_mean_q_greater", float("nan")),
        })
    delta = pd.DataFrame(delta_rows)
    delta_path = out_dir / "v11_observability_q_delta.tsv"
    delta.to_csv(delta_path, sep="\t", index=False)
    print(f"[observability] wrote {audit_path}")
    print("[observability] q-value delta (standard → observability-matched):")
    print(delta[[
        "pathway", "protein_slope_q_standard", "protein_slope_q_observability",
        "signed_mean_q_standard", "signed_mean_q_observability",
    ]].round(4).to_string(index=False))

    # ---- 6. high-coverage subset re-test -----------------------------------
    subset = high_coverage_subset(pool,
                                  min_peptides=args.high_coverage_min_peptides,
                                  max_missing_fraction=args.high_coverage_max_missing)
    print(f"[observability] high-coverage subset: {len(subset)}/{len(pool)} genes")
    if len(subset) >= 30:
        subset = subset.copy()
        subset["match_stratum"] = assign_match_strata(subset).values
        hc = propagation_with_strata(subset, gene_sets, stratum_col="match_stratum",
                                     n_null=args.n_null, seed=args.seed)
        hc = _bh_within_family(hc, "protein_slope_p_two_sided",
                               "protein_slope_q_two_sided")
        hc = _bh_within_family(hc, "signed_mean_p_greater", "signed_mean_q_greater")
        hc["subset_definition"] = (
            f"n_peptides >= {args.high_coverage_min_peptides} & "
            f"missing_fraction <= {args.high_coverage_max_missing}"
        )
        hc_path = out_dir / "v11_observability_high_coverage_propagation.tsv"
        hc.to_csv(hc_path, sep="\t", index=False)
        print(f"[observability] wrote {hc_path}")
    else:
        print("[observability] high-coverage subset too small; skipping re-test")
        hc_path = None

    # ---- 7. NCC/SPAK phosphosite observability check -----------------------
    site_audit = ncc_site_observability(phospho)
    site_path = out_dir / "v11_observability_ncc_site_audit.tsv"
    site_audit.to_csv(site_path, sep="\t", index=False)
    print(f"[observability] NCC site audit:")
    print(site_audit.round(3).to_string(index=False))

    # ---- 8. manifest -------------------------------------------------------
    manifest = {
        "analysis": "v11 Module 3 — proteome observability-bias audit",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "osd462_flight_effects": {"path": str(flight_path), "sha256": _sha256(flight_path)},
            "protein_effects_by_row": {"path": str(by_row_path), "sha256": _sha256(by_row_path)},
            "phospho_all_sites": {"path": str(phospho_path), "sha256": _sha256(phospho_path)},
            "mechanism_gene_sets": {"path": str(args.gene_sets), "sha256": _sha256(args.gene_sets)},
        },
        "outputs": {
            "per_gene_observability": str(out_dir / "v11_observability_per_gene.tsv"),
            "detectability_gradient": str(out_dir / "v11_observability_detectability_gradient.tsv"),
            "matched_propagation": str(audit_path),
            "q_delta": str(delta_path),
            "high_coverage_propagation": str(hc_path) if hc_path else None,
            "ncc_site_audit": str(site_path),
        },
        "parameters": {
            "n_null": args.n_null,
            "seed": args.seed,
            "peptide_filter": args.peptide_filter,
            "high_coverage_min_peptides": args.high_coverage_min_peptides,
            "high_coverage_max_missing_fraction": args.high_coverage_max_missing,
            "standard_null_model": "abundance(5) × peptide(4) matched random gene sets",
            "observability_null_model": "abundance(5) × peptide(4) × missing-fraction(3) matched random gene sets",
        },
        "n_pool": int(len(pool)),
        "n_unique_strata_standard": int(pool["match_stratum"].nunique()),
        "n_unique_strata_observability": int(pool["match_stratum_observability"].nunique()),
        "n_high_coverage_subset": int(len(subset)) if 'subset' in locals() else 0,
    }
    manifest_path = manifests_dir / "v11_observability_audit_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(f"[observability] wrote manifest {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
