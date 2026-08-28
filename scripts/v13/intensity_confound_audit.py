#!/usr/bin/env python3
"""Intensity-confound audit of the v13 continuous phosphosite enrichment.

Motivation
----------
The v13 exact run standardises every parent gene against its own balanced-label
null, which controls gene-specific *variance* but not a systematic
*intensity-dependent shift* in the observed effects. This script measures that
shift directly and re-tests every frozen gene set against an
intensity-stratified competitive null.

It consumes only emitted artefacts of
``run_20260729_v13_continuous_phospho_exact_final`` -- it does not refit the
phosphosite models -- so it is cheap, deterministic, and auditable.

Outputs
-------
``intensity_decile_gradient.tsv``
    Mean gene-level Z by decile of median phosphopeptide signal, per profile.
``intensity_stratified_set_enrichment.tsv``
    Per gene set and profile: raw competitive statistic, statistic after
    within-stratum centring, and an intensity-stratified gene-label
    permutation p-value.
``manifest.json``
    Inputs, digests, parameters, seed.

Usage
-----
    venv/bin/python scripts/v13/intensity_confound_audit.py
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import pathlib
import sys
from typing import Dict, List, Sequence

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
SOURCE_RUN = REPO / "data/results/run_20260729_v13_continuous_phospho_exact_final"
DEFAULT_OUT = REPO / "data/results/run_20260729_v13_intensity_confound_audit"

CALIBRATION = SOURCE_RUN / "parent_gene_null_calibration.tsv"
MEMBERSHIP = SOURCE_RUN / "gene_set_membership_frozen_copy.tsv"

# Frozen primary canonical-axis exclusion, matching the v13 configuration.
PRIMARY_CANONICAL_EXCLUSION = ("Slc12a3", "Stk39", "Oxsr1", "Wnk4")

GENE_SCORE = "median_negative_site_effect"

SET_DISPLAY_ORDER = (
    "DCT2_CNT_transition",
    "ASDN",
    "DCT1",
    "thick_ascending_limb",
    "proximal_tubule",
    "principal_cell",
    "intercalated_cell",
    "podocyte",
    "endothelial",
    "fibroblast",
    "immune",
)

N_STRATA = 20
N_PERMUTATIONS = 20000
SEED = 20260729
MIN_SET_GENES = 3


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def load_gene_sets(path: pathlib.Path) -> Dict[str, set]:
    sets: Dict[str, set] = collections.defaultdict(set)
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sets[row["gene_set"]].add(row["gene_symbol"])
    return dict(sets)


def load_profile(rows: Sequence[dict], profile: str, exclusion: Sequence[str]):
    excl = set(exclusion)
    genes: List[str] = []
    z: List[float] = []
    signal: List[float] = []
    for row in rows:
        if row["profile"] != profile:
            continue
        if row["gene_score"] != GENE_SCORE:
            continue
        if row["fixed_null_universe_eligible"] != "True":
            continue
        if row["gene_symbol"] in excl:
            continue
        try:
            zi = float(row["observed_gene_z"])
            si = float(row["median_log2_signal"])
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(zi) and np.isfinite(si)):
            continue
        genes.append(row["gene_symbol"])
        z.append(zi)
        signal.append(si)
    return np.asarray(genes), np.asarray(z), np.asarray(signal)


def competitive_statistic(z: np.ndarray, total: float, n_total: int, idx: np.ndarray) -> float:
    inside = z[idx].sum()
    n_in = idx.size
    return float(inside / n_in - (total - inside) / (n_total - n_in))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    parser.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument(
        "--profiles",
        nargs="*",
        default=None,
        help="Restrict to these inference profiles (default: all).",
    )
    args = parser.parse_args(argv)

    for path in (CALIBRATION, MEMBERSHIP):
        if not path.exists():
            print(f"missing required input: {path}", file=sys.stderr)
            return 1

    args.out.mkdir(parents=True, exist_ok=True)

    with CALIBRATION.open() as fh:
        calibration_rows = list(csv.DictReader(fh, delimiter="\t"))
    gene_sets = load_gene_sets(MEMBERSHIP)

    profiles = sorted(
        {r["profile"] for r in calibration_rows if r["gene_score"] == GENE_SCORE}
    )
    if args.profiles:
        keep = set(args.profiles)
        profiles = [p for p in profiles if p in keep]

    gradient_rows: List[dict] = []
    set_rows: List[dict] = []
    rng = np.random.default_rng(args.seed)

    for profile in profiles:
        genes, z, signal = load_profile(
            calibration_rows, profile, PRIMARY_CANONICAL_EXCLUSION
        )
        if genes.size < 200:
            continue

        index = {g: i for i, g in enumerate(genes)}
        n_total = genes.size
        total = float(z.sum())
        rho = spearman(z, signal)

        order = np.argsort(signal, kind="stable")
        stratum = np.empty(n_total, dtype=int)
        for bin_id, members in enumerate(np.array_split(order, N_STRATA)):
            stratum[members] = bin_id
        by_stratum = [np.flatnonzero(stratum == b) for b in range(N_STRATA)]

        centred = z.copy()
        for b in range(N_STRATA):
            mask = stratum == b
            centred[mask] = z[mask] - z[mask].mean()
        centred_total = float(centred.sum())

        for decile, members in enumerate(np.array_split(order, 10), start=1):
            gradient_rows.append(
                {
                    "profile": profile,
                    "signal_decile": decile,
                    "n_genes": members.size,
                    "mean_median_log2_signal": round(float(signal[members].mean()), 6),
                    "mean_gene_z": round(float(z[members].mean()), 6),
                    "spearman_z_vs_signal_profile": round(rho, 6),
                }
            )

        for name in SET_DISPLAY_ORDER:
            members = gene_sets.get(name, set())
            idx = np.asarray(
                [index[g] for g in sorted(members) if g in index], dtype=int
            )
            row = {
                "profile": profile,
                "gene_set": name,
                "n_observable_genes": int(idx.size),
                "n_background_genes": int(n_total - idx.size),
                "mean_median_log2_signal_set": "",
                "raw_competitive_statistic": "",
                "within_stratum_centred_statistic": "",
                "intensity_stratified_p_greater": "",
                "n_permutations": args.permutations,
                "evaluation_status": "insufficient_observable_genes",
            }
            if idx.size >= MIN_SET_GENES:
                raw = competitive_statistic(z, total, n_total, idx)
                adjusted = competitive_statistic(
                    centred, centred_total, n_total, idx
                )
                need = collections.Counter(stratum[idx].tolist())
                # Vectorised stratified sampling without replacement: the
                # competitive statistic is a linear function of the sampled
                # sum, so only the per-permutation sum is required.
                sampled_sum = np.zeros(args.permutations)
                for s, c in need.items():
                    pool = z[by_stratum[s]]
                    if c == 1:
                        draw = rng.integers(0, pool.size, size=args.permutations)
                        sampled_sum += pool[draw]
                        continue
                    keys = rng.random((args.permutations, pool.size))
                    take = np.argpartition(keys, c - 1, axis=1)[:, :c]
                    sampled_sum += pool[take].sum(axis=1)
                n_in = idx.size
                null = sampled_sum / n_in - (total - sampled_sum) / (
                    n_total - n_in
                )
                p_value = float(
                    (1 + np.count_nonzero(null >= raw)) / (args.permutations + 1)
                )
                row.update(
                    {
                        "mean_median_log2_signal_set": round(
                            float(signal[idx].mean()), 6
                        ),
                        "raw_competitive_statistic": round(raw, 6),
                        "within_stratum_centred_statistic": round(adjusted, 6),
                        "intensity_stratified_p_greater": round(p_value, 6),
                        "evaluation_status": "evaluated",
                    }
                )
            set_rows.append(row)

    gradient_path = args.out / "intensity_decile_gradient.tsv"
    with gradient_path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=list(gradient_rows[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(gradient_rows)

    set_path = args.out / "intensity_stratified_set_enrichment.tsv"
    with set_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(set_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(set_rows)

    manifest = {
        "purpose": (
            "Measure intensity-dependent bias in v13 gene-level suppression "
            "scores and re-test frozen gene sets against an "
            "intensity-stratified competitive null."
        ),
        "source_run": str(SOURCE_RUN.relative_to(REPO)),
        "inputs": {
            str(CALIBRATION.relative_to(REPO)): sha256(CALIBRATION),
            str(MEMBERSHIP.relative_to(REPO)): sha256(MEMBERSHIP),
        },
        "parameters": {
            "gene_score": GENE_SCORE,
            "primary_canonical_exclusion": list(PRIMARY_CANONICAL_EXCLUSION),
            "n_strata": N_STRATA,
            "n_permutations": args.permutations,
            "min_set_genes": MIN_SET_GENES,
            "seed": args.seed,
        },
        "outputs": {
            gradient_path.name: sha256(gradient_path),
            set_path.name: sha256(set_path),
        },
        "boundary": (
            "This audit does not remove reporter-tag aliasing and does not "
            "assign phosphosites to cells of origin. It only asks whether set "
            "level enrichment survives matching on phosphopeptide signal."
        ),
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print(f"wrote {gradient_path}")
    print(f"wrote {set_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
