#!/usr/bin/env python3
"""Within-block reporter-position diagnostic for OSD-462.

Why this is the clean version
-----------------------------
Condition is perfectly aliased with reporter-tag block (baseline 126-128C,
flight 129N-131N, ground 131C-133C, identical in both plexes), so a
between-block comparison cannot separate biology from tag position.

But *within* a block the five channels hold five biologically exchangeable
animals of the same condition. Any systematic trend across channel position
within a block is therefore a pure tag-position effect with no biological
confound. Six independent estimates are available: 3 blocks x 2 plexes.

If a within-block slope exists, extrapolating it across the block boundary
bounds how much of the flight-minus-ground contrast reporter position alone
could produce.

Outputs
-------
``within_block_position_slopes.tsv``  per block x plex slope and permutation p
``position_effect_bound.tsv``         implied flight-vs-ground positional shift
``channel_profile.tsv``               mean centred log2 signal per channel
``manifest.json``

Usage
-----
    python3 scripts/v13/reporter_position_diagnostic.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import sys

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
PHOSPHO = (
    REPO
    / "data/external/osdr/OSD-462"
    / "GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx"
)
DEFAULT_OUT = REPO / "data/results/run_20260729_osd462_reporter_position"

SHEET = "siteQuant_360"
HEADER_ROW = 3
FIRST_DATA_ROW = 4

# 1-based worksheet columns holding raw signal-to-noise sums.
PLEX_COLUMNS = {1: list(range(11, 26)), 2: list(range(26, 41))}

CHANNELS = [
    "126", "127n", "127c", "128n", "128c",
    "129n", "129c", "130n", "130c", "131",
    "131c", "132n", "132c", "133n", "133c",
]
BLOCKS = [("baseline", 0, 5), ("flight", 5, 10), ("ground", 10, 15)]

N_PERMUTATIONS = 20000
SEED = 20260729


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_matrices(cache: pathlib.Path):
    """Return {plex: (n_sites, 15) array of raw S/N}, cached as .npz."""
    if cache.exists():
        z = np.load(cache)
        return {1: z["plex1"], 2: z["plex2"]}

    import openpyxl

    wb = openpyxl.load_workbook(PHOSPHO, read_only=True, data_only=True)
    ws = wb[SHEET]
    wanted = PLEX_COLUMNS[1] + PLEX_COLUMNS[2]
    lo, hi = min(wanted), max(wanted)
    rows = []
    for row in ws.iter_rows(min_row=FIRST_DATA_ROW, min_col=lo, max_col=hi,
                            values_only=True):
        rows.append([row[c - lo] for c in wanted])
    wb.close()

    arr = np.array(
        [[np.nan if v is None else float(v) for v in r] for r in rows],
        dtype=float,
    )
    plex1, plex2 = arr[:, :15], arr[:, 15:]
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, plex1=plex1, plex2=plex2)
    return {1: plex1, 2: plex2}


def centred_log2(mat: np.ndarray) -> np.ndarray:
    """Keep complete positive rows; log2 then centre each site across channels."""
    ok = np.all(np.isfinite(mat) & (mat > 0), axis=1)
    x = np.log2(mat[ok])
    return x - x.mean(axis=1, keepdims=True)


def slope_vs_position(block: np.ndarray) -> float:
    """Mean per-site OLS slope on within-block position 1..5."""
    pos = np.arange(1, block.shape[1] + 1, dtype=float)
    pos = pos - pos.mean()
    denom = float((pos ** 2).sum())
    return float((block * pos).sum(axis=1).mean() / denom)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    ap.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args(argv)

    if not PHOSPHO.exists():
        print(f"missing workbook: {PHOSPHO}", file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    mats = load_matrices(args.out / "_sn_cache.npz")

    slope_rows, profile_rows = [], []
    block_means = {}

    for plex in (1, 2):
        z = centred_log2(mats[plex])
        n_sites = z.shape[0]

        for j, ch in enumerate(CHANNELS):
            profile_rows.append(
                {
                    "plex": plex,
                    "channel_order": j + 1,
                    "reporter_tag": ch,
                    "block": next(b for b, s, e in BLOCKS if s <= j < e),
                    "within_block_position": j - next(
                        s for b, s, e in BLOCKS if s <= j < e
                    ) + 1,
                    "n_sites": n_sites,
                    "mean_centred_log2": round(float(z[:, j].mean()), 6),
                }
            )

        for name, start, end in BLOCKS:
            block = z[:, start:end]
            block_means[(plex, name)] = float(block.mean())
            observed = slope_vs_position(block)

            null = np.empty(args.permutations)
            for b in range(args.permutations):
                perm = rng.permutation(block.shape[1])
                null[b] = slope_vs_position(block[:, perm])
            p_two = float(
                (1 + np.count_nonzero(np.abs(null) >= abs(observed)))
                / (args.permutations + 1)
            )
            slope_rows.append(
                {
                    "plex": plex,
                    "block": name,
                    "n_sites": n_sites,
                    "slope_log2_per_channel_step": round(observed, 6),
                    "null_sd": round(float(null.std()), 6),
                    "permutation_p_two_sided": round(p_two, 6),
                    "n_permutations": args.permutations,
                }
            )

    # Pooled within-block slope, and what it implies across the FL->GC boundary.
    pooled = float(np.mean([r["slope_log2_per_channel_step"] for r in slope_rows]))
    bound_rows = []
    for plex in (1, 2):
        observed_diff = block_means[(plex, "flight")] - block_means[(plex, "ground")]
        # Flight block centre is channel 8, ground block centre is channel 13.
        predicted = pooled * (8 - 13)
        bound_rows.append(
            {
                "plex": plex,
                "observed_flight_minus_ground_log2": round(observed_diff, 6),
                "pooled_within_block_slope": round(pooled, 6),
                "channel_steps_between_block_centres": 5,
                "positional_prediction_log2": round(predicted, 6),
                "fraction_of_observed_explained": (
                    round(predicted / observed_diff, 4)
                    if observed_diff != 0
                    else ""
                ),
            }
        )

    def write(name, rows):
        path = args.out / name
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return path

    p1 = write("within_block_position_slopes.tsv", slope_rows)
    p2 = write("position_effect_bound.tsv", bound_rows)
    p3 = write("channel_profile.tsv", profile_rows)

    (args.out / "manifest.json").write_text(
        json.dumps(
            {
                "purpose": (
                    "Estimate reporter-tag position effects using within-block "
                    "comparisons among biologically exchangeable animals."
                ),
                "input": {str(PHOSPHO.relative_to(REPO)): sha256(PHOSPHO)},
                "sheet": SHEET,
                "quantification": "raw signal-to-noise sums, log2, site-centred",
                "parameters": {
                    "permutations": args.permutations,
                    "seed": args.seed,
                },
                "outputs": {p.name: sha256(p) for p in (p1, p2, p3)},
                "boundary": (
                    "Within-block slopes are free of biological confounding "
                    "because animals within a block share a condition. The "
                    "extrapolation across block boundaries assumes the "
                    "positional trend continues linearly and is an estimate, "
                    "not a correction."
                ),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {p1}\nwrote {p2}\nwrote {p3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
