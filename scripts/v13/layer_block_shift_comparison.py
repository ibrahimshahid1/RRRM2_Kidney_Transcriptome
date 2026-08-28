#!/usr/bin/env python3
"""Is the flight-block signal dip phospho-specific, or shared with protein?

Question
--------
The OSD-462 channel profile shows a block-shaped dip across the flight
reporter channels (129N-131) in both plexes, recovering immediately at 131C.
A within-block positional slope cannot see a step of that shape, so two
explanations remain:

  (a) a block-level technical effect (sample handling, loading, labelling
      batch) that happens to align with condition; or
  (b) a genuine reduction in phosphopeptide signal in flight animals.

Discriminator
-------------
The protein workbook measures the *same animals* in the *same channel layout*
but without Fe-NTA phosphopeptide enrichment. A handling/loading effect
upstream of the phospho/protein split should appear in both layers at similar
magnitude. Phospho-specific suppression should not.

Two comparisons are made:

1. **Marginal.** Flight-block minus ground-block mean centred log2 signal, per
   layer and plex.
2. **Paired by parent protein.** For every protein quantified in both layers,
   the phosphosite flight-minus-ground effect minus the same protein's own
   flight-minus-ground effect. This removes anything shared by the two layers
   for that protein and is the comparison a "phosphorylation changes without
   abundance changing" claim actually requires.

Boundary
--------
The two layers were labelled in separate reactions (tc882-883 phospho,
tc884-885 protein), so a labelling-batch effect could in principle differ
between them. This test discriminates shared upstream handling effects, not
every possible technical explanation.

Usage
-----
    python3 scripts/v13/layer_block_shift_comparison.py
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
EXT = REPO / "data/external/osdr/OSD-462"
PHOSPHO = EXT / "GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx"
PROTEIN = EXT / "GLDS-462_proteomics_2021-12-31_tc884-885_Protein_WorkUp.xlsx"
DEFAULT_OUT = REPO / "data/results/run_20260729_osd462_layer_block_shift"

# sheet, first data row, gene column, {plex: first summed-S/N column}
LAYERS = {
    "phospho": dict(path=PHOSPHO, sheet="siteQuant_360", first_row=4,
                    gene_col=2, cols={1: 11, 2: 26}),
    "protein": dict(path=PROTEIN, sheet="protein_quant_2721", first_row=4,
                    gene_col=2, cols={1: 7, 2: 22}),
}

FLIGHT = slice(5, 10)    # 129N, 129C, 130N, 130C, 131
GROUND = slice(10, 15)   # 131C, 132N, 132C, 133N, 133C
BASELINE = slice(0, 5)   # 126, 127N, 127C, 128N, 128C

SEED = 20260729


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_layer(spec: dict, cache: pathlib.Path):
    """Return (genes, {plex: (n, 15) raw S/N})."""
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        return list(z["genes"]), {1: z["plex1"], 2: z["plex2"]}

    import openpyxl

    wb = openpyxl.load_workbook(spec["path"], read_only=True, data_only=True)
    ws = wb[spec["sheet"]]
    c1, c2 = spec["cols"][1], spec["cols"][2]
    gcol = spec["gene_col"]
    lo = min(gcol, c1, c2)
    hi = max(c1, c2) + 14
    genes, vals = [], []
    for row in ws.iter_rows(min_row=spec["first_row"], min_col=lo, max_col=hi,
                            values_only=True):
        genes.append(row[gcol - lo])
        block = [row[c1 - lo + k] for k in range(15)] + \
                [row[c2 - lo + k] for k in range(15)]
        vals.append([np.nan if v is None else float(v) for v in block])
    wb.close()

    arr = np.array(vals, dtype=float)
    genes = np.array([("" if g is None else str(g)) for g in genes], dtype=object)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, genes=genes, plex1=arr[:, :15], plex2=arr[:, 15:])
    return list(genes), {1: arr[:, :15], 2: arr[:, 15:]}


def centred(mat: np.ndarray):
    ok = np.all(np.isfinite(mat) & (mat > 0), axis=1)
    x = np.log2(mat[ok])
    return ok, x - x.mean(axis=1, keepdims=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args(argv)

    for spec in LAYERS.values():
        if not spec["path"].exists():
            print(f"missing: {spec['path']}", file=sys.stderr)
            return 1
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    data = {}
    for layer, spec in LAYERS.items():
        genes, mats = load_layer(spec, args.out / f"_{layer}_cache.npz")
        data[layer] = (np.array(genes, dtype=object), mats)

    # ---- 1. marginal block profile per layer and plex ----
    marginal, per_feature = [], {}
    for layer in LAYERS:
        genes, mats = data[layer]
        for plex in (1, 2):
            ok, z = centred(mats[plex])
            fl = z[:, FLIGHT].mean(axis=1)
            gc = z[:, GROUND].mean(axis=1)
            bl = z[:, BASELINE].mean(axis=1)
            marginal.append(
                {
                    "layer": layer,
                    "plex": plex,
                    "n_features": int(z.shape[0]),
                    "baseline_block_mean_log2": round(float(bl.mean()), 6),
                    "flight_block_mean_log2": round(float(fl.mean()), 6),
                    "ground_block_mean_log2": round(float(gc.mean()), 6),
                    "flight_minus_ground_log2": round(float((fl - gc).mean()), 6),
                    "flight_minus_baseline_log2": round(float((fl - bl).mean()), 6),
                    "fraction_features_negative": round(
                        float(np.mean((fl - gc) < 0)), 4
                    ),
                }
            )
            per_feature[(layer, plex)] = (genes[ok], fl - gc)

    # ---- 2. paired by parent protein, averaged over plexes ----
    def gene_effect(layer):
        acc = {}
        for plex in (1, 2):
            g, e = per_feature[(layer, plex)]
            for sym, val in zip(g, e):
                if sym:
                    acc.setdefault(sym, []).append(float(val))
        return {k: float(np.median(v)) for k, v in acc.items()}

    ph, pr = gene_effect("phospho"), gene_effect("protein")
    shared = sorted(set(ph) & set(pr))
    d = np.array([ph[g] - pr[g] for g in shared])

    n_pos = int(np.count_nonzero(d > 0))
    n_neg = int(np.count_nonzero(d < 0))
    # exact two-sided sign test
    from math import comb
    n = n_pos + n_neg
    k = min(n_pos, n_neg)
    tail = sum(comb(n, i) for i in range(k + 1)) / (2 ** n)
    sign_p = min(1.0, 2 * tail)

    boot = np.array([rng.choice(d, d.size, replace=True).mean()
                     for _ in range(10000)])
    paired = [
        {
            "n_shared_parent_proteins": len(shared),
            "mean_phospho_minus_protein_log2": round(float(d.mean()), 6),
            "median_phospho_minus_protein_log2": round(float(np.median(d)), 6),
            "bootstrap_ci_low": round(float(np.percentile(boot, 2.5)), 6),
            "bootstrap_ci_high": round(float(np.percentile(boot, 97.5)), 6),
            "n_more_negative_in_phospho": n_neg,
            "n_more_negative_in_protein": n_pos,
            "exact_sign_test_p_two_sided": float(f"{sign_p:.3e}"),
            "mean_phospho_flight_minus_ground": round(
                float(np.mean([ph[g] for g in shared])), 6),
            "mean_protein_flight_minus_ground": round(
                float(np.mean([pr[g] for g in shared])), 6),
        }
    ]

    def write(name, rows):
        p = args.out / name
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return p

    p1 = write("marginal_block_profile_by_layer.tsv", marginal)
    p2 = write("paired_parent_protein_layer_contrast.tsv", paired)

    (args.out / "manifest.json").write_text(
        json.dumps(
            {
                "purpose": "Discriminate a shared block-level technical shift "
                           "from phospho-specific suppression in OSD-462.",
                "inputs": {
                    str(PHOSPHO.relative_to(REPO)): sha256(PHOSPHO),
                    str(PROTEIN.relative_to(REPO)): sha256(PROTEIN),
                },
                "quantification": "raw summed signal-to-noise, log2, "
                                  "feature-centred within plex",
                "seed": args.seed,
                "outputs": {p.name: sha256(p) for p in (p1, p2)},
                "boundary": "Phospho and protein were labelled in separate "
                            "reactions, so this discriminates shared upstream "
                            "handling effects, not all technical explanations. "
                            "Condition remains aliased with reporter-tag block "
                            "in both layers.",
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {p1}\nwrote {p2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
