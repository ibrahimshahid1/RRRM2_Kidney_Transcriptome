#!/usr/bin/env python3
"""DCT transport score check in the GSE269622 Visium IRI reference."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.v11.spatial_reference_projection import (
    NICHE_MARKERS,
    RUN_ROOT,
    TIME_ORDER,
    VISIUM_DIR,
    dct_adjacent_mask,
    load_visium_sample,
    score_genes,
)


DCT_TRANSPORT_GENES = ["Slc12a3", "Stk39", "Wnk1", "Wnk4", "Oxsr1", "Kcnj10", "Kcnj16", "Clcnkb", "Bsnd", "Calb1"]


def zscore(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    sd = x.std(skipna=True, ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return x * np.nan
    return (x - x.mean(skipna=True)) / sd


def dct_transport_score(adata) -> pd.Series:
    gene_to_idx = {str(g).lower(): i for i, g in enumerate(adata.var_names)}
    present = [gene_to_idx[g.lower()] for g in DCT_TRANSPORT_GENES if g.lower() in gene_to_idx]
    if not present:
        return pd.Series(np.nan, index=adata.obs_names)
    mat = adata.X[:, present]
    arr = np.asarray(mat.toarray() if hasattr(mat, "toarray") else mat, dtype=float)
    z = np.apply_along_axis(lambda v: zscore(pd.Series(v)).to_numpy(dtype=float), 0, arr)
    return pd.Series(np.nanmean(z, axis=1), index=adata.obs_names)


def summarize(values: pd.Series, mask: pd.Series, label: str, timepoint: str) -> dict:
    vals = values.loc[mask.index[mask]].dropna()
    return {
        "timepoint": timepoint,
        "group": label,
        "n_spots": int(len(vals)),
        "mean_dct_transport_score": float(vals.mean()) if len(vals) else np.nan,
        "median_dct_transport_score": float(vals.median()) if len(vals) else np.nan,
        "sd_dct_transport_score": float(vals.std(ddof=1)) if len(vals) > 1 else np.nan,
    }


def run(run_root: Path) -> dict:
    out_dir = run_root / "spatial_reference"
    out_dir.mkdir(parents=True, exist_ok=True)
    extracted_root = out_dir / "visium_extracted_minimal"
    archives = sorted(VISIUM_DIR.glob("GSM*_visium_*_spaceranger.tar.gz"))
    if not archives:
        raise FileNotFoundError(f"No GSE269622 Visium archives found in {VISIUM_DIR}")

    summaries = []
    contrast_rows = []
    spot_rows = []
    by_time_group: dict[tuple[str, str], pd.Series] = {}

    for archive in archives:
        timepoint, adata, positions = load_visium_sample(archive, extracted_root)
        transport = dct_transport_score(adata)
        dct_marker = score_genes(adata, NICHE_MARKERS["dct_enriched"])
        dct_high = dct_marker >= dct_marker.quantile(0.75)
        adjacent = dct_adjacent_mask(adata, dct_marker, positions)
        all_mask = pd.Series(True, index=adata.obs_names)
        masks = {
            "all_spots": all_mask,
            "dct_marker_top_quartile": dct_high,
            "dct_adjacent_spots": adjacent,
            "non_dct_adjacent_spots": ~adjacent,
        }
        for niche, genes in NICHE_MARKERS.items():
            score = score_genes(adata, genes)
            masks[f"{niche}_top_quartile"] = score >= score.quantile(0.75)

        for label, mask in masks.items():
            mask = mask.reindex(adata.obs_names).fillna(False).astype(bool)
            summaries.append(summarize(transport, mask, label, timepoint))
            by_time_group[(timepoint, label)] = transport.loc[mask.index[mask]].dropna()

        spot_rows.append(
            pd.DataFrame(
                {
                    "timepoint": timepoint,
                    "barcode": adata.obs_names.astype(str),
                    "dct_transport_score": transport.to_numpy(dtype=float),
                    "dct_marker_score": dct_marker.reindex(adata.obs_names).to_numpy(dtype=float),
                    "dct_marker_top_quartile": dct_high.reindex(adata.obs_names).fillna(False).to_numpy(dtype=bool),
                    "dct_adjacent_spot": adjacent.reindex(adata.obs_names).fillna(False).to_numpy(dtype=bool),
                }
            )
        )

    summary = pd.DataFrame(summaries)
    summary["timepoint_order"] = summary["timepoint"].map(TIME_ORDER)
    summary = summary.sort_values(["group", "timepoint_order"])
    summary.to_csv(out_dir / "visium_dct_transport_by_niche.tsv", sep="\t", index=False)
    pd.concat(spot_rows, ignore_index=True).to_csv(
        out_dir / "visium_dct_transport_spot_scores.tsv.gz", sep="\t", index=False, compression="gzip"
    )

    for group in summary["group"].unique():
        sham = by_time_group.get(("sham", group), pd.Series(dtype=float))
        for timepoint in ["hour4", "hour12", "day2", "day14", "week6"]:
            vals = by_time_group.get((timepoint, group), pd.Series(dtype=float))
            if len(vals) < 3 or len(sham) < 3:
                continue
            test = stats.ttest_ind(vals, sham, equal_var=False, nan_policy="omit")
            contrast_rows.append(
                {
                    "timepoint": timepoint,
                    "group": group,
                    "comparison": f"{timepoint}_minus_sham_same_group",
                    "mean_difference_vs_sham": float(vals.mean() - sham.mean()),
                    "t_stat": float(test.statistic),
                    "p_value": float(test.pvalue),
                    "n_timepoint_spots": int(len(vals)),
                    "n_sham_spots": int(len(sham)),
                    "statistical_boundary": "spot-level descriptive test; no animal-level spatial replication",
                }
            )

    for timepoint in TIME_ORDER:
        adj = by_time_group.get((timepoint, "dct_adjacent_spots"), pd.Series(dtype=float))
        non = by_time_group.get((timepoint, "non_dct_adjacent_spots"), pd.Series(dtype=float))
        if len(adj) < 3 or len(non) < 3:
            continue
        test = stats.ttest_ind(adj, non, equal_var=False, nan_policy="omit")
        contrast_rows.append(
            {
                "timepoint": timepoint,
                "group": "dct_adjacent_vs_nonadjacent",
                "comparison": "dct_adjacent_minus_nonadjacent_same_timepoint",
                "mean_difference_vs_sham": float(adj.mean() - non.mean()),
                "t_stat": float(test.statistic),
                "p_value": float(test.pvalue),
                "n_timepoint_spots": int(len(adj)),
                "n_sham_spots": int(len(non)),
                "statistical_boundary": "spot-level descriptive test; no animal-level spatial replication",
            }
        )

    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(out_dir / "visium_dct_transport_timepoint_summary.tsv", sep="\t", index=False)

    late = contrasts[
        contrasts["group"].eq("dct_adjacent_spots")
        & contrasts["comparison"].str.endswith("minus_sham_same_group")
        & contrasts["timepoint"].isin(["day14", "week6"])
    ]
    late_negative = bool((late["mean_difference_vs_sham"] < 0).any()) if not late.empty else False
    verdict = {
        "status": "complete",
        "late_iri_dct_adjacent_transport_suppression_present": late_negative,
        "late_iri_dct_adjacent_records": late.to_dict(orient="records"),
        "boundary": (
            "GSE269622 is an external IRI spatial reference. DCT transport scores are spot-level "
            "descriptive summaries, not replicated spatial validation of spaceflight kidney tissue."
        ),
    }
    (out_dir / "visium_dct_transport_verdict.json").write_text(json.dumps(verdict, indent=2))
    return verdict


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-root", type=Path, default=RUN_ROOT)
    args = ap.parse_args()
    verdict = run(args.run_root)
    print(f"[spatial-dct] complete: {args.run_root / 'spatial_reference'}")


if __name__ == "__main__":
    main()
