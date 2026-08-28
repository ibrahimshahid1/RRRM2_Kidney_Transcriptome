#!/usr/bin/env python3
"""Stage 0B: does transcript coverage differ between flight and control?

The fatal-confounding gate. Lai Polo et al. (iScience 2020) showed that
carcass preservation degrades 5' gene-body coverage while leaving RIN high, and
that this distortion can exceed the biological flight effect and completely
change the differential-expression landscape.

The question that decides cohort eligibility is therefore not "is coverage
degraded" -- it is degraded everywhere -- but **whether degradation differs
between the groups being contrasted**. If flight and control samples within a
cohort have systematically different 5'/3' coverage, that cohort's flight
estimate is confounded at the measurement layer and no downstream modelling
repairs it.

Metric: ratio of RSeQC mean gene-body coverage in the 5-20% bin to the 80-95%
bin, taken from GeneLab's published ``*_qc_metrics_*.csv``. A ratio near 1 is
uniform; below 1 indicates 5' loss.

Test: Welch t on the ratio between flight and ground-control samples, plus
Hedges g, per cohort. Reported alongside the same contrast for RIN so the two
can be compared -- RIN is expected to look fine.

Usage
-----
    python3 scripts/stage0/coverage_confounding_gate.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re

REPO = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "data/results/run_20260729_stage0_protocol_inventory"
INVENTORY = DEFAULT_OUT / "sample_integrity_metrics.tsv"

RUNSHEETS = {
    "OSD-102": REPO / "data/external/osdr/OSD-102/GLDS-102_rna_seq_bulkRNASeq_v1_runsheet.csv",
    "OSD-163": REPO / "data/external/osdr/OSD-163/GLDS-163_rna_seq_bulkRNASeq_v1_runsheet.csv",
    "OSD-253": REPO / "data/external/osdr/OSD-253/GLDS-253_rna_seq_bulkRNASeq_v2_runsheet.csv",
    "OSD-513": REPO / "data/external/osdr/OSD-513/GLDS-513_rna_seq_bulkRNASeq_v2_runsheet.csv",
}

FLIGHT_RE = re.compile(r"space\s*flight|^flt$|flight", re.I)
CONTROL_RE = re.compile(r"ground\s*control|^gc$|vivarium|^viv$|basal|^bsl$|control", re.I)


def welch(a: list[float], b: list[float]):
    """Welch t, df, two-sided p (normal approx), and Hedges g."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return None
    ma, mb = sum(a) / na, sum(b) / nb
    va = sum((x - ma) ** 2 for x in a) / (na - 1)
    vb = sum((x - mb) ** 2 for x in b) / (nb - 1)
    se = math.sqrt(va / na + vb / nb)
    if se == 0:
        return None
    t = (ma - mb) / se
    df = (va / na + vb / nb) ** 2 / (
        (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    )
    p = math.erfc(abs(t) / math.sqrt(2))  # normal approximation
    sp = math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    d = (ma - mb) / sp if sp else float("nan")
    J = 1 - 3 / (4 * (na + nb) - 9)
    return dict(mean_flight=ma, mean_control=mb, t=t, df=df, p=p, hedges_g=d * J)


def load_conditions(path: pathlib.Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for r in rows:
        name = r.get("Sample Name") or r.get("Original Sample Name") or ""
        val = r.get("Factor Value[Spaceflight]", "") or ""
        out[name.strip()] = val.strip()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    with INVENTORY.open() as fh:
        samples = list(csv.DictReader(fh, delimiter="\t"))

    results = []
    for cohort, rs in RUNSHEETS.items():
        cond = load_conditions(rs)
        groups = {"flight": {"cov": [], "rin": []},
                  "control": {"cov": [], "rin": []}}
        unmatched = 0
        levels = set()
        for s in samples:
            if s["cohort"] != cohort:
                continue
            c = cond.get(s["sample"].strip(), "")
            if not c:
                unmatched += 1
                continue
            levels.add(c)
            if FLIGHT_RE.search(c):
                key = "flight"
            elif CONTROL_RE.search(c):
                key = "control"
            else:
                continue
            try:
                groups[key]["cov"].append(float(s["cov_5to3_ratio"]))
            except (ValueError, TypeError):
                pass
            try:
                groups[key]["rin"].append(float(s["rin"]))
            except (ValueError, TypeError):
                pass

        row = {
            "cohort": cohort,
            "n_flight": len(groups["flight"]["cov"]),
            "n_control": len(groups["control"]["cov"]),
            "unmatched_samples": unmatched,
            "factor_levels": "|".join(sorted(levels)),
        }
        cov = welch(groups["flight"]["cov"], groups["control"]["cov"])
        rin = welch(groups["flight"]["rin"], groups["control"]["rin"])
        if cov:
            row.update({
                "cov_flight_mean": round(cov["mean_flight"], 4),
                "cov_control_mean": round(cov["mean_control"], 4),
                "cov_hedges_g": round(cov["hedges_g"], 3),
                "cov_p": f"{cov['p']:.3g}",
            })
        if rin:
            row.update({
                "rin_flight_mean": round(rin["mean_flight"], 3),
                "rin_control_mean": round(rin["mean_control"], 3),
                "rin_hedges_g": round(rin["hedges_g"], 3),
                "rin_p": f"{rin['p']:.3g}",
            })
        if cov:
            row["gate"] = ("CONFOUNDED" if cov["p"] < 0.05 and abs(cov["hedges_g"]) >= 0.5
                           else "pass")
        else:
            row["gate"] = "not_evaluable"
        results.append(row)

    path = args.out / "coverage_confounding_gate.tsv"
    keys = list(dict.fromkeys(k for r in results for k in r))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader()
        w.writerows(results)

    (args.out / "coverage_gate_manifest.json").write_text(json.dumps({
        "metric": "RSeQC gene-body coverage 5-20% / 80-95% ratio",
        "gate_rule": "CONFOUNDED if Welch p<0.05 and |Hedges g|>=0.5 between "
                     "flight and control within cohort",
        "rationale_doi": "10.1016/j.isci.2020.101733",
        "note": "p uses a normal approximation to the Welch t; treat borderline "
                "values as indicative and recompute exactly before locking.",
    }, indent=2) + "\n")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
