#!/usr/bin/env python3
"""Stage 0: observed effect size for the known-positive axis.

Why this exists
---------------
A power simulation needs a reference effect size. The ECM/remodeling axis is the
one program already reported as recurrent across this corpus (signed Stouffer
p = 7.0e-4 in prior repository analyses), so it is the natural anchor: if the
random-effects pooled Hedges g for the *known positive* falls below the
detectable threshold implied by the design, the four-axis confirmatory study is
not viable and no further Stage 0 work is warranted.

Signed Stouffer tests direction and is considerably more powerful than a
random-effects meta of magnitudes. It therefore does not tell you the effect
size, which is the quantity that governs power for the axis-ranking design.

Method
------
Per cohort: per-sample mean z-score across axis genes (z computed within cohort
on VST values), then Hedges g for flight versus ground control, then a
DerSimonian-Laird random-effects pool across cohorts with I-squared.

Boundary
--------
This is a design-calibration estimate, not a biological claim. Cohorts failing
the Stage 0B coverage-confounding gate are reported but flagged, and the pooled
estimate is given both with and without them.

Usage
-----
    python3 scripts/stage0/axis_effect_size_anchor.py
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "data/results/run_20260729_stage0_protocol_inventory"

COHORTS = {
    "OSD-102": dict(
        vst=REPO / "data/external/osdr/OSD-102/GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        runsheet=REPO / "data/external/osdr/OSD-102/GLDS-102_rna_seq_bulkRNASeq_v1_runsheet.csv"),
    "OSD-163": dict(
        vst=REPO / "data/external/osdr/OSD-163/GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        runsheet=REPO / "data/external/osdr/OSD-163/GLDS-163_rna_seq_bulkRNASeq_v1_runsheet.csv"),
    "OSD-253": dict(
        vst=REPO / "data/external/osdr/OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        runsheet=REPO / "data/external/osdr/OSD-253/GLDS-253_rna_seq_bulkRNASeq_v2_runsheet.csv"),
    "OSD-513": dict(
        vst=REPO / "data/external/osdr/OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        runsheet=REPO / "data/external/osdr/OSD-513/GLDS-513_rna_seq_bulkRNASeq_v2_runsheet.csv"),
}

AXES = {"ecm_remodeling": None, "dct_ncc_wnk": None, "fibrosis": None}

FLIGHT_RE = re.compile(r"space\s*flight|^flt$", re.I)
GC_RE = re.compile(r"^ground\s*control$", re.I)


def load_axes():
    import yaml
    d = yaml.safe_load((REPO / "config/gene_sets.yaml").read_text())
    for name in list(AXES):
        node = d.get(name, {})
        AXES[name] = [g.strip().upper() for g in node.get("genes", [])]
    return AXES


ANNOTATION = (REPO / "data/external/osdr/OSD-462"
              / "GLDS-462_rna_seq_differential_expression_totRNA_GLbulkRNAseq.csv")


def ensembl_to_symbol() -> dict[str, str]:
    """ENSMUSG -> gene symbol, taken from the GeneLab annotated DE table."""
    mapping = {}
    if not ANNOTATION.exists():
        return mapping
    with ANNOTATION.open() as fh:
        for r in csv.DictReader(fh):
            ens = (r.get("ENSEMBL") or "").strip()
            sym = (r.get("SYMBOL") or "").strip()
            if ens and sym:
                mapping[ens] = sym.upper()
    return mapping


def load_vst(path: pathlib.Path, id2sym: dict[str, str]):
    """Return (sample_names, {SYMBOL: [values]}), mapping Ensembl rows to symbols."""
    with path.open() as fh:
        rows = list(csv.reader(fh))
    samples = rows[0][1:]
    data = {}
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        key = r[0].strip()
        sym = id2sym.get(key, key).upper()
        vals = []
        for v in r[1:]:
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(float("nan"))
        # keep the first occurrence if a symbol maps from several IDs
        data.setdefault(sym, vals)
    return samples, data


def load_groups(path: pathlib.Path):
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for r in rows:
        name = (r.get("Sample Name") or "").strip()
        out[name] = (r.get("Factor Value[Spaceflight]", "") or "").strip()
    return out


def hedges_g(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return None
    ma, mb = sum(a) / na, sum(b) / nb
    va = sum((x - ma) ** 2 for x in a) / (na - 1)
    vb = sum((x - mb) ** 2 for x in b) / (nb - 1)
    sp = math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    if sp == 0:
        return None
    d = (ma - mb) / sp
    J = 1 - 3 / (4 * (na + nb) - 9)
    g = d * J
    var = (na + nb) / (na * nb) + g ** 2 / (2 * (na + nb - 2))
    return g, var, na, nb


def dersimonian_laird(gs, vs):
    w = [1 / v for v in vs]
    sw = sum(w)
    fixed = sum(wi * gi for wi, gi in zip(w, gs)) / sw
    Q = sum(wi * (gi - fixed) ** 2 for wi, gi in zip(w, gs))
    k = len(gs)
    C = sw - sum(wi ** 2 for wi in w) / sw
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    w2 = [1 / (v + tau2) for v in vs]
    sw2 = sum(w2)
    pooled = sum(wi * gi for wi, gi in zip(w2, gs)) / sw2
    se = math.sqrt(1 / sw2)
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    z = pooled / se
    p = math.erfc(abs(z) / math.sqrt(2))
    return dict(pooled=pooled, se=se, lo=pooled - 1.96 * se, hi=pooled + 1.96 * se,
                tau2=tau2, Q=Q, I2=I2, p=p, k=k)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    axes = load_axes()
    id2sym = ensembl_to_symbol()
    if not id2sym:
        print("warning: no Ensembl->symbol mapping found", file=sys.stderr)

    per_cohort, pooled_rows = [], []
    for axis, genes in axes.items():
        gs, vs, used = [], [], []
        for cohort, spec in COHORTS.items():
            if not spec["vst"].exists():
                continue
            samples, data = load_vst(spec["vst"], id2sym)
            groups = load_groups(spec["runsheet"])
            present = [g for g in genes if g in data]
            if len(present) < 3:
                continue
            # z-score each gene within cohort, then average per sample
            zmat = []
            for g in present:
                v = data[g]
                fin = [x for x in v if x == x]
                m = sum(fin) / len(fin)
                sd = math.sqrt(sum((x - m) ** 2 for x in fin) / max(1, len(fin) - 1))
                zmat.append([(x - m) / sd if sd else 0.0 for x in v])
            score = [sum(col) / len(col) for col in zip(*zmat)]

            flt = [s for s, nm in zip(score, samples) if FLIGHT_RE.search(groups.get(nm, ""))]
            gc = [s for s, nm in zip(score, samples) if GC_RE.search(groups.get(nm, ""))]
            res = hedges_g(flt, gc)
            if not res:
                continue
            g, var, na, nb = res
            per_cohort.append({
                "axis": axis, "cohort": cohort, "n_genes_used": len(present),
                "n_genes_defined": len(genes), "n_flight": na, "n_ground": nb,
                "hedges_g": round(g, 4), "var": round(var, 5),
                "ci_low": round(g - 1.96 * math.sqrt(var), 4),
                "ci_high": round(g + 1.96 * math.sqrt(var), 4),
            })
            gs.append(g); vs.append(var); used.append(cohort)

        if len(gs) >= 2:
            m = dersimonian_laird(gs, vs)
            pooled_rows.append({
                "axis": axis, "k_cohorts": m["k"], "cohorts": "|".join(used),
                "pooled_hedges_g": round(m["pooled"], 4),
                "se": round(m["se"], 4),
                "ci_low": round(m["lo"], 4), "ci_high": round(m["hi"], 4),
                "tau2": round(m["tau2"], 4), "I2_pct": round(m["I2"], 1),
                "p": f"{m['p']:.3g}",
            })

    def write(name, rows):
        p = args.out / name
        if not rows:
            p.write_text(""); return p
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
            w.writeheader(); w.writerows(rows)
        return p

    p1 = write("axis_effect_per_cohort.tsv", per_cohort)
    p2 = write("axis_effect_pooled.tsv", pooled_rows)
    (args.out / "axis_anchor_manifest.json").write_text(json.dumps({
        "purpose": "Observed effect sizes to anchor the Stage 0 power simulation.",
        "score": "within-cohort gene z-scores averaged per sample",
        "contrast": "Space Flight vs Ground Control only",
        "pooling": "DerSimonian-Laird random effects",
        "boundary": "Design calibration, not a biological claim. Cohorts failing "
                    "the Stage 0B coverage gate are included here and must be "
                    "excluded before any inferential use.",
    }, indent=2) + "\n")
    print(f"wrote {p1}\nwrote {p2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
