#!/usr/bin/env python3
"""Stage 0A/0C: protocol and transcript-integrity inventory for the OSDR kidney corpus.

Purpose
-------
Before any biological axis is tested, establish for each cohort:

* how the animal was euthanised and where (on-orbit vs post-return);
* whether tissue was dissected immediately or from an intact frozen carcass;
* which library-preparation chemistry was used (polyA vs ribodepletion/total);
* GeneLab processing pipeline version;
* per-sample RNA integrity and, critically, gene-body coverage.

Motivation is Lai Polo et al., iScience 2020 (doi:10.1016/j.isci.2020.101733),
which showed that preservation method can exceed the flight effect in magnitude,
that the distortion is amplified by polyA selection, and that RIN does **not**
predict it -- 5' gene-body coverage does. Choi et al. 2016 (PLOS One
doi:10.1371/journal.pone.0167391) separately found kidney RNA *quality* to be
minimally affected by carcass freezing to 6-7 months, so the open question is
coverage-level distortion, not RIN.

GeneLab's consensus pipeline emits RSeQC gene-body coverage into the published
``*_qc_metrics_*.csv``, so the key metric is downloadable rather than requiring
realignment.

Outputs
-------
``cohort_protocol_inventory.tsv``   one row per cohort
``sample_integrity_metrics.tsv``    one row per sample where QC is available
``eligibility_summary.md``          strata and the frozen eligibility verdict

Eligibility rule (frozen before inspection)
-------------------------------------------
A cohort is *ineligible* if flight and control samples differ systematically in
collection, preservation, or library preparation. Cohorts are then grouped into
preservation x library strata. **If no stratum contains at least three cohorts,
the confirmatory four-axis study does not proceed.**

Usage
-----
    python3 scripts/stage0/protocol_inventory.py
"""

from __future__ import annotations

import argparse
import csv
import glob
import io
import json
import os
import pathlib
import re
import sys
import zipfile

REPO = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "data/results/run_20260729_stage0_protocol_inventory"

# Frozen eligibility threshold: minimum cohorts sharing a stratum.
MIN_COHORTS_PER_STRATUM = 3

COHORTS = {
    "OSD-102": dict(isa=REPO / "data/external/osdr/OSD-102/OSD-102_metadata_OSD-102-ISA.zip",
                    qc=REPO / "data/external/osdr/OSD-102/GLDS-102_rna_seq_qc_metrics_GLbulkRNAseq.csv"),
    "OSD-163": dict(isa=REPO / "data/external/osdr/OSD-163/OSD-163_metadata_GLDS-163-ISA.zip",
                    qc=REPO / "data/external/osdr/OSD-163/GLDS-163_rna_seq_qc_metrics_GLbulkRNAseq.csv"),
    "OSD-253": dict(isa=REPO / "data/external/osdr/OSD-253/OSD-253_metadata_GLDS-253-ISA.zip",
                    qc=REPO / "data/external/osdr/OSD-253/GLDS-253_rna_seq_qc_metrics_GLbulkRNAseq.csv"),
    "OSD-462": dict(isa=REPO / "data/external/osdr/OSD-462/metadata",
                    qc=None),
    "OSD-513": dict(isa=REPO / "data/external/osdr/OSD-513/OSD-513_metadata_OSD-513-ISA.zip",
                    qc=REPO / "data/external/osdr/OSD-513/GLDS-513_rna_seq_qc_metrics_GLbulkRNAseq.csv"),
    "OSD-771": dict(isa=REPO / "data/raw/metadata",
                    qc=None),
}

# Phrases that identify collection/preservation mode in ISA protocol text.
CARCASS_PATTERNS = [
    r"carcass", r"frozen carcass", r"wrapped in aluminum", r"cryochiller",
    r"dissect\w*\s+(?:after|following)\s+return", r"biospecimen sharing",
]
IMMEDIATE_PATTERNS = [
    r"dissected immediately", r"immediately after euthanasia",
    r"immediately (?:frozen|preserved)", r"snap[- ]frozen (?:immediately|on)",
    r"rnalater",
]
ONORBIT_PATTERNS = [r"euthaniz\w+ on[- ]orbit", r"on[- ]orbit", r"euthanized in space",
                    r"aboard the (?:iss|international space station)"]
RETURN_PATTERNS = [r"returned (?:live|alive)", r"live animal return", r"after return to earth",
                   r"post[- ]return", r"returned to earth alive"]


def read_isa_texts(source: pathlib.Path) -> dict[str, str]:
    """Return {filename: text} for ISA files, from a zip or a directory."""
    out: dict[str, str] = {}
    if source is None or not source.exists():
        return out
    if source.is_dir():
        for p in sorted(source.glob("*.txt")):
            out[p.name] = p.read_text(errors="ignore")
    else:
        with zipfile.ZipFile(source) as z:
            for name in z.namelist():
                if name.endswith(".txt"):
                    out[os.path.basename(name)] = z.read(name).decode(
                        "utf-8", errors="ignore"
                    )
    return out


def find_patterns(text: str, patterns: list[str]) -> list[str]:
    low = text.lower()
    return [p for p in patterns if re.search(p, low)]


def classify(texts: dict[str, str]) -> dict:
    """Classify preservation and euthanasia location from ISA free text."""
    blob = "\n".join(texts.values())
    carcass = find_patterns(blob, CARCASS_PATTERNS)
    immediate = find_patterns(blob, IMMEDIATE_PATTERNS)
    onorbit = find_patterns(blob, ONORBIT_PATTERNS)
    ret = find_patterns(blob, RETURN_PATTERNS)

    if carcass and not immediate:
        preservation = "carcass"
    elif immediate and not carcass:
        preservation = "immediate"
    elif carcass and immediate:
        preservation = "mixed_or_ambiguous"
    else:
        preservation = "undetermined"

    if onorbit and ret:
        euth = "both_arms_present"
    elif onorbit:
        euth = "on_orbit"
    elif ret:
        euth = "post_return"
    else:
        euth = "undetermined"

    return {
        "preservation_class": preservation,
        "euthanasia_location": euth,
        "carcass_evidence": ";".join(carcass) or "",
        "immediate_evidence": ";".join(immediate) or "",
    }


def library_from_isa(texts: dict[str, str]) -> str:
    blob = "\n".join(texts.values()).lower()
    hits = []
    if re.search(r"ribo[- ]?deplet|ribo[- ]?zero|total rna", blob):
        hits.append("ribodepletion_or_total")
    if re.search(r"poly\s*\(?a\)?[- ]?(select|enrich)|mrna[- ]?seq|polyadenylated", blob):
        hits.append("polyA")
    return "+".join(hits) if hits else "undetermined"


def qc_summary(path: pathlib.Path) -> tuple[list[dict], dict]:
    """Per-sample integrity rows and a cohort-level summary."""
    if path is None or not path.exists():
        return [], {}
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        return [], {}

    def num(r, k):
        try:
            return float(r.get(k, "") or "nan")
        except ValueError:
            return float("nan")

    per_sample = []
    for r in rows:
        c5 = num(r, "mean_genebody_cov_5_20")
        c3 = num(r, "mean_genebody_cov_80_95")
        ratio = c5 / c3 if c3 and c3 == c3 and c3 != 0 else float("nan")
        per_sample.append(
            {
                "cohort": r.get("osd_num", ""),
                "sample": r.get("sample", ""),
                "tissue": r.get("tissue", ""),
                "library_selection": r.get("library_selection", ""),
                "strandedness": r.get("strandedness", ""),
                "read_length": r.get("read_length", ""),
                "rin": r.get("rin", ""),
                "genebody_cov_5_20": c5,
                "genebody_cov_40_60": num(r, "mean_genebody_cov_40_60"),
                "genebody_cov_80_95": c3,
                "cov_5to3_ratio": round(ratio, 4) if ratio == ratio else "",
            }
        )

    def mean(vals):
        vals = [v for v in vals if isinstance(v, float) and v == v]
        return round(sum(vals) / len(vals), 4) if vals else ""

    libs = sorted({r["library_selection"] for r in per_sample if r["library_selection"]})
    rins = [float(r["rin"]) for r in per_sample if r["rin"] not in ("", None)
            and re.match(r"^[0-9.]+$", str(r["rin"]))]
    summary = {
        "n_samples_qc": len(per_sample),
        "library_selection_qc": "|".join(libs),
        "mean_cov_5to3_ratio": mean([r["cov_5to3_ratio"] for r in per_sample
                                     if isinstance(r["cov_5to3_ratio"], float)]),
        "mean_rin": round(sum(rins) / len(rins), 2) if rins else "",
        "min_rin": min(rins) if rins else "",
    }
    return per_sample, summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    cohort_rows, sample_rows = [], []
    for cohort, spec in COHORTS.items():
        texts = read_isa_texts(spec["isa"])
        cls = classify(texts)
        per_sample, qsum = qc_summary(spec["qc"])
        for r in per_sample:
            r["cohort"] = r["cohort"] or cohort
        sample_rows.extend(per_sample)

        row = {
            "cohort": cohort,
            "isa_files_found": len(texts),
            **cls,
            "library_from_isa": library_from_isa(texts),
            "qc_metrics_available": "yes" if per_sample else "no",
            **{k: qsum.get(k, "") for k in
               ("n_samples_qc", "library_selection_qc", "mean_rin", "min_rin",
                "mean_cov_5to3_ratio")},
        }
        cohort_rows.append(row)

    def write(name, rows):
        p = args.out / name
        if not rows:
            p.write_text("")
            return p
        keys = list(dict.fromkeys(k for r in rows for k in r))
        with p.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        return p

    p1 = write("cohort_protocol_inventory.tsv", cohort_rows)
    p2 = write("sample_integrity_metrics.tsv", sample_rows)

    # Stratify and apply the frozen eligibility rule.
    strata: dict[str, list[str]] = {}
    for r in cohort_rows:
        lib = r["library_selection_qc"] or r["library_from_isa"] or "undetermined"
        key = f"{r['preservation_class']} | {lib}"
        strata.setdefault(key, []).append(r["cohort"])

    largest = max((len(v) for v in strata.values()), default=0)
    verdict = ("PROCEED" if largest >= MIN_COHORTS_PER_STRATUM
               else "DO NOT PROCEED with the confirmatory four-axis study")

    lines = [
        "# Stage 0A/0C eligibility summary",
        "",
        f"Frozen rule: a stratum must contain at least {MIN_COHORTS_PER_STRATUM} "
        "cohorts for the confirmatory study to proceed.",
        "",
        "## Strata (preservation | library selection)",
        "",
    ]
    for k, v in sorted(strata.items(), key=lambda kv: -len(kv[1])):
        lines.append(f"- **{k}** — {len(v)} cohort(s): {', '.join(sorted(v))}")
    lines += [
        "",
        f"Largest stratum: **{largest}** cohort(s).",
        "",
        f"## Verdict: {verdict}",
        "",
        "Automated ISA text classification is a screen, not an adjudication. "
        "Any cohort marked `undetermined` or `mixed_or_ambiguous` must be "
        "resolved by reading the source protocol before the rule is applied.",
    ]
    p3 = args.out / "eligibility_summary.md"
    p3.write_text("\n".join(lines) + "\n")

    (args.out / "manifest.json").write_text(json.dumps({
        "purpose": "Stage 0A/0C protocol and transcript-integrity inventory.",
        "min_cohorts_per_stratum": MIN_COHORTS_PER_STRATUM,
        "cohorts": list(COHORTS),
        "outputs": [p.name for p in (p1, p2, p3)],
        "references": {
            "preservation_confounding": "10.1016/j.isci.2020.101733",
            "kidney_rna_quality_after_carcass_freezing": "10.1371/journal.pone.0167391",
        },
    }, indent=2) + "\n")

    print(f"wrote {p1}\nwrote {p2}\nwrote {p3}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
