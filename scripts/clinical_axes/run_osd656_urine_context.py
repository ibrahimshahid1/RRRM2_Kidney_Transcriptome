#!/usr/bin/env python3
"""Write descriptive OSD-656 urine context for frozen clinical tissue axes."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import yaml


REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from src.clinical_axes.osd656_context import build_osd656_context  # noqa: E402
from src.v11.human_concordance import parse_osd656_submitted  # noqa: E402


DEFAULT_CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"
DEFAULT_INPUT = (
    REPO
    / "data/external/human_spaceflight/OSD-656"
    / "LSDS-64_Multiplex_urine.immune.AlamarPanel_SUBMITTED.xlsx"
)
DEFAULT_OUTPUT = REPO / "data/results/run_20260811_clinical_renal_axes_cross_mission"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run(config_path: Path, input_path: Path, output_dir: Path) -> dict[str, object]:
    config = yaml.safe_load(config_path.read_text())
    long = parse_osd656_submitted(input_path)
    coverage, paired, summary = build_osd656_context(long, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "coverage": output_dir / "osd656_urine_axis_coverage.tsv",
        "paired_changes": output_dir / "osd656_urine_paired_changes.tsv",
        "summary": output_dir / "osd656_urine_context_summary.tsv",
        "manifest": output_dir / "osd656_urine_context_manifest.json",
    }
    coverage.to_csv(paths["coverage"], sep="\t", index=False)
    paired.to_csv(paths["paired_changes"], sep="\t", index=False)
    summary.to_csv(paths["summary"], sep="\t", index=False)

    measured = coverage[coverage["panel_status"].eq("measured_overlap")]
    barrier = coverage[coverage["axis"].eq("glomerular_barrier_identity_loss")]
    havcr1 = coverage[coverage["gene_symbol"].str.casefold().eq("havcr1")]
    manifest: dict[str, object] = {
        "analysis": "OSD-656 descriptive urine context for frozen clinical renal axes",
        "status": "descriptive_context_only",
        "source": "OSD-656 Inspiration4 urine Alamar panel",
        "input": str(input_path),
        "input_sha256": sha256(input_path),
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "n_frozen_primary_genes": int(len(coverage)),
        "n_measured_overlaps": int(len(measured)),
        "measured_overlaps": measured["gene_symbol"].tolist(),
        "n_absent_frozen_genes": int(coverage["panel_status"].eq("absent_from_panel").sum()),
        "barrier_markers_measured": int(barrier["panel_status"].eq("measured_overlap").sum()),
        "havcr1_measured": bool(havcr1["panel_status"].eq("measured_overlap").any()),
        "n_subjects": int(paired["subject"].nunique()),
        "n_r_plus_1_subjects": int(
            paired.loc[paired["timepoint"].eq("R+1"), "subject"].nunique()
        ),
        "interpretation_boundary": (
            "Post-flight human urine context only. These measurements are not inflight "
            "kidney tissue, do not validate the mouse axes, and receive no inferential "
            "weight in the cross-mission analysis."
        ),
        "design_caveats": [
            "four crew members and no unflown control group",
            "no inflight urine; R+1 is the earliest recovery collection",
            "crew member C001 lacks L-3 and R+1 samples",
            "NPQ is a relative assay readout, not a clinical concentration or urine-creatinine-corrected result",
            "urinary protein and kidney-tissue RNA are different biological quantities",
            "the inflammation panel is not kidney-specific",
        ],
        "outputs": {name: str(path) for name, path in paths.items() if name != "manifest"},
    }
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = run(args.config, args.input, args.output_dir)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
