from pathlib import Path

import pandas as pd
import yaml

from src.clinical_axes.analysis import combined_score_design
from src.clinical_axes.data import load_primary_missions
from scripts.clinical_axes.run_compartment_context import (
    BARRIER_CORE_GENES,
    BARRIER_EXPANDED_GENES,
    add_disjoint_podocyte_variants,
)


REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "config/clinical_renal_axes_cross_mission.yaml"


def _load():
    config = yaml.safe_load(CONFIG.read_text())
    return config, load_primary_missions(config, REPO)


def test_frozen_primary_contrasts_resolve_to_expected_biological_animals():
    _, missions = _load()
    expected = {
        "OSD-102": (6, 6, {"all"}),
        "OSD-163": (6, 6, {"all"}),
        "OSD-253": (10, 9, {"day25", "day75"}),
        "OSD-462": (10, 10, {"all"}),
        "OSD-771": (10, 10, {"YNG", "OLD"}),
    }
    assert set(missions) == set(expected)
    for mission, data in missions.items():
        n_flight, n_control, blocks = expected[mission]
        assert int(data.metadata["condition"].eq("FLT").sum()) == n_flight
        assert int(data.metadata["condition"].eq("GC").sum()) == n_control
        assert set(data.metadata["block"]) == blocks
        assert data.metadata.index.is_unique
        assert list(data.expression.columns) == list(data.metadata.index)
        assert list(data.counts.columns) == list(data.metadata.index)


def test_every_frozen_axis_meets_declared_coverage_in_every_primary_mission():
    config, missions = _load()
    scores, design, coverage, _ = combined_score_design(
        missions,
        config["primary_family"],
        cpm_threshold=float(config["eligibility"]["cpm_threshold"]),
    )
    assert scores.shape == (sum(len(x.metadata) for x in missions.values()), 4)
    assert scores.notna().all().all()
    assert set(design["mission"]) == set(missions)
    assert coverage["coverage_pass"].all()
    assert (coverage["n_used"] >= coverage["minimum_required"]).all()


def test_disjoint_podocyte_variants_leave_source_set_unchanged():
    source_genes = {
        gene: 1
        for gene in (
            *BARRIER_EXPANDED_GENES,
            "Nid1",
            "Col4a3",
            "Col4a4",
            "Kirrel1",
            "Myh9",
            "Actn4",
            "Inf2",
            "Nck1",
        )
    }
    family = {
        "podocyte__high_specificity": {
            "role": "secondary_compartment_context",
            "subdomains": {
                "atlas_markers": {
                    "genes": source_genes.copy(),
                    "minimum_present": 8,
                }
            },
        }
    }
    audit = pd.DataFrame(
        [
            {
                "gene_set": "podocyte__high_specificity",
                "report_compartment": "podocyte",
                "tier": "high_specificity",
                "n_defined": len(source_genes),
                "definition_evaluable": True,
            }
        ]
    )

    augmented, augmented_audit = add_disjoint_podocyte_variants(family, audit)

    original = augmented["podocyte__high_specificity"]["subdomains"][
        "atlas_markers"
    ]["genes"]
    core_disjoint = augmented[
        "podocyte__high_specificity_disjoint_barrier_core"
    ]["subdomains"]["atlas_markers"]["genes"]
    expanded_disjoint = augmented[
        "podocyte__high_specificity_disjoint_barrier_expanded"
    ]["subdomains"]["atlas_markers"]["genes"]
    assert set(original) == set(source_genes)
    assert set(core_disjoint) == set(source_genes) - set(BARRIER_CORE_GENES)
    assert set(expanded_disjoint) == set(source_genes) - set(
        BARRIER_EXPANDED_GENES
    )
    assert len(augmented_audit) == 3
