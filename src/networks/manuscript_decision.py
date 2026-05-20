"""Pre-registered manuscript decision tree for the contrast-vector framework."""
from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Mapping


@dataclass
class ManuscriptDecision:
    branch: str
    implication: str
    headline: str
    evidence: dict[str, object]

    def as_dict(self) -> dict:
        return asdict(self)


def decide_manuscript_branch(evidence: Mapping[str, object]) -> ManuscriptDecision:
    """Map observed outcomes to the locked Phase 8 decision tree."""
    within_stable = bool(evidence.get("within_rrrm2_stability_passed", False))
    within_signal = bool(evidence.get("within_rrrm2_projection_significant", False))
    external_recurrence = bool(evidence.get("cross_osdr_recurrence", False))
    external_axis_signal = bool(evidence.get("external_axis_significant", False))
    module_only = bool(evidence.get("module_activity_only_positive", False))
    any_projection = bool(evidence.get("any_projection_signal", False))

    if within_stable and within_signal and external_recurrence:
        return ManuscriptDecision(
            branch="strong_network_aging_paper",
            implication="Strong network-aging paper.",
            headline="Spaceflight modifies kidney aging network architecture.",
            evidence=dict(evidence),
        )
    if external_axis_signal and not external_recurrence:
        return ManuscriptDecision(
            branch="external_aging_axis_alignment",
            implication="External aging-axis flight-network biology paper.",
            headline=(
                "Within-RRRM-2 aging-axis projection was gated off; RRRM-2 flight "
                "effects are interpreted against the independent TMS kidney aging axis."
            ),
            evidence=dict(evidence),
        )
    if not within_signal and external_recurrence:
        return ManuscriptDecision(
            branch="cross_osdr_flight_network_biology",
            implication="Cross-OSDR flight-network biology paper.",
            headline=(
                "RRRM-2 provides the factorial bridge while external cohorts provide "
                "the recurring flight-network direction."
            ),
            evidence=dict(evidence),
        )
    if module_only:
        return ManuscriptDecision(
            branch="modest_module_activity_paper",
            implication="Modest module-activity paper.",
            headline="Module activity shifts are present without a defensible projection signal.",
            evidence=dict(evidence),
        )
    return ManuscriptDecision(
        branch="negative_or_methods_constraint",
        implication="Do not publish as a standalone biology claim.",
        headline="The framework records a negative or methods-constraint result.",
        evidence=dict(evidence),
    )


def write_manuscript_decision(evidence: Mapping[str, object], outpath: str | Path) -> ManuscriptDecision:
    """Write the decision artifact as JSON and return the decision object."""
    decision = decide_manuscript_branch(evidence)
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w") as fh:
        json.dump(decision.as_dict(), fh, indent=2)
    return decision
