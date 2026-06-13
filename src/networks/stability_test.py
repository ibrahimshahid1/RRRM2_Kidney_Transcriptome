"""Guardrail A - bootstrap stability of the control aging vector A^a_GC."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence
import json

import numpy as np
import pandas as pd

from .contrast_vectors import (
    EPS,
    cosine,
    stratified_bootstrap_indices,
)


# Default gates per §2.1 / §4.1
DEFAULT_GATES: dict[str, dict[str, float]] = {
    "edge": {"median": 0.30, "lower": 0.10},
    "gene": {"median": 0.40, "lower": 0.20},
    "lioness_node": {"median": 0.40, "lower": 0.20},
    "pathway": {"median": 0.60, "lower": 0.40},
    "module": {"median": 0.60, "lower": 0.40},
}


def gates_from_config(config: dict) -> dict[str, dict[str, float]]:
    """Extract stability gates from config/contrast_vector_framework.yaml.

    The YAML uses explicit keys (``stability_median`` / ``stability_lower``)
    while the in-memory gate application uses the shorter ``median`` /
    ``lower`` names.
    """
    resolutions = (config or {}).get("resolutions", {})
    if not isinstance(resolutions, dict):
        return DEFAULT_GATES
    gates: dict[str, dict[str, float]] = {}
    for resolution, values in resolutions.items():
        if not isinstance(values, dict):
            continue
        if "stability_median" not in values or "stability_lower" not in values:
            continue
        gates[str(resolution)] = {
            "median": float(values["stability_median"]),
            "lower": float(values["stability_lower"]),
        }
    return gates or DEFAULT_GATES


def resolution_role_from_config(config: dict, resolution: str) -> str:
    """Return the configured role for a resolution."""
    values = ((config or {}).get("resolutions", {}) or {}).get(resolution, {})
    if isinstance(values, dict):
        return str(values.get("role", "primary"))
    return "primary"


def on_fail_from_config(config: dict, resolution: str) -> str:
    """Return the configured failure action for a resolution."""
    values = ((config or {}).get("resolutions", {}) or {}).get(resolution, {})
    if isinstance(values, dict):
        return str(values.get("on_fail", "demote_to_exploratory"))
    return "demote_to_exploratory"


@dataclass
class StabilityReport:
    """Per-(arm, resolution) stability report row."""

    arm: str
    resolution: str
    n_bootstrap: int
    full_vector_norm: float
    cosine_median: float
    cosine_lower: float
    cosine_upper: float

    def as_dict(self) -> dict:
        return {
            "arm": self.arm,
            "resolution": self.resolution,
            "n_bootstrap": self.n_bootstrap,
            "full_vector_norm": self.full_vector_norm,
            "cosine_median": self.cosine_median,
            "cosine_lower": self.cosine_lower,
            "cosine_upper": self.cosine_upper,
        }


@dataclass
class StabilityDecision:
    """Pass/fail decision derived from a StabilityReport + gates config."""

    arm: str
    resolution: str
    median_required: float
    lower_required: float
    cosine_median: float
    cosine_lower: float
    passed: bool
    role: str
    on_fail_action: str

    def as_dict(self) -> dict:
        return {
            "arm": self.arm,
            "resolution": self.resolution,
            "median_required": self.median_required,
            "lower_required": self.lower_required,
            "cosine_median": self.cosine_median,
            "cosine_lower": self.cosine_lower,
            "passed": bool(self.passed),
            "role": self.role,
            "on_fail_action": self.on_fail_action,
        }


def estimate_agc_stability(
    agc_builder: Callable[[np.ndarray | None], np.ndarray],
    strata: Sequence,
    arm: str,
    resolution: str,
    n_iterations: int = 2000,
    rng: np.random.Generator | None = None,
) -> StabilityReport:
    """Compute bootstrap-vs-full-sample angular stability of A^a_GC.

    ``agc_builder(indices)`` should rebuild A^a_GC on the resampled rows.
    Calling it with ``None`` returns the full-sample (no-resample) vector.
    """
    if rng is None:
        rng = np.random.default_rng()

    full_vec = np.asarray(agc_builder(None), dtype=float).ravel()
    if full_vec.size == 0:
        raise ValueError(f"agc_builder returned empty vector at arm={arm}, resolution={resolution}")

    full_norm = float(np.linalg.norm(full_vec))
    if full_norm <= EPS:
        return StabilityReport(
            arm=arm, resolution=resolution, n_bootstrap=int(n_iterations),
            full_vector_norm=full_norm,
            cosine_median=float("nan"),
            cosine_lower=float("nan"),
            cosine_upper=float("nan"),
        )

    cosines = np.empty(int(n_iterations), dtype=float)
    for b in range(int(n_iterations)):
        idx = stratified_bootstrap_indices(strata, rng)
        try:
            vb = np.asarray(agc_builder(idx), dtype=float).ravel()
            cosines[b] = cosine(full_vec, vb)
        except ValueError:
            cosines[b] = np.nan

    finite = cosines[np.isfinite(cosines)]
    if finite.size == 0:
        med, lo, hi = (float("nan"),) * 3
    else:
        med = float(np.median(finite))
        lo = float(np.percentile(finite, 2.5))
        hi = float(np.percentile(finite, 97.5))

    return StabilityReport(
        arm=arm,
        resolution=resolution,
        n_bootstrap=int(n_iterations),
        full_vector_norm=full_norm,
        cosine_median=med,
        cosine_lower=lo,
        cosine_upper=hi,
    )


def apply_stability_gate(
    report: StabilityReport,
    gates: dict[str, dict[str, float]] | None = None,
    role: str = "primary",
    on_fail_action: str = "demote_to_exploratory",
) -> StabilityDecision:
    """Apply the pre-registered gate for a (arm, resolution) report."""
    gates = gates or DEFAULT_GATES
    if report.resolution not in gates:
        raise KeyError(
            f"No gate configured for resolution={report.resolution}. "
            f"Configure one in config/contrast_vector_framework.yaml."
        )
    gate = gates[report.resolution]
    median_required = float(gate["median"])
    lower_required = float(gate["lower"])
    passed = (
        np.isfinite(report.cosine_median)
        and np.isfinite(report.cosine_lower)
        and report.cosine_median >= median_required
        and report.cosine_lower >= lower_required
    )
    return StabilityDecision(
        arm=report.arm,
        resolution=report.resolution,
        median_required=median_required,
        lower_required=lower_required,
        cosine_median=report.cosine_median,
        cosine_lower=report.cosine_lower,
        passed=bool(passed),
        role=role,
        on_fail_action=on_fail_action,
    )


def write_stability_artifacts(
    reports: list[StabilityReport],
    decisions: list[StabilityDecision],
    outdir: Path,
    report_filename: str = "agc_stability_report.tsv",
    decision_filename: str = "agc_stability_decision.json",
    manifest: dict | None = None,
) -> dict[str, Path]:
    """Persist the report + decision artifacts. Returns dict of written paths."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    report_path = outdir / report_filename
    pd.DataFrame([r.as_dict() for r in reports]).to_csv(report_path, sep="\t", index=False)

    decision_path = outdir / decision_filename
    decision_payload: dict[str, object] = {
        "decisions": [d.as_dict() for d in decisions],
        "any_module_failed": any(
            (not d.passed) and d.resolution == "module" for d in decisions
        ),
        "any_pathway_failed": any(
            (not d.passed) and d.resolution == "pathway" for d in decisions
        ),
        "fallback_to_external_axis_only": any(
            (not d.passed)
            and d.resolution in ("module", "pathway")
            and d.on_fail_action in {"fallback_external_axis_only", "external_axis_only"}
            for d in decisions
        ),
    }
    if manifest is not None:
        decision_payload["manifest"] = manifest
    with open(decision_path, "w") as fh:
        json.dump(decision_payload, fh, indent=2)
    return {"report": report_path, "decision": decision_path}


def load_stability_decision(decision_path: Path) -> dict:
    """Re-read the stability decision JSON. Phase 4 must call this before running."""
    decision_path = Path(decision_path)
    if not decision_path.exists():
        raise FileNotFoundError(
            f"Stability decision artifact not found at {decision_path}. "
            f"Phase 3 (Guardrail A) must run before downstream phases."
        )
    with open(decision_path) as fh:
        return json.load(fh)


def resolutions_passing(decision: dict, arm: str | None = None) -> list[str]:
    """Return list of resolutions whose stability gate passed (optionally per arm)."""
    out: list[str] = []
    for row in decision.get("decisions", []):
        if not row.get("passed"):
            continue
        if arm is not None and row.get("arm") != arm:
            continue
        out.append(row["resolution"])
    return sorted(set(out))


def assert_within_projection_allowed(decision: dict, bypass_stability: bool = False) -> None:
    """Enforce Phase 4's stability gate before within-RRRM-2 projection."""
    if bypass_stability:
        return
    if decision.get("fallback_to_external_axis_only"):
        raise RuntimeError(
            "Guardrail A failed at module/pathway level; within-RRRM-2 projection "
            "must not run without --bypass-stability. Use the external-axis path."
        )
