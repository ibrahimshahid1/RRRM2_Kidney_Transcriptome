import pytest
import pandas as pd

from src.validation.external_replication import (
    compare_directional_replication,
    validate_hypothesis_registry,
    validate_study_scope,
)


def test_external_protocol_guards():
    validate_study_scope("OSD-102", "lar_young")
    validate_study_scope("OSD-513", "sex_robustness")
    validate_study_scope("OSD-163", "cross_strain_context")
    validate_study_scope("OSD-253", "rr7_context")
    with pytest.raises(ValueError, match="OSD-568"):
        validate_study_scope("OSD-568", "lar_young")
    with pytest.raises(ValueError, match="lar_young"):
        validate_study_scope("OSD-102", "sex_robustness")


def test_external_registry_rejects_placeholders_and_locks_thresholds():
    bad = pd.DataFrame([{
        "analysis": "lar_young",
        "study": "OSD-102",
        "hypothesis_type": "pathway",
        "feature": "pre_registered_lar_young_pathways",
        "expected_direction": "same",
        "discovery_effect": 1.0,
        "q_threshold": 0.10,
        "claim_boundary": "primary",
    }])
    with pytest.raises(ValueError, match="placeholder"):
        validate_hypothesis_registry(bad)

    good = bad.copy()
    good["feature"] = "PPAR_signaling"
    registry = validate_hypothesis_registry(good)
    assert registry.loc[0, "q_threshold"] == 0.10


def test_external_comparison_requires_registered_features():
    discovery = pd.DataFrame({
        "feature": ["PPAR_signaling", "translation_machinery"],
        "effect": [1.0, 1.0],
        "q_threshold": [0.10, 0.10],
    })
    missing_external = pd.DataFrame({
        "feature": ["PPAR_signaling"],
        "effect": [0.5],
        "q_value": [0.01],
    })
    with pytest.raises(ValueError, match="missing registered features"):
        compare_directional_replication(discovery, missing_external)

    external = pd.DataFrame({
        "feature": ["PPAR_signaling", "translation_machinery"],
        "effect": [0.5, -0.5],
        "q_value": [0.01, 0.01],
    })
    result = compare_directional_replication(discovery, external)
    assert result.loc[result["feature"] == "PPAR_signaling", "replicated"].item()
    assert not result.loc[result["feature"] == "translation_machinery", "replicated"].item()


def test_external_context_rows_are_not_directional_replication_claims():
    discovery = pd.DataFrame({
        "feature": ["ECM_remodeling"],
        "effect": [0.0],
        "q_threshold": [0.10],
    })
    external = pd.DataFrame({
        "feature": ["ECM_remodeling"],
        "effect": [-0.5],
        "q_value": [0.01],
    })
    result = compare_directional_replication(discovery, external)
    assert result.loc[0, "external_pass_q"]
    assert result.loc[0, "context_detected"]
    assert not result.loc[0, "directional_claim"]
    assert not result.loc[0, "replicated"]
