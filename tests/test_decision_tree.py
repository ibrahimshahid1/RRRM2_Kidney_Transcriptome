from src.networks.manuscript_decision import decide_manuscript_branch


def test_decision_tree_strong_network_branch():
    decision = decide_manuscript_branch({
        "within_rrrm2_stability_passed": True,
        "within_rrrm2_projection_significant": True,
        "cross_osdr_recurrence": True,
    })
    assert decision.branch == "strong_network_aging_paper"


def test_decision_tree_external_axis_branch():
    decision = decide_manuscript_branch({
        "within_rrrm2_stability_passed": False,
        "within_rrrm2_projection_significant": False,
        "external_axis_significant": True,
    })
    assert decision.branch == "external_aging_axis_alignment"


def test_decision_tree_negative_branch():
    decision = decide_manuscript_branch({})
    assert decision.branch == "negative_or_methods_constraint"


def test_decision_tree_module_only_branch():
    decision = decide_manuscript_branch({
        "within_rrrm2_stability_passed": True,
        "within_rrrm2_projection_significant": True,
        "module_activity_only_positive": True,
        "any_projection_signal": True,
    })
    assert decision.branch == "modest_module_activity_paper"
