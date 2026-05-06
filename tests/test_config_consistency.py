from pathlib import Path

import yaml

from src.networks.shared_topology import DEFAULT_TOPK
from src.statistics.permutation_bootstrap import DEFAULT_B_BOOT, DEFAULT_K_PERM
from src.validation.cross_validation import DEFAULT_RF_MAX_DEPTH


def test_config_matches_executed_defaults():
    cfg = yaml.safe_load(Path("config/hyperparameters.yaml").read_text())
    assert cfg["network"]["topk"]["k"] == DEFAULT_TOPK
    assert cfg["cross_validation"]["classifier"]["params"]["max_depth"] == DEFAULT_RF_MAX_DEPTH
    assert cfg["statistics"]["permutation"]["n_iterations"] == DEFAULT_K_PERM
    assert cfg["statistics"]["bootstrap"]["n_iterations"] == DEFAULT_B_BOOT
    assert cfg["statistics"]["permutation"]["strict_fdr"] is None
    assert "westfall_young" not in cfg["statistics"]["permutation"].get("focused_testing_modes", [])
