from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.networks.procrustes import load_configured_anchor_indices


def _write_anchor_fixture(tmp_path: Path, n: int = 25):
    cfg = tmp_path / "anchors.yaml"
    genes = "\n".join(f"      - Gene{i}" for i in range(n))
    cfg.write_text(
        "anchors:\n"
        "  stable:\n"
        "    role: primary_alignment\n"
        "    genes:\n"
        f"{genes}\n"
        "validation:\n"
        "  minimum_anchors: 20\n"
        "  warn_threshold: 50\n"
    )
    id_map = tmp_path / "id_map.tsv"
    pd.DataFrame({
        "ensembl_gene_id": [f"ENS{i:05d}" for i in range(n)],
        "mgi_symbol": [f"Gene{i}" for i in range(n)],
        "biotype": ["protein_coding"] * n,
    }).to_csv(id_map, sep="\t", index=False)
    return cfg, id_map


def test_configured_anchors_resolve_and_fail_below_minimum(tmp_path):
    cfg, id_map = _write_anchor_fixture(tmp_path)
    genes = [f"ENS{i:05d}" for i in range(25)]
    anchors = load_configured_anchor_indices(genes, cfg, id_map, tmp_path, min_anchors=20, warn_threshold=50)
    assert anchors.dtype == np.int32
    assert anchors.size == 25
    assert (tmp_path / "anchors_report.tsv").exists()

    with pytest.raises(RuntimeError, match="below minimum"):
        load_configured_anchor_indices(genes[:8], cfg, id_map, tmp_path, min_anchors=20, warn_threshold=50)
