from pathlib import Path

import pandas as pd
import pytest

from src.statistics.full_pipeline_permutation import (
    validate_command_template,
    write_permutation_manifest,
)


def test_full_pipeline_permutation_writes_permuted_metadata_and_uses_it(tmp_path: Path):
    meta = pd.DataFrame([
        {"sample": "s1", "Age": "YNG", "Arm": "LAR", "EnvGroup": "FLT"},
        {"sample": "s2", "Age": "YNG", "Arm": "LAR", "EnvGroup": "FLT"},
        {"sample": "s3", "Age": "YNG", "Arm": "LAR", "EnvGroup": "GC"},
        {"sample": "s4", "Age": "YNG", "Arm": "LAR", "EnvGroup": "GC"},
    ])
    meta_path = tmp_path / "meta.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)

    manifest = write_permutation_manifest(
        outdir=tmp_path / "out",
        top_hits=["ENS1"],
        k_perm=1,
        seed=0,
        command_template="python -m fake --meta={meta} --outdir={outdir}",
        meta_path=meta_path,
        rtech_path=tmp_path / "Rtech.tsv.gz",
        id_map_path=tmp_path / "id_map.tsv",
        anchor_config=tmp_path / "anchors.yaml",
        max_genes=2500,
        topk=80,
        num_seeds=10,
        seed0=0,
        strata_cols=["Age", "Arm"],
    )

    perm_meta = Path(manifest.loc[0, "meta_permuted"])
    assert perm_meta.exists()
    assert str(perm_meta) in manifest.loc[0, "command"]
    assert "--dry-run" not in manifest.loc[0, "command"]
    shuffled = pd.read_csv(perm_meta, sep="\t", compression="gzip")
    assert sorted(shuffled["EnvGroup"].tolist()) == ["FLT", "FLT", "GC", "GC"]


def test_full_pipeline_command_template_must_consume_metadata():
    with pytest.raises(ValueError, match="\\{meta\\}"):
        validate_command_template("python -m fake")
    with pytest.raises(ValueError, match="dry-run"):
        validate_command_template("python -m fake --meta={meta} --dry-run")
