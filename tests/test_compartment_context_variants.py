"""WP-S1: family-wide scoring variants of the compartment-context family
(--method {mean,median}, --gene-universe {per_mission,common}). Synthetic data only."""

from argparse import Namespace
import hashlib
import inspect
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

from src.clinical_axes.data import MissionData
import scripts.clinical_axes.run_compartment_context as rcc
from scripts.clinical_axes.run_compartment_context import (
    load_family,
    prepare_family,
    seed_offset,
    variant_output_dir,
)

CONFIG_SEED = 1000

# (block, n_flight, n_control) per block, one list per synthetic mission.
LAYOUTS = {
    "SYN-A": [("all", 6, 6)],
    "SYN-B": [("d25", 5, 5), ("d75", 5, 4)],
    "SYN-C": [("all", 10, 10)],
}

SET_A = [f"A{i}" for i in range(10)]          # fully common
SET_B = [f"B{i}" for i in range(12)]          # B0-B4 not common -> 7 common genes
SET_C = [f"C{i}" for i in range(9)]           # fully common
SET_SMALL = [f"S{i}" for i in range(6)]       # < 8 defined: never evaluable
B_INELIGIBLE_IN_A = ["B0", "B1"]              # zero counts in SYN-A (CPM-ineligible)
B_MISSING_FROM_C_EXPRESSION = ["B2", "B3", "B4"]  # eligible in SYN-C counts, absent from its VST


def _mission(name, layout, rng):
    condition, block = [], []
    for b, n_flight, n_control in layout:
        condition += ["FLT"] * n_flight + ["GC"] * n_control
        block += [b] * (n_flight + n_control)
    animals = [f"{name}_{i}" for i in range(len(condition))]
    metadata = pd.DataFrame(
        {"condition": condition, "block": block, "source_sample": animals}, index=animals
    )
    genes = SET_A + SET_B + SET_C + SET_SMALL
    flight = (metadata["condition"] == "FLT").to_numpy(float)
    # heavy-tailed gene noise so median and mean scores differ
    values = rng.standard_t(2, size=(len(genes), len(animals))) + 0.6 * flight[None, :] + 7.0
    expression = pd.DataFrame(values, index=pd.Index(genes, name="gene"), columns=animals)
    counts = pd.DataFrame(
        rng.poisson(400, size=(len(genes), len(animals))).astype(float),
        index=pd.Index(genes, name="gene"), columns=animals,
    )
    if name == "SYN-A":
        counts.loc[B_INELIGIBLE_IN_A] = 0.0
    if name == "SYN-C":
        expression = expression.drop(index=B_MISSING_FROM_C_EXPRESSION)
    data = MissionData(
        mission=name, expression=expression, counts=counts, metadata=metadata,
        qc=pd.DataFrame(index=animals),
    )
    data.validate()
    return data


def _missions(seed=0):
    rng = np.random.default_rng(seed)
    return {name: _mission(name, layout, rng) for name, layout in LAYOUTS.items()}


def _gene_sets(tmp_path) -> Path:
    rows = []
    for gene_set, genes, compartment in [
        ("alpha__all", SET_A, "alpha"),
        ("beta__all", SET_B, "beta"),
        ("gamma__all", SET_C, "gamma"),
        ("small__all", SET_SMALL, "small"),
    ]:
        rows += [
            {"gene_set": gene_set, "gene_symbol": g, "final_for_testing": True,
             "report_compartment": compartment, "tier": "all"}
            for g in genes
        ]
    path = tmp_path / "tiers.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def _config(tmp_path) -> Path:
    id_map = tmp_path / "id_map.tsv"
    id_map.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG1\tA0\n")
    config = {
        "seed": CONFIG_SEED,
        "eligibility": {"cpm_threshold": 0.1},
        "gene_mapping": {"path": str(id_map), "annotation_fallback": "unused.csv"},
    }
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump(config))
    return path


def _args(tmp_path, **overrides):
    args = dict(
        config=_config(tmp_path), gene_sets=_gene_sets(tmp_path), results=tmp_path / "run",
        permutations=99, chunk_size=50, add_podocyte_disjoint_variants=False,
        osd253_strain=None, method="mean", gene_universe="per_mission", gene_map=None,
        prereg=tmp_path / "missing_prereg.yaml", allow_uncommitted_prereg=False,
    )
    args.update(overrides)
    return Namespace(**args)


@pytest.fixture
def synthetic(monkeypatch):
    missions = _missions()
    seen = {}

    def fake_loader(config, root):
        seen["gene_map"] = config["gene_mapping"]["path"]
        return dict(missions)

    monkeypatch.setattr(rcc, "load_primary_missions", fake_loader)
    return missions, seen


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _file_digests(folder: Path) -> dict[str, str]:
    return {p.name: _digest(p) for p in folder.glob("compartment_context_*") if p.is_file()}


# --- contracts -------------------------------------------------------------------

def test_default_seed_offset_is_100_and_variants_use_101_to_103():
    assert seed_offset() == 100
    assert seed_offset("mean", "per_mission") == 100
    assert seed_offset("median", "per_mission") == 101
    assert seed_offset("mean", "common") == 102
    assert seed_offset("median", "common") == 103
    with pytest.raises(ValueError):
        seed_offset("trimmed", "per_mission")


def test_variant_output_dirs_are_separate():
    root = Path("run")
    assert variant_output_dir(root) == root
    assert variant_output_dir(root, "median") == root / "compartment_context_median"
    assert variant_output_dir(root, "mean", "common") == root / "compartment_context_common"
    assert variant_output_dir(root, "median", "common") == root / "compartment_context_median_common"
    # registry stages pass the subfolder explicitly; it is not nested twice
    explicit = root / "compartment_context_median"
    assert variant_output_dir(explicit, "median") == explicit
    assert variant_output_dir(explicit, "mean", "common") == explicit / "compartment_context_common"


def test_load_family_signature_unchanged():
    params = inspect.signature(load_family).parameters
    assert list(params) == ["path", "minimum_defined"]
    assert params["minimum_defined"].default == 8


def test_cli_defaults_and_choices(monkeypatch):
    captured = {}
    monkeypatch.setattr(rcc, "run", lambda args: captured.setdefault("args", args))
    monkeypatch.setattr(sys, "argv", ["run_compartment_context.py"])
    rcc.main()
    args = captured["args"]
    assert (args.method, args.gene_universe, args.gene_map) == ("mean", "per_mission", None)
    assert args.allow_uncommitted_prereg is False
    monkeypatch.setattr(sys, "argv", ["x", "--method", "trimmed"])
    with pytest.raises(SystemExit):
        rcc.main()


# --- family construction ---------------------------------------------------------

def test_common_restriction_drops_set_with_fewer_than_8_common_genes(tmp_path):
    missions = _missions()
    tiers = _gene_sets(tmp_path)
    per_mission, _ = prepare_family(missions, tiers, 0.1)
    common, audit = prepare_family(missions, tiers, 0.1, gene_universe="common")

    # per mission, beta is evaluable (>= 8 eligible genes in every mission) ...
    assert list(per_mission) == ["alpha__all", "beta__all", "gamma__all"]
    # ... but only 7 of its genes are eligible everywhere and observed in every VST
    assert list(common) == ["alpha__all", "gamma__all"]
    row = audit.set_index("gene_set").loc["beta__all"]
    assert row["n_common_universe_genes"] == 7
    assert row["n_dropped_common_universe"] == 5
    assert not row["common_universe_kept"]
    assert not row["cross_mission_evaluable"]
    kept = audit.set_index("gene_set").loc["alpha__all"]
    assert kept["n_common_universe_genes"] == 10 and kept["n_dropped_common_universe"] == 0
    assert kept["common_universe_kept"]
    # retained sets require every retained gene in every mission
    for gene_set, spec in common.items():
        sub = spec["subdomains"]["atlas_markers"]
        assert sub["minimum_present"] == len(sub["genes"]) >= 8
    # never-evaluable definitions stay out and are recorded as not kept
    assert not audit.set_index("gene_set").loc["small__all", "common_universe_kept"]


def test_common_universe_requires_eligibility_and_expression_presence():
    universe = rcc.common_gene_universe(_missions(), 0.1)
    assert not set(B_INELIGIBLE_IN_A) & universe          # CPM-ineligible in one mission
    assert not set(B_MISSING_FROM_C_EXPRESSION) & universe  # absent from one VST matrix
    assert set(SET_A) | set(SET_C) | {"B5", "B11"} <= universe


# --- end-to-end runs ---------------------------------------------------------------

def test_median_family_identical_to_mean_family_and_written_separately(tmp_path, synthetic):
    prereg = tmp_path / "prereg.yaml"
    prereg.write_text("analysis_id: test\n")
    mean_args = _args(tmp_path)
    rcc.run(mean_args)
    primary = mean_args.results
    primary_digests = _file_digests(primary)

    median_args = _args(tmp_path, method="median", prereg=prereg, allow_uncommitted_prereg=True)
    rcc.run(median_args)
    median_dir = primary / "compartment_context_median"

    # primary outputs untouched by the variant run
    assert _file_digests(primary) == primary_digests
    assert sorted(p.name for p in median_dir.iterdir()) == sorted(primary_digests)

    mean_meta = pd.read_csv(primary / "compartment_context_meta_results.tsv", sep="\t")
    median_meta = pd.read_csv(median_dir / "compartment_context_meta_results.tsv", sep="\t")
    assert list(mean_meta["axis"]) == list(median_meta["axis"]) == [
        "alpha__all", "beta__all", "gamma__all"
    ]
    assert list(mean_meta.columns) == list(median_meta.columns)
    pd.testing.assert_frame_equal(
        pd.read_csv(primary / "compartment_context_definition_audit.tsv", sep="\t"),
        pd.read_csv(median_dir / "compartment_context_definition_audit.tsv", sep="\t"),
    )
    # the scoring really changed
    assert not np.allclose(mean_meta["estimate"], median_meta["estimate"])

    mean_manifest = json.loads((primary / "compartment_context_manifest.json").read_text())
    median_manifest = json.loads((median_dir / "compartment_context_manifest.json").read_text())
    assert (mean_manifest["scoring_method"], mean_manifest["gene_universe"]) == ("mean", "per_mission")
    assert mean_manifest["seed"] == CONFIG_SEED + 100 and "preregistration" not in mean_manifest
    assert (median_manifest["scoring_method"], median_manifest["gene_universe"]) == ("median", "per_mission")
    assert median_manifest["seed"] == CONFIG_SEED + 101
    assert median_manifest["preregistration"]["prereg_guard"].startswith("bypassed")
    assert median_manifest["preregistration"]["prereg_sha256"] == _digest(prereg)
    assert median_manifest["n_evaluable_sets"] == mean_manifest["n_evaluable_sets"] == 3


def test_common_variant_scores_only_common_genes(tmp_path, synthetic):
    args = _args(tmp_path, gene_universe="common", allow_uncommitted_prereg=True)
    rcc.run(args)
    out = args.results / "compartment_context_common"
    meta = pd.read_csv(out / "compartment_context_meta_results.tsv", sep="\t")
    assert list(meta["axis"]) == ["alpha__all", "gamma__all"]
    coverage = pd.read_csv(out / "compartment_context_gene_coverage.tsv", sep="\t")
    assert (coverage["n_used"] == coverage["n_requested"]).all()
    assert (coverage["n_used"] == coverage["minimum_required"]).all()
    # the same genes in every mission
    assert coverage.groupby("axis")["genes_used"].nunique().eq(1).all()
    audit = pd.read_csv(out / "compartment_context_definition_audit.tsv", sep="\t")
    assert {"n_common_universe_genes", "n_dropped_common_universe", "common_universe_kept"} <= set(audit.columns)
    manifest = json.loads((out / "compartment_context_manifest.json").read_text())
    assert manifest["gene_universe"] == "common" and manifest["seed"] == CONFIG_SEED + 102
    assert manifest["n_evaluable_sets"] == 2

    both = _args(tmp_path, method="median", gene_universe="common", allow_uncommitted_prereg=True)
    rcc.run(both)
    manifest = json.loads(
        (both.results / "compartment_context_median_common" / "compartment_context_manifest.json").read_text()
    )
    assert manifest["seed"] == CONFIG_SEED + 103


def test_variant_refuses_without_committed_preregistration(tmp_path, synthetic):
    with pytest.raises(SystemExit, match="preregistration"):
        rcc.run(_args(tmp_path, method="median"))  # YAML missing
    outside = tmp_path / "prereg.yaml"  # exists but is not a committed repo file
    outside.write_text("analysis_id: test\n")
    with pytest.raises(SystemExit):
        rcc.run(_args(tmp_path, gene_universe="common", prereg=outside))
    assert not (tmp_path / "run" / "compartment_context_common").exists()
    # the primary family needs no preregistration
    rcc.run(_args(tmp_path))


def test_gene_map_override_is_in_memory_only(tmp_path, synthetic):
    _, seen = synthetic
    args = _args(tmp_path)
    alternative = tmp_path / "alt_map.tsv"
    alternative.write_text("ensembl_gene_id\tmgi_symbol\nENSMUSG2\tA1\n")
    before = args.config.read_text()
    args.gene_map = alternative
    rcc.run(args)
    assert seen["gene_map"] == str(alternative.resolve())
    assert args.config.read_text() == before
    manifest = json.loads((args.results / "compartment_context_manifest.json").read_text())
    assert manifest["gene_map"] == str(alternative.resolve())
    assert manifest["gene_map_sha256"] == _digest(alternative)
    assert manifest["gene_map_overridden"] is True
    with pytest.raises(FileNotFoundError):
        rcc.run(_args(tmp_path, gene_map=tmp_path / "nope.tsv"))
