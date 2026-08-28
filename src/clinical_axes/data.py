"""Load and harmonize the frozen cross-mission kidney RNA contrasts.

All public functions in this module are label-transparent: sample membership is
derived only from the frozen configuration and repository metadata.  No effect
estimate is used to select a sample or gene.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Mapping

import numpy as np
import pandas as pd


EPS = 1e-12


@dataclass(frozen=True)
class MissionData:
    """One biological value per animal for a frozen flight-control contrast."""

    mission: str
    expression: pd.DataFrame
    counts: pd.DataFrame
    metadata: pd.DataFrame
    qc: pd.DataFrame
    context: str = "terminal"

    def validate(self) -> None:
        samples = self.metadata.index.astype(str)
        if samples.duplicated().any():
            raise ValueError(f"{self.mission}: duplicate biological animal IDs")
        if list(self.expression.columns.astype(str)) != list(samples):
            raise ValueError(f"{self.mission}: expression/metadata order mismatch")
        if list(self.counts.columns.astype(str)) != list(samples):
            raise ValueError(f"{self.mission}: counts/metadata order mismatch")
        required = {"condition", "block", "source_sample"}
        missing = sorted(required - set(self.metadata.columns))
        if missing:
            raise ValueError(f"{self.mission}: missing metadata columns {missing}")
        levels = set(self.metadata["condition"])
        if levels != {"FLT", "GC"}:
            raise ValueError(f"{self.mission}: expected FLT/GC; found {levels}")
        for block, sub in self.metadata.groupby("block", sort=False):
            counts = sub["condition"].value_counts()
            if min(counts.get("FLT", 0), counts.get("GC", 0)) < 2:
                raise ValueError(
                    f"{self.mission}/{block}: need >=2 animals in FLT and GC"
                )


def _strip_ensembl_version(values: pd.Index) -> pd.Index:
    return values.astype(str).str.replace(r"\.\d+$", "", regex=True)


def load_symbol_map(path: Path, annotation_fallback: Path | None = None) -> dict[str, str]:
    table = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    required = {"ensembl_gene_id", "mgi_symbol"}
    if not required.issubset(table.columns):
        raise ValueError(f"ID map lacks {sorted(required - set(table.columns))}")
    table = table.dropna(subset=["ensembl_gene_id", "mgi_symbol"])
    table["ensembl_gene_id"] = _strip_ensembl_version(
        pd.Index(table["ensembl_gene_id"])
    )
    mapping = dict(zip(table["ensembl_gene_id"], table["mgi_symbol"]))
    if annotation_fallback is not None:
        fallback = pd.read_csv(
            annotation_fallback, usecols=["ENSEMBL", "SYMBOL"], dtype=str
        ).dropna()
        fallback["ENSEMBL"] = _strip_ensembl_version(
            pd.Index(fallback["ENSEMBL"])
        )
        # The OSDR differential-expression table supplies annotations for genes
        # outside the historical 14k repository universe.  The curated map
        # remains authoritative when both sources contain an ID.
        for ens, symbol in zip(fallback["ENSEMBL"], fallback["SYMBOL"]):
            mapping.setdefault(ens, symbol)
    return mapping


def load_matrix(
    path: Path,
    id_to_symbol: Mapping[str, str],
    *,
    aggregate: str,
) -> pd.DataFrame:
    """Read an OSDR gene matrix and collapse Ensembl IDs to MGI symbols."""
    frame = pd.read_csv(path, index_col=0)
    frame.index = _strip_ensembl_version(frame.index)
    mapped = pd.Index([id_to_symbol.get(g, g) for g in frame.index], name="gene")
    frame.index = mapped
    frame = frame.apply(pd.to_numeric, errors="coerce")
    if aggregate == "mean":
        frame = frame.groupby(level=0, sort=False).mean()
    elif aggregate == "sum":
        frame = frame.groupby(level=0, sort=False).sum(min_count=1)
    else:
        raise ValueError(f"Unsupported aggregate: {aggregate}")
    return frame


def _biological_id(sample: str) -> str:
    value = re.sub(r"_techrep\d+$", "", str(sample))
    value = re.sub(r"_(?:totRNA|mRNA|UPX)$", "", value)
    return value


def _collapse_columns(
    frame: pd.DataFrame,
    samples: list[str],
    *,
    aggregate: str,
    column_to_animal: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    missing = sorted(set(samples) - set(frame.columns))
    if missing:
        raise ValueError(f"Matrix missing selected samples: {missing[:8]}")
    selected = frame.loc[:, samples].copy()
    selected.columns = [
        column_to_animal.get(s, _biological_id(s))
        if column_to_animal is not None
        else _biological_id(s)
        for s in selected.columns
    ]
    if aggregate == "mean":
        selected = selected.T.groupby(level=0, sort=False).mean().T
    elif aggregate == "sum":
        selected = selected.T.groupby(level=0, sort=False).sum(min_count=1).T
    else:
        raise ValueError(aggregate)
    return selected


def _collapse_metadata(meta: pd.DataFrame) -> pd.DataFrame:
    meta = meta.copy()
    if "animal" not in meta:
        meta["animal"] = meta["source_sample"].map(_biological_id)
    rows = []
    for animal, sub in meta.groupby("animal", sort=False):
        for col in ("condition", "block"):
            if sub[col].nunique(dropna=False) != 1:
                raise ValueError(f"{animal}: technical replicates disagree on {col}")
        rows.append(
            {
                "animal": animal,
                "condition": sub["condition"].iloc[0],
                "block": sub["block"].iloc[0],
                "source_sample": "|".join(sub["source_sample"].astype(str)),
            }
        )
    return pd.DataFrame(rows).set_index("animal", drop=True)


def _collapse_qc(
    qc: pd.DataFrame,
    eligible_animals: set[str],
    sample_to_animal: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    if qc.empty:
        return pd.DataFrame(index=pd.Index(sorted(eligible_animals)))
    qc = qc.copy()
    if sample_to_animal is None:
        qc["animal"] = qc["sample"].map(_biological_id)
    else:
        qc["animal"] = qc["sample"].map(
            lambda sample: sample_to_animal.get(sample, _biological_id(sample))
        )
    qc = qc[qc["animal"].isin(eligible_animals)]
    numeric = qc.select_dtypes(include=[np.number]).columns.tolist()
    keep = [
        c
        for c in (
            "ratio_genebody_cov_3_to_5",
            "mean_genebody_cov_5_20",
            "mean_genebody_cov_80_95",
            "rin",
            "read_depth",
            "uniquely_mapped_percent",
        )
        if c in numeric
    ]
    return qc.groupby("animal", sort=False)[keep].mean()


def _finalize(
    mission: str,
    expression: pd.DataFrame,
    counts: pd.DataFrame,
    metadata: pd.DataFrame,
    qc: pd.DataFrame,
    *,
    context: str = "terminal",
    qc_sample_to_animal: Mapping[str, str] | None = None,
) -> MissionData:
    order = metadata.index.astype(str).tolist()
    collapsed_qc = _collapse_qc(qc, set(order), qc_sample_to_animal)
    expression = expression.loc[:, order]
    counts = counts.loc[:, order]
    collapsed_qc = collapsed_qc.reindex(order)
    data = MissionData(
        mission=mission,
        expression=expression,
        counts=counts,
        metadata=metadata,
        qc=collapsed_qc,
        context=context,
    )
    data.validate()
    return data


def _load_simple(
    mission: str,
    spec: Mapping[str, object],
    root: Path,
    id_to_symbol: Mapping[str, str],
    *,
    context: str = "terminal",
) -> MissionData:
    expression_all = load_matrix(root / str(spec["vst"]), id_to_symbol, aggregate="mean")
    counts_all = load_matrix(root / str(spec["counts"]), id_to_symbol, aggregate="sum")
    runsheet = pd.read_csv(root / str(spec["runsheet"]), dtype=str)
    label_map = {"Space Flight": "FLT", "Ground Control": "GC"}
    runsheet = runsheet[
        runsheet["Factor Value[Spaceflight]"].isin(label_map)
    ].copy()
    runsheet["condition"] = runsheet["Factor Value[Spaceflight]"].map(label_map)
    runsheet["block"] = "all"
    runsheet["source_sample"] = runsheet["Sample Name"]
    raw_samples = [s for s in runsheet["source_sample"] if s in expression_all.columns]
    runsheet = runsheet[runsheet["source_sample"].isin(raw_samples)]
    metadata = _collapse_metadata(runsheet[["source_sample", "condition", "block"]])
    expression = _collapse_columns(expression_all, raw_samples, aggregate="mean")
    counts = _collapse_columns(counts_all, raw_samples, aggregate="sum")
    qc = pd.read_csv(root / str(spec["qc"]))
    return _finalize(
        mission, expression, counts, metadata, qc, context=context
    )


def _load_osd253(
    spec: Mapping[str, object],
    root: Path,
    id_to_symbol: Mapping[str, str],
    *,
    control: str = "original",
    strain: str = "C57BL/6J",
) -> MissionData:
    expression_all = load_matrix(root / str(spec["vst"]), id_to_symbol, aggregate="mean")
    counts_all = load_matrix(root / str(spec["counts"]), id_to_symbol, aggregate="sum")
    runsheet = pd.read_csv(root / str(spec["runsheet"]), dtype=str)
    runsheet = runsheet[runsheet["Factor Value[Strain]"] == strain].copy()
    if runsheet.empty:
        raise ValueError(f"OSD-253: no samples for strain {strain}")
    if control == "original":
        allowed = {"Space Flight": "FLT", "Ground Control": "GC"}
    elif control == "rerun_white":
        allowed = {"Space Flight": "FLT", "Ground Control Rerun": "GC"}
        runsheet = runsheet[
            (runsheet["Factor Value[Spaceflight]"] != "Ground Control Rerun")
            | (runsheet["Factor Value[Treatment]"] == "White light")
        ]
    else:
        raise ValueError(f"Unknown OSD-253 control: {control}")
    runsheet = runsheet[
        runsheet["Factor Value[Spaceflight]"].isin(allowed)
    ].copy()
    runsheet["condition"] = runsheet["Factor Value[Spaceflight]"].map(allowed)
    runsheet["block"] = (
        runsheet["Factor Value[Duration]"]
        .str.extract(r"(25|75)", expand=False)
        .map({"25": "day25", "75": "day75"})
    )
    runsheet["source_sample"] = runsheet["Sample Name"]
    sample_to_animal: dict[str, str] = {}
    for sample, condition in zip(runsheet["source_sample"], runsheet["condition"]):
        if control == "rerun_white" and condition == "GC":
            match = re.search(
                r"GCrerun_(25|75)days_White_Rep(\d+)_", str(sample)
            )
            if match is None:
                raise ValueError(f"Cannot derive OSD-253 rerun animal from {sample}")
            animal = f"RR7_C57_GCrerun_{match.group(1)}days_White_Rep{match.group(2)}"
        else:
            animal = _biological_id(sample)
        sample_to_animal[str(sample)] = animal
    runsheet["animal"] = runsheet["source_sample"].map(sample_to_animal)
    if runsheet["block"].isna().any():
        raise ValueError("OSD-253: unrecognized duration")
    raw_samples = runsheet["source_sample"].tolist()
    metadata = _collapse_metadata(
        runsheet[["source_sample", "condition", "block", "animal"]]
    )
    expression = _collapse_columns(
        expression_all,
        raw_samples,
        aggregate="mean",
        column_to_animal=sample_to_animal,
    )
    counts = _collapse_columns(
        counts_all,
        raw_samples,
        aggregate="sum",
        column_to_animal=sample_to_animal,
    )
    qc = pd.read_csv(root / str(spec["qc"]))
    mission = "OSD-253" if strain == "C57BL/6J" else f"OSD-253-{strain.split('/')[0]}"
    return _finalize(
        mission,
        expression,
        counts,
        metadata,
        qc,
        qc_sample_to_animal=sample_to_animal,
    )


def _load_osd462(
    spec: Mapping[str, object],
    root: Path,
    id_to_symbol: Mapping[str, str],
    *,
    preparation: str = "totRNA",
) -> MissionData:
    base = root / "data/external/osdr/OSD-462"
    if preparation == "totRNA":
        vst_path = root / str(spec["vst"])
        count_path = root / str(spec["counts"])
        qc_path = root / str(spec["qc"])
    else:
        vst_path = base / f"GLDS-462_rna_seq_VST_Counts_{preparation}_GLbulkRNAseq.csv"
        count_path = base / f"GLDS-462_rna_seq_RSEM_Unnormalized_Counts_{preparation}_GLbulkRNAseq.csv"
        qc_path = base / f"GLDS-462_rna_seq_qc_metrics_{preparation}_GLbulkRNAseq.csv"
    expression_all = load_matrix(vst_path, id_to_symbol, aggregate="mean")
    counts_all = load_matrix(count_path, id_to_symbol, aggregate="sum")
    pattern = re.compile(r"_WT_(FLT|GC)_")
    raw_samples = [s for s in expression_all.columns if pattern.search(s)]
    rows = []
    for sample in raw_samples:
        match = pattern.search(sample)
        rows.append(
            {
                "source_sample": sample,
                "condition": match.group(1),
                "block": "all",
            }
        )
    metadata = _collapse_metadata(pd.DataFrame(rows))
    expression = _collapse_columns(expression_all, raw_samples, aggregate="mean")
    counts = _collapse_columns(counts_all, raw_samples, aggregate="sum")
    qc = pd.read_csv(qc_path)
    return _finalize("OSD-462", expression, counts, metadata, qc)


def _load_osd771(
    spec: Mapping[str, object],
    root: Path,
    id_to_symbol: Mapping[str, str],
    *,
    arm: str = "ISS-T",
) -> MissionData:
    expression_all = load_matrix(root / str(spec["vst"]), id_to_symbol, aggregate="mean")
    counts_all = load_matrix(root / str(spec["counts"]), id_to_symbol, aggregate="sum")
    metadata_all = pd.read_csv(root / str(spec["metadata"]), sep="\t", dtype=str)
    label_map = {"Space Flight": "FLT", "Ground Control": "GC"}
    metadata_all = metadata_all[
        metadata_all["Factor Value[Spaceflight]"].isin(label_map)
        & metadata_all["Sample Name"].str.contains(f"_{arm}_", regex=False)
    ].copy()
    metadata_all["condition"] = metadata_all["Factor Value[Spaceflight]"].map(label_map)
    metadata_all["block"] = metadata_all["Sample Name"].str.extract(
        r"_(YNG|OLD)_", expand=False
    )
    metadata_all["source_sample"] = metadata_all["Sample Name"]
    raw_samples = metadata_all["source_sample"].tolist()
    metadata = _collapse_metadata(
        metadata_all[["source_sample", "condition", "block"]]
    )
    expression = _collapse_columns(expression_all, raw_samples, aggregate="mean")
    counts = _collapse_columns(counts_all, raw_samples, aggregate="sum")
    qc = pd.read_csv(root / str(spec["qc"]))
    context = "terminal" if arm == "ISS-T" else "long_recovery"
    mission = "OSD-771" if arm == "ISS-T" else "OSD-771-LAR"
    return _finalize(mission, expression, counts, metadata, qc, context=context)


def load_primary_missions(
    config: Mapping[str, object], root: Path
) -> dict[str, MissionData]:
    mapping = config["gene_mapping"]
    id_to_symbol = load_symbol_map(
        root / str(mapping["path"]),
        root / str(mapping["annotation_fallback"]),
    )
    specs = config["primary_missions"]
    missions = {
        "OSD-102": _load_simple("OSD-102", specs["OSD-102"], root, id_to_symbol),
        "OSD-163": _load_simple("OSD-163", specs["OSD-163"], root, id_to_symbol),
        "OSD-253": _load_osd253(specs["OSD-253"], root, id_to_symbol),
        "OSD-462": _load_osd462(specs["OSD-462"], root, id_to_symbol),
        "OSD-771": _load_osd771(specs["OSD-771"], root, id_to_symbol),
    }
    return missions


def load_moderator_missions(
    config: Mapping[str, object], root: Path
) -> dict[str, MissionData]:
    mapping = config["gene_mapping"]
    id_to_symbol = load_symbol_map(
        root / str(mapping["path"]),
        root / str(mapping["annotation_fallback"]),
    )
    primary = config["primary_missions"]
    mods = config["moderator_cohorts"]
    return {
        "OSD-513": _load_simple(
            "OSD-513", mods["OSD-513"], root, id_to_symbol, context="live_return"
        ),
        "OSD-771-LAR": _load_osd771(
            primary["OSD-771"], root, id_to_symbol, arm="LAR"
        ),
    }


def load_osd462_preparation(
    config: Mapping[str, object], root: Path, preparation: str
) -> MissionData:
    mapping = config["gene_mapping"]
    id_to_symbol = load_symbol_map(
        root / str(mapping["path"]),
        root / str(mapping["annotation_fallback"]),
    )
    return _load_osd462(
        config["primary_missions"]["OSD-462"],
        root,
        id_to_symbol,
        preparation=preparation,
    )


def load_osd253_control_sensitivity(
    config: Mapping[str, object], root: Path
) -> MissionData:
    mapping = config["gene_mapping"]
    id_to_symbol = load_symbol_map(
        root / str(mapping["path"]),
        root / str(mapping["annotation_fallback"]),
    )
    return _load_osd253(
        config["primary_missions"]["OSD-253"],
        root,
        id_to_symbol,
        control="rerun_white",
    )


def load_osd253_strain_sensitivity(
    config: Mapping[str, object], root: Path, strain: str = "C3H/HeJ"
) -> MissionData:
    """Load the matched original-control OSD-253 contrast for another strain."""
    mapping = config["gene_mapping"]
    id_to_symbol = load_symbol_map(
        root / str(mapping["path"]),
        root / str(mapping["annotation_fallback"]),
    )
    return _load_osd253(
        config["primary_missions"]["OSD-253"],
        root,
        id_to_symbol,
        control="original",
        strain=strain,
    )


def cpm_eligible_genes(data: MissionData, threshold: float = 0.1) -> set[str]:
    """CPM >= threshold in at least half of either contrasted arm."""
    library = data.counts.sum(axis=0)
    if (library <= EPS).any():
        raise ValueError(f"{data.mission}: nonpositive count-library size")
    cpm = data.counts.divide(library, axis=1) * 1_000_000.0
    eligible = pd.Series(False, index=cpm.index)
    for condition in ("FLT", "GC"):
        samples = data.metadata.index[data.metadata["condition"] == condition]
        needed = int(np.ceil(len(samples) / 2))
        eligible |= (cpm.loc[:, samples] >= threshold).sum(axis=1) >= needed
    return set(eligible.index[eligible])
