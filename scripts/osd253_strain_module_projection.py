"""Metadata-aware OSD-253 strain-stratified WGCNA module projection.

This is the decision analysis for the RRRM-2 grey60/TLR4 question:

1. Parse ISA sample and assay metadata, not filename regexes alone.
2. Score RRRM-2 WGCNA modules in OSD-253.
3. Test FLT-minus-control module shifts within Strain x Duration strata.
4. Test whether C3H/HeJ attenuates the C57BL/6J flight effect.
5. Run original-GC and white-light rerun-control sensitivity analyses separately.

The rerun-control analysis fixes the OSD-253 light-wavelength mismatch, but it
introduces sequencing/read-length/run differences. It is therefore a sensitivity
analysis, not a clean replacement for the original GC comparison.
"""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT, bh_fdr


DEFAULT_MODULES = ["grey60", "royalblue", "salmon", "pink", "green", "red", "midnightblue"]


def read_isa_file(metadata_source: Path, prefix: str, assay_contains: str | None = None) -> pd.DataFrame:
    """Read an ISA file from a directory, a direct file path, or an ISA zip."""
    metadata_source = Path(metadata_source)
    if metadata_source.is_file() and metadata_source.suffix == ".zip":
        with zipfile.ZipFile(metadata_source) as zf:
            names = zf.namelist()
            if prefix == "s_":
                candidates = [n for n in names if Path(n).name.startswith("s_") and n.endswith(".txt")]
            elif prefix == "a_":
                candidates = [n for n in names if Path(n).name.startswith("a_") and n.endswith(".txt")]
                if assay_contains:
                    narrowed = [n for n in candidates if assay_contains in Path(n).name]
                    candidates = narrowed or candidates
            else:
                raise ValueError(f"Unsupported ISA prefix: {prefix}")
            if not candidates:
                raise FileNotFoundError(f"No {prefix} ISA file found in {metadata_source}")
            with zf.open(sorted(candidates)[0]) as handle:
                return pd.read_csv(handle, sep="\t")

    if metadata_source.is_file():
        return pd.read_csv(metadata_source, sep="\t")

    if prefix == "s_":
        candidates = sorted(metadata_source.glob("s_*.txt"))
    elif prefix == "a_":
        candidates = sorted(metadata_source.glob("a_*.txt"))
        if assay_contains:
            narrowed = [p for p in candidates if assay_contains in p.name]
            candidates = narrowed or candidates
    else:
        raise ValueError(f"Unsupported ISA prefix: {prefix}")
    if not candidates:
        raise FileNotFoundError(f"No {prefix} ISA file found in {metadata_source}")
    return pd.read_csv(candidates[0], sep="\t")


def load_osd253_metadata(metadata_source: Path) -> pd.DataFrame:
    sample = read_isa_file(metadata_source, "s_")
    assay = read_isa_file(metadata_source, "a_", assay_contains="rna-seq")
    meta = sample.merge(assay, on="Sample Name", how="left", suffixes=("", "_assay"))
    required = [
        "Sample Name",
        "Factor Value[Strain]",
        "Factor Value[Spaceflight]",
        "Factor Value[Duration]",
        "Factor Value[Treatment]",
        "Parameter Value[Read Length]",
    ]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"OSD-253 metadata missing required columns: {missing}")
    return meta


def load_osd771_metadata(metadata_source: Path) -> pd.DataFrame:
    sample = read_isa_file(metadata_source, "s_")
    assay = read_isa_file(metadata_source, "a_", assay_contains="rna-seq")
    meta = sample.merge(assay, on="Sample Name", how="left", suffixes=("", "_assay"))
    required = [
        "Sample Name",
        "Characteristics[Strain]",
        "Factor Value[Spaceflight]",
        "Factor Value[Age]",
        "Factor Value[Euthanasia Location]",
        "Parameter Value[Recovery Duration on Earth]",
        "Parameter Value[duration]",
        "Parameter Value[Read Length]",
        "Parameter Value[Library Batch Number]",
    ]
    missing = [c for c in required if c not in meta.columns]
    if missing:
        raise ValueError(f"OSD-771 metadata missing required columns: {missing}")
    token = meta["Sample Name"].str.extract(r"_(ISS-T|LAR)_")[0]
    meta["Derived Arm"] = token.fillna("unknown")
    return meta


def load_vst(path: Path) -> pd.DataFrame:
    vst = pd.read_csv(path, index_col=0)
    vst.index = vst.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if vst.index.duplicated().any():
        vst = vst.groupby(vst.index).mean()
    return vst


def load_module_gene_lists(wgcna_dir: Path, modules: list[str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for module in modules:
        path = wgcna_dir / "module_gene_lists" / f"{module}.txt"
        if not path.exists():
            raise FileNotFoundError(f"Missing WGCNA module list: {path}")
        genes = [g.strip() for g in path.read_text().splitlines() if g.strip()]
        if len(genes) < 3:
            raise ValueError(f"Module {module} has too few genes: {len(genes)}")
        out[module] = genes
    return out


def module_scores(vst: pd.DataFrame, samples: list[str], modules: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = sorted(set(samples) - set(vst.columns))
    if missing:
        raise ValueError(f"Samples missing from VST matrix: {missing[:10]}")

    rows = {}
    coverage_rows = []
    available = set(vst.index.astype(str))
    for module, genes in modules.items():
        present = [g for g in genes if g in available]
        if len(present) < 3:
            continue
        expr = vst.loc[present, samples].to_numpy(dtype=float)
        mu = expr.mean(axis=1, keepdims=True)
        sd = expr.std(axis=1, ddof=1, keepdims=True)
        sd[sd < 1e-12] = 1.0
        rows[module] = ((expr - mu) / sd).mean(axis=0)
        coverage_rows.append({
            "module": module,
            "n_genes_total": len(genes),
            "n_genes_matched": len(present),
            "coverage": len(present) / len(genes),
        })
    return pd.DataFrame(rows, index=samples), pd.DataFrame(coverage_rows)


def select_osd253_scenario(meta: pd.DataFrame, scenario: str) -> pd.DataFrame:
    base = meta[
        meta["Factor Value[Spaceflight]"].isin(["Space Flight", "Ground Control", "Ground Control Rerun"])
        & meta["Factor Value[Duration]"].isin(["~25", "~75"])
        & meta["Factor Value[Strain]"].isin(["C57BL/6J", "C3H/HeJ"])
    ].copy()

    if scenario == "original_gc_blue":
        keep = (
            ((base["Factor Value[Spaceflight]"] == "Space Flight") & (base["Factor Value[Treatment]"] == "White light"))
            | ((base["Factor Value[Spaceflight]"] == "Ground Control") & (base["Factor Value[Treatment]"] == "Blue Light"))
        )
        label = {"Space Flight": "FLT", "Ground Control": "CTRL"}
        base = base[keep].copy()
        base["control_definition"] = "original_GC_blue_light_read98"
        base["condition"] = base["Factor Value[Spaceflight]"].map(label)
    elif scenario == "rerun_white_gc":
        keep = (
            ((base["Factor Value[Spaceflight]"] == "Space Flight") & (base["Factor Value[Treatment]"] == "White light"))
            | ((base["Factor Value[Spaceflight]"] == "Ground Control Rerun") & (base["Factor Value[Treatment]"] == "White light"))
        )
        label = {"Space Flight": "FLT", "Ground Control Rerun": "CTRL"}
        base = base[keep].copy()
        base["control_definition"] = "rerun_GC_white_light_read151"
        base["condition"] = base["Factor Value[Spaceflight]"].map(label)
    elif scenario == "rerun_blue_gc":
        keep = (
            ((base["Factor Value[Spaceflight]"] == "Space Flight") & (base["Factor Value[Treatment]"] == "White light"))
            | ((base["Factor Value[Spaceflight]"] == "Ground Control Rerun") & (base["Factor Value[Treatment]"] == "Blue Light"))
        )
        label = {"Space Flight": "FLT", "Ground Control Rerun": "CTRL"}
        base = base[keep].copy()
        base["control_definition"] = "rerun_GC_blue_light_read151"
        base["condition"] = base["Factor Value[Spaceflight]"].map(label)
    else:
        raise ValueError(f"Unknown scenario: {scenario}")

    base = base.dropna(subset=["condition"]).copy()
    return base


def average_strain_effect(df: pd.DataFrame, module: str, strain: str) -> float:
    diffs = []
    for duration, sub in df[df["Factor Value[Strain]"] == strain].groupby("Factor Value[Duration]"):
        flt = sub[sub["condition"] == "FLT"][module]
        ctrl = sub[sub["condition"] == "CTRL"][module]
        if len(flt) and len(ctrl):
            diffs.append(float(flt.mean() - ctrl.mean()))
    if not diffs:
        return float("nan")
    return float(np.mean(diffs))


def duration_contrast_rows(df: pd.DataFrame, module: str, scenario: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (strain, duration), sub in df.groupby(["Factor Value[Strain]", "Factor Value[Duration]"]):
        flt = sub[sub["condition"] == "FLT"][module]
        ctrl = sub[sub["condition"] == "CTRL"][module]
        if not len(flt) or not len(ctrl):
            continue
        rows.append({
            "scenario": scenario,
            "module": module,
            "strain": strain,
            "duration": duration,
            "mean_flt": float(flt.mean()),
            "mean_control": float(ctrl.mean()),
            "flight_minus_control": float(flt.mean() - ctrl.mean()),
            "n_flt": int(len(flt)),
            "n_control": int(len(ctrl)),
        })
    return rows


def observed_stats(df: pd.DataFrame, module: str) -> dict[str, float]:
    c57 = average_strain_effect(df, module, "C57BL/6J")
    c3h = average_strain_effect(df, module, "C3H/HeJ")
    return {
        "c57_effect": c57,
        "c3h_effect": c3h,
        "c3h_minus_c57": c3h - c57,
    }


def permute_conditions_within_strata(df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    permuted = df.copy()
    for _, idx in permuted.groupby(["Factor Value[Strain]", "Factor Value[Duration]"]).groups.items():
        labels = permuted.loc[idx, "condition"].to_numpy(copy=True)
        rng.shuffle(labels)
        permuted.loc[idx, "condition"] = labels
    return permuted


def permutation_summary(
    df: pd.DataFrame,
    module: str,
    scenario: str,
    n_perm: int,
    seed: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    obs = observed_stats(df, module)
    exceed_c57 = 0
    exceed_c3h = 0
    exceed_delta = 0
    for _ in range(n_perm):
        permuted = permute_conditions_within_strata(df, rng)
        stat = observed_stats(permuted, module)
        exceed_c57 += int(abs(stat["c57_effect"]) >= abs(obs["c57_effect"]))
        exceed_c3h += int(abs(stat["c3h_effect"]) >= abs(obs["c3h_effect"]))
        exceed_delta += int(abs(stat["c3h_minus_c57"]) >= abs(obs["c3h_minus_c57"]))

    row: dict[str, object] = {
        "scenario": scenario,
        "module": module,
        **obs,
        "abs_c57_effect": abs(obs["c57_effect"]),
        "abs_c3h_effect": abs(obs["c3h_effect"]),
        "same_direction": bool(np.sign(obs["c57_effect"]) == np.sign(obs["c3h_effect"])),
        "attenuation_fraction_vs_c57": np.nan,
        "p_c57_effect": (exceed_c57 + 1) / (n_perm + 1),
        "p_c3h_effect": (exceed_c3h + 1) / (n_perm + 1),
        "p_strain_delta": (exceed_delta + 1) / (n_perm + 1),
        "n_samples": int(len(df)),
        "n_permutations": int(n_perm),
    }
    if np.isfinite(obs["c57_effect"]) and abs(obs["c57_effect"]) > 1e-12:
        row["attenuation_fraction_vs_c57"] = 1.0 - (obs["c3h_effect"] / obs["c57_effect"])
    return row


def metadata_audits(osd253: pd.DataFrame, osd771: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    audit253 = (
        osd253.groupby([
            "Factor Value[Strain]",
            "Factor Value[Spaceflight]",
            "Factor Value[Duration]",
            "Factor Value[Treatment]",
            "Parameter Value[Read Length]",
        ], dropna=False)
        .size()
        .reset_index(name="n")
    )
    audit771 = (
        osd771.groupby([
            "Characteristics[Strain]",
            "Derived Arm",
            "Factor Value[Age]",
            "Factor Value[Spaceflight]",
            "Factor Value[Euthanasia Location]",
            "Parameter Value[Recovery Duration on Earth]",
            "Parameter Value[duration]",
            "Parameter Value[Read Length]",
            "Parameter Value[Library Batch Number]",
        ], dropna=False)
        .size()
        .reset_index(name="n")
    )
    return audit253, audit771


def run(args: argparse.Namespace) -> None:
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    osd253 = load_osd253_metadata(Path(args.osd253_metadata))
    osd771 = load_osd771_metadata(Path(args.osd771_metadata))
    audit253, audit771 = metadata_audits(osd253, osd771)
    audit253.to_csv(outdir / "osd253_metadata_audit.tsv", sep="\t", index=False)
    audit771.to_csv(outdir / "osd771_metadata_audit.tsv", sep="\t", index=False)

    vst = load_vst(Path(args.osd253_vst))
    sample_missing = sorted(set(osd253["Sample Name"]) - set(vst.columns))
    matrix_extra = sorted(set(vst.columns) - set(osd253["Sample Name"]))
    if sample_missing or matrix_extra:
        raise ValueError(
            "OSD-253 metadata/VST sample mismatch: "
            f"metadata_not_matrix={sample_missing[:10]}, matrix_not_metadata={matrix_extra[:10]}"
        )

    modules = [m.strip() for m in args.modules.split(",") if m.strip()]
    module_genes = load_module_gene_lists(Path(args.wgcna_dir), modules)

    summary_rows = []
    duration_rows = []
    coverage_tables = []
    scenario_counts = []
    for scenario in [s.strip() for s in args.scenarios.split(",") if s.strip()]:
        selected = select_osd253_scenario(osd253, scenario)
        samples = selected["Sample Name"].tolist()
        scores, coverage = module_scores(vst, samples, module_genes)
        coverage["scenario"] = scenario
        coverage_tables.append(coverage)

        selected = selected.merge(scores, left_on="Sample Name", right_index=True, how="left")
        counts = (
            selected.groupby([
                "control_definition",
                "Factor Value[Strain]",
                "Factor Value[Duration]",
                "condition",
                "Factor Value[Treatment]",
                "Parameter Value[Read Length]",
            ], dropna=False)
            .size()
            .reset_index(name="n")
        )
        counts.insert(0, "scenario", scenario)
        scenario_counts.append(counts)

        for i, module in enumerate(modules):
            duration_rows.extend(duration_contrast_rows(selected, module, scenario))
            summary_rows.append(
                permutation_summary(
                    selected,
                    module,
                    scenario,
                    n_perm=args.n_perm,
                    seed=args.seed + 1000 * i + len(summary_rows),
                )
            )

    summary = pd.DataFrame(summary_rows)
    for scenario, idx in summary.groupby("scenario").groups.items():
        loc = list(idx)
        summary.loc[loc, "q_c57_effect"] = bh_fdr(summary.loc[loc, "p_c57_effect"].to_numpy())
        summary.loc[loc, "q_c3h_effect"] = bh_fdr(summary.loc[loc, "p_c3h_effect"].to_numpy())
        summary.loc[loc, "q_strain_delta"] = bh_fdr(summary.loc[loc, "p_strain_delta"].to_numpy())

    duration = pd.DataFrame(duration_rows)
    coverage = pd.concat(coverage_tables, ignore_index=True)
    counts = pd.concat(scenario_counts, ignore_index=True)

    summary.to_csv(outdir / "osd253_module_projection_summary.tsv", sep="\t", index=False)
    duration.to_csv(outdir / "osd253_duration_contrasts.tsv", sep="\t", index=False)
    coverage.to_csv(outdir / "osd253_module_gene_coverage.tsv", sep="\t", index=False)
    counts.to_csv(outdir / "osd253_scenario_counts.tsv", sep="\t", index=False)

    manifest = {
        "analysis": "metadata-aware OSD-253 strain-stratified WGCNA module projection",
        "osd253_metadata": str(args.osd253_metadata),
        "osd771_metadata": str(args.osd771_metadata),
        "osd253_vst": str(args.osd253_vst),
        "wgcna_dir": str(args.wgcna_dir),
        "modules": modules,
        "scenarios": [s.strip() for s in args.scenarios.split(",") if s.strip()],
        "n_permutations": args.n_perm,
        "claim_boundary": (
            "Tests module-score flight effects and C3H/HeJ attenuation. "
            "original_gc_blue has light mismatch; rerun_white_gc fixes light but changes read length/run."
        ),
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--osd253_metadata",
        default=str(REPO_ROOT / "data/external/osdr/OSD-253/OSD-253_metadata_GLDS-253-ISA.zip"),
        help="OSD-253 ISA directory, zip, or sample file.",
    )
    ap.add_argument(
        "--osd771_metadata",
        default=str(REPO_ROOT / "data/raw/metadata"),
        help="OSD-771 ISA directory, zip, or sample file.",
    )
    ap.add_argument(
        "--osd253_vst",
        default=str(REPO_ROOT / "data/external/osdr/OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv"),
    )
    ap.add_argument(
        "--wgcna_dir",
        default=str(REPO_ROOT / "data/results/run_20260505_remediated_2500g/wgcna"),
    )
    ap.add_argument(
        "--outdir",
        default=str(REPO_ROOT / "data/results/run_20260517_213205_2500g/osd253_strain_module_projection"),
    )
    ap.add_argument("--modules", default=",".join(DEFAULT_MODULES))
    ap.add_argument("--scenarios", default="original_gc_blue,rerun_white_gc")
    ap.add_argument("--n_perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=771253)
    args = ap.parse_args()
    run(args)


if __name__ == "__main__":
    main()
