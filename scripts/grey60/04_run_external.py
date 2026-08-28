#!/usr/bin/env python3
"""Run cohort-level Grey60 recurrence and random-effects synthesis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import yaml

from scripts.osd253_strain_module_projection import (
    load_osd253_metadata,
    select_osd253_scenario,
)
from src.grey60.adversarial import (
    hedges_g,
    mean_z_score,
    osd462_animal_id,
    random_effects_reml_hk,
)


REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve(path: str) -> Path:
    p = Path(path)
    return p if p.is_absolute() else REPO_ROOT / p


def load_vst(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df


def assert_gene_coverage(
    expression: pd.DataFrame,
    genes: list[str],
    *,
    study: str,
) -> None:
    missing = sorted(set(genes) - set(expression.index))
    if missing or len(set(genes)) != len(genes):
        raise RuntimeError(
            f"{study} does not contain the complete unique frozen signature; "
            f"missing={missing}"
        )


def permutation_p(
    flight: np.ndarray,
    control: np.ndarray,
    *,
    n_perm: int,
    seed: int,
) -> float:
    values = np.concatenate([flight, control]).astype(float)
    nf = len(flight)
    observed = float(np.mean(flight) - np.mean(control))
    rng = np.random.default_rng(seed)
    exceed = 0
    for start in range(0, n_perm, 4096):
        b = min(4096, n_perm - start)
        keys = rng.random((b, len(values)))
        order = np.argsort(keys, axis=1)
        fmean = values[order[:, :nf]].mean(axis=1)
        cmean = values[order[:, nf:]].mean(axis=1)
        exceed += int(np.count_nonzero(np.abs(fmean - cmean) >= abs(observed)))
    return float((exceed + 1) / (n_perm + 1))


def simple_study(
    study: str,
    vst_path: Path,
    genes: list[str],
    *,
    orientation: float,
    n_perm: int,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    vst = load_vst(vst_path)
    assert_gene_coverage(vst, genes, study=study)
    flt = [c for c in vst.columns if "_FLT_" in c]
    gc = [c for c in vst.columns if "_GC_" in c]
    expected_n = {"OSD-102": 6, "OSD-163": 6, "OSD-513": 9}[study]
    if len(flt) != expected_n or len(gc) != expected_n:
        raise RuntimeError(
            f"{study} expected {expected_n}/{expected_n} FLT/GC samples; "
            f"found {len(flt)}/{len(gc)}"
        )
    samples = flt + gc
    score = mean_z_score(vst.loc[:, samples], genes) * orientation
    effect = hedges_g(score.loc[flt], score.loc[gc])
    p = permutation_p(
        score.loc[flt].to_numpy(),
        score.loc[gc].to_numpy(),
        n_perm=n_perm,
        seed=seed,
    )
    row = {
        "study": study,
        "context": "terminal" if study != "OSD-513" else "live_return",
        "estimate_g": effect.estimate,
        "variance_g": effect.variance,
        "se_g": effect.standard_error,
        "ci_low": effect.ci_low,
        "ci_high": effect.ci_high,
        "permutation_p": p,
        "n_flight": effect.n_flight,
        "n_control": effect.n_control,
        "n_genes": len([g for g in genes if g in vst.index]),
        "analysis": "frozen_48gene_mean_z",
    }
    scores = pd.DataFrame(
        {
            "study": study,
            "sample": samples,
            "condition": ["FLT"] * len(flt) + ["GC"] * len(gc),
            "score": score.loc[samples].to_numpy(),
        }
    )
    return row, scores


def fixed_effect_combine(rows: list[dict[str, object]]) -> tuple[float, float]:
    y = np.array([float(row["estimate_g"]) for row in rows])
    v = np.array([float(row["variance_g"]) for row in rows])
    w = 1.0 / v
    return float(np.sum(w * y) / np.sum(w)), float(1.0 / np.sum(w))


def osd253_study(
    vst_path: Path,
    isa_path: Path,
    genes: list[str],
    scenario: str,
    *,
    orientation: float,
    n_perm: int,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame]:
    vst = load_vst(vst_path)
    assert_gene_coverage(vst, genes, study="OSD-253")
    meta = load_osd253_metadata(isa_path)
    selected = select_osd253_scenario(meta, scenario)
    selected = selected[selected["Factor Value[Strain]"] == "C57BL/6J"].copy()
    reference_parts = [
        select_osd253_scenario(meta, reference_scenario)
        for reference_scenario in ("original_gc_blue", "rerun_white_gc")
    ]
    reference = pd.concat(reference_parts, ignore_index=True).drop_duplicates(
        subset=["Sample Name"]
    )
    reference = reference[
        reference["Factor Value[Strain]"] == "C57BL/6J"
    ].copy()
    reference_samples = [
        sample for sample in reference["Sample Name"] if sample in vst.columns
    ]
    if len(reference_samples) != reference["Sample Name"].nunique():
        raise RuntimeError("OSD-253 fixed score reference contains missing samples")
    fixed_score = mean_z_score(vst.loc[:, reference_samples], genes) * orientation
    samples = selected["Sample Name"].tolist()
    if any(sample not in fixed_score.index for sample in samples):
        raise RuntimeError("OSD-253 scenario contains samples outside fixed reference")
    score = fixed_score.loc[samples]
    selected["score"] = selected["Sample Name"].map(score)
    duration_rows = []
    for i, (duration, sub) in enumerate(
        selected.groupby("Factor Value[Duration]", sort=True)
    ):
        xf = sub.loc[sub["condition"] == "FLT", "score"].to_numpy(dtype=float)
        xc = sub.loc[sub["condition"] == "CTRL", "score"].to_numpy(dtype=float)
        g = hedges_g(xf, xc)
        duration_rows.append(
            {
                "study": "OSD-253",
                "scenario": scenario,
                "duration": duration,
                "estimate_g": g.estimate,
                "variance_g": g.variance,
                "ci_low": g.ci_low,
                "ci_high": g.ci_high,
                "permutation_p": permutation_p(
                    xf, xc, n_perm=n_perm, seed=seed + i
                ),
                "n_flight": len(xf),
                "n_control": len(xc),
            }
        )
    estimate, variance = fixed_effect_combine(duration_rows)
    se = float(np.sqrt(variance))
    row = {
        "study": "OSD-253",
        "context": "terminal",
        "estimate_g": estimate,
        "variance_g": variance,
        "se_g": se,
        "ci_low": estimate - 1.96 * se,
        "ci_high": estimate + 1.96 * se,
        "permutation_p": np.nan,
        "n_flight": int((selected["condition"] == "FLT").sum()),
        "n_control": int((selected["condition"] == "CTRL").sum()),
        "n_genes": len([g for g in genes if g in vst.index]),
        "analysis": f"C57_duration_fixed_effect_{scenario}",
    }
    scores = selected[
        [
            "Sample Name",
            "condition",
            "Factor Value[Duration]",
            "Factor Value[Treatment]",
            "score",
        ]
    ].copy()
    scores.insert(0, "study", "OSD-253")
    return row, scores, pd.DataFrame(duration_rows)


def osd462_study(
    study_dir: Path,
    genes: list[str],
    *,
    orientation: float,
    n_perm: int,
    seed: int,
) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prep_frames = []
    prep_rows = []
    for i, prep in enumerate(("UPX", "mRNA", "totRNA")):
        path = study_dir / f"GLDS-462_rna_seq_VST_Counts_{prep}_GLbulkRNAseq.csv"
        vst = load_vst(path)
        assert_gene_coverage(vst, genes, study=f"OSD-462/{prep}")
        flt = [c for c in vst.columns if "_WT_FLT_" in c]
        gc = [c for c in vst.columns if "_WT_GC_" in c]
        if len(flt) != 10 or len(gc) != 10:
            raise RuntimeError(
                f"OSD-462/{prep} expected 10/10 WT FLT/GC samples; "
                f"found {len(flt)}/{len(gc)}"
            )
        samples = flt + gc
        score = mean_z_score(vst.loc[:, samples], genes) * orientation
        # Each preparation gets equal scale before within-animal synthesis.
        score = (score - score.mean()) / score.std(ddof=1)
        xf = score.loc[flt].to_numpy(dtype=float)
        xc = score.loc[gc].to_numpy(dtype=float)
        g = hedges_g(xf, xc)
        prep_rows.append(
            {
                "study": "OSD-462",
                "preparation": prep,
                "mean_difference": float(xf.mean() - xc.mean()),
                "estimate_g": g.estimate,
                "variance_g": g.variance,
                "ci_low": g.ci_low,
                "ci_high": g.ci_high,
                "permutation_p": permutation_p(
                    xf, xc, n_perm=n_perm, seed=seed + i
                ),
                "n_flight": len(xf),
                "n_control": len(xc),
            }
        )
        prep_frames.append(
            pd.DataFrame(
                {
                    "animal": [osd462_animal_id(s) for s in samples],
                    "condition": ["FLT"] * len(flt) + ["GC"] * len(gc),
                    "preparation": prep,
                    "score": score.loc[samples].to_numpy(),
                }
            )
        )
    long = pd.concat(prep_frames, ignore_index=True)
    wide = long.pivot(index=["animal", "condition"], columns="preparation", values="score")
    if wide.isna().any().any() or len(wide) != 20:
        raise RuntimeError("OSD-462 preparations did not align to 20 unique FLT/GC animals")
    wide["score"] = wide[["UPX", "mRNA", "totRNA"]].mean(axis=1)
    wide = wide.reset_index()
    xf = wide.loc[wide["condition"] == "FLT", "score"].to_numpy(dtype=float)
    xc = wide.loc[wide["condition"] == "GC", "score"].to_numpy(dtype=float)
    g = hedges_g(xf, xc)
    row = {
        "study": "OSD-462",
        "context": "terminal",
        "estimate_g": g.estimate,
        "variance_g": g.variance,
        "se_g": g.standard_error,
        "ci_low": g.ci_low,
        "ci_high": g.ci_high,
        "permutation_p": permutation_p(
            xf, xc, n_perm=n_perm, seed=seed + 10
        ),
        "n_flight": len(xf),
        "n_control": len(xc),
        "n_genes": len(genes),
        "analysis": "within_animal_mean_across_UPX_mRNA_totRNA",
    }
    corr = wide[["UPX", "mRNA", "totRNA"]].corr(method="spearman")
    corr_long = (
        corr.rename_axis("preparation_1")
        .reset_index()
        .melt(id_vars="preparation_1", var_name="preparation_2", value_name="spearman_rho")
    )
    corr_long = corr_long[corr_long["preparation_1"] < corr_long["preparation_2"]]
    return row, wide, pd.DataFrame(prep_rows), corr_long


def meta_summary(
    rows: pd.DataFrame,
    *,
    label: str,
) -> tuple[dict[str, object], pd.DataFrame]:
    fit = random_effects_reml_hk(
        rows["estimate_g"], rows["variance_g"], modified=True
    )
    summary = {
        "analysis": label,
        "k": len(rows),
        "estimate": fit.estimate,
        "tau2_reml": fit.tau2,
        "se_hartung_knapp": fit.standard_error_hk,
        "ci_low_hartung_knapp": fit.ci_low_hk,
        "ci_high_hartung_knapp": fit.ci_high_hk,
        "p_hartung_knapp": fit.p_hk,
        "i_squared": fit.i_squared,
        "q": fit.q,
        "maximum_weight": float(fit.weights.max()),
    }
    weights = rows[["study"]].copy()
    weights["random_effect_weight"] = fit.weights
    return summary, weights


def dersimonian_laird_summary(
    rows: pd.DataFrame,
    *,
    label: str,
) -> dict[str, object]:
    y = rows["estimate_g"].to_numpy(dtype=float)
    v = rows["variance_g"].to_numpy(dtype=float)
    fixed_weights = 1.0 / v
    fixed_mean = float(np.sum(fixed_weights * y) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (y - fixed_mean) ** 2))
    df = len(y) - 1
    c = float(
        np.sum(fixed_weights)
        - np.sum(fixed_weights**2) / np.sum(fixed_weights)
    )
    tau2 = float(max(0.0, (q - df) / c))
    random_weights = 1.0 / (v + tau2)
    estimate = float(np.sum(random_weights * y) / np.sum(random_weights))
    se = float(np.sqrt(1.0 / np.sum(random_weights)))
    return {
        "analysis": label,
        "k": len(y),
        "estimate": estimate,
        "tau2": tau2,
        "standard_error": se,
        "ci_low": estimate - 1.96 * se,
        "ci_high": estimate + 1.96 * se,
        "p_normal": float(2 * stats.norm.sf(abs(estimate / se))),
        "q": q,
        "i_squared": float(max(0.0, (q - df) / q) * 100.0) if q > 0 else 0.0,
        "maximum_weight": float((random_weights / random_weights.sum()).max()),
    }


def run(args: argparse.Namespace) -> None:
    config_path = resolve(args.config)
    cfg = yaml.safe_load(config_path.read_text())
    inputs = {k: resolve(v) for k, v in cfg["inputs"].items()}
    root = resolve(args.outdir or cfg["output_dir"])
    outdir = root / "external"
    outdir.mkdir(parents=True, exist_ok=True)
    genes = [
        g.strip()
        for g in inputs["grey60_genes"].read_text().splitlines()
        if g.strip()
    ]
    orientation = float(
        json.loads((root / "internal" / "manifest.json").read_text())[
            "score_orientation"
        ]["mean_z"]["multiplier"]
    )
    n_perm = int(cfg["resampling"]["external_permutations"])
    seed = int(cfg["resampling"]["seed"]) + 500
    external_root = inputs["external_root"]

    rows = []
    score_tables = []
    simple_specs = {
        "OSD-102": external_root
        / "OSD-102/GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "OSD-163": external_root
        / "OSD-163/GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    }
    for i, (study, path) in enumerate(simple_specs.items()):
        row, scores = simple_study(
            study,
            path,
            genes,
            orientation=orientation,
            n_perm=n_perm,
            seed=seed + i,
        )
        rows.append(row)
        score_tables.append(scores)

    osd253_primary, scores253, durations_primary = osd253_study(
        external_root / "OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        external_root / "OSD-253/OSD-253_metadata_GLDS-253-ISA.zip",
        genes,
        cfg["external"]["osd253_primary_control"],
        orientation=orientation,
        n_perm=n_perm,
        seed=seed + 20,
    )
    rows.append(osd253_primary)
    score_tables.append(scores253.rename(columns={"Sample Name": "sample"}))
    osd253_sensitivity, _, durations_sensitivity = osd253_study(
        external_root / "OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        external_root / "OSD-253/OSD-253_metadata_GLDS-253-ISA.zip",
        genes,
        cfg["external"]["osd253_sensitivity_control"],
        orientation=orientation,
        n_perm=n_perm,
        seed=seed + 30,
    )
    osd253_sensitivity["context"] = "terminal_sensitivity_not_meta"

    row462, scores462, prep462, corr462 = osd462_study(
        external_root / "OSD-462",
        genes,
        orientation=orientation,
        n_perm=n_perm,
        seed=seed + 40,
    )
    rows.append(row462)
    scores462_out = scores462.copy()
    scores462_out.insert(0, "study", "OSD-462")
    score_tables.append(scores462_out)

    row513, scores513 = simple_study(
        "OSD-513",
        external_root / "OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        genes,
        orientation=orientation,
        n_perm=n_perm,
        seed=seed + 50,
    )
    rows.append(row513)
    score_tables.append(scores513)

    cohort = pd.DataFrame(rows)
    cohort = pd.concat(
        [cohort, pd.DataFrame([osd253_sensitivity])], ignore_index=True
    )
    cohort.to_csv(outdir / "cohort_effects.tsv", sep="\t", index=False)
    pd.concat(score_tables, ignore_index=True, sort=False).to_csv(
        outdir / "animal_scores.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    pd.concat(
        [
            durations_primary.assign(role="primary"),
            durations_sensitivity.assign(role="sensitivity"),
        ],
        ignore_index=True,
    ).to_csv(outdir / "osd253_duration_effects.tsv", sep="\t", index=False)
    prep462.to_csv(
        outdir / "osd462_preparation_effects.tsv", sep="\t", index=False
    )
    corr462.to_csv(
        outdir / "osd462_preparation_correlations.tsv", sep="\t", index=False
    )

    terminal = cohort[
        cohort["study"].isin(cfg["external"]["terminal_primary"])
        & (cohort["context"] == "terminal")
    ].copy()
    if len(terminal) != 4:
        raise RuntimeError(f"Expected four primary terminal studies; got {len(terminal)}")
    terminal_summary, weights = meta_summary(
        terminal, label="primary_terminal_REML_modified_Hartung_Knapp"
    )
    unmodified_hk = random_effects_reml_hk(
        terminal["estimate_g"], terminal["variance_g"], modified=False
    )
    sensitivity_summaries = [
        {
            "analysis": "sensitivity_terminal_REML_unmodified_Hartung_Knapp",
            "k": len(terminal),
            "estimate": unmodified_hk.estimate,
            "tau2": unmodified_hk.tau2,
            "standard_error": unmodified_hk.standard_error_hk,
            "ci_low": unmodified_hk.ci_low_hk,
            "ci_high": unmodified_hk.ci_high_hk,
            "p": unmodified_hk.p_hk,
            "q": unmodified_hk.q,
            "i_squared": unmodified_hk.i_squared,
            "maximum_weight": float(unmodified_hk.weights.max()),
        },
        dersimonian_laird_summary(
            terminal, label="sensitivity_terminal_Dersimonian_Laird"
        ),
    ]
    pd.DataFrame(sensitivity_summaries).to_csv(
        outdir / "terminal_meta_sensitivities.tsv", sep="\t", index=False
    )
    weights.to_csv(outdir / "terminal_study_weights.tsv", sep="\t", index=False)
    loo_rows = []
    for omitted in terminal["study"]:
        subset = terminal[terminal["study"] != omitted]
        fit, _ = meta_summary(subset, label=f"leave_out_{omitted}")
        fit["omitted_study"] = omitted
        loo_rows.append(fit)
    loo = pd.DataFrame(loo_rows)
    loo.to_csv(outdir / "leave_one_study_out.tsv", sep="\t", index=False)

    thresholds = cfg["go_no_go"]["gate_E"]
    positive_studies = int((terminal["estimate_g"] > 0).sum())
    components = {
        "positive_study_count": positive_studies
        >= thresholds["positive_terminal_studies_gte"],
        "pooled_ci": terminal_summary["ci_low_hartung_knapp"] > 0,
        "heterogeneity": terminal_summary["i_squared"]
        < thresholds["i_squared_lt"],
        "leave_one_study": bool((loo["estimate"] > 0).all()),
        "maximum_weight": terminal_summary["maximum_weight"]
        <= thresholds["maximum_study_weight_lte"],
    }
    gate = {
        "pass": all(components.values()),
        "components": components,
        "positive_terminal_studies": positive_studies,
        "terminal_meta": terminal_summary,
        "leave_one_study_minimum_estimate": float(loo["estimate"].min()),
        "osd513_live_return_effect": float(row513["estimate_g"]),
        "osd513_live_return_p": float(row513["permutation_p"]),
        "boundary": (
            "OSD-513 is a live-return moderator and is excluded from the "
            "primary terminal synthesis."
        ),
    }
    (outdir / "gate_E_status.json").write_text(
        json.dumps(gate, indent=2) + "\n"
    )
    pd.DataFrame([terminal_summary]).to_csv(
        outdir / "terminal_meta_summary.tsv", sep="\t", index=False
    )
    manifest = {
        "analysis_id": cfg["analysis_id"],
        "terminal_studies": cfg["external"]["terminal_primary"],
        "osd253_primary_control": cfg["external"]["osd253_primary_control"],
        "osd253_sensitivity_control": cfg["external"]["osd253_sensitivity_control"],
        "osd462_synthesis": cfg["external"]["osd462_primary_synthesis"],
        "recovery_moderator": cfg["external"]["recovery_moderator"],
        "effect_unit": "animal-level Hedges g",
        "meta_model": "REML with modified Hartung-Knapp uncertainty",
        "n_permutations_per_contrast": n_perm,
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(gate, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", default="config/grey60_adversarial_reanalysis.yaml"
    )
    parser.add_argument("--outdir", default="")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
