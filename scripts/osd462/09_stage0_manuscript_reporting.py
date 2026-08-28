#!/usr/bin/env python3
"""Create manuscript-ready displays from frozen OSD-462 Stage 0 outputs.

This script is deliberately a reporting layer. It reads the Stage 0 TSV
artifacts, validates their internal contracts, and changes no sample, feature,
or assay qualification.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import tempfile
from typing import Iterable

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(tempfile.gettempdir()) / "rrrm2-stage0-reporting-mpl-cache"),
)
os.environ.setdefault(
    "XDG_CACHE_HOME",
    str(Path(tempfile.gettempdir()) / "rrrm2-stage0-reporting-xdg-cache"),
)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
import pandas as pd  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STAGE0_DIR = (
    REPO_ROOT / "data" / "results" / "run_20260728_osd462_stage0"
)
DEFAULT_FIGURE_DIR = REPO_ROOT / "figures" / "v13"

DESIGN_FILE = "osd462_sample_plex_channel_design.tsv"
PROVENANCE_FILE = "osd462_ncc_spak_phosphosite_provenance.tsv"

CONDITION_ORDER = ("BL", "FL", "GC")
CONDITION_COLORS = {
    "BL": "#A7ADB5",
    "FL": "#D55E00",
    "GC": "#0072B2",
}
CONDITION_LABELS = {
    "BL": "BL — baseline context",
    "FL": "FL — primary contrast",
    "GC": "GC — primary contrast",
}
PLEX_LABELS = {
    "Samp1-5": "Plex 1 (Samp1–5)",
    "Samp6-10": "Plex 2 (Samp6–10)",
}


def _read_tsv(path: Path, required: Iterable[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required Stage 0 artifact is missing: {path}")
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing = sorted(set(required) - set(frame.columns))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return frame


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes", "y"}


def _single_value(group: pd.DataFrame, column: str, key: str) -> str:
    values = group[column].drop_duplicates().tolist()
    if len(values) != 1:
        raise ValueError(f"{key} has conflicting {column} values: {values}")
    return str(values[0])


def _ordered_design(design: pd.DataFrame) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Validate and order the identical protein/phosphoprotein channel maps."""
    expected_modalities = {"protein", "phosphoprotein"}
    if set(design["modality"]) != expected_modalities:
        raise ValueError(
            "Expected exactly protein and phosphoprotein modalities; found "
            f"{sorted(set(design['modality']))}"
        )
    if design.duplicated(["modality", "plex", "reporter_tag"]).any():
        raise ValueError("Duplicate modality/plex/reporter-tag rows detected.")

    modality_counts = design.groupby("modality", sort=False).size().to_dict()
    if modality_counts != {"phosphoprotein": 30, "protein": 30}:
        raise ValueError(f"Expected 30 rows per modality; found {modality_counts}")

    phospho = design.loc[design["modality"].eq("phosphoprotein")].copy()
    plexes = phospho["plex"].drop_duplicates().tolist()
    if len(plexes) != 2:
        raise ValueError(f"Expected two plexes; found {plexes}")

    first = phospho.loc[phospho["plex"].eq(plexes[0])].copy()
    first["_row"] = pd.to_numeric(first["workbook_intro_excel_row"], errors="raise")
    first = first.sort_values("_row", kind="stable")
    reporter_tags = first["reporter_tag"].tolist()
    expected_conditions = ["BL"] * 5 + ["FL"] * 5 + ["GC"] * 5
    if first["condition_code"].tolist() != expected_conditions:
        raise ValueError("First plex does not contain contiguous 5/5/5 BL/FL/GC blocks.")

    for modality in ("protein", "phosphoprotein"):
        for plex in plexes:
            current = design.loc[
                design["modality"].eq(modality) & design["plex"].eq(plex)
            ].copy()
            current["_row"] = pd.to_numeric(
                current["workbook_intro_excel_row"], errors="raise"
            )
            current = current.sort_values("_row", kind="stable")
            if current["reporter_tag"].tolist() != reporter_tags:
                raise ValueError(
                    f"Reporter-tag ordering differs for {modality}/{plex}."
                )
            if current["condition_code"].tolist() != expected_conditions:
                raise ValueError(
                    f"Condition-to-tag assignment differs for {modality}/{plex}."
                )

    return design.copy(), reporter_tags, plexes


def build_inclusion_table(design: pd.DataFrame) -> pd.DataFrame:
    """Collapse identical modality maps to one exact row per sample/channel."""
    required = [
        "modality",
        "plex",
        "workbook_sample_name",
        "isa_sample_name",
        "condition_code",
        "condition",
        "animal_number",
        "reporter_tag",
        "reporter_tag_matches_isa",
        "workbook_sample_matches_isa_source",
        "genotype",
        "strain",
        "sex",
        "tissue",
        "analysis_primary_flight_vs_ground",
        "raw_archive",
        "resolved_assay_name",
        "condition_perfectly_confounded_with_reporter_tag",
        "reporter_assignment_swapped_between_plexes",
        "workbook_intro_excel_row",
    ]
    missing = sorted(set(required) - set(design.columns))
    if missing:
        raise ValueError(f"Design table is missing columns: {missing}")

    design, reporter_tags, plexes = _ordered_design(design)
    tag_order = {tag: index + 1 for index, tag in enumerate(reporter_tags)}
    plex_order = {plex: index + 1 for index, plex in enumerate(plexes)}

    shared = [
        "workbook_sample_name",
        "isa_sample_name",
        "condition_code",
        "condition",
        "animal_number",
        "analysis_primary_flight_vs_ground",
        "genotype",
        "strain",
        "sex",
        "tissue",
        "condition_perfectly_confounded_with_reporter_tag",
        "reporter_assignment_swapped_between_plexes",
    ]
    records: list[dict[str, object]] = []
    for (plex, tag), group in design.groupby(
        ["plex", "reporter_tag"], sort=False, dropna=False
    ):
        key = f"{plex}/{tag}"
        if set(group["modality"]) != {"protein", "phosphoprotein"} or len(group) != 2:
            raise ValueError(f"{key} does not have one row in each modality.")
        values = {column: _single_value(group, column, key) for column in shared}
        primary = _truthy(values["analysis_primary_flight_vs_ground"])
        tag_match = group["reporter_tag_matches_isa"].map(_truthy).all()
        sample_match = group["workbook_sample_matches_isa_source"].map(_truthy).all()
        protein = group.loc[group["modality"].eq("protein")].iloc[0]
        phospho = group.loc[group["modality"].eq("phosphoprotein")].iloc[0]
        records.append(
            {
                "plex_order": plex_order[str(plex)],
                "plex": plex,
                "reporter_order": tag_order[str(tag)],
                "reporter_tag": tag,
                "condition_code": values["condition_code"],
                "condition": values["condition"],
                "workbook_sample_name": values["workbook_sample_name"],
                "isa_sample_name": values["isa_sample_name"],
                "animal_number": values["animal_number"],
                "primary_FL_vs_GC_included": primary,
                "inclusion_role": (
                    "included_primary_FL_vs_GC"
                    if primary
                    else "excluded_BL_baseline_context_only"
                ),
                "protein_included": True,
                "phosphoprotein_included": True,
                "reporter_tag_matches_ISA_both_modalities": bool(tag_match),
                "sample_matches_ISA_both_modalities": bool(sample_match),
                "genotype": values["genotype"],
                "strain": values["strain"],
                "sex": values["sex"],
                "tissue": values["tissue"],
                "protein_raw_archive": protein["raw_archive"],
                "phosphoprotein_raw_archive": phospho["raw_archive"],
                "resolved_assay_name": _single_value(
                    group, "resolved_assay_name", key
                ),
                "condition_tag_block_aliased": _truthy(
                    values["condition_perfectly_confounded_with_reporter_tag"]
                ),
                "cross_plex_label_swap": _truthy(
                    values["reporter_assignment_swapped_between_plexes"]
                ),
            }
        )

    out = pd.DataFrame.from_records(records).sort_values(
        ["plex_order", "reporter_order"], kind="stable"
    )
    if len(out) != 30:
        raise ValueError(f"Expected 30 unique sample/channel rows; found {len(out)}")
    if int(out["primary_FL_vs_GC_included"].sum()) != 20:
        raise ValueError("Expected exactly 20 samples in the primary FL-vs-GC contrast.")
    if not out["condition_tag_block_aliased"].all():
        raise ValueError("Stage 0 does not report condition/tag-block aliasing for all rows.")
    if out["cross_plex_label_swap"].any():
        raise ValueError("Stage 0 unexpectedly reports a cross-plex label swap.")
    if not out[
        [
            "protein_included",
            "phosphoprotein_included",
            "reporter_tag_matches_ISA_both_modalities",
            "sample_matches_ISA_both_modalities",
        ]
    ].all(axis=None):
        raise ValueError("One or more sample/channel inclusion or ISA checks failed.")
    return out.reset_index(drop=True)


def _qualification_reason(row: pd.Series) -> str:
    if _truthy(row["isolated_canonical_assay_feature"]):
        return "Yes — isolated strict-canonical feature"
    claim_class = row["canonical_claim_class"]
    if claim_class == "canonical_site_indexed_but_comodified_context_only":
        return (
            "No — canonical-indexed rollup; observed peptide carries "
            "additional modified residue(s)"
        )
    if claim_class == "mixed_composite_not_attributable_to_canonical_site":
        return "No — composite mixes strict-canonical and noncanonical residues"
    if claim_class == "not_in_strict_canonical_set":
        return "No — feature is outside the strict canonical NCC/SPAK site set"
    return f"No — {claim_class.replace('_', ' ')}"


def build_compact_phosphoform_table(provenance: pd.DataFrame) -> pd.DataFrame:
    """Select publication-facing provenance fields without reclassifying rows."""
    required = [
        "gene_symbol",
        "uniprot_accession",
        "source_sheet",
        "site_feature_kind",
        "source_excel_row",
        "workbook_site_position",
        "residue_resolved_site_labels",
        "workbook_max_score_components",
        "peptide_sequence_core",
        "sequence_phosphoform_site_labels",
        "all_localization_scores_ge_13",
        "isolated_canonical_assay_feature",
        "residue_indexed_canonical_but_comodified",
        "n_quantified_peptides_Samp1-5",
        "n_quantified_peptides_Samp6-10",
        "strict_canonical_components",
        "noncanonical_components",
        "canonical_claim_class",
        "n_positive_scaled_BL",
        "n_positive_scaled_FL",
        "n_positive_scaled_GC",
    ]
    missing = sorted(set(required) - set(provenance.columns))
    if missing:
        raise ValueError(f"Provenance table is missing columns: {missing}")

    genes = set(provenance["gene_symbol"])
    if genes != {"Slc12a3", "Stk39"}:
        raise ValueError(f"Expected only Slc12a3 and Stk39; found {sorted(genes)}")

    out = pd.DataFrame(
        {
            "gene_symbol": provenance["gene_symbol"],
            "uniprot_accession": provenance["uniprot_accession"],
            "feature_type": provenance["site_feature_kind"],
            "source_sheet": provenance["source_sheet"],
            "source_excel_row": provenance["source_excel_row"],
            "workbook_position": provenance["workbook_site_position"],
            "residue_resolved_feature": provenance["residue_resolved_site_labels"],
            "observed_peptide_phosphoform": provenance[
                "sequence_phosphoform_site_labels"
            ],
            "modified_peptide": provenance["peptide_sequence_core"],
            "AScore_components": provenance["workbook_max_score_components"],
            "all_AScores_ge_13": provenance["all_localization_scores_ge_13"].map(
                _truthy
            ),
            "quantified_peptides_plex1": provenance[
                "n_quantified_peptides_Samp1-5"
            ],
            "quantified_peptides_plex2": provenance[
                "n_quantified_peptides_Samp6-10"
            ],
            "positive_reporters_BL_FL_GC": provenance.apply(
                lambda row: (
                    f"{row['n_positive_scaled_BL']}/10;"
                    f"{row['n_positive_scaled_FL']}/10;"
                    f"{row['n_positive_scaled_GC']}/10"
                ),
                axis=1,
            ),
            "strict_canonical_components": provenance["strict_canonical_components"],
            "noncanonical_components": provenance["noncanonical_components"],
            "canonical_indexed_but_comodified": provenance[
                "residue_indexed_canonical_but_comodified"
            ].map(_truthy),
            "isolated_canonical_qualified": provenance[
                "isolated_canonical_assay_feature"
            ].map(lambda value: "Yes" if _truthy(value) else "No"),
            "qualification_reason": provenance.apply(_qualification_reason, axis=1),
        }
    )
    out["_gene_order"] = out["gene_symbol"].map({"Slc12a3": 0, "Stk39": 1})
    out["_feature_order"] = out["feature_type"].map(
        {"single_site_rollup": 0, "composite_site_rollup": 1}
    )
    out["_position_order"] = pd.to_numeric(
        out["workbook_position"].str.split(";").str[0], errors="coerce"
    )
    out = out.sort_values(
        ["_gene_order", "_feature_order", "_position_order", "workbook_position"],
        kind="stable",
    ).drop(columns=["_gene_order", "_feature_order", "_position_order"])

    if len(out) != len(provenance):
        raise ValueError("Compact phosphoform reporting dropped provenance rows.")
    if out["isolated_canonical_qualified"].eq("Yes").any():
        raise ValueError(
            "Stage 0 reporting expected zero isolated canonical features but found one."
        )
    return out.reset_index(drop=True)


def plot_reporter_tag_map(
    inclusion: pd.DataFrame,
    figure_dir: Path,
    *,
    dpi: int = 600,
) -> list[Path]:
    """Plot the exact, repeated BL/FL/GC reporter-tag assignment."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    ordered = inclusion.sort_values(
        ["plex_order", "reporter_order"], kind="stable"
    ).copy()
    tags = (
        ordered.loc[ordered["plex_order"].eq(1)]
        .sort_values("reporter_order")["reporter_tag"]
        .tolist()
    )
    plexes = (
        ordered.sort_values("plex_order")["plex"].drop_duplicates().tolist()
    )
    if len(tags) != 15 or len(plexes) != 2:
        raise ValueError("Reporter-tag figure requires 15 tags in each of two plexes.")

    matplotlib.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "svg.hashsalt": "osd462-stage0-reporter-map",
        }
    )
    fig, ax = plt.subplots(figsize=(12.8, 4.6))
    y_for_plex = {plexes[0]: 1.0, plexes[1]: 0.0}

    for _, row in ordered.iterrows():
        x = int(row["reporter_order"]) - 1
        y = y_for_plex[str(row["plex"])]
        condition = str(row["condition_code"])
        ax.add_patch(
            Rectangle(
                (x - 0.46, y - 0.37),
                0.92,
                0.74,
                facecolor=CONDITION_COLORS[condition],
                edgecolor="white",
                linewidth=1.25,
            )
        )
        sample_label = str(row["workbook_sample_name"]).replace("RR-10_", "")
        text_color = "#202124" if condition == "BL" else "white"
        ax.text(
            x,
            y,
            sample_label,
            ha="center",
            va="center",
            fontsize=7.2,
            fontweight="semibold",
            color=text_color,
        )

    block_ranges = {"BL": (0, 4), "FL": (5, 9), "GC": (10, 14)}
    for condition in CONDITION_ORDER:
        start, end = block_ranges[condition]
        ax.plot(
            [start - 0.40, end + 0.40],
            [1.70, 1.70],
            color=CONDITION_COLORS[condition],
            linewidth=4.0,
            solid_capstyle="round",
            clip_on=False,
        )
        ax.text(
            (start + end) / 2,
            1.85,
            CONDITION_LABELS[condition],
            ha="center",
            va="bottom",
            fontsize=9.2,
            fontweight="semibold",
            color="#222222",
        )

    ax.axvline(4.5, color="#565A5E", linestyle=(0, (2, 2)), linewidth=0.8)
    ax.axvline(9.5, color="#565A5E", linestyle=(0, (2, 2)), linewidth=0.8)
    ax.set_xlim(-0.65, 14.65)
    ax.set_ylim(-0.65, 2.15)
    ax.set_xticks(range(15), tags)
    ax.tick_params(axis="x", length=0, pad=4, labelsize=8.6)
    ax.set_xlabel("TMTpro reporter tag", labelpad=8, fontweight="semibold")
    ax.set_yticks(
        [y_for_plex[plex] for plex in plexes],
        [PLEX_LABELS.get(plex, plex) for plex in plexes],
    )
    ax.tick_params(axis="y", length=0, pad=8, labelsize=9.2)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.suptitle(
        "OSD-462 TMTpro reporter-tag assignment",
        x=0.08,
        y=0.985,
        ha="left",
        fontsize=15,
        fontweight="bold",
        color="#17191B",
    )
    fig.text(
        0.08,
        0.925,
        "Identical BL/FL/GC tag blocks in both plexes; no cross-plex label swap",
        ha="left",
        va="top",
        fontsize=10.5,
        color="#4B5055",
    )
    fig.text(
        0.08,
        0.025,
        (
            "Protein and phosphoprotein workbooks use the same channel map. "
            "Condition is perfectly aliased with reporter-tag block; the "
            "FL–GC contrast cannot separate condition from a systematic "
            "tag-block effect."
        ),
        ha="left",
        va="bottom",
        fontsize=8.4,
        color="#3F4448",
    )
    fig.subplots_adjust(left=0.16, right=0.985, top=0.79, bottom=0.23)

    base = figure_dir / "osd462_reporter_tag_map"
    outputs = [base.with_suffix(".png"), base.with_suffix(".pdf"), base.with_suffix(".svg")]
    fig.savefig(
        outputs[0],
        dpi=dpi,
        facecolor="white",
        bbox_inches="tight",
        metadata={
            "Title": "OSD-462 TMTpro reporter-tag assignment",
            "Description": "Exact BL, FL, and GC tag blocks in both OSD-462 plexes.",
        },
    )
    fig.savefig(
        outputs[1],
        facecolor="white",
        bbox_inches="tight",
        metadata={
            "Title": "OSD-462 TMTpro reporter-tag assignment",
            "Subject": "Exact sample-to-reporter-tag design and absence of label swap",
            "Creator": "09_stage0_manuscript_reporting.py",
            "CreationDate": None,
            "ModDate": None,
        },
    )
    fig.savefig(
        outputs[2],
        facecolor="white",
        bbox_inches="tight",
        metadata={
            "Title": "OSD-462 TMTpro reporter-tag assignment",
            "Description": "Exact BL, FL, and GC tag blocks in both OSD-462 plexes.",
            "Date": None,
        },
    )
    plt.close(fig)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Render manuscript-ready sample/channel and phosphoform provenance "
            "outputs from the frozen OSD-462 Stage 0 audit."
        )
    )
    parser.add_argument(
        "--stage0-dir",
        type=Path,
        default=DEFAULT_STAGE0_DIR,
        help=f"Stage 0 input directory (default: {DEFAULT_STAGE0_DIR})",
    )
    parser.add_argument(
        "--reporting-dir",
        type=Path,
        default=None,
        help="Table output directory (default: <stage0-dir>/reporting)",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=DEFAULT_FIGURE_DIR,
        help=f"Figure output directory (default: {DEFAULT_FIGURE_DIR})",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="PNG resolution (default: 600)",
    )
    args = parser.parse_args()

    stage0_dir = args.stage0_dir.resolve()
    reporting_dir = (
        args.reporting_dir.resolve()
        if args.reporting_dir is not None
        else stage0_dir / "reporting"
    )
    figure_dir = args.figure_dir.resolve()
    reporting_dir.mkdir(parents=True, exist_ok=True)

    design = _read_tsv(
        stage0_dir / DESIGN_FILE,
        required=[
            "modality",
            "plex",
            "reporter_tag",
            "condition_code",
            "workbook_sample_name",
        ],
    )
    provenance = _read_tsv(
        stage0_dir / PROVENANCE_FILE,
        required=[
            "gene_symbol",
            "residue_resolved_site_labels",
            "isolated_canonical_assay_feature",
        ],
    )

    inclusion = build_inclusion_table(design)
    compact_provenance = build_compact_phosphoform_table(provenance)

    inclusion_path = (
        reporting_dir / "osd462_exact_sample_channel_inclusion.tsv"
    )
    provenance_path = (
        reporting_dir
        / "osd462_slc12a3_stk39_phosphoform_provenance_compact.tsv"
    )
    inclusion.to_csv(inclusion_path, sep="\t", index=False)
    compact_provenance.to_csv(provenance_path, sep="\t", index=False)
    figure_paths = plot_reporter_tag_map(inclusion, figure_dir, dpi=args.dpi)

    print(f"Wrote {inclusion_path} ({len(inclusion)} rows; 20 primary samples)")
    print(
        f"Wrote {provenance_path} ({len(compact_provenance)} rows; "
        "0 isolated canonical features)"
    )
    for path in figure_paths:
        print(f"Wrote {path}")
    print(
        "Verified: identical protein/phosphoprotein maps; contiguous BL/FL/GC "
        "tag blocks in both plexes; no cross-plex label swap."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
