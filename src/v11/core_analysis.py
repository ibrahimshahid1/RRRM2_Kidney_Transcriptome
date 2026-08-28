#!/usr/bin/env python3
"""Execute the core v11 DCT-subtype-prior phosphoproteome analyses."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
import yaml


RUN_ROOT = Path("data/results/run_20260526_v11_dct1_phospho_mediation")
OSD462 = Path("data/results/run_20260522_osd462_anchor/osd462_anchor")
PHENO = Path("data/results/run_20260522_phenotype_anchor")
CELLTYPE = Path("data/results/run_20260522_celltype_decomposition")
REGULATOR = Path("data/results/run_20260522_regulator_activity")
DCT_PRIOR_DIR = RUN_ROOT / "dct_prior"
PXD_DIR = Path("data/external/phosphoproteomics/PXD001729")
KSEA_NET = Path("data/external/kinase_substrate/renal_kinase_substrate_core.tsv")
MECHANISM_SETS = Path("config/mechanism_gene_sets.yaml")
OSD462_PROTEIN_XLSX = Path("data/external/osdr/OSD-462/GLDS-462_proteomics_2021-12-31_tc884-885_Protein_WorkUp.xlsx")
OSD462_PHOSPHO_XLSX = Path("data/external/osdr/OSD-462/GLDS-462_phosphproteomics_2021-12-31_tc882-883_Pho_WorkUp_JM.xlsx")


ANCHOR_GENES = {"Slc12a3", "Stk39", "Oxsr1", "Wnk1", "Wnk4"}
TRANSPORT_TARGETS = {
    "Slc12a3",
    "Stk39",
    "Oxsr1",
    "Wnk1",
    "Wnk4",
    "Klhl3",
    "Cul3",
    "Nedd4l",
    "Sgk1",
    "Kcnj10",
    "Kcnj16",
    "Trpm6",
    "Pvalb",
    "Calb1",
}

PATHWAY_FAMILIES = {
    "ecm_organization": "matrix_remodeling",
    "fibrosis_tgfb_emt": "matrix_remodeling",
    "integrin_cell_adhesion": "matrix_remodeling",
    "mmp_adam_proteolysis": "matrix_remodeling",
    "dct_ncc_wnk_transport": "distal_tubule_transport",
    "tubular_transport_broad": "distal_tubule_transport",
    "tlr4_innate": "immune_inflammatory",
    "macrophage_inflammation": "immune_inflammatory",
    "oxidative_stress_nrf2": "stress_response",
    "preservation_stress_response": "stress_response",
    "s1p_s1pr3": "s1p_lipid_signaling",
}

FISHER_ALTERNATIVE = "greater"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(pvals: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(pvals), dtype=float)
    out = np.ones_like(p)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    vals = p[ok]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(vals)
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[ok] = restored
    return out


def safe_float(s):
    return pd.to_numeric(s, errors="coerce")


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    if keep.sum() == 0:
        return float("nan")
    a = a[keep]
    b = b[keep]
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return float("nan")
    return float(np.dot(a, b) / denom)


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def ensure_dirs(root: Path):
    for name in [
        "baseline",
        "cross_layer",
        "external_qc",
        "dct_prior",
        "h2_enrichment",
        "h2_pxd",
        "h2_klhl3",
        "h3_mediation",
        "manifests",
    ]:
        (root / name).mkdir(parents=True, exist_ok=True)


def baseline_lock(root: Path):
    summary = json.loads((OSD462 / "results_summary.json").read_text())
    pheno_summary = read_tsv(PHENO / "phenotype_anchor_summary.tsv")
    cell_corr = read_tsv(CELLTYPE / "osd462_compartment_vs_ncc_phospho.tsv")
    ksea = read_tsv(REGULATOR / "osd462_kinase_activity_summary.tsv")

    rows = [
        {
            "component": "OSD462_RNA_recurrence",
            "metric": "pathway_cosine",
            "value": summary["layer4_rna_gate"]["point_cosine"],
            "status": "PASS" if summary["layer4_rna_gate"]["recurrence_pass"] else "FAIL",
            "interpretation": "OSD-462 RNA recurs RRRM-2 ISS-T matrix-high/DCT-low direction",
        },
        {
            "component": "OSD462_protein",
            "metric": "any_target_set_concordant",
            "value": summary["layer1_protein"]["any_set_concordant_in_predicted_direction"],
            "status": "NULL",
            "interpretation": "targeted protein abundance concordance remains null",
        },
        {
            "component": "OSD462_phospho",
            "metric": "isolated_canonical_ncc_spak_coverage",
            "value": 0,
            "status": "UNRESOLVED",
            "interpretation": (
                "Stage 0 found zero isolated canonical assay features; "
                "T53/Y65 and S382/S383 rows are co-modified context only"
            ),
        },
    ]
    for _, row in ksea.iterrows():
        rows.append(
            {
                "component": "KSEA",
                "metric": row["kinase"],
                "value": row["z_score"],
                "status": row["direction"],
                "interpretation": (
                    "kinase activity is not inferred without isolated, "
                    "phosphoform-qualified substrate features"
                ),
            }
        )
    for _, row in pheno_summary.iterrows():
        rows.append(
            {
                "component": "phenotype_anchor",
                "metric": row["comparison"],
                "value": row["spearman_condition_adjusted"],
                "status": row["interpretation"],
                "interpretation": "animal-matched RNA-phospho link remains suggestive/underpowered except non-regulatory wrinkle",
            }
        )
    for _, row in cell_corr.iterrows():
        rows.append(
            {
                "component": "celltype_vs_ncc_phospho",
                "metric": row["panel"],
                "value": np.nan,
                "status": "INVALIDATED_AS_ACTIVITY_EVIDENCE",
                "interpretation": (
                    "historical score used co-modified position-indexed features; "
                    "it is not NCC regulatory-site activity"
                ),
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(root / "baseline" / "v11_baseline_lock_summary.tsv", sep="\t", index=False)
    write_ksea_substrate_table(root)
    write_rna_recurrence_supplements(root)
    write_cross_layer_pathway_matrix(root)
    write_tmt_qc_summaries(root)

    inputs = [
        Path("docs/v11_execution_research_plan.md"),
        MECHANISM_SETS,
        OSD462_PROTEIN_XLSX,
        OSD462_PHOSPHO_XLSX,
        OSD462 / "results_summary.json",
        OSD462 / "osd462_rna_pathway_effects.tsv",
        PHENO / "phenotype_anchor_summary.tsv",
        PHENO / "phenotype_anchor_per_animal.tsv",
        CELLTYPE / "osd462_compartment_scores_per_sample.tsv",
        CELLTYPE / "osd462_compartment_vs_ncc_phospho.tsv",
        REGULATOR / "osd462_kinase_activity_summary.tsv",
    ]
    manifest = {
        "analysis": "v11 baseline lock",
        "run_root": str(root),
        "inputs": [{"path": str(p), "sha256": sha256(p)} for p in inputs if p.exists()],
    }
    (root / "baseline" / "v11_baseline_input_manifest.json").write_text(json.dumps(manifest, indent=2))


def write_ksea_substrate_table(root: Path) -> pd.DataFrame:
    """Write the curated substrate rows behind the targeted KSEA summary."""
    if not KSEA_NET.exists():
        return pd.DataFrame()
    net = read_tsv(KSEA_NET)
    phospho = read_tsv(OSD462 / "phospho_all_sites.tsv")
    phospho["site_position_str"] = phospho["site_position"].astype(str)
    phospho["phospho_q_value_all_sites"] = bh(phospho["phospho_p_value"])
    rows = []
    for _, sub in net.iterrows():
        site = str(sub["substrate_site"]).strip()
        matched = phospho[
            phospho["gene_symbol"].astype(str).eq(str(sub["substrate_gene"]).strip())
            & phospho["site_position_str"].eq(site)
        ].copy()
        if matched.empty:
            rows.append(
                {
                    "kinase": sub["kinase"],
                    "substrate_gene": sub["substrate_gene"],
                    "substrate_site": site,
                    "matched_in_osd462": False,
                    "phospho_effect": np.nan,
                    "phospho_p_value": np.nan,
                    "phospho_q_value_all_sites": np.nan,
                    "n_fl": 0,
                    "n_gc": 0,
                    "n_finite_channels": 0,
                    "evidence": sub.get("evidence", ""),
                }
            )
            continue
        for _, row in matched.iterrows():
            n_fl = pd.to_numeric(row.get("n_fl"), errors="coerce")
            n_gc = pd.to_numeric(row.get("n_gc"), errors="coerce")
            rows.append(
                {
                    "kinase": sub["kinase"],
                    "substrate_gene": sub["substrate_gene"],
                    "substrate_site": site,
                    "matched_in_osd462": True,
                    "phospho_effect": row["phospho_effect"],
                    "phospho_p_value": row["phospho_p_value"],
                    "phospho_q_value_all_sites": row["phospho_q_value_all_sites"],
                    "n_fl": n_fl,
                    "n_gc": n_gc,
                    "n_finite_channels": n_fl + n_gc if np.isfinite(n_fl) and np.isfinite(n_gc) else np.nan,
                    "evidence": sub.get("evidence", ""),
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(root / "baseline" / "osd462_targeted_ksea_substrates.tsv", sep="\t", index=False)
    return out


def write_rna_recurrence_supplements(root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Write transparent gene-set membership and pathway-family sensitivities."""
    if not MECHANISM_SETS.exists():
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    raw = yaml.safe_load(MECHANISM_SETS.read_text()) or {}
    member_rows = []
    for pathway, spec in raw.items():
        genes = list(spec.get("genes") or [])
        protected = {str(g).upper() for g in spec.get("protected_genes") or []}
        family = PATHWAY_FAMILIES.get(pathway, str(spec.get("role", "unassigned")))
        for i, gene in enumerate(genes, start=1):
            member_rows.append(
                {
                    "pathway": pathway,
                    "pathway_family": family,
                    "role": spec.get("role", ""),
                    "description": spec.get("description", ""),
                    "gene_symbol": gene,
                    "gene_rank_in_set": i,
                    "protected_gene": str(gene).upper() in protected,
                }
            )
    members = pd.DataFrame(member_rows)
    members.to_csv(root / "baseline" / "rna_recurrence_gene_set_members.tsv", sep="\t", index=False)

    effects_path = OSD462 / "osd462_rna_pathway_effects.tsv"
    if not effects_path.exists():
        return members, pd.DataFrame(), pd.DataFrame()
    effects = read_tsv(effects_path)
    effects["pathway_family"] = effects["pathway"].map(PATHWAY_FAMILIES).fillna("unassigned")
    full_cosine = cosine(effects["osd462_rna_pathway_effect"], effects["rrrm2_iss_t_pathway_effect"])
    family_rows = []
    for family, part in effects.groupby("pathway_family"):
        keep = effects["pathway_family"].ne(family)
        kept = effects[keep]
        family_rows.append(
            {
                "dropped_pathway_family": family,
                "dropped_pathways": ",".join(part["pathway"].astype(str)),
                "n_pathways_retained": int(len(kept)),
                "full_cosine": full_cosine,
                "cosine_after_drop": cosine(
                    kept["osd462_rna_pathway_effect"],
                    kept["rrrm2_iss_t_pathway_effect"],
                ),
                "boundary": "Pathway-family leave-out over the curated OSD-462/RRRM-2 anchor panel.",
            }
        )
    families = pd.DataFrame(family_rows)
    if not families.empty:
        families["delta_vs_full"] = families["cosine_after_drop"] - families["full_cosine"]
    families.to_csv(root / "baseline" / "osd462_rna_recurrence_leave_one_family.tsv", sep="\t", index=False)

    paired = effects[["osd462_rna_pathway_effect", "rrrm2_iss_t_pathway_effect"]].dropna()
    boot_rows = []
    if len(paired):
        rng = np.random.default_rng(46211)
        x = paired["osd462_rna_pathway_effect"].to_numpy(dtype=float)
        y = paired["rrrm2_iss_t_pathway_effect"].to_numpy(dtype=float)
        for i in range(2000):
            idx = rng.integers(0, len(paired), len(paired))
            boot_rows.append(
                {
                    "bootstrap": i + 1,
                    "n_pathways": int(len(paired)),
                    "cosine": cosine(x[idx], y[idx]),
                    "bootstrap_unit": "paired_curated_pathway",
                }
            )
    boot = pd.DataFrame(boot_rows)
    boot.to_csv(root / "baseline" / "osd462_rna_recurrence_paired_pathway_bootstrap.tsv", sep="\t", index=False)
    return members, families, boot


def _direction(value: float, threshold: float) -> str:
    if not np.isfinite(value) or abs(value) < threshold:
        return "flat"
    return "up" if value > 0 else "down"


def write_cross_layer_pathway_matrix(root: Path) -> pd.DataFrame:
    """Summarize where OSD-462 pathway signals live across RNA, protein, and phosphosite layers."""
    members_path = root / "baseline" / "rna_recurrence_gene_set_members.tsv"
    rna_path = OSD462 / "osd462_rna_pathway_effects.tsv"
    protein_path = OSD462 / "protein_effects_gene_anyplex.tsv"
    phospho_path = OSD462 / "phospho_all_sites.tsv"
    if not (members_path.exists() and rna_path.exists() and protein_path.exists() and phospho_path.exists()):
        return pd.DataFrame()

    members = read_tsv(members_path)
    rna = read_tsv(rna_path)
    proteins = read_tsv(protein_path)
    phospho = read_tsv(phospho_path)

    members = members.copy()
    proteins = proteins.copy()
    phospho = phospho.copy()
    members["gene_upper"] = members["gene_symbol"].astype(str).str.upper()
    proteins["gene_upper"] = proteins["gene_symbol"].astype(str).str.upper()
    phospho["gene_upper"] = phospho["gene_symbol"].astype(str).str.upper()

    protein_effects = proteins.groupby("gene_upper", as_index=False).agg(
        protein_effect=("flight_effect", "mean"),
    )
    phospho_effects = phospho.groupby("gene_upper", as_index=False).agg(
        phospho_effect=("phospho_effect", "mean"),
        n_phospho=("phospho_effect", "count"),
    )

    rows = []
    for pathway, part in members.groupby("pathway", sort=False):
        genes = part["gene_upper"].dropna().astype(str).drop_duplicates()
        rna_row = rna.loc[rna["pathway"].eq(pathway)]
        p = protein_effects.loc[protein_effects["gene_upper"].isin(genes)]
        ph = phospho_effects.loc[phospho_effects["gene_upper"].isin(genes)]
        rows.append(
            {
                "pathway": pathway,
                "n_set": int(len(genes)),
                "rna_effect": float(rna_row["osd462_rna_pathway_effect"].iloc[0]) if not rna_row.empty else np.nan,
                "protein_effect": float(p["protein_effect"].mean()) if len(p) else np.nan,
                "n_protein": int(len(p)),
                "phospho_effect": float(ph["phospho_effect"].mean()) if len(ph) else np.nan,
                "n_phospho": int(len(ph)),
            }
        )

    out = pd.DataFrame(rows)
    for col, zcol in [
        ("rna_effect", "rna_z"),
        ("protein_effect", "protein_z"),
        ("phospho_effect", "phospho_z"),
    ]:
        x = pd.to_numeric(out[col], errors="coerce")
        sd = x.std(ddof=0)
        out[zcol] = (x - x.mean()) / sd if np.isfinite(sd) and sd > 0 else np.nan

    out["rna_dir"] = out["rna_effect"].apply(lambda x: _direction(x, 0.04))
    out["protein_dir"] = out["protein_effect"].apply(lambda x: _direction(x, 0.02))
    out["phospho_dir"] = out["phospho_effect"].apply(lambda x: _direction(x, 0.02))

    def pattern(row: pd.Series) -> str:
        if row["rna_dir"] in {"up", "down"} and row["protein_dir"] in {"up", "down"}:
            if row["rna_dir"] != row["protein_dir"]:
                return "RNA-protein DISCORDANT (opposite)"
            return "RNA-protein concordant"
        if row["rna_dir"] in {"up", "down"} and row["protein_dir"] == "flat":
            return "protein-flat (RNA not propagated)"
        return "RNA-protein concordant"

    out["pattern"] = out.apply(pattern, axis=1)
    out = out.sort_values("rna_effect", ascending=False)
    out.to_csv(root / "cross_layer" / "osd462_cross_layer_pathway_matrix.tsv", sep="\t", index=False)

    summary = (
        out.groupby("pattern", as_index=False)
        .agg(
            n_pathways=("pathway", "size"),
            pathways=("pathway", lambda x: ",".join(x.astype(str))),
        )
        .sort_values("n_pathways", ascending=False)
    )
    summary["n_total_pathways"] = int(len(out))
    summary["fraction_pathways"] = summary["n_pathways"] / summary["n_total_pathways"]
    summary.to_csv(root / "cross_layer" / "osd462_cross_layer_pattern_summary.tsv", sep="\t", index=False)
    return out


def _tmt_qc_from_sheet(xlsx: Path, sheet: str, layer: str) -> pd.DataFrame:
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    if raw.shape[0] < 4:
        return pd.DataFrame()
    section = raw.iloc[0].ffill().astype(str)
    labels = raw.iloc[1].astype(str)
    headers = raw.iloc[2].astype(str)
    data = raw.iloc[3:].reset_index(drop=True)
    rows = []
    for col_idx in range(raw.shape[1]):
        sample = labels.iloc[col_idx].strip()
        if not re.match(r"^(BL|FL|GC)-\d+", sample):
            continue
        values = pd.to_numeric(data.iloc[:, col_idx], errors="coerce")
        header = headers.iloc[col_idx]
        section_label = section.iloc[col_idx].lower()
        if "scaled" in section_label:
            metric = "scaled_signal_to_noise_row_sum_100"
        elif "signal-to-noise" in section_label or "_sn_" in header:
            metric = "summed_signal_to_noise"
        else:
            metric = section_label.replace(" ", "_")
        plex = "Samp1-5" if str(header).startswith("Samp1-5") else "Samp6-10" if str(header).startswith("Samp6-10") else ""
        finite = values[np.isfinite(values)]
        rows.append(
            {
                "layer": layer,
                "source_workbook": str(xlsx),
                "sheet": sheet,
                "metric": metric,
                "plex": plex,
                "sample_label": sample,
                "condition": sample.split("-", 1)[0],
                "n_rows": int(len(values)),
                "n_finite": int(len(finite)),
                "missing_fraction": float(1 - len(finite) / len(values)) if len(values) else np.nan,
                "median": float(finite.median()) if len(finite) else np.nan,
                "mean": float(finite.mean()) if len(finite) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def write_tmt_qc_summaries(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize raw OSD-462 TMT channel medians and missingness by layer."""
    frames = []
    if OSD462_PROTEIN_XLSX.exists():
        frames.append(_tmt_qc_from_sheet(OSD462_PROTEIN_XLSX, "protein_quant_2721", "protein"))
    if OSD462_PHOSPHO_XLSX.exists():
        frames.append(_tmt_qc_from_sheet(OSD462_PHOSPHO_XLSX, "siteQuant_360", "phosphosite_single"))
        frames.append(_tmt_qc_from_sheet(OSD462_PHOSPHO_XLSX, "siteQuant_360_compositeSite", "phosphosite_composite"))
    qc = pd.concat([x for x in frames if not x.empty], ignore_index=True) if frames else pd.DataFrame()
    qc.to_csv(root / "baseline" / "osd462_tmt_channel_qc.tsv", sep="\t", index=False)
    if qc.empty:
        summary = pd.DataFrame()
    else:
        summary = (
            qc.groupby(["layer", "metric", "plex", "condition"], as_index=False)
            .agg(
                n_channels=("sample_label", "nunique"),
                median_channel_median=("median", "median"),
                median_missing_fraction=("missing_fraction", "median"),
                max_missing_fraction=("missing_fraction", "max"),
            )
            .sort_values(["layer", "metric", "plex", "condition"])
        )
    summary.to_csv(root / "baseline" / "osd462_tmt_missingness_by_condition.tsv", sep="\t", index=False)
    return qc, summary


def parse_site(site) -> list[tuple[int, str, str]]:
    if pd.isna(site):
        return []
    text = str(site)
    found = []
    for match in re.finditer(r"(\d+)_([STY])([*?]?)", text):
        found.append((int(match.group(1)), match.group(2), match.group(3) or ""))
    return found


def load_pxd_tables(root: Path):
    xls = PXD_DIR / "41598_2015_BFsrep12829_MOESM6_ESM.xls"
    inventory_rows = []
    for sheet in pd.ExcelFile(xls, engine="xlrd").sheet_names:
        df = pd.read_excel(xls, sheet_name=sheet, engine="xlrd", header=None)
        inventory_rows.append(
            {
                "file": xls.name,
                "sheet": sheet,
                "n_rows_raw": len(df),
                "n_cols_raw": df.shape[1],
            }
        )
    pd.DataFrame(inventory_rows).to_csv(root / "external_qc" / "pxd001729_table_inventory.tsv", sep="\t", index=False)

    total = pd.read_excel(xls, sheet_name="Total phosphosites", engine="xlrd", header=13)
    total = total.rename(columns=lambda c: str(c).strip())
    total["ddavp_effect_log2"] = np.nan
    total["p_neg_log10"] = np.nan
    total["ddavp_direction"] = "identified_total"

    changed_frames = []
    for sheet, direction in [
        ("phos_sites_increased_with_dDAVP", "increased"),
        ("phos_sites_decreased_with_dDAVP", "decreased"),
    ]:
        df = pd.read_excel(xls, sheet_name=sheet, engine="xlrd", header=2)
        df = df.rename(columns=lambda c: str(c).strip())
        df["ddavp_direction"] = direction
        changed_frames.append(df)
    changed = pd.concat(changed_frames, ignore_index=True)

    def normalize(df: pd.DataFrame) -> pd.DataFrame:
        colmap = {
            "Peptide Sequence": "peptide_sequence",
            "Phospho. Sites": "phospho_sites",
            "Phospho. sites": "phospho_sites",
            "Accession No.": "protein_accession",
            "Protein Accession": "protein_accession",
            "Protein Name": "protein_name",
            "Gene Symbol": "gene_symbol",
            "Mean (log2)": "ddavp_effect_log2",
            "p (-log10)": "p_neg_log10",
            "Site Conservation Score": "site_conservation_score",
            "Motif Conservation Score": "motif_conservation_score",
        }
        keep = [c for c in df.columns if c in colmap]
        out = df[keep].rename(columns=colmap).copy()
        for c in ["gene_symbol", "phospho_sites", "peptide_sequence", "protein_accession", "protein_name"]:
            if c in out.columns:
                out[c] = out[c].astype(str).str.strip()
        if "ddavp_effect_log2" not in out.columns:
            out["ddavp_effect_log2"] = np.nan
        if "p_neg_log10" not in out.columns:
            out["p_neg_log10"] = np.nan
        return out

    total_n = normalize(total)
    total_n["source_sheet"] = "Total phosphosites"
    total_n["ddavp_direction"] = "identified_total"
    changed_n = normalize(changed)
    changed_n["source_sheet"] = changed["ddavp_direction"].map(
        {
            "increased": "phos_sites_increased_with_dDAVP",
            "decreased": "phos_sites_decreased_with_dDAVP",
        }
    )
    changed_n["ddavp_direction"] = changed["ddavp_direction"].to_numpy()

    rows = []
    for _, row in pd.concat([total_n, changed_n], ignore_index=True).iterrows():
        for pos, residue, confidence in parse_site(row.get("phospho_sites")):
            d = row.to_dict()
            d["site_position"] = pos
            d["site_residue"] = residue
            d["site_assignment"] = confidence
            d["site_id"] = f"{d.get('gene_symbol')}:{residue}{pos}"
            rows.append(d)
    parsed = pd.DataFrame(rows)
    parsed["ddavp_effect_log2"] = safe_float(parsed["ddavp_effect_log2"])
    parsed["p_neg_log10"] = safe_float(parsed["p_neg_log10"])
    parsed["p_value"] = 10 ** (-parsed["p_neg_log10"])
    parsed.loc[parsed["p_neg_log10"].isna(), "p_value"] = np.nan
    parsed.to_csv(root / "external_qc" / "pxd001729_phosphosite_effects.tsv", sep="\t", index=False)

    target = parsed[parsed["gene_symbol"].isin(TRANSPORT_TARGETS)].copy()
    cov = (
        target.groupby("gene_symbol", dropna=False)
        .agg(
            n_sites=("site_id", "count"),
            n_changed=("ddavp_effect_log2", lambda s: int(s.notna().sum())),
            n_increased=("ddavp_direction", lambda s: int((s == "increased").sum())),
            n_decreased=("ddavp_direction", lambda s: int((s == "decreased").sum())),
        )
        .reset_index()
    )
    cov.to_csv(root / "external_qc" / "pxd001729_target_site_coverage.tsv", sep="\t", index=False)
    changed_n.to_csv(root / "h2_pxd" / "pxd001729_ddavp_direction.tsv", sep="\t", index=False)
    target.to_csv(root / "h2_pxd" / "pxd001729_dct_target_direction.tsv", sep="\t", index=False)
    return parsed


def build_dct_prior_mapping(root: Path):
    de = read_tsv(DCT_PRIOR_DIR / "gse228367_dct1_vs_dct2_de.tsv")
    expressed = (
        (de["pct_detected_dct1"] >= 0.05)
        | (de["pct_detected_dct2"] >= 0.05)
        | (de["mean_expr_dct1"] >= 0.05)
        | (de["mean_expr_dct2"] >= 0.05)
    )
    q75 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.75)
    q25 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.25)
    q90 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.90)
    q10 = de.loc[expressed, "dct1_enrichment_score"].quantile(0.10)
    de["dct1_top_quartile"] = expressed & (de["dct1_enrichment_score"] >= q75)
    de["dct2_bottom_quartile"] = expressed & (de["dct1_enrichment_score"] <= q25)
    de["dct1_top_decile"] = expressed & (de["dct1_enrichment_score"] >= q90)
    de["dct2_bottom_decile"] = expressed & (de["dct1_enrichment_score"] <= q10)
    de["dct1_core_fdr"] = de["dct_expression_class"].eq("DCT1_core")
    de["dct2_core_fdr"] = de["dct_expression_class"].eq("DCT2_core")
    de["dct1_score_percentile"] = np.nan
    de.loc[expressed, "dct1_score_percentile"] = de.loc[expressed, "dct1_enrichment_score"].rank(pct=True)
    de["dct2_leaning_percentile"] = 1 - de["dct1_score_percentile"]

    proteins = read_tsv(OSD462 / "protein_effects_gene_anyplex.tsv")
    phospho = read_tsv(OSD462 / "phospho_all_sites.tsv")
    phospho["phospho_q_value_all_sites"] = bh(phospho["phospho_p_value"])
    phospho["site_position_str"] = phospho["site_position"].astype(str)
    phospho["is_single_site"] = phospho["site_position_str"].str.fullmatch(r"\d+")
    phospho["site_position_int"] = pd.to_numeric(phospho["site_position_str"], errors="coerce")
    phospho["site_id"] = phospho["gene_symbol"].astype(str) + ":" + phospho["site_position_str"]
    phospho["is_anchor_gene"] = phospho["gene_symbol"].isin(ANCHOR_GENES)
    phospho["is_ncc_site"] = phospho["gene_symbol"].eq("Slc12a3")
    phospho["is_suppressed_p05"] = (phospho["phospho_effect"] < 0) & (phospho["phospho_p_value"] < 0.05)
    phospho["is_suppressed_q10"] = (phospho["phospho_effect"] < 0) & (phospho["phospho_q_value_all_sites"] < 0.10)
    effect_q25 = pd.to_numeric(phospho["phospho_effect"], errors="coerce").quantile(0.25)
    phospho["is_effect_bottom_quartile"] = pd.to_numeric(phospho["phospho_effect"], errors="coerce") <= effect_q25

    prior_cols = [
        "gene_symbol",
        "mean_expr_dct1",
        "pct_detected_dct1",
        "mean_expr_dct2",
        "pct_detected_dct2",
        "log2_mean_ratio_dct1_vs_dct2",
        "dct1_enrichment_score",
        "p_value",
        "fdr",
        "dct_expression_class",
        "dct1_top_quartile",
        "dct2_bottom_quartile",
        "dct1_top_decile",
        "dct2_bottom_decile",
        "dct1_core_fdr",
        "dct2_core_fdr",
        "dct1_score_percentile",
        "dct2_leaning_percentile",
    ]
    proteins_prior = proteins.merge(de[prior_cols], on="gene_symbol", how="left", indicator="dct_prior_merge")
    phospho_prior = phospho.merge(de[prior_cols], on="gene_symbol", how="left", indicator="dct_prior_merge")
    phospho_prior = phospho_prior.merge(
        proteins[["gene_symbol", "flight_effect", "n_peptides", "abundance_log2", "plex_coverage"]],
        on="gene_symbol",
        how="left",
        suffixes=("", "_protein"),
    )
    phospho_prior["abundance_bin"] = pd.qcut(phospho_prior["abundance_log2"].rank(method="first"), 4, labels=False, duplicates="drop")
    phospho_prior["peptide_bin"] = pd.qcut(phospho_prior["n_peptides"].rank(method="first"), 4, labels=False, duplicates="drop")
    phospho_prior["match_stratum"] = phospho_prior["abundance_bin"].astype("Int64").astype(str) + "_p" + phospho_prior[
        "peptide_bin"
    ].astype("Int64").astype(str)

    proteins_prior.to_csv(root / "dct_prior" / "osd462_protein_dct1_prior.tsv", sep="\t", index=False)
    phospho_prior.to_csv(root / "dct_prior" / "osd462_phosphosite_dct1_prior.tsv", sep="\t", index=False)
    de.to_csv(root / "dct_prior" / "dct1_enrichment_prior_v1.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(
        [
            {
                "table": "osd462_proteins",
                "n_rows": len(proteins_prior),
                "n_mapped": int(proteins_prior["dct1_enrichment_score"].notna().sum()),
                "fraction_mapped": proteins_prior["dct1_enrichment_score"].notna().mean(),
            },
            {
                "table": "osd462_phosphosites",
                "n_rows": len(phospho_prior),
                "n_mapped": int(phospho_prior["dct1_enrichment_score"].notna().sum()),
                "fraction_mapped": phospho_prior["dct1_enrichment_score"].notna().mean(),
            },
        ]
    )
    target_cov = phospho_prior[phospho_prior["gene_symbol"].isin(TRANSPORT_TARGETS)].groupby("gene_symbol").agg(
        n_phosphosites=("site_id", "count"),
        mapped=("dct1_enrichment_score", lambda s: bool(s.notna().any())),
        dct1_enrichment_score=("dct1_enrichment_score", "first"),
        dct1_score_percentile=("dct1_score_percentile", "first"),
        dct2_leaning_percentile=("dct2_leaning_percentile", "first"),
        dct_expression_class=("dct_expression_class", "first"),
        dct1_top_quartile=("dct1_top_quartile", "first"),
        dct1_top_decile=("dct1_top_decile", "first"),
        dct2_bottom_decile=("dct2_bottom_decile", "first"),
    )
    target_cov = target_cov.reset_index()
    coverage.to_csv(root / "dct_prior" / "osd462_dct1_prior_coverage.tsv", sep="\t", index=False)
    target_cov.to_csv(root / "dct_prior" / "osd462_dct1_prior_target_gene_coverage.tsv", sep="\t", index=False)
    return de, proteins_prior, phospho_prior


def zscore_values(x: pd.Series) -> pd.Series:
    vals = pd.to_numeric(x, errors="coerce")
    sd = vals.std(skipna=True, ddof=0)
    if not np.isfinite(sd) or sd <= 0:
        return vals * np.nan
    return (vals - vals.mean(skipna=True)) / sd


def run_dct_prior_sensitivity(root: Path, phospho_prior: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build leave-one-replicate-out and score-definition DCT prior checks."""
    by_rep_path = root / "dct_prior" / "gse228367_gene_stats_by_rep.tsv"
    if not by_rep_path.exists():
        return pd.DataFrame(), pd.DataFrame()
    stats_by_rep = read_tsv(by_rep_path)
    reps = sorted(stats_by_rep["rep"].dropna().astype(str).unique())
    variants: list[dict[str, object]] = [
        {"variant": "all_reps", "held_out_rep": "none", "keep_reps": reps},
    ]
    variants.extend({"variant": f"leave_one_rep_{r}", "held_out_rep": r, "keep_reps": [x for x in reps if x != r]} for r in reps)

    rows = []
    for variant in variants:
        sub = stats_by_rep[stats_by_rep["rep"].astype(str).isin(variant["keep_reps"])].copy()
        wide = (
            sub.groupby(["gene_symbol", "subtype"], as_index=False)
            .agg(mean_expr=("mean_expr", "mean"), pct_detected=("pct_detected", "mean"))
            .pivot(index="gene_symbol", columns="subtype")
        )
        wide.columns = [f"{metric}_{subtype.lower()}" for metric, subtype in wide.columns]
        wide = wide.reset_index().fillna(0)
        eps = 0.01
        wide["mean_difference"] = wide["mean_expr_dct1"] - wide["mean_expr_dct2"]
        wide["log2_ratio"] = np.log2((wide["mean_expr_dct1"] + eps) / (wide["mean_expr_dct2"] + eps))
        wide["rank_average"] = (
            wide["mean_expr_dct1"].rank(pct=True) - wide["mean_expr_dct2"].rank(pct=True)
        )
        wide["detection_aware"] = (
            zscore_values(wide["mean_difference"]).fillna(0)
            + zscore_values(wide["pct_detected_dct1"] - wide["pct_detected_dct2"]).fillna(0)
        ) / 2
        expressed = (
            (wide["pct_detected_dct1"] >= 0.05)
            | (wide["pct_detected_dct2"] >= 0.05)
            | (wide["mean_expr_dct1"] >= 0.05)
            | (wide["mean_expr_dct2"] >= 0.05)
        )
        for score_def in ["mean_difference", "log2_ratio", "rank_average", "detection_aware"]:
            score = pd.to_numeric(wide[score_def], errors="coerce")
            q90 = score.loc[expressed].quantile(0.90)
            q10 = score.loc[expressed].quantile(0.10)
            pct = pd.Series(np.nan, index=wide.index, dtype=float)
            pct.loc[expressed] = score.loc[expressed].rank(pct=True)
            frame = pd.DataFrame(
                {
                    "prior_variant": variant["variant"],
                    "held_out_rep": variant["held_out_rep"],
                    "score_definition": score_def,
                    "gene_symbol": wide["gene_symbol"],
                    "score": score,
                    "score_percentile": pct,
                    "dct1_top_decile": expressed & (score >= q90),
                    "dct2_bottom_decile": expressed & (score <= q10),
                    "mean_expr_dct1": wide["mean_expr_dct1"],
                    "mean_expr_dct2": wide["mean_expr_dct2"],
                    "pct_detected_dct1": wide["pct_detected_dct1"],
                    "pct_detected_dct2": wide["pct_detected_dct2"],
                }
            )
            rows.append(frame)
    sens = pd.concat(rows, ignore_index=True)
    sens.to_csv(root / "dct_prior" / "gse228367_dct_prior_score_sensitivity.tsv", sep="\t", index=False)

    anchor_genes = sorted(ANCHOR_GENES | {"Pvalb", "Trpm6", "Trpv5", "Calb1", "Klhl3"})
    anchor = sens[sens["gene_symbol"].isin(anchor_genes)].copy()
    anchor.to_csv(root / "dct_prior" / "gse228367_anchor_gene_prior_positions.tsv", sep="\t", index=False)

    parent_genes = set(phospho_prior["gene_symbol"].dropna().astype(str))
    primary = sens[
        sens["prior_variant"].eq("all_reps")
        & sens["score_definition"].eq("mean_difference")
    ][["gene_symbol", "dct1_top_decile", "dct2_bottom_decile"]].copy()
    primary = primary[primary["gene_symbol"].isin(parent_genes)]
    primary_top = set(primary.loc[primary["dct1_top_decile"], "gene_symbol"])
    primary_bottom = set(primary.loc[primary["dct2_bottom_decile"], "gene_symbol"])
    stability_rows = []
    for (variant, score_def), g in sens[sens["gene_symbol"].isin(parent_genes)].groupby(["prior_variant", "score_definition"]):
        top = set(g.loc[g["dct1_top_decile"], "gene_symbol"])
        bottom = set(g.loc[g["dct2_bottom_decile"], "gene_symbol"])
        stability_rows.append(
            {
                "prior_variant": variant,
                "score_definition": score_def,
                "n_osd462_parent_genes_scored": int(g["gene_symbol"].nunique()),
                "n_dct1_top_decile_parent_genes": len(top),
                "n_dct2_bottom_decile_parent_genes": len(bottom),
                "jaccard_vs_primary_dct1_top_decile": len(top & primary_top) / len(top | primary_top) if (top | primary_top) else np.nan,
                "jaccard_vs_primary_dct2_bottom_decile": len(bottom & primary_bottom) / len(bottom | primary_bottom) if (bottom | primary_bottom) else np.nan,
            }
        )
    stability = pd.DataFrame(stability_rows)
    stability.to_csv(root / "dct_prior" / "gse228367_prior_bin_stability.tsv", sep="\t", index=False)
    return sens, anchor


def fisher_table(df: pd.DataFrame, flag: str, suppressed_col: str):
    sub = df[df["dct1_enrichment_score"].notna()].copy()
    tab = pd.crosstab(sub[suppressed_col], sub[flag].fillna(False))
    for r in [False, True]:
        if r not in tab.index:
            tab.loc[r] = 0
    for c in [False, True]:
        if c not in tab.columns:
            tab[c] = 0
    tab = tab.sort_index().sort_index(axis=1)
    arr = np.array([[tab.loc[True, True], tab.loc[True, False]], [tab.loc[False, True], tab.loc[False, False]]])
    odds, p = stats.fisher_exact(arr, alternative=FISHER_ALTERNATIVE)
    frac_sup = arr[0, 0] / max(arr[0].sum(), 1)
    frac_bg = arr[1, 0] / max(arr[1].sum(), 1)
    fold = frac_sup / frac_bg if frac_bg > 0 else np.inf
    return odds, p, fold, arr


def one_representative_site_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    """Select one phosphosite row per parent gene for row-dependence sensitivity.

    The representative row is the most statistically responsive site on the
    parent gene, using phosphosite p value as the primary key and more negative
    flight effect as the tie-breaker. This is distinct from
    ``is_single_site``, which only excludes composite/multi-position rows.
    """
    cols = ["gene_symbol", "phospho_p_value", "phospho_effect", "site_id"]
    sub = df.dropna(subset=["gene_symbol"]).copy()
    return (
        sub.sort_values(cols, ascending=[True, True, True, True])
        .groupby("gene_symbol", as_index=False)
        .head(1)
        .copy()
    )


def one_single_position_representative_site_per_gene(df: pd.DataFrame) -> pd.DataFrame:
    """Select one representative row per gene after excluding composite rows."""
    if "is_single_site" not in df.columns:
        return one_representative_site_per_gene(df)
    single = df[df["is_single_site"].astype(bool)].copy()
    return one_representative_site_per_gene(single)


def build_parent_gene_table(phospho_prior: pd.DataFrame) -> pd.DataFrame:
    """Collapse phosphosite rows to one parent-gene row."""
    sub = phospho_prior[phospho_prior["dct1_enrichment_score"].notna()].copy()
    rows = []
    for gene, g in sub.groupby("gene_symbol", sort=False):
        n_fl = pd.to_numeric(g["n_fl"], errors="coerce") if "n_fl" in g.columns else pd.Series(0, index=g.index)
        n_gc = pd.to_numeric(g["n_gc"], errors="coerce") if "n_gc" in g.columns else pd.Series(0, index=g.index)
        n_valid = n_fl + n_gc
        rows.append(
            {
                "gene_symbol": gene,
                "any_suppressed_p05": bool(g["is_suppressed_p05"].astype(bool).any()),
                "any_suppressed_q10": bool(g["is_suppressed_q10"].astype(bool).any()),
                "any_effect_bottom_quartile": bool(g.get("is_effect_bottom_quartile", pd.Series(False, index=g.index)).astype(bool).any()),
                "min_phospho_p_value": pd.to_numeric(g["phospho_p_value"], errors="coerce").min(),
                "most_negative_phospho_effect": pd.to_numeric(g["phospho_effect"], errors="coerce").min(),
                "mean_phospho_effect": pd.to_numeric(g["phospho_effect"], errors="coerce").mean(),
                "n_quantified_phosphosites": int(len(g)),
                "n_single_position_sites": int(g["is_single_site"].astype(bool).sum()),
                "n_composite_sites": int((~g["is_single_site"].astype(bool)).sum()),
                "mean_n_valid_samples": n_valid.mean(),
                "mean_missing_samples": 20 - n_valid.mean(),
                "dct1_enrichment_score": g["dct1_enrichment_score"].iloc[0],
                "dct1_score_percentile": g.get("dct1_score_percentile", pd.Series(np.nan, index=g.index)).iloc[0],
                "dct2_leaning_percentile": g.get("dct2_leaning_percentile", pd.Series(np.nan, index=g.index)).iloc[0],
                "dct1_top_quartile": bool(g["dct1_top_quartile"].fillna(False).iloc[0]),
                "dct1_top_decile": bool(g["dct1_top_decile"].fillna(False).iloc[0]),
                "dct2_bottom_quartile": bool(g["dct2_bottom_quartile"].fillna(False).iloc[0]),
                "dct2_bottom_decile": bool(g["dct2_bottom_decile"].fillna(False).iloc[0]),
                "dct1_core_fdr": bool(g["dct1_core_fdr"].fillna(False).iloc[0]),
                "dct2_core_fdr": bool(g["dct2_core_fdr"].fillna(False).iloc[0]),
                "parent_protein_effect": pd.to_numeric(g["flight_effect"], errors="coerce").iloc[0],
                "parent_n_peptides": pd.to_numeric(g["n_peptides"], errors="coerce").iloc[0],
                "parent_abundance_log2": pd.to_numeric(g["abundance_log2"], errors="coerce").iloc[0],
            }
        )
    out = pd.DataFrame(rows)
    if len(out):
        for col in ["n_quantified_phosphosites", "parent_n_peptides", "parent_abundance_log2", "mean_missing_samples"]:
            x = pd.to_numeric(out[col], errors="coerce")
            if col.startswith("n_") or col == "parent_n_peptides":
                x = np.log1p(x)
            sd = x.std(skipna=True, ddof=0)
            out[f"{col}_z"] = (x - x.mean(skipna=True)) / sd if np.isfinite(sd) and sd > 0 else np.nan
    return out


def parent_gene_logistic(parent: pd.DataFrame, flag: str, suppressed_col: str) -> dict:
    """Fit a simple parent-gene-level logistic sensitivity model."""
    try:
        import statsmodels.api as sm
    except Exception as exc:  # pragma: no cover - depends on optional dependency
        return {"model_status": f"not_fit: {exc}"}

    covars = [flag, "n_quantified_phosphosites_z", "parent_n_peptides_z", "parent_abundance_log2_z"]
    d = parent.dropna(subset=[suppressed_col] + covars).copy()
    if len(d) < 20 or d[suppressed_col].nunique() < 2 or d[flag].nunique() < 2:
        return {"model_status": "not_fit: insufficient variation"}
    X = sm.add_constant(d[covars].astype(float), has_constant="add")
    y = d[suppressed_col].astype(float)
    try:
        fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
        coef = float(fit.params[flag])
        se = float(fit.bse[flag])
        p = float(fit.pvalues[flag])
        return {
            "model_status": "fit",
            "log_odds_coef": coef,
            "se": se,
            "odds_ratio": float(np.exp(coef)),
            "ci_low": float(np.exp(coef - 1.96 * se)),
            "ci_high": float(np.exp(coef + 1.96 * se)),
            "p_value": p,
            "n_parent_genes_model": int(len(d)),
        }
    except Exception as exc:
        return {"model_status": f"not_fit: {exc}"}


def parent_gene_joint_logistic(parent: pd.DataFrame, suppressed_col: str) -> list[dict]:
    """Fit DCT1-top and DCT2-bottom decile flags in the same parent-gene model."""
    try:
        import statsmodels.api as sm
    except Exception as exc:  # pragma: no cover - depends on optional dependency
        return [{"model_status": f"not_fit: {exc}"}]

    terms = ["dct1_top_decile", "dct2_bottom_decile"]
    covars = terms + [
        "n_quantified_phosphosites_z",
        "parent_n_peptides_z",
        "parent_abundance_log2_z",
        "mean_missing_samples_z",
    ]
    d = parent.dropna(subset=[suppressed_col] + covars).copy()
    if len(d) < 20 or d[suppressed_col].nunique() < 2:
        return [{"model_status": "not_fit: insufficient variation"}]
    X = sm.add_constant(d[covars].astype(float), has_constant="add")
    y = d[suppressed_col].astype(float)
    try:
        fit = sm.GLM(y, X, family=sm.families.Binomial()).fit()
    except Exception as exc:
        return [{"model_status": f"not_fit: {exc}"}]
    rows = []
    for term in terms:
        coef = float(fit.params[term])
        se = float(fit.bse[term])
        rows.append(
            {
                "analysis": suppressed_col,
                "test": f"parent_gene_joint_logistic_{term}",
                "unit": "parent_gene",
                "n_parent_genes": int(len(parent)),
                "n_suppressed_parent_genes": int(parent[suppressed_col].astype(bool).sum()),
                "statistic": float(np.exp(coef)),
                "ci_low": float(np.exp(coef - 1.96 * se)),
                "ci_high": float(np.exp(coef + 1.96 * se)),
                "p_value": float(fit.pvalues[term]),
                "fisher_alternative": np.nan,
                "model_status": "fit",
                "log_odds_coef": coef,
                "se": se,
                "n_parent_genes_model": int(len(d)),
                "covariates": (
                    "dct1_top_decile + dct2_bottom_decile + log1p(n_quantified_phosphosites), "
                    "log1p(parent_n_peptides), parent_abundance_log2, mean_missing_samples"
                ),
            }
        )
    return rows


def _fisher_odds_from_counts(a: int, b: int, c: int, d: int) -> float:
    odds, _ = stats.fisher_exact([[a, b], [c, d]], alternative=FISHER_ALTERNATIVE)
    return float(odds)


def parent_gene_cluster_bootstrap_row_or(
    df: pd.DataFrame,
    flag: str,
    suppressed_col: str,
    n_boot: int = 2000,
    seed: int = 20260526,
) -> dict:
    """Bootstrap row-level odds ratios by resampling parent-gene clusters."""
    sub = df[df["dct1_enrichment_score"].notna()].dropna(subset=["gene_symbol"]).copy()
    if sub.empty:
        return {}
    sub["_flag"] = sub[flag].fillna(False).astype(bool)
    sub["_suppressed"] = sub[suppressed_col].fillna(False).astype(bool)
    counts = (
        sub.groupby("gene_symbol", sort=False)
        .apply(
            lambda g: pd.Series(
                {
                    "a": int((g["_suppressed"] & g["_flag"]).sum()),
                    "b": int((g["_suppressed"] & ~g["_flag"]).sum()),
                    "c": int((~g["_suppressed"] & g["_flag"]).sum()),
                    "d": int((~g["_suppressed"] & ~g["_flag"]).sum()),
                }
            )
        )
        .reset_index(drop=True)
    )
    obs = _fisher_odds_from_counts(int(counts["a"].sum()), int(counts["b"].sum()), int(counts["c"].sum()), int(counts["d"].sum()))
    rng = np.random.default_rng(seed)
    vals = []
    arr = counts[["a", "b", "c", "d"]].to_numpy(dtype=int)
    n = len(arr)
    for _ in range(n_boot):
        draw = arr[rng.integers(0, n, n)]
        total = draw.sum(axis=0)
        vals.append(_fisher_odds_from_counts(int(total[0]), int(total[1]), int(total[2]), int(total[3])))
    vals = np.asarray(vals, dtype=float)
    finite = vals[np.isfinite(vals)]
    return {
        "flag": flag,
        "suppressed_col": suppressed_col,
        "unit": "parent_gene_cluster_bootstrap_of_phosphosite_rows",
        "observed_row_level_odds_ratio": obs,
        "bootstrap_median_odds_ratio": float(np.nanmedian(finite)) if len(finite) else np.nan,
        "ci_low": float(np.nanquantile(finite, 0.025)) if len(finite) else np.nan,
        "ci_high": float(np.nanquantile(finite, 0.975)) if len(finite) else np.nan,
        "n_bootstrap": int(n_boot),
        "n_parent_gene_clusters": int(n),
        "seed": seed,
    }


def _strata_from_parent(parent: pd.DataFrame) -> pd.Series:
    pieces = []
    specs = [
        ("n_quantified_phosphosites", True),
        ("parent_n_peptides", True),
        ("parent_abundance_log2", False),
        ("mean_missing_samples", False),
    ]
    for col, log_first in specs:
        x = pd.to_numeric(parent[col], errors="coerce")
        if log_first:
            x = np.log1p(x)
        try:
            binned = pd.qcut(x.rank(method="first"), q=3, labels=False, duplicates="drop")
        except ValueError:
            binned = pd.Series(0, index=parent.index)
        pieces.append(col + "=" + binned.astype("Int64").astype(str))
    return pd.Series([";".join(vals) for vals in zip(*pieces)], index=parent.index)


def matched_parent_gene_permutation(
    parent: pd.DataFrame,
    flag: str,
    suppressed_col: str,
    n_perm: int = 5000,
    seed: int = 20260526,
) -> dict:
    """Shuffle subtype-prior labels within parent-gene observability strata."""
    use_cols = [flag, suppressed_col, "n_quantified_phosphosites", "parent_n_peptides", "parent_abundance_log2", "mean_missing_samples"]
    d = parent.dropna(subset=use_cols).copy()
    if d.empty:
        return {}
    d["_stratum"] = _strata_from_parent(d)
    y = d[suppressed_col].astype(bool).to_numpy()
    labels = d[flag].astype(bool).to_numpy()
    obs = _fisher_odds_from_counts(
        int((y & labels).sum()),
        int((y & ~labels).sum()),
        int((~y & labels).sum()),
        int((~y & ~labels).sum()),
    )
    rng = np.random.default_rng(seed)
    strata = d["_stratum"].to_numpy()
    hits = 0
    finite = 0
    for _ in range(n_perm):
        perm = labels.copy()
        for stratum in pd.unique(strata):
            idx = np.flatnonzero(strata == stratum)
            if len(idx) > 1:
                perm[idx] = rng.permutation(perm[idx])
        odds = _fisher_odds_from_counts(
            int((y & perm).sum()),
            int((y & ~perm).sum()),
            int((~y & perm).sum()),
            int((~y & ~perm).sum()),
        )
        if np.isfinite(odds):
            finite += 1
            hits += int(odds >= obs)
    return {
        "flag": flag,
        "suppressed_col": suppressed_col,
        "unit": "parent_gene",
        "observed_parent_gene_odds_ratio": obs,
        "n_permutations": int(n_perm),
        "empirical_p_value": float((hits + 1) / (finite + 1)) if finite else np.nan,
        "n_parent_genes": int(len(d)),
        "n_strata": int(d["_stratum"].nunique()),
        "null": "Subtype-prior labels shuffled within parent-gene strata preserving site count, peptide count, parent abundance, and missingness bins.",
        "seed": seed,
    }


def site_count_stratified_permutation(
    df: pd.DataFrame,
    flag: str,
    suppressed_col: str,
    n_perm: int = 2000,
    seed: int = 20260526,
) -> dict:
    """Shuffle subtype-prior flags among parent genes within site-count strata.

    The row-level suppressed/not-suppressed pattern is held fixed; only the
    parent-gene subtype-prior label is permuted within bins of quantified site
    count. This preserves parent-gene site density in the null.
    """
    sub = df[df["dct1_enrichment_score"].notna()].copy()
    sub = sub.dropna(subset=["gene_symbol"])
    if sub.empty:
        return {}
    gene = (
        sub.groupby("gene_symbol", as_index=False)
        .agg(
            flag=(flag, "first"),
            n_sites=("site_id", "count"),
        )
        .copy()
    )
    n_bins = min(5, max(1, gene["n_sites"].nunique()))
    gene["site_count_bin"] = pd.qcut(gene["n_sites"].rank(method="first"), q=n_bins, labels=False, duplicates="drop")
    gene_flag = gene.set_index("gene_symbol")["flag"].astype(bool)
    obs_sub = sub.assign(_perm_flag=sub["gene_symbol"].map(gene_flag).fillna(False).astype(bool))
    obs_odds, _, _, _ = fisher_table(obs_sub.rename(columns={"_perm_flag": "__flag"}), "__flag", suppressed_col)
    rng = np.random.default_rng(seed)
    hits = 0
    finite = 0
    genes = gene["gene_symbol"].to_numpy()
    bins = gene["site_count_bin"].to_numpy()
    flags = gene["flag"].astype(bool).to_numpy()
    gene_index = pd.Series(np.arange(len(gene)), index=gene["gene_symbol"])
    row_gene_idx = sub["gene_symbol"].map(gene_index).to_numpy(dtype=int)
    row_suppressed = sub[suppressed_col].astype(bool).to_numpy()
    for _ in range(n_perm):
        perm_flags = flags.copy()
        for b in pd.unique(bins):
            idx = np.flatnonzero(bins == b)
            perm_flags[idx] = rng.permutation(perm_flags[idx])
        row_flag = perm_flags[row_gene_idx]
        a = int((row_suppressed & row_flag).sum())
        b = int((row_suppressed & ~row_flag).sum())
        c = int((~row_suppressed & row_flag).sum())
        d = int((~row_suppressed & ~row_flag).sum())
        odds, _ = stats.fisher_exact([[a, b], [c, d]], alternative=FISHER_ALTERNATIVE)
        if np.isfinite(odds):
            finite += 1
            hits += int(odds >= obs_odds)
    p_emp = (hits + 1) / (finite + 1) if finite else np.nan
    return {
        "flag": flag,
        "suppressed_col": suppressed_col,
        "observed_row_level_odds_ratio": float(obs_odds),
        "n_permutations": int(n_perm),
        "empirical_p_value": float(p_emp),
        "null": "Subtype-prior flag shuffled among parent genes within quantified-site-count bins",
    }


def matched_null_mean(df: pd.DataFrame, suppressed_col: str, n_draws=5000, seed=20260526):
    rng = np.random.default_rng(seed)
    sub = df[df["dct1_enrichment_score"].notna()].copy()
    sup = sub[sub[suppressed_col]]
    bg = sub[~sub[suppressed_col]]
    obs = sup["dct1_enrichment_score"].mean()
    if len(sup) == 0 or len(bg) == 0:
        return obs, np.nan, np.nan, np.nan, 0
    draws = []
    strata_counts = sup["match_stratum"].fillna("missing").value_counts().to_dict()
    bg_strata = {k: g for k, g in bg.assign(match_stratum=bg["match_stratum"].fillna("missing")).groupby("match_stratum")}
    all_bg = bg
    for _ in range(n_draws):
        sampled = []
        for stratum, n in strata_counts.items():
            pool = bg_strata.get(stratum, all_bg)
            if len(pool) == 0:
                pool = all_bg
            sampled.append(pool.sample(n=n, replace=len(pool) < n, random_state=int(rng.integers(0, 2**31 - 1))))
        draw = pd.concat(sampled, ignore_index=True)
        draws.append(draw["dct1_enrichment_score"].mean())
    draws = np.asarray(draws)
    p_emp = (np.sum(draws >= obs) + 1) / (len(draws) + 1)
    return obs, float(np.median(draws)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975)), float(p_emp)


def run_h2_enrichment(root: Path, phospho_prior: pd.DataFrame):
    rows = []
    one_site = one_representative_site_per_gene(phospho_prior)
    single_position_one_site = one_single_position_representative_site_per_gene(phospho_prior)
    families = [
        ("primary_p05", "is_suppressed_p05", phospho_prior),
        ("strict_q10", "is_suppressed_q10", phospho_prior),
        ("effect_bottom_quartile", "is_effect_bottom_quartile", phospho_prior),
        ("exclude_anchor_genes", "is_suppressed_p05", phospho_prior[~phospho_prior["is_anchor_gene"]]),
        ("exclude_ncc_sites", "is_suppressed_p05", phospho_prior[~phospho_prior["is_ncc_site"]]),
        ("composite_sites_excluded", "is_suppressed_p05", phospho_prior[phospho_prior["is_single_site"]]),
        ("one_site_per_parent_gene", "is_suppressed_p05", one_site),
        ("single_position_one_site_per_parent_gene", "is_suppressed_p05", single_position_one_site),
    ]
    for label, suppressed_col, df in families:
        sub = df[df["dct1_enrichment_score"].notna()].copy()
        sup = sub[sub[suppressed_col]]
        nonsup = sub[~sub[suppressed_col]]
        unit = (
            "single_position_parent_gene_representative_site"
            if label == "single_position_one_site_per_parent_gene"
            else "parent_gene_representative_site"
            if label == "one_site_per_parent_gene"
            else "phosphosite_row"
        )
        mw = stats.mannwhitneyu(
            sup["dct1_enrichment_score"],
            nonsup["dct1_enrichment_score"],
            alternative="greater",
        ) if len(sup) and len(nonsup) else (np.nan, np.nan)
        if hasattr(mw, "statistic"):
            mw_stat, mw_p = mw.statistic, mw.pvalue
        else:
            mw_stat, mw_p = mw
        obs, null_med, null_lo, null_hi, null_p = matched_null_mean(sub, suppressed_col)
        rows.append(
            {
                "analysis": label,
                "test": "mann_whitney_continuous_dct1_score",
                "unit": unit,
                "n_background": len(sub),
                "n_suppressed": len(sup),
                "n_parent_genes": int(sub["gene_symbol"].nunique()),
                "observed_mean_dct1_score_suppressed": obs,
                "background_mean_dct1_score": nonsup["dct1_enrichment_score"].mean(),
                "statistic": mw_stat,
                "p_value": mw_p,
                "fisher_alternative": np.nan,
                "fold_enrichment": np.nan,
                "null_median": np.nan,
                "null_lo": np.nan,
                "null_hi": np.nan,
            }
        )
        rows.append(
            {
                "analysis": label,
                "test": "matched_null_mean_dct1_score",
                "unit": unit,
                "n_background": len(sub),
                "n_suppressed": len(sup),
                "n_parent_genes": int(sub["gene_symbol"].nunique()),
                "observed_mean_dct1_score_suppressed": obs,
                "background_mean_dct1_score": nonsup["dct1_enrichment_score"].mean(),
                "statistic": obs,
                "p_value": null_p,
                "fisher_alternative": np.nan,
                "fold_enrichment": np.nan,
                "null_median": null_med,
                "null_lo": null_lo,
                "null_hi": null_hi,
            }
        )
        for flag in ["dct1_core_fdr", "dct1_top_quartile", "dct1_top_decile", "dct2_bottom_quartile", "dct2_bottom_decile"]:
            odds, p, fold, arr = fisher_table(sub, flag, suppressed_col)
            rows.append(
                {
                    "analysis": label,
                    "test": f"fisher_{flag}",
                    "unit": unit,
                    "n_background": len(sub),
                    "n_suppressed": len(sup),
                    "n_parent_genes": int(sub["gene_symbol"].nunique()),
                    "observed_mean_dct1_score_suppressed": obs,
                    "background_mean_dct1_score": nonsup["dct1_enrichment_score"].mean(),
                    "statistic": odds,
                    "p_value": p,
                    "fisher_alternative": FISHER_ALTERNATIVE,
                    "fold_enrichment": fold,
                    "null_median": np.nan,
                    "null_lo": np.nan,
                    "null_hi": np.nan,
                    "table_suppressed_in_flag": int(arr[0, 0]),
                    "table_suppressed_not_flag": int(arr[0, 1]),
                    "table_background_in_flag": int(arr[1, 0]),
                    "table_background_not_flag": int(arr[1, 1]),
                }
            )
    summary = pd.DataFrame(rows)
    summary["q_value"] = bh(summary["p_value"])
    summary.to_csv(root / "h2_enrichment" / "h2_dct1_phosphosite_enrichment_summary.tsv", sep="\t", index=False)
    phospho_prior.to_csv(root / "h2_enrichment" / "h2_dct1_phosphosite_enrichment_background.tsv", sep="\t", index=False)

    primary_continuous = summary[
        (summary["analysis"] == "primary_p05")
        & (summary["test"].isin(["mann_whitney_continuous_dct1_score", "matched_null_mean_dct1_score"]))
    ]
    primary_binary = summary[
        (summary["analysis"] == "primary_p05")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    anchor_continuous = summary[
        (summary["analysis"] == "exclude_anchor_genes")
        & (summary["test"].isin(["mann_whitney_continuous_dct1_score", "matched_null_mean_dct1_score"]))
    ]
    anchor_binary = summary[
        (summary["analysis"] == "exclude_anchor_genes")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    ncc_binary = summary[
        (summary["analysis"] == "exclude_ncc_sites")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    one_site_binary = summary[
        (summary["analysis"] == "one_site_per_parent_gene")
        & (summary["test"].isin(["fisher_dct1_top_quartile", "fisher_dct1_top_decile"]))
    ]
    dct2_primary_binary = summary[
        (summary["analysis"] == "primary_p05")
        & (summary["test"].isin(["fisher_dct2_bottom_quartile", "fisher_dct2_bottom_decile"]))
    ]
    dct2_one_site_binary = summary[
        (summary["analysis"] == "one_site_per_parent_gene")
        & (summary["test"].isin(["fisher_dct2_bottom_quartile", "fisher_dct2_bottom_decile"]))
    ]
    primary_continuous_pass = bool((primary_continuous["q_value"] <= 0.10).any())
    primary_binary_pass = bool((primary_binary["q_value"] <= 0.10).any())
    survives_anchor_continuous = bool((anchor_continuous["q_value"] <= 0.10).any())
    survives_anchor_binary = bool((anchor_binary["q_value"] <= 0.10).any())
    survives_ncc_binary = bool((ncc_binary["q_value"] <= 0.10).any())
    survives_one_site_per_gene = bool((one_site_binary["q_value"] <= 0.10).any())
    dct2_primary_binary_pass = bool((dct2_primary_binary["q_value"] <= 0.10).any())
    dct2_one_site_per_gene_pass = bool((dct2_one_site_binary["q_value"] <= 0.10).any())
    passes = primary_continuous_pass or primary_binary_pass
    survives_anchor_exclusion = survives_anchor_continuous or survives_anchor_binary
    if primary_binary_pass and dct2_primary_binary_pass:
        interpretation = (
            "Distal-nephron subtype-prior enrichment is supported. The DCT1 top-decile bin is NCC/SPAK/WNK-motivated, "
            "and the DCT2-bottom comparator also passes; do not claim DCT1 exclusivity."
        )
    elif primary_binary_pass and survives_anchor_binary and survives_ncc_binary:
        interpretation = (
            "DCT1-high percentile enrichment is supported and survives anchor/NCC exclusion, "
            "but continuous DCT1-score support is weak; claim as DCT1-prioritized subset enrichment only."
        )
    elif primary_continuous_pass and survives_anchor_continuous:
        interpretation = "Broad continuous DCT1-prior enrichment of suppressed phosphosites is supported."
    elif passes:
        interpretation = "A limited DCT1-prior signal is present, but sensitivity analyses do not support a broad DCT1-specific claim."
    else:
        interpretation = "DCT-subtype-prior enrichment is not supported under the pre-specified tests."
    verdict = {
        "primary_h2_supported_at_fdr_0_10": passes,
        "survives_anchor_gene_exclusion_at_fdr_0_10": survives_anchor_exclusion,
        "primary_continuous_pass": primary_continuous_pass,
        "primary_binary_percentile_pass": primary_binary_pass,
        "survives_anchor_continuous": survives_anchor_continuous,
        "survives_anchor_binary_percentile": survives_anchor_binary,
        "survives_ncc_binary_percentile": survives_ncc_binary,
        "survives_one_site_per_parent_gene": survives_one_site_per_gene,
        "dct2_primary_binary_percentile_pass": dct2_primary_binary_pass,
        "dct2_one_site_per_parent_gene_pass": dct2_one_site_per_gene_pass,
        "interpretation": interpretation,
        "claim_caution": "OSD-462 remains whole-kidney phosphoproteomics; DCT-subtype specificity is reference-prior inference only.",
        "single_position_note": "composite_sites_excluded removes composite/multi-position rows; one_site_per_parent_gene is a one-phosphosite-row-per-parent-gene sensitivity.",
        "fisher_alternative": FISHER_ALTERNATIVE,
    }
    (root / "h2_enrichment" / "h2_dct1_enrichment_verdict.json").write_text(json.dumps(verdict, indent=2))
    summary.to_csv(root / "h2_enrichment" / "h2_dct1_sensitivity_summary.tsv", sep="\t", index=False)
    parent = run_h2_parent_gene_level(root, phospho_prior)
    site_perm = pd.DataFrame(
        [
            site_count_stratified_permutation(phospho_prior, "dct1_top_decile", "is_suppressed_p05"),
            site_count_stratified_permutation(phospho_prior, "dct1_top_quartile", "is_suppressed_p05"),
            site_count_stratified_permutation(phospho_prior, "dct2_bottom_decile", "is_suppressed_p05"),
            site_count_stratified_permutation(phospho_prior, "dct2_bottom_quartile", "is_suppressed_p05"),
        ]
    )
    site_perm.to_csv(root / "h2_enrichment" / "h2_dct1_site_count_stratified_permutation.tsv", sep="\t", index=False)

    def cluster_row(analysis: str, df: pd.DataFrame, flag: str) -> dict:
        row = parent_gene_cluster_bootstrap_row_or(df, flag, "is_suppressed_p05")
        row["analysis"] = analysis
        return row

    cluster_boot = pd.DataFrame(
        [
            cluster_row("all_phosphosite_rows", phospho_prior, "dct1_top_decile"),
            cluster_row("all_phosphosite_rows", phospho_prior, "dct2_bottom_decile"),
            cluster_row("single_position_one_site_per_parent_gene", single_position_one_site, "dct1_top_decile"),
            cluster_row("single_position_one_site_per_parent_gene", single_position_one_site, "dct2_bottom_decile"),
        ]
    )
    cluster_boot.to_csv(root / "h2_enrichment" / "h2_dct_extreme_bin_cluster_bootstrap.tsv", sep="\t", index=False)
    matched_perm = pd.DataFrame(
        [
            matched_parent_gene_permutation(parent, "dct1_top_decile", "any_suppressed_p05"),
            matched_parent_gene_permutation(parent, "dct2_bottom_decile", "any_suppressed_p05"),
            matched_parent_gene_permutation(parent, "dct1_top_decile", "any_effect_bottom_quartile"),
            matched_parent_gene_permutation(parent, "dct2_bottom_decile", "any_effect_bottom_quartile"),
        ]
    )
    matched_perm.to_csv(root / "h2_enrichment" / "h2_dct_matched_parent_gene_permutation.tsv", sep="\t", index=False)
    run_threshold_free_enrichment(root, phospho_prior, parent)
    run_percentile_threshold_sensitivity(root, phospho_prior, parent)
    return summary, verdict


def run_h2_parent_gene_level(root: Path, phospho_prior: pd.DataFrame) -> pd.DataFrame:
    parent = build_parent_gene_table(phospho_prior)
    parent.to_csv(root / "h2_enrichment" / "h2_dct1_parent_gene_background.tsv", sep="\t", index=False)
    rows = []
    for suppressed_col in ["any_suppressed_p05", "any_suppressed_q10"]:
        for flag in ["dct1_top_decile", "dct1_top_quartile", "dct2_bottom_quartile", "dct2_bottom_decile"]:
            odds, p, fold, arr = fisher_table(parent, flag, suppressed_col)
            rows.append(
                {
                    "analysis": suppressed_col,
                    "test": f"parent_gene_fisher_{flag}",
                    "unit": "parent_gene",
                    "n_parent_genes": int(len(parent)),
                    "n_suppressed_parent_genes": int(parent[suppressed_col].astype(bool).sum()),
                    "statistic": odds,
                    "p_value": p,
                    "fisher_alternative": FISHER_ALTERNATIVE,
                    "fold_enrichment": fold,
                    "table_suppressed_in_flag": int(arr[0, 0]),
                    "table_suppressed_not_flag": int(arr[0, 1]),
                    "table_background_in_flag": int(arr[1, 0]),
                    "table_background_not_flag": int(arr[1, 1]),
                }
            )
            logit = parent_gene_logistic(parent, flag, suppressed_col)
            rows.append(
                {
                    "analysis": suppressed_col,
                    "test": f"parent_gene_logistic_{flag}",
                    "unit": "parent_gene",
                    "n_parent_genes": int(len(parent)),
                    "n_suppressed_parent_genes": int(parent[suppressed_col].astype(bool).sum()),
                    "statistic": logit.get("odds_ratio", np.nan),
                    "ci_low": logit.get("ci_low", np.nan),
                    "ci_high": logit.get("ci_high", np.nan),
                    "p_value": logit.get("p_value", np.nan),
                    "fisher_alternative": np.nan,
                    "model_status": logit.get("model_status", "unknown"),
                    "log_odds_coef": logit.get("log_odds_coef", np.nan),
                    "se": logit.get("se", np.nan),
                    "n_parent_genes_model": logit.get("n_parent_genes_model", np.nan),
                    "covariates": "log1p(n_quantified_phosphosites), log1p(parent_n_peptides), parent_abundance_log2",
                }
            )
        rows.extend(parent_gene_joint_logistic(parent, suppressed_col))
    out = pd.DataFrame(rows)
    out["q_value"] = bh(out["p_value"])
    out.to_csv(root / "h2_enrichment" / "h2_dct1_parent_gene_level_summary.tsv", sep="\t", index=False)
    return parent


def _gradient_directional_read(rho: float, p: float, score_label: str) -> str:
    """Sign-aware interpretation of a DCT1-score vs phosphosite-effect Spearman.

    A continuous suppression gradient predicts NEGATIVE rho (more negative effect
    at higher DCT1 prior). The label reports the OBSERVED sign so a positive rho is
    never silently read as supportive of suppression.
    """
    if not np.isfinite(rho):
        return f"insufficient data for Spearman of DCT1 prior vs {score_label}"
    sig = bool(np.isfinite(p) and p < 0.05)
    if rho < 0:
        trend = "higher DCT1 prior tracks a more negative (more suppressed) effect"
        supports = sig
    else:
        trend = ("higher DCT1 prior tracks a less negative (less suppressed) effect, "
                 "opposite the suppression direction")
        supports = False
    return (
        f"observed rho={rho:+.3f} ({'p<0.05' if sig else 'n.s.'}); {trend}; "
        f"continuous gradient {'supported' if supports else 'NOT supported'} -- the "
        f"supported result is extreme-bin enrichment, not a graded DCT1 trend"
    )


def run_threshold_free_enrichment(root: Path, phospho_prior: pd.DataFrame, parent: pd.DataFrame) -> pd.DataFrame:
    """Effect-rank checks that do not define the set by nominal p value."""
    rows = []
    row = phospho_prior.dropna(subset=["dct1_enrichment_score", "phospho_effect"]).copy()
    if len(row):
        rho, p = stats.spearmanr(row["dct1_enrichment_score"], row["phospho_effect"], nan_policy="omit")
        rows.append(
            {
                "analysis": "all_phosphosite_rows",
                "test": "spearman_dct1_score_vs_signed_phosphosite_effect",
                "unit": "phosphosite_row",
                "statistic": rho,
                "p_value": p,
                "directional_read": _gradient_directional_read(rho, p, "signed phosphosite effect"),
                "n_units": int(len(row)),
                "n_parent_genes": int(row["gene_symbol"].nunique()),
            }
        )
        for flag in ["dct1_top_decile", "dct2_bottom_decile", "dct1_top_quartile", "dct2_bottom_quartile"]:
            d = row.dropna(subset=[flag]).copy()
            in_flag = d[d[flag].astype(bool)]["phospho_effect"]
            out_flag = d[~d[flag].astype(bool)]["phospho_effect"]
            mw = stats.mannwhitneyu(in_flag, out_flag, alternative="less") if len(in_flag) and len(out_flag) else (np.nan, np.nan)
            rows.append(
                {
                    "analysis": "all_phosphosite_rows",
                    "test": f"mann_whitney_{flag}_more_negative_effect",
                    "unit": "phosphosite_row",
                    "statistic": mw.statistic if hasattr(mw, "statistic") else mw[0],
                    "p_value": mw.pvalue if hasattr(mw, "pvalue") else mw[1],
                    "mean_effect_in_flag": float(in_flag.mean()) if len(in_flag) else np.nan,
                    "mean_effect_outside_flag": float(out_flag.mean()) if len(out_flag) else np.nan,
                    "directional_read": "one-sided less: subtype-prior bin has more negative signed phosphosite effects",
                    "n_units": int(len(d)),
                    "n_parent_genes": int(d["gene_symbol"].nunique()),
                }
            )

    p = parent.dropna(subset=["dct1_enrichment_score", "most_negative_phospho_effect"]).copy()
    if len(p):
        rho, pval = stats.spearmanr(p["dct1_enrichment_score"], p["most_negative_phospho_effect"], nan_policy="omit")
        rows.append(
            {
                "analysis": "parent_gene_most_negative_site",
                "test": "spearman_dct1_score_vs_most_negative_parent_site_effect",
                "unit": "parent_gene",
                "statistic": rho,
                "p_value": pval,
                "directional_read": _gradient_directional_read(rho, pval, "parent-gene most-negative-site effect"),
                "n_units": int(len(p)),
                "n_parent_genes": int(len(p)),
            }
        )
        for flag in ["dct1_top_decile", "dct2_bottom_decile", "dct1_top_quartile", "dct2_bottom_quartile"]:
            in_flag = p[p[flag].astype(bool)]["most_negative_phospho_effect"]
            out_flag = p[~p[flag].astype(bool)]["most_negative_phospho_effect"]
            mw = stats.mannwhitneyu(in_flag, out_flag, alternative="less") if len(in_flag) and len(out_flag) else (np.nan, np.nan)
            rows.append(
                {
                    "analysis": "parent_gene_most_negative_site",
                    "test": f"mann_whitney_{flag}_more_negative_effect",
                    "unit": "parent_gene",
                    "statistic": mw.statistic if hasattr(mw, "statistic") else mw[0],
                    "p_value": mw.pvalue if hasattr(mw, "pvalue") else mw[1],
                    "mean_effect_in_flag": float(in_flag.mean()) if len(in_flag) else np.nan,
                    "mean_effect_outside_flag": float(out_flag.mean()) if len(out_flag) else np.nan,
                    "directional_read": "one-sided less: subtype-prior bin has more negative parent-gene phosphosite effects",
                    "n_units": int(len(p)),
                    "n_parent_genes": int(len(p)),
                }
            )
    out = pd.DataFrame(rows)
    if len(out):
        out["q_value"] = bh(out["p_value"])
    out.to_csv(root / "h2_enrichment" / "h2_dct_threshold_free_summary.tsv", sep="\t", index=False)
    return out


def run_percentile_threshold_sensitivity(root: Path, phospho_prior: pd.DataFrame, parent: pd.DataFrame) -> pd.DataFrame:
    """Sweep DCT1-top and DCT2-bottom percentile cutoffs at row and parent level."""
    parent_score = parent.dropna(subset=["dct1_enrichment_score"]).copy()
    rows = []
    for pct in [5, 10, 15, 20, 25]:
        hi = parent_score["dct1_enrichment_score"].quantile(1 - pct / 100)
        lo = parent_score["dct1_enrichment_score"].quantile(pct / 100)
        gene_flags = parent_score[["gene_symbol", "dct1_enrichment_score"]].copy()
        gene_flags[f"dct1_top_{pct}pct"] = gene_flags["dct1_enrichment_score"] >= hi
        gene_flags[f"dct2_bottom_{pct}pct"] = gene_flags["dct1_enrichment_score"] <= lo
        row_df = phospho_prior.merge(gene_flags[["gene_symbol", f"dct1_top_{pct}pct", f"dct2_bottom_{pct}pct"]], on="gene_symbol", how="left")
        parent_df = parent.merge(gene_flags[["gene_symbol", f"dct1_top_{pct}pct", f"dct2_bottom_{pct}pct"]], on="gene_symbol", how="left")
        for unit, df, suppressed_col in [
            ("phosphosite_row", row_df, "is_suppressed_p05"),
            ("parent_gene", parent_df, "any_suppressed_p05"),
        ]:
            for direction, flag in [("dct1_top", f"dct1_top_{pct}pct"), ("dct2_bottom", f"dct2_bottom_{pct}pct")]:
                odds, pval, fold, arr = fisher_table(df, flag, suppressed_col)
                rows.append(
                    {
                        "unit": unit,
                        "suppressed_col": suppressed_col,
                        "bin_direction": direction,
                        "percentile_cutoff": pct,
                        "flag": flag,
                        "odds_ratio": odds,
                        "p_value": pval,
                        "fold_enrichment": fold,
                        "table_suppressed_in_flag": int(arr[0, 0]),
                        "table_suppressed_not_flag": int(arr[0, 1]),
                        "table_background_in_flag": int(arr[1, 0]),
                        "table_background_not_flag": int(arr[1, 1]),
                    }
                )
    out = pd.DataFrame(rows)
    out["q_value"] = bh(out["p_value"])
    out.to_csv(root / "h2_enrichment" / "h2_dct_percentile_threshold_sensitivity.tsv", sep="\t", index=False)
    return out


def run_pxd_antialignment(root: Path, pxd: pd.DataFrame, phospho_prior: pd.DataFrame):
    changed = pxd[pxd["ddavp_effect_log2"].notna()].copy()
    changed = changed[changed["gene_symbol"].notna() & changed["site_position"].notna()]
    changed["site_position_int"] = changed["site_position"].astype(int)
    changed_agg = (
        changed.groupby(["gene_symbol", "site_position_int"], as_index=False)
        .agg(
            ddavp_effect_log2=("ddavp_effect_log2", "mean"),
            p_neg_log10=("p_neg_log10", "max"),
            n_pxd_rows=("site_id", "count"),
        )
    )
    osd = phospho_prior[phospho_prior["is_single_site"]].copy()
    osd["site_position_int"] = osd["site_position_int"].astype("Int64")
    shared = osd.merge(changed_agg, on=["gene_symbol", "site_position_int"], how="inner")
    shared["site_id_shared"] = shared["gene_symbol"] + ":" + shared["site_position_int"].astype(str)
    shared.to_csv(root / "h2_pxd" / "pxd001729_osd462_shared_phosphosites.tsv", sep="\t", index=False)
    target_shared = shared[shared["gene_symbol"].isin(TRANSPORT_TARGETS)].copy()

    def summarize(label, df):
        c = cosine(df["ddavp_effect_log2"], df["phospho_effect"])
        rng = np.random.default_rng(20260526)
        boots = []
        if len(df) >= 3:
            for _ in range(2000):
                idx = rng.integers(0, len(df), len(df))
                boots.append(cosine(df["ddavp_effect_log2"].to_numpy()[idx], df["phospho_effect"].to_numpy()[idx]))
        boots = np.asarray(boots, dtype=float)
        return {
            "comparison": label,
            "n_shared_sites": len(df),
            "n_shared_genes": df["gene_symbol"].nunique() if len(df) else 0,
            "cosine_ddavp_vs_spaceflight": c,
            "ci_low": np.nanquantile(boots, 0.025) if len(boots) else np.nan,
            "ci_high": np.nanquantile(boots, 0.975) if len(boots) else np.nan,
            "interpretation": (
                "anti_aligned_with_dct_activation"
                if len(df) >= 3 and c < 0 and np.nanquantile(boots, 0.975) < 0
                else "suggestive_or_not_stable_reference_signal"
            ),
        }

    rows = [
        summarize("all_shared_single_sites", shared),
        summarize("transport_target_shared_sites", target_shared),
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(root / "h2_pxd" / "h2_pxd001729_ddavp_antialignment_summary.tsv", sep="\t", index=False)

    missing_targets = sorted(TRANSPORT_TARGETS - set(target_shared["gene_symbol"]))
    verdict = {
        "n_all_shared_single_sites": int(len(shared)),
        "n_transport_target_shared_sites": int(len(target_shared)),
        "transport_targets_missing_from_shared_map": missing_targets,
        "claim_caution": "PXD001729 is mpkDCT cell-line dDAVP phosphoproteomics, not spaceflight tissue.",
        "verdict": summary.loc[0, "interpretation"] if len(summary) else "not_testable",
    }
    (root / "h2_pxd" / "h2_pxd001729_ddavp_antialignment_verdict.json").write_text(json.dumps(verdict, indent=2))
    failures = phospho_prior[phospho_prior["gene_symbol"].isin(TRANSPORT_TARGETS) & phospho_prior["is_single_site"]].copy()
    failures["present_in_pxd_changed_same_site"] = failures.set_index(["gene_symbol", "site_position_int"]).index.isin(
        changed_agg.set_index(["gene_symbol", "site_position_int"]).index
    )
    failures.to_csv(root / "h2_pxd" / "pxd001729_osd462_mapping_failures.tsv", sep="\t", index=False)
    return summary, verdict


def run_klhl3_check(root: Path, pxd: pd.DataFrame, phospho_prior: pd.DataFrame, proteins_prior: pd.DataFrame):
    genes = ["Klhl3", "Cul3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Slc12a3", "Nedd4l", "Sgk1"]
    osd_sites = phospho_prior[phospho_prior["gene_symbol"].isin(genes)].copy()
    pxd_sites = pxd[pxd["gene_symbol"].isin(genes)].copy()
    prot = proteins_prior[proteins_prior["gene_symbol"].isin(genes)].copy()
    osd_sites.to_csv(root / "h2_klhl3" / "h2_klhl3_cul3_site_coverage.tsv", sep="\t", index=False)
    combined = []
    for g in genes:
        combined.append(
            {
                "gene_symbol": g,
                "osd462_n_phosphosites": int((osd_sites["gene_symbol"] == g).sum()),
                "osd462_min_phospho_p": osd_sites.loc[osd_sites["gene_symbol"] == g, "phospho_p_value"].min(),
                "osd462_mean_phospho_effect": osd_sites.loc[osd_sites["gene_symbol"] == g, "phospho_effect"].mean(),
                "pxd001729_n_phosphosites": int((pxd_sites["gene_symbol"] == g).sum()),
                "pxd001729_n_ddavp_changed": int(pxd_sites.loc[pxd_sites["gene_symbol"] == g, "ddavp_effect_log2"].notna().sum()),
                "osd462_protein_effect": prot.loc[prot["gene_symbol"] == g, "flight_effect"].mean(),
                "osd462_rna_effect": prot.loc[prot["gene_symbol"] == g, "osd462_rna_effect"].mean()
                if "osd462_rna_effect" in prot.columns
                else np.nan,
            }
        )
    effects = pd.DataFrame(combined)
    effects.to_csv(root / "h2_klhl3" / "h2_klhl3_cul3_effects.tsv", sep="\t", index=False)
    has_klhl3_433 = bool(((osd_sites["gene_symbol"] == "Klhl3") & (osd_sites["site_position"].astype(str) == "433")).any())
    interpretation = [
        "# KLHL3/CUL3 exploratory check",
        "",
        f"KLHL3 Ser433 detected in OSD-462 targeted phosphosite table: {has_klhl3_433}.",
        "",
        "This check is exploratory because public OSD-462 does not include ubiquitinomics.",
        "A KLHL3/WNK turnover mechanism cannot be distinguished from ionic/osmotic WNK-SPAK suppression without perturbation or ubiquitin-remnant data.",
    ]
    (root / "h2_klhl3" / "h2_klhl3_cul3_interpretation.md").write_text("\n".join(interpretation) + "\n")


def sample_to_animal(sample: str) -> tuple[str, str] | None:
    m = re.search(r"RR10_KDN_WT_(FLT|GC|BSL|VIV)_([A-Z])(\d+)", sample)
    if not m:
        return None
    cond_code, _, number = m.groups()
    cond = {
        "FLT": "Space Flight",
        "GC": "Ground Control",
        "BSL": "Basal",
        "VIV": "Vivarium",
    }[cond_code]
    return cond, f"{cond}|{int(number)}"


def approximate_bayes_linear(y, X, n_draws=10000, seed=20260526):
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    beta_hat = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta_hat
    df = max(n - p, 1)
    sse = float(np.sum(resid**2))
    xtx_inv = np.linalg.pinv(X.T @ X)
    sigma2 = sse / rng.chisquare(df, size=n_draws)
    draws = np.empty((n_draws, p))
    for k in range(n_draws):
        cov = sigma2[k] * xtx_inv
        draws[k] = rng.multivariate_normal(beta_hat, cov)
    return draws, np.sqrt(sigma2), beta_hat, resid


def run_mediation(root: Path):
    # Stage 0 invalidated the outcome used by the historical mediation:
    # position-indexed T53 and S383 rows are co-modified phosphoforms, and no
    # isolated canonical NCC/SPAK feature qualifies. Fail closed.
    reason = (
        "not_run: zero isolated canonical OSD-462 NCC/SPAK assay features; "
        "historical ncc_activity_score_regulatory is invalid as an activity outcome"
    )
    pd.DataFrame([{"status": "UNRESOLVED", "reason": reason}]).to_csv(
        root / "h3_mediation" / "h3_mediation_input_scores.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(
        columns=[
            "mediator",
            "parameter",
            "posterior_median",
            "ci_low",
            "ci_high",
            "p_less_than_zero",
            "p_greater_than_zero",
            "n_animals",
            "model",
        ]
    ).to_csv(
        root / "h3_mediation" / "h3_mediation_model_summary.tsv",
        sep="\t",
        index=False,
    )
    pd.DataFrame(columns=["status", "reason"]).to_csv(
        root / "h3_mediation" / "h3_mediation_power_simulation.tsv",
        sep="\t",
        index=False,
    )
    (root / "h3_mediation" / "h3_mediation_verdict.json").write_text(
        json.dumps({"status": "UNRESOLVED", "reason": reason}, indent=2)
    )
    return

    pheno = read_tsv(PHENO / "phenotype_anchor_per_animal.tsv")
    scores = pd.read_csv(CELLTYPE / "osd462_compartment_scores_per_sample.tsv", sep="\t", index_col=0).T
    rows = []
    for sample, vals in scores.iterrows():
        parsed = sample_to_animal(sample)
        if parsed is None:
            continue
        condition, animal = parsed
        d = {"sample": sample, "condition": condition, "animal": animal}
        d.update(vals.to_dict())
        rows.append(d)
    comp = pd.DataFrame(rows)
    comp = comp.groupby(["animal", "condition"], as_index=False).mean(numeric_only=True)
    merged = pheno.merge(comp, on=["animal", "condition"], how="inner")
    merged = merged[merged["condition"].isin(["Space Flight", "Ground Control"])].copy()
    merged["flight"] = (merged["condition"] == "Space Flight").astype(float)
    merged["matrix_endothelial_composite"] = merged[[c for c in ["endothelial", "stromal_fibroblast"] if c in merged]].mean(axis=1)
    merged.to_csv(root / "h3_mediation" / "h3_mediation_input_scores.tsv", sep="\t", index=False)

    mediators = [
        "endothelial",
        "stromal_fibroblast",
        "dct_identity",
        "matrix_endothelial_composite",
    ]
    summary_rows = []
    posterior_frames = []
    power_rows = []
    for idx, mediator in enumerate(mediators):
        df = merged[["flight", mediator, "ncc_activity_score_regulatory"]].dropna().copy()
        if len(df) < 10:
            continue
        M = (df[mediator] - df[mediator].mean()) / df[mediator].std(ddof=1)
        Y = (df["ncc_activity_score_regulatory"] - df["ncc_activity_score_regulatory"].mean()) / df[
            "ncc_activity_score_regulatory"
        ].std(ddof=1)
        X_med = np.column_stack([np.ones(len(df)), df["flight"].to_numpy()])
        a_draws, sigma_m, a_hat, resid_m = approximate_bayes_linear(M, X_med, seed=20260526 + idx)
        X_out = np.column_stack([np.ones(len(df)), df["flight"].to_numpy(), M.to_numpy()])
        b_draws, sigma_y, b_hat, resid_y = approximate_bayes_linear(Y, X_out, seed=20260626 + idx)
        a = a_draws[:, 1]
        c_prime = b_draws[:, 1]
        b = b_draws[:, 2]
        indirect = a * b
        total = c_prime + indirect
        for name, draws in [("a", a), ("b", b), ("c_prime", c_prime), ("indirect", indirect), ("total", total)]:
            summary_rows.append(
                {
                    "mediator": mediator,
                    "parameter": name,
                    "posterior_median": np.median(draws),
                    "ci_low": np.quantile(draws, 0.025),
                    "ci_high": np.quantile(draws, 0.975),
                    "p_less_than_zero": np.mean(draws < 0),
                    "p_greater_than_zero": np.mean(draws > 0),
                    "n_animals": len(df),
                    "model": "approximate_bayesian_ols_weak_prior_fallback",
                }
            )
        posterior_frames.append(
            pd.DataFrame(
                {
                    "mediator": mediator,
                    "draw": np.arange(len(indirect)),
                    "a": a,
                    "b": b,
                    "c_prime": c_prime,
                    "indirect": indirect,
                    "total": total,
                }
            )
        )

        # Simple future-n simulation using posterior median paths and residual SDs.
        path_a = float(np.median(a))
        path_b = float(np.median(b))
        path_c = float(np.median(c_prime))
        sd_m = float(np.median(sigma_m))
        sd_y = float(np.median(sigma_y))
        rng = np.random.default_rng(20260726 + idx)
        for n_total in [20, 30, 40, 60, 80, 100, 140]:
            sign_hits = 0
            ci_hits = 0
            sims = 300
            n1 = n_total // 2
            n0 = n_total - n1
            xsim = np.r_[np.zeros(n0), np.ones(n1)]
            for _ in range(sims):
                msim = path_a * xsim + rng.normal(0, sd_m, size=n_total)
                ysim = path_c * xsim + path_b * msim + rng.normal(0, sd_y, size=n_total)
                xm = np.column_stack([np.ones(n_total), xsim])
                xy = np.column_stack([np.ones(n_total), xsim, msim])
                beta_m = np.linalg.lstsq(xm, msim, rcond=None)[0]
                beta_y = np.linalg.lstsq(xy, ysim, rcond=None)[0]
                resid_m_sim = msim - xm @ beta_m
                resid_y_sim = ysim - xy @ beta_y
                cov_m = np.sum(resid_m_sim**2) / max(n_total - xm.shape[1], 1) * np.linalg.pinv(xm.T @ xm)
                cov_y = np.sum(resid_y_sim**2) / max(n_total - xy.shape[1], 1) * np.linalg.pinv(xy.T @ xy)
                ad = beta_m[1]
                bd = beta_y[2]
                se_a = math.sqrt(max(cov_m[1, 1], 0))
                se_b = math.sqrt(max(cov_y[2, 2], 0))
                ind = ad * bd
                se_ind = math.sqrt((bd**2) * (se_a**2) + (ad**2) * (se_b**2))
                lo = ind - 1.96 * se_ind
                hi = ind + 1.96 * se_ind
                sign_hits += int(ind < 0)
                ci_hits += int(ind < 0 and hi < 0)
            power_rows.append(
                {
                    "mediator": mediator,
                    "n_total": n_total,
                    "directional_sign_recovery": sign_hits / sims,
                    "ci_exclusion_power_estimate": ci_hits / sims,
                    "simulation_note": "Sobel-style interval exclusion for posterior-median effect; still approximate and not causal proof",
                }
            )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(root / "h3_mediation" / "h3_mediation_model_summary.tsv", sep="\t", index=False)
    if posterior_frames:
        pd.concat(posterior_frames, ignore_index=True).to_csv(
            root / "h3_mediation" / "h3_mediation_posterior_draws.tsv.gz", sep="\t", index=False, compression="gzip"
        )
    pd.DataFrame(power_rows).to_csv(root / "h3_mediation" / "h3_mediation_power_simulation.tsv", sep="\t", index=False)
    indirect_rows = summary[summary["parameter"] == "indirect"].copy()
    verdict = {
        "model_caveat": "Approximate Bayesian linear-model fallback; cross-sectional n=20 bulk tissue, not causal proof.",
        "mediators": indirect_rows.to_dict(orient="records"),
        "overall_interpretation": "Use as exploratory covariance decomposition. Wide intervals, cross-sectional data, and composition controls prevent mechanistic claims.",
    }
    (root / "h3_mediation" / "h3_mediation_verdict.json").write_text(json.dumps(verdict, indent=2))


def write_manifest(root: Path):
    input_paths = [
        Path("docs/v11_execution_research_plan.md"),
        DCT_PRIOR_DIR / "gse228367_dct1_vs_dct2_de.tsv",
        root / "baseline" / "rna_recurrence_gene_set_members.tsv",
        OSD462 / "osd462_rna_pathway_effects.tsv",
        OSD462 / "phospho_all_sites.tsv",
        OSD462 / "protein_effects_gene_anyplex.tsv",
        PHENO / "phenotype_anchor_per_animal.tsv",
        CELLTYPE / "osd462_compartment_scores_per_sample.tsv",
        PXD_DIR / "41598_2015_BFsrep12829_MOESM6_ESM.xls",
    ]
    manifest = {
        "analysis": "v11 core DCT-subtype-prior phosphoproteome and covariance-decomposition analysis",
        "run_root": str(root),
        "claim_discipline": "GSE228367 and PXD001729 are reference priors/scaffolds, not spaceflight evidence.",
        "inputs": [{"path": str(p), "sha256": sha256(p)} for p in input_paths if p.exists()],
    }
    (root / "manifests" / "v11_core_manifest.json").write_text(json.dumps(manifest, indent=2))


def main():
    global RUN_ROOT, DCT_PRIOR_DIR

    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", default=str(RUN_ROOT))
    args = parser.parse_args()
    root = Path(args.run_root)
    RUN_ROOT = root
    DCT_PRIOR_DIR = RUN_ROOT / "dct_prior"
    ensure_dirs(root)

    baseline_lock(root)
    pxd = load_pxd_tables(root)
    _, proteins_prior, phospho_prior = build_dct_prior_mapping(root)
    run_dct_prior_sensitivity(root, phospho_prior)
    run_h2_enrichment(root, phospho_prior)
    run_pxd_antialignment(root, pxd, phospho_prior)
    run_klhl3_check(root, pxd, phospho_prior, proteins_prior)
    run_mediation(root)
    write_manifest(root)
    print(f"v11 core analysis complete: {root}")


if __name__ == "__main__":
    main()
