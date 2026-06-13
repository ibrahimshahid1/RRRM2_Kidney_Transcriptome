"""Protocol-bound external cohort validation and context mapping."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT, bh_fdr, load_id_map
from src.validation.external_replication import (
    compare_directional_replication,
    registry_for_analysis,
    require_protocol,
    validate_study_scope,
)


PPAR_SYMBOLS = {
    "Ppara", "Ppard", "Pparg", "Rxra", "Rxrb", "Rxrg",
    "Acox1", "Acox2", "Cpt1a", "Cpt1b", "Ehhadh", "Acaa1a", "Acaa1b",
    "Acadl", "Acadm", "Acsl1", "Acsl3", "Acsl4", "Acsl5",
    "Fabp1", "Fabp2", "Fabp3", "Fabp4", "Fabp5",
    "Cd36", "Scd1", "Hmgcs2", "Angptl4", "Plin2", "Lpl", "Adipoq",
    "Slc27a1", "Slc27a2", "Slc27a5",
    "Cyp4a10", "Cyp4a12a", "Cyp4a12b", "Cyp4a14",
}

FEATURE_SYMBOLS = {
    "PPAR_signaling": PPAR_SYMBOLS,
    "cholesterol_biosynthesis": {
        "Hmgcr", "Hmgcs1", "Hmgcs2", "Mvk", "Pmvk", "Mvd", "Idi1", "Idi2",
        "Fdps", "Fdft1", "Sqle", "Lss", "Cyp51", "Msmo1", "Nsdhl",
        "Hsd17b7", "Ebp", "Sc5d", "Dhcr7", "Dhcr24", "Srebf1", "Srebf2",
        "Insig1", "Insig2",
    },
    "ECM_remodeling": {
        "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2", "Col5a1",
        "Col6a1", "Fn1", "Lama1", "Lama2", "Lamb1", "Lamb2", "Mmp2",
        "Mmp9", "Mmp14", "Timp1", "Timp2", "Lox", "Sparc", "Vwf",
        "Fbln5", "Itga1", "Itga5", "Itgb1",
    },
    "EMT_fibrosis": {
        "Tgfb1", "Tgfb2", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3", "Smad4",
        "Smad7", "Acta2", "Vim", "Fn1", "Ctgf", "Lox", "Col1a1",
        "Col1a2", "Col3a1", "Serpine1", "Snai1", "Snai2", "Twist1",
        "Zeb1", "Zeb2",
    },
    "tubular_ion_transport": {
        "Slc12a1", "Slc12a3", "Slc9a3", "Slc4a1", "Slc4a4", "Slc5a2",
        "Slc22a2", "Slc22a6", "Slc22a7", "Slc22a8", "Slc34a1", "Slc34a3",
        "Kcnj1", "Kcnj10", "Kcnj16", "Kcnma1", "Atp1a1", "Atp1b1",
        "Fxyd2", "Clcnkb", "Scnn1a", "Scnn1b", "Scnn1g", "Aqp1", "Aqp2",
        "Aqp3", "Trpv5", "Trpv6", "Calb1", "Pvalb",
    },
    "TGF_beta_Wnt": {
        "Tgfb1", "Tgfb2", "Tgfb3", "Tgfbr1", "Tgfbr2", "Smad2", "Smad3",
        "Smad4", "Smad7", "Wnt4", "Wnt5a", "Wnt7b", "Wnt9b", "Wnt11",
        "Fzd1", "Fzd2", "Fzd4", "Fzd7", "Lrp5", "Lrp6", "Dvl2", "Ctnnb1",
        "Axin2", "Lef1", "Tcf7l2",
    },
    "oxidative_stress": {
        "Sod1", "Sod2", "Sod3", "Cat", "Gpx1", "Gpx3", "Gpx4", "Prdx1",
        "Prdx2", "Prdx3", "Nqo1", "Nfe2l2", "Hmox1", "Keap1", "Txn1",
        "Txnrd1", "Gsr", "Gclc", "Gclm",
    },
}


STUDY_SPECS = {
    "OSD-102": {
        "analysis": "lar_young",
        "vst": "GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "vst_glob": "*VST_Counts*.csv",
        "metadata": "OSD-102_metadata_OSD-102-ISA.zip",
        "metadata_glob": "*ISA.zip",
        "flt_pattern": r"_FLT_",
        "gc_pattern": r"_GC_",
    },
    "OSD-513": {
        "analysis": "sex_robustness",
        "vst": "GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "vst_glob": "*VST_Counts*.csv",
        "metadata": "OSD-513_metadata_OSD-513-ISA.zip",
        "metadata_glob": "*ISA.zip",
        "flt_pattern": r"_FLT_",
        "gc_pattern": r"_GC_",
    },
    "OSD-163": {
        "analysis": "cross_strain_context",
        "vst": "GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "vst_glob": "*VST_Counts*.csv",
        "metadata": "OSD-163_metadata_GLDS-163-ISA.zip",
        "metadata_glob": "*ISA.zip",
        "flt_pattern": r"_FLT_",
        "gc_pattern": r"_GC_",
        "exclude_pattern": r"_BSL_|_VIV_",
    },
    "OSD-253": {
        "analysis": "rr7_context",
        "vst": "GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "vst_glob": "*VST_Counts*.csv",
        "metadata": "GLDS-253_metadata_GLDS-253-ISA.zip",
        "metadata_glob": "*ISA.zip",
        "flt_pattern": r"_FLT_",
        "gc_pattern": r"_GC_",
        "exclude_pattern": r"_BSL_|_VIV_|GCrerun|C3H|C3H-HeJ",
    },
}


def load_vst_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing external VST matrix: {path}")
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    if df.index.duplicated().any():
        df = df.groupby(df.index).mean()
    return df


def resolve_study_file(study_dir: Path, exact_name: str, glob_pattern: str, required: bool = True) -> Path | None:
    exact = study_dir / exact_name
    if exact.exists():
        return exact
    matches = sorted(study_dir.glob(glob_pattern))
    if matches:
        return matches[0]
    if required:
        raise FileNotFoundError(
            f"Missing required external file in {study_dir}: expected {exact_name} or {glob_pattern}"
        )
    return None


def infer_flt_gc_samples(
    columns: list[str],
    flt_pattern: str,
    gc_pattern: str,
    include_pattern: str | None = None,
    exclude_pattern: str | None = None,
) -> tuple[list[str], list[str]]:
    filtered = list(columns)
    if include_pattern:
        include_re = re.compile(include_pattern)
        filtered = [c for c in filtered if include_re.search(c)]
    if exclude_pattern:
        exclude_re = re.compile(exclude_pattern)
        filtered = [c for c in filtered if not exclude_re.search(c)]
    flt_re = re.compile(flt_pattern)
    gc_re = re.compile(gc_pattern)
    flt = [c for c in filtered if flt_re.search(c)]
    gc = [c for c in filtered if gc_re.search(c)]
    if len(flt) < 2 or len(gc) < 2:
        raise ValueError(
            f"Need at least 2 FLT and 2 GC samples after filters; found FLT={len(flt)}, GC={len(gc)}. "
            f"include_pattern={include_pattern!r}, exclude_pattern={exclude_pattern!r}"
        )
    return flt, gc


def feature_gene_ids(feature: str, id_map: pd.DataFrame, expression_genes: set[str]) -> list[str]:
    symbol_to_ids: dict[str, set[str]] = {}
    for _, row in id_map.iterrows():
        symbol_to_ids.setdefault(str(row["mgi_symbol"]).casefold(), set()).add(str(row["ensembl_gene_id"]))

    if feature == "translation_machinery":
        symbols = {
            str(sym)
            for sym in id_map["mgi_symbol"].dropna()
            if re.match(r"^(Rpl|Rps|Mrpl|Mrps)\d", str(sym))
        }
    elif feature in FEATURE_SYMBOLS:
        symbols = FEATURE_SYMBOLS[feature]
    else:
        raise ValueError(f"No external-validation gene-set definition for registered feature: {feature}")

    matched: set[str] = set()
    for symbol in symbols:
        matched |= symbol_to_ids.get(str(symbol).casefold(), set())
    return sorted(matched & expression_genes)


def welch_t(values: np.ndarray, labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return per-gene FLT-GC mean difference and Welch t statistic."""
    flt = labels == 1
    gc = labels == 0
    xf = values[:, flt]
    xg = values[:, gc]
    diff = xf.mean(axis=1) - xg.mean(axis=1)
    vf = xf.var(axis=1, ddof=1)
    vg = xg.var(axis=1, ddof=1)
    denom = np.sqrt(vf / xf.shape[1] + vg / xg.shape[1])
    t = np.divide(diff, denom, out=np.zeros_like(diff), where=denom > 1e-12)
    return diff, t


def feature_gene_sets(
    registry: pd.DataFrame,
    id_map: pd.DataFrame,
    expression_genes: set[str],
    min_genes: int = 3,
) -> dict[str, list[str]]:
    feature_to_genes = {
        row.feature: feature_gene_ids(str(row.feature), id_map, expression_genes)
        for row in registry.itertuples()
    }
    missing = {feature: genes for feature, genes in feature_to_genes.items() if len(genes) < min_genes}
    if missing:
        raise ValueError(f"Registered features with <{min_genes} matched genes in external matrix: {missing}")
    return feature_to_genes


def oriented_statistic(stat: float, discovery_effect: float) -> tuple[float, bool]:
    directional = np.sign(discovery_effect) != 0
    if not directional:
        return abs(stat), False
    sign = 1.0 if discovery_effect >= 0 else -1.0
    return sign * stat, True


def gsea_enrichment_score(scores_desc: np.ndarray, hit_mask_desc: np.ndarray, weight: float = 1.0) -> float:
    """Classic preranked GSEA running-sum statistic on descending scores."""
    hit_mask_desc = hit_mask_desc.astype(bool)
    n_genes = len(scores_desc)
    n_hits = int(hit_mask_desc.sum())
    n_misses = n_genes - n_hits
    if n_hits == 0 or n_misses == 0:
        return 0.0

    hit_weights = np.power(np.abs(scores_desc[hit_mask_desc]), weight)
    hit_norm = float(hit_weights.sum())
    running_steps = np.empty(n_genes, dtype=float)
    if hit_norm > 1e-12:
        running_steps[hit_mask_desc] = np.power(np.abs(scores_desc[hit_mask_desc]), weight) / hit_norm
    else:
        running_steps[hit_mask_desc] = 1.0 / n_hits
    running_steps[~hit_mask_desc] = -1.0 / n_misses
    running = np.cumsum(running_steps)
    max_es = float(running.max())
    min_es = float(running.min())
    return max_es if abs(max_es) >= abs(min_es) else min_es


def mean_t_permutation_table(
    vst: pd.DataFrame,
    flt_samples: list[str],
    gc_samples: list[str],
    registry: pd.DataFrame,
    id_map: pd.DataFrame,
    k_perm: int,
    seed: int,
) -> pd.DataFrame:
    """Build external pathway feature table from mean pathway t statistics."""
    sample_order = flt_samples + gc_samples
    labels = np.array([1] * len(flt_samples) + [0] * len(gc_samples), dtype=int)
    expression_genes = set(vst.index.astype(str))

    feature_to_genes = feature_gene_sets(registry, id_map, expression_genes)

    union_genes = sorted(set().union(*feature_to_genes.values()))
    expr = vst.loc[union_genes, sample_order].values.astype(float)
    gene_index = {g: i for i, g in enumerate(union_genes)}
    feature_masks = {
        feature: np.array([gene_index[g] for g in genes], dtype=int)
        for feature, genes in feature_to_genes.items()
    }

    obs_diff, obs_t = welch_t(expr, labels)
    rng = np.random.default_rng(seed)
    exceed = {feature: 0 for feature in feature_masks}
    obs_oriented = {}
    rows = []

    for row in registry.itertuples():
        feature = str(row.feature)
        idx = feature_masks[feature]
        discovery_effect = float(row.discovery_effect)
        mean_effect = float(obs_diff[idx].mean())
        mean_t = float(obs_t[idx].mean())
        oriented, directional = oriented_statistic(mean_t, discovery_effect)
        obs_oriented[feature] = oriented
        rows.append({
            "feature": feature,
            "effect": mean_effect,
            "mean_t": mean_t,
            "gsea_es": np.nan,
            "oriented_stat": oriented,
            "directional_claim": bool(directional),
            "n_genes_in_feature": int(len(idx)),
            "n_flt": int(len(flt_samples)),
            "n_gc": int(len(gc_samples)),
        })

    for _ in range(k_perm):
        perm = labels.copy()
        rng.shuffle(perm)
        _, perm_t = welch_t(expr, perm)
        for row in registry.itertuples():
            feature = str(row.feature)
            idx = feature_masks[feature]
            discovery_effect = float(row.discovery_effect)
            mean_perm_t = float(perm_t[idx].mean())
            null_stat, _ = oriented_statistic(mean_perm_t, discovery_effect)
            exceed[feature] += int(null_stat >= obs_oriented[feature])

    result = pd.DataFrame(rows)
    result["p_value"] = [
        (1.0 + exceed[feature]) / (k_perm + 1.0)
        for feature in result["feature"]
    ]
    result["q_value"] = bh_fdr(result["p_value"].values)
    result["k_permutations"] = int(k_perm)
    result["pathway_method"] = "mean_t"
    return result


def gsea_permutation_table(
    vst: pd.DataFrame,
    flt_samples: list[str],
    gc_samples: list[str],
    registry: pd.DataFrame,
    id_map: pd.DataFrame,
    k_perm: int,
    seed: int,
    weight: float = 1.0,
) -> pd.DataFrame:
    """Build external pathway feature table from preranked GSEA label permutations."""
    sample_order = flt_samples + gc_samples
    labels = np.array([1] * len(flt_samples) + [0] * len(gc_samples), dtype=int)
    mapped_genes = set(id_map["ensembl_gene_id"].astype(str))
    expression_genes = set(vst.index.astype(str)) & mapped_genes
    feature_to_genes = feature_gene_sets(registry, id_map, expression_genes)

    universe_genes = sorted(expression_genes)
    gene_index = {gene: i for i, gene in enumerate(universe_genes)}
    expr = vst.loc[universe_genes, sample_order].values.astype(float)
    feature_masks_by_gene_order = {
        feature: np.array([gene in set(genes) for gene in universe_genes], dtype=bool)
        for feature, genes in feature_to_genes.items()
    }

    obs_diff, obs_t = welch_t(expr, labels)
    rng = np.random.default_rng(seed)
    exceed = {feature: 0 for feature in feature_masks_by_gene_order}
    obs_oriented = {}
    rows = []

    obs_order = np.argsort(-obs_t, kind="mergesort")
    obs_scores_desc = obs_t[obs_order]

    for row in registry.itertuples():
        feature = str(row.feature)
        discovery_effect = float(row.discovery_effect)
        hit_mask_desc = feature_masks_by_gene_order[feature][obs_order]
        es = gsea_enrichment_score(obs_scores_desc, hit_mask_desc, weight=weight)
        oriented, directional = oriented_statistic(es, discovery_effect)
        genes = feature_to_genes[feature]
        idx = np.array([gene_index[g] for g in genes], dtype=int)
        obs_oriented[feature] = oriented
        rows.append({
            "feature": feature,
            "effect": float(obs_diff[idx].mean()),
            "mean_t": float(obs_t[idx].mean()),
            "gsea_es": float(es),
            "oriented_stat": oriented,
            "directional_claim": bool(directional),
            "n_genes_in_feature": int(len(genes)),
            "n_genes_in_universe": int(len(universe_genes)),
            "n_flt": int(len(flt_samples)),
            "n_gc": int(len(gc_samples)),
        })

    for _ in range(k_perm):
        perm = labels.copy()
        rng.shuffle(perm)
        _, perm_t = welch_t(expr, perm)
        perm_order = np.argsort(-perm_t, kind="mergesort")
        perm_scores_desc = perm_t[perm_order]
        for row in registry.itertuples():
            feature = str(row.feature)
            discovery_effect = float(row.discovery_effect)
            hit_mask_desc = feature_masks_by_gene_order[feature][perm_order]
            es = gsea_enrichment_score(perm_scores_desc, hit_mask_desc, weight=weight)
            null_stat, _ = oriented_statistic(es, discovery_effect)
            exceed[feature] += int(null_stat >= obs_oriented[feature])

    result = pd.DataFrame(rows)
    result["p_value"] = [
        (1.0 + exceed[feature]) / (k_perm + 1.0)
        for feature in result["feature"]
    ]
    result["q_value"] = bh_fdr(result["p_value"].values)
    result["k_permutations"] = int(k_perm)
    result["pathway_method"] = "gsea"
    result["gsea_weight"] = float(weight)
    return result


def pathway_permutation_table(
    vst: pd.DataFrame,
    flt_samples: list[str],
    gc_samples: list[str],
    registry: pd.DataFrame,
    id_map: pd.DataFrame,
    k_perm: int,
    seed: int,
    method: str = "gsea",
) -> pd.DataFrame:
    """Build external pathway feature table using the selected pathway statistic."""
    if method == "mean_t":
        return mean_t_permutation_table(vst, flt_samples, gc_samples, registry, id_map, k_perm, seed)
    if method == "gsea":
        return gsea_permutation_table(vst, flt_samples, gc_samples, registry, id_map, k_perm, seed)
    raise ValueError(f"Unknown external pathway method: {method}")


def run_study(
    study: str,
    external_root: Path,
    protocol_dir: Path,
    id_map_path: Path,
    outdir: Path,
    k_perm: int,
    seed: int,
    method: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = STUDY_SPECS[study]
    analysis = spec["analysis"]
    validate_study_scope(study, analysis)

    _, hypotheses = require_protocol(protocol_dir)
    registry = registry_for_analysis(hypotheses, study, analysis)
    study_dir = external_root / study
    vst_path = resolve_study_file(study_dir, spec["vst"], spec.get("vst_glob", "*VST_Counts*.csv"))
    metadata_path = resolve_study_file(
        study_dir,
        spec["metadata"],
        spec.get("metadata_glob", "*ISA.zip"),
        required=False,
    )
    vst = load_vst_matrix(vst_path)
    flt_samples, gc_samples = infer_flt_gc_samples(
        list(vst.columns),
        spec["flt_pattern"],
        spec["gc_pattern"],
        include_pattern=spec.get("include_pattern"),
        exclude_pattern=spec.get("exclude_pattern"),
    )
    id_map = load_id_map(id_map_path)

    features = pathway_permutation_table(
        vst=vst,
        flt_samples=flt_samples,
        gc_samples=gc_samples,
        registry=registry,
        id_map=id_map,
        k_perm=k_perm,
        seed=seed,
        method=method,
    )
    discovery = registry.rename(columns={"discovery_effect": "effect"})[
        ["feature", "effect", "q_threshold", "expected_direction", "hypothesis_type", "claim_boundary"]
    ].copy()
    replication = compare_directional_replication(discovery, features)
    replication["study"] = study
    replication["analysis"] = analysis
    replication["pathway_method"] = method

    outdir.mkdir(parents=True, exist_ok=True)
    features.to_csv(outdir / f"{study}_{analysis}_external_features.tsv", sep="\t", index=False)
    replication.to_csv(outdir / f"{study}_{analysis}_replication.tsv", sep="\t", index=False)
    manifest = {
        "study": study,
        "analysis": analysis,
        "vst": str(vst_path),
        "metadata": str(metadata_path) if metadata_path is not None else None,
        "n_flt": len(flt_samples),
        "n_gc": len(gc_samples),
        "flt_samples": flt_samples,
        "gc_samples": gc_samples,
        "include_pattern": spec.get("include_pattern"),
        "exclude_pattern": spec.get("exclude_pattern"),
        "k_permutations": k_perm,
        "seed": seed,
        "pathway_method": method,
        "claim_boundary": "independent expression-level pathway validation/context mapping; no cohort pooling",
    }
    (outdir / f"{study}_{analysis}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return features, replication


def main() -> None:
    ap = argparse.ArgumentParser(description="Run protocol-bound OSD external validation/context mapping")
    ap.add_argument("--external_root", default=str(REPO_ROOT / "data/external/osdr"))
    ap.add_argument("--protocol_dir", default=str(REPO_ROOT / "docs/external_replication_protocol"))
    ap.add_argument("--id_map", default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"))
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/external_validation"))
    ap.add_argument("--studies", default="auto",
                    help="Comma-separated studies, or auto for all supported studies with downloaded folders.")
    ap.add_argument("--K_perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--method", choices=["gsea", "mean_t"], default="gsea",
                    help="External pathway statistic. gsea is preranked GSEA over all mapped genes.")
    args = ap.parse_args()

    external_root = Path(args.external_root)
    if args.studies.strip().lower() == "auto":
        studies = [study for study in STUDY_SPECS if (external_root / study).exists()]
    else:
        studies = [s.strip() for s in args.studies.split(",") if s.strip()]
    unknown = sorted(set(studies) - set(STUDY_SPECS))
    if unknown:
        raise ValueError(f"Unknown/unsupported external studies: {unknown}")
    if not studies:
        raise ValueError(f"No supported external study folders found under {external_root}")

    all_replication = []
    for offset, study in enumerate(studies):
        print(f"\n[external] {study}")
        _, replication = run_study(
            study=study,
            external_root=Path(args.external_root),
            protocol_dir=Path(args.protocol_dir),
            id_map_path=Path(args.id_map),
            outdir=Path(args.outdir),
            k_perm=args.K_perm,
            seed=args.seed + offset,
            method=args.method,
        )
        all_replication.append(replication)

    summary = pd.concat(all_replication, ignore_index=True) if all_replication else pd.DataFrame()
    summary.to_csv(Path(args.outdir) / "external_validation_summary.tsv", sep="\t", index=False)
    print(f"\n[OK] External validation outputs written to {args.outdir}")


if __name__ == "__main__":
    main()
