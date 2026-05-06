"""
OSD-102/OSD-513 protocol-bound external validation.

This module builds compact external feature tables from the downloaded GeneLab
processed VST matrices and immediately evaluates them through
src.validation.external_replication. It is intentionally separate from
multi-study pooling: each OSD cohort is analyzed independently.
"""

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


STUDY_SPECS = {
    "OSD-102": {
        "analysis": "lar_young",
        "vst": "GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "metadata": "OSD-102_metadata_OSD-102-ISA.zip",
        "flt_pattern": r"_FLT_",
        "gc_pattern": r"_GC_",
    },
    "OSD-513": {
        "analysis": "sex_robustness",
        "vst": "GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
        "metadata": "OSD-513_metadata_OSD-513-ISA.zip",
        "flt_pattern": r"_FLT_",
        "gc_pattern": r"_GC_",
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


def infer_flt_gc_samples(columns: list[str], flt_pattern: str, gc_pattern: str) -> tuple[list[str], list[str]]:
    flt_re = re.compile(flt_pattern)
    gc_re = re.compile(gc_pattern)
    flt = [c for c in columns if flt_re.search(c)]
    gc = [c for c in columns if gc_re.search(c)]
    if len(flt) < 2 or len(gc) < 2:
        raise ValueError(f"Need at least 2 FLT and 2 GC samples; found FLT={len(flt)}, GC={len(gc)}")
    return flt, gc


def feature_gene_ids(feature: str, id_map: pd.DataFrame, expression_genes: set[str]) -> list[str]:
    symbol_to_ids: dict[str, set[str]] = {}
    for _, row in id_map.iterrows():
        symbol_to_ids.setdefault(str(row["mgi_symbol"]), set()).add(str(row["ensembl_gene_id"]))

    if feature == "PPAR_signaling":
        symbols = PPAR_SYMBOLS
    elif feature == "translation_machinery":
        symbols = {
            str(sym)
            for sym in id_map["mgi_symbol"].dropna()
            if re.match(r"^(Rpl|Rps|Mrpl|Mrps)\d", str(sym))
        }
    else:
        raise ValueError(f"No external-validation gene-set definition for registered feature: {feature}")

    matched: set[str] = set()
    for symbol in symbols:
        matched |= symbol_to_ids.get(symbol, set())
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


def pathway_permutation_table(
    vst: pd.DataFrame,
    flt_samples: list[str],
    gc_samples: list[str],
    registry: pd.DataFrame,
    id_map: pd.DataFrame,
    k_perm: int,
    seed: int,
) -> pd.DataFrame:
    """Build external pathway feature table using sample-label permutations."""
    sample_order = flt_samples + gc_samples
    labels = np.array([1] * len(flt_samples) + [0] * len(gc_samples), dtype=int)
    expression_genes = set(vst.index.astype(str))

    feature_to_genes = {
        row.feature: feature_gene_ids(str(row.feature), id_map, expression_genes)
        for row in registry.itertuples()
    }
    missing = {feature: genes for feature, genes in feature_to_genes.items() if len(genes) < 3}
    if missing:
        raise ValueError(f"Registered features with <3 matched genes in external matrix: {missing}")

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
        sign = 1.0 if float(row.discovery_effect) >= 0 else -1.0
        mean_effect = float(obs_diff[idx].mean())
        mean_t = float(obs_t[idx].mean())
        oriented = sign * mean_t
        obs_oriented[feature] = oriented
        rows.append({
            "feature": feature,
            "effect": mean_effect,
            "mean_t": mean_t,
            "oriented_stat": oriented,
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
            sign = 1.0 if float(row.discovery_effect) >= 0 else -1.0
            null_stat = sign * float(perm_t[idx].mean())
            exceed[feature] += int(null_stat >= obs_oriented[feature])

    result = pd.DataFrame(rows)
    result["p_value"] = [
        (1.0 + exceed[feature]) / (k_perm + 1.0)
        for feature in result["feature"]
    ]
    result["q_value"] = bh_fdr(result["p_value"].values)
    result["k_permutations"] = int(k_perm)
    return result


def run_study(
    study: str,
    external_root: Path,
    protocol_dir: Path,
    id_map_path: Path,
    outdir: Path,
    k_perm: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = STUDY_SPECS[study]
    analysis = spec["analysis"]
    validate_study_scope(study, analysis)

    _, hypotheses = require_protocol(protocol_dir)
    registry = registry_for_analysis(hypotheses, study, analysis)
    study_dir = external_root / study
    vst = load_vst_matrix(study_dir / spec["vst"])
    flt_samples, gc_samples = infer_flt_gc_samples(
        list(vst.columns),
        spec["flt_pattern"],
        spec["gc_pattern"],
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
    )
    discovery = registry.rename(columns={"discovery_effect": "effect"})[
        ["feature", "effect", "q_threshold", "expected_direction", "hypothesis_type", "claim_boundary"]
    ].copy()
    replication = compare_directional_replication(discovery, features)
    replication["study"] = study
    replication["analysis"] = analysis

    outdir.mkdir(parents=True, exist_ok=True)
    features.to_csv(outdir / f"{study}_{analysis}_external_features.tsv", sep="\t", index=False)
    replication.to_csv(outdir / f"{study}_{analysis}_replication.tsv", sep="\t", index=False)
    manifest = {
        "study": study,
        "analysis": analysis,
        "vst": str(study_dir / spec["vst"]),
        "metadata": str(study_dir / spec["metadata"]),
        "n_flt": len(flt_samples),
        "n_gc": len(gc_samples),
        "flt_samples": flt_samples,
        "gc_samples": gc_samples,
        "k_permutations": k_perm,
        "seed": seed,
        "claim_boundary": "independent expression-level pathway validation; no cohort pooling",
    }
    (outdir / f"{study}_{analysis}_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return features, replication


def main() -> None:
    ap = argparse.ArgumentParser(description="Run OSD-102/OSD-513 external validation")
    ap.add_argument("--external_root", default=str(REPO_ROOT / "data/external/osdr"))
    ap.add_argument("--protocol_dir", default=str(REPO_ROOT / "docs/external_replication_protocol"))
    ap.add_argument("--id_map", default=str(REPO_ROOT / "data/processed/resources/id_map.tsv"))
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/external_validation"))
    ap.add_argument("--studies", default="OSD-102,OSD-513")
    ap.add_argument("--K_perm", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    studies = [s.strip() for s in args.studies.split(",") if s.strip()]
    unknown = sorted(set(studies) - set(STUDY_SPECS))
    if unknown:
        raise ValueError(f"Unknown/unsupported external studies: {unknown}")

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
        )
        all_replication.append(replication)

    summary = pd.concat(all_replication, ignore_index=True) if all_replication else pd.DataFrame()
    summary.to_csv(Path(args.outdir) / "external_validation_summary.tsv", sep="\t", index=False)
    print(f"\n[OK] External validation outputs written to {args.outdir}")


if __name__ == "__main__":
    main()
