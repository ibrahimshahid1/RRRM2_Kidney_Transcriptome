#!/usr/bin/env python3
"""Reach E -- hypothesis-generating LINCS/CMap appendix screen.

This is not a treatment-discovery analysis.  It builds a conservative human
L1000-compatible query from the Repair B cross-cohort mouse RNA meta-analysis
and, when the required LINCS metadata are present, computes an approximate local
connectivity score against the Level-5 matrix.

Critical guardrail: local signature scores are not interpretable without
``sig_info`` because the signature id alone does not tell us the perturbagen,
cell line, dose, or time.  If ``sig_info`` is missing, this module writes only
``query_genes.tsv`` plus a verdict marking the local screen as blocked.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

from src.v11.core_analysis import RUN_ROOT, sha256


LINCS_DIR = Path("data/external/lincs_cmap")
GENE_INFO = LINCS_DIR / "GSE92742_Broad_LINCS_gene_info.txt.gz"
SIG_INFO = LINCS_DIR / "GSE92742_Broad_LINCS_sig_info.txt.gz"
PERT_INFO = LINCS_DIR / "GSE92742_Broad_LINCS_pert_info.txt.gz"
LEVEL5 = LINCS_DIR / "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
REPAIR_B_META = Path("data/results/run_20260601_repair_b/cross_osdr_recurrence/recurrence_meta_summary.tsv")

QUERY_N_PER_DIRECTION = 50
DEFAULT_CHUNK_SIZE = 4096
DEFAULT_MAX_SIGNATURES_OUTPUT = 5000


def read_tsv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False)


def load_lincs_gene_info(path: str | Path = GENE_INFO) -> pd.DataFrame:
    gene_info = read_tsv(path)
    required = {"pr_gene_id", "pr_gene_symbol", "pr_is_lm", "pr_is_bing"}
    missing = required - set(gene_info.columns)
    if missing:
        raise ValueError(f"{path} missing required LINCS gene columns: {sorted(missing)}")
    out = gene_info.copy()
    out["human_entrez_id"] = out["pr_gene_id"].astype(str)
    out["human_symbol"] = out["pr_gene_symbol"].astype(str).str.upper()
    out["is_landmark"] = pd.to_numeric(out["pr_is_lm"], errors="coerce").fillna(0).astype(int).astype(bool)
    out["is_bing"] = pd.to_numeric(out["pr_is_bing"], errors="coerce").fillna(0).astype(int).astype(bool)
    out["lincs_space"] = np.where(out["is_landmark"], "landmark", np.where(out["is_bing"], "bing", "aig"))
    return out.drop_duplicates("human_symbol")


def load_ortholog_map(path: str | Path | None) -> dict[str, str]:
    """Load an optional mouse->human symbol map.

    Expected columns are flexible but should contain mouse and human symbol
    names.  If no path is supplied or the file is absent, an empty map is
    returned and the query builder falls back to conservative uppercase symbol
    matching.
    """
    if path is None:
        return {}
    p = Path(path)
    if not p.exists():
        return {}
    tab = read_tsv(p)
    lower = {c.lower(): c for c in tab.columns}
    mouse_col = lower.get("mouse_symbol") or lower.get("mgi_symbol") or lower.get("mouse_gene_symbol")
    human_col = lower.get("human_symbol") or lower.get("hgnc_symbol") or lower.get("human_gene_symbol")
    if not mouse_col or not human_col:
        raise ValueError(f"{p} must contain mouse_symbol/mgi_symbol and human_symbol/hgnc_symbol columns")
    return {
        str(m).strip().lower(): str(h).strip().upper()
        for m, h in zip(tab[mouse_col], tab[human_col])
        if str(m).strip() and str(h).strip()
    }


def mouse_to_human_symbol(mouse_symbol: object, orthologs: Mapping[str, str] | None = None) -> tuple[str, str]:
    sym = str(mouse_symbol).strip()
    if not sym or sym.lower() == "nan":
        return "", "missing_symbol"
    if orthologs and sym.lower() in orthologs:
        return str(orthologs[sym.lower()]).upper(), "ortholog_table"
    return sym.upper(), "uppercase_symbol_assumption"


def build_query_genes(
    meta: pd.DataFrame,
    gene_info: pd.DataFrame,
    *,
    orthologs: Mapping[str, str] | None = None,
    n_per_direction: int = QUERY_N_PER_DIRECTION,
    lincs_space: str = "bing",
) -> pd.DataFrame:
    """Build top up/down human LINCS query genes from Repair B meta-analysis."""
    required = {"gene_symbol", "pooled_effect"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"meta table missing required columns: {sorted(missing)}")
    ginfo = gene_info.copy()
    if lincs_space == "landmark":
        valid = ginfo[ginfo["is_landmark"]]
    elif lincs_space == "bing":
        valid = ginfo[ginfo["is_bing"] | ginfo["is_landmark"]]
    else:
        valid = ginfo
    valid_by_symbol = valid.set_index("human_symbol")

    rows: list[dict[str, object]] = []
    for _, row in meta.iterrows():
        effect = pd.to_numeric(row.get("pooled_effect"), errors="coerce")
        if not np.isfinite(effect) or effect == 0:
            continue
        mouse = row.get("gene_symbol", "")
        human, method = mouse_to_human_symbol(mouse, orthologs)
        direction = "up" if effect > 0 else "down"
        z = pd.to_numeric(row.get("z", np.nan), errors="coerce")
        score = abs(float(z)) if np.isfinite(z) else abs(float(effect))
        in_lincs = human in valid_by_symbol.index
        info = valid_by_symbol.loc[human] if in_lincs else None
        rows.append(
            {
                "mouse_gene_symbol": mouse,
                "human_gene_symbol": human,
                "mapping_method": method,
                "direction": direction,
                "pooled_effect": float(effect),
                "z": None if not np.isfinite(z) else float(z),
                "fdr": row.get("fdr", np.nan),
                "rank_score": score,
                "lincs_valid": bool(in_lincs),
                "human_entrez_id": "" if info is None else str(info["human_entrez_id"]),
                "lincs_space": "" if info is None else str(info["lincs_space"]),
                "is_landmark": False if info is None else bool(info["is_landmark"]),
                "is_bing": False if info is None else bool(info["is_bing"]),
                "included_in_query": False,
                "query_rank": np.nan,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.drop_duplicates(["direction", "human_gene_symbol"], keep="first")
    valid_out = out[out["lincs_valid"]].copy()
    for direction in ["up", "down"]:
        idx = (
            valid_out[valid_out["direction"].eq(direction)]
            .sort_values("rank_score", ascending=False)
            .head(n_per_direction)
            .index
        )
        out.loc[idx, "included_in_query"] = True
        out.loc[idx, "query_rank"] = np.arange(1, len(idx) + 1)
    out["query_limit_per_direction"] = int(n_per_direction)
    return out.sort_values(["included_in_query", "direction", "query_rank"], ascending=[False, True, True])


def decode_bytes(values) -> list[str]:
    return [v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values]


def connectivity_score(values: np.ndarray, up_local: np.ndarray, down_local: np.ndarray) -> np.ndarray:
    """Approximate signed CMap score for signature rows.

    Positive means the signature mimics the mouse flight meta-signature
    (up-query genes high and down-query genes low); negative means reversal.
    """
    up = values[:, up_local]
    down = values[:, down_local]
    return np.nanmean(up, axis=1) - np.nanmean(down, axis=1)


def _matrix_layout(h5) -> dict[str, object]:
    matrix = h5["0/DATA/0/matrix"]
    gene_ids = decode_bytes(h5["0/META/ROW/id"][:])
    sig_ids = decode_bytes(h5["0/META/COL/id"][:])
    if matrix.shape[0] == len(sig_ids) and matrix.shape[1] == len(gene_ids):
        return {"sig_axis": 0, "gene_axis": 1, "matrix": matrix, "gene_ids": gene_ids, "sig_ids": sig_ids}
    if matrix.shape[0] == len(gene_ids) and matrix.shape[1] == len(sig_ids):
        return {"sig_axis": 1, "gene_axis": 0, "matrix": matrix, "gene_ids": gene_ids, "sig_ids": sig_ids}
    raise ValueError(
        f"Cannot reconcile GCTx matrix shape {matrix.shape} with {len(gene_ids)} genes and {len(sig_ids)} signatures"
    )


def score_level5_matrix(
    level5_path: str | Path,
    query: pd.DataFrame,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> pd.DataFrame:
    """Score local LINCS Level-5 signatures against the query."""
    import h5py

    q = query[query["included_in_query"].astype(bool)].copy()
    up_ids = set(q[q["direction"].eq("up")]["human_entrez_id"].astype(str))
    down_ids = set(q[q["direction"].eq("down")]["human_entrez_id"].astype(str))
    if not up_ids or not down_ids:
        raise ValueError("CMap query requires at least one up and one down LINCS gene")

    with h5py.File(level5_path, "r") as h5:
        layout = _matrix_layout(h5)
        gene_to_idx = {gid: i for i, gid in enumerate(layout["gene_ids"])}
        up_idx = [gene_to_idx[g] for g in up_ids if g in gene_to_idx]
        down_idx = [gene_to_idx[g] for g in down_ids if g in gene_to_idx]
        if not up_idx or not down_idx:
            raise ValueError("Query genes are absent from Level-5 matrix gene ids")

        all_idx = sorted(set(up_idx + down_idx))
        local = {gidx: i for i, gidx in enumerate(all_idx)}
        up_local = np.array([local[i] for i in up_idx], dtype=int)
        down_local = np.array([local[i] for i in down_idx], dtype=int)

        matrix = layout["matrix"]
        sig_ids = layout["sig_ids"]
        scores: list[np.ndarray] = []
        sig_chunks: list[list[str]] = []
        n_sig = len(sig_ids)
        for start in range(0, n_sig, chunk_size):
            end = min(start + chunk_size, n_sig)
            if layout["sig_axis"] == 0:
                block = matrix[start:end, all_idx]
            else:
                block = matrix[all_idx, start:end].T
            block = np.asarray(block, dtype=float)
            scores.append(connectivity_score(block, up_local, down_local))
            sig_chunks.append(sig_ids[start:end])

    return pd.DataFrame(
        {
            "sig_id": [s for chunk in sig_chunks for s in chunk],
            "connectivity_score": np.concatenate(scores),
            "n_query_up": len(up_idx),
            "n_query_down": len(down_idx),
        }
    )


def load_sig_info(path: str | Path = SIG_INFO) -> pd.DataFrame:
    return read_tsv(path)


def load_pert_info(path: str | Path = PERT_INFO) -> pd.DataFrame:
    return read_tsv(path)


def merge_signature_metadata(
    scores: pd.DataFrame,
    sig_info: pd.DataFrame,
    pert_info: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if "sig_id" not in sig_info.columns:
        raise ValueError("sig_info is missing sig_id")
    out = scores.merge(sig_info, on="sig_id", how="left")
    if pert_info is not None and "pert_id" in out.columns and "pert_id" in pert_info.columns:
        keep_cols = [c for c in ["pert_id", "moa", "target", "pert_iname", "canonical_smiles"] if c in pert_info.columns]
        if keep_cols:
            suffix_cols = [c for c in keep_cols if c != "pert_id" and c in out.columns]
            out = out.merge(pert_info[keep_cols].drop_duplicates("pert_id"), on="pert_id", how="left", suffixes=("", "_pert"))
            for c in suffix_cols:
                pc = f"{c}_pert"
                if pc in out.columns:
                    out[c] = out[c].where(out[c].notna(), out[pc])
                    out = out.drop(columns=[pc])
    out["connectivity_class"] = np.where(out["connectivity_score"] < 0, "reversal", "mimic")
    return out


def select_extreme_signatures(scored: pd.DataFrame, max_output: int = DEFAULT_MAX_SIGNATURES_OUTPUT) -> pd.DataFrame:
    if len(scored) <= max_output:
        return scored.sort_values("connectivity_score").reset_index(drop=True)
    n_each = max(1, max_output // 2)
    low = scored.nsmallest(n_each, "connectivity_score").index
    high = scored.nlargest(n_each, "connectivity_score").index
    idx = sorted(set(low) | set(high))
    return scored.loc[idx].sort_values("connectivity_score").reset_index(drop=True)


def summarize_perturbagens(scored: pd.DataFrame) -> pd.DataFrame:
    pert_col = "pert_iname" if "pert_iname" in scored.columns else "pert_id" if "pert_id" in scored.columns else None
    if pert_col is None:
        return pd.DataFrame([{"status": "missing_perturbagen_metadata"}])
    group_cols = [pert_col]
    for c in ["pert_id", "pert_type"]:
        if c in scored.columns and c not in group_cols:
            group_cols.append(c)
    agg = (
        scored.groupby(group_cols, dropna=False)["connectivity_score"]
        .agg(n_signatures="count", median_connectivity="median", mean_connectivity="mean",
             best_reversal_score="min", best_mimic_score="max")
        .reset_index()
    )
    agg["appendix_interpretation"] = np.where(
        agg["median_connectivity"] < 0,
        "hypothesis_generating_reversal_class",
        "hypothesis_generating_mimic_class",
    )
    return agg.sort_values(["median_connectivity", "n_signatures"], ascending=[True, False]).reset_index(drop=True)


def summarize_moa(scored: pd.DataFrame) -> pd.DataFrame:
    moa_col = "moa" if "moa" in scored.columns else None
    if moa_col is None:
        return pd.DataFrame([{"status": "moa_metadata_absent", "summary": "Download pert_info or use CLUE MoA annotations."}])
    sub = scored[scored[moa_col].notna() & scored[moa_col].astype(str).ne("")]
    if sub.empty:
        return pd.DataFrame([{"status": "no_moa_annotations_in_scored_subset"}])
    out = (
        sub.groupby(moa_col)["connectivity_score"]
        .agg(n_signatures="count", median_connectivity="median", mean_connectivity="mean",
             best_reversal_score="min", best_mimic_score="max")
        .reset_index()
        .rename(columns={moa_col: "moa"})
        .sort_values(["median_connectivity", "n_signatures"], ascending=[True, False])
        .reset_index(drop=True)
    )
    out["appendix_interpretation"] = np.where(
        out["median_connectivity"] < 0,
        "mechanism_class_reversal_hypothesis",
        "mechanism_class_mimic_hypothesis",
    )
    return out


def _resolve_repair_b_meta(root: Path, override: str | Path | None = None) -> Path:
    if override:
        return Path(override)
    current = root / "cross_osdr_recurrence" / "recurrence_meta_summary.tsv"
    return current if current.exists() else REPAIR_B_META


def _write_verdict(out_dir: Path, verdict: dict[str, object]) -> dict[str, object]:
    summary_path = out_dir / "cmap_summary.tsv"
    scalar = {
        k: v
        for k, v in verdict.items()
        if v is None or isinstance(v, (str, int, float, bool, np.integer, np.floating, np.bool_))
    }
    scalar["analysis_key"] = "cmap_screen"
    pd.DataFrame([scalar]).to_csv(summary_path, sep="\t", index=False)
    outputs = verdict.get("outputs")
    if isinstance(outputs, dict):
        outputs["summary"] = str(summary_path)
    else:
        verdict["outputs"] = {"summary": str(summary_path)}
    (out_dir / "cmap_verdict.json").write_text(json.dumps(verdict, indent=2))
    return verdict


def run_cmap_screen(
    root: str | Path,
    *,
    recurrence_meta: str | Path | None = None,
    gene_info_path: str | Path = GENE_INFO,
    sig_info_path: str | Path = SIG_INFO,
    pert_info_path: str | Path = PERT_INFO,
    level5_path: str | Path = LEVEL5,
    ortholog_map: str | Path | None = None,
    n_per_direction: int = QUERY_N_PER_DIRECTION,
    max_signatures_output: int = DEFAULT_MAX_SIGNATURES_OUTPUT,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    root = Path(root)
    out_dir = root / "cmap_screen"
    out_dir.mkdir(parents=True, exist_ok=True)

    meta_path = _resolve_repair_b_meta(root, recurrence_meta)
    if not meta_path.exists():
        return _write_verdict(out_dir, {
            "analysis": "Reach E -- LINCS/CMap appendix screen",
            "status": "blocked_missing_repair_b_meta",
            "missing": str(meta_path),
        })
    if not Path(gene_info_path).exists():
        return _write_verdict(out_dir, {
            "analysis": "Reach E -- LINCS/CMap appendix screen",
            "status": "blocked_missing_lincs_gene_info",
            "missing": str(gene_info_path),
        })

    meta = read_tsv(meta_path)
    gene_info = load_lincs_gene_info(gene_info_path)
    orthologs = load_ortholog_map(ortholog_map)
    query = build_query_genes(
        meta,
        gene_info,
        orthologs=orthologs,
        n_per_direction=n_per_direction,
        lincs_space="bing",
    )
    query_path = out_dir / "query_genes.tsv"
    query.to_csv(query_path, sep="\t", index=False)
    q_in = query[query.get("included_in_query", pd.Series(dtype=bool)).astype(bool)] if not query.empty else query
    n_up = int((q_in["direction"] == "up").sum()) if not q_in.empty else 0
    n_down = int((q_in["direction"] == "down").sum()) if not q_in.empty else 0

    base_verdict = {
        "analysis": "Reach E -- LINCS/CMap appendix screen",
        "claim_discipline": "Hypothesis-generating appendix only; report perturbagen/MoA classes, not treatments.",
        "recurrence_meta": str(meta_path),
        "query_genes": str(query_path),
        "n_query_up": n_up,
        "n_query_down": n_down,
        "ortholog_mapping": "ortholog_table" if orthologs else "uppercase_symbol_assumption",
        "input_manifest": [
            {"path": str(p), "sha256": sha256(p)}
            for p in [Path(meta_path), Path(gene_info_path)]
            if p.exists()
        ],
    }
    if n_up == 0 or n_down == 0:
        base_verdict["status"] = "query_only_insufficient_up_or_down_genes"
        return _write_verdict(out_dir, base_verdict)
    if not Path(sig_info_path).exists():
        base_verdict["status"] = "query_only_missing_sig_info"
        base_verdict["missing"] = str(sig_info_path)
        base_verdict["interpretation"] = "Local LINCS hits are deliberately not scored/interpreted until sig_info is present."
        return _write_verdict(out_dir, base_verdict)
    if not Path(level5_path).exists():
        base_verdict["status"] = "query_only_missing_level5_matrix"
        base_verdict["missing"] = str(level5_path)
        return _write_verdict(out_dir, base_verdict)

    scores = score_level5_matrix(level5_path, query, chunk_size=chunk_size)
    sig_info = load_sig_info(sig_info_path)
    pert_info = load_pert_info(pert_info_path) if Path(pert_info_path).exists() else None
    scored = merge_signature_metadata(scores, sig_info, pert_info)
    selected = select_extreme_signatures(scored, max_output=max_signatures_output)
    pert = summarize_perturbagens(selected)
    moa = summarize_moa(selected)

    sig_path = out_dir / "signature_connectivity.tsv"
    pert_path = out_dir / "perturbagen_summary.tsv"
    moa_path = out_dir / "moa_summary.tsv"
    selected.to_csv(sig_path, sep="\t", index=False)
    pert.to_csv(pert_path, sep="\t", index=False)
    moa.to_csv(moa_path, sep="\t", index=False)

    base_verdict.update(
        {
            "status": "complete_local_approximate",
            "n_signatures_scored": int(len(scores)),
            "n_signatures_output": int(len(selected)),
            "sig_info": str(sig_info_path),
            "pert_info": str(pert_info_path) if Path(pert_info_path).exists() else None,
            "level5": str(level5_path),
            "outputs": {
                "signature_connectivity": str(sig_path),
                "perturbagen_summary": str(pert_path),
                "moa_summary": str(moa_path),
            },
            "interpretation": (
                "Approximate local connectivity only. Prefer official CLUE query for final appendix numbers "
                "if an API/account workflow is available."
            ),
        }
    )
    return _write_verdict(out_dir, base_verdict)


def main() -> None:
    ap = argparse.ArgumentParser(description="Reach E: LINCS/CMap appendix screen.")
    ap.add_argument("--run-root", default=str(RUN_ROOT))
    ap.add_argument("--recurrence-meta", default=None)
    ap.add_argument("--gene-info", default=str(GENE_INFO))
    ap.add_argument("--sig-info", default=str(SIG_INFO))
    ap.add_argument("--pert-info", default=str(PERT_INFO))
    ap.add_argument("--level5", default=str(LEVEL5))
    ap.add_argument("--ortholog-map", default=None)
    ap.add_argument("--n-per-direction", type=int, default=QUERY_N_PER_DIRECTION)
    ap.add_argument("--max-signatures-output", type=int, default=DEFAULT_MAX_SIGNATURES_OUTPUT)
    ap.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    args = ap.parse_args()
    verdict = run_cmap_screen(
        args.run_root,
        recurrence_meta=args.recurrence_meta,
        gene_info_path=args.gene_info,
        sig_info_path=args.sig_info,
        pert_info_path=args.pert_info,
        level5_path=args.level5,
        ortholog_map=args.ortholog_map,
        n_per_direction=args.n_per_direction,
        max_signatures_output=args.max_signatures_output,
        chunk_size=args.chunk_size,
    )
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
