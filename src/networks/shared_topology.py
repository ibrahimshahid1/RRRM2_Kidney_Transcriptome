# src/networks/shared_topology.py
"""
Phase 2 Step 2.1-2.3: Build Shared Sparse Skeleton E

Cell-standardize expression within each (Age × Arm × EnvGroup) cell,
then build a sparse partial correlation network using Ledoit-Wolf
shrinkage and top-k neighbors per gene.

The skeleton E is fixed for all downstream sample-specific weighting.

Usage:
    python -m src.networks.shared_topology --max_genes 2500 --topk 80

With biotype filtering (recommended):
    python -m src.networks.shared_topology --max_genes 2500 --topk 80 \\
        --id_map data/processed/resources/id_map.tsv \\
        --biotype_filter protein_coding
"""
from __future__ import annotations

import argparse
from itertools import combinations
import re
from pathlib import Path

from src.common import REPO_ROOT, load_anchor_config, resolve_configured_genes
import numpy as np
import pandas as pd
from sklearn.covariance import LedoitWolf

DEFAULT_TOPK = 80


def load_rtech(path: str) -> pd.DataFrame:
    """Load Rtech matrix (genes x samples) from gzipped TSV."""
    return pd.read_csv(path, sep="\t", compression="gzip", index_col=0)


def load_meta(path: str) -> pd.DataFrame:
    """Load metadata from gzipped TSV."""
    return pd.read_csv(path, sep="\t", compression="gzip")


def load_biotype_map(id_map_path: str) -> pd.DataFrame:
    """Load Ensembl biotype annotations from id_map.tsv.

    Returns DataFrame with columns: ensembl_gene_id, mgi_symbol, biotype
    """
    df = pd.read_csv(id_map_path, sep="\t", comment="#")
    cols = ["ensembl_gene_id", "mgi_symbol", "biotype"]
    for c in cols:
        if c not in df.columns:
            raise ValueError(f"id_map.tsv missing required column: {c}")
    return df[cols].copy()


def configured_anchor_ids(
    anchor_config: str | Path,
    id_map: str | Path,
    expression_genes: set[str],
    min_anchors: int | None = None,
    warn_threshold: int | None = None,
) -> tuple[list[str], pd.DataFrame]:
    """Resolve configured anchor symbols to Ensembl IDs present in expression."""
    cfg, records = load_anchor_config(anchor_config)
    symbols = [r["symbol"] for r in records]
    resolved = resolve_configured_genes(symbols, id_map, panel_genes=expression_genes)

    # Attach YAML group/role metadata for auditability.
    rec_df = pd.DataFrame(records).rename(columns={"symbol": "query"})
    if not resolved.empty:
        resolved = resolved.merge(rec_df, on="query", how="left")
    else:
        resolved = rec_df.assign(
            query_type="symbol",
            symbol=rec_df["query"],
            ensembl_gene_id="",
            status="unmapped",
            in_panel=False,
        )

    validation = cfg.get("validation", {}) if isinstance(cfg, dict) else {}
    min_anchors = int(min_anchors if min_anchors is not None else validation.get("minimum_anchors", 20))
    warn_threshold = int(warn_threshold if warn_threshold is not None else validation.get("warn_threshold", 50))

    mapped_present = resolved[
        (resolved["status"] == "mapped") &
        (resolved["ensembl_gene_id"] != "") &
        (resolved["in_panel"] == True)
    ]["ensembl_gene_id"].drop_duplicates().tolist()

    if len(mapped_present) < min_anchors:
        missing = resolved[resolved["in_panel"] != True][["query", "status", "ensembl_gene_id"]].head(20)
        raise RuntimeError(
            f"Only {len(mapped_present)} configured anchors are present in the expression universe; "
            f"minimum_anchors={min_anchors}. First missing/unavailable anchors:\n"
            f"{missing.to_string(index=False)}"
        )
    if len(mapped_present) < warn_threshold:
        print(
            f"  WARNING: {len(mapped_present)} configured anchors are present; "
            f"warn_threshold={warn_threshold}. Alignment may be fragile."
        )
    return mapped_present, resolved


# Regex patterns for noise genes to exclude regardless of biotype
_NOISE_SYMBOL_PATTERNS = [
    re.compile(r"^Gm\d+$"),          # predicted genes (Gm12345)
    re.compile(r"Rik$|Rik\d+$"),     # RIKEN cDNA clones
    re.compile(r"^\d+-[A-Z]"),       # numeric-prefix BAC clones
]


def _is_noise_symbol(symbol: str) -> bool:
    """Check if a gene symbol matches known noise patterns."""
    if not isinstance(symbol, str) or not symbol:
        return True  # unmapped genes are noise
    return any(p.search(symbol) for p in _NOISE_SYMBOL_PATTERNS)


def pick_genes(
    rtech: pd.DataFrame,
    max_genes: int,
    force_include: list[str] | None = None,
    biotype_map: pd.DataFrame | None = None,
    allowed_biotypes: list[str] | None = None,
    exclude_noise_symbols: bool = True,
) -> list[str]:
    """Select top genes by variance, with biotype filtering and force-include.

    Pipeline:
        1. Drop ERCC spike-ins
        2. (Optional) Filter to allowed biotypes via id_map.tsv
        3. (Optional) Exclude Gm-prefix and other noise symbols
        4. Force-include genes get priority slots
        5. HVGs fill remaining slots so total ≤ max_genes

    Args:
        rtech: Expression matrix (genes x samples), index = Ensembl IDs
        max_genes: Maximum total genes in panel
        force_include: Ensembl IDs to force-include (bypass biotype filter)
        biotype_map: DataFrame with ensembl_gene_id, mgi_symbol, biotype
        allowed_biotypes: List of allowed biotypes (e.g. ["protein_coding"])
        exclude_noise_symbols: If True, drop Gm\\d+, Rik, etc. from HVG pool
    """
    keep = ~rtech.index.str.upper().str.startswith("ERCC")
    r = rtech.loc[keep]
    n_before = len(r)

    # ── Biotype filter ───────────────────────────────────────────────────
    if biotype_map is not None and allowed_biotypes:
        bt = biotype_map.set_index("ensembl_gene_id")
        # Genes in allowed biotypes
        bt_pass = set(bt[bt["biotype"].isin(allowed_biotypes)].index)
        # Force-include genes bypass the biotype filter
        force_set = set(force_include) if force_include else set()
        # Keep genes that pass biotype OR are force-included
        bio_mask = r.index.isin(bt_pass | force_set)
        n_removed_biotype = (~bio_mask).sum()
        r = r.loc[bio_mask]
        print(f"  Biotype filter ({', '.join(allowed_biotypes)}):")
        print(f"    Removed {n_removed_biotype} non-{'/'.join(allowed_biotypes)} genes")
        print(f"    Remaining: {len(r)} / {n_before}")
    else:
        bt = None

    # ── Noise-symbol filter (Gm\d+, Rik, etc.) ──────────────────────────
    if exclude_noise_symbols and biotype_map is not None:
        if bt is None:
            bt = biotype_map.set_index("ensembl_gene_id")
        # Build set of noise Ensembl IDs (only from HVG candidates, not forced)
        force_set = set(force_include) if force_include else set()
        noise_ids = set()
        for eid in r.index:
            if eid in force_set:
                continue  # never exclude forced genes
            sym = bt.loc[eid, "mgi_symbol"] if eid in bt.index else ""
            if _is_noise_symbol(sym):
                noise_ids.add(eid)
        if noise_ids:
            n_noise = len(noise_ids)
            r = r.loc[~r.index.isin(noise_ids)]
            print(f"  Noise-symbol filter (Gm\\d+, Rik, unmapped):")
            print(f"    Removed {n_noise} noise-symbol genes")
            print(f"    Remaining: {len(r)}")

    # ── Force-include ────────────────────────────────────────────────────
    forced: set[str] = set()
    if force_include:
        present = set(r.index)
        forced = {g for g in force_include if g in present}
        missing = len(force_include) - len(forced)
        if missing:
            print(f"  Note: {missing}/{len(force_include)} forced genes not in Rtech")
        if len(forced) > max_genes:
            raise ValueError(
                f"Force-include genes ({len(forced)}) exceed max_genes ({max_genes}). "
                f"Increase --max_genes or reduce the marker panel."
            )

    # Fill remaining slots with top-variance genes (excluding already-forced)
    remaining_budget = max(max_genes - len(forced), 0)
    v = r.drop(index=list(forced), errors="ignore").var(axis=1)
    hvg = set(v.sort_values(ascending=False).head(remaining_budget).index.tolist())

    genes = list(forced | hvg)
    print(f"  Gene panel: {len(forced)} forced + {len(hvg)} HVG = {len(genes)} total")
    return genes


def cell_standardize(
    rtech_gxs: pd.DataFrame,
    meta: pd.DataFrame,
    cell_cols: list[str],
    eps: float = 1e-8,
    sd_floor: float = 1e-3,
) -> np.ndarray:
    """
    Standardize within each experimental cell (defined by cell_cols).
    
    Args:
        rtech_gxs: genes x samples DataFrame
        meta: metadata with index = sample IDs
        cell_cols: columns defining experimental cells
        eps: small constant for numerical stability
        sd_floor: minimum SD to avoid division issues
    
    Returns:
        Z: (n_samples x n_genes) cell-standardized matrix
    """
    samples = rtech_gxs.columns.tolist()
    meta_aligned = meta.loc[samples]
    cell_key = meta_aligned[cell_cols].astype(str).agg("|".join, axis=1)
    
    # Work in samples x genes for vectorized operations
    X = rtech_gxs.T.values.astype(np.float64)  # (N x G)
    Z = np.empty_like(X)
    
    for ck in cell_key.unique():
        idx = np.where(cell_key.values == ck)[0]
        Xc = X[idx, :]  # (n_cell x G)
        mu = Xc.mean(axis=0)
        # Fix: ddof=1 blows up for n=1 cells
        n_cell = Xc.shape[0]
        ddof = 1 if n_cell >= 2 else 0
        sd = Xc.std(axis=0, ddof=ddof)
        sd = np.where(np.isfinite(sd), sd, 0.0)  # NaN guard
        sd = np.maximum(sd, sd_floor)
        Z[idx, :] = (Xc - mu) / (sd + eps)
    
    return Z  # (N x G)


def partial_corr_from_precision(P: np.ndarray) -> np.ndarray:
    """Convert precision matrix to partial correlations."""
    d = np.sqrt(np.diag(P))
    pc = -P / np.outer(d, d)
    np.fill_diagonal(pc, 0.0)
    return pc


def topk_skeleton(pc: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Build skeleton by taking top-k neighbors per gene.
    
    The union of all top-k neighbors gives ~G*k edges (not G*k/2).
    """
    G = pc.shape[0]
    abs_pc = np.abs(pc)
    # Fix: guard against k > G-1
    k_eff = min(k, G - 1)
    
    edges = set()
    for i in range(G):
        # Fix: exclude self from consideration
        row = abs_pc[i].copy()
        row[i] = -np.inf
        # Use argpartition for O(n) instead of O(n log n)
        idx = np.argpartition(row, -k_eff)[-k_eff:]
        for j in idx:
            # Store as (min, max) for deduplication
            a, b = (i, j) if i < j else (j, i)
            edges.add((a, b))
    
    ii = np.fromiter((e[0] for e in edges), dtype=np.int32)
    jj = np.fromiter((e[1] for e in edges), dtype=np.int32)
    return ii, jj


def _read_gene_set_yaml(path: str | Path) -> dict:
    try:
        import yaml
    except Exception as exc:  # pragma: no cover
        raise ImportError("PyYAML is required for pathway prior edges") from exc
    path = Path(path)
    if not path.exists():
        return {}
    return yaml.safe_load(path.read_text()) or {}


def _flatten_gene_symbols(value) -> list[str]:
    genes: list[str] = []
    if isinstance(value, str):
        genes.append(value.strip())
    elif isinstance(value, list):
        for item in value:
            genes.extend(_flatten_gene_symbols(item))
    elif isinstance(value, dict):
        for item in value.values():
            genes.extend(_flatten_gene_symbols(item))
    return [g for g in genes if g]


def pathway_prior_edges(
    gene_sets_yaml: str | Path,
    selected_sets: list[str],
    id_map: str | Path,
    expression_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Resolve configured pathway genes and create pre-registered pairwise edges."""
    if not selected_sets:
        return pd.DataFrame(), pd.DataFrame(), []
    cfg = _read_gene_set_yaml(gene_sets_yaml)
    rows: list[dict[str, object]] = []
    reports: list[pd.DataFrame] = []
    force_ids: set[str] = set()
    for set_name in selected_sets:
        if set_name not in cfg or not isinstance(cfg[set_name], dict):
            reports.append(pd.DataFrame([{
                "source": "preregistered_pathway",
                "set": set_name,
                "query": "",
                "ensembl_gene_id": "",
                "status": "missing_gene_set",
                "in_panel": False,
            }]))
            continue
        symbols = _flatten_gene_symbols(cfg[set_name].get("genes", []))
        resolved = resolve_configured_genes(symbols, id_map, panel_genes=expression_genes)
        resolved["source"] = "preregistered_pathway"
        resolved["set"] = set_name
        reports.append(resolved)
        present = sorted(set(resolved.loc[
            (resolved["status"] == "mapped") &
            (resolved["ensembl_gene_id"] != "") &
            (resolved["in_panel"] == True),
            "ensembl_gene_id",
        ]))
        force_ids.update(present)
        for a, b in combinations(present, 2):
            rows.append({
                "gene_i": min(a, b),
                "gene_j": max(a, b),
                "edge_origin": "preregistered_pathway",
                "edge_source_detail": set_name,
                "is_fixed_prior": True,
            })
    report = pd.concat(reports, ignore_index=True) if reports else pd.DataFrame()
    return pd.DataFrame(rows).drop_duplicates(), report, sorted(force_ids)


def external_prior_edges(
    paths: list[str],
    id_map: str | Path,
    expression_genes: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Load external/prior edge tables and resolve gene symbols/IDs to Ensembl IDs."""
    rows: list[dict[str, object]] = []
    reports: list[dict[str, object]] = []
    force_ids: set[str] = set()
    for raw_path in paths:
        path = Path(raw_path)
        if not path.exists():
            reports.append({
                "source": "external_prior",
                "edge_source_detail": str(path),
                "query_i": "",
                "query_j": "",
                "gene_i": "",
                "gene_j": "",
                "status": "missing_file",
            })
            continue
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        table = pd.read_csv(path, sep=sep)
        lower = {c.lower(): c for c in table.columns}
        gi_col = lower.get("gene_i") or lower.get("source") or lower.get("gene1")
        gj_col = lower.get("gene_j") or lower.get("target") or lower.get("gene2")
        if gi_col is None or gj_col is None:
            raise ValueError(f"{path} must contain gene_i/gene_j columns")
        detail_col = lower.get("edge_source_detail") or lower.get("source_detail") or lower.get("evidence")
        for _, row in table.iterrows():
            qi = str(row[gi_col]).strip()
            qj = str(row[gj_col]).strip()
            ri = resolve_configured_genes([qi], id_map, panel_genes=expression_genes)
            rj = resolve_configured_genes([qj], id_map, panel_genes=expression_genes)
            ids_i = sorted(set(ri.loc[
                (ri["status"] == "mapped") & (ri["ensembl_gene_id"] != "") & (ri["in_panel"] == True),
                "ensembl_gene_id",
            ]))
            ids_j = sorted(set(rj.loc[
                (rj["status"] == "mapped") & (rj["ensembl_gene_id"] != "") & (rj["in_panel"] == True),
                "ensembl_gene_id",
            ]))
            detail = str(row[detail_col]).strip() if detail_col else path.name
            if not ids_i or not ids_j:
                reports.append({
                    "source": "external_prior",
                    "edge_source_detail": detail,
                    "query_i": qi,
                    "query_j": qj,
                    "gene_i": ",".join(ids_i),
                    "gene_j": ",".join(ids_j),
                    "status": "unresolved_or_absent_from_expression",
                })
                continue
            for a in ids_i:
                for b in ids_j:
                    if a == b:
                        continue
                    force_ids.update([a, b])
                    rows.append({
                        "gene_i": min(a, b),
                        "gene_j": max(a, b),
                        "edge_origin": "external_prior",
                        "edge_source_detail": detail,
                        "is_fixed_prior": True,
                    })
                    reports.append({
                        "source": "external_prior",
                        "edge_source_detail": detail,
                        "query_i": qi,
                        "query_j": qj,
                        "gene_i": min(a, b),
                        "gene_j": max(a, b),
                        "status": "mapped",
                    })
    edge_df = pd.DataFrame(rows).drop_duplicates() if rows else pd.DataFrame()
    report_df = pd.DataFrame(reports)
    return edge_df, report_df, sorted(force_ids)


def union_skeleton_with_priors(
    genes: list[str],
    pc: np.ndarray,
    ii: np.ndarray,
    jj: np.ndarray,
    prior_edges: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Return a provenance-aware union skeleton and aligned edge index arrays."""
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    records: dict[tuple[int, int], dict[str, object]] = {}

    for a, b in zip(ii.tolist(), jj.tolist()):
        key = (min(a, b), max(a, b))
        records[key] = {
            "gene_i": genes[key[0]],
            "gene_j": genes[key[1]],
            "pcorr": float(pc[key[0], key[1]]),
            "abs_pcorr": float(abs(pc[key[0], key[1]])),
            "i": key[0],
            "j": key[1],
            "edge_origin": "osd771_data",
            "edge_source_detail": "topk_partial_correlation",
            "is_fixed_prior": False,
        }

    if prior_edges is not None and not prior_edges.empty:
        for _, row in prior_edges.iterrows():
            gi = str(row["gene_i"])
            gj = str(row["gene_j"])
            if gi not in gene_to_idx or gj not in gene_to_idx or gi == gj:
                continue
            a, b = sorted([gene_to_idx[gi], gene_to_idx[gj]])
            key = (a, b)
            origin = str(row.get("edge_origin", "external_prior"))
            detail = str(row.get("edge_source_detail", ""))
            if key in records:
                existing_origin = set(str(records[key]["edge_origin"]).split("|"))
                existing_origin.add(origin)
                existing_detail = set(str(records[key]["edge_source_detail"]).split("|"))
                if detail:
                    existing_detail.add(detail)
                records[key]["edge_origin"] = "|".join(sorted(existing_origin))
                records[key]["edge_source_detail"] = "|".join(sorted(existing_detail))
                records[key]["is_fixed_prior"] = True
            else:
                records[key] = {
                    "gene_i": genes[a],
                    "gene_j": genes[b],
                    "pcorr": float(pc[a, b]),
                    "abs_pcorr": float(abs(pc[a, b])),
                    "i": a,
                    "j": b,
                    "edge_origin": origin,
                    "edge_source_detail": detail,
                    "is_fixed_prior": True,
                }

    edge_df = pd.DataFrame(records.values()).sort_values(["i", "j"]).reset_index(drop=True)
    out_i = edge_df["i"].to_numpy(np.int32)
    out_j = edge_df["j"].to_numpy(np.int32)
    return edge_df, out_i, out_j


def main():
    ap = argparse.ArgumentParser(description="Build Phase 2 shared skeleton E")
    ap.add_argument("--rtech", default="data/processed/phase1_residuals/Rtech.tsv.gz",
                    help="Path to Rtech.tsv.gz (genes x samples)")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz",
                    help="Path to meta_phase1.tsv.gz")
    ap.add_argument("--outdir", default="data/processed/networks/phase2",
                    help="Output directory")
    ap.add_argument("--max_genes", type=int, default=2500,
                    help="Maximum genes for skeleton")
    ap.add_argument("--cell_cols", default="Age,Arm,EnvGroup",
                    help="Comma-separated columns defining experimental cells")
    ap.add_argument("--topk", type=int, default=DEFAULT_TOPK,
                    help="Top-k neighbors per gene (~G*k edges)")
    ap.add_argument("--force_include", default="",
                    help="Path to gene list file (one gene per line) to force-include")
    # ── Biotype filter arguments ─────────────────────────────────────────
    ap.add_argument("--id_map", default="data/processed/resources/id_map.tsv",
                    help="Path to id_map.tsv with biotype annotations "
                         "(from build_id_map.py)")
    ap.add_argument("--biotype_filter", default="protein_coding",
                    help="Comma-separated allowed biotypes. "
                         "Set to 'none' to disable. Default: protein_coding")
    ap.add_argument("--no_noise_symbol_filter", action="store_true",
                    help="Disable filtering of Gm-prefix / Rik / unmapped symbols")
    ap.add_argument("--anchor_config", default=str(REPO_ROOT / "config/anchor_genes.yaml"),
                    help="Configured Procrustes anchor YAML. Mapped anchors are force-included.")
    ap.add_argument("--no_anchor_force_include", action="store_true",
                    help="Disable configured-anchor force inclusion (not recommended; for diagnostics only)")
    ap.add_argument("--anchor_report", default="anchor_force_include_report.tsv",
                    help="Filename under outdir for configured-anchor mapping report")
    ap.add_argument("--edge_priors", default="",
                    help="Comma-separated external/prior kidney coexpression edge tables with gene_i/gene_j columns")
    ap.add_argument("--gene_sets", default=str(REPO_ROOT / "config/gene_sets.yaml"),
                    help="Configured gene-set YAML used to create pre-registered pathway prior edges")
    ap.add_argument("--pathway_prior_sets", default="dct_ncc_wnk,ion_transport,calcium_handling",
                    help="Comma-separated gene-set keys to turn into pre-registered pairwise prior edges. "
                         "Set to 'none' to disable.")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 2: Build Shared Sparse Skeleton E")
    print("=" * 60)

    # Load data
    print(f"\nLoading Rtech: {args.rtech}")
    rtech = load_rtech(args.rtech)
    print(f"  Shape: {rtech.shape[0]} genes × {rtech.shape[1]} samples")

    print(f"Loading metadata: {args.meta}")
    meta = load_meta(args.meta)

    # Find sample column
    sample_col = None
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            sample_col = col
            break
    if sample_col is None:
        sample_col = meta.columns[0]
        print(f"  Warning: using first column '{sample_col}' as sample ID")
    meta = meta.set_index(sample_col, drop=False)

    # Align
    common = [s for s in rtech.columns if s in meta.index]
    rtech = rtech[common]
    meta = meta.loc[common]
    print(f"  Aligned: {len(common)} samples")

    # Gene selection
    force_include = []
    if args.force_include and Path(args.force_include).exists():
        force_include.extend([
            line.strip() for line in
            Path(args.force_include).read_text().strip().split("\n")
            if line.strip()
        ])
        print(f"\nLoaded {len(force_include)} force-include genes from {args.force_include}")

    # Load biotype annotations if filtering is enabled
    biotype_map = None
    allowed_biotypes = None
    if args.biotype_filter.lower() != "none":
        id_map_path = Path(args.id_map)
        if not id_map_path.is_absolute():
            id_map_path = REPO_ROOT / id_map_path
        if id_map_path.exists():
            biotype_map = load_biotype_map(str(id_map_path))
            allowed_biotypes = [b.strip() for b in args.biotype_filter.split(",")]
            print(f"\nBiotype filter enabled: {allowed_biotypes}")
            print(f"  Loaded {len(biotype_map)} annotations from {id_map_path}")
        else:
            print(f"\n  WARNING: id_map not found at {id_map_path}, "
                  f"skipping biotype filter. Run build_id_map.py first.")

    prior_edge_tables: list[pd.DataFrame] = []
    prior_reports: list[pd.DataFrame] = []
    prior_force_ids: list[str] = []
    id_map_for_priors = Path(args.id_map)
    if not id_map_for_priors.is_absolute():
        id_map_for_priors = REPO_ROOT / id_map_for_priors
    expression_genes = set(rtech.index)

    if id_map_for_priors.exists():
        pathway_sets = [
            s.strip() for s in args.pathway_prior_sets.split(",")
            if s.strip() and s.strip().lower() != "none"
        ]
        if pathway_sets:
            p_edges, p_report, p_force = pathway_prior_edges(
                gene_sets_yaml=args.gene_sets,
                selected_sets=pathway_sets,
                id_map=id_map_for_priors,
                expression_genes=expression_genes,
            )
            if not p_edges.empty:
                prior_edge_tables.append(p_edges)
            if not p_report.empty:
                prior_reports.append(p_report)
            prior_force_ids.extend(p_force)
            print(f"\nPre-registered pathway prior edges: {len(p_edges)} from sets {pathway_sets}")

        prior_paths = [p.strip() for p in args.edge_priors.split(",") if p.strip()]
        if prior_paths:
            e_edges, e_report, e_force = external_prior_edges(
                paths=prior_paths,
                id_map=id_map_for_priors,
                expression_genes=expression_genes,
            )
            if not e_edges.empty:
                prior_edge_tables.append(e_edges)
            if not e_report.empty:
                prior_reports.append(e_report)
            prior_force_ids.extend(e_force)
            print(f"External/prior kidney edges: {len(e_edges)} from {len(prior_paths)} table(s)")
    else:
        print(f"\n  WARNING: id_map not found at {id_map_for_priors}; skipping prior-edge resolution")

    if prior_force_ids:
        force_include.extend(prior_force_ids)
        print(f"  Prior-edge force-include genes present in expression: {len(set(prior_force_ids))}")

    if not args.no_anchor_force_include:
        anchor_config = Path(args.anchor_config)
        if not anchor_config.is_absolute():
            anchor_config = REPO_ROOT / anchor_config
        id_map_for_anchors = Path(args.id_map)
        if not id_map_for_anchors.is_absolute():
            id_map_for_anchors = REPO_ROOT / id_map_for_anchors
        print(f"\nConfigured-anchor force include: {anchor_config}")
        anchor_ids, anchor_report = configured_anchor_ids(
            anchor_config=anchor_config,
            id_map=id_map_for_anchors,
            expression_genes=set(rtech.index),
        )
        force_include.extend(anchor_ids)
        report_path = outdir / args.anchor_report
        anchor_report.to_csv(report_path, sep="\t", index=False)
        print(f"  Mapped anchors present in expression universe: {len(anchor_ids)}")
        print(f"  Anchor mapping report: {report_path}")

    force_include = sorted(set(force_include))
    if force_include:
        print(f"  Total force-include genes after DCT/anchor merge: {len(force_include)}")
    else:
        force_include = None

    genes = pick_genes(
        rtech,
        args.max_genes,
        force_include=force_include,
        biotype_map=biotype_map,
        allowed_biotypes=allowed_biotypes,
        exclude_noise_symbols=not args.no_noise_symbol_filter,
    )
    (outdir / "phase2_genes.txt").write_text("\n".join(genes) + "\n")
    print(f"\nSelected {len(genes)} genes for skeleton")
    print(f"  → Saved to {outdir / 'phase2_genes.txt'}")

    rtech_gxs = rtech.loc[genes]
    cell_cols = [c.strip() for c in args.cell_cols.split(",") if c.strip()]
    print(f"\nCell columns: {cell_cols}")

    # Step 2.1: Cell-standardize for topology selection
    print("\nStep 2.1: Cell-standardizing within each experimental cell...")
    n_cells = meta[cell_cols].astype(str).agg("|".join, axis=1).nunique()
    print(f"  Found {n_cells} experimental cells (n≈{len(common)//n_cells} per cell)")
    Z = cell_standardize(rtech_gxs, meta, cell_cols=cell_cols)
    print(f"  Cell-standardized matrix: {Z.shape[0]} samples × {Z.shape[1]} genes")

    # Step 2.3: Shrinkage partial correlation
    print("\nStep 2.3: Computing shrinkage partial correlations (Ledoit-Wolf)...")
    lw = LedoitWolf().fit(Z)
    cov = lw.covariance_
    print(f"  Shrinkage: {lw.shrinkage_:.4f}")
    
    # Invert covariance to get precision (with ridge fallback)
    try:
        prec = np.linalg.inv(cov)
    except np.linalg.LinAlgError:
        ridge = 1e-6 * np.eye(cov.shape[0])
        prec = np.linalg.inv(cov + ridge)
        print("  Warning: added tiny ridge for inversion")
    pc = partial_corr_from_precision(prec)
    
    # Build data-driven skeleton
    print(f"\nBuilding skeleton with top-k={args.topk} neighbors per gene...")
    ii_data, jj_data = topk_skeleton(pc, k=args.topk)

    prior_edges = pd.concat(prior_edge_tables, ignore_index=True).drop_duplicates() if prior_edge_tables else pd.DataFrame()
    if prior_reports:
        prior_report = pd.concat(prior_reports, ignore_index=True)
        prior_report.to_csv(outdir / "edge_prior_resolution_report.tsv", sep="\t", index=False)
        print(f"  Edge prior resolution report: {outdir / 'edge_prior_resolution_report.tsv'}")

    edge_df, ii, jj = union_skeleton_with_priors(
        genes=genes,
        pc=pc,
        ii=ii_data,
        jj=jj_data,
        prior_edges=prior_edges,
    )

    # Save edge indices for downstream determinism
    np.save(outdir / "edge_i.npy", ii)
    np.save(outdir / "edge_j.npy", jj)

    edge_df.to_csv(outdir / "skeleton_edges.tsv", sep="\t", index=False)
    provenance_summary = (
        edge_df.groupby(["edge_origin", "is_fixed_prior"], dropna=False)
        .size()
        .reset_index(name="n_edges")
        .sort_values("n_edges", ascending=False)
    )
    provenance_summary.to_csv(outdir / "skeleton_edge_provenance_summary.tsv", sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print("Skeleton E built successfully")
    print(f"{'=' * 60}")
    print(f"  Genes: {len(genes)}")
    print(f"  Edges: {len(edge_df)} (data-driven target: ~{len(genes)*args.topk//2} to ~{len(genes)*args.topk})")
    if not prior_edges.empty:
        print(f"  Prior edges requested: {len(prior_edges)}; fixed-prior edges in skeleton: {int(edge_df['is_fixed_prior'].sum())}")
    print(f"\nOutputs in {outdir}:")
    print(f"  - phase2_genes.txt")
    print(f"  - skeleton_edges.tsv (with pcorr weights + provenance)")
    print(f"  - edge_i.npy, edge_j.npy (indices)")
    print(f"  - skeleton_edge_provenance_summary.tsv")


if __name__ == "__main__":
    main()
