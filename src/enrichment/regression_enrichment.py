# src/enrichment/regression_enrichment.py
"""
Phase 7 supplement: Run gene-set enrichment on regression-significant genes
(from Phase 6 full_regression) instead of rewiring-top-decile.

This tests: "Are genes with significant arm×flight or age×flight interaction
enriched for DCT/NCC-WNK or other predefined pathways?"

Usage:
    python -m src.enrichment.regression_enrichment \
        --reg_dir data/results/<run>/phase6_regression \
        --outdir data/results/<run>/phase7_regression_enrichment \
        --map data/processed/resources/id_map.tsv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

REPO_ROOT = Path(__file__).resolve().parents[2]

from src.enrichment.gene_set_loader import load_gene_sets


def load_id_map(map_path: Path) -> tuple[dict[str, str], dict[str, set[str]]]:
    """Load Ensembl→Symbol mapping. Returns (ens_to_sym, sym_to_ens)."""
    ens_to_sym = {}
    sym_to_ens: dict[str, set[str]] = {}
    if not map_path.exists():
        return ens_to_sym, sym_to_ens

    m = pd.read_csv(map_path, sep="\t", comment="#")
    col_map = {}
    for c in m.columns:
        cl = c.lower().strip()
        if "ensembl" in cl:
            col_map[c] = "ensembl_gene_id"
        elif "symbol" in cl or "mgi" in cl:
            col_map[c] = "mgi_symbol"
    m = m.rename(columns=col_map)

    if "ensembl_gene_id" in m.columns and "mgi_symbol" in m.columns:
        for _, row in m.iterrows():
            eid = str(row["ensembl_gene_id"]).strip()
            sym = str(row["mgi_symbol"]).strip()
            if eid and sym and sym != "nan":
                ens_to_sym[eid] = sym
                sym_to_ens.setdefault(sym.lower(), set()).add(eid)
    return ens_to_sym, sym_to_ens


def run_enrichment(sig_genes: set[str], universe: set[str],
                   gene_sets: dict[str, list[str]],
                   ens_to_sym: dict[str, str],
                   sym_to_ens: dict[str, set[str]],
                   is_ensembl: bool,
                   set_to_library: dict[str, str] | None = None) -> pd.DataFrame:
    """Fisher's exact test for each gene set."""
    results = []

    for gs_name, gs_symbols in gene_sets.items():
        # Map gene set symbols to universe IDs
        if is_ensembl:
            gs_ids = set()
            missing = []
            for sym in gs_symbols:
                ens_ids = sym_to_ens.get(sym.lower(), set())
                matched = ens_ids & universe
                if matched:
                    gs_ids.update(matched)
                else:
                    missing.append(sym)
            if missing:
                print(f"  {gs_name}: {len(missing)} symbols not in universe ({', '.join(missing[:3])})")
        else:
            gs_ids = {s for s in gs_symbols if s in universe}

        overlap = gs_ids & universe
        hits = gs_ids & sig_genes

        # 2x2 contingency table
        a = len(hits)
        b = len(sig_genes - gs_ids)
        c = len(overlap - sig_genes)
        d = len(universe - sig_genes - overlap)

        if len(overlap) == 0:
            results.append({
                "gene_set": gs_name,
                "library": set_to_library.get(gs_name, "") if set_to_library else "",
                "overlap": 0, "hits": 0,
                "odds_ratio": 0, "p_fisher": 1.0, "sig_genes_in_set": ""
            })
            continue

        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        # Get symbol names for hits
        if is_ensembl:
            hit_symbols = sorted([ens_to_sym.get(g, g) for g in hits])
        else:
            hit_symbols = sorted(hits)

        results.append({
            "gene_set": gs_name,
            "library": set_to_library.get(gs_name, "") if set_to_library else "",
            "overlap": len(overlap),
            "hits": a,
            "odds_ratio": round(odds, 2),
            "p_fisher": pval,
            "sig_genes_in_set": "; ".join(hit_symbols),
        })

    df = pd.DataFrame(results)
    if len(df) > 0:
        # BH correction
        df = df.sort_values("p_fisher").reset_index(drop=True)
        m = len(df)
        df["q_BH"] = np.minimum(1.0, df["p_fisher"] * m / (np.arange(1, m + 1)))
        df["q_BH"] = df["q_BH"][::-1].cummin()[::-1]
    return df


def main():
    ap = argparse.ArgumentParser(
        description="Gene-set enrichment on regression-significant genes"
    )
    ap.add_argument("--reg_dir", required=True,
                    help="Directory with Phase 6 regression outputs")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/results/phase7_regression_enrichment"))
    ap.add_argument("--map", default="",
                    help="TSV: ensembl_gene_id, mgi_symbol")
    ap.add_argument("--q_threshold", type=float, default=0.05,
                    help="FDR threshold for 'significant' genes (default: 0.05)")
    ap.add_argument("--libraries", default="",
                    help="Comma-separated Enrichr library names (default: KEGG_2019_Mouse,"
                         "WikiPathway_2023_Mouse,Reactome_2022,MSigDB_Hallmark_2020)")
    ap.add_argument("--gmt", default="",
                    help="Comma-separated paths to .gmt gene set files")
    ap.add_argument("--min_set_size", type=int, default=5,
                    help="Minimum gene set size (default: 5)")
    ap.add_argument("--max_set_size", type=int, default=500,
                    help="Maximum gene set size (default: 500)")
    ap.add_argument("--include_curated", action="store_true",
                    help="Also include legacy curated pre-registered gene sets")
    ap.add_argument("--curated_only", action="store_true",
                    help="Use ONLY curated sets (skip database fetch)")
    args = ap.parse_args()

    reg_dir = Path(args.reg_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load ID map
    ens_to_sym, sym_to_ens = {}, {}
    if args.map:
        ens_to_sym, sym_to_ens = load_id_map(Path(args.map))
        print(f"Loaded {len(ens_to_sym)} ID mappings")

    # Find all regression result files
    reg_files = sorted(reg_dir.glob("gene_*.tsv"))
    if not reg_files:
        print(f"No regression result files found in {reg_dir}")
        return

    print(f"Found {len(reg_files)} regression result files")

    # Load gene sets from database / curated
    libraries = None
    if args.curated_only:
        libraries = []
    elif args.libraries:
        libraries = [s.strip() for s in args.libraries.split(",") if s.strip()]

    gmt_files = [s.strip() for s in args.gmt.split(",") if s.strip()] if args.gmt else None
    include_curated = args.include_curated or args.curated_only

    loaded_sets, set_to_library = load_gene_sets(
        libraries=libraries,
        gmt_files=gmt_files,
        min_size=args.min_set_size,
        max_size=args.max_set_size,
        include_curated=include_curated,
    )

    for reg_file in reg_files:
        effect_name = reg_file.stem  # e.g. "gene_flight_effect"
        print(f"\n{'='*60}")
        print(f"Enrichment for: {effect_name}")
        print(f"{'='*60}")

        df = pd.read_csv(reg_file, sep="\t")
        if "gene" not in df.columns or "q_BH" not in df.columns:
            print(f"  Skipping — missing gene or q_BH column")
            continue

        universe = set(df["gene"].astype(str))
        sig_mask = df["q_BH"] < args.q_threshold
        sig_genes = set(df.loc[sig_mask, "gene"].astype(str))
        n_sig = len(sig_genes)

        print(f"  Universe: {len(universe)} genes")
        print(f"  Significant (q < {args.q_threshold}): {n_sig} genes")

        if n_sig == 0:
            print(f"  No significant genes — skipping enrichment")
            continue

        is_ensembl = sum(1 for g in universe if g.startswith("ENSMUSG")) / len(universe) > 0.5

        enrich = run_enrichment(sig_genes, universe, loaded_sets,
                                ens_to_sym, sym_to_ens, is_ensembl,
                                set_to_library=set_to_library)

        out_path = outdir / f"{effect_name}_enrichment.tsv"
        enrich.to_csv(out_path, sep="\t", index=False)

        # Save significant subset
        sig = enrich[enrich["q_BH"] < 0.05]
        if len(sig) > 0:
            sig_path = outdir / f"{effect_name}_enrichment_significant.tsv"
            sig.to_csv(sig_path, sep="\t", index=False)
            print(f"  {len(sig)} significant at FDR < 0.05")

        # Print top results (up to 20)
        show = enrich.head(20)
        for _, row in show.iterrows():
            flag = " *" if row["q_BH"] < 0.05 else "  "
            genes_str = f" [{row['sig_genes_in_set']}]" if row["hits"] > 0 else ""
            print(f"{flag} {row['gene_set']:55s} overlap={row['overlap']:3d}  "
                  f"hits={row['hits']}  OR={row['odds_ratio']:.2f}  "
                  f"q={row['q_BH']:.3f}{genes_str}")

        print(f"  Wrote: {out_path}")

    print(f"\n[OK] All outputs in: {outdir}")


if __name__ == "__main__":
    main()
