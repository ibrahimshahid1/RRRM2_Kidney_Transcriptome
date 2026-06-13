# src/enrichment/biological_grounding.py
"""
Phase 7: Fast biological grounding (poster-friendly).

1) Test enrichment of pre-registered gene sets among high-Δ genes (top decile)
   using Fisher's exact test.

Optional:
2) Cluster a reference embedding (k-means) and test cluster enrichment of high-Δ genes.

Inputs:
  - rewiring agg table (has gene + rewiring_mean)
  - optional embedding npy for module clustering
  - optional gene mapping TSV (Ensembl -> Symbol) for readability

Outputs:
  data/results/phase7_grounding/
    gene_set_enrichment.tsv
    cluster_enrichment.tsv (if embedding provided)
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from src.common import REPO_ROOT, bh_fdr

from src.enrichment.gene_set_loader import load_gene_sets


def fisher_exact(a, b, c, d):
    """Fisher's exact test with scipy fallback."""
    try:
        from scipy.stats import fisher_exact as fe
        OR, p = fe([[a, b], [c, d]], alternative="greater")
        return OR, p
    except ImportError:
        # Fallback: compute odds ratio only
        OR = (a * d + 1e-9) / (b * c + 1e-9)
        return OR, np.nan





def main():
    ap = argparse.ArgumentParser(description="Fast biological grounding for rewiring")
    ap.add_argument("--rewiring", required=True,
                    help="Path to rewiring agg table (TSV with gene, rewiring_mean)")
    ap.add_argument("--outdir", 
                    default=str(REPO_ROOT / "data/results/phase7_grounding"),
                    help="Output directory")
    ap.add_argument("--top_quantile", type=float, default=0.9,
                    help="Quantile threshold for 'high rewiring' (default: top decile = 0.9)")
    ap.add_argument("--map", default="", 
                    help="TSV with columns: ensembl_gene_id, mgi_symbol (REQUIRED when "
                         "universe uses Ensembl IDs and gene sets use symbols)")
    ap.add_argument("--embed", default="", 
                    help="Optional embedding .npy for clustering (genes x dim)")
    ap.add_argument("--k", type=int, default=12,
                    help="Number of clusters for k-means (if --embed provided)")
    ap.add_argument("--gene_de", default="",
                    help="Path to gene-level DE TSV (e.g. from Phase 6 regression). "
                         "If provided, expands the gene universe beyond network genes "
                         "to include all expressed genes for enrichment testing.")
    ap.add_argument("--libraries", default="",
                    help="Comma-separated Enrichr library names to fetch gene sets from. "
                         "Default: KEGG_2019_Mouse,WikiPathway_2023_Mouse,"
                         "Reactome_2022,MSigDB_Hallmark_2020")
    ap.add_argument("--gmt", default="",
                    help="Comma-separated paths to .gmt gene set files")
    ap.add_argument("--min_set_size", type=int, default=5,
                    help="Minimum gene set size to include (default: 5)")
    ap.add_argument("--max_set_size", type=int, default=500,
                    help="Maximum gene set size to include (default: 500)")
    ap.add_argument("--include_curated", action="store_true",
                    help="Also include legacy curated pre-registered gene sets")
    ap.add_argument("--curated_only", action="store_true",
                    help="Use ONLY curated sets (skip database fetch)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load rewiring table
    rw = pd.read_csv(args.rewiring, sep="\t")
    if "rewiring_mean" not in rw.columns:
        for c in ["delta_mean", "mean", "rewiring"]:
            if c in rw.columns:
                rw = rw.rename(columns={c: "rewiring_mean"})
                break
    
    if "gene" not in rw.columns or "rewiring_mean" not in rw.columns:
        raise ValueError(f"rewiring table must have 'gene' and 'rewiring_mean'. Found: {list(rw.columns)}")
    
    genes = rw["gene"].astype(str).tolist()
    print(f"Loaded {len(genes)} genes from {args.rewiring}")

    # Optionally expand universe with gene-level DE results (all expressed genes)
    extra_universe_genes = set()
    if args.gene_de and Path(args.gene_de).exists():
        de_df = pd.read_csv(args.gene_de, sep="\t")
        # Find gene column
        gene_col = None
        for c in ["gene", "gene_id", "ensembl_gene_id", "Gene"]:
            if c in de_df.columns:
                gene_col = c
                break
        if gene_col is None and de_df.index.dtype == object:
            extra_universe_genes = set(de_df.index.astype(str))
        elif gene_col:
            extra_universe_genes = set(de_df[gene_col].astype(str))
        if extra_universe_genes:
            n_new = len(extra_universe_genes - set(genes))
            print(f"Expanded universe with {n_new} additional genes from DE results "
                  f"(total universe: {len(set(genes) | extra_universe_genes)})")

    # Detect ID type
    frac_ens = sum(1 for g in genes if g.startswith("ENSMUSG")) / max(len(genes), 1)
    is_ensembl = frac_ens > 0.5
    print(f"ID type: {'Ensembl' if is_ensembl else 'Symbol'} ({frac_ens:.0%} ENSMUSG)")

    # Load mapping file
    sym_to_ens: dict[str, set[str]] = {}  # symbol (lowercase) → set of Ensembl IDs
    ens_to_sym: dict[str, str] = {}        # Ensembl ID → symbol

    if args.map and Path(args.map).exists():
        m = pd.read_csv(args.map, sep="\t", comment="#")
        # Normalise column names
        col_map = {}
        for c in m.columns:
            cl = c.lower().strip()
            if "ensembl" in cl:
                col_map[c] = "ensembl_gene_id"
            elif "symbol" in cl or "mgi" in cl:
                col_map[c] = "mgi_symbol"
            elif "in_universe" in cl:
                col_map[c] = "in_universe"
        m = m.rename(columns=col_map)

        # Filter to universe rows only (extras are for reporting, not enrichment)
        if "in_universe" in m.columns:
            n_before = len(m)
            m = m[m["in_universe"].astype(str).str.lower().isin(["true", "1"])]
            n_filtered = n_before - len(m)
            if n_filtered > 0:
                print(f"  Filtered out {n_filtered} non-universe rows")

        if "ensembl_gene_id" in m.columns and "mgi_symbol" in m.columns:
            for _, row in m.iterrows():
                eid = str(row["ensembl_gene_id"]).strip()
                sym = str(row["mgi_symbol"]).strip()
                if eid and sym and sym != "nan" and sym != "":
                    ens_to_sym[eid] = sym
                    # Reverse map: symbol → set of Ensembl IDs (handles 1:many)
                    sym_to_ens.setdefault(sym.lower(), set()).add(eid)
            print(f"Loaded mapping: {len(ens_to_sym)} Ensembl → Symbol (universe only)")
            print(f"  Unique symbols: {len(sym_to_ens)}")
        else:
            print(f"WARNING: Map file columns not recognized: {list(m.columns)}")

        # Optionally load extras for annotation/reporting
        extras_path = Path(args.map).parent / "id_map_extras.tsv"
        if extras_path.exists():
            me = pd.read_csv(extras_path, sep="\t", comment="#")
            print(f"  Loaded {len(me)} extra symbol mappings for annotation")
    elif is_ensembl:
        print("WARNING: Universe uses Ensembl IDs but no --map provided. "
              "Gene set overlap will likely be 0.")

    # Define high-rewiring genes
    thr = rw["rewiring_mean"].quantile(args.top_quantile)
    high = set(rw.loc[rw["rewiring_mean"] >= thr, "gene"].astype(str))
    print(f"High rewiring threshold (q={args.top_quantile}): {thr:.4f} ({len(high)} genes)")

    # Load gene sets from database / curated
    libraries = None  # use defaults
    if args.curated_only:
        libraries = []  # skip Enrichr
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

    gene_universe = set(genes) | extra_universe_genes
    results = []
    total_unresolved = 0

    # Resolve gene sets to universe IDs and run Fisher's exact test
    for setname, symbols in loaded_sets.items():
        set_genes = set()
        unresolved = []

        for s in symbols:
            s_lower = s.lower()
            if sym_to_ens and is_ensembl:
                ens_ids = sym_to_ens.get(s_lower, set())
                matched = ens_ids & gene_universe
                if matched:
                    set_genes |= matched
                else:
                    unresolved.append(s)
            else:
                for g in gene_universe:
                    if g.lower() == s_lower:
                        set_genes.add(g)
                        break
                else:
                    unresolved.append(s)

        total_unresolved += len(unresolved)

        A = len(high & set_genes)
        B = len(high - set_genes)
        C = len((gene_universe - high) & set_genes)
        D = len((gene_universe - high) - set_genes)
        OR, p = fisher_exact(A, B, C, D)

        matched_syms = sorted(ens_to_sym.get(g, g) for g in (high & set_genes))

        results.append({
            "gene_set": setname,
            "library": set_to_library.get(setname, "unknown"),
            "set_size_in_universe": len(set_genes),
            "hits_in_top_decile": A,
            "odds_ratio": round(OR, 4),
            "p": p,
            "hit_symbols": "; ".join(matched_syms) if matched_syms else "",
        })

    if total_unresolved:
        print(f"\nNote: {total_unresolved} total symbol→ID lookups unresolved "
              f"across {len(loaded_sets)} gene sets")

    res = pd.DataFrame(results)
    res["q_BH"] = bh_fdr(res["p"].fillna(1.0).values)
    res = res.sort_values("p").reset_index(drop=True)

    # Warn if ALL gene sets have 0 overlap (ID mismatch)
    total_overlap = res["set_size_in_universe"].sum()
    if total_overlap == 0:
        msg = ("ERROR: All gene sets have 0 overlap with the gene universe. "
               "This almost certainly means an ID type mismatch.\n")
        if is_ensembl and not sym_to_ens:
            msg += ("Universe uses Ensembl IDs but no --map file was loaded. "
                    "Run: python scripts/build_id_map.py --genes <gene_list> "
                    "--outdir data/processed/resources")
        print(msg)
        raise RuntimeError(msg)

    # Save full results
    enrich_path = outdir / "gene_set_enrichment.tsv"
    res.to_csv(enrich_path, sep="\t", index=False)
    print(f"\nWrote {enrich_path} ({len(res)} gene sets tested)")

    # Save significant subset
    sig = res[res["q_BH"] < 0.05]
    if len(sig) > 0:
        sig_path = outdir / "gene_set_enrichment_significant.tsv"
        sig.to_csv(sig_path, sep="\t", index=False)
        print(f"Wrote {sig_path} ({len(sig)} significant at FDR < 0.05)")

    # Print top results (up to 30)
    show = res.head(30)
    print(f"\nTop results ({len(show)} of {len(res)}):")
    for _, row in show.iterrows():
        marker = "*" if row["q_BH"] < 0.05 else " "
        print(f"  {marker} {row['gene_set']:55s}  "
              f"overlap={row['set_size_in_universe']:3d}  "
              f"hits={row['hits_in_top_decile']}  "
              f"OR={row['odds_ratio']:.2f}  "
              f"q={row['q_BH']:.4f}")

    # Optional: embedding module clustering
    if args.embed and Path(args.embed).exists():
        X = np.load(args.embed)
        if X.shape[0] != len(genes):
            raise ValueError(f"Embedding rows ({X.shape[0]}) must match #genes ({len(genes)})")
        
        try:
            from sklearn.cluster import KMeans
            km = KMeans(n_clusters=args.k, n_init=10, random_state=0)
            cl = km.fit_predict(X)
            print(f"Clustered into {args.k} clusters")
        except ImportError:
            print("[WARN] sklearn not available, skipping cluster enrichment")
            cl = None

        if cl is not None:
            dfc = pd.DataFrame({"gene": genes, "cluster": cl})
            dfc["is_high"] = dfc["gene"].isin(high)

            rows = []
            for c in sorted(dfc["cluster"].unique()):
                sub = dfc[dfc["cluster"] == c]
                A = int(sub["is_high"].sum())
                B = int(dfc["is_high"].sum() - A)
                C = int(len(sub) - A)
                D = int(len(dfc) - len(sub) - B)
                OR, p = fisher_exact(A, B, C, D)
                rows.append([c, len(sub), A, OR, p])

            ce = pd.DataFrame(rows, columns=["cluster", "cluster_size", "hits_in_top_decile", "odds_ratio", "p"])
            ce["q_BH"] = bh_fdr(ce["p"].fillna(1.0).values)
            cluster_path = outdir / "cluster_enrichment.tsv"
            ce.to_csv(cluster_path, sep="\t", index=False)
            print(f"Wrote {cluster_path}")

    print(f"\n[OK] Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
