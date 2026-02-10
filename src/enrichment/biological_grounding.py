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

# Repository root (2 levels up from src/enrichment/)
REPO_ROOT = Path(__file__).resolve().parents[2]


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


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    p = np.asarray(p, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


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
        m = m.rename(columns=col_map)

        if "ensembl_gene_id" in m.columns and "mgi_symbol" in m.columns:
            for _, row in m.iterrows():
                eid = str(row["ensembl_gene_id"]).strip()
                sym = str(row["mgi_symbol"]).strip()
                if eid and sym and sym != "nan" and sym != "":
                    ens_to_sym[eid] = sym
                    # Reverse map: symbol → set of Ensembl IDs (handles 1:many)
                    sym_to_ens.setdefault(sym.lower(), set()).add(eid)
            print(f"Loaded mapping: {len(ens_to_sym)} Ensembl → Symbol")
            print(f"  Unique symbols: {len(sym_to_ens)}")
        else:
            print(f"WARNING: Map file columns not recognized: {list(m.columns)}")
    elif is_ensembl:
        print("WARNING: Universe uses Ensembl IDs but no --map provided. "
              "Gene set overlap will likely be 0.")

    # Define high-rewiring genes
    thr = rw["rewiring_mean"].quantile(args.top_quantile)
    high = set(rw.loc[rw["rewiring_mean"] >= thr, "gene"].astype(str))
    print(f"High rewiring threshold (q={args.top_quantile}): {thr:.4f} ({len(high)} genes)")

    # Pre-registered gene sets — using MOUSE symbol style
    # (case-insensitive matching handles WNK1 vs Wnk1 vs wnk1)
    prereg = {
        "DCT_NCC_WNK_axis": [
            "Wnk1", "Wnk4", "Stk39", "Slc12a3", "Kcnj10", "Kcnj16", 
            "Scnn1a", "Scnn1b", "Scnn1g"
        ],
        "Oxidative_stress": [
            "Nfe2l2", "Sod1", "Sod2", "Cat", "Gpx1", "Prdx1", "Prdx2", "Hmox1"
        ],
        "ECM_remodeling": [
            "Col1a1", "Col1a2", "Col3a1", "Fn1", "Mmp2", "Mmp9", "Timp1", "Sparc"
        ],
        "Lipid_metabolism": [
            "Ppara", "Pparg", "Srebf1", "Srebf2", "Acaca", "Fasn", "Cpt1a"
        ],
        "Calcium_handling": [
            "Atp2b1", "Atp2b4", "Itpr1", "Itpr2", "Ryr1", "Ryr2", "Camk2a", "Camk2d"
        ],
        "Inflammation": [
            "Il6", "Tnf", "Il1b", "Ccl2", "Nfkb1", "Nfkb2", "Rela"
        ],
        "Apoptosis": [
            "Bcl2", "Bax", "Casp3", "Casp9", "Trp53", "Cycs"
        ],
    }

    gene_universe = set(genes)
    results = []

    # Resolve gene sets to universe IDs
    for setname, symbols in prereg.items():
        set_genes = set()
        unresolved = []

        for s in symbols:
            s_lower = s.lower()
            if sym_to_ens and is_ensembl:
                # Look up symbol → Ensembl IDs, intersect with universe
                ens_ids = sym_to_ens.get(s_lower, set())
                matched = ens_ids & gene_universe
                if matched:
                    set_genes |= matched
                else:
                    unresolved.append(s)
            else:
                # Direct symbol matching (case-insensitive)
                for g in gene_universe:
                    if g.lower() == s_lower:
                        set_genes.add(g)
                        break
                else:
                    unresolved.append(s)

        if unresolved:
            print(f"  {setname}: {len(unresolved)} symbols not in universe "
                  f"({', '.join(unresolved[:5])}{'...' if len(unresolved) > 5 else ''})")

        A = len(high & set_genes)         # high & in set
        B = len(high - set_genes)         # high & not in set
        C = len((gene_universe - high) & set_genes)
        D = len((gene_universe - high) - set_genes)
        OR, p = fisher_exact(A, B, C, D)

        # Include matched symbol names for readability
        matched_syms = []
        for g in (high & set_genes):
            matched_syms.append(ens_to_sym.get(g, g))

        results.append([setname, len(set_genes), A, OR, p,
                        "; ".join(sorted(matched_syms)) if matched_syms else ""])

    res = pd.DataFrame(results, columns=[
        "gene_set", "set_size_in_universe", "hits_in_top_decile",
        "odds_ratio", "p", "hit_symbols"
    ])
    res["q_BH"] = bh_fdr(res["p"].fillna(1.0).values)
    res = res.sort_values("p").reset_index(drop=True)
    
    # Hard-fail if ALL gene sets have 0 overlap (ID mismatch)
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

    enrich_path = outdir / "gene_set_enrichment.tsv"
    res.to_csv(enrich_path, sep="\t", index=False)
    print(f"Wrote {enrich_path}")
    print(f"\nResults:")
    for _, row in res.iterrows():
        marker = "*" if row["q_BH"] < 0.05 else " "
        print(f"  {marker} {row['gene_set']:25s}  "
              f"overlap={row['set_size_in_universe']:2d}  "
              f"hits={row['hits_in_top_decile']}  "
              f"OR={row['odds_ratio']:.2f}  "
              f"q={row['q_BH']:.3f}")

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
