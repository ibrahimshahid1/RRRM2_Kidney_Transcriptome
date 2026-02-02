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
                    help="Optional TSV with columns: gene, symbol")
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

    # Optional gene symbol mapping
    sym = None
    if args.map and Path(args.map).exists():
        m = pd.read_csv(args.map, sep=None, engine="python")
        m = m.rename(columns={"Gene": "gene", "SYMBOL": "symbol", "Symbol": "symbol"})
        if "gene" in m.columns and "symbol" in m.columns:
            sym = dict(zip(m["gene"].astype(str), m["symbol"].astype(str)))
            print(f"Loaded gene->symbol mapping for {len(sym)} genes")

    # Define high-rewiring genes
    thr = rw["rewiring_mean"].quantile(args.top_quantile)
    high = set(rw.loc[rw["rewiring_mean"] >= thr, "gene"].astype(str))
    print(f"High rewiring threshold (q={args.top_quantile}): {thr:.4f} ({len(high)} genes)")

    # Pre-registered gene sets (SYMBOLS)
    # These are example kidney-relevant pathways - adjust as needed
    prereg = {
        "DCT_NCC_WNK_axis": [
            "WNK1", "WNK4", "STK39", "SLC12A3", "KCNJ10", "KCNJ16", 
            "SCNN1A", "SCNN1B", "SCNN1G"
        ],
        "Oxidative_stress": [
            "NFE2L2", "SOD1", "SOD2", "CAT", "GPX1", "PRDX1", "PRDX2", "HMOX1"
        ],
        "ECM_remodeling": [
            "COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2", "MMP9", "TIMP1", "SPARC"
        ],
        "Lipid_metabolism": [
            "PPARA", "PPARG", "SREBF1", "SREBF2", "ACACA", "FASN", "CPT1A"
        ],
        "Calcium_handling": [
            "ATP2B1", "ATP2B4", "ITPR1", "ITPR2", "RYR1", "RYR2", "CAMK2A", "CAMK2D"
        ],
        "Inflammation": [
            "IL6", "TNF", "IL1B", "CCL2", "NFKB1", "NFKB2", "RELA"
        ],
        "Apoptosis": [
            "BCL2", "BAX", "CASP3", "CASP9", "TP53", "CYCS"
        ],
    }

    gene_universe = set(genes)
    results = []

    # If we have symbol mapping, convert symbols to Ensembl IDs
    if sym is not None:
        # Build reverse map: symbol -> list of gene IDs
        rev = {}
        for g, s in sym.items():
            rev.setdefault(s.upper(), []).append(g)

        for setname, symbols in prereg.items():
            set_genes = set()
            for s in symbols:
                set_genes |= set(rev.get(s.upper(), []))
            set_genes &= gene_universe

            A = len(high & set_genes)         # high & in set
            B = len(high - set_genes)         # high & not in set
            C = len((gene_universe - high) & set_genes)
            D = len((gene_universe - high) - set_genes)
            OR, p = fisher_exact(A, B, C, D)
            results.append([setname, len(set_genes), A, OR, p])
    else:
        # Assume genes are already symbols
        for setname, symbols in prereg.items():
            set_genes = set(s.upper() for s in symbols) & set(g.upper() for g in gene_universe)
            # Map back to original case
            set_genes_orig = {g for g in gene_universe if g.upper() in set_genes}
            high_upper = {g.upper() for g in high}

            A = len({g for g in set_genes_orig if g.upper() in high_upper})
            B = len(high) - A
            C = len(set_genes_orig) - A
            D = len(gene_universe) - len(set_genes_orig) - B
            OR, p = fisher_exact(A, B, C, D)
            results.append([setname, len(set_genes_orig), A, OR, p])

    res = pd.DataFrame(results, columns=["gene_set", "set_size_in_universe", "hits_in_top_decile", "odds_ratio", "p"])
    res["q_BH"] = bh_fdr(res["p"].fillna(1.0).values)
    res = res.sort_values("p").reset_index(drop=True)
    
    enrich_path = outdir / "gene_set_enrichment.tsv"
    res.to_csv(enrich_path, sep="\t", index=False)
    print(f"Wrote {enrich_path}")

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
