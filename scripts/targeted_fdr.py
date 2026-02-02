# scripts/targeted_fdr.py
import argparse
import pandas as pd
from statsmodels.stats.multitest import multipletests

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm", required=True, help="phase6 *_perm_pvals.tsv")
    ap.add_argument("--targets", required=True, help="txt file with one Ensembl gene ID per line")
    ap.add_argument("--out", required=True, help="output TSV")
    args = ap.parse_args()

    perm = pd.read_csv(args.perm, sep="\t")
    targets = [x.strip() for x in open(args.targets) if x.strip() and not x.startswith("#")]
    targ = perm[perm["gene"].isin(targets)].copy()

    if targ.empty:
        raise SystemExit("No target genes matched the perm file. Check ID format (Ensembl versions?).")

    # BH only within target list
    targ["q_BH_targets"] = multipletests(targ["p_perm"].values, method="fdr_bh")[1]

    # Useful ranks
    perm_sorted = perm.sort_values("rewiring_abs_obs", ascending=False).reset_index(drop=True)
    perm_sorted["rank_all"] = perm_sorted.index + 1
    targ = targ.merge(perm_sorted[["gene", "rank_all"]], on="gene", how="left")

    targ = targ.sort_values(["q_BH_targets", "p_perm", "rank_all"])
    targ.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] Wrote {args.out} ({len(targ)} targets)")

if __name__ == "__main__":
    main()
