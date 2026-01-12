from __future__ import annotations
from pathlib import Path
import re
import pandas as pd

IN_DIR = Path("tms_kidney_female")  # adjust if needed
OUT = IN_DIR / "gene_id_map_master.csv"

# Candidate columns that might contain symbols in cellxgene exports / AnnData var
SYMBOL_COL_CANDIDATES = [
    "gene_symbol", "gene_symbols", "symbol", "symbols",
    "feature_name", "gene_name", "name"
]

ENSM_RE = re.compile(r"^ENSMUSG\d+")

def looks_ensembl(s: str) -> bool:
    return bool(ENSM_RE.match(str(s)))

def choose_symbol_column(df: pd.DataFrame) -> str | None:
    cols = [c for c in df.columns if c in SYMBOL_COL_CANDIDATES]
    # Prefer a column that looks like gene symbols (not Ensembl)
    for c in cols:
        sample = df[c].dropna().astype(str).head(200)
        if len(sample) == 0:
            continue
        frac_ens = (sample.str.match(ENSM_RE)).mean()
        # if most values are NOT Ensembl, it's probably a symbol/name column
        if frac_ens < 0.2:
            return c
    # fallback: first candidate if exists
    return cols[0] if cols else None

def main():
    var_files = sorted(IN_DIR.glob("var_metadata_*.csv"))
    if not var_files:
        raise SystemExit(f"No var_metadata_*.csv found in {IN_DIR.resolve()}")

    rows = []
    for vf in var_files:
        df = pd.read_csv(vf, index_col=0)
        df.index = df.index.astype(str)

        # Ensembl IDs are usually the index in your case
        if not looks_ensembl(df.index[0]):
            # Sometimes Ensembl is in a column; try to find it
            ens_col = None
            for c in df.columns:
                sample = df[c].dropna().astype(str).head(200)
                if len(sample) and (sample.str.match(ENSM_RE).mean() > 0.8):
                    ens_col = c
                    break
            if ens_col is None:
                print(f"[WARN] Could not find Ensembl IDs in {vf.name}. Skipping.")
                continue
            ensembl = df[ens_col].astype(str)
        else:
            ensembl = pd.Series(df.index, index=df.index, name="ensembl_id")

        sym_col = choose_symbol_column(df)
        if sym_col is None:
            print(f"[WARN] No symbol-like column found in {vf.name}. Skipping symbols.")
            symbols = pd.Series([pd.NA] * len(df), index=df.index, dtype="string")
        else:
            symbols = df[sym_col].astype("string")

        tmp = pd.DataFrame({
            "ensembl_id": ensembl.values,
            "symbol": symbols.values,
            "source_file": vf.name,
        })

        # clean
        tmp["ensembl_id"] = tmp["ensembl_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
        tmp["symbol"] = tmp["symbol"].astype("string").str.strip()

        # keep only plausible Ensembl
        tmp = tmp[tmp["ensembl_id"].str.match(ENSM_RE, na=False)].copy()

        rows.append(tmp)

    m = pd.concat(rows, ignore_index=True)

    # Drop empty symbols
    m = m[~m["symbol"].isna() & (m["symbol"] != "")].copy()

    # Deduplicate: symbol -> most common Ensembl (or first)
    # (There are rare 1-to-many mappings; for markers, "best-effort" is fine.)
    sym_to_ens = (
        m.groupby(["symbol", "ensembl_id"])
         .size()
         .reset_index(name="n")
         .sort_values(["symbol", "n"], ascending=[True, False])
         .drop_duplicates("symbol")
         .drop(columns="n")
    )

    # Also store reverse mapping (Ensembl -> a symbol)
    ens_to_sym = (
        m.groupby(["ensembl_id", "symbol"])
         .size()
         .reset_index(name="n")
         .sort_values(["ensembl_id", "n"], ascending=[True, False])
         .drop_duplicates("ensembl_id")
         .drop(columns="n")
    )

    # Merge into one table for convenience
    merged = sym_to_ens.merge(ens_to_sym, on="ensembl_id", how="outer", suffixes=("_symbol_to_ens", "_ens_to_symbol"))
    merged.rename(columns={"symbol_symbol_to_ens": "symbol"}, inplace=True)

    # Final cleanup
    merged["symbol"] = merged["symbol"].astype("string")
    merged["symbol_ens_to_symbol"] = merged["symbol_ens_to_symbol"].astype("string")
    merged.to_csv(OUT, index=False)

    print(f"Saved master mapping: {OUT.resolve()}")
    print(f"Unique symbols mapped: {merged['symbol'].dropna().nunique()}")
    print(f"Unique Ensembl mapped: {merged['ensembl_id'].dropna().nunique()}")

if __name__ == "__main__":
    main()
