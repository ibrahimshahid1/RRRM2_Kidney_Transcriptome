#!/usr/bin/env python3
"""
Build Ensembl → Symbol ID Mapping
==================================

Queries the Ensembl REST API to map mouse Ensembl gene IDs to MGI symbols.
Generates a cached TSV for use by Phase 7 biological grounding.

Uses two strategies:
  1) POST /lookup/id (batch, fast) — maps Ensembl → display_name
  2) GET /xrefs/symbol/mus_musculus/<symbol> (fallback for unmapped symbols)

Usage:
    python scripts/build_id_map.py \\
        --genes data/processed/networks/run_xxx/phase2_genes.txt \\
        --outdir data/processed/resources
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from datetime import datetime

import requests
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
ENSEMBL_REST = "https://rest.ensembl.org"


def ensembl_batch_lookup(ids: list[str], chunk: int = 500,
                         sleep: float = 0.3, retries: int = 5) -> dict:
    """
    Batch POST to /lookup/id.
    Returns dict: {ensembl_id: {display_name, biotype, description, ...}}
    """
    out = {}
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    total_chunks = (len(ids) + chunk - 1) // chunk
    for c_idx in range(0, len(ids), chunk):
        batch = ids[c_idx:c_idx + chunk]
        chunk_num = c_idx // chunk + 1
        print(f"  Querying chunk {chunk_num}/{total_chunks} "
              f"({len(batch)} IDs)...", end="", flush=True)

        for attempt in range(retries):
            try:
                r = requests.post(
                    f"{ENSEMBL_REST}/lookup/id",
                    headers=headers,
                    data=json.dumps({"ids": batch}),
                    timeout=60
                )
                if r.status_code == 429:
                    wait = float(r.headers.get("Retry-After", 2.0 + attempt))
                    print(f" rate-limited, waiting {wait}s...", end="")
                    time.sleep(wait)
                    continue
                r.raise_for_status()
                data = r.json()
                out.update(data)
                print(f" OK ({len(data)} mapped)")
                break
            except (requests.RequestException, json.JSONDecodeError) as e:
                print(f" error: {e}, retry {attempt + 1}/{retries}")
                time.sleep(2.0 + attempt)
        else:
            print(f" FAILED after {retries} retries")

        time.sleep(sleep)

    return out


def ensembl_xref_symbol(symbol: str, species: str = "mus_musculus",
                        retries: int = 3) -> list[str]:
    """
    Fallback: look up a gene symbol → Ensembl IDs via xrefs.
    Returns list of Ensembl gene IDs matching the symbol.
    """
    headers = {"Accept": "application/json"}
    url = f"{ENSEMBL_REST}/xrefs/symbol/{species}/{symbol}"

    for attempt in range(retries):
        try:
            r = requests.get(url, headers=headers, timeout=30)
            if r.status_code == 429:
                time.sleep(float(r.headers.get("Retry-After", 2.0)))
                continue
            if r.status_code == 404:
                return []
            r.raise_for_status()
            data = r.json()
            return [entry["id"] for entry in data
                    if entry.get("type") == "gene"]
        except (requests.RequestException, json.JSONDecodeError):
            time.sleep(1.0 + attempt)

    return []


def main():
    ap = argparse.ArgumentParser(
        description="Build Ensembl → Symbol ID mapping via REST API"
    )
    ap.add_argument("--genes", required=True,
                    help="Path to gene list (one Ensembl ID per line)")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/processed/resources"),
                    help="Output directory for id_map.tsv")
    ap.add_argument("--extra_symbols", default="",
                    help="Comma-separated extra symbols to look up via xrefs "
                         "(e.g. Wnk1,Slc12a3,Kcnj10)")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Build Ensembl → Symbol ID Mapping")
    print("=" * 60)

    # Load gene IDs, strip Ensembl versions
    gene_ids = [
        line.strip().split(".")[0] for line in
        Path(args.genes).read_text().strip().split("\n")
        if line.strip()
    ]
    # Deduplicate while preserving order
    seen = set()
    gene_ids = [g for g in gene_ids if g not in seen and not seen.add(g)]
    print(f"\nLoaded {len(gene_ids)} Ensembl gene IDs from {args.genes}")

    # Query Ensembl REST release for reproducibility
    rest_release = "unknown"
    try:
        r = requests.get(f"{ENSEMBL_REST}/info/data/",
                         headers={"Accept": "application/json"}, timeout=10)
        if r.ok:
            releases = r.json().get("releases", [])
            if releases:
                rest_release = str(releases[0])
                print(f"Ensembl REST release: {rest_release}")
    except Exception:
        pass

    # Step 1: Batch lookup via REST
    print("\nStep 1: Batch POST /lookup/id ...")
    lookup = ensembl_batch_lookup(gene_ids)

    rows = []
    for gid in gene_ids:
        info = lookup.get(gid, {})
        if isinstance(info, dict) and "display_name" in info:
            rows.append({
                "ensembl_gene_id": gid,
                "mgi_symbol": info.get("display_name", ""),
                "biotype": info.get("biotype", ""),
                "description": info.get("description", ""),
            })
        else:
            rows.append({
                "ensembl_gene_id": gid,
                "mgi_symbol": "",
                "biotype": "",
                "description": "",
            })

    df = pd.DataFrame(rows)
    n_mapped = (df["mgi_symbol"] != "").sum()
    print(f"\n  Mapped: {n_mapped}/{len(gene_ids)} genes have symbols")

    # Step 2: Fallback xref for extra symbols (pathway genes that might
    # not be in the universe but we still want their Ensembl IDs)
    if args.extra_symbols:
        extra = [s.strip() for s in args.extra_symbols.split(",") if s.strip()]
        print(f"\nStep 2: xref lookup for {len(extra)} extra symbols...")
        xref_rows = []
        for sym in extra:
            ens_ids = ensembl_xref_symbol(sym)
            for eid in ens_ids:
                if eid not in set(df["ensembl_gene_id"]):
                    xref_rows.append({
                        "ensembl_gene_id": eid,
                        "mgi_symbol": sym,
                        "biotype": "",
                        "description": f"xref lookup for {sym}",
                    })
            time.sleep(0.2)  # rate limit
        if xref_rows:
            df = pd.concat([df, pd.DataFrame(xref_rows)], ignore_index=True)
            print(f"  Added {len(xref_rows)} extra mappings from xrefs")

    # Save
    map_path = outdir / "id_map.tsv"
    # Write header comment with timestamp for reproducibility
    with open(map_path, "w") as f:
        f.write(f"# Generated {datetime.now().isoformat()}\n")
        f.write(f"# Source: {ENSEMBL_REST}\n")
        f.write(f"# Ensembl release: {rest_release}\n")
        f.write(f"# Genes: {len(gene_ids)}\n")
    df.to_csv(map_path, sep="\t", index=False, mode="a")
    print(f"\nWrote: {map_path} ({len(df)} entries)")

    # Also save as JSON for programmatic use
    json_path = outdir / "id_map.json"
    mapping = {row["ensembl_gene_id"]: row["mgi_symbol"]
               for _, row in df.iterrows() if row["mgi_symbol"]}
    with open(json_path, "w") as f:
        json.dump(mapping, f, indent=2)
    print(f"Wrote: {json_path}")

    print(f"\n[OK] ID mapping complete → {outdir}")


if __name__ == "__main__":
    main()
