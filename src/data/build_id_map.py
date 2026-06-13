#!/usr/bin/env python3
"""Build Ensembl to Symbol ID Mapping"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from datetime import datetime

import yaml
import requests
import pandas as pd

from src.common import REPO_ROOT

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
    Look up a gene symbol → Ensembl gene IDs via xrefs.
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


def write_tsv_with_header(df: pd.DataFrame, path: Path,
                          header_lines: list[str]) -> None:
    """Write a TSV with # comment header lines for provenance."""
    with open(path, "w") as f:
        for line in header_lines:
            f.write(f"# {line}\n")
    df.to_csv(path, sep="\t", index=False, mode="a")


def main():
    ap = argparse.ArgumentParser(
        description="Build Ensembl → Symbol ID mapping via REST API"
    )
    ap.add_argument("--genes", required=True,
                    help="Path to gene list (one Ensembl ID per line) — "
                         "this defines the enrichment UNIVERSE")
    ap.add_argument("--outdir",
                    default=str(REPO_ROOT / "data/processed/resources"),
                    help="Output directory for id_map.tsv")
    ap.add_argument("--extra_symbols", default="",
                    help="Comma-separated extra symbols for annotation/reporting "
                         "(NOT added to enrichment universe). E.g. Wnk1,Slc12a3")
    ap.add_argument("--curated_yaml", default="",
                    help="Path to curated gene sets YAML. All symbols in the YAML "
                         "will be auto-resolved and added to the extras map. "
                         "(default: config/gene_sets.yaml if it exists)")
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
    seen: set[str] = set()
    gene_ids = [g for g in gene_ids if g not in seen and not seen.add(g)]
    gene_list_hash = hashlib.sha1("\n".join(gene_ids).encode()).hexdigest()[:12]
    print(f"\nLoaded {len(gene_ids)} Ensembl gene IDs from {args.genes}")
    print(f"  SHA1 (first 12): {gene_list_hash}")

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

    # Step 1: Batch lookup universe IDs
    print("\nStep 1: Batch POST /lookup/id (universe) ...")
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
                "source": "universe",
                "in_universe": True,
            })
        else:
            rows.append({
                "ensembl_gene_id": gid,
                "mgi_symbol": "",
                "biotype": "",
                "description": "",
                "source": "universe",
                "in_universe": True,
            })

    df_universe = pd.DataFrame(rows)
    n_mapped = (df_universe["mgi_symbol"] != "").sum()
    print(f"\n  Mapped: {n_mapped}/{len(gene_ids)} universe genes have symbols")

    # Step 2: Extra symbols for reporting (NOT in universe)
    # Collect extra symbols from --extra_symbols flag AND --curated_yaml
    extra_syms: list[str] = []
    if args.extra_symbols:
        extra_syms.extend(s.strip() for s in args.extra_symbols.split(",") if s.strip())

    # Auto-resolve curated YAML symbols
    curated_yaml_path = Path(args.curated_yaml) if args.curated_yaml else (REPO_ROOT / "config" / "gene_sets.yaml")
    if curated_yaml_path.exists():
        print(f"\nLoading curated gene sets from: {curated_yaml_path}")
        with open(curated_yaml_path) as f:
            curated_config = yaml.safe_load(f)
        for key, val in curated_config.items():
            if isinstance(val, dict) and "genes" in val:
                genes_raw = val["genes"]
                if isinstance(genes_raw, list):
                    for item in genes_raw:
                        if isinstance(item, str):
                            extra_syms.append(item.strip())
                        elif isinstance(item, list):
                            extra_syms.extend(s.strip() for s in item if isinstance(s, str))
                        elif isinstance(item, dict):
                            for sublist in item.values():
                                if isinstance(sublist, list):
                                    extra_syms.extend(s.strip() for s in sublist if isinstance(s, str))
        # Deduplicate
        extra_syms = list(dict.fromkeys(extra_syms))
        print(f"  Collected {len(extra_syms)} unique symbols from YAML + CLI")

    df_extras = pd.DataFrame()
    if extra_syms:
        print(f"\nStep 2: xref lookup for {len(extra_syms)} extra symbols...")

        present = set(df_universe["ensembl_gene_id"])
        xref_ens_ids = []
        xref_requested: dict[str, str] = {}  # ens_id → requested symbol

        for sym in extra_syms:
            ens_ids = ensembl_xref_symbol(sym)
            for eid in ens_ids:
                if eid not in present:
                    xref_ens_ids.append(eid)
                    xref_requested[eid] = sym
                    present.add(eid)
                else:
                    print(f"    {sym} → {eid} (already in universe)")
            time.sleep(0.2)  # rate limit

        # Re-lookup xref IDs to get canonical display_name
        if xref_ens_ids:
            print(f"  Re-looking up {len(xref_ens_ids)} xref IDs for canonical names...")
            xref_lookup = ensembl_batch_lookup(xref_ens_ids)

            xref_rows = []
            for eid in xref_ens_ids:
                info = xref_lookup.get(eid, {})
                canonical = info.get("display_name", "") if isinstance(info, dict) else ""
                requested = xref_requested[eid]
                xref_rows.append({
                    "ensembl_gene_id": eid,
                    "mgi_symbol": canonical or requested,  # prefer canonical
                    "requested_symbol": requested,
                    "biotype": info.get("biotype", "") if isinstance(info, dict) else "",
                    "description": info.get("description", "") if isinstance(info, dict) else "",
                    "source": "extra_symbol",
                    "in_universe": False,
                })

            df_extras = pd.DataFrame(xref_rows)
            print(f"  Added {len(df_extras)} extra mappings (for reporting only)")

    # Dedup universe: prefer non-empty symbol, drop true duplicates
    df_universe = df_universe.sort_values(
        ["ensembl_gene_id", "mgi_symbol"],
        key=lambda s: s.apply(lambda x: (0 if x else 1) if isinstance(x, str) else x)
    )
    df_universe = df_universe.drop_duplicates(subset="ensembl_gene_id", keep="first")

    # Provenance header
    header = [
        f"Generated {datetime.now().isoformat()}",
        f"Source: {ENSEMBL_REST}",
        f"Ensembl release: {rest_release}",
        f"Universe genes: {len(df_universe)}",
        f"Gene list: {args.genes}",
        f"Gene list SHA1: {gene_list_hash}",
    ]

    # Save universe map
    map_path = outdir / "id_map.tsv"
    write_tsv_with_header(df_universe, map_path, header)
    print(f"\nWrote: {map_path} ({len(df_universe)} universe entries)")

    # Save extras separately
    if len(df_extras) > 0:
        extras_path = outdir / "id_map_extras.tsv"
        write_tsv_with_header(df_extras, extras_path, header)
        print(f"Wrote: {extras_path} ({len(df_extras)} extra entries)")

    # JSON for programmatic use (universe only)
    json_path = outdir / "id_map.json"
    mapping = {row["ensembl_gene_id"]: row["mgi_symbol"]
               for _, row in df_universe.iterrows() if row["mgi_symbol"]}
    with open(json_path, "w") as f:
        json.dump(mapping, f, indent=2)
    print(f"Wrote: {json_path}")

    print(f"\n[OK] ID mapping complete → {outdir}")
    print(f"     Universe (enrichment denominator): id_map.tsv")
    if len(df_extras) > 0:
        print(f"     Extra symbols (reporting only):    id_map_extras.tsv")


if __name__ == "__main__":
    main()
