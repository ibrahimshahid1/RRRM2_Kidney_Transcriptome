# src/enrichment/gene_set_loader.py
"""
Gene set loader for enrichment analysis.

Fetches gene set collections from Enrichr (via gseapy), caches locally
for reproducibility, and supports offline .gmt files.

Symbol resolution uses the Ensembl-backed id_map.json (built by
src.data.build_id_map via the Ensembl REST API) instead of the
old title-casing heuristic.

Usage:
    from src.enrichment.gene_set_loader import load_gene_sets

    sets, libs = load_gene_sets()
    # sets : dict[set_name → list[mouse_symbol]]
    # libs : dict[set_name → library_name]
"""
from __future__ import annotations

import json
from pathlib import Path

import yaml

from src.common import REPO_ROOT

DEFAULT_CACHE = REPO_ROOT / "data" / "processed" / "gene_sets"
DEFAULT_CURATED_YAML = REPO_ROOT / "config" / "gene_sets.yaml"
DEFAULT_ID_MAP = REPO_ROOT / "data" / "processed" / "resources" / "id_map.json"
DEFAULT_SEGMENT_MARKERS = REPO_ROOT / "data" / "processed" / "segment_markers"

# Default Enrichr libraries (mouse-specific first)
DEFAULT_LIBRARIES: list[str] = [
    "KEGG_2019_Mouse",
    "WikiPathways_2024_Mouse",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]


# Ensembl-backed symbol resolution

def load_id_map(path: Path | None = None) -> dict[str, str]:
    """Load the Ensembl→Symbol JSON map and build a case-insensitive
    symbol → canonical symbol lookup.

    Returns dict: lowercase_symbol → canonical_mouse_symbol
    """
    if path is None:
        path = DEFAULT_ID_MAP

    lookup: dict[str, str] = {}

    if not path.exists():
        print(f"  WARNING: id_map not found at {path}")
        print("  Run `python -m src.data.build_id_map` first.")
        print("  Symbol resolution will fall back to raw strings.")
        return lookup

    with open(path) as f:
        raw = json.load(f)  # {ensembl_id: symbol}

    # Build case-insensitive symbol → canonical symbol
    for _ens_id, sym in raw.items():
        if sym:
            lookup[sym.lower()] = sym
            # Also map the Ensembl ID itself
            lookup[_ens_id.lower()] = sym

    return lookup


def resolve_symbol(sym: str, lookup: dict[str, str]) -> str:
    """Resolve a gene symbol to its canonical mouse form using the
    Ensembl-backed id_map.

    Falls back to the original string if not found.
    """
    canonical = lookup.get(sym.lower())
    if canonical:
        return canonical
    return sym


# Curated sets from YAML

def load_curated_from_yaml(
    yaml_path: Path | None = None,
    id_map_path: Path | None = None,
) -> dict[str, list[str]]:
    """Load curated gene sets from config/gene_sets.yaml and resolve
    each symbol against the Ensembl-backed id_map.

    Returns: dict[set_name → list[resolved_mouse_symbol]]
    """
    if yaml_path is None:
        yaml_path = DEFAULT_CURATED_YAML
    if id_map_path is None:
        id_map_path = DEFAULT_ID_MAP

    if not yaml_path.exists():
        print(f"  WARNING: Curated gene sets YAML not found: {yaml_path}")
        return {}

    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    lookup = load_id_map(id_map_path)

    sets: dict[str, list[str]] = {}
    skipped = []

    for key, val in config.items():
        # Skip non-gene-set entries (enrichment_databases, etc.)
        if not isinstance(val, dict) or "genes" not in val:
            continue

        genes_raw = val["genes"]
        if not isinstance(genes_raw, list):
            continue

        # Flatten: could be a flat list of strings, or a dict of sublists
        flat: list[str] = []
        for item in genes_raw:
            if isinstance(item, str):
                flat.append(item.strip())
            elif isinstance(item, list):
                flat.extend(s.strip() for s in item if isinstance(s, str))
            elif isinstance(item, dict):
                # nested subcategory like {kinases: [...]}
                for sublist in item.values():
                    if isinstance(sublist, list):
                        flat.extend(s.strip() for s in sublist if isinstance(s, str))

        # Resolve symbols
        resolved = []
        for sym in flat:
            canonical = resolve_symbol(sym, lookup)
            resolved.append(canonical)

        if resolved:
            sets[key] = resolved
        else:
            skipped.append(key)

    if skipped:
        print(f"  Skipped empty gene sets: {skipped}")

    n_resolved = sum(1 for s in sets.values() for g in s
                     if g.lower() in lookup)
    n_total = sum(len(s) for s in sets.values())
    if lookup:
        print(f"  Curated sets: {len(sets)} loaded, "
              f"{n_resolved}/{n_total} symbols resolved via Ensembl id_map")
    else:
        print(f"  Curated sets: {len(sets)} loaded (no id_map — using raw symbols)")

    return sets


# Data-driven segment marker panels

def load_segment_marker_panels(
    marker_dir: Path | None = None,
) -> dict[str, list[str]]:
    """Load data-driven marker panels from Phase 1.5b (discover_markers.py).

    Each segment has a <SEGMENT>_marker_panel.txt file containing one gene
    per line.  These are loaded as 'segment_markers::SEGMENT' gene sets.

    Also checks the run-specific results directory.
    """
    import os

    sets: dict[str, list[str]] = {}
    dirs_to_check: list[Path] = []

    if marker_dir is not None:
        dirs_to_check.append(marker_dir)

    # Check current run's segment_markers dir
    results_dir = os.environ.get("RRRM_RESULTS_DIR")
    if results_dir:
        dirs_to_check.append(Path(results_dir) / "segment_markers")

    # Check default location
    dirs_to_check.append(DEFAULT_SEGMENT_MARKERS)

    found_dir = None
    for d in dirs_to_check:
        if d.exists() and list(d.glob("*_marker_panel.txt")):
            found_dir = d
            break

    if found_dir is None:
        return sets

    for panel_file in sorted(found_dir.glob("*_marker_panel.txt")):
        segment = panel_file.stem.replace("_marker_panel", "")
        genes = [g.strip() for g in panel_file.read_text().strip().split("\n")
                 if g.strip()]
        if genes:
            sets[segment] = genes

    return sets


# Enrichr fetcher

def fetch_enrichr_library(name: str, cache_dir: Path) -> dict[str, list[str]]:
    """Fetch an Enrichr gene-set library; cache as JSON for offline reuse."""
    cache_file = cache_dir / f"{name}.json"

    if cache_file.exists():
        print(f"  Loading from cache: {cache_file.name}")
        with open(cache_file) as f:
            return json.load(f)

    try:
        import gseapy as gp
        lib: dict[str, list[str]] = gp.get_library(name=name)
        cache_dir.mkdir(parents=True, exist_ok=True)
        with open(cache_file, "w") as f:
            json.dump(lib, f, indent=1)
        print(f"  Fetched and cached: {len(lib)} gene sets")
        return lib
    except ImportError:
        print("  WARNING: gseapy not installed — run: pip install gseapy")
        return {}
    except Exception as e:
        print(f"  WARNING: Could not fetch '{name}': {e}")
        return {}


# GMT file loader

def load_gmt(path: Path) -> dict[str, list[str]]:
    """Parse a standard .gmt file (MSigDB / GSEA format)."""
    gene_sets: dict[str, list[str]] = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g.strip() for g in parts[2:] if g.strip()]
            if genes:
                gene_sets[name] = genes
    return gene_sets


# Main entry point

def load_gene_sets(
    libraries: list[str] | None = None,
    gmt_files: list[str | Path] | None = None,
    cache_dir: Path | None = None,
    id_map_path: Path | None = None,
    curated_yaml: Path | None = None,
    min_size: int = 5,
    max_size: int = 500,
    include_curated: bool = True,
    include_segment_markers: bool = True,
    segment_marker_dir: Path | None = None,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """
    Load gene sets from Enrichr libraries, .gmt files, and/or curated YAML.

    All symbols are resolved to canonical mouse form using the Ensembl-backed
    id_map.json (if available).  Falls back to raw strings if unavailable.

    Parameters
    ----------
    libraries : list of Enrichr library names (None → DEFAULT_LIBRARIES)
    gmt_files : list of paths to .gmt files
    cache_dir : directory for caching Enrichr downloads
    id_map_path : path to id_map.json (None → default location)
    curated_yaml : path to curated gene sets YAML (None → config/gene_sets.yaml)
    min_size, max_size : gene-set size filter (inclusive)
    include_curated : also include curated sets from YAML

    Returns
    -------
    gene_sets : dict[str, list[str]]
        set_name → list of mouse gene symbols
    set_to_library : dict[str, str]
        set_name → source library name
    """
    if cache_dir is None:
        cache_dir = DEFAULT_CACHE
    cache_dir.mkdir(parents=True, exist_ok=True)

    if libraries is None:
        libraries = list(DEFAULT_LIBRARIES)

    # Load symbol resolver
    lookup = load_id_map(id_map_path)

    gene_sets: dict[str, list[str]] = {}
    set_to_library: dict[str, str] = {}

    # 1) Enrichr libraries
    for lib_name in libraries:
        print(f"Gene set library: {lib_name}")
        raw = fetch_enrichr_library(lib_name, cache_dir)
        if not raw:
            continue

        # Resolve symbols via Ensembl id_map (no more title-casing hack)
        normed: dict[str, list[str]] = {}
        for k, v in raw.items():
            normed[k] = [resolve_symbol(s, lookup) for s in v]

        n_kept = 0
        for set_name, genes in normed.items():
            if min_size <= len(genes) <= max_size:
                key = f"{lib_name}::{set_name}"
                gene_sets[key] = genes
                set_to_library[key] = lib_name
                n_kept += 1

        print(f"  Kept {n_kept} / {len(normed)} sets (size {min_size}–{max_size})")

    # 2) .gmt files
    if gmt_files:
        for gmt_path in gmt_files:
            p = Path(gmt_path)
            if not p.exists():
                print(f"  WARNING: GMT not found: {p}")
                continue
            raw = load_gmt(p)
            lib_name = p.stem
            n_kept = 0
            for set_name, genes in raw.items():
                genes = [resolve_symbol(s, lookup) for s in genes]
                if min_size <= len(genes) <= max_size:
                    key = f"{lib_name}::{set_name}"
                    gene_sets[key] = genes
                    set_to_library[key] = lib_name
                    n_kept += 1
            print(f"Gene set library: {lib_name} — {n_kept} / {len(raw)} sets from .gmt")

    # 3) Curated sets from YAML
    if include_curated:
        curated = load_curated_from_yaml(curated_yaml, id_map_path)
        n_cur = 0
        for set_name, genes in curated.items():
            key = f"curated::{set_name}"
            if min_size <= len(genes) <= max_size:
                gene_sets[key] = genes
                set_to_library[key] = "curated"
                n_cur += 1
        print(f"Gene set library: curated — {n_cur} sets from YAML")

    # 4) Data-driven segment marker panels
    if include_segment_markers:
        panels = load_segment_marker_panels(segment_marker_dir)
        n_seg = 0
        for seg_name, genes in panels.items():
            key = f"segment_markers::{seg_name}"
            if min_size <= len(genes) <= max_size:
                gene_sets[key] = genes
                set_to_library[key] = "segment_markers"
                n_seg += 1
        if panels:
            print(f"Gene set library: segment_markers — {n_seg} data-driven panels")

    print(f"\nTotal gene sets loaded: {len(gene_sets)}")
    return gene_sets, set_to_library


# CLI: list available Enrichr libraries
if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="List or fetch Enrichr gene set libraries")
    ap.add_argument("--list", action="store_true", help="List available Enrichr libraries")
    ap.add_argument("--fetch", nargs="*", default=None,
                    help="Fetch and cache these libraries (default: DEFAULT_LIBRARIES)")
    ap.add_argument("--cache", default=str(DEFAULT_CACHE))
    args = ap.parse_args()

    if args.list:
        try:
            import gseapy as gp
            names = gp.get_library_name()
            mouse = [n for n in names if "mouse" in n.lower()]
            print(f"--- Mouse-specific ({len(mouse)}) ---")
            for n in sorted(mouse):
                print(f"  {n}")
            print(f"\n--- All ({len(names)}) ---")
            for n in sorted(names):
                print(f"  {n}")
        except ImportError:
            print("gseapy not installed — run: pip install gseapy")

    if args.fetch is not None:
        libs = args.fetch if args.fetch else DEFAULT_LIBRARIES
        load_gene_sets(libraries=libs, cache_dir=Path(args.cache))
