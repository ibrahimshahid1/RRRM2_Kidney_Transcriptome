# src/enrichment/gene_set_loader.py
"""
Gene set loader for enrichment analysis.

Fetches gene set collections from Enrichr (via gseapy), caches locally
for reproducibility, and supports offline .gmt files. All gene symbols
are normalised to mouse convention (title case).

Usage:
    from src.enrichment.gene_set_loader import load_gene_sets

    sets, libs = load_gene_sets()
    # sets : dict[set_name → list[mouse_symbol]]
    # libs : dict[set_name → library_name]
"""
from __future__ import annotations

import json
from pathlib import Path

from src.common import REPO_ROOT
DEFAULT_CACHE = REPO_ROOT / "data" / "processed" / "gene_sets"

# ── Default Enrichr libraries (mouse-specific first) ──────────────────────
DEFAULT_LIBRARIES: list[str] = [
    "KEGG_2019_Mouse",
    "WikiPathways_2024_Mouse",
    "Reactome_2022",
    "MSigDB_Hallmark_2020",
]

# ── Curated pre-registered sets (expanded, best-effort accurate) ──────────
CURATED_SETS: dict[str, list[str]] = {
    "DCT_NCC_WNK_axis": [
        "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Slc12a3",
        "Kcnj10", "Kcnj16", "Scnn1a", "Scnn1b", "Scnn1g",
        "Nedd4l", "Sgk1", "Klhl3", "Cul3",
    ],
    "Oxidative_stress": [
        "Nfe2l2", "Sod1", "Sod2", "Cat", "Gpx1", "Gpx4",
        "Prdx1", "Prdx2", "Hmox1", "Keap1",
    ],
    "ECM_remodeling": [
        "Col1a1", "Col1a2", "Col3a1", "Col4a1", "Col4a2",
        "Fn1", "Mmp2", "Mmp9", "Mmp14", "Timp1", "Timp2",
        "Lox", "Sparc", "Tgfb1", "Lama2", "Lamb1",
    ],
    "Lipid_metabolism": [
        "Ppara", "Pparg", "Acox1", "Cpt1a", "Fabp1",
        "Hmgcr", "Srebf1", "Srebf2", "Fasn", "Acaca", "Scd1",
    ],
    "Calcium_handling": [
        "Trpv5", "Trpv6", "Calb1", "S100g", "Atp2b1", "Atp2b4",
        "Ryr1", "Ryr2", "Itpr1", "Itpr2", "Camk2a", "Camk2d",
    ],
    "Inflammation": [
        "Il6", "Tnf", "Ccl2", "Il1b", "Nfkb1", "Nfkb2",
        "Rela", "Tlr4", "Cxcl1", "Il10",
    ],
    "Apoptosis": [
        "Bax", "Bcl2", "Casp3", "Casp9", "Trp53", "Fas", "Cycs",
    ],
    "Ion_transport": [
        "Slc12a3", "Slc12a1", "Kcnj1", "Kcnma1", "Atp1a1",
        "Atp1b1", "Slc4a1", "Slc9a3", "Slc22a6", "Kcnj10", "Kcnj16",
    ],
    "Fibrosis": [
        "Tgfb1", "Tgfb2", "Acta2", "Col1a1", "Col3a1",
        "Vim", "Fn1", "Ctgf", "Lox",
    ],
}


# ── Symbol normalisation ──────────────────────────────────────────────────
def _to_mouse_symbol(sym: str) -> str:
    """Best-effort human (UPPER) → mouse (Title) symbol conversion."""
    if not sym or len(sym) < 2:
        return sym
    # Already mouse-style (first upper, rest has lowercase)
    if sym[0].isupper() and any(c.islower() for c in sym[1:]):
        return sym
    # ALL-CAPS → Title case:  COL1A1 → Col1a1,  SLC12A3 → Slc12a3
    return sym[0].upper() + sym[1:].lower()


# ── Enrichr fetcher ───────────────────────────────────────────────────────
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


# ── GMT file loader ───────────────────────────────────────────────────────
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


# ── Main entry point ──────────────────────────────────────────────────────
def load_gene_sets(
    libraries: list[str] | None = None,
    gmt_files: list[str | Path] | None = None,
    cache_dir: Path | None = None,
    min_size: int = 5,
    max_size: int = 500,
    include_curated: bool = True,
) -> tuple[dict[str, list[str]], dict[str, str]]:
    """
    Load gene sets from Enrichr libraries, .gmt files, and/or curated lists.

    Parameters
    ----------
    libraries : list of Enrichr library names (None → DEFAULT_LIBRARIES)
    gmt_files : list of paths to .gmt files
    cache_dir : directory for caching Enrichr downloads
    min_size, max_size : gene-set size filter (inclusive)
    include_curated : also include the CURATED_SETS dict

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

    gene_sets: dict[str, list[str]] = {}
    set_to_library: dict[str, str] = {}

    # ── 1) Enrichr libraries ──────────────────────────────────────────
    for lib_name in libraries:
        print(f"Gene set library: {lib_name}")
        raw = fetch_enrichr_library(lib_name, cache_dir)
        if not raw:
            continue

        # Always normalise to mouse-style symbols (Title case)
        # Even "Mouse" libraries on Enrichr use UPPERCASE human symbols
        normed: dict[str, list[str]] = {}
        for k, v in raw.items():
            normed[k] = [_to_mouse_symbol(s) for s in v]

        n_kept = 0
        for set_name, genes in normed.items():
            if min_size <= len(genes) <= max_size:
                key = f"{lib_name}::{set_name}"
                gene_sets[key] = genes
                set_to_library[key] = lib_name
                n_kept += 1

        print(f"  Kept {n_kept} / {len(normed)} sets (size {min_size}–{max_size})")

    # ── 2) .gmt files ────────────────────────────────────────────────
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
                genes = [_to_mouse_symbol(s) for s in genes]
                if min_size <= len(genes) <= max_size:
                    key = f"{lib_name}::{set_name}"
                    gene_sets[key] = genes
                    set_to_library[key] = lib_name
                    n_kept += 1
            print(f"Gene set library: {lib_name} — {n_kept} / {len(raw)} sets from .gmt")

    # ── 3) Curated pre-registered sets ───────────────────────────────
    if include_curated:
        n_cur = 0
        for set_name, genes in CURATED_SETS.items():
            key = f"curated::{set_name}"
            if min_size <= len(genes) <= max_size:
                gene_sets[key] = genes
                set_to_library[key] = "curated"
                n_cur += 1
        print(f"Gene set library: curated — {n_cur} pre-registered sets")

    print(f"\nTotal gene sets loaded: {len(gene_sets)}")
    return gene_sets, set_to_library


# ── CLI: list available Enrichr libraries ─────────────────────────────────
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
