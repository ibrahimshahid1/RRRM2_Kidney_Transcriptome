"""Frozen gene-panel registry with offline resolution and API-backed verification.

Panels are read from config/panels/*.yaml at analysis time — never from the
network — so a run stays deterministic and its provenance hash stays verifiable.
The Ensembl REST calls here exist only to *verify* that stored symbols still
resolve and to *refresh* a panel into a new dated snapshot under explicit
operator control. They are never on the analysis path.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator

import yaml

from .paths import PANELS_DIR, REPO_ROOT, relative

ENSEMBL_REST = "https://rest.ensembl.org"
SPECIES = {"mouse": "mus_musculus", "human": "homo_sapiens"}


class PanelError(ValueError):
    """A panel file is malformed, or a requested panel does not exist."""


@dataclass(frozen=True)
class Panel:
    """One named gene panel with direction weights and provenance."""

    id: str
    description: str
    species: str
    genes: dict[str, int]
    minimum_present: int | None
    provenance: dict[str, Any]
    source: Path

    @property
    def symbols(self) -> list[str]:
        return list(self.genes)

    @property
    def frozen(self) -> bool:
        return bool(self.provenance.get("frozen", False))

    def digest(self) -> str:
        """Stable SHA-256 over membership and direction, independent of file layout."""
        payload = json.dumps(
            {"id": self.id, "species": self.species,
             "genes": dict(sorted(self.genes.items()))},
            sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()


def iter_panel_files(directory: Path | None = None) -> Iterator[Path]:
    """Yield panel files, skipping underscore-prefixed helpers."""
    directory = directory or PANELS_DIR
    if not directory.exists():
        return
    for path in sorted(directory.glob("*.yaml")):
        if not path.name.startswith("_"):
            yield path


def load_panels(directory: Path | None = None) -> dict[str, Panel]:
    """Load every panel from the registry, keyed by panel id."""
    out: dict[str, Panel] = {}
    for path in iter_panel_files(directory):
        raw = yaml.safe_load(path.read_text()) or {}
        for pid, body in (raw.get("panels") or {}).items():
            if pid in out:
                raise PanelError(
                    f"duplicate panel id {pid!r} in {relative(path)} "
                    f"(already defined in {relative(out[pid].source)})"
                )
            genes = body.get("genes") or {}
            if isinstance(genes, list):
                genes = {g: 1 for g in genes}
            if not genes:
                raise PanelError(f"{relative(path)}: panel {pid!r} has no genes")
            out[pid] = Panel(
                id=pid,
                description=(body.get("description") or "").strip(),
                species=body.get("species") or "mouse",
                genes={str(k): int(v) for k, v in genes.items()},
                minimum_present=body.get("minimum_present"),
                provenance=body.get("provenance") or {},
                source=path,
            )
    return out


def get_panel(panel_id: str, directory: Path | None = None) -> Panel:
    """Resolve one panel by id. This is the function analysis code should call."""
    panels = load_panels(directory)
    if panel_id not in panels:
        raise PanelError(
            f"unknown panel {panel_id!r}; known panels: {', '.join(sorted(panels))}"
        )
    return panels[panel_id]


def resolve(panel_id: str, *, directory: Path | None = None,
            signed: bool = False) -> list[str] | dict[str, int]:
    """Return a panel's symbols, or its symbol->direction map when ``signed``."""
    panel = get_panel(panel_id, directory)
    return dict(panel.genes) if signed else panel.symbols


# --- API-backed verification and refresh (never on the analysis path) --------

def _lookup_symbols(symbols: list[str], species: str, *,
                    timeout: int = 60, chunk: int = 200) -> dict[str, Any]:
    """POST /lookup/symbol for a species; returns the raw Ensembl payload."""
    import requests  # imported lazily so offline analysis never needs it

    sp = SPECIES.get(species, species)
    out: dict[str, Any] = {}
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    for i in range(0, len(symbols), chunk):
        batch = symbols[i:i + chunk]
        resp = requests.post(
            f"{ENSEMBL_REST}/lookup/symbol/{sp}",
            headers=headers, data=json.dumps({"symbols": batch}), timeout=timeout)
        resp.raise_for_status()
        out.update(resp.json())
    return out


def verify_panel(panel: Panel, *, timeout: int = 60) -> dict[str, Any]:
    """Check every symbol against Ensembl; report unresolved and renamed symbols."""
    payload = _lookup_symbols(panel.symbols, panel.species, timeout=timeout)
    resolved, unresolved, renamed = {}, [], []
    for symbol in panel.symbols:
        hit = payload.get(symbol)
        if not hit:
            unresolved.append(symbol)
            continue
        display = hit.get("display_name") or symbol
        resolved[symbol] = {
            "ensembl_id": hit.get("id"),
            "display_name": display,
            "biotype": hit.get("biotype"),
        }
        if display != symbol:
            renamed.append({"stored": symbol, "current": display})
    return {
        "panel": panel.id,
        "species": panel.species,
        "frozen": panel.frozen,
        "digest": panel.digest(),
        "n_genes": len(panel.symbols),
        "n_resolved": len(resolved),
        "unresolved": unresolved,
        "renamed": renamed,
        "resolved": resolved,
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "source": ENSEMBL_REST,
    }


def snapshot_path(panel: Panel, directory: Path | None = None) -> Path:
    """Dated snapshot path for a refreshed panel; refresh never overwrites in place."""
    directory = directory or (PANELS_DIR / "_snapshots")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d")
    return directory / f"{panel.id}.{stamp}.yaml"


def write_snapshot(panel: Panel, verification: dict[str, Any],
                   directory: Path | None = None) -> Path:
    """Write a dated, hashed snapshot carrying the API resolution alongside membership."""
    path = snapshot_path(panel, directory)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = {
        "schema_version": 1,
        "snapshot_of": panel.id,
        "snapshot_utc": verification["checked_utc"],
        "source_file": relative(panel.source),
        "membership_digest": panel.digest(),
        "panels": {
            panel.id: {
                "description": panel.description,
                "species": panel.species,
                "provenance": {**panel.provenance,
                               "verified_against": ENSEMBL_REST,
                               "verified_utc": verification["checked_utc"]},
                "genes": dict(panel.genes),
                "minimum_present": panel.minimum_present,
                "ensembl": verification["resolved"],
            }
        },
    }
    path.write_text(yaml.safe_dump(body, sort_keys=False, width=88))
    return path


def diff_panels(left: Panel, right: Panel) -> dict[str, Any]:
    """Set and direction differences between two panels."""
    ls, rs = set(left.genes), set(right.genes)
    flipped = sorted(g for g in ls & rs if left.genes[g] != right.genes[g])
    return {
        "left": left.id,
        "right": right.id,
        "shared": sorted(ls & rs),
        "left_only": sorted(ls - rs),
        "right_only": sorted(rs - ls),
        "direction_conflicts": flipped,
        "identical": not (ls ^ rs) and not flipped,
    }
