#!/usr/bin/env python3
"""Download analysis inputs listed in a manifest, verify sizes, and enforce registration gates.

Manifests are TSVs with columns ``target_path``, ``url`` (or ``osdr:<OSD number>:<file name>``,
resolved through the OSDR file API), ``size_bytes`` (blank = unchecked) and, optionally, ``group``
and ``gate``. A row with a ``gate`` path is downloaded only if that file is tracked and committed
in git, so directional predictions are registered before unseen data are fetched.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OSDR_FILES_API = "https://osdr.nasa.gov/osdr/data/osd/files/{number}"
OSD771_METADATA_DIR = "data/raw/metadata"


def resolve(url: str, cache: dict[str, dict[str, str]]) -> str:
    """Turn ``osdr:<number>:<file name>`` into a download URL via the OSDR file API."""
    if not url.startswith("osdr:"):
        return url
    _, number, name = url.split(":", 2)
    if number not in cache:
        with urllib.request.urlopen(OSDR_FILES_API.format(number=number), timeout=120) as response:
            listing = json.load(response)
        cache[number] = {
            item["file_name"]: "https://osdr.nasa.gov" + item["remote_url"]
            for study in listing["studies"].values()
            for item in study["study_files"]
        }
    if name not in cache[number]:
        raise FileNotFoundError(f"OSD-{number} has no file named {name}")
    return cache[number][name]


def gate_committed(root: Path, gate: str) -> bool:
    tracked = subprocess.run(
        ["git", "ls-files", "--error-unmatch", gate], cwd=root, capture_output=True, text=True
    )
    dirty = subprocess.run(
        ["git", "status", "--porcelain", "--", gate], cwd=root, capture_output=True, text=True
    )
    return tracked.returncode == 0 and dirty.stdout.strip() == ""


def download(url: str, dest: Path, attempts: int = 4) -> None:
    """Fetch ``url`` to ``dest`` atomically, retrying transient failures with exponential backoff."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=600) as response, tmp.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1 << 20)
            tmp.replace(dest)
            return
        except Exception:
            if attempt == attempts - 1:
                raise
            time.sleep(2 ** (attempt + 1))


def extract_isa(root: Path, archive: Path) -> None:
    """Unpack an ISA bundle next to it; OSD-771's sample/assay tables also go to data/raw/metadata."""
    out = archive.parent / "metadata"
    out.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.infolist():
            name = Path(member.filename).name
            if member.is_dir() or not name or name.startswith("."):
                continue
            (out / name).write_bytes(bundle.read(member))
    if "OSD-771" in archive.name:
        legacy = root / OSD771_METADATA_DIR
        legacy.mkdir(parents=True, exist_ok=True)
        for table in out.iterdir():
            if table.name.startswith(("s_", "a_", "i_")):
                shutil.copyfile(table, legacy / table.name)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=REPO / "resources/osdr_inputs.tsv")
    parser.add_argument("--group", action="append", help="only rows whose group matches (repeatable)")
    parser.add_argument("--root", type=Path, default=REPO)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--include-supplementary", action="store_true",
                        help="also fetch rows whose role is 'supplementary' (MultiQC bundles etc.)")
    args = parser.parse_args()

    with args.manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if args.group:
        rows = [row for row in rows if row.get("group") in set(args.group)]
    if not args.include_supplementary:
        rows = [row for row in rows if row.get("role") != "supplementary"]

    cache: dict[str, dict[str, str]] = {}
    failures = 0
    for row in rows:
        target = args.root / row["target_path"]
        size = int(row["size_bytes"]) if row.get("size_bytes") else None
        gate = (row.get("gate") or "").strip()
        if gate and not gate_committed(args.root, gate):
            print(f"GATED        {row['target_path']}: commit {gate} first")
            failures += 1
            continue
        if target.exists() and (size is None or target.stat().st_size == size):
            print(f"ok (present) {row['target_path']}")
            continue
        if args.dry_run:
            print(f"would fetch  {row['target_path']}")
            continue
        try:
            download(resolve(row["url"], cache), target)
        except Exception as error:  # network or API failure: report and continue
            print(f"FAILED       {row['target_path']}: {error}")
            failures += 1
            continue
        if size is not None and target.stat().st_size != size:
            print(f"SIZE MISMATCH {row['target_path']}: {target.stat().st_size} != {size}")
            failures += 1
            continue
        if target.name.endswith("ISA.zip"):
            extract_isa(args.root, target)
        print(f"fetched      {row['target_path']}")
    print(f"done: {len(rows) - failures}/{len(rows)} rows ok")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
