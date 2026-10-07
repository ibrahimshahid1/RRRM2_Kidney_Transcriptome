#!/usr/bin/env python3
"""Restore the frozen derived resources (gene map, tiers, KEGG, atlas pseudobulk) and verify sha256."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_BUNDLE = REPO / "resources/frozen/2026-10-06"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=DEFAULT_BUNDLE)
    parser.add_argument("--root", type=Path, default=REPO, help="repository root to restore into")
    parser.add_argument("--force", action="store_true", help="overwrite a target whose sha256 differs")
    args = parser.parse_args()

    rows = (args.bundle / "MANIFEST.tsv").read_text().splitlines()[1:]
    failures = 0
    for row in rows:
        packed, target, compression, expected = row.split("\t")
        dest = args.root / target
        if dest.exists():
            if sha256(dest) == expected:
                print(f"ok (present)  {target}")
                continue
            if not args.force:
                print(f"REFUSED      {target}: exists with a different sha256 (use --force)")
                failures += 1
                continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        source = args.bundle / packed
        tmp = dest.with_name(dest.name + ".restoring")
        if compression == "gzip":
            with gzip.open(source, "rb") as fi, tmp.open("wb") as fo:
                shutil.copyfileobj(fi, fo)
        else:
            shutil.copyfile(source, tmp)
        if sha256(tmp) != expected:
            tmp.unlink()
            print(f"FAILED       {target}: restored sha256 does not match MANIFEST")
            failures += 1
            continue
        tmp.replace(dest)
        print(f"restored     {target}")
    print(f"{len(rows) - failures}/{len(rows)} resources verified")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
