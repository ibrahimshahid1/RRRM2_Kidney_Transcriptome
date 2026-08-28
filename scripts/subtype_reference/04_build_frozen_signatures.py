#!/usr/bin/env python3
"""CLI for flight-blind distal-nephron signature construction."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from src.subtype_reference.reference_builder import build_signatures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--discovery")
    parser.add_argument("--validation")
    parser.add_argument("--segment-validation")
    parser.add_argument("--atlas")
    args = parser.parse_args()
    gate = build_signatures(
        config_path=args.config,
        output_dir=args.output_dir,
        discovery_path=args.discovery,
        validation_path=args.validation,
        segment_validation_path=args.segment_validation,
        atlas_path=args.atlas,
    )
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
