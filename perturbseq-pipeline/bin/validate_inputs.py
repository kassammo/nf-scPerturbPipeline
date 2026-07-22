#!/usr/bin/env python3
"""Validate the samplesheet and analysis configuration.

Author: Mohamed Kassam
Date: 2026-07-19
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

import yaml

BASE_REQUIRED = {"sample_id", "h5_path"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--h5-list", required=True, type=Path)
    args = parser.parse_args()

    with args.config.open(encoding="utf-8") as handle:
        config = yaml.safe_load(handle) or {}

    configured_columns = set(config.get("metadata", {}).get("required_columns", []))
    required = BASE_REQUIRED | configured_columns

    with args.samplesheet.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    missing = required - set(fieldnames)
    if missing:
        sys.exit(f"Missing columns: {sorted(missing)}")
    if not rows:
        sys.exit("Samplesheet is empty")

    ids = [row["sample_id"].strip() for row in rows]
    if any(not sample_id for sample_id in ids):
        sys.exit("sample_id values cannot be empty")
    if len(ids) != len(set(ids)):
        sys.exit("sample_id values must be unique")

    for row in rows:
        h5_path = os.path.expandvars(os.path.expanduser(row["h5_path"].strip()))
        if not os.path.isfile(h5_path):
            sys.exit(f"Missing h5_path: {h5_path}")
        row["h5_path"] = os.path.abspath(h5_path)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    args.h5_list.parent.mkdir(parents=True, exist_ok=True)
    with args.h5_list.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(row["h5_path"] for row in rows) + "\n")

    print(f"Validated {len(rows)} samples")


if __name__ == "__main__":
    main()
