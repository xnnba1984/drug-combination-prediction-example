#!/usr/bin/env python3
"""Stream the DrugComb v1.5 summary file and retain one study."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import sys
import urllib.request
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--url")
    source.add_argument("--input", type=Path)
    parser.add_argument("--study", default="ONEIL")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metadata.parent.mkdir(parents=True, exist_ok=True)

    request = urllib.request.Request(args.url, headers={"User-Agent": "Section3WorkedExample/1.0"}) if args.url else None
    rows_read = 0
    rows_written = 0

    with (urllib.request.urlopen(request, timeout=120) if request else args.input.open("rb")) as response:
        text_stream = io.TextIOWrapper(response, encoding="utf-8", newline="")
        reader = csv.DictReader(text_stream)
        if reader.fieldnames is None or "study_name" not in reader.fieldnames:
            raise ValueError("DrugComb file does not contain study_name")

        with args.output.open("w", encoding="utf-8", newline="") as output_handle:
            writer = csv.DictWriter(output_handle, fieldnames=reader.fieldnames)
            writer.writeheader()
            for row in reader:
                rows_read += 1
                if row.get("study_name", "").strip().upper() == args.study.upper():
                    writer.writerow(row)
                    rows_written += 1
                if rows_read % 1_000_000 == 0:
                    print(
                        f"Read {rows_read:,} rows; retained {rows_written:,}",
                        file=sys.stderr,
                        flush=True,
                    )

    digest = hashlib.sha256()
    with args.output.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)

    metadata = {
        "source_url": args.url or "https://zenodo.org/api/records/15235991/files/summary_v_1_5.csv/content",
        "study": args.study,
        "rows_read": rows_read,
        "rows_written": rows_written,
        "output": str(args.output),
        "output_sha256": digest.hexdigest(),
    }
    args.metadata.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
