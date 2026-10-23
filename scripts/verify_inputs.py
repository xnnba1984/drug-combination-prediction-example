"""Verify every downloaded and extracted input before fitting models."""

import json

from acquire_sources import ROOT, checksum


def main():
    manifest = json.loads((ROOT / "data/input_checksums.json").read_text())
    for relative, expected in manifest["files"].items():
        path = ROOT / relative
        if not path.is_file():
            raise SystemExit(f"Missing {relative}; run scripts/acquire_sources.py first.")
        if checksum(path, manifest["algorithm"]) != expected:
            raise SystemExit(f"Checksum mismatch: {relative}")
        print(f"Verified {relative}")


if __name__ == "__main__":
    main()
