"""Audit a clean public release tree against its file manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

import pymupdf

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "release_manifest.json"
LOCAL_DIRECTORIES = {".git", ".Rlib", ".venv", "__pycache__"}
PRIVATE_PATTERNS = [
    re.compile(r"/(?:Users|home)/[A-Za-z0-9_.-]+/"),
    re.compile(r"/private/(?:tmp)/"),
    re.compile(r"/var/(?:folders)/"),
    re.compile(r"[A-Za-z0-9._%+-]+@(?:ucla\.edu|gmail\.com)"),
    re.compile(r"gh[pousr]_[A-Za-z0-9]{20,}"),
    re.compile(r"(?i)-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    manifest = json.loads(MANIFEST.read_text())
    expected = manifest["files"]
    actual = {
        path.relative_to(ROOT).as_posix(): path
        for path in ROOT.rglob("*")
        if path.is_file() and path != MANIFEST
        and not (set(path.relative_to(ROOT).parts) & LOCAL_DIRECTORIES)
    }
    if set(actual) != set(expected):
        raise SystemExit(f"Release inventory differs: {sorted(set(actual) ^ set(expected))}")
    for relative, path in actual.items():
        assert not path.is_symlink(), f"Symlink in release: {relative}"
        assert digest(path) == expected[relative]["sha256"], f"Changed file: {relative}"
        assert path.stat().st_size == expected[relative]["bytes"], f"Size mismatch: {relative}"
        assert path.suffix.lower() not in {".docx", ".rds", ".xlsx", ".gz", ".zip"}, relative
        assert not relative.startswith(("data/raw/", "data/interim/", "data/processed/", "outputs/", "logs/")), relative
        if path.suffix in {".md", ".R", ".py", ".json", ".lock", ".txt", ".html", ".cff", ".csv"} or not path.suffix:
            text = path.read_text(encoding="utf-8")
            for pattern in PRIVATE_PATTERNS:
                assert not pattern.search(text), f"Review a sensitive-path pattern in {relative}"
        if path.suffix == ".pdf":
            with pymupdf.open(path) as document:
                assert not document.get_xml_metadata(), f"XMP remains in {relative}"
                assert document.xref_get_key(-1, "Info")[0] == "null", f"Info remains in {relative}"
                assert document.xref_get_key(-1, "ID")[0] == "null", f"Document ID remains in {relative}"
                assert all(not value for key, value in document.metadata.items()
                           if key not in {"format", "encryption"}), f"Properties remain in {relative}"
    print(f"Public release audit passed: {len(actual)} files; hashes, inventory, privacy patterns, and PDF metadata checked.")


if __name__ == "__main__":
    main()
