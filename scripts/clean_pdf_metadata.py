"""Clear nonessential PDF properties and verify unchanged page content."""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys

import pymupdf


def page_signature(document):
    return [
        (page.get_text(), tuple(page.rect), hashlib.sha256(page.get_pixmap().samples).hexdigest())
        for page in document
    ]


def clean(path):
    temporary = path.with_suffix(".metadata-clean.tmp.pdf")
    with pymupdf.open(path) as document:
        before = page_signature(document)
        document.set_metadata({})
        document.del_xml_metadata()
        document.xref_set_key(-1, "Info", "null")
        document.xref_set_key(-1, "ID", "null")
        document.save(temporary, garbage=4, deflate=True, no_new_id=True)
    with pymupdf.open(temporary) as document:
        assert page_signature(document) == before, f"PDF content changed: {path}"
        assert not document.get_xml_metadata(), f"XMP remains: {path}"
        assert document.xref_get_key(-1, "Info")[0] == "null", f"Info remains: {path}"
        assert document.xref_get_key(-1, "ID")[0] == "null", f"Document ID remains: {path}"
        assert all(not value for key, value in document.metadata.items()
                   if key not in {"format", "encryption"}), f"Properties remain: {path}"
    temporary.replace(path)
    print(f"PDF metadata cleared; page text and pixels unchanged: {path.name}")


def main():
    paths = []
    for name in sys.argv[1:]:
        path = Path(name)
        paths.extend(sorted(path.glob("*.pdf")) if path.is_dir() else [path])
    if not paths:
        raise SystemExit("Provide a PDF file or a directory containing PDFs.")
    for path in paths:
        clean(path)


if __name__ == "__main__":
    main()
