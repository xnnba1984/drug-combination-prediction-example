import contextlib
import csv
import hashlib
import io
import json
from pathlib import Path
import pickle
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from acquire_sources import checksum
from clean_pdf_metadata import clean, page_signature
from extract_deepsynergy_features import RestrictedUnpickler
import verify_inputs


class HelperTests(unittest.TestCase):
    def test_checksum(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            path.write_bytes(b"abc")
            self.assertEqual(checksum(path), hashlib.sha256(b"abc").hexdigest())

    def test_restricted_numpy_pickle(self):
        array = np.array([[1.5, 2.5]])
        restored = RestrictedUnpickler(io.BytesIO(pickle.dumps(array, protocol=4))).load()
        np.testing.assert_array_equal(restored, array)

    def test_restricted_pickle_rejects_other_globals(self):
        with self.assertRaises(pickle.UnpicklingError):
            RestrictedUnpickler(io.BytesIO(pickle.dumps(eval, protocol=4))).load()

    def test_source_filter_preserves_header_and_selected_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.csv"
            source.write_text("study_name,value\nONEIL,7\nOTHER,9\nONEIL,11\n")
            output, metadata = root / "filtered.csv", root / "filter.json"
            subprocess.run([sys.executable, str(ROOT / "scripts/filter_drugcomb_v15.py"),
                            "--input", str(source), "--study", "ONEIL", "--output", str(output),
                            "--metadata", str(metadata)], check=True, capture_output=True, text=True)
            with output.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows, [{"study_name": "ONEIL", "value": "7"},
                                    {"study_name": "ONEIL", "value": "11"}])
            record = json.loads(metadata.read_text())
            self.assertEqual(record["rows_read"], 3)
            self.assertEqual(record["rows_written"], 2)
            self.assertEqual(record["output_sha256"], checksum(output))

    def test_input_verification_rejects_modified_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "data").mkdir()
            path = root / "data/input.csv"
            path.write_bytes(b"original")
            (root / "data/input_checksums.json").write_text(json.dumps({
                "algorithm": "sha256", "files": {"data/input.csv": checksum(path)}}))
            with patch.object(verify_inputs, "ROOT", root), contextlib.redirect_stdout(io.StringIO()):
                verify_inputs.main()
                path.write_bytes(b"modified")
                with self.assertRaisesRegex(SystemExit, "Checksum mismatch"):
                    verify_inputs.main()

    def test_pdf_cleanup_preserves_content(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "figure.pdf"
            with pymupdf.open() as document:
                page = document.new_page(width=240, height=160)
                page.insert_text((25, 60), "Prediction results")
                document.set_metadata({"title": "Temporary property", "producer": "Temporary generator"})
                document.set_xml_metadata('<x:xmpmeta xmlns:x="adobe:ns:meta/"><value>temporary</value></x:xmpmeta>')
                document.save(path)
            with pymupdf.open(path) as document:
                before = page_signature(document)
            with contextlib.redirect_stdout(io.StringIO()):
                clean(path)
            with pymupdf.open(path) as document:
                self.assertEqual(page_signature(document), before)
                self.assertEqual(document.get_xml_metadata(), "")
                self.assertEqual(document.xref_get_key(-1, "Info")[0], "null")
                self.assertEqual(document.xref_get_key(-1, "ID")[0], "null")


if __name__ == "__main__":
    unittest.main()
