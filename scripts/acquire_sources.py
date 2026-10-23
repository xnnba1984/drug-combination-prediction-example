"""Download pinned source files, check their hashes, and extract analysis inputs."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys
import time
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
LABEL_URL = "https://raw.githubusercontent.com/tastanlab/ohe/4d3368d0c004949e319842c56f7900fe3ccd5f4b/DeepSynergy/labels.csv"
FEATURE_URL = "https://huggingface.co/datasets/ebcandir/DeepSynergy/resolve/2cc715df1e91d0a42269cb2f58ab32fce30ad2dd/X.p.gz"
DRUGCOMB_URL = "https://zenodo.org/api/records/15235991/files/summary_v_1_5.csv/content"
IDENTIFIER_URL = "https://zenodo.org/api/records/18449193/files/DrugComb_drug_identifiers.xlsx/content"
LABEL_SHA = "4ef40678fff51faa84fa4de63344c5b51e2e4abe7c28a417af24ec2b50fa135a"
FEATURE_SHA = "e399f6cdb40184048ac34b370f89cc6e4582bec72893cc2c98d90c1634517899"
SUBSET_SHA = "6bf6e7e434e4866d48c4e455f5a1039ed4bc4cf355c1f8022bcef71c9ea96349"


def checksum(path, algorithm="sha256"):
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url, relative, expected, algorithm="sha256"):
    path = ROOT / relative
    if path.exists() and checksum(path, algorithm) == expected:
        print(f"Verified {relative}", flush=True)
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".download")
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "Section3WorkedExample/1.0"})
            with urllib.request.urlopen(request, timeout=120) as response, temporary.open("wb") as stream:
                received = 0
                last_update = time.monotonic()
                for block in iter(lambda: response.read(1024 * 1024), b""):
                    stream.write(block)
                    received += len(block)
                    if time.monotonic() - last_update >= 10:
                        print(f"Downloading {path.name}: {received / 1024**2:.1f} MiB", flush=True)
                        last_update = time.monotonic()
            if checksum(temporary, algorithm) != expected:
                raise ValueError(f"{algorithm} mismatch for {relative}")
            temporary.replace(path)
            print(f"Downloaded and verified {relative}", flush=True)
            return path
        except Exception:
            temporary.unlink(missing_ok=True)
            if attempt == 2:
                raise
            print(f"Retrying {path.name}", flush=True)
            time.sleep(2)


def helper(script, *args):
    subprocess.run([sys.executable, f"scripts/{script}", *args], cwd=ROOT, check=True)


def main():
    for directory in ("data/raw", "data/interim", "logs"):
        (ROOT / directory).mkdir(parents=True, exist_ok=True)
    download(LABEL_URL, "data/raw/deepsynergy_labels.csv", LABEL_SHA)
    download(FEATURE_URL, "data/raw/deepsynergy_X.p.gz", FEATURE_SHA)
    subset = ROOT / "data/raw/drugcomb_v15_oneil.csv"
    if not subset.exists() or checksum(subset) != SUBSET_SHA:
        archive = download(DRUGCOMB_URL, "data/raw/summary_v_1_5.csv", "cd3e24758c10ffaea0124b6e80626c32", "md5")
        helper("filter_drugcomb_v15.py", "--input", str(archive), "--study", "ONEIL", "--output", "data/raw/drugcomb_v15_oneil.csv", "--metadata", "logs/drugcomb_v15_filter_metadata.json")
        if checksum(subset) != SUBSET_SHA:
            raise ValueError("Filtered DrugComb checksum differs from the reference subset")
        archive.unlink()
    else:
        print("Verified data/raw/drugcomb_v15_oneil.csv", flush=True)
    download(IDENTIFIER_URL, "data/raw/DrugComb_drug_identifiers_v1.4.xlsx", "06894e63e5e06baa23d2bd59d6bd8808", "md5")
    helper("extract_deepsynergy_features.py", "--labels", "data/raw/deepsynergy_labels.csv", "--matrix", "data/raw/deepsynergy_X.p.gz", "--drug-output", "data/interim/drug_features.csv", "--cell-output", "data/interim/cell_features.csv", "--qc-output", "logs/feature_extraction_qc.json")
    helper("extract_drugcomb_identifiers.py", "--workbook", "data/raw/DrugComb_drug_identifiers_v1.4.xlsx", "--oneil", "data/raw/drugcomb_v15_oneil.csv", "--output", "data/interim/drugcomb_oneil_identifiers.csv")
    print("All inputs prepared from the pinned sources.", flush=True)


if __name__ == "__main__":
    main()
