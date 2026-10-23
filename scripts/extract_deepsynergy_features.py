#!/usr/bin/env python3
"""Extract compact entity-level features from the published DeepSynergy matrix."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import pickle
from pathlib import Path

import numpy as np


N_DRUG_FEATURES = 4_387
N_CELL_FEATURES = 3_984
EXPECTED_COLUMNS = 2 * N_DRUG_FEATURES + N_CELL_FEATURES


class RestrictedUnpickler(pickle.Unpickler):
    """Allow only the NumPy objects required by the published array."""

    ALLOWED = {
        ("numpy.core.multiarray", "_reconstruct"): np.core.multiarray._reconstruct,
        ("numpy", "ndarray"): np.ndarray,
        ("numpy", "dtype"): np.dtype,
    }

    def find_class(self, module: str, name: str):
        key = (module, name)
        if key not in self.ALLOWED:
            raise pickle.UnpicklingError(f"Blocked pickle global: {module}.{name}")
        return self.ALLOWED[key]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--labels", required=True, type=Path)
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--drug-output", required=True, type=Path)
    parser.add_argument("--cell-output", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    return parser.parse_args()


def read_labels(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"drug_a_name", "drug_b_name", "cell_line", "synergy", "fold"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("Label file is missing required columns")
        return list(reader)


def max_entity_deviation(
    labels: list[dict[str, str]],
    matrix: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], float, float]:
    drugs: dict[str, np.ndarray] = {}
    cells: dict[str, np.ndarray] = {}
    max_drug_deviation = 0.0
    max_cell_deviation = 0.0

    for index, row in enumerate(labels):
        drug_a = matrix[index, :N_DRUG_FEATURES]
        drug_b = matrix[index, N_DRUG_FEATURES : 2 * N_DRUG_FEATURES]
        cell = matrix[index, 2 * N_DRUG_FEATURES :]

        for name, values in ((row["drug_a_name"], drug_a), (row["drug_b_name"], drug_b)):
            if name in drugs:
                max_drug_deviation = max(
                    max_drug_deviation,
                    float(np.max(np.abs(values - drugs[name]))),
                )
            else:
                drugs[name] = values.copy()

        cell_name = row["cell_line"]
        if cell_name in cells:
            max_cell_deviation = max(
                max_cell_deviation,
                float(np.max(np.abs(cell - cells[cell_name]))),
            )
        else:
            cells[cell_name] = cell.copy()

    return drugs, cells, max_drug_deviation, max_cell_deviation


def write_feature_table(
    path: Path,
    id_name: str,
    entities: dict[str, np.ndarray],
    prefix: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = sorted(entities)
    feature_count = len(entities[names[0]])
    header = [id_name] + [f"{prefix}_{index:04d}" for index in range(1, feature_count + 1)]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for name in names:
            writer.writerow([name] + [f"{value:.12g}" for value in entities[name]])


def main() -> int:
    args = parse_args()
    labels = read_labels(args.labels)

    with gzip.open(args.matrix, "rb") as handle:
        matrix = RestrictedUnpickler(handle).load()

    expected_shape = (2 * len(labels), EXPECTED_COLUMNS)
    if not isinstance(matrix, np.ndarray) or matrix.shape != expected_shape:
        raise ValueError(f"Expected matrix shape {expected_shape}; found {getattr(matrix, 'shape', None)}")
    if not np.isfinite(matrix).all():
        raise ValueError("Feature matrix contains non-finite values")

    n_rows = len(labels)
    first = matrix[:n_rows]
    second = matrix[n_rows:]
    swap_a_to_b = float(
        np.max(
            np.abs(
                first[:, :N_DRUG_FEATURES]
                - second[:, N_DRUG_FEATURES : 2 * N_DRUG_FEATURES]
            )
        )
    )
    swap_b_to_a = float(
        np.max(
            np.abs(
                first[:, N_DRUG_FEATURES : 2 * N_DRUG_FEATURES]
                - second[:, :N_DRUG_FEATURES]
            )
        )
    )
    cell_copy_deviation = float(
        np.max(np.abs(first[:, 2 * N_DRUG_FEATURES :] - second[:, 2 * N_DRUG_FEATURES :]))
    )
    if max(swap_a_to_b, swap_b_to_a, cell_copy_deviation) > 1e-12:
        raise ValueError("The second matrix half is not the expected Drug 2-Drug 1 copy")

    drugs, cells, drug_deviation, cell_deviation = max_entity_deviation(labels, first)
    if drug_deviation > 1e-12 or cell_deviation > 1e-12:
        raise ValueError("Entity features vary across repeated observations")

    write_feature_table(args.drug_output, "drug_name", drugs, "drug_feature")
    write_feature_table(args.cell_output, "cell_line", cells, "cell_feature")

    qc = {
        "matrix_shape": list(matrix.shape),
        "label_rows": n_rows,
        "drug_count": len(drugs),
        "cell_line_count": len(cells),
        "drug_feature_count": N_DRUG_FEATURES,
        "cell_feature_count": N_CELL_FEATURES,
        "swap_a_to_b_max_abs_deviation": swap_a_to_b,
        "swap_b_to_a_max_abs_deviation": swap_b_to_a,
        "cell_copy_max_abs_deviation": cell_copy_deviation,
        "within_drug_max_abs_deviation": drug_deviation,
        "within_cell_max_abs_deviation": cell_deviation,
    }
    args.qc_output.parent.mkdir(parents=True, exist_ok=True)
    args.qc_output.write_text(json.dumps(qc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(qc, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
