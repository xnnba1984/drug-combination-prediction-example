# Data Sources

The primary outcome is the published DeepSynergy Loewe synergy score. DrugComb supplies combination sensitivity score (CSS) records and identifiers. Its recalculated Loewe score is retained for a reconciliation check; the models use the published DeepSynergy outcome throughout.

## Fixed Sources

| Input | Fixed version | Access |
| :--- | :--- | :--- |
| DeepSynergy labels and five folds | `tastanlab/ohe`, commit `4d3368d0c004949e319842c56f7900fe3ccd5f4b` | [labels.csv](https://raw.githubusercontent.com/tastanlab/ohe/4d3368d0c004949e319842c56f7900fe3ccd5f4b/DeepSynergy/labels.csv) |
| Published DeepSynergy feature matrix | `ebcandir/DeepSynergy`, revision `2cc715df1e91d0a42269cb2f58ab32fce30ad2dd` | [X.p.gz](https://huggingface.co/datasets/ebcandir/DeepSynergy/resolve/2cc715df1e91d0a42269cb2f58ab32fce30ad2dd/X.p.gz) |
| DrugComb summary | Version 1.5, Zenodo record 15235991 | [Source record](https://zenodo.org/records/15235991) |
| DrugComb drug identifiers | Version 1.4, Zenodo record 18449193 | [Source record](https://zenodo.org/records/18449193) |

The acquisition script verifies the full DrugComb archive against its published MD5 value, `cd3e24758c10ffaea0124b6e80626c32`. It extracts records with `study_name` equal to `ONEIL`, checks the subset's SHA-256 hash, then removes its downloaded full archive. The identifier workbook is verified against MD5 `06894e63e5e06baa23d2bd59d6bd8808`.

The feature matrix is a NumPy pickle compressed with gzip. Its SHA-256 hash is checked before it is opened with a restricted unpickler. Extraction checks the matrix shape, the drug-order copies, and the consistency of repeated profiles. [input_checksums.json](input_checksums.json) records SHA-256 hashes for all seven inputs to the statistical analysis.

## Reuse And Redistribution

Both DrugComb records declare **CC BY 4.0** in their Zenodo metadata. Attribute DrugComb and the underlying studies when reusing these data. The analysis filters the v1.5 archive, reconciles identifiers, and averages matched CSS records for each analysis triplet.

No explicit redistribution license was found for the pinned label repository or the feature-matrix mirror during the source check on October 23, 2026. This package therefore provides access links and reconstruction scripts, with raw files and record-level derived data kept out of the release. Users obtain those inputs from the source hosts under the applicable terms. The code license applies to this repository's original code and documentation; it does not grant rights to third-party data.

## Generated Locally

- `raw/`: downloaded labels, the compressed feature matrix, the filtered O'Neil records, and the identifier workbook.
- `interim/`: separate drug and cell line feature tables, harmonized identifiers, and model checkpoints.
- `processed/`: the 22,737 analysis triplets, feature matrices, and out-of-fold predictions.

These directories are ignored by Git. The released summaries in `reference_results/` report model performance and figure-level aggregates.

## Key Fields

| Field | Meaning |
| :--- | :--- |
| `triplet_id` | Canonical unordered drug pair and cell line |
| `pair_id` | Canonical unordered drug pair |
| `outer_fold` | Published pair-level fold, numbered 0 through 4 |
| `loewe` | Mean of the published DeepSynergy labels for the triplet |
| `replicate_count` | Number of label records contributing to that mean |
| `css` | Mean matched DrugComb CSS, used for retrospective evaluation |
| `drugcomb_loewe` | Recalculated DrugComb Loewe score, used in a reconciliation diagnostic |
