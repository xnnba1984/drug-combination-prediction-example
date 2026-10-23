# Reproduction

## Entry Points

`run_figures.R` reads the released summary tables and writes the two-panel chapter figure and separate panels. It uses no raw inputs. This is useful for checking labels, sizing, and plotting code.

`run_all.R` rebuilds the analysis from locally acquired inputs. It reruns data preparation, PCA, nested model fitting, evaluation, bootstrap intervals, and pair-block permutations. Figures and checks are then regenerated. Existing model checkpoints are overwritten by new fits and are never used to skip a fold.

Run `scripts/acquire_sources.py` before the full analysis. Every source download is checked against a fixed hash. The full runner independently verifies all seven source and extracted input files before fitting models.

## Environment

- R 4.6.1; 30 locked packages, including `glmnet` 5.0, `ranger` 0.18.0, `data.table` 1.18.4, and `ggplot2` 4.0.3.
- Python 3.12; exact package versions are listed in `requirements.txt`.
- Eight random forest threads by default. Set `SECTION3_THREADS` before starting R to use a different count.

The clean test used macOS on Apple silicon. Linux and Windows follow the same relative-path entry points; separate execution on those operating systems has not been tested. Numerical details may vary with the compiler, numerical libraries, or platform. The manuscript checks expose changes that affect the reported values.

`setup.R` restores packages into `.Rlib/` without using another project's package library. It recognizes an installed Homebrew `gettext` on macOS if a source build needs `libintl.h`. A missing compiler or system library must be installed through the tools appropriate to the operating system. R packages that ship with R must match the lockfile; using R 4.6.1 is the simplest way to reproduce the tested setup.

## Checks

| Stage | Checks |
| :--- | :--- |
| Acquisition | Fixed source hashes, full DrugComb archive checksum, and filtered subset hash |
| Feature extraction | Array dimensions, finite values, swapped drug blocks, and repeated entity profiles |
| Data preparation | Canonical pair identifiers, duplicate aggregation, pair-level folds, and complete CSS matching |
| Features | Drug-order invariance and coverage of every drug and cell line in each training fold |
| Fitting | Aligned input rows, expected design dimensions, and no outcome or CSS column in the design matrices |
| Final QC | Sixteen checks of counts, folds, predictions, bootstrap sizes, permutations, and candidate lists |
| Manuscript | Forty-four checks of reported numbers and prediction-design constraints |
| PDF output | Empty Info and XMP metadata, with identical page text and rendered pixels after cleanup |

`outputs/tables/manuscript_checks.csv` compares results against the values reported in the chapter. Tolerances reflect the displayed precision. Counts and structural constraints must match exactly. Failure stops the run, preserving the check table for inspection.

Small helper tests cover checksum validation, source filtering, restricted pickle loading, and PDF cleanup:

```sh
.venv/bin/python -m unittest discover -s tests -v
```

The released tables retain numerical precision. The README and chapter round values for reading. The weighted fold enrichment is checked directly; it is not reconstructed by dividing the displayed percentages. The 8.55 mean difference is calculated over all 60 selected pairs.

## Outputs

Full analysis results are written under `outputs/`, with processed data under `data/processed/`. `logs/run_all.log` records elapsed time for each stage, and `logs/environment.txt` records R and package versions without personal file paths. The scripts stop on failed assertions.

The release contains a small, fixed set of reference summaries. Model-level metrics and the 60 plotted pair aggregates support the chapter figure. Raw outcomes, individual predictions, and extracted entity features are generated locally and excluded from the public release.

## Verified Run

On October 23, 2026, the full analysis was run in a separate directory after fresh source acquisition and package restoration. The environment used R 4.6.1 and Python 3.12.14 on macOS with Apple silicon. Statistical analysis took 546.3 seconds, excluding package installation and data acquisition. The R process reached about 4 GB of resident memory.

All five sets of predictions matched the reference analysis exactly, and 26 existing result tables were byte-identical. The 16 final analysis checks and 44 manuscript checks passed. Six helper tests also passed. A separate environment redrew the chapter figure using only the released summaries; its PNG pixels matched the full-run figure exactly.

The HTML page was checked at 390, 1,280, and 1,920 pixel widths. The published figure PDF was checked after metadata removal. Its text and rendered pixels were unchanged.

To audit the public package itself, run:

```sh
.venv/bin/python scripts/check_release.py
```

This verifies the release manifest and checks file types, sensitive-path patterns, and the metadata of the released figure PDF. Run it on the clean release tree before adding local data or environment directories.
