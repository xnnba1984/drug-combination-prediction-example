# Third-Party Sources And Dependencies

The MIT license in this repository covers its original analysis code and documentation. It does not replace the terms of downloaded datasets or installed software dependencies.

## Data

DrugComb v1.5 and the v1.4 identifier record declare CC BY 4.0. The pinned DeepSynergy label repository and feature mirror had no explicit redistribution license in the checked source records. Input files are acquired from those hosts and excluded from this release. See [data/README.md](data/README.md) for versions, attribution, and processing details.

## Software

Dependencies are installed from their distribution sources; their code and binary packages are not bundled here. The following license strings were checked in the installed package metadata:

| Package | Version | Declared license |
| :--- | :--- | :--- |
| glmnet | 5.0 | GPL-2 |
| ranger | 0.18.0 | GPL-3 |
| data.table | 1.18.4 | MPL-2.0 or the package LICENSE file |
| ggplot2 | 4.0.3 | MIT with package LICENSE file |
| renv | 1.2.4 | MIT with package LICENSE file |
| PyMuPDF | 1.27.1 | GNU Affero GPL 3.0 or an Artifex commercial license |

PyMuPDF is used to clean and check figure PDFs. NumPy, pandas, and openpyxl support input preparation. Consult each dependency's distributed license for its full terms when reusing or redistributing software built around this workflow.
