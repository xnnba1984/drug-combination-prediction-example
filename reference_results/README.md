# Reference Results

These summaries accompany the worked example in Section 3.6. They support Figure 3.1 and provide a reference for a complete rerun.

## Main Files

| File | Contents |
| :--- | :--- |
| `figures/worked_example_results.png` | Chapter figure, 6.5 by 2.65 inches at 320 dpi |
| `figures/worked_example_results.pdf` | Vector version with nonessential PDF metadata cleared |
| `tables/model_performance_overall.csv` | RMSE, mean absolute error, and correlations for five models |
| `tables/model_performance_pair_bootstrap_ci.csv` | Estimates and 95% pair-bootstrap intervals |
| `tables/chapter_pair_points.csv` | The 60 pair-level aggregates shown in panel B |
| `tables/pair_calibration_summary.csv` | Count above prediction and mean difference over all 60 pairs |
| `tables/triplet_top_decile_overall.csv` | Ranking performance and weighted fold enrichment |
| `tables/pair_candidate_drug_frequency.csv` | Drug frequencies among the 60 selected pairs |
| `tables/identity_ridge_negative_control_summary.csv` | Summary of 20 pair-block permutations |
| `tables/final_qc_checks.csv` | Sixteen analysis checks |
| `tables/manuscript_checks.csv` | Forty-four checks of manuscript values and design constraints |

Other tables document fold coverage, PCA, fold-specific performance, and endpoint reconciliation. Computed summaries are included here; raw labels, full prediction tables, and individual entity features are reconstructed locally.

## Panel B Fields

Each row of `chapter_pair_points.csv` represents one selected pair. `pair_prediction` is the mean prediction for its three highest-ranked cell lines. `selected_observed_loewe` and `selected_css` summarize the same three cell lines after ranking. `observed_minus_predicted` is the difference reported in the calibration summary. The mean difference is taken over all 60 rows.

The quick figure command reads these values without refitting models. A full run produces new results under `outputs/`; it leaves this reference directory unchanged.
