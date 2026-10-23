suppressPackageStartupMessages(library(data.table))
source("R/utils.R")

outcomes <- readRDS("data/processed/oneil_analysis_outcomes.rds")
performance <- fread("outputs/tables/model_performance_overall.csv")
intervals <- fread("outputs/tables/model_performance_pair_bootstrap_ci.csv")
pca <- fread("outputs/tables/pca_summary.csv")
ranking <- fread("outputs/tables/triplet_top_decile_overall.csv")[model == "Identity ridge"]
points <- fread("outputs/tables/chapter_pair_points.csv")
calibration <- fread("outputs/tables/pair_calibration_summary.csv")
permutations <- fread("outputs/tables/identity_ridge_pair_permutation_control.csv")
frequency <- fread("outputs/tables/pair_candidate_drug_frequency.csv")
coverage <- fread("outputs/tables/outer_fold_coverage.csv")
predictors <- fread("outputs/tables/model_predictor_columns.csv")
identity <- performance[model == "Identity ridge"]
identity_ci <- intervals[model == "Identity ridge" & metric == "rmse"]
example <- points[drug_1 == "BEZ-235" & drug_2 == "DASATINIB"]
assert_true(nrow(example) == 1L, "Expected example pair is missing or duplicated")

checks <- list()
check_value <- function(name, observed, expected, tolerance = 0) {
  passed <- length(observed) == 1L && is.finite(observed) &&
    abs(observed - expected) <= tolerance + 1e-12
  checks[[length(checks) + 1L]] <<- data.table(
    check = name, observed = observed, expected = expected,
    tolerance = tolerance, status = if (passed) "PASS" else "FAIL"
  )
}

check_value("Published label rows", sum(outcomes$replicate_count), 23052)
check_value("Unique triplets", nrow(outcomes), 22737)
check_value("Repeated triplets", sum(outcomes$replicate_count == 2L), 315)
check_value("Unique pairs", uniqueN(outcomes$pair_id), 583)
check_value("Drugs", uniqueN(c(outcomes$drug_1, outcomes$drug_2)), 38)
check_value("Cell lines", uniqueN(outcomes$cell_line), 39)
check_value("Triplets with CSS", sum(is.finite(outcomes$css)), 22737)
check_value("Drug features", pca[entity == "Drug", original_features], 4387)
check_value("Cell line features", pca[entity == "Cell line", original_features], 3984)
check_value("Drug components", pca[entity == "Drug", retained_components], 23)
check_value("Cell line components", pca[entity == "Cell line", retained_components], 14)
check_value("Drug variance percent", 100 * pca[entity == "Drug", retained_variance], 80.7, 0.05)
check_value("Cell line variance percent", 100 * pca[entity == "Cell line", retained_variance], 80.9, 0.05)
check_value("Identity ridge RMSE", identity$rmse, 16.476, 0.0005)
check_value("Identity ridge RMSE lower 95% limit", identity_ci$lower_95, 15.608, 0.0005)
check_value("Identity ridge RMSE upper 95% limit", identity_ci$upper_95, 17.36, 0.005)
check_value("Identity ridge MAE", identity$mae, 10.552, 0.0005)
check_value("Identity ridge Pearson correlation", identity$pearson, 0.688, 0.0005)
check_value("Identity ridge Spearman correlation", identity$spearman, 0.693, 0.0005)
check_value("Training mean RMSE", performance[model == "Training mean", rmse], 22.739, 0.0005)
check_value("Feature elastic net RMSE", performance[model == "Feature elastic net", rmse], 19.815, 0.0005)
check_value("Feature random forest RMSE", performance[model == "Feature random forest", rmse], 19.637, 0.0005)
check_value("Hybrid elastic net RMSE", performance[model == "Hybrid elastic net", rmse], 16.551, 0.0005)
check_value("Selected triplet positive percent", 100 * ranking$positive_rate, 91, 0.05)
check_value("Baseline triplet positive percent", 100 * ranking$baseline_positive_rate, 34.4, 0.05)
check_value("Weighted fold enrichment", ranking$positive_enrichment, 2.68, 0.005)
check_value("Example pair prediction", example$pair_prediction, 81.6, 0.05)
check_value("Example observed mean Loewe", example$selected_observed_loewe, 114, 0.5)
check_value("Example observed mean CSS", example$selected_css, 69, 0.5)
check_value("Selected pairs", nrow(points), 60)
check_value("Observed above predicted pairs", sum(points$observed_minus_predicted > 0), 43)
check_value("Mean difference over all 60 pairs", mean(points$observed_minus_predicted), 8.55, 0.005)
check_value("Mean difference summary matches figure points", calibration$mean_observed_minus_predicted,
  mean(points$selected_observed_loewe - points$pair_prediction), 1e-10)
check_value("Permutation replicates", nrow(permutations), 20)
check_value("Permutation RMSE minimum", min(permutations$rmse), 22.526, 0.0005)
check_value("Permutation RMSE maximum", max(permutations$rmse), 22.715, 0.0005)
check_value("Permutation correlation minimum", min(permutations$pearson), 0.015, 0.0005)
check_value("Permutation correlation maximum", max(permutations$pearson), 0.126, 0.0005)
check_value("BEZ-235 candidate pairs", frequency[drug == "BEZ-235", candidate_pairs], 34)
check_value("DASATINIB candidate pairs", frequency[drug == "DASATINIB", candidate_pairs], 23)
check_value("Train-test pair overlap", sum(coverage$overlapping_pairs), 0)
check_value("Training drug coverage", sum(coverage$training_drugs == 38L), 5)
check_value("Training cell line coverage", sum(coverage$training_cell_lines == 39L), 5)
check_value("Outcome names in predictor matrices",
  sum(grepl("css|loewe|synergy|observed|outcome", predictors$predictor, ignore.case = TRUE)), 0)

result <- rbindlist(checks)
fwrite(result, "outputs/tables/manuscript_checks.csv")
assert_true(all(result$status == "PASS"), "A manuscript-number check failed; inspect outputs/tables/manuscript_checks.csv")
message(sprintf("Manuscript and predictor checks: %d of %d passed.", sum(result$status == "PASS"), nrow(result)))
