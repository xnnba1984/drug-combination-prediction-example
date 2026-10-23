suppressPackageStartupMessages(library(data.table))
source("R/utils.R")
ensure_directories()

required_files <- c(
  "data/processed/oneil_analysis_outcomes.rds",
  "data/processed/model_features.rds",
  "data/processed/out_of_fold_predictions.rds",
  "outputs/tables/model_performance_overall.csv",
  "outputs/tables/model_performance_pair_bootstrap_ci.csv",
  "outputs/tables/candidate_shortlist_identity_ridge.csv",
  "outputs/tables/identity_ridge_pair_permutation_control.csv",
  "outputs/figures/worked_example_results.png",
  "outputs/figures/worked_example_results.pdf"
)
assert_true(all(file.exists(required_files)), "One or more required outputs are missing")
assert_true(all(file.info(required_files)$size > 0), "One or more required outputs are empty")

outcomes <- readRDS("data/processed/oneil_analysis_outcomes.rds")
features <- readRDS("data/processed/model_features.rds")
predictions <- readRDS("data/processed/out_of_fold_predictions.rds")
performance <- fread("outputs/tables/model_performance_overall.csv")
bootstrap_ci <- fread("outputs/tables/model_performance_pair_bootstrap_ci.csv")
shortlist <- fread("outputs/tables/candidate_shortlist_identity_ridge.csv")
permutation_control <- fread("outputs/tables/identity_ridge_pair_permutation_control.csv")
pair_rankings <- fread("outputs/tables/pair_candidate_rankings.csv")

model_columns <- c(
  "pred_training_mean",
  "pred_identity_ridge",
  "pred_feature_elastic_net",
  "pred_hybrid_elastic_net",
  "pred_feature_forest"
)

checks <- data.table(
  check = c(
    "Unique analysis triplets",
    "Unique drug pairs",
    "Cell lines",
    "Outer folds",
    "One fold per pair",
    "Complete 39-cell-line pair profiles",
    "Feature rows match outcomes",
    "Prediction rows match outcomes",
    "Unique out-of-fold predictions",
    "Complete predictions for all models",
    "Five performance rows",
    "Bootstrap intervals use 1000 pair samples",
    "Ten shortlisted pairs",
    "Two shortlisted pairs per fold",
    "Twenty pair-block permutations",
    "Sixty top-decile identity-ridge pairs"
  ),
  observed = c(
    as.character(uniqueN(outcomes$triplet_id)),
    as.character(uniqueN(outcomes$pair_id)),
    as.character(uniqueN(outcomes$cell_line)),
    paste(sort(unique(outcomes$outer_fold)), collapse = ","),
    as.character(outcomes[, uniqueN(outer_fold), by = pair_id][, max(V1)]),
    paste(sort(unique(outcomes[, .N, by = pair_id]$N)), collapse = ","),
    as.character(nrow(features$base_feature_matrix)),
    as.character(nrow(predictions)),
    as.character(uniqueN(predictions$triplet_id)),
    as.character(sum(is.na(predictions[, ..model_columns]))),
    as.character(nrow(performance)),
    paste(sort(unique(bootstrap_ci$bootstrap_replicates)), collapse = ","),
    as.character(nrow(shortlist)),
    paste(shortlist[, .N, by = outer_fold][order(outer_fold), N], collapse = ","),
    as.character(nrow(permutation_control)),
    as.character(pair_rankings[
      model == "Identity ridge",
      .SD[fold_rank <= ceiling(0.10 * .N)],
      by = outer_fold
    ][, .N])
  ),
  expected = c(
    "22737", "583", "39", "0,1,2,3,4", "1", "39", "22737", "22737",
    "22737", "0", "5", "1000", "10", "2,2,2,2,2", "20", "60"
  )
)
checks[, status := fifelse(observed == expected, "PASS", "FAIL")]
fwrite(checks, "outputs/tables/final_qc_checks.csv")

assert_true(all(checks$status == "PASS"), "One or more final quality-control checks failed")
assert_true(
  identical(sort(outcomes$triplet_id), sort(predictions$triplet_id)),
  "Outcome and prediction triplet identifiers do not match"
)
assert_true(
  identical(outcomes$triplet_id, features$triplet_id),
  "Outcome and feature triplet identifiers do not match"
)
assert_true(
  performance[which.min(rmse), model] == "Identity ridge",
  "Identity ridge is no longer the model with the lowest cross-validated RMSE"
)

log_lines <- c(
  sprintf("Final quality-control checks passed: %d of %d.", nrow(checks), nrow(checks)),
  sprintf("Validated %s unique triplets and %s drug pairs.", checks[1, observed], checks[2, observed]),
  "Out-of-fold predictions are complete and pair-disjoint.",
  "The candidate shortlist and pair-block permutation control have the frozen sizes.",
  "The primary PNG and vector PDF figures are present and nonempty."
)
writeLines(log_lines, "logs/07_final_qc.log")
cat(paste(log_lines, collapse = "\n"), "\n")
