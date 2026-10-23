suppressPackageStartupMessages(library(data.table))
source("R/utils.R")
ensure_directories()

predictions <- readRDS("data/processed/out_of_fold_predictions.rds")
model_columns <- c(
  "pred_training_mean",
  "pred_identity_ridge",
  "pred_feature_elastic_net",
  "pred_hybrid_elastic_net",
  "pred_feature_forest"
)
model_labels <- c(
  pred_training_mean = "Training mean",
  pred_identity_ridge = "Identity ridge",
  pred_feature_elastic_net = "Feature elastic net",
  pred_hybrid_elastic_net = "Hybrid elastic net",
  pred_feature_forest = "Feature random forest"
)

assert_true(all(model_columns %in% names(predictions)), "Prediction table is missing a model")
assert_true(!anyNA(predictions[, ..model_columns]), "Prediction table contains missing values")
assert_true(uniqueN(predictions$triplet_id) == nrow(predictions), "Prediction table has duplicate triplets")
assert_true(uniqueN(predictions$pair_id) == 583L, "Prediction table does not contain 583 pairs")

safe_cor <- function(x, y, method = "pearson") {
  if (sd(x) == 0 || sd(y) == 0) NA_real_ else cor(x, y, method = method)
}

metric_row <- function(data, prediction_column, fold = NA_integer_) {
  observed <- data$observed_loewe
  predicted <- data[[prediction_column]]
  data.table(
    outer_fold = fold,
    model = unname(model_labels[prediction_column]),
    n = nrow(data),
    rmse = rmse(observed, predicted),
    mae = mae(observed, predicted),
    pearson = safe_cor(observed, predicted, "pearson"),
    spearman = safe_cor(observed, predicted, "spearman")
  )
}

overall_performance <- rbindlist(lapply(model_columns, function(model) {
  metric_row(predictions, model)
}))
fold_performance <- rbindlist(lapply(0:4, function(fold) {
  fold_data <- predictions[outer_fold == fold]
  rbindlist(lapply(model_columns, function(model) metric_row(fold_data, model, fold)))
}))
fwrite(overall_performance, "outputs/tables/model_performance_overall.csv")
fwrite(fold_performance, "outputs/tables/model_performance_by_fold.csv")

pair_rows <- split(seq_len(nrow(predictions)), predictions$pair_id)
pair_names <- names(pair_rows)
bootstrap_replicates <- 1000L
set.seed(20260809L)
bootstrap_array <- array(
  NA_real_,
  dim = c(bootstrap_replicates, length(model_columns), 3L),
  dimnames = list(NULL, model_columns, c("rmse", "mae", "pearson"))
)

for (replicate_index in seq_len(bootstrap_replicates)) {
  sampled_pairs <- sample(pair_names, length(pair_names), replace = TRUE)
  sampled_rows <- unlist(pair_rows[sampled_pairs], use.names = FALSE)
  observed <- predictions$observed_loewe[sampled_rows]
  for (model_index in seq_along(model_columns)) {
    predicted <- predictions[[model_columns[model_index]]][sampled_rows]
    bootstrap_array[replicate_index, model_index, "rmse"] <- rmse(observed, predicted)
    bootstrap_array[replicate_index, model_index, "mae"] <- mae(observed, predicted)
    bootstrap_array[replicate_index, model_index, "pearson"] <- safe_cor(observed, predicted)
  }
}

bootstrap_ci <- rbindlist(lapply(seq_along(model_columns), function(model_index) {
  rbindlist(lapply(c("rmse", "mae", "pearson"), function(metric) {
    values <- bootstrap_array[, model_index, metric]
    estimate <- overall_performance[model == model_labels[model_columns[model_index]]][[metric]]
    data.table(
      model = unname(model_labels[model_columns[model_index]]),
      metric = metric,
      estimate = estimate,
      lower_95 = quantile(values, 0.025, na.rm = TRUE, names = FALSE),
      upper_95 = quantile(values, 0.975, na.rm = TRUE, names = FALSE),
      bootstrap_replicates = bootstrap_replicates
    )
  }))
}))
fwrite(bootstrap_ci, "outputs/tables/model_performance_pair_bootstrap_ci.csv")

identity_index <- match("pred_identity_ridge", model_columns)
rmse_difference <- rbindlist(lapply(seq_along(model_columns), function(model_index) {
  if (model_index == identity_index) return(NULL)
  differences <- bootstrap_array[, model_index, "rmse"] - bootstrap_array[, identity_index, "rmse"]
  data.table(
    comparison = paste0(model_labels[model_columns[model_index]], " minus Identity ridge"),
    estimate = overall_performance[model == model_labels[model_columns[model_index]], rmse] -
      overall_performance[model == "Identity ridge", rmse],
    lower_95 = quantile(differences, 0.025, names = FALSE),
    upper_95 = quantile(differences, 0.975, names = FALSE),
    probability_difference_below_zero = mean(differences < 0)
  )
}))
fwrite(rmse_difference, "outputs/tables/rmse_differences_vs_identity_ridge.csv")

synergy_threshold <- 10
top_fraction <- 0.10
ranking_models <- setdiff(model_columns, "pred_training_mean")

ranking_by_fold <- rbindlist(lapply(0:4, function(fold) {
  fold_data <- predictions[outer_fold == fold]
  baseline_positive_rate <- mean(fold_data$observed_loewe > synergy_threshold)
  rbindlist(lapply(ranking_models, function(model) {
    top_n <- ceiling(top_fraction * nrow(fold_data))
    selected <- fold_data[order(-get(model))][seq_len(top_n)]
    positive_rate <- mean(selected$observed_loewe > synergy_threshold)
    data.table(
      outer_fold = fold,
      model = unname(model_labels[model]),
      n_candidates = top_n,
      mean_observed_loewe = mean(selected$observed_loewe),
      median_observed_loewe = median(selected$observed_loewe),
      positive_rate = positive_rate,
      fold_positive_rate = baseline_positive_rate,
      positive_enrichment = positive_rate / baseline_positive_rate,
      mean_css = mean(selected$css),
      median_css = median(selected$css)
    )
  }))
}))
fwrite(ranking_by_fold, "outputs/tables/triplet_top_decile_by_fold.csv")

ranking_overall <- ranking_by_fold[, .(
  folds = .N,
  candidates = sum(n_candidates),
  mean_observed_loewe = weighted.mean(mean_observed_loewe, n_candidates),
  median_fold_observed_loewe = median(median_observed_loewe),
  positive_rate = weighted.mean(positive_rate, n_candidates),
  baseline_positive_rate = weighted.mean(fold_positive_rate, n_candidates),
  positive_enrichment = weighted.mean(positive_enrichment, n_candidates),
  mean_css = weighted.mean(mean_css, n_candidates),
  median_fold_css = median(median_css)
), by = model]
fwrite(ranking_overall, "outputs/tables/triplet_top_decile_overall.csv")

make_pair_rankings <- function(model_column) {
  pair_results <- predictions[, {
    predicted <- get(model_column)
    selected_index <- order(predicted, decreasing = TRUE)[seq_len(3L)]
    best_index <- selected_index[1L]
    .(
      pair_prediction = mean(predicted[selected_index]),
      selected_observed_loewe = mean(observed_loewe[selected_index]),
      selected_css = mean(css[selected_index]),
      best_cell_line = cell_line[best_index],
      best_cell_prediction = predicted[best_index],
      best_cell_observed_loewe = observed_loewe[best_index],
      best_cell_css = css[best_index],
      selected_cell_lines = paste(cell_line[selected_index], collapse = "; "),
      selected_cell_observed_loewe = paste(sprintf("%.2f", observed_loewe[selected_index]), collapse = "; "),
      selected_cell_css = paste(sprintf("%.2f", css[selected_index]), collapse = "; ")
    )
  }, by = .(pair_id, drug_1, drug_2, outer_fold)]
  pair_results[, model := unname(model_labels[model_column])]
  pair_results[, fold_rank := frank(-pair_prediction, ties.method = "first"), by = outer_fold]
  pair_results[]
}

pair_rankings <- rbindlist(lapply(ranking_models, make_pair_rankings), use.names = TRUE)
setcolorder(pair_rankings, c("model", setdiff(names(pair_rankings), "model")))
fwrite(pair_rankings, "outputs/tables/pair_candidate_rankings.csv")

pair_top_decile <- pair_rankings[, {
  selected <- .SD[fold_rank <= ceiling(top_fraction * .N)]
  .(
    candidate_pairs = nrow(selected),
    mean_confirmation_loewe = mean(selected$selected_observed_loewe),
    median_confirmation_loewe = median(selected$selected_observed_loewe),
    confirmation_rate_above_10 = mean(selected$selected_observed_loewe > synergy_threshold),
    mean_confirmation_css = mean(selected$selected_css),
    median_confirmation_css = median(selected$selected_css)
  )
}, by = .(model, outer_fold)]
fwrite(pair_top_decile, "outputs/tables/pair_top_decile_by_fold.csv")

selected_model <- overall_performance[which.min(rmse), model]
assert_true(selected_model == "Identity ridge", "The expected primary model was not selected")
candidate_shortlist <- pair_rankings[
  model == selected_model & fold_rank <= 2L
][order(outer_fold, fold_rank)]
fwrite(candidate_shortlist, "outputs/tables/candidate_shortlist_identity_ridge.csv")

identity_pair_top_decile <- pair_rankings[
  model == selected_model,
  .SD[fold_rank <= ceiling(top_fraction * .N)],
  by = outer_fold
]
chapter_pair_points <- identity_pair_top_decile[, .(
  outer_fold, fold_rank, drug_1, drug_2,
  pair_prediction, selected_observed_loewe, selected_css,
  observed_minus_predicted = selected_observed_loewe - pair_prediction
)]
fwrite(chapter_pair_points, "outputs/tables/chapter_pair_points.csv")
pair_calibration_summary <- chapter_pair_points[, .(
  selected_pairs = .N,
  observed_above_predicted = sum(observed_minus_predicted > 0),
  mean_observed_minus_predicted = mean(observed_minus_predicted)
)]
fwrite(pair_calibration_summary, "outputs/tables/pair_calibration_summary.csv")
pair_candidate_drug_frequency <- data.table(
  drug = c(identity_pair_top_decile$drug_1, identity_pair_top_decile$drug_2)
)[, .(candidate_pairs = .N), by = drug][order(-candidate_pairs, drug)]
pair_candidate_drug_frequency[, proportion_of_60_pairs := candidate_pairs / nrow(identity_pair_top_decile)]
fwrite(pair_candidate_drug_frequency, "outputs/tables/pair_candidate_drug_frequency.csv")

identity_triplet_top_decile <- rbindlist(lapply(0:4, function(fold) {
  fold_data <- predictions[outer_fold == fold][order(-pred_identity_ridge)]
  fold_data[seq_len(ceiling(top_fraction * nrow(fold_data)))]
}))
triplet_candidate_drug_frequency <- data.table(
  drug = c(identity_triplet_top_decile$drug_1, identity_triplet_top_decile$drug_2)
)[, .(candidate_triplets = .N), by = drug][order(-candidate_triplets, drug)]
triplet_candidate_drug_frequency[, proportion_of_selected_drug_slots :=
  candidate_triplets / (2L * nrow(identity_triplet_top_decile))]
fwrite(triplet_candidate_drug_frequency, "outputs/tables/triplet_candidate_drug_frequency.csv")

triplet_candidate_cell_frequency <- identity_triplet_top_decile[
  , .(candidate_triplets = .N), by = cell_line
][order(-candidate_triplets, cell_line)]
triplet_candidate_cell_frequency[, proportion_of_selected_triplets :=
  candidate_triplets / nrow(identity_triplet_top_decile)]
fwrite(triplet_candidate_cell_frequency, "outputs/tables/triplet_candidate_cell_line_frequency.csv")

central_lower <- quantile(predictions$observed_loewe, 0.01)
central_upper <- quantile(predictions$observed_loewe, 0.99)
robustness_groups <- list(
  "All triplets" = rep(TRUE, nrow(predictions)),
  "Single-measurement triplets" = predictions$replicate_count == 1L,
  "Repeated triplets after averaging" = predictions$replicate_count == 2L,
  "Central 98% of observed Loewe scores" =
    predictions$observed_loewe >= central_lower & predictions$observed_loewe <= central_upper
)
robustness_metrics <- rbindlist(lapply(names(robustness_groups), function(group_name) {
  subset_data <- predictions[robustness_groups[[group_name]]]
  rows <- rbindlist(lapply(model_columns, function(model) metric_row(subset_data, model)))
  rows[, group := group_name]
  rows
}))
setcolorder(robustness_metrics, c("group", setdiff(names(robustness_metrics), "group")))
fwrite(robustness_metrics, "outputs/tables/model_robustness_metrics.csv")

endpoint_sensitivity <- data.table(
  comparison = c(
    "Published DeepSynergy Loewe versus DrugComb recalculated Loewe",
    "Published DeepSynergy Loewe versus CSS",
    "DrugComb recalculated Loewe versus CSS"
  ),
  pearson = c(
    cor(readRDS("data/processed/oneil_analysis_outcomes.rds")$loewe,
      readRDS("data/processed/oneil_analysis_outcomes.rds")$drugcomb_loewe),
    cor(readRDS("data/processed/oneil_analysis_outcomes.rds")$loewe,
      readRDS("data/processed/oneil_analysis_outcomes.rds")$css),
    cor(readRDS("data/processed/oneil_analysis_outcomes.rds")$drugcomb_loewe,
      readRDS("data/processed/oneil_analysis_outcomes.rds")$css)
  ),
  spearman = c(
    cor(readRDS("data/processed/oneil_analysis_outcomes.rds")$loewe,
      readRDS("data/processed/oneil_analysis_outcomes.rds")$drugcomb_loewe, method = "spearman"),
    cor(readRDS("data/processed/oneil_analysis_outcomes.rds")$loewe,
      readRDS("data/processed/oneil_analysis_outcomes.rds")$css, method = "spearman"),
    cor(readRDS("data/processed/oneil_analysis_outcomes.rds")$drugcomb_loewe,
      readRDS("data/processed/oneil_analysis_outcomes.rds")$css, method = "spearman")
  )
)
fwrite(endpoint_sensitivity, "outputs/tables/endpoint_relationships.csv")

summary_lines <- c(
  sprintf("Selected model: %s.", selected_model),
  sprintf("Identity ridge RMSE %.3f (pair-bootstrap 95%% CI %.3f to %.3f).",
    overall_performance[model == selected_model, rmse],
    bootstrap_ci[model == selected_model & metric == "rmse", lower_95],
    bootstrap_ci[model == selected_model & metric == "rmse", upper_95]),
  sprintf("Identity ridge Pearson correlation %.3f and Spearman correlation %.3f.",
    overall_performance[model == selected_model, pearson],
    overall_performance[model == selected_model, spearman]),
  sprintf("Identity ridge top-decile triplet positive rate %.1f%% versus %.1f%% overall.",
    100 * ranking_overall[model == selected_model, positive_rate],
    100 * ranking_overall[model == selected_model, baseline_positive_rate]),
  "Candidate shortlist contains the top two held-out pairs from each outer fold.",
  "Observed Loewe and CSS values were revealed only after ranking the out-of-fold predictions."
)
writeLines(summary_lines, "logs/04_evaluate_models.log")
cat(paste(summary_lines, collapse = "\n"), "\n")
