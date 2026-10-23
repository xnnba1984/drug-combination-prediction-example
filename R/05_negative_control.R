suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(glmnet)
})
source("R/utils.R")
ensure_directories()

analysis_data <- readRDS("data/processed/oneil_analysis_outcomes.rds")
observed_predictions <- readRDS("data/processed/out_of_fold_predictions.rds")

drugs <- sort(unique(c(analysis_data$drug_1, analysis_data$drug_2)))
cell_lines <- sort(unique(analysis_data$cell_line))
drug_index <- setNames(seq_along(drugs), drugs)
cell_index <- setNames(seq_along(cell_lines), cell_lines)
n <- nrow(analysis_data)
n_drugs <- length(drugs)
n_cells <- length(cell_lines)
row_index <- seq_len(n)
drug_1_index <- unname(drug_index[analysis_data$drug_1])
drug_2_index <- unname(drug_index[analysis_data$drug_2])
cell_row_index <- unname(cell_index[analysis_data$cell_line])

identity_matrix <- sparseMatrix(
  i = rep(row_index, 5L),
  j = c(
    drug_1_index,
    drug_2_index,
    n_drugs + cell_row_index,
    n_drugs + n_cells + (drug_1_index - 1L) * n_cells + cell_row_index,
    n_drugs + n_cells + (drug_2_index - 1L) * n_cells + cell_row_index
  ),
  x = 1,
  dims = c(n, n_drugs + n_cells + n_drugs * n_cells)
)

pair_cell_counts <- analysis_data[, .N, by = pair_id]
assert_true(all(pair_cell_counts$N == 39L), "Pair-block permutation requires 39 rows per pair")

fit_oof_ridge <- function(outcome) {
  prediction <- rep(NA_real_, n)
  for (test_fold in 0:4) {
    train_index <- which(analysis_data$outer_fold != test_fold)
    test_index <- which(analysis_data$outer_fold == test_fold)
    inner_folds <- match(
      analysis_data$outer_fold[train_index],
      sort(unique(analysis_data$outer_fold[train_index]))
    )
    fit <- cv.glmnet(
      x = identity_matrix[train_index, , drop = FALSE],
      y = outcome[train_index],
      family = "gaussian",
      alpha = 0,
      foldid = inner_folds,
      type.measure = "mse",
      standardize = TRUE,
      nlambda = 60
    )
    prediction[test_index] <- as.numeric(predict(
      fit,
      identity_matrix[test_index, , drop = FALSE],
      s = "lambda.min"
    ))
  }
  prediction
}

pair_ids <- sort(unique(analysis_data$pair_id))
outcome_lookup <- analysis_data[, .(pair_id, cell_line, loewe)]
setkey(outcome_lookup, pair_id, cell_line)
permutation_count <- 20L
set.seed(20260809L)
permutation_results <- vector("list", permutation_count)

for (permutation_index in seq_len(permutation_count)) {
  donor_pairs <- sample(pair_ids, length(pair_ids), replace = FALSE)
  donor_map <- setNames(donor_pairs, pair_ids)
  permuted_outcome <- outcome_lookup[
    data.table(
      pair_id = unname(donor_map[analysis_data$pair_id]),
      cell_line = analysis_data$cell_line
    ),
    on = .(pair_id, cell_line),
    x.loewe
  ]
  assert_true(length(permuted_outcome) == n && all(is.finite(permuted_outcome)),
    "Pair-block permutation failed")
  permuted_prediction <- fit_oof_ridge(permuted_outcome)
  permutation_results[[permutation_index]] <- data.table(
    permutation = permutation_index,
    rmse = rmse(permuted_outcome, permuted_prediction),
    mae = mae(permuted_outcome, permuted_prediction),
    pearson = cor(permuted_outcome, permuted_prediction),
    spearman = cor(permuted_outcome, permuted_prediction, method = "spearman")
  )
  message(sprintf(
    "Negative-control permutation %s/%s: RMSE %.3f; Pearson %.3f",
    permutation_index,
    permutation_count,
    permutation_results[[permutation_index]]$rmse,
    permutation_results[[permutation_index]]$pearson
  ))
}

permutation_table <- rbindlist(permutation_results)
fwrite(permutation_table, "outputs/tables/identity_ridge_pair_permutation_control.csv")

actual <- data.table(
  rmse = rmse(observed_predictions$observed_loewe, observed_predictions$pred_identity_ridge),
  mae = mae(observed_predictions$observed_loewe, observed_predictions$pred_identity_ridge),
  pearson = cor(observed_predictions$observed_loewe, observed_predictions$pred_identity_ridge),
  spearman = cor(
    observed_predictions$observed_loewe,
    observed_predictions$pred_identity_ridge,
    method = "spearman"
  )
)
control_summary <- data.table(
  metric = c("RMSE", "MAE", "Pearson", "Spearman"),
  actual = unlist(actual[1L]),
  permutation_mean = c(
    mean(permutation_table$rmse),
    mean(permutation_table$mae),
    mean(permutation_table$pearson),
    mean(permutation_table$spearman)
  ),
  permutation_min = c(
    min(permutation_table$rmse),
    min(permutation_table$mae),
    min(permutation_table$pearson),
    min(permutation_table$spearman)
  ),
  permutation_max = c(
    max(permutation_table$rmse),
    max(permutation_table$mae),
    max(permutation_table$pearson),
    max(permutation_table$spearman)
  ),
  permutations = permutation_count
)
fwrite(control_summary, "outputs/tables/identity_ridge_negative_control_summary.csv")

summary_lines <- c(
  sprintf("Completed %s pair-block outcome permutations.", permutation_count),
  sprintf("Actual identity-ridge RMSE %.3f versus permutation range %.3f to %.3f.",
    actual$rmse, min(permutation_table$rmse), max(permutation_table$rmse)),
  sprintf("Actual identity-ridge Pearson correlation %.3f versus permutation range %.3f to %.3f.",
    actual$pearson, min(permutation_table$pearson), max(permutation_table$pearson)),
  "Each permutation retained the donor pair's 39-cell-line outcome profile while breaking its link to drug identity."
)
writeLines(summary_lines, "logs/05_negative_control.log")
cat(paste(summary_lines, collapse = "\n"), "\n")
