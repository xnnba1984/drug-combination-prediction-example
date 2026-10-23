suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(glmnet)
  library(ranger)
})
source("R/utils.R")
ensure_directories()
threads <- as.integer(Sys.getenv("SECTION3_THREADS", "8"))
assert_true(length(threads) == 1L && is.finite(threads) && threads >= 1L, "SECTION3_THREADS must be a positive integer")

analysis_data <- readRDS("data/processed/oneil_analysis_outcomes.rds")
feature_objects <- readRDS("data/processed/model_features.rds")
assert_true(identical(analysis_data$triplet_id, feature_objects$triplet_id),
  "Outcome and feature rows are not aligned")

y <- analysis_data$loewe
outer_fold <- analysis_data$outer_fold
base_features <- feature_objects$base_feature_matrix
linear_features <- feature_objects$linear_feature_matrix

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

identity_i <- c(
  row_index,
  row_index,
  row_index,
  row_index,
  row_index
)
identity_j <- c(
  drug_1_index,
  drug_2_index,
  n_drugs + cell_row_index,
  n_drugs + n_cells + (drug_1_index - 1L) * n_cells + cell_row_index,
  n_drugs + n_cells + (drug_2_index - 1L) * n_cells + cell_row_index
)
identity_matrix <- sparseMatrix(
  i = identity_i,
  j = identity_j,
  x = 1,
  dims = c(n, n_drugs + n_cells + n_drugs * n_cells)
)
colnames(identity_matrix) <- c(
  paste0("drug_present_", drugs),
  paste0("cell_", cell_lines),
  unlist(lapply(drugs, function(drug) paste0("drug_cell_", drug, "__", cell_lines)))
)
assert_true(ncol(identity_matrix) == 1559L, "Unexpected identity feature count")
assert_true(all(Matrix::rowSums(identity_matrix) == 5), "Identity matrix rows have unexpected nonzero counts")

hybrid_matrix <- cbind(identity_matrix, Matrix(linear_features, sparse = TRUE))
colnames(hybrid_matrix) <- c(colnames(identity_matrix), paste0("external_", colnames(linear_features)))
assert_true(nrow(hybrid_matrix) == n && ncol(hybrid_matrix) == 2263L,
  "Unexpected hybrid feature dimensions")

predictor_columns <- rbindlist(list(
  data.table(design = "identity", predictor = colnames(identity_matrix)),
  data.table(design = "base", predictor = colnames(base_features)),
  data.table(design = "linear", predictor = colnames(linear_features)),
  data.table(design = "hybrid", predictor = colnames(hybrid_matrix))
))
assert_true(!any(grepl("css|loewe|synergy|observed|outcome", predictor_columns$predictor,
  ignore.case = TRUE)), "An outcome column entered the predictor matrices")
fwrite(predictor_columns, "outputs/tables/model_predictor_columns.csv")

recode_inner_folds <- function(folds) {
  match(folds, sort(unique(folds)))
}

fit_identity_ridge <- function(train_index, test_index) {
  inner_fold_id <- recode_inner_folds(outer_fold[train_index])
  fit <- cv.glmnet(
    x = identity_matrix[train_index, , drop = FALSE],
    y = y[train_index],
    family = "gaussian",
    alpha = 0,
    foldid = inner_fold_id,
    type.measure = "mse",
    standardize = TRUE,
    nlambda = 80,
    keep = FALSE
  )
  list(
    prediction = as.numeric(predict(fit, identity_matrix[test_index, , drop = FALSE], s = "lambda.min")),
    lambda = fit$lambda.min,
    inner_rmse = sqrt(min(fit$cvm))
  )
}

fit_feature_elastic_net <- function(train_index, test_index) {
  inner_fold_id <- recode_inner_folds(outer_fold[train_index])
  alpha_candidates <- c(0, 0.5, 1)
  candidate_fits <- lapply(alpha_candidates, function(alpha_value) {
    cv.glmnet(
      x = linear_features[train_index, , drop = FALSE],
      y = y[train_index],
      family = "gaussian",
      alpha = alpha_value,
      foldid = inner_fold_id,
      type.measure = "mse",
      standardize = TRUE,
      nlambda = 80,
      keep = FALSE
    )
  })
  candidate_rmse <- vapply(candidate_fits, function(fit) sqrt(min(fit$cvm)), numeric(1L))
  selected <- which.min(candidate_rmse)
  selected_fit <- candidate_fits[[selected]]
  list(
    prediction = as.numeric(predict(
      selected_fit,
      linear_features[test_index, , drop = FALSE],
      s = "lambda.min"
    )),
    alpha = alpha_candidates[selected],
    lambda = selected_fit$lambda.min,
    inner_rmse = candidate_rmse[selected],
    alpha_candidates = alpha_candidates,
    candidate_rmse = candidate_rmse
  )
}

fit_hybrid_elastic_net <- function(train_index, test_index) {
  inner_fold_id <- recode_inner_folds(outer_fold[train_index])
  alpha_candidates <- c(0, 0.5, 1)
  candidate_fits <- lapply(alpha_candidates, function(alpha_value) {
    cv.glmnet(
      x = hybrid_matrix[train_index, , drop = FALSE],
      y = y[train_index],
      family = "gaussian",
      alpha = alpha_value,
      foldid = inner_fold_id,
      type.measure = "mse",
      standardize = TRUE,
      nlambda = 80,
      keep = FALSE
    )
  })
  candidate_rmse <- vapply(candidate_fits, function(fit) sqrt(min(fit$cvm)), numeric(1L))
  selected <- which.min(candidate_rmse)
  selected_fit <- candidate_fits[[selected]]
  list(
    prediction = as.numeric(predict(
      selected_fit,
      hybrid_matrix[test_index, , drop = FALSE],
      s = "lambda.min"
    )),
    alpha = alpha_candidates[selected],
    lambda = selected_fit$lambda.min,
    inner_rmse = candidate_rmse[selected],
    alpha_candidates = alpha_candidates,
    candidate_rmse = candidate_rmse
  )
}

rf_grid <- unique(data.table(
  mtry = rep(c(round(sqrt(ncol(base_features))), round(ncol(base_features) / 3)), each = 2L),
  min_node_size = rep(c(5L, 20L), times = 2L)
))

fit_feature_forest <- function(train_index, test_index, test_fold) {
  inner_folds <- sort(unique(outer_fold[train_index]))
  tuning_rows <- vector("list", nrow(rf_grid))

  for (grid_index in seq_len(nrow(rf_grid))) {
    squared_error <- 0
    validation_count <- 0L
    for (inner_fold in inner_folds) {
      inner_validation <- train_index[outer_fold[train_index] == inner_fold]
      inner_training <- train_index[outer_fold[train_index] != inner_fold]
      fit <- ranger(
        x = as.data.frame(base_features[inner_training, , drop = FALSE]),
        y = y[inner_training],
        num.trees = 250,
        mtry = rf_grid$mtry[grid_index],
        min.node.size = rf_grid$min_node_size[grid_index],
        seed = 20260809L + 1000L * test_fold + 100L * grid_index + inner_fold,
        num.threads = threads,
        write.forest = TRUE
      )
      prediction <- predict(
        fit,
        data = as.data.frame(base_features[inner_validation, , drop = FALSE]),
        num.threads = threads
      )$predictions
      squared_error <- squared_error + sum((y[inner_validation] - prediction)^2)
      validation_count <- validation_count + length(inner_validation)
    }
    tuning_rows[[grid_index]] <- data.table(
      test_fold = test_fold,
      mtry = rf_grid$mtry[grid_index],
      min_node_size = rf_grid$min_node_size[grid_index],
      inner_rmse = sqrt(squared_error / validation_count)
    )
  }

  tuning <- rbindlist(tuning_rows)
  selected <- tuning[which.min(inner_rmse)]
  final_fit <- ranger(
    x = as.data.frame(base_features[train_index, , drop = FALSE]),
    y = y[train_index],
    num.trees = 750,
    mtry = selected$mtry,
    min.node.size = selected$min_node_size,
    seed = 20260809L + 1000L * test_fold + 999L,
    num.threads = threads,
    write.forest = TRUE
  )
  prediction <- predict(
    final_fit,
    data = as.data.frame(base_features[test_index, , drop = FALSE]),
    num.threads = threads
  )$predictions
  list(prediction = prediction, selected = selected, tuning = tuning)
}

prediction_table <- analysis_data[, .(
  triplet_id,
  pair_id,
  drug_1,
  drug_2,
  cell_line,
  outer_fold,
  observed_loewe = loewe,
  css,
  replicate_count,
  loewe_range
)]
prediction_table[, `:=`(
  pred_training_mean = NA_real_,
  pred_identity_ridge = NA_real_,
  pred_feature_elastic_net = NA_real_,
  pred_hybrid_elastic_net = NA_real_,
  pred_feature_forest = NA_real_
)]

tuning_results <- list()
model_log <- character()

for (test_fold in 0:4) {
  test_index <- which(outer_fold == test_fold)
  train_index <- which(outer_fold != test_fold)
  prediction_table$pred_training_mean[test_index] <- mean(y[train_index])

  message(sprintf("Outer fold %s: fitting identity ridge", test_fold))
  ridge_result <- fit_identity_ridge(train_index, test_index)
  prediction_table$pred_identity_ridge[test_index] <- ridge_result$prediction
  tuning_results[[length(tuning_results) + 1L]] <- data.table(
    test_fold = test_fold,
    model = "Identity ridge",
    alpha = 0,
    lambda = ridge_result$lambda,
    mtry = NA_integer_,
    min_node_size = NA_integer_,
    inner_rmse = ridge_result$inner_rmse
  )

  message(sprintf("Outer fold %s: fitting feature elastic net", test_fold))
  elastic_result <- fit_feature_elastic_net(train_index, test_index)
  prediction_table$pred_feature_elastic_net[test_index] <- elastic_result$prediction
  tuning_results[[length(tuning_results) + 1L]] <- data.table(
    test_fold = test_fold,
    model = "Feature elastic net",
    alpha = elastic_result$alpha,
    lambda = elastic_result$lambda,
    mtry = NA_integer_,
    min_node_size = NA_integer_,
    inner_rmse = elastic_result$inner_rmse
  )
  tuning_results[[length(tuning_results) + 1L]] <- data.table(
    test_fold = test_fold,
    model = paste0("Feature elastic net alpha candidate ", elastic_result$alpha_candidates),
    alpha = elastic_result$alpha_candidates,
    lambda = NA_real_,
    mtry = NA_integer_,
    min_node_size = NA_integer_,
    inner_rmse = elastic_result$candidate_rmse
  )

  message(sprintf("Outer fold %s: fitting hybrid identity-plus-feature elastic net", test_fold))
  hybrid_result <- fit_hybrid_elastic_net(train_index, test_index)
  prediction_table$pred_hybrid_elastic_net[test_index] <- hybrid_result$prediction
  tuning_results[[length(tuning_results) + 1L]] <- data.table(
    test_fold = test_fold,
    model = "Hybrid elastic net",
    alpha = hybrid_result$alpha,
    lambda = hybrid_result$lambda,
    mtry = NA_integer_,
    min_node_size = NA_integer_,
    inner_rmse = hybrid_result$inner_rmse
  )
  tuning_results[[length(tuning_results) + 1L]] <- data.table(
    test_fold = test_fold,
    model = paste0("Hybrid elastic net alpha candidate ", hybrid_result$alpha_candidates),
    alpha = hybrid_result$alpha_candidates,
    lambda = NA_real_,
    mtry = NA_integer_,
    min_node_size = NA_integer_,
    inner_rmse = hybrid_result$candidate_rmse
  )

  message(sprintf("Outer fold %s: tuning and fitting feature random forest", test_fold))
  forest_result <- fit_feature_forest(train_index, test_index, test_fold)
  prediction_table$pred_feature_forest[test_index] <- forest_result$prediction
  tuning_results[[length(tuning_results) + 1L]] <- forest_result$tuning[, .(
    test_fold,
    model = "Feature random forest candidate",
    alpha = NA_real_,
    lambda = NA_real_,
    mtry,
    min_node_size,
    inner_rmse
  )]
  tuning_results[[length(tuning_results) + 1L]] <- forest_result$selected[, .(
    test_fold,
    model = "Feature random forest selected",
    alpha = NA_real_,
    lambda = NA_real_,
    mtry,
    min_node_size,
    inner_rmse
  )]

  fold_checkpoint <- list(
    test_fold = test_fold,
    predictions = prediction_table[outer_fold == test_fold],
    tuning = rbindlist(tuning_results, fill = TRUE)
  )
  saveRDS(fold_checkpoint, sprintf("data/interim/model_checkpoint_fold_%s.rds", test_fold), compress = "xz")

  fold_message <- sprintf(
    "Completed outer fold %s: %s training rows, %s test rows; ridge RMSE %.3f; feature elastic-net RMSE %.3f; hybrid RMSE %.3f; forest RMSE %.3f.",
    test_fold,
    length(train_index),
    length(test_index),
    rmse(y[test_index], ridge_result$prediction),
    rmse(y[test_index], elastic_result$prediction),
    rmse(y[test_index], hybrid_result$prediction),
    rmse(y[test_index], forest_result$prediction)
  )
  message(fold_message)
  model_log <- c(model_log, fold_message)
  writeLines(model_log, "logs/03_fit_models.log")
}

prediction_columns <- grep("^pred_", names(prediction_table), value = TRUE)
assert_true(!anyNA(prediction_table[, ..prediction_columns]), "At least one out-of-fold prediction is missing")
assert_true(all(is.finite(as.matrix(prediction_table[, ..prediction_columns]))),
  "At least one out-of-fold prediction is non-finite")

tuning_table <- rbindlist(tuning_results, fill = TRUE)
fwrite(prediction_table, "data/processed/out_of_fold_predictions.csv")
saveRDS(prediction_table, "data/processed/out_of_fold_predictions.rds", compress = "xz")
fwrite(tuning_table, "outputs/tables/model_tuning_results.csv")
writeLines(c(model_log, capture.output(sessionInfo())), "logs/03_fit_models_session_info.txt")

cat(paste(model_log, collapse = "\n"), "\n")
