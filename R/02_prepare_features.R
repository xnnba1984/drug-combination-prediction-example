suppressPackageStartupMessages(library(data.table))
source("R/utils.R")
ensure_directories()

outcome_path <- "data/processed/oneil_analysis_outcomes.rds"
drug_feature_path <- "data/interim/drug_features.csv"
cell_feature_path <- "data/interim/cell_features.csv"
for (path in c(outcome_path, drug_feature_path, cell_feature_path)) {
  assert_true(file.exists(path), paste("Required file is missing:", path))
}

analysis_data <- readRDS(outcome_path)
drug_features <- fread(drug_feature_path)
cell_features <- fread(cell_feature_path)

assert_true(nrow(drug_features) == 38L && ncol(drug_features) == 4388L,
  "Unexpected drug feature dimensions")
assert_true(nrow(cell_features) == 39L && ncol(cell_features) == 3985L,
  "Unexpected cell-line feature dimensions")
assert_true(setequal(drug_features$drug_name, unique(c(analysis_data$drug_1, analysis_data$drug_2))),
  "Drug features do not cover the analysis drugs")
assert_true(setequal(cell_features$cell_line, unique(analysis_data$cell_line)),
  "Cell-line features do not cover the analysis cell lines")

drug_matrix <- as.matrix(drug_features[, -"drug_name"])
rownames(drug_matrix) <- drug_features$drug_name
cell_matrix <- as.matrix(cell_features[, -"cell_line"])
rownames(cell_matrix) <- cell_features$cell_line

assert_true(all(is.finite(drug_matrix)), "Drug feature matrix contains non-finite values")
assert_true(all(is.finite(cell_matrix)), "Cell-line feature matrix contains non-finite values")

drug_sds <- apply(drug_matrix, 2L, sd)
cell_sds <- apply(cell_matrix, 2L, sd)
drug_variable <- is.finite(drug_sds) & drug_sds > 0
cell_variable <- is.finite(cell_sds) & cell_sds > 0

drug_pca <- prcomp(drug_matrix[, drug_variable, drop = FALSE], center = TRUE, scale. = TRUE)
cell_pca <- prcomp(cell_matrix[, cell_variable, drop = FALSE], center = TRUE, scale. = TRUE)

drug_cumulative_variance <- cumsum(drug_pca$sdev^2) / sum(drug_pca$sdev^2)
cell_cumulative_variance <- cumsum(cell_pca$sdev^2) / sum(cell_pca$sdev^2)
variance_threshold <- 0.80
drug_pc_count <- which(drug_cumulative_variance >= variance_threshold)[1L]
cell_pc_count <- which(cell_cumulative_variance >= variance_threshold)[1L]

drug_scores <- drug_pca$x[, seq_len(drug_pc_count), drop = FALSE]
cell_scores <- cell_pca$x[, seq_len(cell_pc_count), drop = FALSE]
colnames(drug_scores) <- sprintf("drug_pc_%02d", seq_len(drug_pc_count))
colnames(cell_scores) <- sprintf("cell_pc_%02d", seq_len(cell_pc_count))

drug_1_scores <- drug_scores[analysis_data$drug_1, , drop = FALSE]
drug_2_scores <- drug_scores[analysis_data$drug_2, , drop = FALSE]
cell_observation_scores <- cell_scores[analysis_data$cell_line, , drop = FALSE]

drug_mean <- (drug_1_scores + drug_2_scores) / 2
drug_absdiff <- abs(drug_1_scores - drug_2_scores)
colnames(drug_mean) <- sprintf("drug_mean_pc_%02d", seq_len(drug_pc_count))
colnames(drug_absdiff) <- sprintf("drug_absdiff_pc_%02d", seq_len(drug_pc_count))

base_feature_matrix <- cbind(drug_mean, drug_absdiff, cell_observation_scores)
storage.mode(base_feature_matrix) <- "double"
assert_true(nrow(base_feature_matrix) == nrow(analysis_data), "Feature rows do not match outcome rows")
assert_true(all(is.finite(base_feature_matrix)), "Base feature matrix contains non-finite values")

swapped_mean <- (drug_2_scores + drug_1_scores) / 2
swapped_absdiff <- abs(drug_2_scores - drug_1_scores)
assert_true(max(abs(drug_mean - swapped_mean)) == 0, "Drug-mean features are not order invariant")
assert_true(max(abs(drug_absdiff - swapped_absdiff)) == 0, "Drug-difference features are not order invariant")

pair_feature_count <- 2L * drug_pc_count
cell_start <- pair_feature_count + 1L
pair_features <- base_feature_matrix[, seq_len(pair_feature_count), drop = FALSE]
cell_features_for_rows <- base_feature_matrix[, cell_start:ncol(base_feature_matrix), drop = FALSE]

interaction_columns <- vector("list", ncol(cell_features_for_rows))
for (cell_index in seq_len(ncol(cell_features_for_rows))) {
  block <- pair_features * cell_features_for_rows[, cell_index]
  colnames(block) <- paste0(
    colnames(pair_features),
    "_x_",
    colnames(cell_features_for_rows)[cell_index]
  )
  interaction_columns[[cell_index]] <- block
}
linear_feature_matrix <- cbind(base_feature_matrix, do.call(cbind, interaction_columns))
storage.mode(linear_feature_matrix) <- "double"
assert_true(all(is.finite(linear_feature_matrix)), "Linear-model feature matrix contains non-finite values")

feature_objects <- list(
  triplet_id = analysis_data$triplet_id,
  base_feature_matrix = base_feature_matrix,
  linear_feature_matrix = linear_feature_matrix,
  drug_scores = drug_scores,
  cell_scores = cell_scores,
  variance_threshold = variance_threshold,
  drug_pc_count = drug_pc_count,
  cell_pc_count = cell_pc_count,
  drug_pca = drug_pca,
  cell_pca = cell_pca,
  drug_variable = drug_variable,
  cell_variable = cell_variable
)
saveRDS(feature_objects, "data/processed/model_features.rds", compress = "xz")

fwrite(
  data.table(drug_name = rownames(drug_scores), drug_scores),
  "data/processed/drug_pca_scores.csv"
)
fwrite(
  data.table(cell_line = rownames(cell_scores), cell_scores),
  "data/processed/cell_line_pca_scores.csv"
)

pca_summary <- rbindlist(list(
  data.table(
    entity = "Drug",
    original_features = ncol(drug_matrix),
    nonconstant_features = sum(drug_variable),
    retained_components = drug_pc_count,
    retained_variance = drug_cumulative_variance[drug_pc_count]
  ),
  data.table(
    entity = "Cell line",
    original_features = ncol(cell_matrix),
    nonconstant_features = sum(cell_variable),
    retained_components = cell_pc_count,
    retained_variance = cell_cumulative_variance[cell_pc_count]
  )
))
fwrite(pca_summary, "outputs/tables/pca_summary.csv")

drug_feature_groups <- data.table(
  feature_group = c("ECFP6 fingerprint", "Physicochemical properties", "Toxicophore indicators"),
  first_column = c(1L, 1310L, 2112L),
  last_column = c(1309L, 2111L, 4387L)
)
drug_feature_groups[, original_features := last_column - first_column + 1L]
drug_feature_groups[, nonconstant_features := vapply(
  seq_len(.N),
  function(index) sum(drug_variable[first_column[index]:last_column[index]]),
  integer(1L)
)]
fwrite(drug_feature_groups, "outputs/tables/drug_feature_group_audit.csv")

fold_coverage <- rbindlist(lapply(0:4, function(test_fold) {
  training <- analysis_data[outer_fold != test_fold]
  testing <- analysis_data[outer_fold == test_fold]
  data.table(
    test_fold = test_fold,
    training_pairs = uniqueN(training$pair_id),
    test_pairs = uniqueN(testing$pair_id),
    overlapping_pairs = length(intersect(unique(training$pair_id), unique(testing$pair_id))),
    training_drugs = uniqueN(c(training$drug_1, training$drug_2)),
    test_drugs = uniqueN(c(testing$drug_1, testing$drug_2)),
    training_cell_lines = uniqueN(training$cell_line),
    test_cell_lines = uniqueN(testing$cell_line)
  )
}))
assert_true(all(fold_coverage$overlapping_pairs == 0L), "Outer folds share a drug pair")
assert_true(all(fold_coverage$training_drugs == 38L), "An outer training fold lacks a drug")
assert_true(all(fold_coverage$training_cell_lines == 39L), "An outer training fold lacks a cell line")
fwrite(fold_coverage, "outputs/tables/outer_fold_coverage.csv")

log_lines <- c(
  sprintf("Removed %s constant drug features and %s constant cell-line features.", sum(!drug_variable), sum(!cell_variable)),
  sprintf("Retained %s drug PCs (%.1f%% variance) and %s cell-line PCs (%.1f%% variance).",
    drug_pc_count, 100 * drug_cumulative_variance[drug_pc_count],
    cell_pc_count, 100 * cell_cumulative_variance[cell_pc_count]),
  sprintf("Created %s base features and %s linear-model features for %s triplets.",
    ncol(base_feature_matrix), ncol(linear_feature_matrix), nrow(base_feature_matrix)),
  "Verified exact feature invariance to Drug 1-Drug 2 exchange.",
  "Verified that every outer training fold contains all 38 drugs and all 39 cell lines."
)
writeLines(log_lines, "logs/02_prepare_features.log")
cat(paste(log_lines, collapse = "\n"), "\n")
