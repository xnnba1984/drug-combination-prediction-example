assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

normalize_identifier <- function(x) {
  toupper(gsub("[^A-Za-z0-9]", "", x))
}

canonical_pair_id <- function(drug_1, drug_2) {
  paste(pmin(drug_1, drug_2), pmax(drug_1, drug_2), sep = "||")
}

triplet_id <- function(drug_1, drug_2, cell_line) {
  paste(canonical_pair_id(drug_1, drug_2), cell_line, sep = "||")
}

finite_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

finite_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

rmse <- function(observed, predicted) {
  sqrt(mean((observed - predicted)^2))
}

mae <- function(observed, predicted) {
  mean(abs(observed - predicted))
}

ensure_directories <- function() {
  directories <- c(
    "data/raw",
    "data/interim",
    "data/processed",
    "outputs/tables",
    "outputs/figures",
    "logs"
  )
  invisible(lapply(directories, dir.create, recursive = TRUE, showWarnings = FALSE))
}
