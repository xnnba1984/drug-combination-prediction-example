.libPaths(c(file.path(getwd(), ".Rlib"), .Library))

check_environment <- function() {
  if (!requireNamespace("renv", quietly = TRUE)) {
    stop("Run Rscript --vanilla setup.R before the analysis.", call. = FALSE)
  }
  lock <- renv::lockfile_read("renv.lock")
  wanted <- vapply(lock$Packages, function(x) x$Version, character(1L))
  actual <- vapply(names(wanted), function(package) {
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    packageDescription(package, fields = "Version")
  }, character(1L))
  wrong <- names(wanted)[is.na(actual) | actual != wanted]
  if (length(wrong)) stop(paste("Locked package versions are missing:", paste(wrong, collapse = ", ")), call. = FALSE)
  if (as.character(getRversion()) != lock$R$Version) {
    warning(sprintf("The reference analysis used R %s; this session uses R %s.", lock$R$Version, getRversion()))
  }
  invisible(actual)
}

write_environment_record <- function(path) {
  versions <- check_environment()
  record <- c(
    paste("R", as.character(getRversion())),
    paste("Platform", R.version$platform),
    paste(names(versions), versions)
  )
  writeLines(record, path)
}
