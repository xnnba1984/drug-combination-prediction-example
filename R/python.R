run_python <- function(script, ...) {
  default <- if (.Platform$OS.type == "windows") ".venv/Scripts/python.exe" else ".venv/bin/python"
  python <- Sys.getenv("SECTION3_PYTHON", default)
  if (!file.exists(python)) python <- Sys.which(python)
  if (!nzchar(python)) stop("Create .venv first or set SECTION3_PYTHON to a Python executable.", call. = FALSE)
  status <- system2(python, c(shQuote(script), vapply(list(...), shQuote, character(1L))))
  if (status != 0L) stop(sprintf("Python step failed: %s", script), call. = FALSE)
}
