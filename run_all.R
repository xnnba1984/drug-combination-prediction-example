if (!file.exists("renv.lock")) stop("Run this script from the repository root.", call. = FALSE)
source("R/environment.R")
source("R/python.R")
check_environment()
run_python("scripts/verify_inputs.py")

analysis_scripts <- sprintf("R/%s", c(
  "01_prepare_data.R", "02_prepare_features.R", "03_fit_models.R",
  "04_evaluate_models.R", "05_negative_control.R", "06_make_figures.R",
  "07_final_qc.R", "08_manuscript_checks.R"
))
dir.create("logs", showWarnings = FALSE)
started <- Sys.time()
run_log <- character()
for (script in analysis_scripts) {
  message(sprintf("\nRunning %s", script))
  step_started <- Sys.time()
  source(script, local = new.env(parent = globalenv()))
  elapsed <- as.numeric(difftime(Sys.time(), step_started, units = "secs"))
  run_log <- c(run_log, sprintf("Completed %s in %.1f seconds", script, elapsed))
  writeLines(run_log, "logs/run_all.log")
}
run_python("scripts/clean_pdf_metadata.py", "outputs/figures")
write_environment_record("logs/environment.txt")
run_log <- c(run_log, sprintf("Total elapsed time: %.1f seconds",
  as.numeric(difftime(Sys.time(), started, units = "secs"))))
writeLines(run_log, "logs/run_all.log")
cat(paste(run_log, collapse = "\n"), "\n")
