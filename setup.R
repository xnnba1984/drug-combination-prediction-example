local_library <- file.path(getwd(), ".Rlib")
dir.create(local_library, showWarnings = FALSE)
.libPaths(c(local_library, .Library))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
Sys.setenv(RENV_CONFIG_CACHE_ENABLED = "FALSE")

# Source builds on some macOS installations need the existing gettext headers.
if (Sys.info()[["sysname"]] == "Darwin") {
  prefixes <- c(Sys.getenv("HOMEBREW_PREFIX"), "/opt/homebrew", "/usr/local")
  gettext <- file.path(prefixes[nzchar(prefixes)], "opt", "gettext")
  gettext <- gettext[file.exists(file.path(gettext, "include", "libintl.h"))]
  if (length(gettext)) {
    prepend_path <- function(value, existing) {
      paste(c(value, existing[nzchar(existing)]), collapse = .Platform$path.sep)
    }
    Sys.setenv(
      CPATH = prepend_path(file.path(gettext[1L], "include"), Sys.getenv("CPATH")),
      LIBRARY_PATH = prepend_path(file.path(gettext[1L], "lib"), Sys.getenv("LIBRARY_PATH"))
    )
    message("Using installed gettext headers and libraries for macOS source builds.")
  }
}

renv_version <- "1.2.4"
installed <- installed.packages(lib.loc = local_library)
if (!"renv" %in% rownames(installed) || installed["renv", "Version"] != renv_version) {
  archive <- tempfile(fileext = ".tar.gz")
  urls <- c(
    sprintf("https://cloud.r-project.org/src/contrib/renv_%s.tar.gz", renv_version),
    sprintf("https://cloud.r-project.org/src/contrib/Archive/renv/renv_%s.tar.gz", renv_version)
  )
  downloaded <- FALSE
  for (url in urls) {
    status <- suppressWarnings(try(download.file(url, archive, mode = "wb"), silent = TRUE))
    if (!inherits(status, "try-error") && status == 0L) {
      downloaded <- TRUE
      break
    }
  }
  if (!downloaded) stop("Could not download the locked renv version.", call. = FALSE)
  install.packages(archive, repos = NULL, type = "source", lib = local_library)
  unlink(archive)
}
renv::restore(project = getwd(), library = local_library, lockfile = "renv.lock", prompt = FALSE)
source("R/environment.R")
check_environment()
message("R environment restored and checked.")
