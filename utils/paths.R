# Project directories. Run scripts from the repository root; set SPG_ROOT to keep data and results elsewhere.
REPO      <- normalizePath(".")
ROOT      <- Sys.getenv("SPG_ROOT", REPO)
DATA      <- file.path(ROOT, "data")
RAW       <- file.path(DATA, "raw")
METADATA  <- file.path(DATA, "metadata")
PROCESSED <- file.path(DATA, "processed")
EXTERNAL  <- file.path(DATA, "external")
RESULTS   <- file.path(ROOT, "results")
FIGURES   <- file.path(RESULTS, "figures")
VENDOR    <- file.path(ROOT, "vendor")
RESOURCES <- file.path(REPO, "resources")

out_dir <- function(...) {
  d <- file.path(RESULTS, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
