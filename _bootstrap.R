# _bootstrap.R -- load config.R and expose common path aliases.
#
# Rscript does not expose the script's own path via sys.frame() (verified on
# R 4.4.1), so we read the process command line instead. This is the only
# reliable method under `Rscript script.R`.
#
# Usage -- first executable line of any script in scripts/ :
#     source(file.path(dirname(dirname(normalizePath(
#       sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
#     ))), "_bootstrap.R"))
#
# After sourcing this file you have: ROOT, DATA_DIR, RES_DIR, LOG_DIR
# plus the lowercase aliases data_dir / res_dir / log_dir.

local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f) == 0) stop("_bootstrap.R: no --file= in commandArgs; run with Rscript")
  repo_root <- dirname(dirname(normalizePath(sub("^--file=", "", f[1]))))
  cfg <- file.path(repo_root, "config.R")
  if (!file.exists(cfg)) stop(sprintf("_bootstrap.R: config.R not found at %s", repo_root))
  source(cfg, local = FALSE)
})

data_dir <- DATA_DIR
res_dir  <- RES_DIR
log_dir  <- LOG_DIR
