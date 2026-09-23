# config.R -- portable path resolution for the OSA-HF MR project.
#
# Works whether a script is run with:
#     Rscript scripts/14_definitive_analysis.R      (from repo root)
#     Rscript /abs/path/to/repo/scripts/xx.R        (from anywhere)
# or sourced interactively.
#
# Resolution order:
#   1. environment variable OSA_HF_ROOT, if set and valid
#   2. derived from the running script's own path (strip /scripts/, /R/)
#   3. walk up from the working directory looking for a `data` or `results` dir
#   4. the working directory itself
#
# Defines: ROOT, DATA_DIR, RES_DIR, LOG_DIR   (and creates the latter three)

.osahf_script_dir <- function() {
  # --file=... is how Rscript records the script path
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grepl("^--file=", a)]
  if (length(f)) {
    return(dirname(normalizePath(sub("^--file=", "", f[1]))))
  }
  # sourced interactively: try the ofile of the calling frame
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of) && is.character(of) && nzchar(of)) return(dirname(normalizePath(of)))
  NA_character_
}

.osahf_root <- function() {
  # 1) explicit override
  env <- Sys.getenv("OSA_HF_ROOT", "")
  if (nzchar(env) && dir.exists(env)) return(normalizePath(env))

  # 2) from the script's own location
  sd <- .osahf_script_dir()
  if (!is.na(sd)) {
    # if the script lives in scripts/ or R/, the repo root is one level up
    cand <- if (basename(sd) %in% c("scripts", "R", "code")) dirname(sd) else sd
    for (p in unique(c(cand, sd))) {
      if (dir.exists(file.path(p, "data")) || dir.exists(file.path(p, "results")))
        return(normalizePath(p))
    }
    return(normalizePath(cand))
  }

  # 3) walk up from the working directory
  wd <- normalizePath(getwd())
  for (i in 1:6) {
    if (dir.exists(file.path(wd, "data")) || dir.exists(file.path(wd, "results")))
      return(wd)
    parent <- dirname(wd)
    if (identical(parent, wd)) break
    wd <- parent
  }

  # 4) last resort
  normalizePath(getwd())
}

ROOT     <- .osahf_root()
DATA_DIR <- file.path(ROOT, "data")
RES_DIR  <- file.path(ROOT, "results")
LOG_DIR  <- file.path(ROOT, "logs")

for (d in c(DATA_DIR, RES_DIR, LOG_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
