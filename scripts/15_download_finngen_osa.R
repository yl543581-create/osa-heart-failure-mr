# 15_download_finngen_osa.R
# Download FinnGen R13 OSA + heart failure summary statistics.
# Purpose: (a) meta-analyse with MVP to strengthen the OSA instrument
#          (b) independent replication cohort for OSA -> HF

suppressMessages({ library(curl) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

base <- "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
targets <- c(
  G6_SLEEPAPNO_INCLAVO = "finnngen_R13_G6_SLEEPAPNO_INCLAVO.gz",
  I9_HEARTFAIL         = "finngen_R13_I9_HEARTFAIL.gz"
)

dl <- function(pheno, fname) {
  url <- paste0(base, "finngen_R13_", pheno, ".gz")
  dest <- file.path(data_dir, fname)
  cat("\n=== ", pheno, " ===\n", sep = "")
  cat("url:", url, "\n")
  if (file.exists(dest) && file.size(dest) > 1e8) {
    cat("already present:", round(file.size(dest)/1e6, 1), "MB\n"); return(invisible(dest))
  }
  h <- curl::new_handle()
  curl::handle_setopt(h, timeout = 3600L, connecttimeout = 60L, followlocation = TRUE)
  t0 <- Sys.time()
  tryCatch({
    curl::curl_download(url, dest, handle = h, quiet = TRUE)
    cat("downloaded", round(file.size(dest)/1e6, 1), "MB in",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
  }, error = function(e) cat("DOWNLOAD FAILED:", conditionMessage(e), "\n"))
  invisible(dest)
}

for (p in names(targets)) dl(p, targets[[p]])

# peek at the FinnGen format
f <- file.path(data_dir, "finngen_R13_G6_SLEEPAPNO_INCLAVO.gz")
if (file.exists(f) && file.size(f) > 1e6) {
  cat("\n=== FinnGen OSA format (first rows) ===\n")
  con <- gzfile(f, "rt"); hdr <- readLines(con, n = 3); close(con)
  cat(paste(hdr, collapse = "\n"), "\n")
}
cat("\n=== done ===\n")
