# 17_download_finngen_v2.R
# Restart the FinnGen download with resume + progress logging.
suppressMessages(library(curl))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR

targets <- c(
  "finnngen_R13_G6_SLEEPAPNO_INCLAVO.gz" = "G6_SLEEPAPNO_INCLAVO",
  "finngen_R13_I9_HEARTFAIL.gz"          = "I9_HEARTFAIL"
)
base <- "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"

# clean up the stalled temp files
for (f in list.files(data_dir, pattern="curltmp$", full.names=TRUE)) {
  cat("removing stalled temp:", f, "\n"); file.remove(f)
}

for (dest_name in names(targets)) {
  pheno <- targets[[dest_name]]
  url <- paste0(base, "finngen_R13_", pheno, ".gz")
  dest <- file.path(data_dir, dest_name)
  cat("\n=== ", pheno, " ===\n", sep="")
  if (file.exists(dest) && file.size(dest) > 1e8) {
    cat("already complete:", round(file.size(dest)/1e6,1), "MB\n"); next
  }
  t0 <- Sys.time()
  ok <- tryCatch({
    curl::curl_download(url, dest, quiet = FALSE, mode = "wb")
    TRUE
  }, error = function(e) { cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("downloaded", round(file.size(dest)/1e6,1), "MB in",
              round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
}

cat("\n=== final data dir ===\n")
for (f in list.files(data_dir, full.names=TRUE)) {
  cat(sprintf("%-45s %8.1f MB\n", basename(f), file.size(f)/1e6))
}
