# 41_download_overweight_hf.R
# Download the two stratified heart-failure endpoints:
#   I9_HEARTFAIL_AND_OVERWEIGHT  (HF + BMI>=25)   <- adiposity-enriched
#   I9_HEARTFAIL_AND_CHD         (HF + CHD)       <- non-adiposity comparator
suppressMessages(library(curl))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
base <- "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/"
targets <- c(
  "finngen_R13_I9_HEARTFAIL_AND_OVERWEIGHT.gz" = "I9_HEARTFAIL_AND_OVERWEIGHT",
  "finngen_R13_I9_HEARTFAIL_AND_CHD.gz"        = "I9_HEARTFAIL_AND_CHD"
)
for (dest_name in names(targets)) {
  pheno <- targets[[dest_name]]
  url <- paste0(base, "finngen_R13_", pheno, ".gz")
  dest <- file.path(data_dir, dest_name)
  cat("\n=== ", pheno, " ===\n", sep="")
  if (file.exists(dest) && file.size(dest) > 1e8) {
    cat("already present:", round(file.size(dest)/1e6,1), "MB\n"); next
  }
  t0 <- Sys.time()
  ok <- tryCatch({ curl::curl_download(url, dest, quiet=TRUE, mode="wb"); TRUE },
                 error=function(e){ cat("FAILED:", conditionMessage(e), "\n"); FALSE })
  if (ok) cat("downloaded", round(file.size(dest)/1e6,1), "MB in",
              round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
}
cat("\n=== data dir ===\n")
for (f in list.files(data_dir, pattern="^finngen_R13", full.names=TRUE))
  cat(sprintf("%-48s %8.1f MB\n", basename(f), file.size(f)/1e6))
