# 18_diagnose_finngen.R
suppressMessages(library(curl))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
u <- "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_I9_HEARTFAIL.gz"

cat("=== test 1: range request (first 2 MB) ===\n")
h <- curl::new_handle()
curl::handle_setopt(h, range = "0-2000000", timeout = 120, connecttimeout = 30)
r <- tryCatch(curl::curl_fetch_memory(u, h), error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
if (!is.null(r)) cat("status:", r$status_code, " bytes:", length(r$content), "\n")

cat("\n=== test 2: timed download of 5 MB via curl_download to temp ===\n")
tf <- tempfile()
t0 <- Sys.time()
res <- tryCatch({
  h2 <- curl::new_handle()
  curl::handle_setopt(h2, range = "0-5000000", timeout = 120, connecttimeout = 30,
                      followlocation = TRUE, useragent = "Mozilla/5.0")
  curl::curl_download(u, tf, handle = h2, quiet = TRUE)
  "ok"
}, error=function(e) paste("FAIL:", conditionMessage(e)))
el <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat("result:", res, "\n")
if (file.exists(tf)) cat("bytes:", file.size(tf), " in ", round(el,1), "s  =>",
                          round(file.size(tf)/1e6/max(el,0.1), 2), "MB/s\n")

cat("\n=== test 3: d4file / alternative mirrors ===\n")
alts <- c(
  "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv",
  "https://finngen-public-data-r13.storage.googleapis.com/summary_stats/finngen_R13_manifest.tsv"
)
for (a in alts) {
  rr <- tryCatch({ hh<-curl::new_handle(); curl::handle_setopt(hh, timeout=60, nobody=TRUE)
                   curl::curl_fetch_memory(a, hh) }, error=function(e) NULL)
  cat(a, "->", ifelse(is.null(rr), "FAIL", rr$status_code), "\n")
}
