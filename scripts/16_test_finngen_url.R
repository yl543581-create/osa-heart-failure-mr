# 16_test_finngen_url.R
suppressMessages(library(curl))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
urls <- c(
  "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_G6_SLEEPAPNO_INCLAVO.gz",
  "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_I9_HEARTFAIL.gz"
)
for (u in urls) {
  cat("\n===", basename(u), "===\n")
  h <- curl::new_handle()
  curl::handle_setopt(h, nobody = TRUE, timeout = 60, connecttimeout = 30)
  r <- tryCatch(curl::curl_fetch_memory(u, h), error = function(e) { cat("FAIL:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) {
    cat("status:", r$status_code, "\n")
    cat("headers of interest:\n")
    for (nm in intersect(names(r$headers), c("content-length","content-type","location","HTTP/2 200"))) {
      cat("  ", nm, ":", r$headers[[nm]], "\n", sep="")
    }
  }
}
cat("\n=== alternative: R13 manifest (small) ===\n")
h2 <- curl::new_handle(); curl::handle_setopt(h2, timeout=60)
r2 <- tryCatch(curl::curl_fetch_memory(
  "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv", h2),
  error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
if (!is.null(r2)) cat("manifest status:", r2$status_code, " bytes:", length(r2$content), "\n")
