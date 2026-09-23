# 33_get_manifest.R
# Download the FinnGen R13 manifest to find usable endpoints for a PheWAS
# and to document the exact OSA endpoint definition.
suppressMessages(library(curl))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
dest <- file.path(data_dir, "finngen_R13_manifest.tsv")
if (!file.exists(dest) || file.size(dest) < 1e5) {
  curl::curl_download(
    "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv",
    dest, quiet = TRUE)
}
cat("manifest size:", file.size(dest), "bytes\n")
m <- read.delim(dest, sep="\t", stringsAsFactors=FALSE)
cat("rows:", nrow(m), " cols:", paste(names(m), collapse=", "), "\n\n")
# the endpoints we care about
key <- m[grepl("SLEEPAPNO|HEARTFAIL", m$phenocode, ignore.case=TRUE), ]
cat("=== OSA / heart failure endpoints ===\n")
print(key[, intersect(c("phenocode","phenostring","n_case","n_control","num_gw_signif"), names(m))])
cat("\n=== sample of columns for one row ===\n")
if (nrow(key) > 0) { str(key[1, ]) }
