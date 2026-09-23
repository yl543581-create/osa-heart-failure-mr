# 13_clump_proper.R
# Replace the 10 Mb distance-pruning proxy with REAL LD-based clumping (r2<0.001, 10 Mb).
# Strategy: try the OpenGWAS LD clumping API first (no local reference needed).
# Fall back to a local 1000G EUR PLINK reference if the API is unavailable.

suppressMessages({ library(data.table); library(ieugwasr); library(TwoSampleMR) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR

sig <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))
cat("significant variants loaded:", nrow(sig), "\n")
cat("columns:", paste(names(sig), collapse=", "), "\n")

# rebuild self-consistent beta/se
sig[, beta2 := log(or)]
sig[, se2   := abs(beta2) / qnorm(p/2, lower.tail = FALSE)]
sig <- sig[is.finite(beta2) & is.finite(se2) & se2 > 0]
cat("usable:", nrow(sig), "\n")
cat("p range:", range(sig$p), "\n\n")

# ---------- attempt API clumping ----------
cat("=== attempting OpenGWAS API clumping (r2<0.001, kb=10000, EUR) ===\n")
clumped <- tryCatch({
  ieugwasr::ld_clump(
    dat    = data.frame(rsid = sig$rsid, pval = sig$p),
    clump_r2 = 0.001,
    clump_kb = 10000,
    pop    = "EUR"
  )
}, error = function(e) { cat("API clump FAILED:", conditionMessage(e), "\n"); NULL })

if (!is.null(clumped) && nrow(clumped) > 0) {
  cat("API clumping SUCCEEDED ->", nrow(clumped), "independent SNPs\n")
  ivs <- sig[rsid %in% clumped$rsid]
  cat("matched back:", nrow(ivs), "\n")
  used_api <- TRUE
} else {
  cat("\nAPI unavailable. Downloading local 1000G EUR reference instead...\n")
  used_api <- FALSE
  ivs <- NULL
}

if (!used_api) {
  # MRC IEU 1000G EUR reference (PLINK format)
  ref_dir <- file.path(ROOT, "ldref")
  dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)
  base <- "https://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz"
  tgz  <- file.path(ref_dir, "1kg.v3.tgz")
  if (!file.exists(tgz) || file.size(tgz) < 1e6) {
    cat("downloading", base, "...\n")
    try(utils::download.file(base, tgz, mode = "wb", quiet = FALSE))
  }
  cat("ldref archive size:", file.size(tgz), "bytes\n")
}

cat("\n=== done stage 1 ===\n")
cat("used_api:", used_api, "\n")
if (!is.null(ivs)) {
  fwrite(ivs, file.path(res_dir, "mvp_osa_instruments_LDclumped.tsv"), sep = "\t")
  cat("saved LD-clumped instruments:", nrow(ivs), "\n")
}
