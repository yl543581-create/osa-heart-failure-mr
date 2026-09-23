# 04_download_and_extract_mvp_osa.R
# Stream-download the MVP OSA GWAS (548 MB) and keep only genome-wide
# significant variants (p < 5e-8) -> small, fast, reusable instrument set.

suppressMessages({ library(data.table); library(curl) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
              "GCST90475001-GCST90476000/GCST90475824/GCST90475824.tsv.gz")
gz  <- file.path(data_dir, "GCST90475824.tsv.gz")

cat("=== downloading MVP OSA (GCST90475824) ===\n")
cat("target:", gz, "\n")

if (!file.exists(gz) || file.size(gz) < 5e8) {
  h <- curl::new_handle()
  curl::handle_setopt(h, timeout = 3600L, connecttimeout = 60L,
                      followlocation = TRUE, failonerror = FALSE)
  t0 <- Sys.time()
  curl::curl_download(url, gz, handle = h, quiet = FALSE)
  cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
} else {
  cat("already present, size =", file.size(gz), "\n")
}
cat("final size (bytes):", file.size(gz), "\n")

cat("\n=== streaming parse, keeping p < 5e-8 ===\n")
# read in chunks; select only needed columns to save memory
cols <- c("chromosome","base_pair_location","effect_allele","other_allele",
          "odds_ratio","standard_error","effect_allele_frequency","p_value",
          "rsid","num_cases","num_controls")

keep <- list()
chunk_size <- 2e6
n_read <- 0L
con <- gzfile(gz, "rt")
header <- strsplit(readLines(con, n = 1), "\t")[[1]]
cat("header:", paste(header, collapse = " | "), "\n")
idx <- match(cols, header)
cat("column indices:", paste(idx, collapse = ","), "\n")

repeat {
  chunk <- tryCatch(
    fread(text = readLines(con, n = chunk_size), sep = "\t",
          header = FALSE, select = idx, col.names = cols,
          na.strings = c("#NA","NA",""), showProgress = FALSE),
    error = function(e) NULL)
  if (is.null(chunk) || nrow(chunk) == 0) break
  n_read <- n_read + nrow(chunk)
  sub <- chunk[!is.na(p_value) & p_value < 5e-8]
  if (nrow(sub) > 0) keep[[length(keep) + 1]] <- sub
  if (n_read %% 2e7 < chunk_size) cat("  rows read:", n_read,
                                      " | sig so far:", sum(sapply(keep, nrow)), "\n")
}
close(con)

cat("total rows read:", n_read, "\n")
sig <- rbindlist(keep)
cat("genome-wide significant variants (p<5e-8):", nrow(sig), "\n")

sig[, beta := log(odds_ratio)]
sig <- sig[is.finite(beta) & !is.na(standard_error) & standard_error > 0]
cat("usable after beta/se filter:", nrow(sig), "\n")

# save ALL significant variants (for clumping) 
fwrite(sig, file.path(data_dir, "GCST90475824_p5e8.tsv"), sep = "\t")

# quick look: top independent-ish signals by p, thinned to 1 Mb windows
setorder(sig, p_value)
thin <- sig[, .SD[1], by = .(chromosome, window = floor(base_pair_location / 1e6))]
setorder(thin, p_value)
cat("\n=== top signals (1 Mb thinned), first 30 ===\n")
print(thin[1:min(30, .N), .(rsid, chromosome, base_pair_location,
                            effect_allele, other_allele, eaf = effect_allele_frequency,
                            or = odds_ratio, p_value, num_cases)])

cat("\n=== summary ===\n")
cat("MVP OSA: n_case =", sig$num_cases[1], " n_control =", sig$num_controls[1], "\n")
cat("Saved:", file.path(data_dir, "GCST90475824_p5e8.tsv"), "\n")
cat("=== done ===\n")
