# 07_extract_mvp_fast.R
# Fast extraction of genome-wide significant variants from the MVP OSA GWAS.
# Strategy: decompress in C (R's gzfile is slow for line-by-line), read big
# chunks with data.table::fread on a pipe-free connection, filter, discard.

suppressMessages({ library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
gz  <- file.path(data_dir, "GCST90475824.tsv.gz")
out <- file.path(data_dir, "GCST90475824_p5e8.tsv")
cat("input:", gz, " size:", file.size(gz), "bytes\n")

# --- Pass 1: decompress to a temporary plain file? 548MB gz -> ~7GB plain, too big.
# --- Better: fread supports reading directly from a .gz file, streaming.
cat("\nstreaming with fread (gz-aware)...\n")
t0 <- Sys.time()

# Read only the columns we need; fread handles gzip natively.
sig <- fread(
  gz,
  select = c("chromosome","base_pair_location","effect_allele","other_allele",
             "odds_ratio","standard_error","effect_allele_frequency","p_value",
             "rsid","num_cases","num_controls"),
  sep = "\t", header = TRUE, showProgress = TRUE, nThread = 4
)
cat("read rows:", nrow(sig), " in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

cat("\ncolumn check:\n"); print(head(sig, 3))
cat("\np_value class:", class(sig$p_value), "\n")
cat("NA in p_value:", sum(is.na(sig$p_value)), "\n")

sig <- sig[!is.na(p_value) & p_value < 5e-8]
cat("significant (p<5e-8):", nrow(sig), "\n")

sig[, beta := log(odds_ratio)]
sig <- sig[is.finite(beta) & !is.na(standard_error) & standard_error > 0]
cat("usable:", nrow(sig), "\n")

fwrite(sig, out, sep = "\t")
cat("saved:", out, "\n")

# thin to 1 Mb windows for a first look
setorder(sig, p_value)
thin <- sig[, .SD[1], by = .(chromosome, window = floor(base_pair_location / 1e6))]
setorder(thin, p_value)
cat("\n=== top 40 signals (1 Mb thinned) ===\n")
print(thin[1:min(40, .N), .(rsid, chromosome, base_pair_location,
                            effect_allele, other_allele,
                            eaf = effect_allele_frequency,
                            or = odds_ratio, p_value)])
cat("\nn_case:", sig$num_cases[1], " n_control:", sig$num_controls[1], "\n")
cat("=== done ===\n")
