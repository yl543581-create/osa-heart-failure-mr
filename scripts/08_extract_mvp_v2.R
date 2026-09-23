# 08_extract_mvp_v2.R
# FIX: the MVP file's `standard_error` column is entirely "#NA".
# The SE is recoverable from the CI columns:
#     se = (log(ci_upper) - log(ci_lower)) / (2 * 1.96)
# Re-extract significant variants with the CI columns included.

suppressMessages({ library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
gz  <- file.path(data_dir, "GCST90475824.tsv.gz")
out <- file.path(data_dir, "GCST90475824_p5e8.tsv")

cat("streaming with fread...\n")
t0 <- Sys.time()
sig <- fread(
  gz,
  select = c("chromosome","base_pair_location","effect_allele","other_allele",
             "odds_ratio","ci_upper","ci_lower","effect_allele_frequency",
             "p_value","rsid","num_cases","num_controls"),
  sep = "\t", header = TRUE, showProgress = TRUE, nThread = 4
)
cat("read", nrow(sig), "rows in",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

cat("\n--- structure ---\n"); str(sig)
cat("\n--- head ---\n"); print(head(sig, 3))

sig <- sig[!is.na(p_value) & p_value < 5e-8]
cat("\nsignificant (p<5e-8):", nrow(sig), "\n")

cat("\n--- class of ci_upper/ci_lower in significant set ---\n")
cat("ci_upper:", class(sig$ci_upper), " | ci_lower:", class(sig$ci_lower), "\n")
cat("example values:\n"); print(head(sig[, .(rsid, odds_ratio, ci_lower, ci_upper, p_value)], 5))

# numeric coercion (they may have been read as character because of "#NA")
sig[, `:=`(or_n  = suppressWarnings(as.numeric(odds_ratio)),
           cu    = suppressWarnings(as.numeric(ci_upper)),
           cl    = suppressWarnings(as.numeric(ci_lower)),
           eaf_n = suppressWarnings(as.numeric(effect_allele_frequency)),
           p_n   = suppressWarnings(as.numeric(p_value)))]
cat("\nnon-numeric odds_ratio:", sum(is.na(sig$or_n)), "\n")
cat("non-numeric ci_upper  :", sum(is.na(sig$cu)), "\n")
cat("non-numeric ci_lower  :", sum(is.na(sig$cl)), "\n")

# derive beta and se
sig <- sig[!is.na(or_n) & or_n > 0]
sig[, beta := log(or_n)]
sig[, se_ci := (log(cu) - log(cl)) / (2 * 1.96)]
sig <- sig[is.finite(beta) & is.finite(se_ci) & se_ci > 0]
cat("\nusable (beta + se from CI):", nrow(sig), "\n")

# sanity: |beta/se| should approximate the z from p
sig[, z_from_p := qnorm(p_n / 2, lower.tail = FALSE)]
sig[, z_from_ci := abs(beta / se_ci)]
cat("correlation of z(p) vs z(se_ci):",
    round(cor(sig$z_from_p, sig$z_from_ci, use = "complete.obs"), 4), "\n")
cat("median ratio z_ci/z_p:",
    round(median(sig$z_from_ci / sig$z_from_p, na.rm = TRUE), 3), "\n")

keep <- c("rsid","chromosome","base_pair_location","effect_allele","other_allele",
          "eaf_n","or_n","beta","se_ci","p_n","num_cases","num_controls")
res <- sig[, ..keep]
setnames(res, c("rsid","chromosome","bp","effect_allele","other_allele",
                "eaf","or","beta","se","p","num_cases","num_controls"))
setorder(res, p)
fwrite(res, out, sep = "\t")
cat("\nsaved:", out, " rows:", nrow(res), "\n")

# thin to 1 Mb for an overview
thin <- res[, .SD[1], by = .(chromosome, window = floor(bp / 1e6))]
setorder(thin, p)
cat("\n=== top 40 signals (1 Mb thinned) ===\n")
print(thin[1:min(40, .N), .(rsid, chromosome, bp, effect_allele, other_allele,
                            eaf, or, p, num_cases)])
cat("\nn_case:", res$num_cases[1], " n_control:", res$num_controls[1], "\n")
cat("=== done ===\n")
