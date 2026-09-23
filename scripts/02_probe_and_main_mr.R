# 02_probe_exposure_file.R
# Two jobs:
#   A) probe the MVP OSA summary-statistics file format using an HTTP Range request
#      (no need to download all 548 MB just to learn the columns)
#   B) run the FIRST REAL MR:  OSA (FinnGen) -> Heart failure (HERMES)
#      with the full battery of sensitivity analyses.

suppressMessages({
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
  library(TwoSampleMR); library(data.table); library(curl)
})

dir.create(file.path(ROOT, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(ROOT, "logs"),    showWarnings = FALSE, recursive = TRUE)

cat("################################################################\n")
cat("# PART A -- probe MVP OSA file format via HTTP Range request\n")
cat("################################################################\n")

url <- paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
              "GCST90475001-GCST90476000/GCST90475824/GCST90475824.tsv.gz")

h <- curl::new_handle()
curl::handle_setopt(h, range = "0-400000", timeout = 120)
raw <- tryCatch(curl::curl_fetch_memory(url, h), error = function(e) NULL)

if (is.null(raw)) {
  cat("Range request FAILED\n")
} else {
  cat("HTTP status :", raw$status_code, "\n")
  cat("bytes       :", length(raw$content), "\n")
  # decompress the partial gzip stream
  con <- gzcon(rawConnection(raw$content, "rb"))
  txt <- tryCatch(readLines(con, n = 8, warn = FALSE), error = function(e)
                    paste("decompress error:", conditionMessage(e)))
  close(con)
  cat("\n--- first lines ---\n")
  cat(paste(txt, collapse = "\n"), "\n")
}

cat("\n################################################################\n")
cat("# PART B -- FIRST REAL MR: OSA -> Heart failure\n")
cat("################################################################\n")

osa_id <- "finn-b-G6_SLEEPAPNO"     # FinnGen OSA (validated above)
hf_id  <- "ebi-a-GCST009541"        # HERMES heart failure (47,309 cases)

exp_dat <- extract_instruments(outcomes = osa_id, p1 = 5e-8,
                               clump = TRUE, r2 = 0.001, kb = 10000)
cat("\nOSA instruments (clumped, r2<0.001, 10Mb):", nrow(exp_dat), "\n")
print(as.data.table(exp_dat)[, .(SNP, effect_allele.exposure, other_allele.exposure,
                                eaf.exposure, beta.exposure, se.exposure,
                                pval.exposure)])

# instrument strength
exp_dat$F <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
cat("\nF statistics:\n"); print(summary(exp_dat$F))
cat("min F =", round(min(exp_dat$F), 1), " | n with F>10:",
    sum(exp_dat$F > 10), "/", nrow(exp_dat), "\n")
# total R2 (approx, for binary exposure use get_r_from_lor-free approximation)
cat("mean F =", round(mean(exp_dat$F), 1), "\n")

out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = hf_id)
cat("\nOutcome SNPs retrieved:", nrow(out_dat), "of", nrow(exp_dat), "\n")

dat <- harmonise_data(exp_dat, out_dat, action = 2)
cat("After harmonisation:", nrow(dat), "SNPs\n")
cat("Palindromic/ambiguous removed:\n")
print(as.data.table(dat)[, .(SNP, mr_keep)])

methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median",
             "mr_weighted_mode")
res <- mr(dat, method_list = methods)
cat("\n=== MR RESULTS: OSA -> Heart failure ===\n")
print(as.data.table(res)[, .(method, nsnp, b, se, pval,
                            or = exp(b), or_lci = exp(b - 1.96*se),
                            or_uci = exp(b + 1.96*se))])

cat("\n--- OR scale (per 1 log-odds genetic liability to OSA) ---\n")
res_or <- generate_odds_ratios(res)
print(as.data.table(res_or)[, .(method, nsnp, or, or_lci95, or_uci95, pval)])

cat("\n=== Sensitivity analyses ===\n")
het <- tryCatch(mr_heterogeneity(dat), error = function(e) paste("FAIL:", conditionMessage(e)))
cat("\n-- Cochran Q / heterogeneity --\n"); print(het)

plei <- tryCatch(mr_pleiotropy_test(dat), error = function(e) paste("FAIL:", conditionMessage(e)))
cat("\n-- MR-Egger intercept --\n"); print(plei)

loo <- tryCatch(mr_leaveoneout(dat), error = function(e) paste("FAIL:", conditionMessage(e)))
if (is.data.frame(loo)) {
  cat("\n-- Leave-one-out --\n")
  print(as.data.table(loo)[, .(SNP, b, se, p)])
}

presso <- tryCatch({
  if (requireNamespace("MRPRESSO", quietly = TRUE)) {
    MRPRESSO::mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
                        SdOutcome = "se.outcome", SdExposure = "se.exposure",
                        OUTCOME = "outcome", EXPOSURE = "exposure",
                        data = as.data.frame(dat), NbDistribution = 3000,
                        SignifThreshold = 0.05)
  } else "MRPRESSO not installed"
}, error = function(e) paste("FAIL:", conditionMessage(e)))
cat("\n-- MR-PRESSO --\n")
if (is.list(presso)) {
  print(presso$`Main MR results`)
  cat("\nGlobal test p-value:", presso$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
} else print(presso)

cat("\n=== Steiger directionality ===\n")
st <- tryCatch(directionality_test(dat), error = function(e) paste("FAIL:", conditionMessage(e)))
print(st)

fwrite(as.data.table(res_or), file.path(ROOT, "results/main_OSA_HF.csv"))
cat("\nSaved: main_OSA_HF.csv\n")
cat("\n=== done ===\n")
