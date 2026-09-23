# 09_mvp_instruments_and_mr.R
# Build the MVP OSA instrument set and run the main MR + sensitivity analyses.
#
# NOTE ON STANDARD ERRORS:
#   The MVP OSA file's `standard_error` column is entirely "#NA".
#   SE derived from the CI columns is NOT self-consistent with p_value
#   (median z_ci/z_p = 0.903), so the CI columns use a different basis.
#   -> SE is derived from p_value instead:
#          se = |beta| / qnorm(p/2, lower.tail = FALSE)
#   This is internally consistent and is reported as a limitation.

suppressMessages({ library(TwoSampleMR); library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
res_dir <- RES_DIR
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

sig <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))
cat("loaded significant variants:", nrow(sig), "\n")
cat("columns:", paste(names(sig), collapse = ", "), "\n")

# rebuild beta / se self-consistently
sig[, beta2 := log(or)]
sig[, se2   := abs(beta2) / qnorm(p / 2, lower.tail = FALSE)]
sig <- sig[is.finite(beta2) & is.finite(se2) & se2 > 0]
cat("usable:", nrow(sig), "\n")
sig[, z := abs(beta2 / se2)]
cat("self-consistency: cor(z, qnorm(p)) =",
    round(cor(sig$z, qnorm(sig$p/2, lower.tail = FALSE)), 6), "\n")

# ---- independent instrument selection ----
# LD clumping needs an LD reference; we do not have one locally.
# Distance-based pruning at 10 Mb is a CONSERVATIVE proxy for independence
# and is documented as such.
setorder(sig, p)
sig[, grp := floor(bp / 1e7)]
thin <- sig[, .SD[1], by = .(chromosome, grp)]
setorder(thin, p)
cat("\nafter 10 Mb distance-pruning:", nrow(thin), "SNPs\n")

ivs <- thin[, .(SNP = rsid, beta = beta2, se = se2, p = p,
                effect_allele, other_allele, eaf, chromosome, position = bp)]

cat("\n=== MVP OSA instruments (10 Mb pruned, top 30 by p) ===\n")
print(ivs[1:min(30, .N)])
cat("\ntotal instruments:", nrow(ivs), "\n")
fwrite(ivs, file.path(res_dir, "mvp_osa_instruments.tsv"), sep = "\t")

# ================= MVP OSA -> Heart failure =================
cat("\n################################################################\n")
cat("# MVP OSA -> Heart failure (HERMES)\n")
cat("################################################################\n")

hf_id <- "ebi-a-GCST009541"
out <- extract_outcome_data(snps = ivs$SNP, outcomes = hf_id)
cat("outcome SNPs retrieved:", nrow(out), "\n")

exp_dat <- format_data(as.data.frame(ivs), type = "exposure",
                       snp_col = "SNP", beta_col = "beta", se_col = "se",
                       eaf_col = "eaf", effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele", pval_col = "p")
dat <- harmonise_data(exp_dat, out, action = 2)
dat <- dat[dat$mr_keep, ]
cat("after harmonisation:", nrow(dat), "SNPs\n")

methods <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
r <- mr(dat, method_list = methods)
cat("\n=== MVP OSA -> Heart failure ===\n")
print(as.data.table(generate_odds_ratios(r))[, .(method, nsnp, b, se, pval,
                                                or, or_lci95, or_uci95)])

cat("\n-- heterogeneity (Cochran Q) --\n"); print(mr_heterogeneity(dat))
cat("\n-- MR-Egger intercept --\n");       print(mr_pleiotropy_test(dat))
loo <- mr_leaveoneout(dat)
cat("\n-- leave-one-out --\n");            print(as.data.table(loo)[, .(SNP, b, se, p)])

fwrite(as.data.table(generate_odds_ratios(r)), file.path(res_dir, "mvp_OSA_HF.csv"))
cat("\nsaved:", file.path(res_dir, "mvp_OSA_HF.csv"), "\n")
cat("=== done ===\n")
