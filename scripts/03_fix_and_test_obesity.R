# 03_fix_and_test_obesity.R
# (a) fix MR-PRESSO + Steiger
# (b) test the decisive question: IS THE OSA INSTRUMENT JUST PROXYING BMI?
#     - what is the OSA instrument's causal effect ON BMI?
#     - does a pure BMI instrument reproduce the same HF effect?

suppressMessages({
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
  library(TwoSampleMR); library(data.table); library(MRPRESSO)
})
options(warn = 1)

osa_id <- "finn-b-G6_SLEEPAPNO"
hf_id  <- "ebi-a-GCST009541"
bmi_id <- "ieu-b-40"          # GIANT + UKB BMI

cat("################################################################\n")
cat("# (a) OSA -> Heart failure, with MR-PRESSO fixed\n")
cat("################################################################\n")

exp_dat <- extract_instruments(outcomes = osa_id, p1 = 5e-8, clump = TRUE,
                               r2 = 0.001, kb = 10000)
out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = hf_id)
dat <- harmonise_data(exp_dat, out_dat, action = 2)
dat <- dat[dat$mr_keep, ]

pr <- MRPRESSO::mr_presso(
  BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
  SdOutcome   = "se.outcome",   SdExposure   = "se.exposure",
  data = as.data.frame(dat), NbDistribution = 5000, SignifThreshold = 0.05)
cat("\n-- MR-PRESSO main results --\n"); print(pr$`Main MR results`)
cat("\n-- Global heterogeneity test p --\n")
cat(pr$`MR-PRESSO results`$`Global Test`$Pvalue, "\n")
cat("\n-- Outlier test --\n")
print(pr$`MR-PRESSO results`$`Outlier Test`)

cat("\n\n################################################################\n")
cat("# (b) DECISIVE TEST 1: OSA -> BMI  (is the instrument proxying obesity?)\n")
cat("################################################################\n")

bmi_out <- extract_outcome_data(snps = exp_dat$SNP, outcomes = bmi_id)
d2 <- harmonise_data(exp_dat, bmi_out, action = 2)
d2 <- d2[d2$mr_keep, ]
r2 <- mr(d2, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
cat("\n=== OSA -> BMI ===\n")
print(as.data.table(generate_odds_ratios(r2))[, .(method, nsnp, b, se, pval, or, or_lci95, or_uci95)])
cat("\n(interpretation: a large positive b here means the OSA instrument strongly\n")
cat(" raises BMI -> the OSA->HF estimate is confounded/proxied by adiposity)\n")

cat("\n\n################################################################\n")
cat("# (b) DECISIVE TEST 2: pure BMI -> Heart failure (positive control)\n")
cat("################################################################\n")

bmi_exp <- extract_instruments(outcomes = bmi_id, p1 = 5e-8, clump = TRUE,
                               r2 = 0.001, kb = 10000)
cat("BMI instruments:", nrow(bmi_exp), "\n")
bmi_exp$F <- (bmi_exp$beta.exposure / bmi_exp$se.exposure)^2
bmi_exp <- bmi_exp[bmi_exp$F > 10, ]
cat("BMI instruments with F>10:", nrow(bmi_exp), "\n")

hf_for_bmi <- extract_outcome_data(snps = bmi_exp$SNP, outcomes = hf_id)
d3 <- harmonise_data(bmi_exp, hf_for_bmi, action = 2)
d3 <- d3[d3$mr_keep, ]
r3 <- mr(d3, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
cat("\n=== BMI -> Heart failure (POSITIVE CONTROL) ===\n")
print(as.data.table(generate_odds_ratios(r3))[, .(method, nsnp, b, se, pval, or, or_lci95, or_uci95)])

cat("\n\n################################################################\n")
cat("# (b) DECISIVE TEST 3: MVMR  OSA + BMI -> Heart failure\n")
cat("################################################################\n")

if (!requireNamespace("MVMR", quietly = TRUE)) {
  cat("MVMR package NOT installed -- installing from mrcieu r-universe...\n")
  try(install.packages("MVMR", repos = c("https://mrcieu.r-universe.dev",
                                        "https://cloud.r-project.org")), silent = FALSE)
}
if (requireNamespace("MVMR", quietly = TRUE)) {
  cat("MVMR version:", as.character(packageVersion("MVMR")), "\n")
  library(MVMR)
  # build combined instrument set: union of OSA and BMI instruments
  iv <- unique(c(exp_dat$SNP, bmi_exp$SNP))
  cat("Union instrument set:", length(iv), "SNPs\n")
  ex <- extract_instruments(outcomes = c(osa_id, bmi_id))
  # keep only union IVs
  ex <- ex[ex$SNP %in% iv, ]
  ex <- format_data(ex, type = "exposure")
  out <- extract_outcome_data(snps = unique(ex$SNP), outcomes = hf_id)
  mvmr_dat <- mv_harmonise_data(ex, out)
  res_mv <- tryCatch(
    ivw_mvmr(mvmr_dat),
    error = function(e) paste("ivw_mvmr FAIL:", conditionMessage(e)))
  cat("\n-- MVMR-IVW (OSA + BMI -> HF) --\n")
  print(res_mv)
  # conditional F statistics
  cf <- tryCatch({
    Fmat <- mvmr_conditional_F(mvmr_dat)
    Fmat
  }, error = function(e) paste("conditional F FAIL:", conditionMessage(e)))
  cat("\n-- Conditional F statistics --\n"); print(cf)
} else {
  cat("MVMR unavailable -- skipping (network install may be blocked)\n")
}

cat("\n=== done ===\n")
