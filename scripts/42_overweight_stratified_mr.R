# 42_overweight_stratified_mr.R
# Natural-experiment test of the adiposity-mediation model.
#
# If OSA's heart-failure association is adiposity-mediated, then stratifying the
# outcome by adiposity should CHANGE the effect size:
#
#   I9_HEARTFAIL_AND_OVERWEIGHT (HF + BMI>=25)  -> adiposity-enriched
#   I9_HEARTFAIL_AND_CHD        (HF + CHD)      -> adiposity-neutral comparator
#   I9_HEARTFAIL                (strict HF)     -> reference
#   I9_HEARTFAIL_AND_HYPERTCARDIOM              -> non-adiposity comparator
#
# CAVEAT: stratified endpoints differ in case/control ratio and ascertainment,
# so a naive comparison of log-ORs is not strictly like-for-like. We report it
# as a descriptive, hypothesis-consistent pattern and say so.

suppressMessages({ library(data.table); library(TwoSampleMR); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

# ---- instrument ----
meta <- fread(file.path(res_dir,"meta_osa_snps.tsv"))
mv <- fread(file.path(data_dir,"GCST90475824_p5e8.tsv"))[
        , .(SNP=rsid, effect_allele, other_allele, eaf)]
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP, pval=meta$p_meta),
             clump_r2=0.001, clump_kb=10000, pop="EUR"), error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
ivs <- merge(meta[, .(SNP, beta=beta_meta, se=se_meta, p=p_meta)], mv, by="SNP")
cat("instrument SNPs:", nrow(ivs), "\n")

ex <- format_data(as.data.frame(ivs), type="exposure", snp_col="SNP",
      beta_col="beta", se_col="se", eaf_col="eaf",
      effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")

# ---- endpoint files ----
endpoints <- data.table(
  label = c("Strict HF", "HF + BMI>=25", "HF + CHD", "HF + hypertrophic CM"),
  file  = c("finngen_R13_I9_HEARTFAIL.gz",
            "finngen_R13_I9_HEARTFAIL_AND_OVERWEIGHT.gz",
            "finngen_R13_I9_HEARTFAIL_AND_CHD.gz",
            "finngen_R13_I9_HEARTFAIL_AND_HYPERTCARDIOM.gz"),
  ncase = c(41591, 25129, 26459, 622),
  nctrl = c(458595, 204333, 409472, 498666)
)

cols <- c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt")
res <- list()
for (i in seq_len(nrow(endpoints))) {
  f <- file.path(data_dir, endpoints$file[i])
  if (!file.exists(f)) { cat(sprintf("%-22s file missing\n", endpoints$label[i])); next }
  cat(sprintf("\n=== %s (%d cases) ===\n", endpoints$label[i], endpoints$ncase[i]))
  d <- fread(f, select=cols, sep="\t", header=TRUE, nThread=4, showProgress=FALSE)
  setnames(d, c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
              c("chr","pos","other_allele.outcome","effect_allele.outcome","SNP",
                "pval.outcome","beta.outcome","se.outcome","eaf.outcome"))
  d[, SNP := tstrsplit(SNP, ",", fixed=TRUE)[[1]]]
  d <- d[!is.na(SNP) & !is.na(beta.outcome) & !is.na(se.outcome) & se.outcome>0]
  d <- d[!duplicated(SNP)]
  o <- d[SNP %in% ivs$SNP]
  cat("  instrument SNPs in file:", nrow(o), "\n")
  if (nrow(o) < 3) next
  o <- as.data.frame(o)
  o$id.outcome <- endpoints$label[i]
  o$outcome <- endpoints$label[i]
  o$mr_keep.outcome <- TRUE
  h <- tryCatch(harmonise_data(ex, o, action=2), error=function(e){cat("  err:",conditionMessage(e),"\n");NULL})
  if (is.null(h)) next
  h <- h[h$mr_keep, ]
  cat("  harmonised:", nrow(h), "\n")
  if (nrow(h) < 3) next
  r <- generate_odds_ratios(mr(h, method_list=c("mr_ivw","mr_weighted_median")))
  ri <- r[r$method=="Inverse variance weighted", ]
  # convert to per-SD-ish comparable scale
  res[[length(res)+1]] <- data.table(
    endpoint=endpoints$label[i], ncase=endpoints$ncase[i], nctrl=endpoints$nctrl[i],
    nsnp=ri$nsnp, b=ri$b, se=ri$se, pval=ri$pval,
    OR=ri$or, lo=ri$or_lci95, hi=ri$or_uci95)
  cat(sprintf("  IVW: OR=%.3f (%.3f-%.3f) p=%.3g  [b=%.4f]\n",
              ri$or, ri$or_lci95, ri$or_uci95, ri$pval, ri$b))
}

if (length(res) > 0) {
  R <- rbindlist(res)
  cat("\n\n========= STRATIFIED COMPARISON =========\n")
  print(R[, .(endpoint, ncase, nsnp, b=round(b,4), OR=round(OR,3),
              lo=round(lo,3), hi=round(hi,3), pval=signif(pval,3))])
  fwrite(R, file.path(res_dir,"STRATIFIED_HF_comparison.csv"))

  cat("\n--- adiposity-stratification contrast ---\n")
  ov <- R[endpoint=="HF + BMI>=25"]; st <- R[endpoint=="Strict HF"]
  if (nrow(ov)>0 && nrow(st)>0) {
    cat(sprintf("HF+BMI>=25 : b = %+.4f\n", ov$b))
    cat(sprintf("Strict HF  : b = %+.4f\n", st$b))
    cat(sprintf("difference : %+.4f\n", ov$b - st$b))
    # crude z for the difference (assumes independence; approximate)
    z <- (ov$b - st$b)/sqrt(ov$se^2 + st$se^2)
    cat(sprintf("z = %.2f, p = %.3g (approximate, assumes independence)\n",
                z, 2*pnorm(-abs(z))))
  }
} else cat("\nno results\n")
cat("\n=== done ===\n")
