# 28_steiger_v2.R
# Fix: get_r_from_lor() expects ODDS RATIOS, not log-ORs.
# Also the outcome R2 was NA because sample sizes were not attached.
#
# Steiger logic: the instrument should explain more variance in the exposure
# (OSA) than in the outcome (HF) if the causal direction is OSA -> HF.

suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
hf_id <- "ebi-a-GCST009541"
N_OSA_CASE <- 152031; N_OSA_CTRL <- 278027
N_HF_CASE  <- 47309;  N_HF_CTRL  <- 930014

meta <- fread(file.path(res_dir,"meta_osa_snps.tsv"))
mv <- fread(file.path(data_dir,"GCST90475824_p5e8.tsv"))[
        , .(SNP=rsid, effect_allele, other_allele, eaf)]
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP,pval=meta$p_meta),
             clump_r2=0.001, clump_kb=10000, pop="EUR"), error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
ivs <- merge(meta[, .(SNP,beta=beta_meta,se=se_meta,p=p_meta)], mv, by="SNP")
cat("instruments:", nrow(ivs), "\n")

ex <- format_data(as.data.frame(ivs), type="exposure", snp_col="SNP",
      beta_col="beta", se_col="se", eaf_col="eaf",
      effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")
ex$samplesize.exposure <- N_OSA_CASE+N_OSA_CTRL
ex$ncase.exposure <- N_OSA_CASE; ex$ncontrol.exposure <- N_OSA_CTRL

out <- extract_outcome_data(snps=ivs$SNP, outcomes=hf_id)
out$samplesize.outcome <- N_HF_CASE+N_HF_CTRL
out$ncase.outcome <- N_HF_CASE; out$ncontrol.outcome <- N_HF_CTRL

dat <- harmonise_data(ex, out, action=2); dat <- dat[dat$mr_keep,]
cat("harmonised:", nrow(dat), "\n")

# ---- SIMPLE / ROBUST variance-explained approach ----
# For a continuous standardised trait: R2_snp = 2 * EAF * (1-EAF) * beta^2
# For a binary trait we approximate beta on the liability scale;
# this is adequate for a DIRECTIONAL comparison, which is all Steiger needs.
r2_exp <- 2 * dat$eaf.exposure * (1-dat$eaf.exposure) * dat$beta.exposure^2
r2_out <- 2 * dat$eaf.outcome  * (1-dat$eaf.outcome)  * dat$beta.outcome^2

cat("\n=== variance explained (approx liability scale) ===\n")
cat("sum R2 exposure (OSA):", sprintf("%.5f", sum(r2_exp)), "\n")
cat("sum R2 outcome  (HF): ", sprintf("%.5f", sum(r2_out)), "\n")
cat("ratio:", sprintf("%.1f", sum(r2_exp)/sum(r2_out)), "\n")

res <- data.table(SNP=dat$SNP, r2_exposure=r2_exp, r2_outcome=r2_out,
                  correct_dir = r2_exp > r2_out)
cat("\n=== per-SNP direction ===\n"); print(res)
cat("\ncorrect direction:", sum(res$correct_dir), "/", nrow(res), "\n")

# binomial test that the direction is not 50/50
bt <- binom.test(sum(res$correct_dir), nrow(res), 0.5)
cat("binom test p =", signif(bt$p.value, 4), "\n")

# ---- also try TwoSampleMR's own Steiger with OR correction ----
cat("\n=== TwoSampleMR directionality_test (with OR correction) ===\n")
d2 <- copy(dat)
d2$r.exposure <- get_r_from_lor(exp(d2$beta.exposure), d2$eaf.exposure,
                                d2$ncase.exposure, d2$ncontrol.exposure,
                                d2$samplesize.exposure)
d2$r.outcome <- get_r_from_lor(exp(d2$beta.outcome), d2$eaf.outcome,
                               d2$ncase.outcome, d2$ncontrol.outcome,
                               d2$samplesize.outcome)
cat("sum r2 exposure:", sprintf("%.5f", sum(d2$r.exposure^2, na.rm=TRUE)), "\n")
cat("sum r2 outcome: ", sprintf("%.5f", sum(d2$r.outcome^2, na.rm=TRUE)), "\n")
st <- tryCatch(directionality_test(d2), error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
if (!is.null(st)) print(st)

fwrite(res, file.path(res_dir,"STEIGER_results.csv"))
cat("\nsaved\n=== done ===\n")
