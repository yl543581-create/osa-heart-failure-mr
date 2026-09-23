# 26_steiger.R
# Steiger directionality filtering: does the genetic instrument explain more
# variance in OSA (exposure) than in heart failure (outcome)?
# If yes, the causal direction OSA -> HF is supported over reverse causation.

suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
hf_id <- "ebi-a-GCST009541"

meta <- fread(file.path(res_dir, "meta_osa_snps.tsv"))
mv <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))[
        , .(SNP=rsid, effect_allele, other_allele, eaf)]
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP, pval=meta$p_meta),
                                  clump_r2=0.001, clump_kb=10000, pop="EUR"), error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
ivs <- merge(meta[, .(SNP, beta=beta_meta, se=se_meta, p=p_meta)], mv, by="SNP")
cat("instruments:", nrow(ivs), "\n")

N_OSA_CASE <- 152031; N_OSA_CTRL <- 278027
N_HF_CASE  <- 47309;  N_HF_CTRL  <- 930014

ex <- format_data(as.data.frame(ivs), type="exposure",
      snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
      effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")
ex$samplesize.exposure <- N_OSA_CASE + N_OSA_CTRL
ex$ncase.exposure <- N_OSA_CASE
ex$ncontrol.exposure <- N_OSA_CTRL

out <- extract_outcome_data(snps = ivs$SNP, outcomes = hf_id)
out$samplesize.outcome <- N_HF_CASE + N_HF_CTRL
out$ncase.outcome <- N_HF_CASE
out$ncontrol.outcome <- N_HF_CTRL

dat <- harmonise_data(ex, out, action = 2)
dat <- dat[dat$mr_keep, ]
cat("harmonised:", nrow(dat), "\n")

# ---- Steiger: variance explained in exposure vs outcome ----
dat$r.exposure <- get_r_from_lor(
  dat$beta.exposure, dat$eaf.exposure, dat$ncase.exposure,
  dat$ncontrol.exposure, dat$samplesize.exposure)
dat$r.outcome <- get_r_from_lor(
  dat$beta.outcome, dat$eaf.outcome, dat$ncase.outcome,
  dat$ncontrol.outcome, dat$samplesize.outcome)

cat("\n=== variance explained ===\n")
cat("R2 exposure (OSA):  ", sprintf("%.4f%%", 100*sum(dat$r.exposure^2)), "\n")
cat("R2 outcome  (HF):   ", sprintf("%.4f%%", 100*sum(dat$r.outcome^2)), "\n")
cat("ratio R2_OSA/R2_HF:", sprintf("%.1f", sum(dat$r.exposure^2)/sum(dat$r.outcome^2)), "\n")

cat("\n=== directionality test (Steiger) ===\n")
st <- tryCatch(directionality_test(dat), error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
if (!is.null(st)) print(st)

cat("\n=== per-SNP correct direction ===\n")
dat$correct_dir <- dat$r.exposure^2 > dat$r.outcome^2
cat("correctly oriented:", sum(dat$correct_dir), "/", nrow(dat), "\n")
print(dat[, c("SNP","r.exposure","r.outcome","correct_dir")])

fwrite(as.data.table(dat)[, .(SNP, r.exposure, r.outcome, correct_dir)],
       file.path(res_dir, "STEIGER_results.csv"))
cat("\nsaved\n=== done ===\n")
