# 29_coloc_steiger_final.R
# Final fixes:
#  - coloc: drop SNPs whose MAF is unavailable (coloc rejects NA in MAF)
#  - steiger: drop SNPs with missing outcome EAF, then test direction

suppressMessages({ library(data.table); library(coloc); library(ieugwasr)
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
                   library(TwoSampleMR) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
hf_id <- "ebi-a-GCST009541"
N_OSA_CASE <- 152031; N_OSA_CTRL <- 278027
N_HF_CASE  <- 47309;  N_HF_CTRL  <- 930014

ivs <- fread(file.path(res_dir,"META_mvmr_input.tsv"))
comp <- fread(file.path(data_dir,"GCST90475824_p5e8.tsv"))[
          , .(SNP=rsid, chr=chromosome, pos=bp, ea=effect_allele,
              oa=other_allele, eaf, beta=log(or), p)]
ivs <- merge(ivs, comp[, .(SNP,chr,pos)], by="SNP", all.x=TRUE)
ivs <- ivs[!is.na(pos)]
cat("loci:", nrow(ivs), "\n")

# ================= COLOC =================
cat("\n############ COLOCALISATION ############\n")
run_coloc <- function(chr, pos, snp) {
  lo <- pos-5e5; hi <- pos+5e5
  reg <- comp[chr==get("chr") & pos>=lo & pos<=hi]
  if (nrow(reg) < 20) { cat(sprintf("\n%s: only %d local SNPs\n", snp, nrow(reg))); return(NULL) }
  get <- function(id) {
    r <- tryCatch(ieugwasr::associations(variants=reg$SNP, id=id), error=function(e) NULL)
    if (is.null(r)||nrow(r)==0) return(NULL)
    d <- as.data.table(r)
    d <- d[!is.na(beta) & !is.na(se) & se>0 & !is.na(eaf)]
    d[!duplicated(rsid)]
  }
  a <- get("finn-b-G6_SLEEPAPNO"); b <- get(hf_id)
  if (is.null(a)||is.null(b)) { cat(sprintf("\n%s: missing data\n", snp)); return(NULL) }
  m <- merge(a[, .(SNP=rsid, ea, nea, eaf, beta, se)],
             b[, .(SNP=rsid, ea2=ea, nea2=nea, eaf2=eaf, beta2=beta, se2=se)], by="SNP")
  m[, flip := ifelse(ea==ea2 & nea==nea2, 1, ifelse(ea==nea2 & nea==ea2, -1, NA_integer_))]
  m <- m[!is.na(flip)]
  m[, beta2a := beta2*flip]
  m <- m[!is.na(eaf) & !is.na(eaf2) & eaf>0 & eaf<1 & eaf2>0 & eaf2<1]
  cat(sprintf("\n%s (chr%d %.2f-%.2f Mb): %d SNPs\n", snp, chr, lo/1e6, hi/1e6, nrow(m)))
  if (nrow(m) < 30) { cat("  too few\n"); return(NULL) }
  d1 <- list(snp=m$SNP, beta=m$beta, varbeta=m$se^2, type="cc",
             N=N_OSA_CASE+N_OSA_CTRL, MAF=pmin(m$eaf,1-m$eaf),
             s=N_OSA_CASE/(N_OSA_CASE+N_OSA_CTRL))
  d2 <- list(snp=m$SNP, beta=m$beta2a, varbeta=m$se2^2, type="cc",
             N=N_HF_CASE+N_HF_CTRL, MAF=pmin(m$eaf2,1-m$eaf2),
             s=N_HF_CASE/(N_HF_CASE+N_HF_CTRL))
  r <- tryCatch(coloc.abf(d1,d2), error=function(e){cat("  err:",conditionMessage(e),"\n");NULL})
  if (is.null(r)) return(NULL)
  pp <- r$summary
  cat(sprintf("  PP.H3=%.3f  PP.H4=%.3f\n", pp[["PP.H3.abf"]], pp[["PP.H4.abf"]]))
  data.table(SNP=snp, chr=chr, pos=pos, nsnp=nrow(m),
             PP.H0=pp[["PP.H0.abf"]], PP.H1=pp[["PP.H1.abf"]],
             PP.H2=pp[["PP.H2.abf"]], PP.H3=pp[["PP.H3.abf"]], PP.H4=pp[["PP.H4.abf"]])
}
out <- list()
for (i in seq_len(nrow(ivs))) {
  r <- tryCatch(run_coloc(ivs$chr[i], ivs$pos[i], ivs$SNP[i]), error=function(e) NULL)
  if (!is.null(r)) out[[length(out)+1]] <- r
}
if (length(out)>0) {
  R <- rbindlist(out)
  cat("\n===== COLOC SUMMARY =====\n")
  print(R[, .(SNP, nsnp, PP.H3=round(PP.H3,3), PP.H4=round(PP.H4,3))])
  fwrite(R, file.path(res_dir,"COLOC_results.csv"))
  cat("\nPP.H4>=0.8:", sum(R$PP.H4>=0.8), "of", nrow(R),
      " | PP.H3>=0.8:", sum(R$PP.H3>=0.8), "\n")
} else cat("\nnone colocalised\n")

# ================= STEIGER =================
cat("\n\n############ STEIGER ############\n")
meta <- fread(file.path(res_dir,"meta_osa_snps.tsv"))
mv <- fread(file.path(data_dir,"GCST90475824_p5e8.tsv"))[
        , .(SNP=rsid, effect_allele, other_allele, eaf)]
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP,pval=meta$p_meta),
             clump_r2=0.001, clump_kb=10000, pop="EUR"), error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
ivs2 <- merge(meta[, .(SNP,beta=beta_meta,se=se_meta,p=p_meta)], mv, by="SNP")

ex <- format_data(as.data.frame(ivs2), type="exposure", snp_col="SNP",
      beta_col="beta", se_col="se", eaf_col="eaf",
      effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")
out2 <- extract_outcome_data(snps=ivs2$SNP, outcomes=hf_id)
dat <- harmonise_data(ex, out2, action=2); dat <- dat[dat$mr_keep,]
cat("harmonised:", nrow(dat), "\n")

r2e <- 2*dat$eaf.exposure*(1-dat$eaf.exposure)*dat$beta.exposure^2
r2o <- 2*dat$eaf.outcome *(1-dat$eaf.outcome )*dat$beta.outcome^2
keep <- is.finite(r2e) & is.finite(r2o) & r2e>0 & r2o>0
cat("SNPs with both R2 computable:", sum(keep), "of", nrow(dat), "\n")

if (sum(keep) >= 3) {
  res <- data.table(SNP=dat$SNP[keep], r2_exposure=r2e[keep],
                    r2_outcome=r2o[keep], ratio=r2e[keep]/r2o[keep],
                    correct_dir = r2e[keep] > r2o[keep])
  cat("\nsum R2 OSA:", sprintf("%.5f", sum(res$r2_exposure)), "\n")
  cat("sum R2 HF :", sprintf("%.5f", sum(res$r2_outcome)), "\n")
  cat("ratio    :", sprintf("%.1f", sum(res$r2_exposure)/sum(res$r2_outcome)), "\n")
  print(res)
  cat("\ncorrect direction:", sum(res$correct_dir), "/", nrow(res), "\n")
  bt <- binom.test(sum(res$correct_dir), nrow(res), 0.5)
  cat("binom p =", signif(bt$p.value,4), "\n")
  fwrite(res, file.path(res_dir,"STEIGER_results.csv"))
} else cat("too few SNPs for Steiger\n")
cat("\n=== done ===\n")
