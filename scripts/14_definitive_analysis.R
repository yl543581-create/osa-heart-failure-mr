# 14_definitive_analysis.R
# Definitive analysis on the PROPERLY LD-CLUMPED instrument set (r2<0.001, 10 Mb).
#   (1) main MR: OSA -> heart failure
#   (2) OSA -> BMI
#   (3) MVMR: OSA + BMI -> HF  with Sanderson conditional F
#   (4) FTO-region sensitivity
#   (5) colocalisation-ready outputs

suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
hf_id   <- "ebi-a-GCST009541"
bmi_id  <- "ieu-b-40"

ivs <- fread(file.path(res_dir, "mvp_osa_instruments_LDclumped.tsv"))
ivs <- ivs[, .(SNP = rsid, beta = beta2, se = se2, p,
               effect_allele, other_allele, eaf, chromosome, position = bp)]
cat("LD-clumped instruments:", nrow(ivs), "\n")
cat("mean F =", round(mean((ivs$beta/ivs$se)^2), 1),
    " | min F =", round(min((ivs$beta/ivs$se)^2), 1), "\n\n")

fmt_exp <- function(dt) format_data(as.data.frame(dt), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
  effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")

# ================= (1) main MR =================
cat("################ (1) OSA -> Heart failure ################\n")
out <- extract_outcome_data(snps = ivs$SNP, outcomes = hf_id)
dat <- harmonise_data(fmt_exp(ivs), out, action = 2); dat <- dat[dat$mr_keep, ]
r <- generate_odds_ratios(mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode")))
for (i in seq_len(nrow(r)))
  cat(sprintf("%-22s nsnp=%3d OR=%.3f (%.3f-%.3f) p=%.3g\n",
      r$method[i], r$nsnp[i], r$or[i], r$or_lci95[i], r$or_uci95[i], r$pval[i]))
het <- mr_heterogeneity(dat); ple <- mr_pleiotropy_test(dat)
cat("\nCochran Q p =", signif(het$Q_pval[het$method=="Inverse variance weighted"],3), "\n")
cat("Egger intercept =", signif(ple$egger_intercept,4), " p =", signif(ple$pval,3), "\n")
fwrite(as.data.table(r), file.path(res_dir,"FINAL_main_OSA_HF.csv"))
saveRDS(dat, file.path(res_dir,"FINAL_harmonised_dat.rds"))

# ================= (2) OSA -> BMI =================
cat("\n################ (2) OSA -> BMI ################\n")
outb <- extract_outcome_data(snps = ivs$SNP, outcomes = bmi_id)
db <- harmonise_data(fmt_exp(ivs), outb, action = 2); db <- db[db$mr_keep, ]
rb <- mr(db, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median"))
for (i in seq_len(nrow(rb)))
  cat(sprintf("%-22s nsnp=%3d beta=%+.4f (%.4f-%.4f) p=%.3g\n",
      rb$method[i], rb$nsnp[i], rb$b[i], rb$b[i]-1.96*rb$se[i],
      rb$b[i]+1.96*rb$se[i], rb$pval[i]))

# ================= (3) MVMR =================
cat("\n################ (3) MVMR: OSA + BMI -> HF ################\n")
get_assoc <- function(id) {
  rr <- tryCatch(ieugwasr::associations(variants = ivs$SNP, id = id), error=function(e) NULL)
  if (is.null(rr)) return(NULL)
  d <- as.data.table(rr)
  d <- d[!is.na(beta) & !is.na(se) & !is.na(ea) & !is.na(nea)]
  d[, .(SNP = rsid, ea, oa = nea, beta, se)]
}
a_osa <- get_assoc("finn-b-G6_SLEEPAPNO")
a_bmi <- get_assoc(bmi_id)
o2 <- as.data.table(out)[, .(SNP, ea_o=effect_allele.outcome,
                            oa_o=other_allele.outcome, by=beta.outcome, sy=se.outcome)]
m <- merge(merge(a_osa, a_bmi, by="SNP", suffixes=c(".osa",".bmi")), o2, by="SNP")
al <- function(e1,o1,e2,o2v) ifelse(e1==e2 & o1==o2v, 1L, ifelse(e1==o2v & o1==e2, -1L, NA_integer_))
f1 <- al(m$ea.osa,m$oa.osa,m$ea_o,m$oa_o); f2 <- al(m$ea.bmi,m$oa.bmi,m$ea_o,m$oa_o)
k <- !is.na(f1) & !is.na(f2)
D <- data.table(SNP=m$SNP[k], bx1=m$beta.osa[k]*f1[k], bx2=m$beta.bmi[k]*f2[k],
                sx1=m$se.osa[k], sx2=m$se.bmi[k], by=m$by[k], sy=m$sy[k])
D <- D[is.finite(bx1)&is.finite(bx2)&is.finite(by)&sy>0]
cat("MVMR SNPs:", nrow(D), "\n")

w <- 1/D$sy^2; X <- cbind(1,D$bx1,D$bx2); XtWX <- t(X)%*%(w*X)
fit <- solve(XtWX, t(X)%*%(w*D$by)); se <- sqrt(diag(solve(XtWX)))
res <- data.table()
for (j in 1:2) {
  b <- fit[j+1]; s <- se[j+1]
  res <- rbind(res, data.table(exposure=c("OSA (direct)","BMI (direct)")[j], beta=b, se=s,
              pval=2*pnorm(-abs(b/s)), or=exp(b), or_l=exp(b-1.96*s), or_u=exp(b+1.96*s)))
}
for (i in seq_len(nrow(res)))
  cat(sprintf("%-14s beta=%+.4f (%.4f to %.4f) p=%.3g | OR=%.3f (%.3f-%.3f)\n",
      res$exposure[i], res$beta[i], res$beta[i]-1.96*res$se[i], res$beta[i]+1.96*res$se[i],
      res$pval[i], res$or[i], res$or_l[i], res$or_u[i]))

K <- nrow(D)
cf <- function(j) { y<-D[[paste0("bx",j)]]; s<-D[[paste0("sx",j)]]; o<-D[[paste0("bx",3-j)]]
  wj<-1/s^2; g<-sum(wj*y*o)/sum(wj*o^2); Q<-sum(wj*(y-g*o)^2); (Q-(K-1))/(K-2) }
cat(sprintf("\nConditional F (Sanderson 2019):\n  F(OSA|BMI) = %.1f\n  F(BMI|OSA) = %.1f\n", cf(1), cf(2)))

uv <- function(bx,lab){ b<-sum(w*bx*D$by)/sum(w*bx^2); s<-sqrt(1/sum(w*bx^2))
  cat(sprintf("  %-18s beta=%+.4f p=%.3g OR=%.3f\n", lab, b, 2*pnorm(-abs(b/s)), exp(b))) }
cat("\nUnivariable on same SNPs:\n"); uv(D$bx1,"OSA univariable"); uv(D$bx2,"BMI univariable")
fwrite(res, file.path(res_dir,"FINAL_mvmr_osa_bmi_hf.csv"))
fwrite(D, file.path(res_dir,"FINAL_mvmr_input.tsv"), sep="\t")

# ================= (4) FTO sensitivity =================
cat("\n################ (4) FTO-region sensitivity ################\n")
nofto <- ivs[!(chromosome==16 & position>=53.2e6 & position<=54.3e6)]
cat("removed", nrow(ivs)-nrow(nofto), "FTO-region SNP(s)\n")
out2 <- extract_outcome_data(snps=nofto$SNP, outcomes=hf_id)
d2 <- harmonise_data(fmt_exp(nofto), out2, action=2); d2 <- d2[d2$mr_keep,]
r2 <- generate_odds_ratios(mr(d2, method_list=c("mr_ivw","mr_weighted_median")))
for (i in seq_len(nrow(r2)))
  cat(sprintf("%-22s nsnp=%3d OR=%.3f (%.3f-%.3f) p=%.3g\n",
      r2$method[i], r2$nsnp[i], r2$or[i], r2$or_lci95[i], r2$or_uci95[i], r2$pval[i]))

cat("\n=== ALL DONE ===\n")
