# 10_final_mvmr_and_summary.R
# Decisive analyses on the properly-powered MVP OSA instrument set.
#   (1) exact main results
#   (2) MVP OSA -> BMI
#   (3) MVMR: OSA + BMI -> Heart failure  (manual, with Sanderson conditional F)

suppressMessages({ library(TwoSampleMR); library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR

hf_id  <- "ebi-a-GCST009541"
bmi_id <- "ieu-b-40"

ivs <- fread(file.path(res_dir, "mvp_osa_instruments.tsv"))
cat("MVP OSA instruments loaded:", nrow(ivs), "\n")

# ---------- refit main model to print exact CIs ----------
out <- extract_outcome_data(snps = ivs$SNP, outcomes = hf_id)
exp_dat <- format_data(as.data.frame(ivs), type = "exposure",
                       snp_col = "SNP", beta_col = "beta", se_col = "se",
                       eaf_col = "eaf", effect_allele_col = "effect_allele",
                       other_allele_col = "other_allele", pval_col = "p")
dat <- harmonise_data(exp_dat, out, action = 2)
dat <- dat[dat$mr_keep, ]
r <- mr(dat, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
r <- generate_odds_ratios(r)
cat("\n================ MVP OSA -> Heart failure ================\n")
for (i in seq_len(nrow(r))) {
  cat(sprintf("%-22s nsnp=%3d  %s\n", r$method[i], r$nsnp[i],
              sprintf("OR=%.3f (%.3f-%.3f) p=%.3g", r$or[i], r$or_lci95[i], r$or_uci95[i], r$pval[i])))
}

# ---------- OSA -> BMI ----------
cat("\n================ MVP OSA -> BMI ================\n")
out_bmi <- extract_outcome_data(snps = ivs$SNP, outcomes = bmi_id)
d_bmi <- harmonise_data(exp_dat, out_bmi, action = 2)
d_bmi <- d_bmi[d_bmi$mr_keep, ]
rb <- mr(d_bmi, method_list = c("mr_ivw","mr_egger_regression","mr_weighted_median"))
for (i in seq_len(nrow(rb))) {
  cat(sprintf("%-22s nsnp=%3d  beta=%.4f (%.4f-%.4f) p=%.3g\n",
              rb$method[i], rb$nsnp[i], rb$b[i], rb$b[i]-1.96*rb$se[i],
              rb$b[i]+1.96*rb$se[i], rb$pval[i]))
}
cat("(SD units of BMI per log-odds genetic liability to OSA)\n")

# ---------- MVMR: OSA + BMI -> HF ----------
cat("\n================ MVMR: OSA + BMI -> Heart failure ================\n")

get_assoc <- function(id) {
  r <- tryCatch(ieugwasr::associations(variants = ivs$SNP, id = id),
                error = function(e) { cat("assoc failed for", id, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  as.data.table(r)[, .(SNP = rsid, ea, oa, beta, se)][!is.na(beta) & !is.na(se)]
}
a_osa <- get_assoc("finn-b-G6_SLEEPAPNO")   # OSA effects at MVP IVs (independent source)
a_bmi <- get_assoc(bmi_id)
out2  <- as.data.table(out)[, .(SNP, oa_o = other_allele.outcome,
                               ea_o = effect_allele.outcome,
                               by = beta.outcome, sy = se.outcome)]
cat("assoc rows -- OSA:", ifelse(is.null(a_osa),0,nrow(a_osa)),
    " BMI:", ifelse(is.null(a_bmi),0,nrow(a_bmi)),
    " outcome:", nrow(out2), "\n")

m <- merge(a_osa, a_bmi, by = "SNP", suffixes = c(".osa",".bmi"))
m <- merge(m, out2, by = "SNP")
cat("merged:", nrow(m), "\n")

al <- function(ea_e, oa_e, ea_o, oa_o)
  ifelse(ea_e==ea_o & oa_e==oa_o, 1L, ifelse(ea_e==oa_o & oa_e==ea_o, -1L, NA_integer_))
f1 <- al(m$ea.osa, m$oa.osa, m$ea_o, m$oa_o)
f2 <- al(m$ea.bmi, m$oa.bmi, m$ea_o, m$oa_o)
k  <- !is.na(f1) & !is.na(f2)
cat("harmonisable:", sum(k), "of", nrow(m), "\n")

D <- data.table(SNP=m$SNP[k], bx1=m$beta.osa[k]*f1[k], bx2=m$beta.bmi[k]*f2[k],
                sx1=m$se.osa[k], sx2=m$se.bmi[k], by=m$by[k], sy=m$sy[k])

if (nrow(D) >= 3) {
  w <- 1/D$sy^2
  X <- cbind(1, D$bx1, D$bx2)
  XtWX <- t(X) %*% (w*X)
  fit <- solve(XtWX, t(X) %*% (w*D$by))
  se  <- sqrt(diag(solve(XtWX)))
  cat("\n-- MVMR-IVW --\n")
  for (j in 1:2) {
    lab <- c("OSA (direct)","BMI (direct)")[j]
    b <- fit[j+1]; s <- se[j+1]
    cat(sprintf("%-14s beta=%+.4f (%.4f-%.4f) p=%.3g   OR=%.3f (%.3f-%.3f)\n",
                lab, b, b-1.96*s, b+1.96*s, 2*pnorm(-abs(b/s)),
                exp(b), exp(b-1.96*s), exp(b+1.96*s)))
  }

  # ---- Sanderson conditional F ----
  # Q_j = sum_k w_jk (bx_jk - g_j * bx_other,k)^2 , w_jk = 1/sx_jk^2
  # F_j = (Q_j - (K-1)) / (K - 2)
  K <- nrow(D)
  cf <- function(j) {
    y <- D[[paste0("bx",j)]]; s <- D[[paste0("sx",j)]]
    o <- D[[paste0("bx",3-j)]]
    wj <- 1/s^2
    g <- sum(wj*y*o)/sum(wj*o^2)
    Q <- sum(wj*(y - g*o)^2)
    (Q - (K-1))/(K-2)
  }
  cat(sprintf("\n-- Conditional F (Sanderson 2019) --\n"))
  cat(sprintf("F(OSA | BMI) = %.1f\n", cf(1)))
  cat(sprintf("F(BMI | OSA) = %.1f\n", cf(2)))
  cat("(>10 = adequate)\n")
} else cat("too few SNPs for MVMR\n")

cat("\n=== done ===\n")
