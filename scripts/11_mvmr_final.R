# 11_mvmr_final.R
# MVMR: OSA + BMI -> Heart failure, with robust handling of OpenGWAS columns.

suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
hf_id  <- "ebi-a-GCST009541"
bmi_id <- "ieu-b-40"

ivs <- fread(file.path(res_dir, "mvp_osa_instruments.tsv"))
cat("IVs:", nrow(ivs), "\n")

# --- robust association fetch ---
get_assoc <- function(id) {
  r <- tryCatch(ieugwasr::associations(variants = ivs$SNP, id = id),
                error = function(e) { cat("FAIL", id, ":", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) return(NULL)
  d <- as.data.table(r)
  cat("\ncolumns for", id, ":", paste(names(d), collapse = ", "), "\n")
  # normalise column names
  nm <- names(d)
  pick <- function(cands) { h <- intersect(cands, nm); if (length(h)) h[1] else NA_character_ }
  snp_c <- pick(c("rsid","SNP","snp"))
  ea_c  <- pick(c("ea","effect_allele","A1"))
  oa_c  <- pick(c("nea","other_allele","oa","A2"))
  b_c   <- pick(c("beta","b"))
  se_c  <- pick(c("se","standard_error"))
  out <- data.table(SNP = d[[snp_c]], ea = d[[ea_c]], oa = d[[oa_c]],
                    beta = suppressWarnings(as.numeric(d[[b_c]])),
                    se   = suppressWarnings(as.numeric(d[[se_c]])))
  out[!is.na(beta) & !is.na(se) & !is.na(ea) & !is.na(oa)]
}

a_osa <- get_assoc("finn-b-G6_SLEEPAPNO")
a_bmi <- get_assoc(bmi_id)
cat("\nrows -- OSA:", ifelse(is.null(a_osa),0,nrow(a_osa)),
    " BMI:", ifelse(is.null(a_bmi),0,nrow(a_bmi)), "\n")

out <- extract_outcome_data(snps = ivs$SNP, outcomes = hf_id)
out <- as.data.table(out)[, .(SNP, oa_o = other_allele.outcome, ea_o = effect_allele.outcome,
                              by = beta.outcome, sy = se.outcome)]
cat("outcome rows:", nrow(out), "\n")

m <- merge(a_osa, a_bmi, by = "SNP", suffixes = c(".osa",".bmi"))
m <- merge(m, out, by = "SNP")
cat("merged:", nrow(m), "\n")

al <- function(ea_e, oa_e, ea_o, oa_o)
  ifelse(ea_e==ea_o & oa_e==oa_o, 1L, ifelse(ea_e==oa_o & oa_e==ea_o, -1L, NA_integer_))
f1 <- al(m$ea.osa, m$oa.osa, m$ea_o, m$oa_o)
f2 <- al(m$ea.bmi, m$oa.bmi, m$ea_o, m$oa_o)
k  <- !is.na(f1) & !is.na(f2)
cat("harmonisable:", sum(k), "of", nrow(m), "\n")

D <- data.table(SNP = m$SNP[k], bx1 = m$beta.osa[k]*f1[k], bx2 = m$beta.bmi[k]*f2[k],
                sx1 = m$se.osa[k], sx2 = m$se.bmi[k], by = m$by[k], sy = m$sy[k])
D <- D[is.finite(bx1) & is.finite(bx2) & is.finite(by) & sy > 0]
cat("final MVMR SNPs:", nrow(D), "\n")

if (nrow(D) >= 5) {
  w <- 1/D$sy^2
  X <- cbind(1, D$bx1, D$bx2)
  XtWX <- t(X) %*% (w*X)
  fit <- solve(XtWX, t(X) %*% (w*D$by))
  se  <- sqrt(diag(solve(XtWX)))

  cat("\n================ MVMR-IVW: OSA + BMI -> Heart failure ================\n")
  res <- data.table()
  for (j in 1:2) {
    lab <- c("OSA (direct)","BMI (direct)")[j]
    b <- fit[j+1]; s <- se[j+1]
    res <- rbind(res, data.table(exposure=lab, beta=b, se=s,
                pval=2*pnorm(-abs(b/s)), or=exp(b),
                or_l=exp(b-1.96*s), or_u=exp(b+1.96*s)))
  }
  for (i in seq_len(nrow(res)))
    cat(sprintf("%-14s beta=%+.4f (%.4f to %.4f) p=%.3g | OR=%.3f (%.3f-%.3f)\n",
        res$exposure[i], res$beta[i], res$beta[i]-1.96*res$se[i],
        res$beta[i]+1.96*res$se[i], res$pval[i], res$or[i], res$or_l[i], res$or_u[i]))

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

  # univariable on the same SNP set for comparison
  cat("\n-- Univariable on the SAME SNP set --\n")
  uv <- function(bx, lab) {
    b <- sum(w*bx*D$by)/sum(w*bx^2); s <- sqrt(1/sum(w*bx^2))
    cat(sprintf("%-20s beta=%+.4f p=%.3g OR=%.3f\n", lab, b, 2*pnorm(-abs(b/s)), exp(b)))
  }
  uv(D$bx1, "OSA univariable"); uv(D$bx2, "BMI univariable")

  fwrite(res, file.path(res_dir, "mvmr_osa_bmi_hf.csv"))
  fwrite(D,  file.path(res_dir, "mvmr_input_snps.tsv"), sep = "\t")
  cat("\nsaved results\n")
} else cat("too few SNPs\n")
cat("=== done ===\n")
