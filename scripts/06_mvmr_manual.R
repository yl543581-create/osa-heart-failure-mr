# 06_mvmr_manual.R
# MVMR package is unavailable (github times out, r-universe 403s).
# Implement MVMR-IVW + conditional F directly from summary statistics.
#
#   MVMR-IVW: WLS of beta.outcome on (beta.exposure1, beta.exposure2),
#             weights = 1 / se.outcome^2
#   Conditional F (Sanderson, Davey Smith, Windmeijer & Bowden 2019, Stat Med):
#             uses Q_j = sum_k w_k * (beta_Xj,k - (sum_j' gamma_jj' beta_Xj',k))^2
#             approximated here via the standard summary-data formulation.

suppressMessages({ library(TwoSampleMR); library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

osa_id <- "finn-b-G6_SLEEPAPNO"
bmi_id <- "ieu-b-40"
hf_id  <- "ebi-a-GCST009541"

cat("=== instrument sets ===\n")
ex_osa <- extract_instruments(outcomes = osa_id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
ex_bmi <- extract_instruments(outcomes = bmi_id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
cat("OSA IVs:", nrow(ex_osa), " | BMI IVs:", nrow(ex_bmi), "\n")
ivs <- unique(c(ex_osa$SNP, ex_bmi$SNP))
cat("union IVs:", length(ivs), "\n\n")

# effects of BOTH exposures at ALL union IVs
assoc <- function(id) {
  r <- tryCatch(ieugwasr::associations(variants = ivs, id = id), error = function(e) NULL)
  if (is.null(r) || nrow(r) == 0) return(NULL)
  as.data.table(r)[, .(SNP = rsid, ea, oa, eaf, beta, se)][!is.na(beta) & !is.na(se)]
}
a_osa <- assoc(osa_id)
a_bmi <- assoc(bmi_id)
if (is.null(a_osa) || is.null(a_bmi)) stop("failed to fetch exposure associations")
cat("assoc OSA:", nrow(a_osa), " | assoc BMI:", nrow(a_bmi), "\n")

out <- extract_outcome_data(snps = ivs, outcomes = hf_id)
out <- as.data.table(out)[, .(SNP, oa_out = other_allele.outcome,
                              ea_out = effect_allele.outcome,
                              beta_out = beta.outcome, se_out = se.outcome)]
cat("outcome rows:", nrow(out), "\n")

m <- merge(a_osa, a_bmi, by = "SNP", suffixes = c(".osa", ".bmi"))
m <- merge(m, out, by = "SNP")
cat("merged:", nrow(m), "SNPs\n")

align <- function(ea_e, oa_e, ea_o, oa_o) {
  ifelse(ea_e == ea_o & oa_e == oa_o,  1L,
  ifelse(ea_e == oa_o & oa_e == ea_o, -1L, NA_integer_))
}
f1 <- align(m$ea.osa, m$oa.osa, m$ea_out, m$oa_out)
f2 <- align(m$ea.bmi, m$oa.bmi, m$ea_out, m$oa_out)
keep <- !is.na(f1) & !is.na(f2)
cat("harmonisables:", sum(keep), "\n\n")

D <- data.table(
  SNP = m$SNP[keep],
  bx1 = m$beta.osa[keep] * f1[keep],      # OSA
  bx2 = m$beta.bmi[keep] * f2[keep],      # BMI
  sx1 = m$se.osa[keep],
  sx2 = m$se.bmi[keep],
  by  = m$beta_out[keep],
  sy  = m$se_out[keep]
)
print(D)

# ---------- MVMR-IVW ----------
w  <- 1 / D$sy^2
X  <- cbind(1, D$bx1, D$bx2)
XtWX <- t(X) %*% (w * X)
fit  <- solve(XtWX, t(X) %*% (w * D$by))
se   <- sqrt(diag(solve(XtWX)))

res <- data.table(exposure = c("OSA (direct effect)", "BMI (direct effect)"),
                  beta = as.numeric(fit[2:3]), se = as.numeric(se[2:3]))
res[, `:=`(pval = 2 * pnorm(-abs(beta / se)), or = exp(beta),
           or_l = exp(beta - 1.96 * se), or_u = exp(beta + 1.96 * se))]

cat("\n=== MVMR-IVW: OSA + BMI -> Heart failure ===\n")
print(res)

# ---------- univariable on the SAME SNP set, for comparability ----------
uv <- function(bx, sx, lab) {
  wi <- 1 / D$sy^2
  b  <- sum(wi * bx * D$by) / sum(wi * bx^2)
  s  <- sqrt(1 / sum(wi * bx^2))
  data.table(exposure = lab, beta = b, se = s, pval = 2 * pnorm(-abs(b / s)),
             or = exp(b), or_l = exp(b - 1.96 * s), or_u = exp(b + 1.96 * s))
}
cat("\n=== Univariable on the SAME SNPs (comparability) ===\n")
print(rbind(uv(D$bx1, D$sx1, "OSA univariable"),
            uv(D$bx2, D$sx2, "BMI univariable")))

# ---------- Conditional F (Sanderson 2019, summary-data form) ----------
# Q_j = sum_k [ (bx_jk - sum_{j'} g_jj' bx_j'k) / sx_jk ]^2
# where g solves the weighted regression of exposure j on the other exposures.
# Conditional F_j = (Q_j - (K - 1)) / (K - n_exposures) , K = #instruments
condF <- function(j) {
  y   <- D[[paste0("bx", j)]]
  s   <- D[[paste0("sx", j)]]
  oth <- D[[paste0("bx", 3 - j)]]
  wt  <- 1 / s^2
  # estimate the coefficient of `oth` in the conditional model
  g   <- sum(wt * y * oth) / sum(wt * oth^2)
  Q   <- sum(wt * (y - g * oth)^2)
  K   <- nrow(D)
  (Q - (K - 1)) / (K - 2)
}
F1 <- condF(1); F2 <- condF(2)
cat("\n=== Conditional F statistics (Sanderson 2019) ===\n")
cat("conditional F (OSA | BMI):", round(F1, 2), "\n")
cat("conditional F (BMI | OSA):", round(F2, 2), "\n")
cat("(rule of thumb: conditional F > 10 = adequate instrument strength)\n")

# ---------- Instrument strength (unconditional) ----------
D[, `:=`(F1 = (bx1 / sx1)^2, F2 = (bx2 / sx2)^2)]
cat("\nUnconditional mean F -- OSA:", round(mean(D$F1), 1),
    " | BMI:", round(mean(D$F2), 1), "\n")

fwrite(res, file.path(ROOT, "results/mvmr_osa_bmi_hf.csv"))
cat("\nSaved: mvmr_osa_bmi_hf.csv\n=== done ===\n")
