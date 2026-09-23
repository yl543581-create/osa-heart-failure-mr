# 22_mvmr_meta.R
# FINAL MVMR with a meta-analysed (stronger) OSA instrument.
#
# Logic: the weak conditional F(OSA|BMI)=1.3 reflects genuine genetic overlap
# between OSA and BMI. A stronger instrument is the only principled fix.
# Here we compare, on the SAME SNPs:
#   Instrument A: MVP-only OSA effects
#   Instrument B: MVP+FinnGen meta-analysed OSA effects (smaller SE -> stronger)
# If conditional F does not improve, the weakness is biological, not statistical.

suppressMessages({ library(TwoSampleMR); library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
hf_id <- "ebi-a-GCST009541"; bmi_id <- "ieu-b-40"

meta <- fread(file.path(res_dir, "meta_osa_snps.tsv"))
cat("meta SNPs:", nrow(meta), "\n")

# ---------- clump the meta SNPs into independent instruments ----------
cl <- tryCatch(ieugwasr::ld_clump(
        data.frame(rsid = meta$SNP, pval = meta$p_meta),
        clump_r2 = 0.001, clump_kb = 10000, pop = "EUR"),
      error = function(e) { cat("clump failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(cl) && nrow(cl) > 0) {
  cat("independent meta instruments:", nrow(cl), "\n")
  meta <- meta[SNP %in% cl$rsid]
} else cat("using all meta SNPs\n")
cat("final meta instrument SNPs:", nrow(meta), "\n")

# ================= Part 1: meta OSA -> HF (univariable) =================
cat("\n########## PART 1: meta OSA -> Heart failure ##########\n")
ivs_meta <- meta[, .(SNP, beta = beta_meta, se = se_meta, p = p_meta)]
# need alleles: take from MVP extraction
mv <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))
mv <- mv[, .(SNP = rsid, effect_allele, other_allele, eaf)]
ivs_meta <- merge(ivs_meta, mv, by = "SNP")
cat("with alleles:", nrow(ivs_meta), "\n")

fmt <- function(dt) format_data(as.data.frame(dt), type="exposure",
  snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
  effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")

outH <- extract_outcome_data(snps = ivs_meta$SNP, outcomes = hf_id)
dH <- harmonise_data(fmt(ivs_meta), outH, action = 2); dH <- dH[dH$mr_keep, ]
rH <- generate_odds_ratios(mr(dH, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")))
for (i in seq_len(nrow(rH)))
  cat(sprintf("%-22s nsnp=%3d OR=%.3f (%.3f-%.3f) p=%.3g\n",
      rH$method[i], rH$nsnp[i], rH$or[i], rH$or_lci95[i], rH$or_uci95[i], rH$pval[i]))
cat("mean F =", round(mean((ivs_meta$beta/ivs_meta$se)^2),1), "\n")
fwrite(as.data.table(rH), file.path(res_dir,"META_OSA_HF.csv"))

# ================= Part 2: MVMR comparison =================
cat("\n########## PART 2: MVMR on identical SNPs ##########\n")

ga <- function(id) {
  r <- tryCatch(ieugwasr::associations(variants = meta$SNP, id = id), error=function(e) NULL)
  if (is.null(r)) return(NULL)
  d <- as.data.table(r)
  d <- d[!is.na(beta) & !is.na(se) & !is.na(ea) & !is.na(nea)]
  d <- d[!duplicated(rsid)]
  d[, .(SNP = rsid, ea, oa = nea, beta, se)]
}
a_osa <- ga("finn-b-G6_SLEEPAPNO")
a_bmi <- ga(bmi_id)
outA  <- as.data.table(outH)[, .(SNP, ea_o=effect_allele.outcome,
                                oa_o=other_allele.outcome, by=beta.outcome, sy=se.outcome)]
outA <- outA[!duplicated(SNP)]
cat("assoc rows  OSA:", ifelse(is.null(a_osa),0,nrow(a_osa)),
    " BMI:", ifelse(is.null(a_bmi),0,nrow(a_bmi)),
    " HF:", nrow(outA), "\n")

m <- merge(merge(a_osa, a_bmi, by="SNP", suffixes=c(".o",".b")), outA, by="SNP")
m <- merge(m, meta[, .(SNP, beta_meta, se_meta)], by="SNP")
cat("merged:", nrow(m), "\n")

al <- function(e1,o1,e2,o2v) ifelse(e1==e2 & o1==o2v, 1L, ifelse(e1==o2v & o1==e2, -1L, NA_integer_))
f1 <- al(m$ea.o, m$oa.o, m$ea_o, m$oa_o)     # FinnGen-OSA  -> outcome
f2 <- al(m$ea.b, m$oa.b, m$ea_o, m$oa_o)     # BMI          -> outcome
fm <- al(m$ea.o, m$oa.o, m$ea_o, m$oa_o)     # meta-OSA alleles vs outcome (same ea as FinnGen OSA? check)
# meta effects are aligned to MVP alleles; MVP ea/oa may differ from FinnGen's.
# Use the MVP allele column we merged earlier instead:
m2 <- merge(m, mv, by="SNP")
fm2 <- al(m2$effect_allele, m2$other_allele, m2$ea_o, m2$oa_o)
cat("alignable (osa_fg / bmi / meta):", sum(!is.na(f1)), sum(!is.na(f2)), sum(!is.na(fm2)), "\n")

k <- !is.na(f1) & !is.na(f2) & !is.na(fm2)
D <- data.table(
  SNP  = m$SNP[k],
  bx1_fg   = m$beta.o[k]*f1[k],   sx1_fg = m$se.o[k],
  bx1_meta = m2$beta_meta[k]*fm2[k], sx1_meta = m2$se_meta[k],
  bx2_bmi  = m$beta.b[k]*f2[k],   sx2    = m$se.b[k],
  by       = m$by[k],             sy     = m$sy[k]
)
D <- D[is.finite(bx1_fg)&is.finite(bx1_meta)&is.finite(bx2_bmi)&is.finite(by)&sy>0]
cat("final MVMR SNPs:", nrow(D), "\n")

run_mvmr <- function(bx1, sx1, label) {
  w <- 1/D$sy^2
  X <- cbind(1, bx1, D$bx2_bmi)
  XtWX <- t(X) %*% (w*X)
  f <- tryCatch(solve(XtWX), error=function(e) NULL)
  if (is.null(f)) { cat(label, ": singular\n"); return(invisible()) }
  fit <- f %*% (t(X) %*% (w*D$by)); se <- sqrt(diag(f))
  cat("\n---", label, "---\n")
  for (j in 1:2) {
    b <- fit[j+1]; s <- se[j+1]
    cat(sprintf("  %-13s OR=%.3f (%.3f-%.3f) p=%.3g\n",
        c("OSA (direct)","BMI (direct)")[j], exp(b),
        exp(b-1.96*s), exp(b+1.96*s), 2*pnorm(-abs(b/s))))
  }
  # Sanderson conditional F_j
  K <- nrow(D)
  cf <- function(y, s, o) { wj <- 1/s^2; g <- sum(wj*y*o)/sum(wj*o^2)
                            Q <- sum(wj*(y-g*o)^2); (Q-(K-1))/(K-2) }
  cat(sprintf("  conditional F(OSA|BMI) = %.1f\n", cf(bx1, sx1, D$bx2_bmi)))
  cat(sprintf("  conditional F(BMI|OSA) = %.1f\n", cf(D$bx2_bmi, D$sx2, bx1)))
}

run_mvmr(D$bx1_fg,   D$sx1_fg,   "Instrument A: FinnGen-only OSA")
run_mvmr(D$bx1_meta, D$sx1_meta, "Instrument B: MVP+FinnGen META OSA")

fwrite(D, file.path(res_dir,"META_mvmr_input.tsv"), sep="\t")
cat("\n=== done ===\n")
