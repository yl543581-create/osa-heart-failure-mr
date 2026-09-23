# 32_steiger_local.R
# Steiger directionality WITHOUT OpenGWAS (which lacked outcome eaf).
# Uses the LOCAL FinnGen files, which carry af_alt.
#
# Method: convert each SNP's log-OR to a liability-scale R2 for both traits on a
# COMMON scale, then test whether the instrument explains more variance in the
# exposure (OSA) than in the outcome (HF).
#
# Liability-scale conversion (Lee et al. 2011 / standard MR practice):
#   For a binary trait with case proportion s, let
#     z = qnorm(1 - s)          (threshold)
#     i = dnorm(z)              (mean liability of cases)
#   Then  beta_liability = beta_logOR * s * (1 - s) / i
#   and   R2_snp = 2 * MAF * (1 - MAF) * beta_liability^2
#
# Both traits are converted on the SAME (liability) scale, so the comparison
# is not confounded by differing case proportions or sample sizes.

suppressMessages({ library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

ivs <- fread(file.path(res_dir, "META_mvmr_input.tsv"))
comp <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))[
          , .(SNP=rsid, chr=chromosome, pos=bp)]
ivs <- merge(ivs, comp, by="SNP", all.x=TRUE); ivs <- ivs[!is.na(pos)]
cat("meta instrument loci:", nrow(ivs), "\n")

cols <- c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt")
rd <- function(f) {
  d <- fread(f, select=cols, sep="\t", header=TRUE, nThread=4, showProgress=FALSE)
  setnames(d, c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
              c("chr","pos","oa","ea","SNP","p","b","se","eaf"))
  d[, SNP := tstrsplit(SNP, ",", fixed=TRUE)[[1]]]
  d[!is.na(SNP) & !is.na(b) & !is.na(se) & se>0 & !is.na(eaf) & eaf>0 & eaf<1]
}
cat("reading FinnGen OSA...\n")
osa <- rd(file.path(data_dir, "finnngen_R13_G6_SLEEPAPNO_INCLAVO.gz"))
cat("reading FinnGen HF...\n")
hf  <- rd(file.path(data_dir, "finngen_R13_I9_HEARTFAIL.gz"))

# restrict to instrument loci (+/- 1 Mb) to speed the merge
setkey(osa, SNP); setkey(hf, SNP)
want <- unique(c(ivs$SNP))
# also include any SNP in the loci so we get several SNPs per locus where available
sel <- c()
for (i in seq_len(nrow(ivs))) {
  lo <- ivs$pos[i]-1e6; hi <- ivs$pos[i]+1e6
  sel <- c(sel, osa[chr==ivs$chr[i] & pos>=lo & pos<=hi, SNP])
}
sel <- unique(sel[!is.na(sel)])
cat("SNPs in instrument loci:", length(sel), "\n")

m <- merge(osa[SNP %in% sel, .(SNP, ea, oa, eaf, b_osa=b, se_osa=se, p_osa=p)],
           hf [SNP %in% sel, .(SNP, ea2=ea, oa2=oa, eaf2=eaf, b_hf=b, se_hf=se, p_hf=p)],
           by="SNP")
m[, flip := ifelse(ea==ea2 & oa==oa2, 1, ifelse(ea==oa2 & oa==ea2, -1, NA_integer_))]
m <- m[!is.na(flip)]
m[, b_hf2 := b_hf * flip]
m <- m[order(p_osa)][!duplicated(SNP)]
cat("overlapping SNPs with alleles+MAF:", nrow(m), "\n")

# ---- liability-scale conversion ----
liab <- function(beta, maf, s) {
  z <- qnorm(1 - s); i <- dnorm(z)
  bl <- beta * s * (1 - s) / i
  2 * maf * (1 - maf) * bl^2
}
S_OSA <- 74697/(74697+430000)
S_HF  <- 41591/(41591+458595)
cat(sprintf("case proportions -- OSA: %.3f  HF: %.3f\n", S_OSA, S_HF))

m[, r2_osa := liab(b_osa, pmin(eaf,1-eaf), S_OSA)]
m[, r2_hf  := liab(b_hf2, pmin(eaf2,1-eaf2), S_HF)]
m <- m[is.finite(r2_osa) & is.finite(r2_hf)]

cat("\n=== Steiger (liability scale, same cohort) ===\n")
cat("sum R2 OSA :", sprintf("%.6f", sum(m$r2_osa)), "\n")
cat("sum R2 HF  :", sprintf("%.6f", sum(m$r2_hf)),  "\n")
cat("ratio OSA/HF:", sprintf("%.2f", sum(m$r2_osa)/sum(m$r2_hf)), "\n")
m[, correct := r2_osa > r2_hf]
cat("\ncorrect direction:", sum(m$correct), "/", nrow(m), "\n")
if (nrow(m) >= 5) {
  bt <- binom.test(sum(m$correct), nrow(m), 0.5)
  cat("binom p =", signif(bt$p.value,4), "\n")
}

# also restricted to the 10 independent instrument SNPs
m10 <- m[SNP %in% ivs$SNP]
cat("\n--- restricted to the", nrow(m10), "independent instrument SNPs ---\n")
if (nrow(m10) > 0) {
  cat("sum R2 OSA:", sprintf("%.6f", sum(m10$r2_osa)),
      "| sum R2 HF:", sprintf("%.6f", sum(m10$r2_hf)), "\n")
  cat("correct:", sum(m10$r2_osa > m10$r2_hf), "/", nrow(m10), "\n")
  print(m10[, .(SNP, r2_osa=signif(r2_osa,6), r2_hf=signif(r2_hf,6),
                correct = r2_osa > r2_hf)])
}
fwrite(m[, .(SNP, ea, oa, eaf, b_osa, b_hf=b_hf2, r2_osa, r2_hf, correct)],
       file.path(res_dir,"STEIGER_local.csv"))
cat("\nsaved: results/STEIGER_local.csv\n=== done ===\n")
