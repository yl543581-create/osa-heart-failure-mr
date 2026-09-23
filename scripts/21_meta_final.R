# 21_meta_final.R
# Meta-analysis, with all arithmetic precomputed as plain vectors
# (avoids data.table's column-creation ordering issue).

suppressMessages({ library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
res_dir <- RES_DIR

mvp <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))
mvp[, beta := log(or)]
mvp[, se   := abs(beta) / qnorm(p/2, lower.tail = FALSE)]
mvp <- mvp[is.finite(beta) & is.finite(se) & se > 0]
setnames(mvp, c("rsid","bp","effect_allele","other_allele"), c("SNP","pos","ea","oa"))
mvp <- mvp[, .(SNP, pos, ea, oa, eaf, beta, se, p)]

fg <- fread(file.path(data_dir, "finngen_osa_p5e8.tsv"))
fg[, SNP := tstrsplit(SNP, ",", fixed = TRUE)[[1]]]
fg <- fg[!is.na(SNP) & SNP != ""]
fg <- fg[is.finite(beta) & is.finite(se) & se > 0]

m <- merge(
  mvp[, .(SNP, ea_m=ea, oa_m=oa, beta_m=beta, se_m=se, p_m=p, eaf_m=eaf)],
  fg[,  .(SNP, ea_f=ea, oa_f=oa, beta_f=beta, se_f=se, p_f=p)],
  by = "SNP")
m <- m[!is.na(SNP)]
m[, flip := ifelse(ea_f==ea_m & oa_f==oa_m, 1L,
            ifelse(ea_f==oa_m & oa_f==ea_m, -1L, NA_integer_))]
m <- m[!is.na(flip)]
cat("alignable overlapping SNPs:", nrow(m), "\n")

# ---- plain vectors ----
bm <- m$beta_m;  sm <- m$se_m
bf <- m$beta_f * m$flip
sf <- m$se_f

cat("\n=== VALIDATION ===\n")
cat("cor(beta):", round(cor(bm, bf), 4), "\n")
cat("same sign:", sum(sign(bm)==sign(bf)), "/", length(bm), "\n")
cat("cor(z):", round(cor(bm/sm, bf/sf), 4), "\n")

w1 <- 1/sm^2
w2 <- 1/sf^2
bmeta <- (w1*bm + w2*bf)/(w1+w2)
smeta <- sqrt(1/(w1+w2))
pmeta <- 2*pnorm(-abs(bmeta/smeta))
Q <- w1*(bm-bmeta)^2 + w2*(bf-bmeta)^2

Out <- data.table(SNP=m$SNP, beta_mvp=bm, se_mvp=sm, beta_fg=bf, se_fg=sf,
                  beta_meta=bmeta, se_meta=smeta, p_meta=pmeta, Q=Q)
setorder(Out, p_meta)

cat("\n=== META-ANALYSIS (IVW fixed effects) ===\n")
cat("SNPs:", nrow(Out), "\n")
Qtot <- sum(Q); df <- length(Q)
cat(sprintf("Cochran Q = %.1f (df=%d), p = %.3g\n", Qtot, df, pchisq(Qtot, df, lower.tail=FALSE)))
cat(sprintf("I2 = %.1f%%\n", max(0,(Qtot-df)/Qtot*100)))
cat("meta p<5e-8:", sum(Out$p_meta < 5e-8), "of", nrow(Out), "\n")
cat("median |z|:", round(median(abs(bmeta/smeta)),2), "\n")
cat("\n=== top 25 meta signals ===\n")
print(Out[1:min(25,.N), .(SNP, beta_mvp, beta_fg, beta_meta, se_meta, p_meta)])
fwrite(Out, file.path(res_dir,"meta_osa_snps.tsv"), sep="\t")
cat("\nsaved\n=== done ===\n")
