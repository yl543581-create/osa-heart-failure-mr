# 20_meta_analysis_v2.R
# Fix the data.table scoping bug and add an essential validation:
# do MVP and FinnGen actually AGREE in direction? (they are very different cohorts)

suppressMessages({ library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
res_dir <- RES_DIR

mvp <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))
mvp[, beta := log(or)]
mvp[, se   := abs(beta) / qnorm(p/2, lower.tail = FALSE)]
mvp <- mvp[is.finite(beta) & is.finite(se) & se > 0]
setnames(mvp, c("rsid","bp","effect_allele","other_allele"),
              c("SNP","pos","ea","oa"))
mvp <- mvp[, .(SNP, pos, ea, oa, eaf, beta, se, p)]
cat("MVP significant:", nrow(mvp), "\n")

fg <- fread(file.path(data_dir, "finngen_osa_p5e8.tsv"))
fg[, SNP := tstrsplit(SNP, ",", fixed = TRUE)[[1]]]
fg <- fg[!is.na(SNP) & SNP != ""]
fg <- fg[is.finite(beta) & is.finite(se) & se > 0]
cat("FinnGen significant:", nrow(fg), "\n")

m <- merge(
  mvp[, .(SNP, ea_m = ea, oa_m = oa, beta_m = beta, se_m = se, p_m = p, eaf_m = eaf)],
  fg[,  .(SNP, ea_f = ea, oa_f = oa, beta_f = beta, se_f = se, p_f = p)],
  by = "SNP")
cat("overlapping significant SNPs:", nrow(m), "\n")

m[, flip := ifelse(ea_f == ea_m & oa_f == oa_m,  1,
            ifelse(ea_f == oa_m & oa_f == ea_m, -1, NA_integer_))]
m <- m[!is.na(flip)]
m[, beta_f2 := beta_f * flip]
cat("alignable:", nrow(m), "\n")

# ---- VALIDATION: do the two cohorts agree? ----
cat("\n=== VALIDATION: MVP vs FinnGen agreement ===\n")
cat("correlation of beta (aligned):",
    round(cor(m$beta_m, m$beta_f2), 3), "\n")
cat("same sign:", sum(sign(m$beta_m) == sign(m$beta_f2)), "/", nrow(m),
    sprintf("(%.1f%%)\n", 100*mean(sign(m$beta_m) == sign(m$beta_f2))))
# compare z-scores (independent of scale differences)
m[, z_m := beta_m/se_m][, z_f := beta_f2/se_f]
cat("correlation of z-scores:", round(cor(m$z_m, m$z_f), 3), "\n")
# regression slope of FinnGen on MVP
sl <- coef(lm(beta_f2 ~ beta_m, data = m))
cat("slope (FinnGen ~ MVP):", round(sl[2], 3), "\n")

# meta-analysis
m[, w1 := 1/se_m^2]
m[, w2 := 1/se_f2^2]
m[, beta_meta := (w1*beta_m + w2*beta_f2)/(w1+w2)]
m[, se_meta   := sqrt(1/(w1+w2))]
m[, p_meta    := 2*pnorm(-abs(beta_meta/se_meta))]
m[, Q := w1*(beta_m-beta_meta)^2 + w2*(beta_f2-beta_meta)^2]

cat("\n=== META-ANALYSIS (IVW fixed effects) ===\n")
cat("SNPs meta-analysed:", nrow(m), "\n")
Qtot <- sum(m$Q); df <- nrow(m)
cat(sprintf("Cochran Q = %.1f (df=%d), p = %.3g\n", Qtot, df,
            pchisq(Qtot, df, lower.tail = FALSE)))
cat("I2 approx =", sprintf("%.1f%%\n", max(0, (Qtot-df)/Qtot*100)))
cat("median |z_meta| =", round(median(abs(m$beta_meta/m$se_meta)), 2), "\n")
cat("meta p<5e-8:", sum(m$p_meta < 5e-8), "\n")

setorder(m, p_meta)
cat("\n=== top 25 meta signals ===\n")
print(m[1:min(25,.N), .(SNP, beta_m, beta_f2, beta_meta, se_meta, p_meta)])

fwrite(m, file.path(res_dir, "meta_osa_snps.tsv"), sep = "\t")
cat("\nsaved:", file.path(res_dir, "meta_osa_snps.tsv"), "\n")
cat("=== done ===\n")
