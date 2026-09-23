# 19_meta_analysis_osa.R
# Fix the weak-instrument problem by META-ANALYSING the two independent OSA GWAS:
#   (a) MVP  -- GCST90475824, 152,031 cases / 278,027 controls, European
#   (b) FinnGen R13 -- G6_SLEEPAPNO_INCLAVO, ~74.7k cases (incl. outpatient)
# A stronger OSA instrument should raise the conditional F in the MVMR.
#
# Meta-analysis: inverse-variance fixed effects on log-OR, aligned by rsID+alleles.

suppressMessages({ library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
res_dir <- RES_DIR

# ---------- (a) MVP significant variants (already extracted) ----------
mvp <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))
mvp[, beta := log(or)]
mvp[, se   := abs(beta) / qnorm(p/2, lower.tail = FALSE)]
mvp <- mvp[is.finite(beta) & is.finite(se) & se > 0]
setnames(mvp, c("rsid","chromosome","bp","effect_allele","other_allele"),
              c("SNP","chr","pos","ea","oa"))
mvp <- mvp[, .(SNP, chr, pos, ea, oa, eaf, beta, se, p)]
cat("MVP significant variants:", nrow(mvp), "\n")

# ---------- (b) FinnGen: stream, keep p < 5e-8 ----------
fgz <- file.path(data_dir, "finnngen_R13_G6_SLEEPAPNO_INCLAVO.gz")
cat("FinnGen file size:", round(file.size(fgz)/1e6,1), "MB\n")

# peek at header first
con <- gzfile(fgz, "rt"); hdr <- readLines(con, n = 1); close(con)
cat("FinnGen header:", hdr, "\n")
cols <- strsplit(hdr, "\t")[[1]]

keepcols <- intersect(c("#chrom","chrom","pos","ref","alt","rsids","rsid",
                        "af_alt","af_alt_cases","af_alt_controls",
                        "beta","sebeta","pval","mlogp","n_hom_cases",
                        "n_hom_ref","n_het","n_case","n_control","odds_ratio"),
                      cols)
cat("available cols:", paste(cols, collapse=", "), "\n")
cat("keeping:", paste(keepcols, collapse=", "), "\n")

cat("\nreading FinnGen (streaming)...\n")
t0 <- Sys.time()
fg <- fread(fgz, select = keepcols, sep = "\t", header = TRUE,
            showProgress = TRUE, nThread = 4)
cat("read", nrow(fg), "rows in",
    round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")

# normalise names
setnames(fg, old = intersect(c("#chrom","chrom"), names(fg)), new = "chr", skip_absent = TRUE)
setnames(fg, old = intersect(c("rsids","rsid"), names(fg)), new = "SNP", skip_absent = TRUE)
setnames(fg, old = intersect(c("pval"), names(fg)), new = "p", skip_absent = TRUE)
setnames(fg, old = intersect(c("sebeta"), names(fg)), new = "se", skip_absent = TRUE)
setnames(fg, old = intersect(c("alt"), names(fg)), new = "ea", skip_absent = TRUE)
setnames(fg, old = intersect(c("ref"), names(fg)), new = "oa", skip_absent = TRUE)
setnames(fg, old = intersect(c("af_alt"), names(fg)), new = "eaf", skip_absent = TRUE)
cat("columns after rename:", paste(names(fg), collapse=", "), "\n")

fg <- fg[!is.na(p) & p < 5e-8]
cat("FinnGen significant:", nrow(fg), "\n")
if ("beta" %in% names(fg)) {
  fg[, beta := as.numeric(beta)]
  fg[, se   := as.numeric(se)]
  cat("using reported beta/se\n")
} else if ("odds_ratio" %in% names(fg)) {
  fg[, beta := log(as.numeric(odds_ratio))]
  fg[, se   := abs(beta)/qnorm(p/2, lower.tail=FALSE)]
  cat("derived beta/se from odds_ratio\n")
}
fg <- fg[is.finite(beta) & is.finite(se) & se > 0]
cat("FinnGen usable:", nrow(fg), "\n")

fwrite(fg, file.path(data_dir, "finngen_osa_p5e8.tsv"), sep="\t")
cat("saved FinnGen significant variants\n")

# ---------- (c) meta-analysis on the union ----------
cat("\n=== meta-analysis (inverse variance, fixed effects) ===\n")
# handle FinnGen multi-rsid fields (comma separated) -> take first
fg2 <- copy(fg)
fg2[, SNP := tstrsplit(SNP, ",", fixed=TRUE)[[1]]]
fg2 <- fg2[!is.na(SNP) & SNP != ""]
cat("FinnGen rows with a usable rsID:", nrow(fg2), "\n")

m <- merge(
  mvp[, .(SNP, ea_m = ea, oa_m = oa, beta_m = beta, se_m = se, p_m = p, eaf_m = eaf)],
  fg2[, .(SNP, ea_f = ea, oa_f = oa, beta_f = beta, se_f = se, p_f = p)],
  by = "SNP")
cat("overlapping significant SNPs:", nrow(m), "\n")

# align FinnGen to MVP alleles
al <- ifelse(m$ea_f == m$ea_m & m$oa_f == m$oa_m,  1,
      ifelse(m$ea_f == m$oa_m & m$oa_f == m$ea_m, -1, NA))
m[, flip := al]
m <- m[!is.na(flip)]
m[, beta_f2 := beta_f * flip]
cat("alignable:", nrow(m), "\n")

if (nrow(m) > 0) {
  m[, w1 := 1/se_m^2][, w2 := 1/se_f2^2]
  m[, beta_meta := (w1*beta_m + w2*beta_f2)/(w1+w2)]
  m[, se_meta   := sqrt(1/(w1+w2))]
  m[, p_meta    := 2*pnorm(-abs(beta_meta/se_meta))]
  m[, z_meta    := beta_meta/se_meta]
  # heterogeneity
  m[, Q := w1*(beta_m-beta_meta)^2 + w2*(beta_f2-beta_meta)^2]
  cat("\nmeta-analysed SNPs:", nrow(m), "\n")
  cat("Cochran Q total:", round(sum(m$Q),1), " df:", nrow(m),
      " p:", signif(pchisq(sum(m$Q), nrow(m), lower.tail=FALSE),3), "\n")
  cat("\nsignificant in meta (p<5e-8):", sum(m$p_meta < 5e-8), "\n")
  cat("median |z|:", round(median(abs(m$z_meta)),2), "\n")
  print(head(m[order(p_meta), .(SNP, beta_m, beta_f2, beta_meta, se_meta, p_meta)], 20))
  fwrite(m, file.path(res_dir, "meta_osa_snps.tsv"), sep="\t")
}

cat("\n=== done ===\n")
