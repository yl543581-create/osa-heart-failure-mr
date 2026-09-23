# 25_coloc.R
# Colocalisation (coloc.abf) at each independent OSA instrument locus.
#
# QUESTION: is the OSA signal and the heart-failure signal driven by the SAME
# causal variant (PP.H4 high), or by distinct variants in LD (PP.H3)?
#
# Given our finding that OSA->HF is adiposity-mediated, we EXPECT low PP.H4 and
# relatively high PP.H3 at obesity loci -- which is itself the informative result.

suppressMessages({ library(data.table); library(coloc); library(ieugwasr); library(TwoSampleMR) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
hf_id   <- "ebi-a-GCST009541"

ivs <- fread(file.path(res_dir, "META_mvmr_input.tsv"))
cat("loci from MVMR input:", nrow(ivs), "\n")

# need position for regional extraction; get from the meta file + MVP extraction
comp <- fread(file.path(ROOT, "data/GCST90475824_p5e8.tsv"))[
          , .(SNP = rsid, chr = chromosome, pos = bp)]
ivs <- merge(ivs, comp, by = "SNP", all.x = TRUE)
ivs <- ivs[!is.na(pos)]
cat("with positions:", nrow(ivs), "\n")

# sample sizes
N_OSA_EUR <- 152031 + 278027                       # MVP European
N_HF      <- 47309 + 930014                        # HERMES

window <- 5e5    # +/- 500 kb

run_coloc <- function(chr, pos, snp) {
  lo <- pos - window; hi <- pos + window
  cat("\n=== locus", snp, "chr", chr, sprintf("%.2f-%.2f Mb", lo/1e6, hi/1e6), "===\n")
  reg <- sprintf("%s:%d-%d", chr, lo, hi)

  get <- function(id) {
    r <- tryCatch(ieugwasr::associations(variants = reg, id = id), error = function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    d <- as.data.table(r)
    if (!all(c("rsid","beta","se","ea","nea","eaf") %in% names(d))) return(NULL)
    d <- d[!is.na(beta) & !is.na(se) & se > 0 & !is.na(eaf)]
    d[!duplicated(rsid)]
  }
  a <- get("finn-b-G6_SLEEPAPNO")
  b <- get(hf_id)
  if (is.null(a) || is.null(b)) { cat("  no data\n"); return(NULL) }
  m <- merge(a[, .(SNP=rsid, ea, nea, eaf, beta, se)],
             b[, .(SNP=rsid, ea2=ea, nea2=nea, eaf2=eaf, beta2=beta, se2=se)],
             by = "SNP")
  cat("  overlapping SNPs:", nrow(m), "\n")
  if (nrow(m) < 30) { cat("  too few SNPs for coloc\n"); return(NULL) }

  # align to a common effect allele
  flip <- ifelse(m$ea == m$ea2 & m$nea == m$nea2,  1,
          ifelse(m$ea == m$nea2 & m$nea == m$ea2, -1, NA))
  m[, flip := flip]
  m <- m[!is.na(flip)]
  m[, beta2a := beta2 * flip]

  d1 <- list(snp = m$SNP, beta = m$beta, varbeta = m$se^2,
             type = "cc", N = N_OSA_EUR, MAF = pmin(m$eaf, 1-m$eaf),
             s = 152031/(152031+278027))
  d2 <- list(snp = m$SNP, beta = m$beta2a, varbeta = m$se2^2,
             type = "cc", N = N_HF, MAF = pmin(m$eaf2, 1-m$eaf2),
             s = 47309/(47309+930014))
  res <- tryCatch(coloc.abf(d1, d2), error = function(e) { cat("  coloc error:", conditionMessage(e), "\n"); NULL })
  if (is.null(res)) return(NULL)
  pp <- res$summary
  cat(sprintf("  PP.H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f\n",
              pp[["PP.H0.abf"]], pp[["PP.H1.abf"]], pp[["PP.H2.abf"]],
              pp[["PP.H3.abf"]], pp[["PP.H4.abf"]]))
  data.table(SNP = snp, chr = chr, pos = pos, nsnp = nrow(m),
             PP.H0 = pp[["PP.H0.abf"]], PP.H1 = pp[["PP.H1.abf"]],
             PP.H2 = pp[["PP.H2.abf"]], PP.H3 = pp[["PP.H3.abf"]],
             PP.H4 = pp[["PP.H4.abf"]])
}

out <- list()
for (i in seq_len(nrow(ivs))) {
  r <- tryCatch(run_coloc(ivs$chr[i], ivs$pos[i], ivs$SNP[i]),
                error = function(e) { cat("  FAIL:", conditionMessage(e), "\n"); NULL })
  if (!is.null(r)) out[[length(out)+1]] <- r
}

if (length(out) > 0) {
  R <- rbindlist(out)
  cat("\n\n================ COLOC SUMMARY ================\n")
  print(R[, .(SNP, nsnp, PP.H3 = round(PP.H3,3), PP.H4 = round(PP.H4,3))])
  fwrite(R, file.path(res_dir, "COLOC_results.csv"))
  cat("\nsaved:", file.path(res_dir,"COLOC_results.csv"), "\n")
  cat("\nloci with PP.H4 >= 0.8:", sum(R$PP.H4 >= 0.8), "of", nrow(R), "\n")
} else cat("\nno loci successfully colocalised\n")
cat("=== done ===\n")
