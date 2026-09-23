# 27_coloc_v2.R
# Fix: ieugwasr::associations with a "chr:start-end" string returned 0 SNPs.
# Instead, fetch the full significant-variant lists we already have locally
# (MVP OSA p<5e-8, FinnGen OSA p<5e-8) and query OpenGWAS for those exact rsIDs.

suppressMessages({ library(data.table); library(coloc); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
hf_id <- "ebi-a-GCST009541"

N_OSA_CASE <- 152031; N_OSA_CTRL <- 278027
N_HF_CASE  <- 47309;  N_HF_CTRL  <- 930014

ivs <- fread(file.path(res_dir, "META_mvmr_input.tsv"))
comp <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))[
          , .(SNP = rsid, chr = chromosome, pos = bp)]
ivs <- merge(ivs, comp, by = "SNP", all.x = TRUE)
ivs <- ivs[!is.na(pos)]
cat("loci:", nrow(ivs), "\n")

# local OSA significant variants (for regional SNP sets)
osa_local <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))[
   , .(SNP = rsid, chr = chromosome, pos = bp, ea = effect_allele,
       oa = other_allele, eaf, beta = log(or), p)]
setkey(osa_local, chr, pos)
cat("local OSA significant variants:", nrow(osa_local), "\n")

run_coloc <- function(chr, pos, snp) {
  lo <- pos - 5e5; hi <- pos + 5e5
  cat(sprintf("\n=== %s (chr%d, %.2f-%.2f Mb) ===\n", snp, chr, lo/1e6, hi/1e6))
  reg <- osa_local[chr == get("chr") & pos >= lo & pos <= hi]
  cat("  local OSA SNPs in region:", nrow(reg), "\n")
  if (nrow(reg) < 20) { cat("  too few local SNPs\n"); return(NULL) }

  # get OSA + HF effects for these SNPs from OpenGWAS
  get <- function(id) {
    r <- tryCatch(ieugwasr::associations(variants = reg$SNP, id = id), error=function(e) NULL)
    if (is.null(r) || nrow(r) == 0) return(NULL)
    d <- as.data.table(r)
    d <- d[!is.na(beta) & !is.na(se) & se > 0]
    d[!duplicated(rsid)]
  }
  a <- get("finn-b-G6_SLEEPAPNO")
  b <- get(hf_id)
  if (is.null(a) || is.null(b)) { cat("  missing data\n"); return(NULL) }
  cat("  OSA SNPs:", nrow(a), " HF SNPs:", nrow(b), "\n")

  m <- merge(a[, .(SNP=rsid, ea, nea, eaf, beta, se)],
             b[, .(SNP=rsid, ea2=ea, nea2=nea, eaf2=eaf, beta2=beta, se2=se)],
             by = "SNP")
  cat("  overlapping:", nrow(m), "\n")
  if (nrow(m) < 30) { cat("  too few overlapping\n"); return(NULL) }

  m[, flip := ifelse(ea==ea2 & nea==nea2, 1, ifelse(ea==nea2 & nea==ea2, -1, NA_integer_))]
  m <- m[!is.na(flip)]
  m[, beta2a := beta2 * flip]

  d1 <- list(snp=m$SNP, beta=m$beta, varbeta=m$se^2, type="cc",
             N=N_OSA_CASE+N_OSA_CTRL, MAF=pmin(m$eaf,1-m$eaf),
             s=N_OSA_CASE/(N_OSA_CASE+N_OSA_CTRL))
  d2 <- list(snp=m$SNP, beta=m$beta2a, varbeta=m$se2^2, type="cc",
             N=N_HF_CASE+N_HF_CTRL, MAF=pmin(m$eaf2,1-m$eaf2),
             s=N_HF_CASE/(N_HF_CASE+N_HF_CTRL))
  res <- tryCatch(coloc.abf(d1, d2), error=function(e){cat("  coloc err:",conditionMessage(e),"\n");NULL})
  if (is.null(res)) return(NULL)
  pp <- res$summary
  cat(sprintf("  PP.H0=%.3f H1=%.3f H2=%.3f H3=%.3f H4=%.3f\n",
      pp[["PP.H0.abf"]],pp[["PP.H1.abf"]],pp[["PP.H2.abf"]],pp[["PP.H3.abf"]],pp[["PP.H4.abf"]]))
  data.table(SNP=snp, chr=chr, pos=pos, nsnp=nrow(m),
             PP.H0=pp[["PP.H0.abf"]], PP.H1=pp[["PP.H1.abf"]],
             PP.H2=pp[["PP.H2.abf"]], PP.H3=pp[["PP.H3.abf"]], PP.H4=pp[["PP.H4.abf"]])
}

out <- list()
for (i in seq_len(nrow(ivs))) {
  r <- tryCatch(run_coloc(ivs$chr[i], ivs$pos[i], ivs$SNP[i]),
                error=function(e){cat("  FAIL:",conditionMessage(e),"\n");NULL})
  if (!is.null(r)) out[[length(out)+1]] <- r
}
if (length(out) > 0) {
  R <- rbindlist(out)
  cat("\n\n============ COLOC SUMMARY ============\n")
  print(R[, .(SNP, nsnp, PP.H3=round(PP.H3,3), PP.H4=round(PP.H4,3))])
  fwrite(R, file.path(res_dir,"COLOC_results.csv"))
  cat("\nPP.H4 >= 0.8:", sum(R$PP.H4>=0.8), "of", nrow(R), "\n")
  cat("PP.H3 >= 0.8:", sum(R$PP.H3>=0.8), "of", nrow(R), "\n")
} else cat("\nnone colocalised\n")
cat("=== done ===\n")
