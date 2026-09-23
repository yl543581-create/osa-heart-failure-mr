# 23_finngen_replication.R
# Independent replication in FinnGen R13 heart failure (NOT part of HERMES).
# Single streaming pass over the local FinnGen HF file.

suppressMessages({ library(data.table); library(TwoSampleMR); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR
res_dir <- RES_DIR

meta <- fread(file.path(res_dir, "meta_osa_snps.tsv"))
mv   <- fread(file.path(data_dir, "GCST90475824_p5e8.tsv"))[
          , .(SNP=rsid, effect_allele, other_allele, eaf)]
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP, pval=meta$p_meta),
                                  clump_r2=0.001, clump_kb=10000, pop="EUR"),
               error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
ivs <- merge(meta[, .(SNP, beta=beta_meta, se=se_meta, p=p_meta)], mv, by="SNP")
cat("instrument SNPs:", nrow(ivs), "\n")
print(ivs)

fgz <- file.path(data_dir, "finngen_R13_I9_HEARTFAIL.gz")
want <- ivs$SNP

cat("\nsingle streaming pass over FinnGen HF...\n")
con <- gzfile(fgz, "rt")
hdr <- strsplit(readLines(con, n=1), "\t")[[1]]
cat("header:", paste(hdr, collapse=" | "), "\n")
ci <- function(nm) { i <- which(hdr == nm); if(length(i)) i[1] else NA_integer_ }
rs_c  <- ci("rsids"); p_c <- ci("pval"); b_c <- ci("beta"); s_c <- ci("sebeta")
ref_c <- ci("ref");   alt_c <- ci("alt")

hits <- list(); n <- 0L
repeat {
  ln <- readLines(con, n = 200000L)
  if (length(ln) == 0) break
  n <- n + length(ln)
  f <- strsplit(ln, "\t", fixed = TRUE)
  rs <- vapply(f, function(x) if (length(x) >= rs_c) x[rs_c] else NA_character_, "")
  sel <- which(rs %in% want)
  if (length(sel) > 0) {
    hits[[length(hits)+1]] <- data.table(
      rsids = rs[sel],
      ref   = vapply(f[sel], function(x) x[ref_c], ""),
      alt   = vapply(f[sel], function(x) x[alt_c], ""),
      pval  = as.numeric(vapply(f[sel], function(x) x[p_c], "")),
      beta  = as.numeric(vapply(f[sel], function(x) x[b_c], "")),
      sebeta= as.numeric(vapply(f[sel], function(x) x[s_c], "")))
  }
}
close(con)
cat("rows scanned:", n, "\n")

if (length(hits) == 0) stop("no instrument SNPs found in FinnGen HF")
ou <- rbindlist(hits)
ou[, rsids := tstrsplit(rsids, ",", fixed=TRUE)[[1]]]
ou <- ou[!duplicated(rsids)]
ou <- ou[!is.na(beta) & !is.na(sebeta) & sebeta > 0]
setnames(ou, c("rsids","ref","alt","pval","beta","sebeta"),
             c("SNP","other_allele.outcome","effect_allele.outcome",
               "pval.outcome","beta.outcome","se.outcome"))
cat("FinnGen HF instrument SNPs recovered:", nrow(ou), "\n")
fwrite(ou, file.path(res_dir,"finngen_hf_instrument_snps.tsv"), sep="\t")

# ---------- replication MR ----------
ex <- format_data(as.data.frame(ivs), type="exposure",
        snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
        effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")
dat <- harmonise_data(ex, as.data.frame(ou), action=2)
cat("\nharmonised rows:", nrow(dat), " | kept:", sum(dat$mr_keep), "\n")
dat <- dat[dat$mr_keep, ]
if (nrow(dat) < 2) {
  cat("\n*** too few SNPs survived harmonisation - inspecting ***\n")
  print(dat[, c("SNP","effect_allele.exposure","other_allele.exposure",
                "effect_allele.outcome","other_allele.outcome","mr_keep")])
} else {
  r <- generate_odds_ratios(mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")))
  cat("\n=== REPLICATION: meta-OSA -> FinnGen heart failure ===\n")
  for (i in seq_len(nrow(r)))
    cat(sprintf("%-22s nsnp=%3d OR=%.3f (%.3f-%.3f) p=%.3g\n",
        r$method[i], r$nsnp[i], r$or[i], r$or_lci95[i], r$or_uci95[i], r$pval[i]))
  fwrite(as.data.table(r), file.path(res_dir,"REPLICATION_finngen_HF.csv"))
}
cat("\n=== done ===\n")
