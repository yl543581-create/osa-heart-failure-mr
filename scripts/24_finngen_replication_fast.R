# 24_finngen_replication_fast.R
# Fast replication: read only the needed columns from the FinnGen HF file with fread.

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

fgz <- file.path(data_dir, "finngen_R13_I9_HEARTFAIL.gz")
cat("\nreading FinnGen HF (selected columns only)...\n")
t0 <- Sys.time()
fg <- fread(fgz, select = c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
            sep="\t", header=TRUE, showProgress=TRUE, nThread=4)
cat("read", nrow(fg), "rows in",
    round(as.numeric(difftime(Sys.time(), t0, units="mins")),1), "min\n")
setnames(fg, c("#chrom","rsids","pval","beta","sebeta","af_alt"),
              c("chr","SNP","pval.outcome","beta.outcome","se.outcome","eaf.outcome"))
fg[, SNP := tstrsplit(SNP, ",", fixed=TRUE)[[1]]]
fg <- fg[!is.na(SNP)]

ou <- fg[SNP %in% ivs$SNP]
ou <- ou[!duplicated(SNP)]
setnames(ou, c("ref","alt"), c("other_allele.outcome","effect_allele.outcome"))
ou <- ou[!is.na(beta.outcome) & !is.na(se.outcome) & se.outcome > 0]
# harmonise_data requires outcome metadata columns
ou[, id.outcome := "FinnGen_R13_I9_HEARTFAIL"]
ou[, outcome := "Heart failure (FinnGen R13)"]
ou[, mr_keep.outcome := TRUE]
cat("FinnGen HF instrument SNPs recovered:", nrow(ou), "of", nrow(ivs), "\n")

if (nrow(ou) < 2) { cat("too few - aborting\n"); quit(save="no") }
fwrite(ou, file.path(res_dir,"finngen_hf_instrument_snps.tsv"), sep="\t")

ex <- format_data(as.data.frame(ivs), type="exposure",
        snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
        effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")
dat <- harmonise_data(ex, as.data.frame(ou), action=2)
cat("harmonised:", nrow(dat), " kept:", sum(dat$mr_keep), "\n")
if (sum(dat$mr_keep) >= 2) {
  dat2 <- dat[dat$mr_keep, ]
  r <- generate_odds_ratios(mr(dat2, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median")))
  cat("\n=== REPLICATION: meta-OSA -> FinnGen heart failure ===\n")
  for (i in seq_len(nrow(r)))
    cat(sprintf("%-22s nsnp=%3d OR=%.3f (%.3f-%.3f) p=%.3g\n",
        r$method[i], r$nsnp[i], r$or[i], r$or_lci95[i], r$or_uci95[i], r$pval[i]))
  fwrite(as.data.table(r), file.path(res_dir,"REPLICATION_finngen_HF.csv"))
} else {
  cat("too few kept; showing harmonisation:\n")
  print(dat[, c("SNP","effect_allele.exposure","other_allele.exposure",
                "effect_allele.outcome","other_allele.outcome","mr_keep")])
}
cat("\n=== done ===\n")
