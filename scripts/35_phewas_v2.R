# 35_phewas_v2.R
# Bug fix: the previous version overwrote the instrument set with a FinnGen-only
# extraction (5 SNPs). Now the exposure is built once, properly, from the
# meta-analysed instrument + MVP alleles.

suppressMessages({ library(data.table); library(TwoSampleMR); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

# ---------- exact endpoint counts ----------
mani <- fread(file.path(data_dir,"finngen_R13_manifest.tsv"))
key <- mani[grepl("SLEEPAPNO|HEARTFAIL|OVERWEIGHT|OBESITY", phenocode, ignore.case=TRUE),
            .(phenocode, phenotype, num_cases, num_controls)]
cat("=== FinnGen R13 endpoint counts ===\n"); print(key)
fwrite(key, file.path(res_dir,"finngen_endpoint_counts.csv"))

# ---------- build the instrument set ONCE ----------
meta <- fread(file.path(res_dir,"meta_osa_snps.tsv"))
mv   <- fread(file.path(data_dir,"GCST90475824_p5e8.tsv"))[
          , .(SNP=rsid, effect_allele, other_allele, eaf)]
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP, pval=meta$p_meta),
             clump_r2=0.001, clump_kb=10000, pop="EUR"), error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
ivs <- merge(meta[, .(SNP, beta=beta_meta, se=se_meta, p=p_meta)], mv, by="SNP")
cat("\ninstrument SNPs:", nrow(ivs), "\n"); print(ivs[, .(SNP, beta, p)])

ex <- format_data(as.data.frame(ivs), type="exposure", snp_col="SNP",
      beta_col="beta", se_col="se", eaf_col="eaf",
      effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="p")

# ---------- PheWAS panel ----------
positive <- c("ieu-b-40","ebi-a-GCST006867","ieu-a-2","finn-b-E4_OBESITY",
              "ieu-b-109","ebi-a-GCST005179","ieu-b-107","finn-b-E4_DM2",
              "ieu-b-4952","ieu-a-300","ieu-a-302","finn-b-I9_HYPTENS")
negative <- c("ieu-b-4965","ieu-b-4970","ieu-b-4961","finn-b-C3_COLORECTAL_EXALLC",
              "finn-b-M13_ARTHROSIS","finn-b-H7_CATARACT","finn-b-C3_PROSTATE_EXALLC")
targets <- data.table(id=c(positive,negative),
                      expect=c(rep("positive",length(positive)),
                               rep("negative",length(negative))))
targets <- unique(targets, by="id")
cat("\nPheWAS targets:", nrow(targets), "\n")

res <- list()
for (i in seq_len(nrow(targets))) {
  id <- targets$id[i]
  o <- tryCatch(extract_outcome_data(snps = ivs$SNP, outcomes = id), error=function(e) NULL)
  if (is.null(o) || nrow(o) < 3) { cat(sprintf("%-34s skipped\n", id)); next }
  d <- tryCatch(harmonise_data(ex, o, action=2), error=function(e) NULL)
  if (is.null(d)) { cat(sprintf("%-34s harmonise failed\n", id)); next }
  d <- d[d$mr_keep, ]
  if (nrow(d) < 3) { cat(sprintf("%-34s <3 SNPs\n", id)); next }
  r <- tryCatch(mr(d, method_list=c("mr_ivw","mr_weighted_median")), error=function(e) NULL)
  if (is.null(r)) next
  ri <- r[r$method=="Inverse variance weighted", ]
  res[[length(res)+1]] <- data.table(id=id, expect=targets$expect[i], nsnp=ri$nsnp,
                                     b=ri$b, se=ri$se, pval=ri$pval)
  cat(sprintf("%-34s nsnp=%3d beta=%+.4f p=%.3g\n", id, ri$nsnp, ri$b, ri$pval))
}

if (length(res) > 0) {
  R <- rbindlist(res); R[, fdr := p.adjust(pval, method="BH")]; R[, sig := fdr<0.05]
  cat("\n\n===== PheWAS SUMMARY (sorted by p) =====\n")
  print(R[order(pval), .(id, expect, nsnp, beta=round(b,4), pval=signif(pval,3),
                          fdr=signif(fdr,3), sig)])
  fwrite(R, file.path(res_dir,"PHEWAS_results.csv"))
  cat("\n--- concordance ---\n")
  cat("positive significant:", R[expect=="positive" & sig, .N], "/", R[expect=="positive", .N], "\n")
  cat("negative significant:", R[expect=="negative" & sig, .N], "/", R[expect=="negative", .N], "\n")
} else cat("\nno results\n")
cat("\n=== done ===\n")
