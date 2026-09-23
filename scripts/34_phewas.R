# 34_phewas.R
# Instrument PheWAS: do the OSA instruments hit expected traits and miss
# unrelated ones? Also extract exact endpoint counts for the methods section.

suppressMessages({ library(data.table); library(TwoSampleMR); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

# ---------- exact counts from the manifest ----------
m <- fread(file.path(data_dir,"finngen_R13_manifest.tsv"))
key <- m[grepl("SLEEPAPNO|HEARTFAIL|OVERWEIGHT|OBESITY", phenocode, ignore.case=TRUE),
         .(phenocode, phenotype, num_cases, num_controls)]
cat("=== FinnGen R13 endpoint counts ===\n"); print(key)
fwrite(key, file.path(res_dir,"finngen_endpoint_counts.csv"))

# ---------- instrument set ----------
meta <- fread(file.path(res_dir,"meta_osa_snps.tsv"))
cl <- tryCatch(ieugwasr::ld_clump(data.frame(rsid=meta$SNP, pval=meta$p_meta),
             clump_r2=0.001, clump_kb=10000, pop="EUR"), error=function(e) NULL)
if (!is.null(cl)) meta <- meta[SNP %in% cl$rsid]
snps <- meta$SNP
cat("\ninstrument SNPs for PheWAS:", length(snps), "\n")

# ---------- PheWAS panel ----------
panel <- data.table(
  id = c(
    # expected POSITIVE (obesity / metabolic / sleep-related)
    "ieu-b-40",            # BMI
    "ebi-a-GCST006867",    # BMI (GIANT)
    "ieu-a-2",             # BMI
    "finn-b-E4_OBESITY",   # obesity
    "ieu-b-109",           # waist-hip ratio
    "ebi-a-GCST005179",    # waist-hip ratio adj BMI
    "ieu-b-107",           # body fat percentage
    "finn-b-E4_DM2",       # type 2 diabetes
    "ieu-b-4952",          # fasting glucose
    "ieu-a-300",           # HDL cholesterol
    "ieu-a-302",           # triglycerides
    "finn-b-I9_HYPTENS",   # hypertension
    "ieu-b-109"            # (dup-safe)
  ),
  expect = c(rep("positive", 13))
)
# de-duplicate ids
panel <- unique(panel, by="id")

# expected NEGATIVE controls: should NOT be causally affected by OSA
neg <- data.table(
  id = c("ieu-b-4965","ieu-b-4970","ieu-b-4961","ieu-b-4972",
         "finn-b-C3_COLORECTAL_EXALLC","finn-b-M13_ARTHROSIS",
         "finn-b-H7_CATARACT","finn-b-C3_PROSTATE_EXALLC"),
  expect = "negative"
)
panel <- rbind(panel, neg)
cat("PheWAS targets:", nrow(panel), "\n")

# ---------- run ----------
res <- list()
for (i in seq_len(nrow(panel))) {
  id <- panel$id[i]
  o <- tryCatch(extract_outcome_data(snps = snps, outcomes = id), error=function(e) NULL)
  if (is.null(o) || nrow(o) < 3) { cat(sprintf("%-32s skipped (no data)\n", id)); next }
  ex <- format_data(as.data.frame(meta[, .(SNP, beta=beta_meta, se=se_meta, p=p_meta)]),
                    type="exposure", snp_col="SNP", beta_col="beta", se_col="se",
                    effect_allele_col="effect_allele", other_allele_col="other_allele",
                    eaf_col="eaf", pval_col="p")
  ex <- extract_instruments(outcomes="finn-b-G6_SLEEPAPNO", p1=5e-8, clump=TRUE)
  ex <- ex[ex$SNP %in% snps, ]
  if (nrow(ex) < 3) { cat(sprintf("%-32s skipped (too few IVs)\n", id)); next }
  d <- tryCatch(harmonise_data(ex, o, action=2), error=function(e) NULL)
  if (is.null(d)) { cat(sprintf("%-32s harmonise failed\n", id)); next }
  d <- d[d$mr_keep, ]
  if (nrow(d) < 3) { cat(sprintf("%-32s <3 SNPs after harmonise\n", id)); next }
  r <- tryCatch(mr(d, method_list="mr_ivw"), error=function(e) NULL)
  if (is.null(r)) next
  res[[length(res)+1]] <- data.table(id=id, expect=panel$expect[i], nsnp=r$nsnp[1],
                                     b=r$b[1], se=r$se[1], pval=r$pval[1])
  cat(sprintf("%-32s nsnp=%3d beta=%+.4f p=%.3g\n", id, r$nsnp[1], r$b[1], r$pval[1]))
}

if (length(res) > 0) {
  R <- rbindlist(res)
  R[, fdr := p.adjust(pval, method="BH")]
  R[, sig := fdr < 0.05]
  cat("\n\n===== PheWAS SUMMARY =====\n")
  print(R[order(pval), .(id, expect, nsnp, beta=round(b,4), pval=signif(pval,3),
                          fdr=signif(fdr,3), sig)])
  fwrite(R, file.path(res_dir,"PHEWAS_results.csv"))
  cat("\n--- concordance with expectation ---\n")
  cat("positive controls significant:", R[expect=="positive" & sig, .N],
      "/", R[expect=="positive", .N], "\n")
  cat("negative controls significant:", R[expect=="negative" & sig, .N],
      "/", R[expect=="negative", .N], "\n")
} else cat("\nno PheWAS results\n")
cat("\n=== done ===\n")
