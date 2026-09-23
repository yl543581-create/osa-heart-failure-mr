# 12_sensitivity_final.R
# (a) FTO sensitivity: does the univariable OSA->HF effect survive removal of
#     the chr16 FTO region (+/-500 kb)?
# (b) leave-one-out summary already obtained separately
# (c) FinnGen OSA definition variants for comparison

suppressMessages({ library(TwoSampleMR); library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR
hf_id  <- "ebi-a-GCST009541"

ivs <- fread(file.path(res_dir, "mvp_osa_instruments.tsv"))
cat("IVs:", nrow(ivs), "\n")

run_mr <- function(dt, label) {
  if (nrow(dt) < 3) { cat(label, ": too few SNPs\n"); return(invisible(NULL)) }
  out <- extract_outcome_data(snps = dt$SNP, outcomes = hf_id)
  ex <- format_data(as.data.frame(dt), type = "exposure",
                    snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf",
                    effect_allele_col="effect_allele", other_allele_col="other_allele",
                    pval_col="p")
  d <- harmonise_data(ex, out, action = 2); d <- d[d$mr_keep, ]
  r <- generate_odds_ratios(mr(d, method_list=c("mr_ivw","mr_weighted_median")))
  cat("\n=== ", label, " (n=", nrow(d), " SNPs) ===\n", sep="")
  for (i in seq_len(nrow(r)))
    cat(sprintf("  %-20s OR=%.3f (%.3f-%.3f) p=%.3g\n",
                r$method[i], r$or[i], r$or_lci95[i], r$or_uci95[i], r$pval[i]))
  invisible(r)
}

# (a) full set
run_mr(ivs, "ALL instruments")

# (b) exclude FTO region (chr16: 53.2-54.3 Mb covers rs1421085 @ 53,767,042)
fto_lo <- 53.2e6; fto_hi <- 54.3e6
no_fto <- ivs[!(chromosome == 16 & position >= fto_lo & position <= fto_hi)]
cat("\nFTO-region SNPs removed:", nrow(ivs) - nrow(no_fto), "\n")
if (nrow(ivs) - nrow(no_fto) > 0)
  print(ivs[chromosome == 16 & position >= fto_lo & position <= fto_hi, .(SNP, position, beta, se, p)])
run_mr(no_fto, "EXCLUDING FTO region (chr16 53.2-54.3 Mb)")

# (c) only genome-wide "strong" instruments
strong <- ivs[p < 5e-10]
run_mr(strong, "p < 5e-10 only")

# (d) BMI-proxy check: how many of the 110 IVs sit near known obesity loci?
obesity_loci <- data.table(
  gene = c("FTO","MC4R","TMEM18","BDNF","SEC16B","GNPDA2","SH2B1","FAIM2",
           "NPC1","KCTD15","MTCH2","TFAP2B","ETV5","MAP2K5","NRXN3"),
  chr  = c(16,18,2,11,1,4,16,12,18,19,11,6,3,15,14),
  mb   = c(53.8,58.0,0.6,27.6,177.9,45.2,28.8,50.2,21.1,34.0,47.6,50.8,185.8,67.8,79.0)
)
ivs[, near_obesity := FALSE]
for (i in seq_len(nrow(obesity_loci))) {
  ivs[chromosome == obesity_loci$chr[i] &
      abs(position/1e6 - obesity_loci$mb[i]) < 1,
      near_obesity := TRUE]
}
cat("\nIVs within 1 Mb of a known obesity locus:",
    sum(ivs$near_obesity), "of", nrow(ivs), "\n")
print(ivs[near_obesity == TRUE, .(SNP, chromosome, position, beta, p)])

cat("\n=== done ===\n")
