# 39_finalize_outputs.R
# Trim the Steiger output to the interpretable subset and consolidate
# the endpoint counts for the methods section.

suppressMessages(library(data.table))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
res_dir <- RES_DIR

# ---- Steiger: keep the instrument SNPs, plus a summary of the regional set ----
s <- fread(file.path(res_dir, "STEIGER_local.csv"))
cat("full Steiger rows:", nrow(s), "\n")

ivs <- fread(file.path(res_dir, "META_mvmr_input.tsv"))
inst <- s[SNP %in% ivs$SNP]
cat("instrument SNP rows:", nrow(inst), "\n")

# summary rows: overall regional + instrument-only
summ <- data.table(
  set = c("all SNPs in instrument loci", "independent instrument SNPs only"),
  n_snp = c(nrow(s), nrow(inst)),
  sum_R2_OSA = c(sum(s$r2_osa), sum(inst$r2_osa)),
  sum_R2_HF  = c(sum(s$r2_hf),  sum(inst$r2_hf)),
  ratio = c(sum(s$r2_osa)/sum(s$r2_hf), sum(inst$r2_osa)/sum(inst$r2_hf)),
  n_correct = c(sum(s$correct), sum(inst$correct))
)
summ[, pct_correct := round(100*n_correct/n_snp, 1)]
summ[, R2_ratio := round(ratio, 2)]

cat("\n=== Steiger summary ===\n"); print(summ[, .(set, n_snp, R2_ratio, pct_correct)])
cat("\n=== instrument-level detail ===\n")
print(inst[, .(SNP, r2_OSA=signif(r2_osa,6), r2_HF=signif(r2_hf,6),
               ratio=round(r2_osa/r2_hf,1), correct)])

fwrite(summ, file.path(res_dir, "STEIGER_summary.csv"))
fwrite(inst, file.path(res_dir, "STEIGER_instruments.csv"))
file.remove(file.path(res_dir, "STEIGER_local.csv"))   # 12.5 MB, superseded
cat("\nremoved oversized STEIGER_local.csv; wrote summary + instrument files\n")
cat("=== done ===\n")
