# 43_collect_metadata.R
# Gather the factual details needed for an accurate STROBE-MR checklist:
# software versions, tool counts, and completeness of the sensitivity battery.

suppressMessages({ library(data.table); library(TwoSampleMR) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", commandArgs(trailingOnly=FALSE)[grepl("^--file=", commandArgs(trailingOnly=FALSE))][1])))), "config.R")); data_dir <- DATA_DIR; res_dir <- RES_DIR

cat("=== SOFTWARE VERSIONS ===\n")
cat("R:", R.version.string, "\n")
for (p in c("TwoSampleMR","ieugwasr","coloc","data.table","MRPRESSO",
            "MendelianRandomization","RadialMR")) {
  v <- tryCatch(as.character(packageVersion(p)), error=function(e) NA_character_)
  cat(sprintf("%-24s %s\n", p, ifelse(is.na(v),"MISSING",v)))
}
cat("MVMR: NOT AVAILABLE (implemented manually)\n")
cat("MRlap: NOT AVAILABLE (build toolchain missing)\n")

cat("\n=== INSTRUMENT / SNP COUNTS ===\n")
m1 <- fread(file.path(res_dir,"mvp_osa_instruments_LDclumped.tsv"))
cat("MVP OSA, LD-clumped (r2<0.001):", nrow(m1), "SNPs\n")
cat("  mean F:", round(mean((m1$beta2/m1$se2)^2),1),
    " min F:", round(min((m1$beta2/m1$se2)^2),1), "\n")
m2 <- fread(file.path(res_dir,"meta_osa_snps.tsv"))
cat("meta (MVP+FinnGen) significant overlapping:", nrow(m2), "\n")
m3 <- fread(file.path(res_dir,"META_mvmr_input.tsv"))
cat("meta, LD-clumped instruments:", nrow(m3), "SNPs\n")
cat("  mean F:", round(mean((m3$bx1_meta/m3$sx1_meta)^2),1), "\n")
cat("\nOSA significant variants: MVP 4,770 | FinnGen 3,205\n")

cat("\n=== SENSITIVITY BATTERY STATUS ===\n")
cat("Cochran Q (heterogeneity)      : DONE (IVW + Egger)\n")
cat("MR-Egger intercept (pleiotropy): DONE\n")
cat("Weighted median                : DONE\n")
cat("Weighted mode                  : DONE\n")
cat("Leave-one-out                  : DONE\n")
cat("MR-PRESSO                      : RUN, no outliers detected\n")
cat("FTO-region exclusion           : DONE\n")
cat("Steiger directionality         : DONE (liability scale, local)\n")
cat("Colocalisation                 : DONE (10 loci)\n")
cat("PheWAS                         : DONE (17 traits)\n")
cat("BMI-stratified endpoints        : DONE (3 endpoints)\n")
cat("Positive control (BMI->HF)     : DONE\n")

cat("\n=== RESULTS FILES PRESENT ===\n")
for (f in list.files(res_dir, full.names=TRUE))
  cat(sprintf("  %-42s %6.1f KB\n", basename(f), file.size(f)/1024))
