# 30_diagnose_coloc.R
suppressMessages({ library(data.table); library(ieugwasr) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
# pick the chr1 locus
snps <- c("rs10875017","rs11672660","rs62048402")
for (id in c("finn-b-G6_SLEEPAPNO","ebi-a-GCST009541")) {
  cat("\n=====", id, "=====\n")
  r <- tryCatch(ieugwasr::associations(variants=snps, id=id), error=function(e){cat("FAIL:",conditionMessage(e),"\n");NULL})
  if (!is.null(r)) { d <- as.data.table(r); cat("cols:", paste(names(d),collapse=", "), "\n"); print(d) }
}
cat("\n=== regional query test: does a string region work? ===\n")
for (v in c("1:95980000-96980000", "chr1:95980000-96980000")) {
  r <- tryCatch(ieugwasr::associations(variants=v, id="ebi-a-GCST009541"), error=function(e){cat("err:",conditionMessage(e),"\n");NULL})
  cat(v, "->", ifelse(is.null(r), "NULL", nrow(r)), "rows\n")
}
cat("\n=== how many of the 99 local SNPs does HF return? ===\n")
comp <- fread(file.path(ROOT, "data/GCST90475824_p5e8.tsv"))
loc <- comp[chromosome==1 & bp>=95980000 & bp<=96980000]
cat("local SNPs in region:", nrow(loc), "\n")
r <- tryCatch(ieugwasr::associations(variants=loc$rsid, id="ebi-a-GCST009541"), error=function(e) NULL)
if (!is.null(r)) {
  d <- as.data.table(r)
  cat("HF returned:", nrow(d), "\n")
  cat("with non-NA eaf:", sum(!is.na(d$eaf)), "\n")
  cat("with non-NA beta:", sum(!is.na(d$beta)), "\n")
  print(head(d, 5))
}
