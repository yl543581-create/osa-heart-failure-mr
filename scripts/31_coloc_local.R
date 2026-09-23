# 31_coloc_local.R
# Colocalisation using ONLY local data (no slow API calls):
#   OSA  : FinnGen R13 G6_SLEEPAPNO_INCLAVO (local, 772 MB)
#   HF   : FinnGen R13 I9_HEARTFAIL          (local, 771 MB)
# Both files carry ref/alt/af_alt, so MAF is available.
#
# NOTE: this coloc is FinnGen-vs-FinnGen, so it tests whether OSA and HF share a
# causal variant WITHIN the same cohort. Sample overlap means this is a
# "same-cohort" coloc -- interpretable, but stated as such.

suppressMessages({ library(data.table); library(coloc) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

ivs <- fread(file.path(res_dir,"META_mvmr_input.tsv"))
comp <- fread(file.path(data_dir,"GCST90475824_p5e8.tsv"))[
          , .(SNP=rsid, chr=chromosome, pos=bp)]
ivs <- merge(ivs, comp, by="SNP", all.x=TRUE); ivs <- ivs[!is.na(pos)]
cat("loci:", nrow(ivs), "\n")

cols <- c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt")
cat("\nreading OSA (selected cols)...\n")
osa <- fread(file.path(data_dir,"finnngen_R13_G6_SLEEPAPNO_INCLAVO.gz"),
             select=cols, sep="\t", header=TRUE, nThread=4, showProgress=FALSE)
setnames(osa, c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
              c("chr","pos","oa","ea","SNP","p_osa","b_osa","se_osa","eaf_osa"))
osa[, SNP := tstrsplit(SNP,",",fixed=TRUE)[[1]]]
osa <- osa[!is.na(SNP) & !is.na(b_osa) & !is.na(se_osa) & se_osa>0]
cat("OSA rows:", nrow(osa), "\n")

cat("reading HF (selected cols)...\n")
hf <- fread(file.path(data_dir,"finngen_R13_I9_HEARTFAIL.gz"),
            select=cols, sep="\t", header=TRUE, nThread=4, showProgress=FALSE)
setnames(hf, c("#chrom","pos","ref","alt","rsids","pval","beta","sebeta","af_alt"),
             c("chr","pos","oa2","ea2","SNP","p_hf","b_hf","se_hf","eaf_hf"))
hf[, SNP := tstrsplit(SNP,",",fixed=TRUE)[[1]]]
hf <- hf[!is.na(SNP) & !is.na(b_hf) & !is.na(se_hf) & se_hf>0]
cat("HF rows:", nrow(hf), "\n")

setkey(osa, SNP); setkey(hf, SNP)

run <- function(chr, pos, snp) {
  lo <- pos-5e5; hi <- pos+5e5
  a <- osa[chr==get("chr") & pos>=lo & pos<=hi]
  b <- hf [chr==get("chr") & pos>=lo & pos<=hi]
  m <- merge(a, b, by="SNP", suffixes=c("",".h"))
  m <- m[!is.na(eaf_osa) & !is.na(eaf_hf) & eaf_osa>0 & eaf_osa<1 & eaf_hf>0 & eaf_hf<1]
  m <- m[!is.na(b_osa) & !is.na(b_hf)]
  # drop duplicate rsIDs (multiallelic / multi-position) -- keep lowest OSA p
  m <- m[order(p_osa)]
  m <- m[!duplicated(SNP)]
  cat(sprintf("\n%s chr%d %d SNPs", snp, chr, nrow(m)))
  if (nrow(m) < 30) { cat(" - too few\n"); return(NULL) }
  flip <- ifelse(m$ea==m$ea2 & m$oa==m$oa2, 1,
          ifelse(m$ea==m$oa2 & m$oa==m$ea2, -1, NA_integer_))
  m[, flip := flip]; m <- m[!is.na(flip)]
  m[, b_hf2 := b_hf*flip]
  if (nrow(m) < 30) { cat(" - too few after align\n"); return(NULL) }
  d1 <- list(snp=m$SNP, beta=m$b_osa, varbeta=m$se_osa^2, type="cc",
             N=74697+430000, MAF=pmin(m$eaf_osa,1-m$eaf_osa), s=74697/(74697+430000))
  d2 <- list(snp=m$SNP, beta=m$b_hf2, varbeta=m$se_hf^2, type="cc",
             N=41591+458595, MAF=pmin(m$eaf_hf,1-m$eaf_hf), s=41591/(41591+458595))
  r <- tryCatch(coloc.abf(d1,d2), error=function(e){cat(" err:",conditionMessage(e),"\n");NULL})
  if (is.null(r)) return(NULL)
  pp <- r$summary
  cat(sprintf(" | PP.H3=%.3f PP.H4=%.3f\n", pp[["PP.H3.abf"]], pp[["PP.H4.abf"]]))
  data.table(SNP=snp, chr=chr, pos=pos, nsnp=nrow(m),
             PP.H0=pp[["PP.H0.abf"]],PP.H1=pp[["PP.H1.abf"]],PP.H2=pp[["PP.H2.abf"]],
             PP.H3=pp[["PP.H3.abf"]],PP.H4=pp[["PP.H4.abf"]])
}
out <- list()
for (i in seq_len(nrow(ivs))) {
  r <- tryCatch(run(ivs$chr[i], ivs$pos[i], ivs$SNP[i]), error=function(e){cat(" FAIL:",conditionMessage(e),"\n");NULL})
  if (!is.null(r)) out[[length(out)+1]] <- r
}
if (length(out)>0) {
  R <- rbindlist(out)
  cat("\n\n===== COLOC SUMMARY (FinnGen OSA vs FinnGen HF) =====\n")
  print(R[, .(SNP, nsnp, PP.H3=round(PP.H3,3), PP.H4=round(PP.H4,3))])
  fwrite(R, file.path(res_dir,"COLOC_results.csv"))
  cat("\nPP.H4>=0.8:", sum(R$PP.H4>=0.8), "of", nrow(R),
      " | PP.H3>=0.8:", sum(R$PP.H3>=0.8), "\n")
} else cat("\nnone\n")
cat("=== done ===\n")
