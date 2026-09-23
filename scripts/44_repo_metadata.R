# Generate repository metadata: sessionInfo, file manifest, and checksums.
# Run from the repository root or anywhere (paths resolve via config.R).

source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

si_path <- file.path(RES_DIR, "sessionInfo.txt")
con <- file(si_path, "w")
writeLines(c(
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("Repository root:", ROOT),
  ""
), con)
writeLines(capture.output(sessionInfo()), con)
close(con)
cat("wrote:", si_path, "\n")

# --- results file manifest with sizes and md5 ---
files <- list.files(RES_DIR, full.names = TRUE)
files <- files[file.info(files)$isdir == FALSE]
man <- data.frame(
  file = basename(files),
  bytes = file.size(files),
  md5 = vapply(files, function(f) tools::md5sum(f)[[1]], character(1)),
  stringsAsFactors = FALSE
)
man <- man[order(man$file), ]
out <- file.path(RES_DIR, "RESULTS_MANIFEST.tsv")
write.table(man, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("wrote:", out, " (", nrow(man), " files )\n", sep = "")

# --- data file provenance (no checksums for multi-GB files; sizes only) ---
if (dir.exists(DATA_DIR)) {
  dfiles <- list.files(DATA_DIR, full.names = TRUE)
  dfiles <- dfiles[file.info(dfiles)$isdir == FALSE]
  dman <- data.frame(
    file = basename(dfiles),
    MB = round(file.size(dfiles)/1e6, 1),
    stringsAsFactors = FALSE
  )
  dman <- dman[order(dman$file), ]
  dout <- file.path(DATA_DIR, "DATA_MANIFEST.tsv")
  write.table(dman, dout, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("wrote:", dout, " (", nrow(dman), " files )\n", sep = "")
  print(dman)
}

cat("\n=== done ===\n")
