# 40_check_overweight_endpoint.R
# Before downloading 771 MB, check what I9_HEARTFAIL_AND_OVERWEIGHT actually is
# and how it compares with the plain heart-failure endpoint.

suppressMessages(library(data.table))
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))
mani <- fread(file.path(ROOT, "data/finngen_R13_manifest.tsv"))

cols <- intersect(c("phenocode","phenotype","category","num_cases",
                    "num_controls","path_https"), names(mani))
key <- mani[grepl("HEARTFAIL|OVERWEIGHT|OBESITY|SLEEPAPNO", phenocode, ignore.case=TRUE), ..cols]
setorder(key, -num_cases)
cat("=== all heart-failure / overweight / OSA endpoints ===\n")
print(key)

cat("\n=== specifically the comparison we want ===\n")
want <- c("I9_HEARTFAIL", "I9_HEARTFAIL_AND_OVERWEIGHT",
          "I9_HEARTFAIL_AND_CHD", "I9_HEARTFAIL_AND_HYPERTCARDIOM")
cmp <- mani[phenocode %in% want, ..cols]
print(cmp)

cat("\n=== interpretation check ===\n")
cat("If I9_HEARTFAIL_AND_OVERWEIGHT has FEWER cases than I9_HEARTFAIL,\n")
cat("it is a SUBSET (HF patients who are also overweight) -> a stratified\n")
cat("endpoint. Under the adiposity-mediated model we would predict the OSA\n")
cat("effect to be AT LEAST as strong here as for all HF.\n")

# size estimate: FINNGEN files are ~770 MB each
cat("\n=== download URL ===\n")
u <- mani[phenocode=="I9_HEARTFAIL_AND_OVERWEIGHT", path_https]
cat(u, "\n")
