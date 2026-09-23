# 45_lookup_trait_names.R
# Resolve the true trait names for the PheWAS panel, so the OSF/paper text
# describes the traits accurately (some IDs were mislabelled in the draft).
suppressMessages({ library(httr); library(jsonlite); library(data.table) })
jwt <- Sys.getenv("OPENGWAS_JWT")
ids <- c("ieu-b-40","ebi-a-GCST006867","ieu-a-2","finn-b-E4_OBESITY",
         "ieu-b-109","ebi-a-GCST005179","ieu-b-107","finn-b-E4_DM2",
         "ieu-a-300","ieu-a-302","finn-b-I9_HYPTENS",
         "ieu-b-4965","ieu-b-4970","ieu-b-4961",
         "finn-b-C3_COLORECTAL_EXALLC","finn-b-M13_ARTHROSIS",
         "finn-b-C3_PROSTATE_EXALLC")
out <- list()
for (id in ids) {
  r <- tryCatch(httr::GET(paste0("https://api.opengwas.io/api/gwasinfo/", id),
                httr::add_headers(Authorization = paste("Bearer", jwt)),
                httr::timeout(60)), error = function(e) NULL)
  if (is.null(r) || httr::status_code(r) != 200) {
    cat(sprintf("%-32s LOOKUP_FAILED\n", id)); next
  }
  j <- fromJSON(httr::content(r, "text", encoding="UTF-8"))
  out[[length(out)+1]] <- data.table(
    id = id,
    trait = if (!is.null(j$trait)) j$trait else NA_character_,
    ncase = if (!is.null(j$ncase)) j$ncase else NA,
    ncontrol = if (!is.null(j$ncontrol)) j$ncontrol else NA,
    sample_size = if (!is.null(j$sample_size)) j$sample_size else NA,
    population = if (!is.null(j$population)) j$population else NA)
  cat(sprintf("%-32s %s\n", id, out[[length(out)]]$trait))
}
if (length(out)) {
  R <- rbindlist(out)
  cat("\n=== FULL TABLE ===\n"); print(R)
  fwrite(R, "D:/harness/osa_hf/results/PHEWAS_trait_names.csv", sep="\t")
  cat("\nsaved: results/PHEWAS_trait_names.csv\n")
}
