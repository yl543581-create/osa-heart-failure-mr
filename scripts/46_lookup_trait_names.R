# 46_lookup_trait_names_v2.R
# Use the POST /api/gwasinfo endpoint (known to work in this environment).
suppressMessages({ library(httr); library(jsonlite); library(data.table) })
jwt <- Sys.getenv("OPENGWAS_JWT")
ids <- c("ieu-b-40","ebi-a-GCST006867","ieu-a-2","finn-b-E4_OBESITY",
         "ieu-b-109","ebi-a-GCST005179","ieu-b-107","finn-b-E4_DM2",
         "ieu-a-300","ieu-a-302","finn-b-I9_HYPTENS",
         "ieu-b-4965","ieu-b-4970","ieu-b-4961",
         "finn-b-C3_COLORECTAL_EXALLC","finn-b-M13_ARTHROSIS",
         "finn-b-C3_PROSTATE_EXALLC")

r <- httr::POST("https://api.opengwas.io/api/gwasinfo",
                httr::add_headers(Authorization = paste("Bearer", jwt)),
                body = jsonlite::toJSON(list(id = ids), auto_unbox = TRUE),
                httr::content_type_json(), httr::timeout(180))
cat("HTTP:", httr::status_code(r), "\n")
if (httr::status_code(r) == 200) {
  x <- fromJSON(httr::content(r, "text", encoding="UTF-8"), flatten = TRUE)
  d <- as.data.table(x)
  cat("returned rows:", nrow(d), "\n")
  cat("columns:", paste(names(d), collapse=", "), "\n\n")
  keep <- intersect(c("id","trait","ncase","ncontrol","sample_size","population","sex","year"), names(d))
  print(d[, ..keep])
  fwrite(d[, ..keep], "D:/harness/osa_hf/results/PHEWAS_trait_names.csv", sep="\t")
  cat("\nsaved: results/PHEWAS_trait_names.csv\n")
}
