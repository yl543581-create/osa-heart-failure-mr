# 01_search_opengwas.R
# Find usable OSA / heart failure / cardiac-MRI datasets on OpenGWAS,
# then report the REAL instrument counts (the "閻㈢喐顒村Λ鈧? numbers).

suppressMessages({ library(TwoSampleMR); library(jsonlite); library(data.table) })
source(file.path(dirname(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))), "_bootstrap.R"))

`%||%` <- function(a, b) if (is.null(a)) b else a

cat("=== OpenGWAS connectivity (authenticated) ===\n")
jwt <- Sys.getenv("OPENGWAS_JWT")
if (!nzchar(jwt)) stop("OPENGWAS_JWT not set")

api <- function(path) {
  r <- httr::GET(paste0("https://api.opengwas.io/api/", path),
                 httr::add_headers(Authorization = paste("Bearer", jwt)),
                 httr::timeout(60))
  if (httr::status_code(r) != 200) {
    stop(sprintf("HTTP %s for %s :: %s", httr::status_code(r), path,
                 substr(httr::content(r, "text", encoding = "UTF-8"), 1, 300)))
  }
  httr::content(r, "text", encoding = "UTF-8")
}

cat("status:", substr(api("status"), 1, 200), "\n\n")

query_gwas <- function(q, pagesize = 100) {
  body <- list(query = q, pagesize = pagesize)
  r <- httr::POST("https://api.opengwas.io/api/gwasinfo",
                  httr::add_headers(Authorization = paste("Bearer", jwt)),
                  body = jsonlite::toJSON(body, auto_unbox = TRUE),
                  httr::content_type_json(), httr::timeout(120))
  if (httr::status_code(r) != 200) {
    stop(sprintf("HTTP %s :: %s", httr::status_code(r),
                 substr(httr::content(r, "text", encoding = "UTF-8"), 1, 300)))
  }
  x <- fromJSON(httr::content(r, "text", encoding = "UTF-8"), flatten = TRUE)
  if (!is.data.frame(x) || nrow(x) == 0) return(NULL)
  x
}

show <- function(x, label) {
  cat("\n--- ", label, " (n=", ifelse(is.null(x), 0, nrow(x)), ") ---\n", sep = "")
  if (is.null(x)) { cat("(no results)\n"); return(invisible(NULL)) }
  keep <- intersect(c("id","trait","sample_size","ncase","ncontrol",
                      "population","year","author","sex","nsnp","build"),
                    names(x))
  print(as.data.table(x)[, ..keep][order(-sample_size)][1:min(15, .N)])
  invisible(x)
}

osa <- query_gwas("sleep apnoea OR sleep apnea OR apnoea OR apnea")
show(osa, "OSA-ish datasets")

hf <- query_gwas("heart failure")
show(hf, "heart failure datasets")

cmr <- query_gwas("cardiac magnetic resonance OR ventricular")
show(cmr, "cardiac MRI / ventricular datasets")

cat("\n\n=== INSTRUMENT COUNTS (P < 5e-8) ===\n")
test_ids <- c("finn-b-G6_SLEEPAPNO", "finn-b-G6_SLEEPAPNO_INCLAVO",
              "finn-b-I9_HEARTFAIL", "ieu-b-108", "ebi-a-GCST009541")
for (id in test_ids) {
  out <- tryCatch({
    d <- extract_instruments(outcomes = id, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
    sprintf("OK  n_SNP=%d", nrow(d))
  }, error = function(e) paste("FAIL:", conditionMessage(e)))
  cat(sprintf("%-32s %s\n", id, out))
}

cat("\n=== done ===\n")
