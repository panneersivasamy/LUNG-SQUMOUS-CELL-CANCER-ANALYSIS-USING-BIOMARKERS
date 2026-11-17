library(TCGAbiolinks)
library(SummarizedExperiment)
query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
library(TCGAbiolinks)
options(timeout = 600)
httr::set_config(httr::config(connecttimeout = 120, timeout = 600, http_version = 0))

query <- NULL
for (i in 1:30) {
  query <- try(GDCquery(
    project = "TCGA-LUSC",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  ), silent = TRUE)
  if (!inherits(query, "try-error")) { message("Query OK!"); break }
  message("Query retry ", i, " ...")
  Sys.sleep(10)
}
saveRDS(query, "query_LUSC.rds")
for (i in 1:30) {
  res <- try(GDCdownload(query, method = "api", files.per.chunk = 3), silent = TRUE)
  if (!inherits(res, "try-error")) { message("Done!"); break }
  message("Download retry ", i, " ...")
  Sys.sleep(15)
}
length(list.files("GDCdata/TCGA-LUSC", recursive = TRUE, pattern = "\\.tsv$"))
options(download.file.method = "wininet")
options(timeout = 600)
con <- url("https://api.gdc.cancer.gov/status")
readLines(con, warn = FALSE)
close(con)
options(download.file.method = "wininet")
readLines("https://api.gdc.cancer.gov/status", warn = FALSE)
readLines("https://api.gdc.cancer.gov/status", warn = FALSE)
length(list.files("GDCdata/TCGA-LUSC", recursive = TRUE, pattern = "\\.tsv$"))
library(SummarizedExperiment)
data <- GDCprepare(query)
saveRDS(data, "TCGA_LUSC_Data.rds")
data
length(list.files("GDCdata/TCGA-LUSC", recursive = TRUE, pattern = "\\.tsv$"))
library(TCGAbiolinks)
library(SummarizedExperiment)
query <- readRDS("query_LUSC.rds")

