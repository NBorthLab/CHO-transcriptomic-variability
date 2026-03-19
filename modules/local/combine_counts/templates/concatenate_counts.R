#!/usr/bin/env Rscript --vanilla

library(tximport)

files <- list.files("counts",
  pattern = "*.genes.results",
  recursive = TRUE,
  full.names = TRUE
)
filenames <- basename(files)
samples <- gsub("(\\\\w+)\\\\.genes\\\\.results", "\\\\1", filenames)
names(files) <- samples

all_gene_counts <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
saveRDS(all_gene_counts, "all_gene_counts.rds")

r.version <- strsplit(version[["version.string"]], " ")[[1]][3]
tximport.version <- as.character(packageVersion("tximport"))

writeLines(
  c(
    '"${task.process}":',
    paste("    bioconductor-tximport:", tximport.version)
  ),
  "versions.yml"
)
