build_workflowr <- function(full = FALSE) {
  if (full) {
    rmd_files <- list.files("analysis", pattern = "*.Rmd", full.names = TRUE)
    workflowr::wflow_build(rmd_files)
  } else {
    workflowr::wflow_build()
  }
}

get_assay <- function(se, assay) {
  matrix <- assay(se, assay)
  SummarizedExperiment(
    assay = matrix,
    colData = colData(se)
  )
}

mtx2tibble <- function(mtx) {
  mtx %>%
    as_tibble(rownames = "gene") %>%
    pivot_longer(-gene, names_to = "samples", values_to = "values")
}
