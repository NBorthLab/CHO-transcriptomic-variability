suppressPackageStartupMessages({
  library(tibble)
  library(targets)
  library(tarchetypes)
  library(tidyverse)
  library(crew)
  library(SummarizedExperiment)
})


controller <- crew_controller_local(
  name = "my_controller",
  workers = 8
)

global_packages <- c(
  "tidyverse",
  "conflicted",
  "SummarizedExperiment",
  "tidyplots",
  "ggpubr"
)

tar_option_set(
  packages = global_packages,
  controller = controller
)

tar_source()

params <- yaml::read_yaml("targets_config.yaml")


# Get all R package versions


list(
  tar_target(
    r_version,
    command = paste0(R.Version()$major, ".", R.Version()$minor)
  ),
  tar_target(
    session_info,
    command = sessioninfo::session_info(
      to_file = here::here("results/R_session_info.txt")
    ),
    cue = tar_cue(mode = "always"),
    packages = c(
      "DESeq2",
      "edgeR",
      "clusterProfiler",
      "AnnotationHub",
      "GenomicDistributions",
      "moments"
    )
  ),
  tar_target(
    package_versions,
    command = session_info$packages %>% pull(loadedversion, name = package)
  ),
  tar_target(
    versions_RData,
    command = {
      file <- here::here("results/versions.RData")
      save(r_version, package_versions, file = file)
      file
    },
    format = "file"
  )
)
