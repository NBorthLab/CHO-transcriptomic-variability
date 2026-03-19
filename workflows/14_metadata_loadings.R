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

plotdir <- "plots/metadata_loadings"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}


list(
  tar_quarto(
    report,
    "analysis/14_metadata_loadings.qmd"
  )
)
