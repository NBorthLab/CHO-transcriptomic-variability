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

plotdir <- "plots/motifs"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}
resultsdir <- "results/motifs"
if (!dir.exists(resultsdir)) {
  dir.create(resultsdir, recursive = TRUE)
}


# ================================================
#   Inputs
# ================================================

studies <- readr::read_csv("resources/datasets.csv",
  progress = FALSE,
  show_col_types = FALSE
)

tar_load(first_pc, store = "results/R/principal_components")


# ================================================
#   Pipeline
# ================================================

list(
  tar_map(
    values = list(motif = c("TATA", "Inr")),
    per_decile_motifs <- tar_map(
      values = list(decile = 1:10),
      tar_target(
        best_site_np,
        paste0(
          resultsdir, "/fimo/", decile, ".", motif, ".best_site.narrowPeak"
        ),
        format = "file"
      ),
      tar_target(
        best_site,
        plyranges::read_narrowpeaks(best_site_np)
      )
    ),
    tar_combine(
      decile_motifs,
      per_decile_motifs[["best_site"]],
      command = GenomicRanges::GRangesList(!!!.x)
    )
  )
)
