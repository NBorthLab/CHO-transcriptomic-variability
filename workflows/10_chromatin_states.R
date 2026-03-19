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

plotdir <- "plots/chromatin_states"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
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

#   NOTE:
#   This analysis requires additional dependencies.
#   Install the conda environments under `envs/chromhmm.yml` and
#   `envs/python.yaml`.


combinations_dec_tp <- expand.grid(tp = paste0("Tp", 0:17), dec = 1:10)

list(

  #
  # Overlap enrichment
  #
  per_decile_timepoint <- tar_map(
    values = combinations_dec_tp,
    names = "dec",
    tar_target(
      overlap,
      command = load_overlap_results(dec, tp)
    )
  ),
  tar_combine(
    overlaps,
    per_decile_timepoint[["overlap"]],
    command = bind_rows(!!!.x)
  ),

  #
  # Coverage
  #
  per_decile <- tar_map(
    values = list(dec = 1:10),
    tar_target(
      sub_annotation,
      command = {
        decile_str <- paste0("decile", str_pad(dec, 2, pad = "0"))
        str_glue("results/chromatin_states/{decile_str}_genes.gtf")
      }
    ),
    per_timepoint <- tar_map(
      values = list(tp = paste0("Tp", 0:17)),
      tar_target(
        coverage,
        command = calculate_coverage(dec, tp)
      )
    ),
    tar_combine(
      decile_coverages,
      per_timepoint[["coverage"]],
      command = bind_rows(!!!.x)
    )
  ),
  tar_combine(
    coverages,
    per_decile[["decile_coverages"]],
    command = load_coverages(!!!.x)
  ),
  tar_target(
    emissions,
    command = load_emissions()
  ),
  tar_quarto(
    report,
    "analysis/10_chromatin_states.qmd"
  )
) %>%
  tar_hook_before(
    hook = conflicted::conflicts_prefer(
      dplyr::filter,
      dplyr::rename
    )
  )
