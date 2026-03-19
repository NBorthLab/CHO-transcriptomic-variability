suppressPackageStartupMessages({
  library(tibble)
  library(targets)
  library(tarchetypes)
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
  "tidyplots"
)

tar_option_set(
  packages = global_packages,
  controller = controller
)

tar_source()

params <- yaml::read_yaml("targets_config.yaml")


plotdir <- "plots/loess_residuals"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}
resultsdir <- "results/loess_residuals"
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

tar_load(coding_counts, store = "results/R/preprocessing")


# ================================================
#   Pipeline
# ================================================


list(
  per_study_metric <- tar_map(
    values = studies,
    names = "dataset",
    tar_target(
      dds_rds,
      command = str_c("results/R/transformation/objects/dds_", dataset)
    ),
    tar_target(
      dds,
      command = readRDS(dds_rds)
    ),
    tar_target(
      expression_stats,
      command = get_expression_statistics(dds)
    ),
    tar_target(
      loess_fit,
      command = fit_loess(expression_stats)
    ),
    tar_target(
      loess_residuals,
      command = get_loess_residuals(loess_fit)
    ),
    tar_target(
      plot_residuals,
      command = plot_loess_residuals(loess_residuals),
      format = "file"
    )
  ),
  tar_combine(
    loess_resids,
    per_study_metric[["loess_residuals"]],
    command = bind_rows(!!!.x)
  ),
  tar_target(
    loess_matrix,
    command = filter_metric(loess_resids)
  ),
  tar_target(
    rank_matrix,
    command = apply(loess_matrix, 2L, rank)
  ),
  tar_target(
    cor,
    command = correlate_metric(rank_matrix)
  ),
  tar_target(
    eigen,
    command = calculate_eigenvectors(cor)
  ),
  tar_target(
    principal_components,
    command = get_principal_components(rank_matrix, eigen$vectors)
  ),
  tar_target(
    first_pc,
    command = get_first_component(principal_components)
  )
) %>%
  tar_hook_before(
    hook = conflicted::conflicts_prefer(
      dplyr::filter,
      dplyr::rename
    )
  )
