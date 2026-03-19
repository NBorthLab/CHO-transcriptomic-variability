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
  "tidyplots",
  "ggpubr"
)

tar_option_set(
  packages = global_packages,
  controller = controller
)

tar_source()

params <- yaml::read_yaml("targets_config.yaml")

plotdir <- "plots/principal_components"
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


# ================================================
#   Pipeline
# ================================================

list(
  per_study_metric <- tar_map(
    values = studies,
    names = "dataset",
    tar_target(
      std_rds,
      command = str_c(
        "results/R/transformation/objects/std_", blindness, "_", dataset
      ),
      format = "file"
    ),
    tar_target(
      std,
      command = readRDS(std_rds)
    )
  ),
  tar_combine(
    stds,
    per_study_metric[["std"]],
    command = bind_rows(!!!.x)
  ),
  tar_target(
    stds_matrix,
    command = filter_metric(stds)
  ),
  tar_target(
    rank_matrix,
    command = apply(stds_matrix, 2L, rank)
  ),
  tar_target(
    cor,
    command = correlate_metric(stds_matrix)
  ),
  tar_target(
    cor_plot,
    command = plot_correlation_matrix(cor),
    format = "file"
  ),
  tar_target(
    cor_hist,
    command = plot_correlation_histogram(cor),
    format = "file"
  ),

  #
  # PCA
  #
  # WARN: This part should be jackknifed
  tar_target(
    eigen,
    command = calculate_eigenvectors(cor)
  ),
  tar_target(
    principal_components,
    command = get_principal_components(rank_matrix, eigen$vectors)
  ),
  tar_target(
    loadings,
    command = calculate_loadings(eigen)
  ),
  tar_target(
    first_pc,
    command = get_first_component(principal_components)
  ),
  tar_target(
    scree_plot,
    command = plot_scree(eigen),
    format = "file"
  ),
  tar_target(
    pca_plot,
    command = plot_pc(principal_components, eigen),
    format = "file"
  ),
  tar_target(
    global_variability_genes,
    command = {
      file <- "results/principal_components/global_variability_genes.csv"
      write_csv(first_pc, file)
      file
    },
    format = "file"
  ),

  #
  # Subsets of first component
  #
  tar_target(
    variable_genes,
    command = filter(first_pc, decile == 1)
  ),
  tar_target(
    stable_genes,
    command = filter(first_pc, decile == 10)
  ),
  # Used for chromatin states analysis
  # tar_target(dec, 1:10),
  # tar_target(
  #   dec_genes,
  #   command = extract_genes(first_pc, dec),
  #   pattern = map(dec)
  # ),
  tar_map(
    values = list(dec = 1:10),
    tar_target(
      dec_genes,
      command = extract_genes(first_pc, dec),
      format = "file"
    )
  ),

  #
  tar_quarto(
    report,
    "analysis/04_principal_components.qmd"
  ),
  NULL
) %>%
  tar_hook_before(
    hook = conflicted::conflicts_prefer(
      dplyr::filter,
      dplyr::rename
    )
  )
