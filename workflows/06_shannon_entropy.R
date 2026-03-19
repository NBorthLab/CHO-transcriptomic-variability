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
  "SummarizedExperiment"
)

tar_option_set(
  packages = global_packages,
  controller = controller
)

tar_source()

params <- yaml::read_yaml("targets_config.yaml")

plotdir <- "plots/shannon_entropy"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}
resultsdir <- "results/shannon_entropy"
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

# inputs <- list(
#   tar_target(
#     raw_counts_rds,
#     "results/R/preprocessing/objects/raw_counts",
#     format = "file"
#   ),
#   tar_target(
#     raw_counts,
#     command = readRDS(raw_counts_rds)
#   )
# )
#

# ================================================
#   Pipeline
# ================================================

list(
  per_study_entropy <- tar_map(
    values = studies,
    names = "dataset",
    tar_target(
      counts_filt_rds,
      command = str_c(
        "results/R/transformation/objects/counts_filt_", dataset
      )
    ),
    tar_target(
      counts_filt,
      command = readRDS(counts_filt_rds)
    ),
    tar_target(
      rpk,
      command = calculate_rpk2(counts_filt)
    ),
    tar_target(
      dge_norm,
      command = normalize_tmm(rpk)
    ),
    tar_target(
      cpm,
      command = get_cpm(dge_norm)
    ),
    tar_target(
      logcpm,
      command = get_cpm(dge_norm, log = TRUE)
    ),
    tar_target(
      gene_specificity,
      command = calculate_entropy(cpm)
    )
  ),
  tar_combine(
    norm_counts,
    per_study_entropy[["dge_norm"]],
    command = list(!!!.x)
  ),
  tar_combine(
    cpms,
    per_study_entropy[["cpm"]],
    command = list(!!!.x)
  ),
  tar_combine(
    logcpms,
    per_study_entropy[["logcpm"]],
    command = list(!!!.x)
  ),
  tar_combine(
    entropies,
    per_study_entropy[["gene_specificity"]],
    command = bind_rows(!!!.x)
  ),

  #
  # PCA
  #
  tar_target(
    entropies_matrix,
    command = filter_metric(entropies)
  ),
  tar_target(
    rank_matrix,
    command = apply(entropies_matrix, 2L, rank)
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
