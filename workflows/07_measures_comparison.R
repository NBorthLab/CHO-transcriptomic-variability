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

params <- yaml::read_yaml("targets_config.yaml")

# ================================================
#   Functions
# ================================================


combine_ranks <- function(std, entropy, loess) {
  combined_ranks <- std %>%
    select(gene, rank) %>%
    full_join(
      loess %>% select(gene, rank),
      by = "gene",
      suffix = c("_std", "_loess")
    ) %>%
    full_join(
      entropy %>% select(gene, rank),
      by = "gene"
    ) %>%
    rename(
      `Entropy` = rank,
      `SD` = rank_std,
      `LOESS` = rank_loess
    ) %>%
    column_to_rownames("gene")
  return(combined_ranks)
}



# ================================================
#   Inputs
# ================================================

studies <- readr::read_csv("resources/datasets.csv",
  progress = FALSE,
  show_col_types = FALSE
)


# Standard deviation
ranks_std <- tar_read(first_pc,
  store = "results/R/principal_components"
)

# Shannon Entropy
ranks_entropy <- tar_read(first_pc,
  store = "results/R/shannon_entropy"
)

# LOESS Residuals
ranks_loess <- tar_read(first_pc,
  store = "results/R/loess_residuals"
)

# ================================================
#   Pipeline
# ================================================

list(
  tar_target(
    combined_df,
    command = combine_ranks(ranks_std, ranks_entropy, ranks_loess)
  ),
  tar_target(
    ranks_corr,
    command = cor(combined_df, method = "spearman")
  ),
  tar_quarto(
    report,
    "analysis/07_measures_comparison.qmd"
  )
) %>%
  tar_hook_before(
    hook = conflicted::conflicts_prefer(
      dplyr::filter,
      base::setdiff,
      dplyr::rename
    )
  )
