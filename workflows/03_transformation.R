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

plotdir <- "plots/transformation"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}
resultsdir <- "results/transformation"
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
subset_dataset(coding_counts, "Dhiman")

list(
  per_study_metric <- tar_map(
    values = studies,
    names = "dataset",
    tar_target(
      dataset_counts,
      command = subset_dataset(coding_counts, dataset)
    ),
    tar_target(
      design,
      command = make_design(dataset_counts)
    ),
    tar_target(
      counts_filt,
      command = filter_genes(dataset_counts, design)
    ),
    tar_target(
      dds,
      command = normalize_deseq2(counts_filt, design)
    ),

    # Run blind and non-blind VST
    tar_map(
      values = list(blind = c("blind", "nonblind")),
      tar_target(
        vsd,
        command = transform_vst(counts_filt, design, blind)
      ),
      tar_target(
        std,
        command = get_std(vsd)
      )
    ),

    # Optionally, also try out rlog transformation
    # !! Very time consuming computation !!
    if (params$run_rlog) {
      list(
        tar_target(
          rld,
          command = transform_rlog(dds, design)
        ),
        tar_target(
          std_rld,
          command = get_std(rld)
        )
      )
    },
    NULL
  ),

  # Normalized counts
  tar_combine(
    ddss,
    per_study_metric[["dds"]],
    command = list(!!!.x)
  ),

  # Transformed counts
  tar_combine(
    vsds_blind,
    per_study_metric[["vsd_blind"]],
    command = list(!!!.x)
  ),
  tar_combine(
    vsds_nonblind,
    per_study_metric[["vsd_nonblind"]],
    command = list(!!!.x)
  ),

  # Standard deviations
  tar_combine(
    stds_blind,
    per_study_metric[["std_blind"]],
    command = bind_rows(!!!.x)
  ),
  tar_combine(
    stds_nonblind,
    per_study_metric[["std_nonblind"]],
    command = bind_rows(!!!.x)
  ),

  # Plots
  tar_target(
    dispersions_plot,
    command = plot_dispersions(ddss),
    format = "file"
  ),
  tar_target(
    mean_sd_plot,
    command = plot_mean_sds(vsds_blind, vsds_nonblind),
    format = "file"
  ),

  # Regularized log
  if (params$run_rlog) {
    list(
      tar_combine(
        rlds,
        per_study_metric[["rld"]],
        command = list(!!!.x)
      ),
      tar_combine(
        standard_deviations_rld,
        per_study_metric[["std_rld"]],
        command = bind_rows(!!!.x)
      )
    )
  },
  tar_quarto(
    report,
    "analysis/03_transformation.qmd",
    quiet = FALSE
  )
)
