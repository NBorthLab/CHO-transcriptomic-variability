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

plotdir <- "plots/gene_sets"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}


tar_load(first_pc, store = "results/R/principal_components")
tar_load(principal_components, store = "results/R/principal_components")


list(
  tar_target(
    housekeeping_genes_csv,
    command = "resources/gene_sets/housekeeping_genes_Joshi.csv",
    format = "file"
  ),
  tar_target(
    housekeeping_genes,
    command = get_housekeeping_genes(housekeeping_genes_csv)
  ),
  tar_target(
    brown_housekeeping_genes, # et al., not the colour
    command = c(
      "Prkar1a", "Fkbp1a", "Mmadhc", "Gnb1", "Tmed2", "Actb", "Pgam1", "Gapdh"
    )
  ),
  tar_target(
    essential_genes_csv,
    command = "resources/gene_sets/xiong_2021_essential_genes.csv",
    format = "file"
  ),
  tar_target(
    essential_genes,
    command = get_essential_genes(essential_genes_csv)
  ),
  tar_target(
    geneset_data,
    command = get_geneset_data(
      first_pc, housekeeping_genes, essential_genes, brown_housekeeping_genes
    )
  ),

  # Skewness
  tar_target(
    hk_skewness,
    command = calculate_skewness(first_pc, housekeeping_genes)
  ),
  tar_target(
    ess_skewness,
    command = calculate_skewness(first_pc, essential_genes)
  ),

  # Consensus housekeepers
  tar_target(
    consensus_hk,
    command = geneset_data %>%
      filter(housekeeper & decile == 10) %>%
      pull(gene)
  ),
  tar_target(
    consensus_hk_gobp,
    command = enrich_GO_terms(consensus_hk, first_pc$gene, "BP")
  ),
  tar_target(
    consensus_hk_csv,
    command = {
      file <- "results/gene_sets/consensus_housekeepinggenes.csv"
      dir.create("results/gene_sets")
      write_csv(data.frame(gene = consensus_hk), file)
      file
    },
    format = "file"
  ),

  # Variable essential genes
  tar_target(
    essential_variable,
    command = geneset_data %>%
      filter(essential & decile == 1) %>%
      pull(gene)
  ),
  tar_target(
    essential_genes_go,
    command = enrich_GO_terms(essential_variable, first_pc$gene, "BP")
  ),

  # Report and Figures
  tar_quarto(
    report,
    "analysis/08_gene_sets.qmd"
  )
) %>%
  tar_hook_before(
    hook = conflicted::conflicts_prefer(
      dplyr::filter
    )
  )
