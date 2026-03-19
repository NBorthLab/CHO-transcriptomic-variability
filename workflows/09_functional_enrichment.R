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

plotdir <- "plots/functional_enrichment"
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

# orthologs <- get_all_orthologs(first_pc$gene)
# mouse_genes <- convert_orthologs(stable_genes, orthologs)

list(
  tar_target(
    go_generic_slim_file,
    command = "resources/gene_ontology/goslim_generic.tsv",
    format = "file"
  ),
  tar_target(
    secRecon_xlsx,
    command = "resources/secRecon/Supplementary_Tables_S1-S3.xlsx",
    format = "file"
  ),
  tar_target(
    gene_list,
    command = first_pc %>% pull(PC1, name = gene) %>% sort(decreasing = TRUE)
  ),

  #
  # Full GO
  #
  tar_target(
    go_gsea,
    command = run_gseGO(gene_list)
  ),
  # KEGG
  tar_target(
    kegg_gsea,
    command = run_gseKEGG(gene_list)
  ),
  tar_target(
    kegg_gsea_df,
    command = as.data.frame(kegg_gsea)
  ),
  tar_target(
    kegg_gsea_csv,
    command = {
      file <- "results/functional_enrichment/gsea_kegg.csv"
      write_csv(kegg_gsea_df, file)
      file
    },
    format = "file"
  ),

  #
  # GO Slim
  #
  tar_target(
    go_slim,
    command = get_slim_annotation(go_generic_slim_file, first_pc$gene)
  ),
  tar_target(
    go_slim_gsea,
    command = run_gsea(gene_list, go_slim$term2gene, go_slim$term2name)
  ),
  tar_target(
    go_slim_gsea_df,
    command = as.data.frame(go_slim_gsea)
  ),
  tar_target(
    go_slim_gsea_csv,
    command = {
      file <- "results/functional_enrichment/gsea_goslim.csv"
      write_csv(go_slim_gsea_df, file)
      file
    },
    format = "file"
  ),

  #
  # secRecon
  #
  tar_target(
    secRecon,
    command = get_secRecon_ontology(secRecon_xlsx)
  ),
  tar_target(
    secRecon_gsea,
    command = run_gsea(gene_list, secRecon)
  ),
  tar_target(
    secRecon_gsea_df,
    command = as.data.frame(secRecon_gsea)
  ),
  tar_target(
    secRecon_gsea_csv,
    command = {
      file <- "results/functional_enrichment/gsea_secRecon.csv"
      write_csv(secRecon_gsea_df, file)
      file
    },
    format = "file"
  ),
  tar_quarto(
    report,
    "analysis/09_functional_enrichment.qmd"
  )
) %>%
  tar_hook_before(
    hook = conflicted::conflicts_prefer(
      dplyr::filter,
      dplyr::select,
      dplyr::rename
    )
  )
