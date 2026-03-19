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

plotdir <- "plots/gene_charactersistics"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}
resultsdir <- "results/gene_charactersistics"
if (!dir.exists(resultsdir)) {
  dir.create(resultsdir, recursive = TRUE)
}


# ================================================
#   Inputs
# ================================================

tar_load(first_pc, store = "results/R/principal_components")
tar_load(gtf_granges, store = "results/R/gene_features")
tar_load(txdb, store = "results/R/gene_features")


subset_decile_genes <- function(all_genes, first_pc) {
  #' Subset the GRanges genes object to genes in our analysis (in deciles)

  mcols(all_genes) <- left_join(
    data.frame(gene = all_genes$gene_id),
    first_pc,
    by = join_by(gene)
  )

  genes <- all_genes[!is.na(all_genes$decile)]

  genes
}

get_gene_flanks <- function(genes) {
  #' Extend gene bodies +- 1000bp

  genes_flanked <- genes %>%
    plyranges::anchor_center() %>%
    plyranges::stretch(2 * 1000)

  genes_flanked
}

get_decile_genes <- function(dec, genes) {
  genes[genes$decile == dec & !is.na(genes$decile)]
}


# ================================================
#   Pipeline
# ================================================


list(
  tar_target(
    genes,
    subset_decile_genes(txdb$genes, first_pc)
  ),
  tar_target(
    genes_flanked,
    get_gene_flanks(genes)
  ),
  per_decile_calc <- tar_map(
    values = list(decile = 1:10),

    # Get decile GRanges objects
    tar_target(
      genes_decile,
      get_decile_genes(decile, genes)
    ),
    tar_target(
      genes_flanked_dec,
      get_decile_genes(decile, genes_flanked)
    ),

    # Gene Length and neighbor distance
    tar_target(
      gene_length,
      GenomicDistributions::calcWidth(genes_decile)
    ),
    tar_target(
      gene_neigh_dist,
      GenomicDistributions::calcNeighborDist(genes_decile)
    ),

    # GC Content in gene body +- 1000bp
    tar_target(
      gc_content,
      GenomicDistributions::calcGCContent(
        genes_flanked_dec,
        BSgenome.Cgriseus.NCBI.CriGriPICRH1.0::BSgenome.Cgriseus.NCBI.CriGriPICRH1.0
      )
    )
  ),
  tar_combine(
    gene_lengths,
    per_decile_calc[["gene_length"]],
    command = list(!!!.x)
  ),
  tar_combine(
    gene_dists,
    per_decile_calc[["gene_neigh_dist"]],
    command = list(!!!.x)
  ),
  tar_combine(
    gc_contents,
    per_decile_calc[["gc_content"]],
    command = list(!!!.x)
  ),
  tar_quarto(
    report,
    "analysis/13_gene_characteristics.qmd"
  )
)
