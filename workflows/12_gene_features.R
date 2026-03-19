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

plotdir <- "plots/gene_features"
if (!dir.exists(plotdir)) {
  dir.create(plotdir, recursive = TRUE)
}
resultsdir <- "results/gene_features"
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

get_gtf_annotation <- function(gtf) {
  rtracklayer::import(gtf)
}

get_genefeatures <- function(gtf_granges) {
  txdb <- txdbmaker::makeTxDbFromGRanges(gtf_granges)
  txdb_genes <- GenomicFeatures::genes(txdb)
  # txdb_promoters <- promoters(txdb, columns = c("gene_id"))
  txdb_promoters <- GenomicFeatures::promoters(txdb_genes)
  txdb_exons <- GenomicFeatures::exonsBy(txdb, by = "gene")

  # Subset to only primary transcript of a gene
  # i.e. transcript with most cpm expressoin
  tx2gene <- get_primary_transcript()
  txdb_transcripts <- GenomicFeatures::transcripts(txdb,
    columns = c("gene_id", "tx_name")
  ) %>%
    plyranges::filter(tx_name %in% tx2gene$transcript_id)

  list(
    genes = txdb_genes,
    promoters = txdb_promoters,
    exons = txdb_exons,
    transcripts = txdb_transcripts
  )
}

get_chromosome_alias <- function(alias_tsv) {
  chromosome_id_df <- read_tsv(alias_tsv,
    col_names = FALSE
  ) %>%
    column_to_rownames("X1")
  return(chromosome_id_df)
}

read_revised_tss <- function(revised_tss_tsv, chromosome_alias) {
  tss <- read_tsv(revised_tss_tsv,
    col_names = c(
      "chr", "start", "end", "name", "stat", "strand", "tss_id", "gene",
      "confident", "open_chromatin_CHO", "open_chromatin_tissue"
    )
  ) %>%
    dplyr::filter(stat > 0) %>%
    # filter(confident) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    plyranges::mutate(start = end(.)) %>%
    plyranges::mutate(
      tss_number = str_replace(name, "p(.)@.*", "\\1"),
      tss_tx = str_replace(name, "p.@[A-Za-z0-9]+_(.*)", "\\1")
    )

  GenomeInfoDb::seqlevels(tss) <-
    chromosome_alias[GenomeInfoDb::seqlevels(tss), ]

  return(tss)
}
# revised_tss <- read_revised_tss(revised_tss_tsv, chromosome_alias)

# FIXME: Make dependent on isoform.results input
get_primary_transcript <- function() {
  transcript_counts <- read_tsv(
    "results/counts/Stor/Stor01.isoforms.results"
  )
  tx2gene <- transcript_counts %>%
    group_by(gene_id) %>%
    arrange(desc(TPM), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::select(gene_id, transcript_id) %>%
    ungroup() %>%
    mutate(gene_id = str_remove(gene_id, "_[0-9]"))
  return(tx2gene)
}

subset_tss <- function(tss) {
  # Revised TSS is from PICR, need to convert transcript IDs to PICRH
  compare_prev <- read_tsv(
    "resources/genome/GCF_003668045.3_CriGri-PICRH-1.0_compare_prev.txt.gz"
  )
  picr2picrh <- compare_prev %>%
    dplyr::select(
      `previous transcript accession`, `current transcript accession`
    ) %>%
    dplyr::rename(
      PICR = `previous transcript accession`,
      PICRH = `current transcript accession`
    ) %>%
    dplyr::filter(!is.na(PICR)) %>%
    pull(PICRH, name = PICR)

  tss$picrh_txid <- picr2picrh[tss$tss_tx]

  # Subset to only primary transcript of a gene
  primary_tx <- get_primary_transcript()

  filtered_tss <- tss %>%
    `[`(tss$picrh_txid %in% primary_tx$transcript_id)

  return(filtered_tss)
}
# filtered_tss <- subset_tss(revised_tss)

get_decile_tss <- function(dec, revised_tss, genes) {
  gene_strands_df <- data.frame(
    gene = as.character(genes$gene_id),
    strand = GenomicRanges::strand(genes)
  ) %>%
    column_to_rownames("gene")

  #' Overwrite TSS strand information with gene's strand
  correct_strand <- function(tss) {
    strand(tss) <- gene_strands_df[tss$gene, ]
    return(tss)
  }

  dec_genes <- first_pc %>%
    dplyr::filter(decile == dec) %>%
    pull(gene)
  decile_tss <- revised_tss[revised_tss$gene %in% dec_genes] %>%
    correct_strand() %>%
    plyranges::filter(tss_number == 1)

  return(decile_tss)
}

# get_decile_tss(1, filtered_tss, txdb_genes)

# Homer BED file

# tar_load(decile_tss)
# decile_tss[[1]] %>%
#   plyranges::mutate(start = end(.) - 1) %>%
#   plyranges::select(name) %>%
#   write_bed("/dev/shm/homer/tss_var.bed")

write_decile_bed <- function(dec_tss, dec) {
  file <- str_glue("{resultsdir}/{dec}_tss.bed")
  dec_tss %>%
    # BED files automatically add that or something
    # plyranges::mutate(start = end(.) - 1) %>%
    plyranges::select(name) %>%
    plyranges::write_bed(file)
  return(file)
}

# END

get_tss_flanks <- function(tss) {
  tss %>%
    plyranges::anchor_center() %>%
    plyranges::stretch(100)
}

get_sequence <- function(regions) {
  sequences <- memes::get_sequence(
    regions,
    BSgenome.Cgriseus.NCBI.CriGriPICRH1.0::BSgenome.Cgriseus.NCBI.CriGriPICRH1.0
  )
  old_names <- names(sequences)
  new_names <- paste0(regions$name, "__", old_names)
  names(sequences) <- new_names
  return(sequences)
}

write_fasta <- function(seqs, decile) {
  file <- stringr::str_glue("{resultsdir}/{decile}_tss_flanks.fa")
  memes::write_fasta(seqs, file)
  return(file)
}

list(
  tar_target(
    gtf,
    "resources/chromatin_states/annotation_aliased.gtf",
    format = "file"
  ),
  tar_target(
    alias_tsv,
    "resources/chromatin_states/PICRH_aliases_rev.txt",
    format = "file"
  ),
  tar_target(
    revised_tss_tsv,
    "resources/exact_tss/SD2_revisedTSS_RefSeq_PICRH_sort.bed",
    format = "file"
  ),
  tar_target(
    gtf_granges,
    get_gtf_annotation(gtf)
  ),
  tar_target(
    chromosome_alias,
    get_chromosome_alias(alias_tsv)
  ),
  tar_target(
    revised_tss,
    read_revised_tss(revised_tss_tsv, chromosome_alias)
  ),
  tar_target(
    filtered_tss,
    subset_tss(revised_tss)
  ),
  tar_target(
    txdb,
    get_genefeatures(gtf_granges)
  ),
  per_decile_tss <- tar_map(
    values = list(decile = 1:10),
    tar_target(
      tss,
      get_decile_tss(decile, revised_tss, txdb$genes)
    ),
    tar_target(
      tss_bed,
      write_decile_bed(tss, decile),
      format = "file"
    ),
    tar_target(
      tss_flanks,
      get_tss_flanks(tss)
    ),
    tar_target(
      tss_flanks_seqs,
      get_sequence(tss_flanks)
    ),
    tar_target(
      tss_flanks_seqs_fasta,
      write_fasta(tss_flanks_seqs, decile),
      format = "file"
    )
  ),
  tar_combine(
    decile_tss_flanks,
    per_decile_tss[["tss_flanks"]],
    command = GenomicRanges::GRangesList(!!!.x)
  ),
  tar_combine(
    decile_tss,
    per_decile_tss[["tss"]],
    command = GenomicRanges::GRangesList(!!!.x)
  ),
  tar_target(
    all_tss,
    unlist(decile_tss)
  ),
  tar_target(
    revised_tss_bed,
    write_decile_bed(revised_tss, "all")
  )
)
