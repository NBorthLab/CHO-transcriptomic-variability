#' Get only genes from the annotation
#'
read_annotation <- function(gtf) {
  genome_annotation <- rtracklayer::readGFF(gtf)
  is_relevant_transcript <- genome_annotation$gene_biotype %in% c(
    "protein_coding"
  )
  genome_annotation[is_relevant_transcript, ]
}

read_counts_tables <- function(counts_rds) {
  counts <- readRDS(counts_rds)
  counts$countsFromAbundance <- NULL
  return(counts)
}

filter_samples <- function(counts) {
  # Filter out samples that were low quality from MultiQC reports
  # Samples are described in Lab notebook entry from 11.08.2025

  to_remove <- c(
    "Chiang04", "Chiang10", "Malm11", "Malm15", "Malm13", "NovakLnc04",
    "Orellana09", "Orellana08", "Rucker01"
  )

  counts <- counts[, !colnames(counts) %in% to_remove]

  return(counts)
}

filter_protein_coding_genes <- function(counts, annotation) {
  protein_coding_genes <- annotation$gene
  cleaned_gene_symbols <- str_replace(rownames(counts), "(.*)_1", "\\1")
  rownames(counts) <- cleaned_gene_symbols

  is_protein_coding <- rownames(counts) %in% protein_coding_genes
  counts <- counts[is_protein_coding, ]

  return(counts)
}

get_metadata <- function(metadata_csv) {
  conflicts_prefer(dplyr::filter)
  read_csv(metadata_csv) %>%
    suppressMessages() %>%
    # Condense technical replicates
    select(-RunAccession) %>%
    unique() %>%
    # Alphabetical sorting by dataset
    arrange(SampleName) %>%
    # This sample had an error in the RNA-seq pipeline
    filter(SampleName != "Malm19") %>%
    # Numeric metadata to factors
    mutate(
      Temperature = factor(Temperature),
      GlutamineConcentration = factor(GlutamineConcentration),
      Serum = factor(Serum)
    )
}

perform_batch_correction <- function(raw_counts) {
  # raw_counts <- tar_read(summarized_experiment)
  counts_matrix <- assay(raw_counts, "counts")
  metadata <- colData(raw_counts)
  batches <- factor(metadata$DatasetName)

  corr_counts <- sva::ComBat_seq(
    counts = counts_matrix,
    batch = batches,
    full_mod = FALSE
  )

  SummarizedExperiment(
    assay = corr_counts,
    colData = metadata
  )
}

remove_batch_effect <- function(norm_logcpm) {
  # norm_logcpm <- tar_read(getmm_logcpm_counts_uncorr)
  expr_matrix <- assay(norm_logcpm)
  metadata <- colData(norm_logcpm)
  batches <- factor(metadata$LibraryStrategy)
  batches2 <- factor(metadata$DatasetName)
  expr_corr <- removeBatchEffect(expr_matrix,
    batch = batches
  )
  SummarizedExperiment(
    assay = expr_corr,
    colData = metadata
  )
}

subset_dataset <- function(se, dataset) {
  dataset_samples <- se$Dataset == dataset
  return(se[, dataset_samples])
}

make_summarized_experiment <- function(counts_tables, metadata) {
  SummarizedExperiment(
    assay = list(
      counts = counts_tables$counts,
      tpm = counts_tables$abundance,
      length = counts_tables$length
    ),
    colData = metadata
  )
}

filter_genes5 <- function(se) {
  is_expressed <- assay(se) %>%
    `>`(10) %>%
    rowSums() %>%
    `>`(100)

  se[is_expressed, ]
}

filter_genes4 <- function(se) {
  assay(se) %>%
    `>`(10) %>%
    rowSums() %>%
    `>`(0.5 * ncol(se))
}

filter_genes3 <- function(se) {
  # Filter low expressed genes out
  cpm <- edgeR::cpm(assay(se))

  # Min 0.1 CPM in all samples and min average of 0.5 CPM
  pass_min_expr <- rowSums(cpm < 1) == 0
  # pass_avg_expr <- rowMeans(cpm) > 0.5

  # table(pass_min_expr)
  # table(pass_min_expr & pass_avg_expr)

  se <- se[pass_min_expr, ] # & pass_avg_expr, ]
  return(se)
}

filter_genes2 <- function(se) {
  # Filter low expressed genes out
  is_kept <- edgeR::filterByExpr(se)
  se <- se[is_kept, ]
  return(se)
}

filter_genes1 <- function(se, annotation) {
  # Filter out duplicated genes (only keeping genes with "_1") based on the GTF
  # annotation
  rownames(se) <- str_replace(rownames(se), "(.*)_1", "\\1")
  protein_coding_genes <- annotation$gene %>%
    `[`(annotation$gene_biotype == "protein_coding")
  is_relevant <- rownames(se) %in% protein_coding_genes
  se <- se[is_relevant, ]

  # Filter low expressed genes out
  is_kept <- filterByExpr(
    se,
    min.count = 10,
    large.n = 10,
    min.prop = 0.05
  )

  se <- se[is_kept, ]

  # Filter by counts
  # is_expressed <- SummarizedExperiment::assays(summarized_experiment) %>%
  #   `$`(counts) %>%
  #   `>`(10) %>%
  #   rowSums() %>%
  #   `>`(0)
  # summarized_experiment[is_expressed, ]
  return(se)
}
