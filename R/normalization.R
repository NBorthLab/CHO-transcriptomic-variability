#' Normalize by DESeq2's Median of Ratios method followed by
#' Variance-stabilizing transformation.
#'
transform_vst <- function(counts, design, blind = c("blind", "nonblind")) {
  dds <- make_dds(counts, design)
  # Run VST
  vsd <- DESeq2::vst(dds, blind = blind == "blind")

  # Optionally, run batch effect correction
  if (any(!is.na(dds$Batch))) {
    corr_mat <- limma::removeBatchEffect(
      assay(vsd),
      batch = dds$Batch, design = design
    )
    assay(vsd) <- corr_mat
  }

  return(vsd)
}

transform_rlog <- function(dds, design, blind = c("blind", "nonblind")) {
  # Run VST
  rld <- DESeq2::rlog(dds, blind = blind == "blind")

  # Optionally, run batch effect correction
  if (any(!is.na(dds$Batch))) {
    corr_mat <- limma::removeBatchEffect(
      assay(rld),
      batch = dds$Batch, design = design
    )
    assay(rld) <- corr_mat
  }

  return(rld)
}

make_design <- function(counts) {
  if (any(!is.na(counts$Design))) {
    metadata <- colData(counts) %>% as.data.frame()
    metadata <- mutate(
      metadata,
      Temperature = droplevels(Temperature),
      GlutamineConcentration = droplevels(GlutamineConcentration),
      Serum = droplevels(Serum)
    )
    design <- model.matrix(
      paste("~ 0 +", unique(counts$Design)) %>% as.formula(),
      metadata
    )
    # if (!limma::is.fullrank(design)) {
    #   stop(str_glue("Design matrix not fullrank"))
    # }
  } else {
    design <- model.matrix(~1, data = colData(counts))
  }

  return(design)
}

filter_genes <- function(se, design) {
  # Filter low expressed genes out
  is_kept <- edgeR::filterByExpr(se, design)
  se <- se[is_kept, ]
  return(se)
}


filter_cooks <- function(counts, design) {
  norm_counts <- normalize_deseq2(counts, design)
  cutoff <- qf(0.7, ncol(design), ncol(norm_counts) - ncol(design))
  is_outlier <- assay(norm_counts, "cooks") > cutoff
  outlier_genes <- which(is_outlier, arr.ind = TRUE) %>% rownames()
  if (length(outlier_genes) != 0) {
    counts <- counts[!rownames(counts) %in% outlier_genes, ]
  }
  return(counts)
}

make_dds <- function(counts, design) {
  # Create DESeqDataSet
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = assay(counts) %>% round(),
    colData = colData(counts),
    design = design
  )
}

normalize_deseq2 <- function(counts, design, fit = "parametric") {
  dds <- make_dds(counts, design)
  dds <- DESeq2::DESeq(dds, fitType = fit)
  return(dds)
}

calculate_rpk2 <- function(se) {
  # counts <- tar_read(counts_corrected)
  # lengths <- tar_read(effective_gene_lengths)
  cnts <- assay(se, "counts")
  lens <- assay(se, "length")

  rpk <- cnts / (lens / 1000)

  NA_gene <- which(is.na(rpk) | is.infinite(rpk), arr.ind = TRUE) %>%
    `[`(, "row")

  if (length(NA_gene) != 0) {
    rpk <- rpk[-NA_gene, ]
  }

  SummarizedExperiment(
    assay = rpk,
    colData = colData(se)
  )
}


calculate_rpk <- function(counts, lengths) {
  # counts <- tar_read(counts_corrected)
  # lengths <- tar_read(effective_gene_lengths)
  cnts <- assay(counts)
  lens <- assay(lengths)

  # Some genes have zero lengths in samples. Probably due to short transcript
  # lengths than the fragment size of the RNA-seq experiment. Those are only a
  # few genes (~300) so we filter them out.
  has_zero_len <- rowSums(lens == 0) > 0

  cnts <- cnts[!has_zero_len, ]
  lens <- lens[!has_zero_len, ]

  rpk <- cnts / (lens / 1000)

  SummarizedExperiment(
    assay = rpk,
    colData = colData(counts)
  )
}

#' Normalize by Trimmed Means of M-values
#'
normalize_tmm <- function(se) {
  counts_mtx <- assay(se)

  dge <- edgeR::DGEList(counts = counts_mtx, samples = colData(se))
  norm_dge <- edgeR::calcNormFactors(dge)

  return(norm_dge)
}

get_cpm <- function(dge, log = FALSE) {
  norm_cpm <- edgeR::cpm(dge, log = log)

  SummarizedExperiment(
    assay = norm_cpm,
    colData = dge$samples
  )
}

# get_norm_cpm <- function(norm_counts) {
#
# }


#' Normalize by Gene-length corrected Trimmed Means of M-values.
#' Using the (!!!)Effective Length(!!!) as computed by RSEM.
#'
normalize_getmm <- function(summarized_experiment) {
  count_matrix <- assay(summarized_experiment, "counts") %>% round()
  length_matrix <- assay(summarized_experiment, "length")

  # Some genes have zero lengths in samples. Probably due to short transcript
  # lengths than the fragment size of the RNA-seq experiment. Those are only a
  # few genes (~300) so we filter them out.
  is_zero_len <- rowSums(length_matrix == 0) > 0

  cnts <- count_matrix[!is_zero_len, ]
  lens <- length_matrix[!is_zero_len, ]

  # Calculate Read per Kilobase of the (effective) read length
  rpk_matrix <- cnts / (lens / 1000)

  rpk_se <- SummarizedExperiment(
    assay = list(counts = rpk_matrix),
    colData = colData(summarized_experiment)
  )

  rpk_dge <- SE2DGEList(rpk_se)

  # TMM normalization from gene-length corrected values
  norm_cnts <- calcNormFactors(rpk_dge)

  # Calculate CPM and log2CPM
  norm_cpm <- cpm(norm_cnts, log = FALSE)
  norm_logcpm <- cpm(norm_cnts, log = TRUE)

  results_se <- SummarizedExperiment(
    assay = list(counts = norm_cnts, cpm = norm_cpm, logcpm = norm_logcpm),
    colData = colData(summarized_experiment)
  )

  return(results_se)
}
