calculate_entropies <- function(exprs) {
  expr_matrix <- SummarizedExperiment::assay(exprs)

  specialization <- BioQC::sampleSpecialization(expr_matrix) # , norm = TRUE)
  SummarizedExperiment::colData(exprs)$specialization <- specialization

  specificity <- BioQC::entropySpecificity(expr_matrix) # , norm = TRUE)
  SummarizedExperiment::rowData(exprs)$specificity <- specificity

  diversity <- BioQC::entropyDiversity(expr_matrix) # , norm = TRUE)
  SummarizedExperiment::colData(exprs)$diversity <- diversity

  return(exprs)
}


calculate_entropy <- function(cpm) {
  cpm_matrix <- assay(cpm)
  specificity <- BioQC::entropySpecificity(cpm_matrix)
  dataset <- unique(colData(cpm)$Dataset)

  entropy <- data.frame(
    gene = names(specificity),
    metric = specificity,
    dataset = dataset
  ) %>%
    mutate(metric_rank = rank(specificity))

  return(entropy)
}

# source("scratch/47_weighted_entropy/RNentropy/R/RNentropy-internal.R")

weighted_relative_entropy <- function(cpm) {
  cpm_matrix <- assay(cpm)
  test <- apply(cpm_matrix, 1L, function(gene_expr) {
    cum_expr <- sum(gene_expr)
    n_samples <- length(gene_expr)
    mean_expr <- cum_expr / n_samples
    g_statistic <- sum(gene_expr * log(gene_expr / mean_expr))
  })

  test
}
