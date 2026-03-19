get_std <- function(vsd) {
  mtx <- assay(vsd)
  sds <- rowSds(mtx)
  dataset <- unique(colData(vsd)$Dataset)

  gene_ranks <- data.frame(
    gene = names(sds),
    metric = sds,
    dataset = dataset
  ) %>%
    mutate(
      metric_rank = rank(metric)
    )

  return(gene_ranks)
}
