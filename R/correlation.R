# stds <- tar_read(entropies, store = "results/R/entropy")

filter_metric <- function(metric) {
  metric_matrix <- metric %>%
    select(gene, dataset, metric) %>%
    pivot_wider(names_from = "dataset", values_from = "metric") %>%
    column_to_rownames("gene") %>%
    as.matrix()

  has_NA <- is.na(metric_matrix) %>%
    rowSums() %>%
    `>`(0)
  metric_matrix <- metric_matrix[!has_NA, ]

  metric_matrix
}

correlate_metric <- function(metric) {
  corr_matrix <- cor(
    metric,
    use = "complete.obs",
    method = "spearman"
  )
  return(corr_matrix)
}

plot_correlation_matrix <- function(corr_matrix) {
  filename <- file.path(plotdir, "sd_rank_correlations.pdf")

  pdf(filename)
  corrplot::corrplot(
    corr_matrix,
    method = "ellipse",
    type = "lower",
    diag = FALSE
  )
  dev.off()

  # hist(corr_matrix, breaks = 50)

  filename
}

# corr_matrix <- tar_read(cor)
plot_correlation_histogram <- function(corr_matrix) {
  filename <- file.path(plotdir, "correlation_histogram.pdf")

  data <- corr_matrix %>%
    mtx2tibble() %>%
    rename(dataset1 = "gene", dataset2 = "samples") %>%
    filter(values < 0.99)

  p <- data %>%
    gghistogram(
      "values",
      bins = 40,
      add = "median",
      xlab = "Spearman correlation",
      rug = TRUE
    )

  ggsave(filename, plot = p)

  filename
}
