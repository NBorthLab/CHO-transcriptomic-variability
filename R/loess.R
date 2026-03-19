get_expression_statistics <- function(normalized_dds) {
  expressions <- DESeq2::counts(normalized_dds, normalized = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression")

  dataset <- unique(colData(normalized_dds)$Dataset)

  expression_statistics <- expressions %>%
    group_by(gene) %>%
    summarize(
      dataset = dataset,
      mean = mean(expression),
      sd = sd(expression),
      median = median(expression),
      cv = sd(expression) / mean(expression)
    )

  return(expression_statistics)
}

# ------------------------------------------------

fit_loess <- function(expression_statistics) {
  fit <- loess(
    formula = log(cv) ~ log(mean),
    data = expression_statistics,
    degree = 1,
    span = 0.75
  )

  results <- list(
    expression_statistics = expression_statistics,
    fit = fit
  )

  return(results)
}

# ------------------------------------------------

get_loess_residuals <- function(fit_data) {
  expression_statistics <- fit_data$expression_statistics
  fit <- fit_data$fit

  residuals <- expression_statistics %>%
    mutate(
      metric = resid(fit),
      rank = rank(metric)
    )

  return(residuals)
}

# ------------------------------------------------

plot_loess_residuals <- function(data) {
  dataset <- unique(data$dataset)
  filename <- file.path(plotdir, str_glue("loess_residuals_{dataset}.pdf"))

  p1 <- data %>%
    tidyplot(x = mean, y = cv) %>%
    add_data_points(alpha = 0.2, size = 0.5) %>%
    add(geom_smooth(
      formula = y ~ x,
      color = "red", method = "loess", span = 0.75,
      method.args = list(degree = 1)
    )) %>%
    adjust_x_axis(transform = "log10", breaks = c(10, 100, 1000)) %>%
    adjust_y_axis(transform = "log2")

  p2 <- data %>%
    tidyplot(x = mean, y = metric) %>%
    add_data_points(alpha = 0.2, size = 0.5) %>%
    add(geom_smooth(
      formula = y ~ x,
      color = "red",
      method = "loess", span = 0.75,
      method.args = list(degree = 1)
    )) %>%
    adjust_x_axis(transform = "log10", breaks = c(10, 100, 1000))

  ggsave(filename,
    plot = p1 - p2
  )
  filename
}
