#' Calculate eigen decomposition for a correlation matrix
#'
#' Compute the eigen decomposition of a square correlation matrix and
#' normalize the sign/orientation of the eigenvectors in a deterministic
#' way. This function returns the result of `eigen()`; if the entries of
#' the first eigenvector meet a simple orientation check the set of
#' eigenvectors is multiplied by -1 (to keep the sign convention
#' consistent across datasets).
#'
#' @param correlation_matrix Numeric square matrix. A correlation or
#'   covariance-like matrix for which the eigen decomposition should be
#'   computed.
#' @return A list with components `values` and `vectors` as returned by
#'   `eigen()`; the `vectors` component may have been multiplied by
#'   -1 to enforce a consistent orientation.
#' @examples
#' m <- diag(3)
#' calculate_eigenvectors(m)
calculate_eigenvectors <- function(correlation_matrix) {
  corr_eigen <- eigen(correlation_matrix)
  if (all(corr_eigen$vectors[, 1] < 1)) {
    corr_eigen$vectors <- -corr_eigen$vectors
  }
  return(corr_eigen)
}

#' Project centered ranks onto eigenvectors to obtain principal components
#'
#' Center the input rank matrix by subtracting column means (no scaling)
#' and project the centered data onto the supplied eigenvectors to obtain
#' principal component scores.
#'
#' @param metric_ranks Numeric matrix-like object with rows representing
#'   observations (genes) and columns representing variables (datasets).
#'   Row names will be preserved in the returned object when applicable.
#' @param eigenvectors Numeric matrix of eigenvectors (columns are
#'   eigenvectors) produced by `eigen()` or `calculate_eigenvectors()`.
#' @return A numeric matrix of principal component scores (observations
#'   x components).
#' @examples
#' mat <- matrix(rnorm(30), nrow = 10)
#' eig <- eigen(cov(mat))$vectors
#' get_principal_components(mat, eig)
get_principal_components <- function(metric_ranks, eigenvectors) {
  # metric_ranks <- metric_ranks %>%
  #   select(gene, dataset, rank) %>%
  #   pivot_wider(names_from = "dataset", values_from = "rank") %>%
  #   column_to_rownames("gene") %>%
  #   as.matrix()

  centered_ranks <- scale(metric_ranks, center = TRUE, scale = FALSE)
  principal_components <- centered_ranks %*% eigenvectors

  return(principal_components)
}

#' Calculate PCA loadings from eigen decomposition
#'
#' Given an eigen decomposition (as returned by `eigen()`), compute the
#' loadings matrix for principal component analysis. Loadings are
#' calculated as the eigenvectors multiplied by the square root of the
#' corresponding eigenvalues.
#'
#' @param eigen A list with components `vectors` and `values` (an
#'   eigen decomposition object).
#' @return A numeric matrix of loadings (variables x components).
#' @examples
#' e <- eigen(diag(3))
#' calculate_loadings(e)
calculate_loadings <- function(eigen) {
  eigen$vectors %*% diag(sqrt(eigen$values))
}

#' Extract and summarize the first principal component
#'
#' Convert a principal components matrix into a tibble containing the
#' first principal component (`PC1`) per observation (gene), and add a
#' descending rank and decile grouping. Rows with `NA` in `PC1` are
#' removed.
#'
#' @param principal_components Numeric matrix or object coercible to a
#'   data frame where the first component corresponds to column `V1`.
#'   Row names should contain gene identifiers.
#' @return A tibble with columns `gene`, `PC1`, `rank` (descending), and
#'   `decile` (1..10).
#' @examples
#' pc <- matrix(rnorm(30), nrow = 10)
#' rownames(pc) <- paste0("gene", seq_len(nrow(pc)))
#' get_first_component(pc)
get_first_component <- function(principal_components) {
  first_pc <- principal_components %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::select(gene, V1) %>%
    dplyr::rename(PC1 = "V1") %>%
    dplyr::filter(!is.na(PC1)) %>%
    dplyr::arrange(PC1) %>%
    dplyr::mutate(
      rank = rank(-PC1),
      decile = ntile(rank, 10)
    )

  return(first_pc)
}

#' Create a scree plot PDF for eigenvalues
#'
#' Produce and save a scree plot (percentage variance explained per
#' principal component) to a PDF file in the global `plotdir`. The
#' function expects `plotdir` to be defined in the environment where it
#' is called.
#'
#' @param eigen An eigen decomposition list with component `values`.
#' @return The full path to the saved PDF file (invisible side effect
#'   is the saved plot).
#' @note This function relies on plotting helper functions from the
#'   project (e.g. `tidyplot`, `add_sum_bar`, `save_plot`) and on the
#'   global `plotdir` variable.
#' @examples
#' # requires plotdir and plotting helpers to be available
#' # plot_scree(eigen(cov(matrix(rnorm(100), ncol=5))))
plot_scree <- function(eigen) {
  filename <- file.path(plotdir, "scree_plot.pdf")
  pdata <- data.frame(
    pc = seq_along(eigen$values),
    vars = (eigen$values / sum(eigen$values)) * 100
  ) %>%
    mutate(
      cum = cumsum(vars)
    )

  pdata %>%
    tidyplot(x = pc, y = vars) %>%
    add_sum_bar(saturation = 0.7) %>%
    add_sum_line(linewidth = 1) %>%
    adjust_x_axis(
      breaks = seq(1, 21)
    ) %>%
    adjust_y_axis(
      breaks = seq(0, 30, 5),
      limits = c(0, 30),
      padding = c(0, 0.05)
    ) %>%
    adjust_size(120, 100) %>%
    adjust_font(fontsize = 9) %>%
    save_plot(filename)

  filename
}

#' Plot principal component scores (PC1 vs PC2) and save to PDF
#'
#' Create a scatter plot of the first two principal components and save
#' the figure as a PDF in `plotdir`. Axis titles include the percent
#' variance explained by each component derived from `eigen$values`.
#'
#' @param principal_components Numeric matrix or data frame where the
#'   first two components are in columns `V1` and `V2`. Row names are
#'   treated as observation labels (genes).
#' @param eigen Eigen decomposition list (with `values`) used to
#'   compute percent variance annotations for axis titles.
#' @return The full path to the saved PDF file.
#' @examples
#' # requires plotdir and plotting helpers to be available
#' # plot_pc(pc_matrix, eigen(cov(matrix(rnorm(100), ncol=5))))
plot_pc <- function(principal_components, eigen) {
  filename <- file.path(plotdir, "ranks_pca.pdf")

  var <- round((eigen$values / sum(eigen$values)) * 100, 2)
  principal_components %>%
    as_tibble(rownames = "gene") %>%
    tidyplot(x = V1, y = V2) %>%
    add_data_points(alpha = 0.3, size = 0.4) %>%
    add(coord_fixed()) %>%
    adjust_x_axis(title = str_glue("PC1 ({var[1]}%)")) %>%
    adjust_y_axis(title = str_glue("PC2 ({var[2]}%)")) %>%
    adjust_size(120, 100) %>%
    adjust_font(fontsize = 9) %>%
    save_plot(filename)

  filename
}
