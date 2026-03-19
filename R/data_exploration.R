plot_dispersions <- function(dds_list) {
  plotfile <- file.path(plotdir, "dispersions.pdf")

  datasets <- names(dds_list) %>% str_remove("dds_")

  pdf(plotfile)
  map2(dds_list, datasets, function(x, y) {
    plotDispEsts(
      x,
      main = y, xlim = c(1e1, 1e6), ylim = c(1e-8, 1e1)
    )
  })
  dev.off()

  return(plotfile)
}


plot_mean_sds <- function(vsds_blind, vsds_nonblind) {
  filename <- file.path(plotdir, "mean_sd_plot.pdf")
  # plots <- map2(vsd_list, datasets, function(x, y) {
  #   p <- vsn::meanSdPlot(assay(x), plot = FALSE)
  #   p$gg <- p$gg + labs(title = y, subtitle = x$Design)
  #   p$gg
  # })
  #
  # ggsave(filename,
  #   plot = gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)
  # )
  datasets <- names(vsds_blind) %>% str_remove("vsd_blind_")

  test <- pmap(list(vsds_blind, vsds_nonblind, datasets), function(b, nb, d) {
    pb <- vsn::meanSdPlot(assay(b), plot = FALSE)
    pb$gg <- pb$gg + labs(title = d, subtitle = "blind")
    pnb <- vsn::meanSdPlot(assay(nb), plot = FALSE)
    pnb$gg <- pnb$gg + labs(title = d, subtitle = "non-blind")
    list(pb$gg, pnb$gg)
  }) %>%
    flatten()

  ggsave(filename,
    plot = gridExtra::marrangeGrob(test, nrow = 1, ncol = 2),
    width = 12
  )

  return(filename)
}

plot_overall_pca <- function(raw_counts) {
  filename <- "plots/pca_overall.pdf"

  filtered_counts <- filter_genes2(raw_counts)
  vsd <- DESeq2::vst(
    DESeq2::DESeqDataSetFromMatrix(
      countData = assay(filtered_counts) %>% round(),
      colData = colData(filtered_counts),
      design = ~1
    ),
    blind = TRUE
  )
  data <- DESeq2::plotPCA(vsd,
    intgroup = "SampleName",
    returnData = TRUE
  )
  pca_data <- data %>%
    select(PC1, PC2) %>%
    cbind(colData(vsd))
  pca_data %>%
    mutate(celllib = interaction(CellLine, LibraryStrategy)) %>%
    tidyplot(x = PC1, y = PC2, color = Dataset) %>%
    add_data_points() %>%
    add(ggforce::geom_mark_ellipse(
      data = pca_data,
      mapping = aes(x = PC1, y = PC2, fill = Dataset, label = Dataset),
      expand = unit(2, "mm"),
      label.fontsize = 8,
      con.cap = 0
    )) %>%
    adjust_size(NA, NA) %>%
    save_plot(filename, width = 400, height = 300)
}

plot_pca <- function(vsd_list, filename) {
  plot_data <- map2_dfr(vsd_list, datasets, function(x, y) {
    # intgroup <- case_when(
    #   y == "Hofer" ~ c("Temperature", "CulturePhase")
    # )
    d <- DESeq2::plotPCA(x,
      intgroup = "FullSampleName",
      returnData = TRUE
    )
    vars <- attr(d, "percentVar")

    d <- mutate(d,
      dataset = y,
      var1 = vars[1],
      var2 = vars[2]
    )
    d
  })

  plot_data %>%
    tidyplot(x = PC1, y = PC2) %>%
    add_data_points() %>%
    add_data_labels_repel(label = FullSampleName) %>%
    add(coord_fixed()) %>%
    split_plot(dataset, nrow = 1, ncol = 1, widths = 150, heights = 150) %>%
    save_plot(filename, bg = "white")

  # ggsave(
  #   plotfile,
  #   plot = gridExtra::marrangeGrob(plots, nrow = 1, ncol = 1)
  # )

  return(filename)
}

plot_pairwise_sd <- function(sd_blind, sd_nonblind) {
  filename <- "plots/standard_deviations_blind_nonblind.pdf"

  full_join(sd_blind, sd_nonblind, join_by(gene, dataset)) %>%
    tidyplot(x = sd.x, y = sd.y) %>%
    add_data_points(rasterize_dpi = 100) %>%
    split_plot(dataset, width = 50, heights = 40) %>%
    save_plot(filename, bg = "white", view_plot = FALSE)

  return(filename)
}

plot_pairwise_sd_ranks <- function(sd_blind, sd_nonblind) {
  filename <- "plots/standard_deviations_ranks_blind_nonblind.pdf"

  full_join(sd_blind, sd_nonblind, join_by(gene, dataset)) %>%
    tidyplot(x = sd_rank.x, y = sd_rank.y) %>%
    add_data_points(rasterize_dpi = 100) %>%
    split_plot(dataset, width = 50, heights = 40) %>%
    save_plot(filename, bg = "white", view_plot = FALSE)

  return(filename)
}
