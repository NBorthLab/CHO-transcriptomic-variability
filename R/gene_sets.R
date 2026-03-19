get_housekeeping_genes <- function(csv) {
  hk_genes <- read_csv(csv)$Genes

  # WARN: This dataset has some issues:
  #   - The genes were saved in an Excel which lead to the common problem of
  #     date conversion...
  #   - The gene symbols are in all caps which is AGAINST the Chinese hamster
  #     gene nomenclature.
  #  Infuriating... really...

  hk_genes <- hk_genes %>%
    str_to_sentence()

  return(hk_genes)
}

plot_housekeeper_histogram <- function(first_pc, hk_genes) {
  filename <- file.path(plotdir, "housekeepers_histogram.pdf")

  data <- first_pc %>%
    arrange(rank) %>%
    mutate(housekeeper = toupper(gene) %in% hk_genes)

  break_pos <- cumsum(ntile(data$rank, 10) %>% table())

  data %>%
    filter(housekeeper) %>%
    tidyplot(x = rank) %>%
    add(geom_histogram(
      breaks = c(1, break_pos),
      alpha = 0.5,
      color = "black"
    )) %>%
    add(scale_x_reverse(
      limits = c(10000, 1),
      breaks = c(1, 2500, 5000, 7500, 10000)
    )) %>%
    adjust_x_axis_title("Gene expression variability rank") %>%
    adjust_y_axis_title("Count") %>%
    adjust_size(120, 120) %>%
    adjust_font(fontsize = 14) %>%
    save_plot(filename)

  filename
}

plot_housekeeper_pca <- function(pcs, hk_genes) {
  filename <- file.path(plotdir, "housekeepers_pca.pdf")

  principal_components %>%
    as_tibble(rownames = "gene") %>%
    select(gene, V1, V2) %>%
    mutate(housekeeper = toupper(gene) %in% hk_genes) %>%
    tidyplot(x = V1, y = V2) %>%
    add_reference_lines(x = 0, y = 0, color = "grey40") %>%
    add(
      geom_density2d(color = "grey70")
    ) %>%
    add(
      geom_density2d(
        data = filter_rows(housekeeper),
        color = "red"
      )
    ) %>%
    adjust_x_axis_title("PC1") %>%
    adjust_y_axis_title("PC2") %>%
    add(coord_fixed(xlim = c(-17000, 20000), ylim = c(-14000, 14000))) %>%
    remove_x_axis_labels() %>%
    remove_y_axis_labels() %>%
    adjust_size(120, 120) %>%
    adjust_font(fontsize = 14) %>%
    save_plot(filename)

  filename
}


get_geneset_data <- function(
    first_pc, housekeeping_genes, essential_genes, brown_housekeeping_genes) {
  first_pc %>%
    mutate(gene = gene) %>%
    left_join(
      data.frame(gene = housekeeping_genes, housekeeper = TRUE),
      by = join_by(gene)
    ) %>%
    left_join(
      data.frame(gene = essential_genes, essential = TRUE),
      by = join_by(gene)
    ) %>%
    left_join(
      data.frame(gene = brown_housekeeping_genes, brown_hk = TRUE),
      by = join_by(gene)
    ) %>%
    mutate(
      hk_rank = ifelse(housekeeper, rank, NA),
      ess_rank = ifelse(essential, rank, NA),
      brown_hk_rank = ifelse(brown_hk, rank, NA)
    )
}

calculate_skewness <- function(first_pc, gene_set) {
  first_pc %>%
    mutate(gene = toupper(gene)) %>%
    right_join(data.frame(gene = toupper({{ gene_set }})), by = "gene") %>%
    filter(!is.na(rank)) %>%
    pull(decile) %>%
    moments::skewness()
}


get_essential_genes <- function(essential_genes_csv) {
  read_tsv(essential_genes_csv) %>%
    mutate(gene = str_remove(gene, "__[0-9]")) %>%
    pull(gene) %>%
    str_split("\\|") %>%
    unlist()
}
