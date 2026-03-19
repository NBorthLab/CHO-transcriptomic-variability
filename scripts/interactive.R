library(tidyverse)
library(targets)
library(tarchetypes)
library(SummarizedExperiment)
library(shiny)
library(plotly)
library(esquisse)


studies <- readr::read_csv("resources/datasets.csv",
  progress = FALSE,
  show_col_types = FALSE
)

# ================================================
#   Entropies vs Std deviation
# ================================================

# Entropy vs std
tar_load(entropies, store = "results/R/entropy")
tar_load(stds, store = "results/R/variance")

combined_df <- full_join(entropies, stds, join_by(gene, dataset))

tar_load(cpms, store = "results/R/entropy/")

cpm_df <- cpms %>%
  map(assay) %>%
  map(as_tibble, rownames = "gene") %>%
  map(pivot_longer, -gene, names_to = "sample", values_to = "cpm") %>%
  map2(str_remove(names(cpms), "cpm_"), (\(x, y) mutate(x, dataset = y))) %>%
  bind_rows()

# LOESS vs std
tar_load(stds_nonblind, store = "results/R/transformation")
loess_resids <- tar_read(loess_resids, store = "results/R/loess_residuals") %>%
  dplyr::rename(
    resid_rank = "rank", resid = "metric", sd_counts = "sd"
  )

combined_df <- full_join(stds_nonblind, loess_resids, join_by(gene, dataset))

tar_load(ddss, store = "results/R/transformation")

counts_df <- ddss %>%
  map(counts, normalized = TRUE) %>%
  map(as_tibble, rownames = "gene") %>%
  map(pivot_longer, -gene, names_to = "sample", values_to = "counts") %>%
  map2(str_remove(names(ddss), "dds_"), (\(x, y) mutate(x, dataset = y))) %>%
  bind_rows()


ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectizeInput("gene", "Gene", choices = NULL),
      selectInput("dataset", "Dataset", studies$dataset, "Dhiman")
    ),
    mainPanel(
      plotOutput("counts_plot"),
      plotlyOutput("rank_plot"),
      plotlyOutput("values_plot")
    )
  )
)

server <- function(input, output, session) {
  updateSelectizeInput(session, "gene",
    choices = unique(counts_df$gene), server = TRUE
  )

  output$counts_plot <- renderPlot({
    counts_df %>%
      filter(dataset == input$dataset) %>%
      filter(gene %in% input$gene) %>%
      ggplot(aes(x = sample, y = counts)) +
      geom_point()
  })

  # output$cpm_plot <- renderPlot({
  #   cpm_df %>%
  #     filter(dataset == input$dataset) %>%
  #     filter(gene %in% input$gene) %>%
  #     ggplot(aes(x = sample, y = cpm)) +
  #     geom_point()
  # })

  output$values_plot <- renderPlotly({
    highligh <- combined_df %>%
      filter(gene == input$gene) %>%
      filter(dataset == input$dataset)
    p <- combined_df %>%
      # mutate(tooltip = str_c(gene, sd, specificity, g, sep = "_")) %>%
      filter(dataset == input$dataset) %>%
      ggplot(aes(x = sd, y = resid, label = gene)) +
      geom_point(alpha = 0.5) +
      geom_point(
        data = highligh,
        mapping = aes(sd, resid),
        color = "red",
        size = 3
      )

    ggplotly(p)
  })

  output$rank_plot <- renderPlotly({
    highligh <- combined_df %>%
      filter(gene == input$gene) %>%
      filter(dataset == input$dataset)
    p <- combined_df %>%
      filter(dataset == input$dataset) %>%
      ggplot(aes(x = sd_rank, y = resid_rank, label = gene)) +
      geom_point(alpha = 0.5) +
      geom_point(
        data = highligh,
        mapping = aes(sd_rank, resid_rank),
        color = "red",
        size = 2
      )

    ggplotly(p)
  })
}

shinyApp(ui, server)


# ================================================
#   AnnotationHub interactive
# ================================================

BiocHubsShiny::BiocHubsShiny()
