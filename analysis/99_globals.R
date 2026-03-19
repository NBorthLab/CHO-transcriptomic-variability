base_fontsize <- 8
tag_fontsize <- 10

theme_tidy <- function() {
  theme(
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
    plot.background = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.key = element_rect(fill = NA, colour = NA),
    strip.background = element_rect(fill = NA, colour = NA),
    panel.background = element_rect(fill = NA, colour = NA),
    # panel.border = element_rect(
    #   fill = NA,
    #   colour = "black",
    #   linewidth = 0.5
    # ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.ticks = element_line(colour = "black", linewidth = 0.25),
    panel.border = element_blank(),
    axis.line = element_line(linewidth = 0.25, color = "black"),
    axis.ticks = element_line(linewidth = 0.25, color = "black")
  )
}

theme_border <- function() {
  theme_tidy() +
    theme(
      panel.border = element_rect(
        fill = NA,
        colour = "black",
        linewidth = 0.5
      )
    )
}
