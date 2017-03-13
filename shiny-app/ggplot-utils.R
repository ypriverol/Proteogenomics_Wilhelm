PLOT_POINT_SIZE <- ggplot2::rel(2)
PLOT_LINE_WIDTH <- ggplot2::rel(.7)
PLOT_TEXT_SIZE <- ggplot2::rel(2)

ggplot_theme <- function(base_size = 12, base_family = "") {
    requireNamespace("ggplot2", quietly = TRUE)

    ggplot2::theme_bw(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = ggplot2::rel(1)),
            plot.margin = grid::unit(c(0.2, 0.4, 0.5, 0.2), "lines"),
            panel.background = ggplot2::element_rect(fill = NA, color = NA),
            plot.background = ggplot2::element_rect(fill = NA, color = NA),
            legend.title = ggplot2::element_text(size = ggplot2::rel(0.7)),
            legend.text = ggplot2::element_text(size = ggplot2::rel(0.7)),
            axis.title = ggplot2::element_text(size = ggplot2::rel(0.75)),
            axis.text = ggplot2::element_text(size = ggplot2::rel(0.7)),
            axis.title.y = ggplot2::element_text(vjust = 1),
            axis.title.x = ggplot2::element_text(vjust = 0),
            panel.grid.major = ggplot2::element_line(color = "gray30", size = ggplot2::rel(0.5),
                                                     linetype="dotted"),
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_rect(fill = "#ffffff", color = "gray50",
                                                     size = 0.3),
            strip.text = ggplot2::element_text(size = ggplot2::rel(0.75)),
            panel.border = ggplot2::element_rect(color = "gray50", size = 0.3),
            legend.background = ggplot2::element_blank()
        )
}
