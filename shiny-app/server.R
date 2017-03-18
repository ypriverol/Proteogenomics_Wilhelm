
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(dplyr)
source("ggplot-utils.R")

OUT_IMG_PPI <- 144
OUT_IMG_BASESIZE <- 12
PLOT_POINT_SIZE <- ggplot2::rel(3)
PLOT_LINE_WIDTH <- ggplot2::rel(.7)
PLOT_TEXT_SIZE <- ggplot2::rel(3)


## Load individual data
# gene.name.map <- readRDS("data-cache/gene-names.rds")
gene.individual.data <- readRDS("data-cache/gene-individual-data.rds")

## Make gene search more efficient
gene.search.choices <- gene.summary.data %>%
    arrange(desc(avail.both), desc(cor)) %>%
    mutate(
        avail.both = as.character(avail.both),
        cor = sprintf("%.2f", cor)
    ) %>%
    rename(
        availBoth = avail.both,
        geneName = gene.name
    ) %>%
    as.data.frame()

##
## Color-scheme
##
color.palette <- c(
    "Testis" = "#56B4E9",
    "Stomach" = "#0072B2",
    "Salivary gland" = "#001833",
    "Kidney" = "#F0E442",
    "Oesophagus" = "#E69F00",
    "Adrenal gland" = "#D55E00",
    "Spleen" = "#747311",
    "Prostate" = "#009E73",
    "Pancreas" = "#2D7411",
    "Ovary" = "#2A4113",
    "Uterus" = "#CC79A7",
    "Thyroid gland" = "#742343"
)

## Pre-render legend
legend.plot <- data.frame(
    tissue = names(color.palette),
    col = rep(c(0, 1), times = 6),
    row = rev(rep(seq_len(length(color.palette) %/% 2L), each = 2))
) %>%
    ggplot(aes(x = col, y = row, color = tissue, label = tissue)) +
    geom_point(size = PLOT_POINT_SIZE, shape = 1) +
    geom_point(size = PLOT_POINT_SIZE, alpha = .4) +
    geom_text(size = PLOT_TEXT_SIZE, hjust = 0, nudge_x = .1) +
    scale_color_manual(values = color.palette, guide = "none") +
    # scale_x_continuous(expand = c(0, 0), limits = c(-.8, 2)) +
    coord_cartesian(xlim = c(-0.5, 2), ylim = c(-1, 8)) +
    theme_void() +
    theme(
        plot.background = ggplot2::element_rect(fill = NA, color = NA)
    )
# print(legend.plot)
# ggsave("www/legend.png", width = 3.5, height = 3, dpi = 300, bg = "transparent")

##
##
##
renderIndividualPlots <- function(selected.gene, output) {
    output$gene.name <- renderText({selected.gene})

    output$gene.stat = renderText({
        indiv.stat <- gene.individual.data[[selected.gene]]
        sprintf("Correlation = %.2f (based on %d available pairs)", indiv.stat$info["cor"],
                indiv.stat$info["avail.both"])
    })

    output$gene.cor.1 <- output$gene.cor.2 <- renderText({
        indiv.stat <- gene.individual.data[[selected.gene]]
        sprintf("%.2f", indiv.stat$info["cor"])
    })

    output$scatter <- renderPlot({
        indiv.stat <- gene.individual.data[[selected.gene]]
        max.xy <- with(indiv.stat$data, c(max(mrna, na.rm = TRUE), max(prot, na.rm = TRUE)))

        ggplot(indiv.stat$data, aes(x = mrna, y = prot, color = tissue)) +
            geom_abline(intercept = 0, slope = indiv.stat$info["medr"],
                        linetype = "dashed", size = PLOT_LINE_WIDTH,
                        color = "#0072B2") +
            geom_point(size = PLOT_POINT_SIZE, shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, alpha = .4) +
            scale_color_manual(values = color.palette, guide = "none") +
            xlab("miRNA Expression") +
            ylab("Protein Expression") +
            coord_cartesian(xlim = c(0, max.xy[1]), ylim = c(0, max.xy[2])) +
            ggplot_theme(base_size = OUT_IMG_BASESIZE)
    }, res = OUT_IMG_PPI)

    output$ratio <- renderPlot({
        indiv.stat <- gene.individual.data[[selected.gene]]

        ggplot(indiv.stat$data, aes(x = tissue, y = ratio, color = tissue)) +
            geom_abline(intercept = indiv.stat$info["medr"], slope = 0,
                        linetype = "dashed", size = PLOT_LINE_WIDTH,
                        color = "#0072B2") +
            geom_point(size = PLOT_POINT_SIZE, shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, alpha = .4) +
            scale_color_manual(values = color.palette, guide = "none") +
            xlab("Tissue") +
            ylab("Ratio protein/miRNA") +
            ggplot_theme(base_size = OUT_IMG_BASESIZE) +
            theme(
                axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
            )
    }, res = OUT_IMG_PPI)

    output$prediction <- renderPlot({
        indiv.stat <- gene.individual.data[[selected.gene]]
        good.points <- with(indiv.stat$data, !is.na(prot) & !is.na(pred.prot))
        max.xy <- with(indiv.stat$data[good.points, ], max(c(prot, pred.prot), na.rm = TRUE))

        ggplot(indiv.stat$data, aes(y = prot, x = pred.prot, color = tissue)) +
            geom_abline(slope = 1,
                        linetype = "dotted", size = PLOT_LINE_WIDTH,
                        color = "#333333") +
            geom_point(size = PLOT_POINT_SIZE, shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, alpha = .4) +
            scale_color_manual(values = color.palette, guide = "none") +
            ylab("Protein Expression") +
            scale_x_continuous("Predicted Protein Expression") +
            coord_fixed(xlim = c(0, max.xy), ylim = c(0, max.xy)) +
            ggplot_theme(base_size = OUT_IMG_BASESIZE)
    }, res = OUT_IMG_PPI)
}

##
## Reactive server
##

shinyServer(function(input, output, session) {
    updateSelectizeInput(
        session,
        'lookup-gene',
        server = TRUE,
        choices = gene.search.choices
    )

    ## Reacte to link-clicks
    reactive.vals <- reactiveValues(
        sel.gene = gene.search.choices$gene[1L]
    )
    lapply(selected.example.genes, function (gene) {
        observeEvent(input[[sprintf("gene-sel-%s", gene)]], {
            reactive.vals$sel.gene <- gene
            updateSelectizeInput(
                session,
                'lookup-gene',
                selected = NULL
            )
        })
    })

    ## Determine the gene selected in the search box or by a click on a link
    observe({
        dd.sel <- input$`lookup-gene`
        if (!is.null(dd.sel) & nchar(dd.sel) > 0L) {
            reactive.vals$sel.gene <- dd.sel
        }
    })
    # selected.gene <- reactive({
    #     dd.sel <- input$`lookup-gene`
    #     if (nchar(dd.sel) == 0L || is.null(dd.sel)) {
    #         dd.sel <- gene.search.choices$gene[1L]
    #     }
    #     return(dd.sel)
    # })

    ##
    observe({
        # sel.gene <- selected.gene()
        sel.gene <- reactive.vals$sel.gene
        renderIndividualPlots(sel.gene, output)
    })

    output$legend <- renderPlot(legend.plot, res = OUT_IMG_PPI)

    output$cors.hist <- renderPlot({
        fd <- filter(gene.summary.data,
               avail.both >= input$`summary-avail`[1L],
               avail.both <= input$`summary-avail`[2L])

        ggplot(fd, aes(x = cor)) +
            geom_histogram(bins = nclass.FD(na.omit(fd$cor)),
                           fill = "#0072B255", color = "#0072B2FF") +
            scale_x_continuous("Correlation", expand = c(.03, 0)) +
            scale_y_continuous("Count", expand = c(0, 0)) +
            coord_cartesian(xlim = c(-1, 1)) +
            ggplot_theme(base_size = OUT_IMG_BASESIZE)
    }, res = OUT_IMG_PPI)

    ## The first gene that is shown is picked at random
    renderIndividualPlots(sample(gene.summary.data$gene, 1), output)
})
