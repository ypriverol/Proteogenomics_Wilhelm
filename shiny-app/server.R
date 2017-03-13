
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(dplyr)
source("ggplot-utils.R")

## Load data
prot.summary.data <- readRDS("data-cache/prot-summary-data.rds")
prot.individual.data <- readRDS("data-cache/prot-individual-data.rds")

## Pre-compute orders
avail.ord <- with(prot.summary.data, order(avail.both, -cor, prot))

prots.by.avail <- with(prot.summary.data,
                       split(prot[avail.ord],
                             factor(avail.both[avail.ord],
                                    levels = rev(seq(min(avail.both), max(avail.both))))))
prots.by.avail <- with(prot.summary.data,
                       split(prot[avail.ord],
                             factor(avail.both[avail.ord],
                                    levels = rev(seq(min(avail.both), max(avail.both))))))

cor.split <- cut(prot.summary.data$cor, breaks = 8, right = FALSE, include.lowest = T)
cor.split <- factor(cor.split, levels = rev(levels(cor.split)))
levels(cor.split) <- sub("^[\\(\\[]([-\\.0-9]{1,5}).+$", "≥ \\1", levels(cor.split))
cor.ord <- with(prot.summary.data, order(cor.split, -avail.both, -cor))
prots.by.cor <- split(prot.summary.data$prot[cor.ord], cor.split[cor.ord])

##
## Reactive server
##

shinyServer(function(input, output, session) {
    selected.prot <- reactive({
        if (!is.null(input$prot)) {
            input$prot
        } else {
            names(prot.individual.data)[[1]]
        }
    })

    observe({
        sorted.choices <- switch(input$`prot-sort`,
                                 "nr-samples" = prots.by.avail,
                                 "cor" = prots.by.cor,
                                 prot.summary.data$prot)
        updateSelectInput(session, "prot",
                          choices = sorted.choices,
                          selected = selected.prot())
    })

    output$gene.stat = renderText({
        indiv.stat <- prot.individual.data[[selected.prot()]]
        sprintf("Correlation = %.2f (based on %d available pairs)", indiv.stat$info["cor"],
                indiv.stat$info["avail.both"])
    })

    output$cors.hist <- renderPlot({
        fd <- filter(prot.summary.data,
               avail.both >= input$`summary-avail`[1L],
               avail.both <= input$`summary-avail`[2L])

        ggplot(fd, aes(x = cor)) +
            geom_histogram(bins = nclass.FD(na.omit(fd$cor)),
                           fill = "#0072B255", color = "#0072B2FF") +
            scale_x_continuous("Correlation", expand = c(.03, 0)) +
            scale_y_continuous("Count", expand = c(0, 0)) +
            coord_cartesian(xlim = c(-1, 1)) +
            ggplot_theme(base_size = 16)
    })

    output$prot.name <- renderText({selected.prot()})

    output$scatter <- renderPlot({
        indiv.stat <- prot.individual.data[[selected.prot()]]
        max.xy <- with(indiv.stat$data, c(max(mrna, na.rm = TRUE), max(prot, na.rm = TRUE)))

        ggplot(indiv.stat$data, aes(x = mrna, y = prot)) +
            geom_abline(intercept = 0, slope = indiv.stat$info["medr"],
                        linetype = "dashed", size = PLOT_LINE_WIDTH,
                        color = "#0072B2") +
            geom_point(size = PLOT_POINT_SIZE, color = "#0072B2", shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, color = "#0072B2", alpha = .4) +
            xlab("miRNA Expression") +
            ylab("Protein Expression") +
            coord_cartesian(xlim = c(0, max.xy[1]), ylim = c(0, max.xy[2])) +
            ggplot_theme(base_size = 16)
    })

    output$ratio <- renderPlot({
        indiv.stat <- prot.individual.data[[selected.prot()]]

        ggplot(indiv.stat$data, aes(x = tissue, y = ratio)) +
            geom_abline(intercept = indiv.stat$info["medr"], slope = 0,
                        linetype = "dashed", size = PLOT_LINE_WIDTH,
                        color = "#0072B2") +
            geom_point(size = PLOT_POINT_SIZE, color = "#0072B2", shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, color = "#0072B2", alpha = .4) +
            xlab("Tissue") +
            ylab("Ratio protein/miRNA") +
            ggplot_theme(base_size = 16) +
            theme(
                axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1)
            )
    })

    output$prediction <- renderPlot({
        indiv.stat <- prot.individual.data[[selected.prot()]]
        good.points <- with(indiv.stat$data, !is.na(prot) & !is.na(pred.prot))
        max.xy <- with(indiv.stat$data[good.points, ], max(c(prot, pred.prot), na.rm = TRUE))

        ggplot(indiv.stat$data, aes(x = prot, y = pred.prot)) +
            geom_abline(slope = 1,
                        linetype = "dotted", size = PLOT_LINE_WIDTH,
                        color = "#333333") +
            geom_point(size = PLOT_POINT_SIZE, color = "#D55E00", shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, color = "#D55E00", alpha = .4) +
            ylab("Predicted Protein Expression") +
            scale_x_continuous("Protein Expression") +
            coord_fixed(xlim = c(0, max.xy), ylim = c(0, max.xy)) +
            ggplot_theme(base_size = 16)
    })
})
