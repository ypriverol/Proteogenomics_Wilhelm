
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)
library(ggplot2)
library(dplyr)
source("ggplot-utils.R")

## Load individual data
# prot.summary.data <- readRDS("data-cache/prot-summary-data.rds")
prot.individual.data <- readRDS("data-cache/prot-individual-data.rds")

## Create list of genes, split by # of available observations
avail.ord <- with(prot.summary.data, order(avail.both, -cor, prot))

prots.by.avail <- with(prot.summary.data,
                       split(prot[avail.ord],
                             factor(avail.both[avail.ord],
                                    levels = rev(seq(min(avail.both), max(avail.both))))))

prots.by.avail$`0` <- NULL
prots.by.avail$`1` <- NULL

## Create list of genes, split by common prefix
prots.by.alpha <- split(prot.summary.data$prot, substr(prot.summary.data$prot, 1, 11))

## Pre-compute the links groups
prots.by.avail.links <- create.group.containers(prots.by.avail)
prots.by.alpha.links <- create.group.containers(prots.by.alpha)

## States of current gene list
MAX_NR_GENE_GROUPS <- 30
current.grouping <- prots.by.cor.links
current.grouping.plain <- prots.by.cor
current.group.observers <- NULL

##
##
##
renderIndividualPlots <- function(selected.gene, output) {
    output$prot.name <- renderText({selected.gene})

    output$gene.stat = renderText({
        indiv.stat <- prot.individual.data[[selected.gene]]
        sprintf("Correlation = %.2f (based on %d available pairs)", indiv.stat$info["cor"],
                indiv.stat$info["avail.both"])
    })

    output$scatter <- renderPlot({
        indiv.stat <- prot.individual.data[[selected.gene]]
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
        indiv.stat <- prot.individual.data[[selected.gene]]

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
        indiv.stat <- prot.individual.data[[selected.gene]]
        good.points <- with(indiv.stat$data, !is.na(prot) & !is.na(pred.prot))
        max.xy <- with(indiv.stat$data[good.points, ], max(c(prot, pred.prot), na.rm = TRUE))

        ggplot(indiv.stat$data, aes(y = prot, x = pred.prot)) +
            geom_abline(slope = 1,
                        linetype = "dotted", size = PLOT_LINE_WIDTH,
                        color = "#333333") +
            geom_point(size = PLOT_POINT_SIZE, color = "#D55E00", shape = 1) +
            geom_point(size = PLOT_POINT_SIZE, color = "#D55E00", alpha = .4) +
            ylab("Protein Expression") +
            scale_x_continuous("Predicted Protein Expression") +
            coord_fixed(xlim = c(0, max.xy), ylim = c(0, max.xy)) +
            ggplot_theme(base_size = 16)
    })
}

##
## Reactive server
##

shinyServer(function(input, output, session) {
    lapply(seq_len(MAX_NR_GENE_GROUPS), function (i) {
        observeEvent(input[[sprintf("gene-group-%d", i)]], {
            group.attribs <- current.grouping$containers$children[[i]]$children[[2]]$attribs

            # Toggle collapsed state
            current.grouping$containers$children[[i]]$children[[2]]$attribs$collapsed <<- !current.grouping$containers$children[[i]]$children[[2]]$attribs$collapsed
            if (!isTRUE(group.attribs$loaded)) {
                current.grouping$containers$children[[i]]$children[[2]]$attribs$loaded <<- TRUE
                insertUI(
                    sprintf("#%s", group.attribs$id),
                    where = "afterBegin",
                    ui = current.grouping$gene.links[[i]]
                )
            }

            if (!isTRUE(current.grouping$containers$children[[i]]$children[[2]]$attribs$collapsed)) {
                if (!is.null(current.group.observers)) {
                    lapply(current.group.observers, function(obs) {
                        obs$destroy()
                    })
                }
                current.group.observers <<- lapply(seq_along(current.grouping$gene.links[[i]]$children), function(j) {
                    observeEvent(input[[sprintf("gene-id-%d-%d", i, j)]], {
                        session$sendCustomMessage(type = "scrollCallback", 1)
                        renderIndividualPlots(current.grouping.plain[[i]][j], output)
                    })
                })
            }
        })
    })

    observe({
        current.grouping <<- switch(input$`gene-sort`,
                                    "nr-samples" = prots.by.avail.links,
                                    "cor" = prots.by.cor.links,
                                    prots.by.alpha.links)


        current.grouping.plain <<- switch(input$`gene-sort`,
                                          "nr-samples" = prots.by.avail,
                                          "cor" = prots.by.cor,
                                          prots.by.alpha)

        # Set all containers to "not loaded" and "collapsed"
        current.grouping$containers$children <- lapply(current.grouping$containers$children, function (node) {
            node$children[[2L]]$attribs$loaded <- FALSE
            node$children[[2L]]$attribs$collapsed <- TRUE
            return(node)
        })


        removeUI("#gene-links", session = session, immediate = TRUE)
        insertUI(
            "#gene-links-container",
            where = "afterBegin",
            ui = current.grouping$containers
        )
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

    ## The first gene that is shown is picked at random
    renderIndividualPlots(sample(prot.summary.data$prot, 1), output)
})
