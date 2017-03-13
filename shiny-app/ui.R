
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

## Load data
prot.summary.data <- readRDS("data-cache/prot-summary-data.rds")
cor.split <- cut(prot.summary.data$cor, breaks = 8, right = FALSE, include.lowest = T)
cor.split <- factor(cor.split, levels = rev(levels(cor.split)))
levels(cor.split) <- sub("^[\\(\\[]([-\\.0-9]{1,5}).+$", "≥ \\1", levels(cor.split))
cor.ord <- with(prot.summary.data, order(cor.split, -avail.both, -cor))
prots.by.cor <- split(prot.summary.data$prot[cor.ord], cor.split[cor.ord])


shinyUI(fluidPage(
    titlePanel("Alleluia"),

    tabsetPanel(
        tabPanel("Summary",
            h3("Summary of all Genes"),
            p("Correlations between the Protein Expression and the predicted Protein Expression"),
            plotOutput("cors.hist", width = "75%"),
            sliderInput("summary-avail", "Number of Available Observations per Gene:",
                        min = 2, max = 12, value = c(8, 12),
                        width = "40%"
            )
        ),
        tabPanel("Protein Specific", sidebarLayout(
            sidebarPanel(
                width = 0.25 * 12,
                h4("Gene Selection"),
                selectInput(
                    "prot-sort",
                    label = "Sorting:",
                    choices = c(
                        "Correlation" = "cor",
                        "# of samples" = "nr-samples",
                        "Alphanumeric" = "abc"
                    )
                ),
                selectInput(
                    "prot",
                    label = NULL,
                    size = 30,
                    selectize = FALSE,
                    choices = prots.by.cor
                )
            ),
            mainPanel(
                width = 0.75 * 12,
                h3("Gene ", textOutput("prot.name", container = em)),
                helpText(textOutput("gene.stat")),
                splitLayout(
                    div(
                        h4(HTML("Ratio of Protein to miRNA")),
                        withMathJax(helpText("The line depicts the median ratio: $$\\hat{r} = \\operatorname{Med}_{i = 1, \\dotsc, n} (\\text{protein}_i / \\text{miRNA}_i)$$")),
                        plotOutput("ratio", width = "100%")
                    ),
                    div(
                        h4("Scatterplot of Measured Expression Levels"),
                        helpText("The line depicts the predicted protein expression by $$\\hat{\\text{protein}}_i = \\hat{r} \\cdot \\text{miRNA}_i$$"),
                        plotOutput("scatter", width = "100%")
                    ),
                    cellArgs = list(style = "padding: 1em")
                ),
                splitLayout(
                    div(
                        h4("Predicted vs. Measured Protein Expression"),
                        plotOutput("prediction", width = "100%"),
                        style = "text-align: center;"
                    ),
                    cellArgs = list(style = "padding: 1em")
                )
            )
        ))
    )
))
