
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(
            '
            Shiny.addCustomMessageHandler("scrollCallback", function (msg) {
                window.scrollTo(0, 0);
            });'
        )
    ),

    titlePanel("Alleluia"),

    tabsetPanel(
        tabPanel("Summary",
            h3("Summary of all Genes"),
            p("Correlations between the measured protein expression and the predicted protein expression"),
            plotOutput("cors.hist", width = "75%"),
            sliderInput("summary-avail", "Number of available observations per gene:",
                        min = 2, max = 12, value = c(8, 12),
                        width = "40%"
            )
        ),
        tabPanel("Protein Specific", sidebarLayout(
            sidebarPanel(
                width = 0.25 * 12,
                h4("Gene Selection"),
                selectizeInput(
                    "gene-sort",
                    label = "Group by:",
                    choices = c(
                        "Correlation" = "cor",
                        "# of samples" = "nr-samples",
                        "Alphanumeric" = "abc"
                    )
                ),
                div(
                    id = "gene-links-container",
                    prots.by.cor.links$containers
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
                        h4("Measured vs. Predicted Protein Expression"),
                        plotOutput("prediction", width = "100%"),
                        style = "text-align: center;"
                    ),
                    cellArgs = list(style = "padding: 1em")
                )
            )
        ))
    )
))
