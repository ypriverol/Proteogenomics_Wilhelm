
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(type = "text/javascript", src = "search-genes.js")
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
        tabPanel("Example Genes", sidebarLayout(
            sidebarPanel(
                width = 0.25 * 12,
                h4("Select Gene"),
                div(
                    class = "list-group",
                    id = "select-gene-nav",
                    lapply(selected.example.genes, function (gene) {
                        actionLink(
                            sprintf("gene-sel-%s", gene),
                            label = gene,
                            icon = icon("arrow-circle-right", class = "fa-lg"),
                            class = "list-group-item"
                        )
                    })
                ),
                selectizeInput(
                    "lookup-gene",
                    label = NULL,
                    width = "100%",
                    options = list(
                        placeholder = "Search Genes",
                        diacritics = FALSE,
                        maxOptions = 20,
                        searchField = c("gene", "geneName", "ensembl"),
                        valueField = "gene",
                        labelField = "gene",
                        render = I('{option: renderGene}')
                    ),
                    choices = NULL
                )
            ),
            mainPanel(
                width = 0.75 * 12,
                h3("Gene ", textOutput("gene.name", container = em)),
                helpText(textOutput("gene.stat")),
                splitLayout(
                    div(
                        h4("Scatterplot of Measured Expression Levels"),
                        helpText("The line depicts the predicted protein expression by $$\\hat{\\text{protein}}_i = \\hat{r} \\cdot \\text{miRNA}_i$$"),
                        plotOutput("scatter", width = "100%")
                    ),
                    div(
                        h4("Measured vs. Predicted Protein Expression"),
                        helpText("The line depicts the predicted protein expression by $$\\hat{\\text{protein}}_i = \\hat{r} \\cdot \\text{miRNA}_i$$"),
                        plotOutput("prediction", width = "100%"),
                        style = "text-align: center;",
                        class = "hide-help"
                    ),
                    cellArgs = list(style = "padding: 1em")
                ),
                splitLayout(
                    div(
                        h4(HTML("Ratio of Protein to miRNA")),
                        withMathJax(helpText("The line depicts the median ratio: $$\\hat{r} = \\operatorname{Med}_{i = 1, \\dotsc, n} (\\text{protein}_i / \\text{miRNA}_i)$$")),
                        plotOutput("ratio", width = "100%")
                    ),
                    div(
                        h4(" "),
                        plotOutput("legend", width = "100%"),
                        style = "text-align: center;"
                    ),
                    cellArgs = list(style = "padding: 1em")
                )
            )
        ))
    )
))
