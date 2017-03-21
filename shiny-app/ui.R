
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

tag.class <- function (add.class, tag = span) {
    return(function(..., class = "") {
        el <- tag(
            ...,
            class = paste(add.class, class)
        )
        return(el)
    })
}


shinyUI(fluidPage(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(type = "text/javascript", src = "search-genes.js")
    ),

    titlePanel("Central Dogma"),

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
                    mapply(function (gene.desc, gene) {
                        actionLink(
                            sprintf("gene-sel-%s", gene),
                            label = HTML(gene.desc),
                            icon = icon("arrow-circle-right", class = "fa-lg"),
                            class = "list-group-item"
                        )
                    }, selected.example.genes.with.cor, selected.example.genes,
                    SIMPLIFY = FALSE)
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
                h3(
                    "Gene ",
                    textOutput("gene.name", container = tag.class("gene-name-title", em)),
                    textOutput("gene.accnr", container = tag.class("acc-nr"))),
                helpText(textOutput("gene.stat")),
                splitLayout(
                    div(
                        h4(
                            "Measured Levels",
                            textOutput("gene.cor.1", container = tag.class("cor-output"))
                        ),
                        helpText("The line depicts the predicted protein expression by $$\\hat{\\text{protein}}_i = \\hat{r} \\cdot \\text{miRNA}_i$$"),
                        plotOutput("scatter", width = "100%")
                    ),
                    div(
                        h4(
                            "Predictions",
                            textOutput("gene.cor.2", container = tag.class("cor-output"))
                        ),
                        helpText("The line depicts the predicted protein expression by $$\\hat{\\text{protein}}_i = \\hat{r} \\cdot \\text{mRNA}_i$$"),
                        plotOutput("prediction", width = "100%"),
                        style = "text-align: center;",
                        class = "hide-help"
                    ),
                    cellArgs = list(style = "padding: 1em", class = "plot-panel")
                ),
                splitLayout(
                    div(
                        h4("Ratio of Protein to mRNA"),
                        withMathJax(helpText("The line depicts the median ratio: $$\\hat{r} = \\operatorname{Med}_{i = 1, \\dotsc, n} (\\text{protein}_i / \\text{mRNA}_i)$$")),
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
        )),
        tabPanel("About", div(id = "about-tab",
                 h3("About"),
                 p("…"),
                 h4("Authors"),
                 tags$ul(
                     tags$li("Gabriela Cohen-Freue",
                             span(a("gcohen at stat.ubc.ca", href = "mailto:gcohen@stat.ubc.ca"))),
                     tags$li("…")
                 ),
                 h4("Acknowledgement"),
                 p(HTML(
                     'This shiny app is maintained and developed by ',
                     '<a href="mailto:d.kepplinger@stat.ubc.ca">David Kepplinger</a>,',
                     'using some parts and inspirations from ',
                     '<a href="http://neuroexpresso.org">neuroexpresso.org</a>',
                     'by <a href="https://github.com/oganm/">Ogan Mancarci</a>.'
                 ))
        ))
    )
))
