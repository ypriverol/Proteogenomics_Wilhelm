library(shiny)

## Load summary data
gene.summary.data <- readRDS("data-cache/gene-summary-data.rds")

selected.example.genes <- c(
    "DIABLO", #cor = 0.17
    "PPP1CC", #cor = 0.07
    "PPIA", # cor = -0.02
    "UFC1", # cor = -0.04
    "FN3K", # cor = 0.89
    "HNRNPH2" # cor = -0.75
    # "BTF3", # cor = 0.38
    # "TBCD", # cor = 0.03
    # "PRDX5", # cor = 0
    # "SRSF1", # cor = -0.41
)

selected.example.genes.with.cor <- with(gene.summary.data, {
    with(gene.summary.data[match(selected.example.genes, gene), ],
         sprintf("%s <span>Correlation = %.2f</span>", gene, cor))
})

