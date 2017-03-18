library(shiny)
source("create.gene.links.R", local = TRUE)

## Load summary data
gene.summary.data <- readRDS("data-cache/gene-summary-data.rds")

selected.example.genes <- c(
    "FN3K", # cor = 0.89
    "CRYL1", # cor = 0.41
    "TBCD", # cor = 0.03
    "PRDX5", # cor = 0
    "SRSF1", # cor = -0.41
    "HNRNPH2" # cor = -0.75
)



