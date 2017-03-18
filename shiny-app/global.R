library(shiny)
source("create.gene.links.R", local = TRUE)

## Load summary data
gene.summary.data <- readRDS("data-cache/gene-summary-data.rds")

selected.example.genes <- c(
    "FN3K",
    "SAMHD1",
    "ALDH2",
    "CORO1B",
    "PHYHD1",
    "IDH2"
)



