library(shiny)
source("create.gene.links.R", local = TRUE)

## Load summary data
prot.summary.data <- readRDS("data-cache/prot-summary-data.rds")

## Create list of genes, split by correlation
cor.split <- cut(prot.summary.data$cor, breaks = 20, right = FALSE, include.lowest = T)
# table(cor.split)
cor.split <- factor(cor.split, levels = rev(levels(cor.split)))
levels(cor.split) <- sub("^[\\(\\[]([-\\.0-9]{1,5}).+$", "≥ \\1", levels(cor.split))
cor.ord <- with(prot.summary.data, order(cor.split, -avail.both, -cor))
prots.by.cor <- split(prot.summary.data$prot[cor.ord], cor.split[cor.ord])


## Pre-compute the links groups
prots.by.cor.links <- create.group.containers(prots.by.cor)
