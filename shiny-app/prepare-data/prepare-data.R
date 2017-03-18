################################################################################
##
## Pre-compute the values for the graphical output
##
################################################################################
library(dplyr)
library(tibble)

load("../../final/3_Mats.RData", envir = environment())
gene.name.map <- readRDS("../data-cache/gene-names.rds")

tissue.names <- c(
    testis = "Testis",
    stomach = "Stomach",
    salivary.gland = "Salivary gland",
    kidney = "Kidney",
    esophagus = "Oesophagus",
    adrenal.gland = "Adrenal gland",
    spleen = "Spleen",
    prostate = "Prostate",
    pancreas = "Pancreas",
    ovary = "Ovary",
    uterus = "Uterus",
    thyroid.gland = "Thyroid gland"
)
tissue.names <- factor(tissue.names, levels = tissue.names)

gene.names.for.prots <- with(gene.name.map, SYMBOL[match(proteins, ENSEMBL)])

indiv.data <- lapply(proteins[!is.na(gene.names.for.prots)], function (p) {
    ind <- match(p, proteins)

    ratio <- with(mats, prots[ind, ] / mRNAs[ind, ])
    medr <- median(ratio, na.rm = TRUE)

    ret.dat <- with(mats, data.frame(
        tissue = tissue.names[tissues],
        mrna = mRNAs[ind, ],
        prot = prots[ind, ],
        ratio = ratio,
        pred.prot = mRNAs[ind, ] * medr
    ))

    list(
        data = ret.dat,
        info = with(ret.dat, c(
            medr = medr,
            cor = cor(mrna, prot,
                      use = "pairwise.complete.obs",
                      method = "spearman"),
            avail.both = sum(is.finite(mrna) & is.finite(prot)),
            avail.mrna = sum(is.finite(mrna)),
            avail.prot = sum(is.finite(prot))
        ))
    )
})


names(indiv.data) <- na.omit(gene.names.for.prots)

saveRDS(indiv.data, "../data-cache/gene-individual-data.rds")

summary.data <- sapply(indiv.data, function (id) {
    c(id$info["avail.both"], id$info["cor"])
})

summary.data <- rownames_to_column(as.data.frame(t(summary.data)), var = "gene") %>%
    left_join(select(
        gene.name.map,
        gene = SYMBOL,
        ensembl = ENSEMBL,
        gene.name = GENENAME
    ), by = "gene")

saveRDS(summary.data, "../data-cache/gene-summary-data.rds")
