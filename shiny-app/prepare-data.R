################################################################################
##
## Pre-compute the values for the graphical output
##
################################################################################
library(tibble)

load("../data//3_Mats.RData", envir = environment())

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

indiv.data <- lapply(proteins, function (p) {
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
            # cor.pred = with(ret.dat, cor(pred.prot, prot,
            #                              use = "pairwise.complete.obs",
            #                              method = "spearman")),
            cor = cor(mrna, prot,
                      use = "pairwise.complete.obs",
                      method = "spearman"),
            avail.both = sum(is.finite(mrna) & is.finite(prot)),
            avail.mrna = sum(is.finite(mrna)),
            avail.prot = sum(is.finite(prot))
        ))
    )
})

names(indiv.data) <- proteins

saveRDS(indiv.data, "data-cache/prot-individual-data.rds")


summary.data <- sapply(indiv.data, function (id) {
    c(id$info["avail.both"], id$info["cor"])
})

summary.data <- rownames_to_column(as.data.frame(t(summary.data)), var = "prot")

saveRDS(summary.data, "data-cache/prot-summary-data.rds")
