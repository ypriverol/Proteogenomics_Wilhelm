################################################################################
##
## Get gene symobls & names for ENSEMBL gene record IDs
##
################################################################################
library(org.Hs.eg.db)

all.genes <- readRDS("data-cache/prot-summary-data.rds")

matches <- select(
    org.Hs.eg.db,
    keys = all.genes$prot,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "GENENAME", "SYMBOL")
)

# Filter unmatched gene records
matches <- matches[!is.na(matches$SYMBOL), ]

# Handle duplicates
# --> Remove the gene records with less data
matches$ENSEMBL[with(matches, SYMBOL %in% SYMBOL[duplicated(SYMBOL)])]
remove.dupl <- c(
    "ENSG00000150526",
    "ENSG00000185414",
    "ENSG00000256671"
)

matches <- matches[!c(matches$ENSEMBL %in% remove.dupl), ]
saveRDS(matches, file = "data-cache/gene-names.rds")
