# DATA IMPORT -------------------------------------------------------------
str(prot <- data.matrix(read.csv("../data/proteins.csv", row.names = 1)))
str(gene <- data.matrix(read.csv("../data/transcripts.csv", row.names = 1)))
str(ratio <- data.matrix(read.csv("../data/ratios.csv", row.names = 1)))

library(gplots)
venn(list(p = row.names(prot), g= row.names(gene), r=row.names(ratio)))

colnames(prot)
colnames(gene)
colnames(gene) <- colnames(prot)
stopifnot(all(colnames(prot) == colnames(gene)))

stopifnot(all(rownames(prot) == rownames(gene)))
stopifnot(all(rownames(prot) == rownames(ratio)))

# THAT's WHAT THEY DID
geneUN <- (10**gene)/10**6
protUN <- 10**(prot-10)

# Some controls (do i get their ratio?):
stopifnot(
  round(log10(
    median(protUN[2,]/geneUN[2,])
    ),7) == round(ratio[2,1],7))
# just checking all this logging
stopifnot(
  round(log10(
    median(10**(prot[2,]-10)/(10**gene[2,]/10**6))
    ),7) == round(ratio[2,1],7))

# GET THE INFORMATION I FINALLY USE ----------------------------------------------------------------

tissues <- colnames(geneUN)
table(apply(!is.na(geneUN), 1, sum), apply(!is.na(protUN), 1, sum))
str(proteins <- rownames(geneUN)) # using all genes

mats <- list(
  mRNAs = geneUN[proteins, tissues],
  prots = protUN[proteins, tissues]
)

save(mats, tissues, proteins, file="1_Mats.RData")