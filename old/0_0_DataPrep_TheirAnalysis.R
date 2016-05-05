setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")


# DATA IMPORT -------------------------------------------------------------
str(prot <- data.matrix(read.csv("data/proteins.csv", row.names = 1)))
str(gene <- data.matrix(read.csv("data/transcripts.csv", row.names = 1)))
str(ratio <- data.matrix(read.csv("data/ratios.csv", row.names = 1)))
library(gplots)
venn(list(p = row.names(prot), g= row.names(gene), r=row.names(ratio)))

colnames(prot)
colnames(gene)
colnames(gene) <- colnames(prot)
stopifnot(all(colnames(prot) == colnames(gene)))

stopifnot(all(rownames(prot) == rownames(gene)))
stopifnot(all(rownames(prot) == rownames(ratio)))

#########################
# ORIGINAL VALUES ---------------------------------------------------------
#########################
plot(density(apply(prot,1,function(x){return(mean(x, na.rm = T))}), na.rm=T))
sum(apply(prot,1,function(x){return(mean(x, na.rm = T))}) > 8.5, na.rm = T)


# Getting ratios ----------------------------------------------------------
log10(10**(prot[2,1]-10)/(10**gene[2,1]/10**6))
ratio[2,1]
##
log10(median(10**(prot[2,]-10)/(10**gene[2,]/10**6)))
##
geneUN <- (10**gene)/10**6
protUN <- 10**(prot-10)
round(log10(median(protUN[2,]/geneUN[2,])),7)== round(ratio[2,1],7)
round(log10(median(10**(prot[2,]-10)/(10**gene[2,]/10**6))),7) == round(ratio[2,1],7)

str(ratio <- data.frame(ratio))
str(ratio$nik <- log10(apply(protUN/geneUN, 1, function(x) median(x, na.rm=T))))
head(ratio)
table(ratio$tr <- (round(ratio$nik,5) == round(ratio[,1],5)))
subset(ratio, is.na(tr))
sum(is.na(ratio$tr))
sum(is.na(ratio$protein.mRNA.ratio))
ratio$rat <- ratio$protein.mRNA.ratio
ratio$rat[is.na(ratio$rat)] <- ratio$X[is.na(ratio$rat)]
head(ratio)

# Median protein vs ratio --------------------------------------------------
str(pDat <- data.frame(ratio=ratio$nik, prot=apply(prot,1,function(x){return(median(x, na.rm = T))})))
library(ggplot2)
ggplot(pDat, aes(x=prot, y=ratio)) + geom_point()
cor(pDat, use="pairwise.complete.obs")

# Predicting protein from RNA ---------------------------------------------
str(protPred <- geneUN * (10**ratio$protein.mRNA.ratio))
str(protPred <- log10(protPred)+10)

# Reproducing published results -------------------------------------------
tissueMats <- lapply(colnames(prot), function(tis) data.frame(p=prot[,tis], pr=protPred[,tis]))
names(tissueMats) <- colnames(prot)
# compare with figure SF7 or fig 5a bottom right (Salvary gland)
ggplot(tissueMats$uterus, aes(x=pr, y=p)) + geom_point() + theme_bw()
round(claimed <- sapply(tissueMats, function(mat) cor(mat, method="spearman", use="pairwise.complete.obs")[2,1]),2)

#########################
# WITHOUT OVERFITTING (same tissue not used to predict) ---------------------------------------------------------
#########################

# Matrices each missing one tissue ----------------------------------------
nikGenes <- lapply(colnames(geneUN), function(tis) geneUN[,-which(colnames(geneUN) == tis)])
names(nikGenes) <- colnames(geneUN)
nikProteins <- lapply(colnames(protUN), function(tis) protUN[,-which(colnames(protUN) == tis)])
names(nikProteins) <- colnames(protUN)
nikRatios <- lapply(colnames(geneUN), function(tis) apply((nikProteins[[tis]]/nikGenes[[tis]]), 1, function(x) median(x, na.rm=T)))
names(nikRatios) <- colnames(geneUN)
# Test this another way
str(ratioMat <- protUN/geneUN)
nikRatios2 <- lapply(colnames(ratioMat), function(tis) apply(ratioMat[,-which(colnames(ratioMat) == tis)], 1, function(x) median(x, na.rm=T)))
names(nikRatios2) <- colnames(ratioMat)
for(i in names(nikRatios)){ stopifnot(all(nikRatios[[i]] == nikRatios2[[i]], na.rm=T)) }

# Do one example to test it worked -----------------------------------------
str(nikGenesUterus <- geneUN[,-which(colnames(geneUN) == "uterus")])
str(nikProteinUterus <- protUN[,-which(colnames(protUN) == "uterus")])
colnames(nikGenesUterus)
stopifnot(all(colnames(nikGenesUterus) == colnames(nikProteinUterus)))
str(nikRatiosUterusMat <- nikProteinUterus/nikGenesUterus)
stopifnot(nikRatiosUterusMat[21, "prostate"] == protUN[21, "prostate"]/geneUN[21, "prostate"])
str(nikRatiosUterusVec <- apply(nikRatiosUterusMat, 1, function(row){return(median(row, na.rm=T))}))
stopifnot(all(nikRatiosUterusVec == nikRatios$uterus, na.rm=T))
stopifnot(all(nikRatiosUterusVec == nikRatios2$uterus, na.rm=T))
(x <- subset(data.frame(r = nikRatiosUterusVec, r1 = nikRatios$uterus, r2 = nikRatios$uterus, same = nikRatiosUterusVec == nikRatios$uterus), is.na(same)))
stopifnot(all(is.na(x$r1)))
stopifnot(all(is.na(x$r2)))

# Predict proteins --------------------------------------------------------
sapply(nikRatios, str)
str(nikRatiosMat <- do.call(cbind, nikRatios))
all(colnames(nikRatiosMat) == colnames(geneUN))
str(nikProtPred <- log10(geneUN * nikRatiosMat)+10)

# creating the same numbers and figures
nikTissueMats <- lapply(colnames(geneUN), function(tis) data.frame(p=prot[,tis], pr=nikProtPred[,tis]))
names(nikTissueMats) <- colnames(geneUN)
str(nikTissueMats)
round(sapply(nikTissueMats, function(mat) cor(mat, method="spearman", use="pairwise.complete.obs")[2,1]),2)
ggplot(nikTissueMats$uterus, aes(x=pr, y=p)) + geom_point() + theme_bw()

boxplot(list(
  claimed=claimed,
  overfit_corrected=sapply(nikTissueMats, function(mat) cor(mat, method="spearman", use="pairwise.complete.obs")[2,1])
  ), las = 1)

save(gene, geneUN, prot, protUN, protPred, nikRatiosMat, nikProtPred, claimed, nikTissueMats, file="0_0_data.RData")
