library(gdata)
library(ggplot2)
library(reshape)
library(reshape2)

load("1_Mats.RData")

str(mDat <- melt(mats$mRNAs))
str(pDat <- melt(mats$prots))
str(modelDat <- merge(mDat, pDat, by=c("X1", "X2")))
colnames(modelDat) <- c("Gene", "Tissue", "mRNA", "Protein")
genes <- proteins

# Get NAs
table(na.count <- sapply(genes, function(p){
  return(sum(is.na(mats$prots[p,]) | is.na(mats$mRNA[p,])))
}))

# which genes to look at?
str(genes <- genes[which(na.count <= 9)])
# str(genes <- genes[1:150]) # for a first test, just a subset of genes
str(modelDat <- drop.levels(subset(modelDat, Gene %in% genes)))
str(modelDat <- drop.levels(subset(modelDat, !is.na(Protein) & !is.na(mRNA))))

# # SLOPE ONLY --------------------------------------------------------------
# lmRatio <- lm(Protein ~ Gene:mRNA-1, data=modelDat)
# summary(lmRatio)$r.squared
# head(lmRatio$coefficients)
# length(lmRatio$coefficients)

# INTERCEPT only ---------------------------------------------------------------
lmIntercept <- lm(Protein ~ Gene, data=modelDat)
summary(lmIntercept)$r.squared
head(lmIntercept$coefficients)
length(lmIntercept$coefficients)

# FULL MODEL --------------------------------------------------------------
lmFull <- lm(Protein ~ Gene*mRNA, data=modelDat)
summary(lmFull)$r.squared
head(lmFull$coefficients)
length(lmFull$coefficients)

# ANOVA -------------------------------------------------------------------
# this is what they should have done:
# anova(lmRatio, lmIntercept)
# anova(lmIntercept, lmRatio)
# this is our way:
anova(lmFull, lmIntercept)
anova(lmIntercept, lmFull)
