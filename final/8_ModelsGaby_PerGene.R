source("../0_SETUP.R")
setwd("final/")

library(gdata)
library(ggplot2)
library(reshape)
library(reshape2)

load("1_Mats.RData")

genes <- proteins

# Get NAs
table(values.av <- sapply(genes, function(p){
  return(sum(!is.na(mats$prots[p,]) & !is.na(mats$mRNA[p,])))
}))
ggplot(data.frame(Values=factor(values.av)), aes(x=Values)) + geom_histogram() + theme_bw(24)
ggsave("8_Values.pdf", width=7, height=5)

names(values.av) <- genes
str(genes <- genes[which(values.av >= 3)])

# ONE MODEL (slope only) FOR EACH GENE -------------------------------------------------------
pModels <- lapply(genes, function(p){ return(lm(mats$prots[p,] ~ mats$mRNAs[p,] - 1)) })
names(pModels) <- genes

# get R2
r2 <- sapply(pModels[which(!is.na(pModels))], function(m){
  return(summary(m)$r.squared)
})

# show one gene
g1 <- names(r2[which(r2 == max(r2))])
plot(mats$mRNA[g1,], mats$prots[g1,])

# SPearman correlation
spearCor <- sapply(genes, function(p){ return(cor(mats$prots[p,], mats$mRNAs[p,], use="pairwise.complete.obs", method="spearman"))})
names(spearCor) <- genes

str(pDat <- data.frame(R2 = r2, Values = factor(values.av[genes]), Correlation=spearCor))

ggplot(data=pDat, aes(x=Values, y=Correlation)) + geom_boxplot() + theme_bw(24)
ggsave("8_Correlation.pdf", height=5, width=7)
round(sapply(split(pDat$Correlation, pDat$Values), median),2)

ggplot(data=pDat, aes(x=Values, y=R2)) + geom_boxplot() + theme_bw(24)
ggsave("8_R2.pdf", height=5, width =7 )
round(sapply(split(pDat$R2, pDat$Values), median),2)