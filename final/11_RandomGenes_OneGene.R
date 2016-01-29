library(gdata)
library(ggplot2)
library(reshape)
library(reshape2)

load("3_Mats.RData")

genes <- proteins
table(values.av <- sapply(genes, function(p){
  return(sum(!is.na(mats$prots[p,]) & !is.na(mats$mRNA[p,])))
}))
# ggplot(data.frame(Values=factor(values.av)), aes(x=Values)) + geom_histogram() + theme_bw(24)

names(values.av) <- genes

str(aDat <- data.frame(Values = values.av[genes]))
rownames(aDat) <- genes

set.seed(18828)

# PREDICT PROTEINS FROM ALL OTHER GENES (Average after) -------------------------------------------------------------------
library(doMC)
gns.idces <- 1:length(proteins)
randomVecs <- foreach(p.idx=1:length(proteins)) %do% {
  # make prediction from one gene
  g.idx <- sample(gns.idces[-p.idx],1)
  # get ratio for this gene
  ratio <- median(mats$prots[p.idx, ]/mats$mRNAs[g.idx, ], na.rm=T)
  # predict proteins
  return(mats$mRNAs[g.idx, ] * ratio)
}

str(randomMat <- do.call(rbind, randomVecs))
save(randomMat, file="11_RandomMat.RData")

# ORGANIZE THIS FOR COMPARISON --------------------------------------------
mats2 <- list(
  Original = mats$Ratio,
  Random = randomMat
)
for(matNam in names(mats2)){
  rownames(mats2[[matNam]]) <- proteins
  colnames(mats2[[matNam]]) <- tissues
}

# CORRELATION ACROSS PROTEINS PER TISSUE ----------------------------------
corTissues <- list() # this is where the correlations are stored finally
for(matNam in names(mats2)){
  corTissues[[matNam]] <- sapply(tissues, function(t1){cor(mats$prots[,t1], mats2[[matNam]][,t1], method="spearman", use="pairwise.complete.obs")})
}
# BOXPLOT
matriXes <- c("Original", "Random")
names(corTissues)
sapply(corTissues, length)
pDat <- melt.list(corTissues[matriXes])
pDat$predictionMat <- factor(pDat$L1, levels=names(corTissues)) # sorting of names
ggplot(pDat, aes(x=predictionMat, y=value)) + theme_bw(24) + 
  xlab("Method") + ylab("Correlation (R)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_boxplot(fill="lightgrey", colour="grey", outlier.colour = "grey") + 
  geom_jitter(position = position_jitter(width = .1))
ggsave(paste0("10_Random_AcrossGenesPerTissue.pdf"), width=5, height=5)



# CORRELATING ACROSS TISSUES PER GENE  -------------------
corGenes <- list()
for(matNam in names(mats2)){
  corGenes[[matNam]] <- sapply(proteins, function(p1){cor(mats$prots[p1,], mats2[[matNam]][p1,], method="spearman", use="pairwise.complete.obs")})
}
# BOXPLOT
matriXes <- c("Original", "Random")
names(corGenes)
sapply(corGenes, length)
str(pDat <- melt.list(corGenes[matriXes]))
pDat$predictionMat <- factor(pDat$L1, levels=matriXes) # sorting of names
ggplot(pDat, aes(x=value, colour=predictionMat, group=predictionMat)) + theme_bw(24) + 
  xlab("Correlation (R)") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_density(fill=NA, size=2)
ggsave(paste0("10_Random_AcrossTissuesPerGene.pdf"), width=5, height=5)
