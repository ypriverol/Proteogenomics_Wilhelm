library(gdata)
library(ggplot2)
library(reshape)
library(reshape2)
library(plyr)

load("3_Mats.RData")

g <- c("ENSG00000102309", # Peptidyl-prolyl cis-trans isomerase NIMA-interacting 4
"ENSG00000072042",        # Retinol dehydrogenase 11
"ENSG00000184047",        # DIABLO
"ENSG00000143222",        # Ubiquitin-fold modifier-conjugating enzyme 1
"ENSG00000115241"         # Protein phosphatase 1 (?)
)

str(pDat <- merge(melt(mats$prots[g,]), melt(mats$mRNAs[g,]), by=c("X1", "X2")))
colnames(pDat) <- c("Gene", "Tissue", "Protein", "mRNA")
pDat$Ratio <- pDat$mRNA/pDat$Protein

ggplot(pDat, aes(x=Gene, y=Ratio)) + geom_point() + 
  theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pDat2 <- merge(melt(mats$prots), melt(mats$mRNAs), by=c("X1", "X2"))
colnames(pDat2) <- c("Gene", "Tissue", "Protein", "mRNA")
pDat2$Ratio <- pDat2$mRNA/pDat2$Protein
str(pDat3 <- ddply(pDat2, .(Gene), summarize, sd_ratio = sd(Ratio, na.rm=T), mean_ratio=mean(Ratio, na.rm=T)))
ggplot(pDat3, aes(x=sd_ratio)) + geom_histogram() + xlim(0, 25)
quantile(pDat3$sd_ratio, na.rm=T)

ggplot(pDat3, aes(x=log2(sd_ratio/mean_ratio))) + geom_histogram() + xlim(-4, 4)
quantile(pDat3$sd_ratio/pDat3$mean_ratio, na.rm=T)
quantile(log2(pDat3$sd_ratio/pDat3$mean_ratio), na.rm=T)
ggplot(pDat3, aes(x=sd_ratio/mean_ratio)) + stat_ecdf() + theme_bw(24) + 
  xlab("Relative standard deviation") + ylab("Fraction")
ggsave("9_relSD_ecdf.pdf", width = 5, height = 5)

ggplot(pDat3, aes(x=sd_ratio)) + stat_ecdf() + xlim(0, 25)
quantile(pDat3$sd_ratio, na.rm=T)
sum(pDat3$sd_ratio > 2, na.rm=T)


ggplot(subset(pDat, Gene=="ENSG00000102309"), aes(x=mRNA, y=Protein)) + geom_point() + theme_bw(24)
cor(subset(pDat, Gene=="ENSG00000102309", select = c(Protein, mRNA)), 
    use='pairwise.complete.obs', method="spearman")
subset(pDat3, Gene=="ENSG00000102309")
