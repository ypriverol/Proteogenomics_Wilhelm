source("../0_SETUP.R")
setwd("final/")

library(reshape)
library(ggplot2)
library(gdata)
library(grid)

(load(file="3_Mats.RData"))
names(mats)

# CORRELATIONS IN SAMPLES ACROSS PROTEINS------------------------------------------------------------
stopifnot(all(rownames(mats$prots) == proteins))
str(pCors <- sapply(1:length(proteins), function(p.idx){cor(mats$prots[p.idx,], mats$mRNAs[p.idx,], method="spearman", use="pairwise.complete.obs")}))
hist(pCors)

allExamples <- proteins[which(pCors > -0.01 & pCors < 0.01)]

rownames(mats$Ratio) = rownames(mats$prots)
colnames(mats$Ratio) = colnames(mats$prots)

str(pDat <- data.frame(
  melt(mats$prots[allExamples,]),
  mRNA = melt(mats$mRNAs[allExamples,])$value,
  Ratio= melt(mats$Ratio[allExamples,])$value
  ))
colnames(pDat)[1:3] <- c("Gene", "Tissue", "Protein")
pDat$kidney <- pDat$Tissue == "kidney"
head(pDat)


# Three examples
examplesProteins <- proteins[which(pCors > 0 & pCors < 0.01)[c(1,3,4,7)]]
examplesProteins %in% allExamples
str(pDat2 <- drop.levels(subset(pDat, Gene %in% examplesProteins)))
# ggplot(pDat2, aes(y=Protein, colour=Gene)) + #geom_point(aes(x=mRNA)) + 
#   geom_abline(slope=1, size = 2, colour = "grey") +
#   geom_segment(aes(x=mRNA, xend=Ratio, y=Protein, yend=Protein, colour=Gene), arrow=arrow(angle=15, length = unit(0.15, "inches")), size = 0.8) + 
#   theme_bw(24) +
#   scale_x_log10() + scale_y_log10() + 
#   xlab("mRNA to predicted protein") + ylab("Measured protein")
# ggsave("5_Example.pdf", width=9, height=5)

with(subset(pDat2, Tissue=="kidney"), cor(mRNA, Protein, method="spearman", use="pairwise.complete.obs"))
with(subset(pDat2, Tissue=="kidney"), cor(Ratio, Protein, method="spearman", use="pairwise.complete.obs"))
ggplot(pDat2, aes(y=Protein)) + 
  geom_point(aes(x=mRNA, colour=Gene), shape = 1) + 
  geom_point(aes(x=Ratio, colour=Gene), shape = 16) + 
  geom_abline(slope=1, size = 2, colour = "grey") +
  geom_segment(data = subset(pDat2, Tissue=="kidney"), 
               aes(x=mRNA, xend=Ratio, y=Protein, yend=Protein, colour=Gene), 
               arrow=arrow(angle=25, length = unit(0.15, "inches")), size = 0.8) + 
  theme_bw(24) +
  scale_x_log10() + scale_y_log10() + 
  theme(legend.position="none") + 
  xlab("mRNA / predicted protein") + ylab("Observed protein")
ggsave("5_Example.pdf", width=5, height=5)
