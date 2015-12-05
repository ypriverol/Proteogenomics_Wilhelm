setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

library(reshape)
library(ggplot2)

(load(file="4_mats.RData"))
names(mats)

# CORRELATIONS IN SAMPLES ACROSS PROTEINS------------------------------------------------------------
stopifnot(all(rownames(mats$prots) == proteins))
str(pCors <- sapply(1:length(proteins), function(p.idx){cor(mats$prots[p.idx,], mats$mRNAs[p.idx,], method="spearman", use="pairwise.complete.obs")}))
hist(pCors)

allExamples <- proteins[which(pCors > -0.01 & pCors < 0.01)]

rownames(mats$ratio_LOO) = rownames(mats$prots)
colnames(mats$ratio_LOO) = colnames(mats$prots)

str(pDat <- data.frame(
  melt(mats$prots[allExamples,]),
  mRNA = melt(mats$mRNAs[allExamples,])$value,
  ProteinPrediction= melt(mats$ratio_LOO[allExamples,])$value
  ))
colnames(pDat)[1:3] <- c("Gene", "Tissue", "Protein")
pDat$kidney <- pDat$Tissue == "kidney"
head(pDat)


# Three examples
examplesProteins <- proteins[which(pCors > 0 & pCors < 0.01)[c(1,3,4,7)]]
examplesProteins %in% allExamples
str(pDat2 <- drop.levels(subset(pDat, Gene %in% examplesProteins)))
ggplot(pDat2, aes(y=Protein, colour=Gene)) + #geom_point(aes(x=mRNA)) + 
  geom_abline(slope=1, size = 2, colour = "grey") +
  geom_segment(aes(x=mRNA, xend=ProteinPrediction, y=Protein, yend=Protein, colour=Gene), arrow=arrow(angle=15, length = unit(0.15, "inches")), size = 0.8) + 
  theme_bw(24) +
  scale_x_log10() + scale_y_log10() + 
  xlab("mRNA to Predicted Protein")
ggsave("4_3_Example.pdf", width=9, height=5)

with(subset(pDat2, Tissue=="kidney"), cor(mRNA, Protein, method="spearman", use="pairwise.complete.obs"))
with(subset(pDat2, Tissue=="kidney"), cor(ProteinPrediction, Protein, method="spearman", use="pairwise.complete.obs"))
ggplot(pDat2, aes(y=Protein)) + 
  geom_point(aes(x=mRNA, colour=Gene), shape = 1) + 
  geom_point(aes(x=ProteinPrediction, colour=Gene), shape = 16) + 
  geom_abline(slope=1, size = 2, colour = "grey") +
  geom_segment(data = subset(pDat2, Tissue=="kidney"), 
               aes(x=mRNA, xend=ProteinPrediction, y=Protein, yend=Protein, colour=Gene), 
               arrow=arrow(angle=25, length = unit(0.15, "inches")), size = 0.8) + 
  theme_bw(24) +
  scale_x_log10() + scale_y_log10() + 
  xlab("mRNA to Predicted Protein")
ggsave("4_3_Example.pdf", width=9, height=5)
