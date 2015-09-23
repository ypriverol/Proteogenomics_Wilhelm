setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

# Evaluation of the predicted proteins by mean and median squared error

(load("2_mats.RData"))

(names(mats) <- gsub("LOO", "zLOO", gsub("NI", "zNI", gsub("TFs", "zTFs", names(mats)))))
lapply(mats, str)
sapply(is.na(mats), sum)
sapply(mats, function(m) sum(m==0))

library(reshape)
library(ggplot2)

# MSE by Tissue--------------------------------------------------------------

# MEAN
mses <- lapply(mats, function(m) return((m-mats$prot)**2))
tissueMeanSEs <- lapply(mses, function(m) return(apply(m, 2, mean)))
str(pDat <- melt.list(tissueMeanSEs))
pDat$data <- factor(pDat$L1)
table(pDat$L1)
ggplot(pDat, aes(x=L1, y=value)) + geom_boxplot(aes(colour = L1)) + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Mean squared error") + scale_y_log10()
ggsave("3_MeanSquaredError.pdf", height=8)

# MEDIAN
tissueMedianSEs <- lapply(mses, function(m) return(apply(m, 2, median)))
str(pDat <- melt.list(tissueMedianSEs))
pDat$data <- factor(pDat$L1)
table(pDat$L1)
ggplot(pDat, aes(x=L1, y=value)) + geom_boxplot(aes(colour = L1)) + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Median squared error") + scale_y_log10()
ggsave("3_MedianSquaredError.pdf", height=8)

# MSE by protein--------------------------------------------------------------

# MEAN
proteinMeanSEs <- lapply(mses, function(m) return(apply(m, 1, mean)))
str(pDat <- melt.list(proteinMeanSEs))
pDat$data <- factor(pDat$L1)
table(pDat$L1)
ggplot(pDat, aes(x=L1, y=value)) + geom_boxplot(aes(colour = L1)) + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Mean squared error") + scale_y_log10()
ggsave("3_MeanSquaredError_Proteins.pdf", height=8)

