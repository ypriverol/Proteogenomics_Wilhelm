setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

(load("2_mats.RData"))

(names(mats) <- gsub("LOO", "zLOO", gsub("NI", "zNI", gsub("TFs", "zTFs", names(mats)))))
lapply(mats, str)
sapply(is.na(mats), sum)
sapply(mats, function(m) sum(m==0))

library(reshape)
library(ggplot2)

# CORRELATIONS ------------------------------------------------------------
library(robustbase)
corMats <- list()
for(matNam in names(mats[-c(2)])){
  rownames(mats[[matNam]]) <- proteins
  colnames(mats[[matNam]]) <- tissues
  corMats[[matNam]] <- matrix(NA, length(tissues), length(tissues))
  colnames(corMats[[matNam]]) <- tissues
  rownames(corMats[[matNam]]) <- tissues
  for(t1 in tissues){
    for(t2 in tissues){
      #       corMats[[matNam]][t1, t2] <- cor(mats$prots[,t1], mats[[matNam]][,t2], method="spearman", use="pairwise.complete.obs")
      try(corMats[[matNam]][t1, t2] <- covMcd(data.frame(mats$prots[,t1], mats[[matNam]][,t2]) , cor=T)$cor[1,2])
    }
  }
}

library(gdata)
diags <- lapply(corMats, diag)
names(diags) <- paste0(names(diags), "_diag")
rest <- lapply(corMats, function(mat){
  diag(mat) <- NA
  return(as.vector(mat))
})
names(rest) <- paste0(names(rest), "_rest")
library(reshape)
library(ggplot2)
str(pDat <- melt.list(c(diags, rest)))
pDat$data <- factor(pDat$L1)
pDat$method <- factor(sapply(strsplit(pDat$L1, "_"),  function(v) return(v[[1]][1])))
ggplot(pDat, aes(x=data, y=value)) + geom_boxplot(aes(colour=method)) + theme_bw(24) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("3_Correlation.pdf", height=8)