setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

(load("0_1_data.RData"))

library(reshape)
library(ggplot2)

# EXPLORE THESE DATA IN RELATIONSHIP --------------------------------------
average <- apply(protUN, 1, function(row){mean(row, na.rm=T)})
xx <- data.frame(
  claimed_MSE = apply((protUN-protPred)**2, 1, function(row){return(mean(row, na.rm=T))})/average**2,
  overfit_cor_MSE = apply((protUN-nikProtPred)**2, 1, function(row){return(mean(row, na.rm=T))})/average**2,
  var_prot = apply(protUN, 1, function(row){var(row, na.rm=T)})/average**2,
#   var_gene = apply(geneUN, 1, function(row){var(row, na.rm=T)})/average,
  cnt_prot = apply(protUN, 1, function(r) sum(!is.na(r))),
  cnt_gene = apply(geneUN, 1, function(r) sum(!is.na(r))),
  cnt = apply(!is.na(geneUN) & !is.na(protUN), 1 ,sum)
)
quantile(xx$overfit_effect <- (xx$overfit_cor_MSE/xx$claimed_MSE), na.rm=T)
table(xx$cnt_prot, xx$cnt_gene)
round(cor(xx, use="pairwise.complete.obs", method="spearman"),2)
heatmap(round(cor(xx, use="pairwise.complete.obs", method="spearman"),2), scale="none")



apply(is.na(xx), 2,sum)
apply(xx, 2, function(col) quantile(col, na.rm=T))
# ggplot(melt(xx), aes(x=value, colour = variable)) + stat_ecdf()
# ggplot(melt(xx), aes(x=value, colour = variable)) + stat_ecdf() + xlim(c(0,3))
# ggplot(melt(xx), aes(x=value, colour = variable)) + stat_ecdf() + scale_y_log10()

library(hexbin)
hexplom(xx)

# average
ggplot(xx, aes(y=claimed_MSE, x=average)) + geom_point() + theme_bw(24) + scale_y_log10() + scale_x_log10()
ggplot(xx, aes(y=overfit_cor_MSE, x=average)) + geom_point() + theme_bw(24) + scale_y_log10() + scale_x_log10()
ggplot(xx, aes(y=var_prot, x=average)) + geom_point() + theme_bw(24) + scale_y_log10() + scale_x_log10()

# count
table(xx$cnt)
ggplot(xx, aes(y=claimed_MSE, group=cnt, x=factor(cnt))) + geom_boxplot() + theme_bw(24)
ggplot(xx, aes(y=claimed_MSE, group=cnt, x=factor(cnt))) + geom_boxplot() + theme_bw(24) + scale_y_log10()
dev.print(pdf, "0_0_1_Count_ClaimedMSE.pdf")
ggplot(xx, aes(y=overfit_cor_MSE, group=cnt, x=factor(cnt))) + geom_boxplot() + theme_bw(24) + scale_y_log10()
dev.print(pdf, "0_0_1_Count_OverfitCorMSE.pdf")
ggplot(xx, aes(y=overfit_effect, group=cnt, x=factor(cnt))) + geom_boxplot() + theme_bw(24) + scale_y_log10()
quantile(subset(xx, cnt==2)$overfit_effect, na.rm=T)
head(subset(xx, cnt==2, select=c(claimed_MSE, overfit_cor_MSE, overfit_effect)))
sapply(split(xx$overfit_effect, factor(xx$cnt)), function(v) return(median(v, na.rm=T)))

# why are there non zeros at cnt==1
head(nonZero <- subset(xx, cnt==1 & claimed_MSE != 0, select = claimed_MSE))
(x <- rownames(nonZero)[1])
geneUN[x,]
protUN[x,]
(ratios <- protUN[x,]/geneUN[x,])
(medianRatio <- median(ratios, na.rm = T))
medianRatio * geneUN[x,]
protUN[x,]
medianRatio * geneUN[x,]-protUN[x,]

# variability
# low variability has lower MSE are lower overfit effect
ggplot(xx, aes(y=claimed_MSE, group=var_prot, x=var_prot)) + geom_point(alpha=0.3, colour="blue") + theme_bw(24)+ scale_y_log10() + scale_x_log10()
ggplot(xx, aes(y=overfit_cor_MSE, group=var_prot, x=var_prot)) + geom_point(alpha=0.3, colour="blue") + theme_bw(24) + scale_y_log10() + scale_x_log10()
ggplot(xx, aes(y=overfit_effect, group=var_prot, x=var_prot)) + geom_point(alpha=0.3, colour="blue") + theme_bw(24) + scale_y_log10() + scale_x_log10()


# MAKE GROUPS BASED ON VARIABILITY ----------------------------------------
plot(density(xx$var_prot, na.rm=T))
subGroups <- split(rownames(xx), cut(xx$var_prot, breaks=quantile(xx$var_prot, na.rm=T)))
# plot(density(log2(xx$var_prot), na.rm=T))
# logVar <- log2(xx$var_prot)
# subGroups <- split(rownames(xx), cut(logVar, floor(min(logVar, na.rm=T)):ceiling(max(logVar, na.rm=T))))
# sapply(subGroups, length)
# subGroups <- subGroups[which(sapply(subGroups, length) > 50)]
# sapply(subGroups, length)
# subGroups <- lapply(subGroups, function(v) return(sample(v, 50)))
sapply(subGroups, length)
names(subGroups) <- sub("^\\((\\-?\\d),(\\-?\\d)]", "\\1", names(subGroups))
subProtUN <- lapply(subGroups, function(p){return(protUN[p,])})
subGeneUN <- lapply(subGroups, function(p){return(geneUN[p,])})
subRatioMat <- list()
subProtPred <- list()
# predict proteins with leave one out
for(group in names(subGroups)){
  ratioVec <- list()
  for(tissue.idx in 1:length(tissues)){
    tissue <- tissues[tissue.idx]
    ratioVec[[tissue]] <- apply(subProtUN[[group]][,tissues[-tissue.idx]]/subGeneUN[[group]][,tissues[-tissue.idx]], 1, function(row) median(row, na.rm=T))
  }  
  subRatioMat[[group]] <- do.call(cbind, ratioVec)
  subProtPred[[group]] <- subGeneUN[[group]][,tissues] * subRatioMat[[group]][,tissues]
}
sapply(subRatioMat, dim)
# Get correlations
subCorrelations <- list()
for(group in names(subGroups)){
  subCorrelations[[group]] <- sapply(tissues, function(tissue){return(cor(subProtPred[[group]][,tissue], subProtUN[[group]][,tissue], method="spearman", use="pairwise.complete.obs"))})
}
str(pDat <- melt.list(subCorrelations))
pDat$log2Bin <- factor(pDat$L1, levels=unique(pDat$L1))
ggplot(pDat, aes(x=log2Bin, y=value)) + geom_boxplot() + theme_bw(24)
dev.print(pdf, "0_0_1_Variability_Correlation.pdf")

# WRONG Tissues
subCorrelationsWrongTissue <- list()
for(group in names(subGroups)){
  subCorrelationsWrongTissue[[group]] <- c()
  for(t1 in tissues){
    for(t2 in tissues){
      if(t1 != t2){
        subCorrelationsWrongTissue[[group]] <- c(subCorrelationsWrongTissue[[group]], cor(subProtPred[[group]][,t1], log10(subProtUN[[group]][,t2]), method="spearman", use="pairwise.complete.obs"))
      }
      
    }
  }  
}

names(subCorrelationsWrongTissue) <- paste0(c("1_low", "2_low_mid", "3_high_mid", "4_high"), "_rest")
names(subCorrelations) <- paste0(c("1_low", "2_low_mid", "3_high_mid", "4_high"), "_diag")
str(pDat <- melt.list(c(subCorrelations, subCorrelationsWrongTissue)))
table(pDat$L1)
ggplot(pDat, aes(x=L1, y=value)) + geom_boxplot() + theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.print(pdf, "0_0_1_Variability_Correlation_WrongTissue.pdf")