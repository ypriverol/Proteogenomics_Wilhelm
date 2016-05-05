setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

# evaluating the predictions by getting correlation (across proteins per sample as done in the original publication)
# also evaluating the effect of missing data and variability

library(reshape)
library(ggplot2)

(load(file="5_mats.RData"))
names(mats)

# QUICK LOOK --------------------------------------------------------------
# (names(mats) <- gsub("LOO", "zLOO", gsub("NI", "zNI", gsub("TFs", "zTFs", names(mats)))))
lapply(mats, str)
sapply(is.na(mats), sum)
sapply(mats, function(m) sum(m==0, na.rm=T))

# CORRELATIONS IN SAMPLES ACROSS PROTEINS------------------------------------------------------------
corMats <- list() # this is where the correlations are stored finally
for(matNam in names(mats)[c(1,3,4)]){
  rownames(mats[[matNam]]) <- proteins
  colnames(mats[[matNam]]) <- tissues
  corMats[[matNam]] <- sapply(tissues, function(t1){cor(mats$prots[,t1], mats[[matNam]][,t1], method="spearman", use="pairwise.complete.obs")})
}

# Predicting just with intercept 
predMat <- matrix(NA, nrow=length(proteins), ncol=length(tissues))
for(t.idx in 1:length(tissues)){
  predMat[,t.idx] <- apply(mats$prots[, -t.idx], 1, mean, na.rm=T)
}
corMats[["intercept_LOO"]] <- sapply(1:length(tissues), function(t1){cor(mats$prots[,t1], predMat[,t1], method="spearman", use="pairwise.complete.obs")})

# Correlation by number of missing samples ------------------------
# for(proteinsPresent in c(3, 6, 9, 12)){
#   str(proteins2 <- rownames(mats$mRNAs)[which(apply(!is.na(mats$mRNAs), 1, sum) >= proteinsPresent & apply(!is.na(mats$prots), 1, sum) >= proteinsPresent)])
#   predMat <- matrix(NA, nrow=length(proteins2), ncol=length(tissues))
#   for(p.idx in 1:length(proteins2)){
#     for(t.idx in 1:length(tissues)){
#       ratio <- median(mats$prots[proteins2[p.idx], -t.idx]/mats$mRNAs[proteins2[p.idx], -t.idx], na.rm=T)
#       predMat[p.idx, t.idx] <- mats$mRNAs[proteins2[p.idx], t.idx] * ratio
#     }
#   }
#   corMats[[paste0("ratio_LOO_", proteinsPresent)]] <- sapply(1:length(tissues), function(t1){cor(mats$prots[proteins2,t1], predMat[,t1], method="spearman", use="pairwise.complete.obs")})  
# }


# GROUPS BY VARIANCE (norm to mean) ---------------------------------------
# variance is normalized to the mean expression ** 2 otherwise i get a strange effect.
# averageProt <- apply(mats$prots, 1, function(row){mean(row, na.rm=T)})
# varProt <- apply(mats$prots, 1, function(row){var(row, na.rm=T)})/averageProt**2
# subGroups <- split(proteins, cut(varProt, breaks=quantile(varProt, na.rm=T)))
# sapply(subGroups, length)
# names(subGroups) <- c("low", "midLow", "midHigh", "high")

# predicting proteins by LOO ratio and getting correlation for each sample
# now that I look at this, it looks like i'm repredicting protein, by ratio_LOO which is not necessarys
# for(subGroupNam in names(subGroups)){
#   proteins2 <- subGroups[[subGroupNam]]
#   predMat <- matrix(NA, nrow=length(proteins2), ncol=length(tissues))
#   for(p.idx in 1:length(proteins2)){
#     for(t.idx in 1:length(tissues)){
#       ratio <- median(mats$prots[proteins2[p.idx], -t.idx]/mats$mRNAs[proteins2[p.idx], -t.idx], na.rm=T)
#       predMat[p.idx, t.idx] <- mats$mRNAs[proteins2[p.idx], t.idx] * ratio
#     }
#   }
#   corMats[[paste0("ratio_LOO_", subGroupNam, "Var")]] <- sapply(1:length(tissues), function(t1){cor(mats$prots[proteins2,t1], predMat[,t1], method="spearman", use="pairwise.complete.obs")})
# }




# MAKE FIGURE -------------------------------------------------------------
names(corMats[1:3])
sapply(corMats[1:3], length)
pDat <- melt.list(corMats[1:3])
pDat$predictionMat <- factor(pDat$L1, levels=names(corMats)) # sorting of names
ggplot(pDat, aes(x=predictionMat, y=value)) + theme_bw(24) + 
  ylab("Spearman correlation") + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_boxplot(fill="lightgrey", colour="grey", outlier.colour = "grey") + geom_jitter(position = position_jitter(width = .1))
ggsave("5_1_FinalFigure_Normalized.pdf", width=7, height=5)



# CORRELATIONS BY PROTEIN ACROSS SAMPLES ------------------------------------
pCorMats <- list()
for(matNam in names(mats)[c(1,3,4)]){
  pCorMats[[matNam]] <- sapply(1:length(proteins), function(p.idx){cor(mats$prots[p.idx,], mats[[matNam]][p.idx,], method="spearman", use="pairwise.complete.obs")})
}
sapply(pCorMats, function(vec) sum(!is.na(vec)))
str(pDat <- melt(pCorMats))
ggplot(pDat, aes(x=L1, y=value)) + geom_violin(fill = "lightgrey") + 
  xlab("Method") + ylab("Spearman correlation") + 
  theme_bw(24)
ggsave("5_2_AcrossSamples.pdf", width=7, height=5)


