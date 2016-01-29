library(reshape)
library(ggplot2)
library(gdata)

(load(file="3_Mats.RData"))
label = ""
names(mats)
mats$RandomGenes <- randomMat

# QUICK LOOK --------------------------------------------------------------
names(mats)
lapply(mats, str)
sapply(mats, function(m) sum(is.na(m), na.rm=T))
table(apply(is.na(mats$prots), 1, sum))

# CORRELATIONS IN SAMPLES ACROSS PROTEINS------------------------------------------------------------
load("11_RandomMat.RData")
mats$rGenes <- randomMat
corMats <- list() # this is where the correlations are stored finally
for(matNam in names(mats)){
  rownames(mats[[matNam]]) <- proteins
  colnames(mats[[matNam]]) <- tissues
  corMats[[matNam]] <- sapply(tissues, function(t1){cor(mats$prots[,t1], mats[[matNam]][,t1], method="spearman", use="pairwise.complete.obs")})
}
# write to file
write.csv(do.call(cbind, corMats), file=paste0("4_AcrossGenes",label,".csv"))
round(sapply(corMats, median),2)


# MAKE Boxplots ACROSS GENES PER TISSUE  -------------------------------------------------------------
matriXes <- c("mRNAs", "Ratio", "Control", "rGenes")
names(corMats)
sapply(corMats, length)
pDat <- melt.list(corMats[matriXes])
pDat$predictionMat <- factor(pDat$L1, levels=names(corMats)) # sorting of names
ggplot(pDat, aes(x=predictionMat, y=value)) + theme_bw(24) + 
  xlab("Method") + ylab("Correlation\nacross genes") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_boxplot(fill="lightgrey", colour="grey", outlier.colour = "grey") + 
  geom_jitter(position = position_jitter(width = .1))
ggsave(paste0("4_AcrossGenes",label,".pdf"), width=5, height=5)


# The same with OVERFITTING  -------------------------------------------------------------
matriXes2 <- c("mRNAs", "Ratio", "RatioLOO", "Control", "ControlLOO")
names(corMats)
sapply(corMats, length)
pDat <- melt.list(corMats[matriXes2])
pDat$predictionMat <- factor(pDat$L1, levels=matriXes2) # sorting of names
ggplot(pDat, aes(x=predictionMat, y=value)) + theme_bw(24) + 
  xlab("") + ylab("Correlation\nacross genes") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_boxplot(fill="lightgrey", colour="grey", outlier.colour = "grey") + 
  geom_jitter(position = position_jitter(width = .1))
ggsave(paste0("4_AcrossGenes_OverfitCor",label,".pdf"), width=7, height=5)
