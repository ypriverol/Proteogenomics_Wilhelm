source("../SETUP.R")
setwd("final/")

library(reshape)
library(ggplot2)
library(gdata)

(load(file="3_Mats_norm.RData"))
label = "_norm"


# QUICK LOOK --------------------------------------------------------------
names(mats)
lapply(mats, str)
sapply(mats, function(m) sum(is.na(m), na.rm=T))
table(apply(is.na(mats$prots), 1, sum))
77*12
77*12+248

# CORRELATIONS IN SAMPLES ACROSS PROTEINS------------------------------------------------------------
corMats <- list() # this is where the correlations are stored finally
for(matNam in names(mats)){
  rownames(mats[[matNam]]) <- proteins
  colnames(mats[[matNam]]) <- tissues
  corMats[[matNam]] <- sapply(tissues, function(t1){cor(mats$prots[,t1], mats[[matNam]][,t1], method="spearman", use="pairwise.complete.obs")})
}
# write to file
write.csv(do.call(cbind, corMats), file=paste0("4_AcrossGenes",label,".csv"))


# MAKE Boxplots ACROSS GENES -------------------------------------------------------------
matriXes <- c("mRNAs", "Ratio", "RatioLOO")
names(corMats)
sapply(corMats, length)
pDat <- melt.list(corMats[matriXes])
pDat$predictionMat <- factor(pDat$L1, levels=names(corMats)) # sorting of names
ggplot(pDat, aes(x=predictionMat, y=value)) + theme_bw(24) + 
  xlab("Method") + ylab("Correlation (R)") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  geom_boxplot(fill="lightgrey", colour="grey", outlier.colour = "grey") + 
  geom_jitter(position = position_jitter(width = .1))
ggsave(paste0("4_AcrossGenes",label,".pdf"), width=5, height=5)


# # CORRELATIONS BY PROTEIN ACROSS SAMPLES ------------------------------------
# # pCorMats <- list()
# # for(matNam in c("mRNAs", "Ratio")){
# #   pCorMats[[matNam]] <- sapply(1:length(proteins), function(p.idx){cor(mats$prots[p.idx,], mats[[matNam]][p.idx,], method="spearman", use="pairwise.complete.obs")})
# # }
# # sapply(pCorMats, function(vec) sum(!is.na(vec)))
# # str(pDat <- melt(pCorMats))
# # ggplot(pDat, aes(x=L1, y=value)) + geom_violin(fill = "lightgrey") + 
# #   xlab("Method") + ylab("Correlation (R)") + 
# #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
# #   theme_bw(24)
# # ggsave(paste0("4_AcrossSamples",label,".pdf"), width=5, height=5)
# 
# 
# # OVERFIT EFFECT ----------------------------------------------------------
# xx <- data.frame(
#   claimed_MSE = apply((mats$prots-mats$Ratio)**2, 1, function(row){return(mean(row, na.rm=T))}),
#   overfit_cor_MSE = apply((mats$prots-mats$RatioLOO)**2, 1, function(row){return(mean(row, na.rm=T))}),
#   cnt_prot = apply(mats$prots, 1, function(r) sum(!is.na(r))),
#   cnt_gene = apply(mats$mRNAs, 1, function(r) sum(!is.na(r))),
#   cnt = apply(!is.na(mats$mRNAs) & !is.na(mats$prots), 1 ,sum)
# )
# quantile(xx$overfit_effect <- (xx$overfit_cor_MSE/xx$claimed_MSE), na.rm=T)
# 
# ggplot(xx, aes(y=overfit_effect, group=cnt, x=factor(cnt))) + geom_boxplot() + 
#   scale_y_log10() + 
#   xlab("Available values per gene") + ylab("Overfit") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme_bw(24)
# ggsave(paste0("4_OverfitEffect",label,".pdf"), height=5, width=7)
# quantile(subset(xx, cnt==2)$overfit_effect, na.rm=T)
# head(subset(xx, cnt==2, select=c(claimed_MSE, overfit_cor_MSE, overfit_effect)))
# sapply(split(xx$overfit_effect, factor(xx$cnt)), function(v) return(median(v, na.rm=T)))
# 
# 
# # MSE ----------------------------------------------------------
# names(mats)
# mses <- lapply(mats, function(m){
#   x <- (m - mats$prots)**2
#   return(apply(x, 1, mean, na.rm=T))
# })
# str(mses)
# names(mses)
# # str(msesDif <- lapply(mses[1:3], function(v){
# #   return(log2(v/mses$intercept_LOO))
# # }))
# 
# matriXes2 <- c("Ratio", "Control")
# str(msesPlot <- melt(mses[matriXes2]))
# ggplot(msesPlot, aes(x=L1, y=value)) + geom_boxplot(fill="lightgrey") + 
#   scale_y_log10() + theme_bw(24) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ylab("MSE") + xlab("Method")
# ggsave(paste0("4_MSE",label,".pdf"), height=5, width=5)
# 
# ggplot(data.frame(do.call(cbind, mses)), aes(y=Ratio, x=Control)) + 
#   geom_point() + 
#   scale_y_log10() + scale_x_log10() + 
#   theme_bw(24) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ylab("Ratio") + xlab("Control")
# 
# matriXes3 <- c("RatioLOO", "ControlLOO")
# str(msesPlot <- melt(mses[matriXes3]))
# ggplot(msesPlot, aes(x=L1, y=value)) + geom_violin(fill="lightgrey") + 
#   scale_y_log10() + theme_bw(24) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   ylab("MSE") + xlab("Method")
# ggsave(paste0("4_MSE_LOO",label,".pdf"), height=5, width=5)
