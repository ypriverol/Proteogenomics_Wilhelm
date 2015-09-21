setwd("~/OtherAnalysis/2015_07_03_mRNAvsProtein/DraftMapOfHumanProteome/")

(load("2_mats.RData"))

(names(mats) <- gsub("LOO", "zLOO", gsub("NI", "zNI", gsub("TFs", "zTFs", names(mats)))))
lapply(mats, str)
sapply(is.na(mats), sum)
sapply(mats, function(m) sum(m==0))

library(reshape)
library(ggplot2)

mses <- lapply(mats, function(m) return((m-mats$prot)**2))

# FUNCTION TO PLOT ONE PROTEIN
plotProtein <- function(p, pred){
  if(is.character(p)){p <- which(proteins == p)}
  pDat <- data.frame(
    tissue = tissues, 
    mRNA=mats$mRNAs[p,tissues], 
    protein=mats$prots[p,tissues])
  for(x in pred){
    pDat[,x] <- mats[[x]][p,]
  }
  plot <- ggplot(melt(pDat), aes(x=tissue, y=value, colour=variable, group=variable)) + geom_line() + theme_bw(24) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  return(plot)
}


# COMPARISON BETWEEN METHODS ----------------------------------------------

str(mses2 <- do.call(rbind, lapply(mses, function(mat) apply(mat, 1, mean))))
str(mses2 <- data.frame(t(mses2)))
cor(mses2, method="spearman")

# Linear model without intercept is pretty much like ratio!!
str(dif <- mses2$lm_zNI - mses2$ratio_zLOO)
plotProtein(which(dif == tail(sort(dif), n = 1))[1], c("lm_zNI", "ratio_zLOO")) + scale_y_log10()
plotProtein(which(dif == head(sort(dif), n = 1))[1], c("lm_zNI", "ratio_zLOO")) + scale_y_log10()

# Linear model is different from ratio but mostly because of the intercept!!
str(dif <- mses2$lm - mses2$ratio_zLOO)
plotProtein(which(dif == tail(sort(dif), n = 1))[1], c("lm", "ratio_zLOO")) + scale_y_log10()
plotProtein(which(dif == head(sort(dif), n = 1))[1], c("lm", "ratio_zLOO")) + scale_y_log10()

# same is true for the elastic net
str(dif <- mses2$elastic_zTFs_zNI - mses2$ratio_zLOO)
plotProtein(which(dif == tail(sort(dif), n = 1))[1], c("elastic_zTFs_zNI", "ratio_zLOO")) + scale_y_log10()
plotProtein(which(dif == head(sort(dif), n = 1))[1], c("elastic_zTFs_zNI", "ratio_zLOO")) + scale_y_log10()



# SOME EXAMPLES? --------------------------------------------------------

# HIGH MSE in RATIO
maxMSE <- apply(mses$ratio, 1, max)
tail(sort(maxMSE))
plotProtein(which(maxMSE %in% tail(sort(maxMSE), n=1)), c("lm", "ratio")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-1]), c("lm", "ratio")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-2]), c("lm", "ratio")) + scale_y_log10()

# HIGH MSE in LM
maxMSE <- apply(mses$lm, 1, max)
tail(sort(maxMSE))
plotProtein(which(maxMSE %in% tail(sort(maxMSE), n=1)), c("lm", "ratio")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-1]), c("lm", "ratio")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-2]), c("lm", "ratio")) + scale_y_log10()

# HIGH MSE in LM without zNI
maxMSE <- apply(mses$lm_zNI, 1, max)
tail(sort(maxMSE))
plotProtein(which(maxMSE %in% tail(sort(maxMSE), n=1)), c("lm", "ratio", "lm_zNI")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-1]), c("lm", "ratio", "lm_zNI")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-2]), c("lm", "ratio", "lm_zNI")) + scale_y_log10()

# HIGH MSE in elastic without zNI
maxMSE <- apply(mses$elastic_zTFs_zNI, 1, max)
tail(sort(maxMSE))
plotProtein(which(maxMSE %in% tail(sort(maxMSE), n=1)), c("lm", "ratio", "elastic_zTFs_zNI")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-1]), c("lm", "ratio", "elastic_zTFs_zNI")) + scale_y_log10()
plotProtein(which(maxMSE %in% sort(maxMSE)[length(maxMSE)-2]), c("lm", "ratio", "elastic_zTFs_zNI")) + scale_y_log10()



